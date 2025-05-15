#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import multiprocessing as mp
import concurrent.futures
import numpy as np
import editdistance
import time
import traceback
from functools import partial

from utils.basic import loadBlatOutput
from utils.coverage import average_edit_distance
from utils.align.maf_th import MultipleAlignmentReader
from utils.rna import RNAcofold2, containCommonSNP, containGquad2, countCpG, GCcontent

from jklib.genome import locus
from jklib.bioDB import CommonSNP


# ──────────────────────────────────────────────────────────────
# 1) 워커 초기화 & 쿼리 함수  ───────────────────────────────────
# ──────────────────────────────────────────────────────────────
_reader = None           # 워커 프로세스 안에서만 쓰는 전역 객체

def _init_worker(maf_dir, ref_assembly, query_assembly):
    """
    각 워커가 시작될 때 한 번만 실행.
    무거운 MultipleAlignmentReader 를 전역으로 만들어 둔다.
    """
    global _reader
    _reader = MultipleAlignmentReader(
        ref_assembly=ref_assembly,
        query_assembly=query_assembly,
        maf_dir=maf_dir,
    )

def _query_region(region: str):
    """
    워커가 실제로 수행할 일.
    region 하나를 받아서 reader.query() 결과를 반환.
    """
    return _reader.query(region, verbose=False)

# ──────────────────────────────────────────────────────────────
# 2) ASOdesign 클래스  ─────────────────────────────────────────
# ──────────────────────────────────────────────────────────────
class ASOdesign:
    def __init__(self,
                 transid="NM_002415",
                 refFlat_path="/Users/dowonkim/Dropbox/data/UCSC/hg38/refFlat/refFlat_200817.txt",
                 maf_dir='/Users/dowonkim/Dropbox/data/offtarget_test/maf',
                 query_assembly=["mm39"],      # tuple/리스트 허용
                 ref_assembly="hg38",
                 tile_length=17):
        print(f"[ASOdesign] transid={transid}")
        self.refFlat      = loadBlatOutput(refFlat_path, by='transID')
        self.cSNP = CommonSNP()
        self.maf_dir      = maf_dir
        self.query_asm    = query_assembly
        self.ref_asm      = ref_assembly
        self.transid      = transid
        self.tile_length  = tile_length

        self.transInfo    = self._get_transInfo()
        self.txn_tiles    = self._tile_region()
        self.txn_tile_seq = self._tile_txn_seq()

    # ───── 내부 헬퍼 ──────────────────────────────────────────
    def _get_transInfo(self):
        return [t for t in self.refFlat[self.transid]
                if len(t['chrom'].split('_')) == 1][0]

    def _tile_region(self):
        anti = '+' if self.transInfo['strand'] == '-' else '-'
        chrom = self.transInfo['chrom']
        start = self.transInfo['txnSta']
        end   = self.transInfo['txnEnd']
        L     = self.tile_length
        return [locus(f"{chrom}:{i}-{i+L-1}{anti}")   # inclusive end = i+L-1
            for i in range(start, end-L+1)]

    def _tile_txn_seq(self):
        anti = '+' if self.transInfo['strand'] == '-' else '-'
        chrom = self.transInfo['chrom']
        start = self.transInfo['txnSta']
        end   = self.transInfo['txnEnd']
        L     = self.tile_length
        txn_seq = locus(f"{chrom}:{start}-{end-1}{anti}").twoBitFrag().upper()
        result_array = [ txn_seq[i:i+L] for i in range(0, end-start-L+1)]
        if anti == '-':
            result_array = result_array[::-1]
        return result_array

    def _chunks(self, seq, n):
        """seq 를 n 등분하여 순차적으로 yield"""
        k = max(1, len(seq)//n)
        for i in range(0, len(seq), k):
            yield seq[i:i+k]

    # ───── 퍼블릭 메서드 ─────────────────────────────────────
    def process_main(self, chunk_division=3, max_workers=3):
        print(f"#tiles={len(self.txn_tiles)}, assemblies={self.query_asm}")

        
        result_dict = {}
        for asm in self.query_asm:                      # 어셈블리별로 별도 풀
            print(f"[{asm}]")
            with concurrent.futures.ProcessPoolExecutor(
                    max_workers=max_workers,
                    initializer=_init_worker,
                    initargs=(self.maf_dir, self.ref_asm, asm),
            ) as ex:
                all_results_locInfo = []
                all_results_maf, all_editdist = [], []
                for chunk_loc, chunk_seq in zip(self._chunks(self.txn_tiles, chunk_division), self._chunks(self.txn_tile_seq, chunk_division)):
                    # region 문자열만 워커에 전달
                    #chunk_locInfo = list(ex.map(self.getlocInfo, chunk))
                    chunk_locInfo = [self.getlocInfo(tile_loc, tile_seq) for tile_loc, tile_seq in zip(chunk_loc, chunk_seq)]
                    chunk_maf_results = list(ex.map(_query_region, chunk_loc))
                    dists   = [self._editdistance_safe(r) for r in chunk_maf_results]
                    all_results_locInfo.extend(chunk_locInfo)
                    all_results_maf.extend(chunk_maf_results)
                    all_editdist.extend(dists)
                result_dict[asm] = {"maf_seq": all_results_maf, "coverage": all_editdist}

        print("finished. first 3 results:", result_dict[self.query_asm[0]]["maf_seq"][:3])
        print("finished. first 3 results:", result_dict[self.query_asm[0]]["coverage"][:3])
        return result_dict

    @staticmethod
    def _editdistance_safe(res):
        try:
            human, maf = res.values()
            return average_edit_distance([human.upper(), maf.upper()])
            #return editdistance.eval(human.upper(), maf.upper())
        except Exception:
            return None
    
    
    def getlocInfo(self, loc, sequence):
        try:
            #sequence = loc.twoBitFrag().upper()
            locStr   = loc.toString()
            loc_tmp = locus(f'{loc.chrom}:{loc.chrSta+1}-{loc.chrEnd-1}{loc.strand}')
            flagL = [ e[3] for e in loc_tmp.regionType()]
            flag = []
            for a in flagL:
                #print(a[0][0])
                if (('cds' in a[0][0]) or ('utr' in a[0][0])): #exon flag ['cds','utr']
                    flag.append('e') 
                elif 'int' in a[0][0] : #intron
                    flag.append('i')

            regionT = '/'.join(map(str, loc.regionType()))
            flag = ':'.join(flag)
            homo, mono = RNAcofold2(sequence)
            #'type','gene','transcriptID','locus','sequence','length','regionType','commonSNP','Gquad','CpG','GC_content','Homo_dimer','Monomer'
            return [flag,
                    self.transInfo['transName'],
                    self.transInfo['transID'],
                    locStr, sequence, len(sequence), regionT,
                    containCommonSNP(loc, self.cSNP),
                    containGquad2(sequence),
                    countCpG(sequence),
                    GCcontent(sequence),
                    homo, mono]
        except Exception as e:
            print(traceback.format_exc())
            return e.args

# ──────────────────────────────────────────────────────────────
# 3) 실행 예시  ────────────────────────────────────────────────
if __name__ == "__main__":
    mp.set_start_method("spawn", force=True)   # macOS/Linux 안전

    t0   = time.time()
    aso  = ASOdesign(transid="NM_002415",
                     query_assembly=["mm39", "mm10"],
                     ref_assembly="hg38",
                     tile_length=17)

    aso.process_main(chunk_division=5, max_workers=4)
    print("Elapsed:", time.time() - t0, "sec")