#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import multiprocessing as mp

import concurrent.futures
import numpy as np
import pandas as pd
from collections import defaultdict

import os
import time
import traceback
#import editdistance
#from functools import partial
from cyvcf2 import VCF   # pip install cyvcf2

from asopipe.utils.basic import loadBlatOutput
from asopipe.utils.coverage import average_edit_distance
from asopipe.utils.csv import save_csv_std, save_csv_pyarrow, save_csv_polars
from asopipe.utils.align.maf_th import check_wobble
from asopipe.utils.align.maf_th import MultipleAlignmentReader
from asopipe.utils.rna import RNAcofold2, containCommonSNP, containGquad2, countCpG, GCcontent
from asopipe.pipeline.gapmer import gapmer
from jklib.genome import locus, getRegionType
#from jklib.bioDB import CommonSNP


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
                 dbsnp_path ="/Users/dowonkim/Dropbox/data/VCF/dbsnp.bcf",
                 dbsnp_index_path ="/Users/dowonkim/Dropbox/data/VCF/dbsnp.bcf.csi",  # .csi 인덱스 자동 사용
                 query_assembly=["mm39"],      # tuple/리스트 허용
                 ref_assembly="hg38",
                 tile_length=17):
        print(f"[ASOdesign] transid={transid}")
        self.refFlat      = loadBlatOutput(refFlat_path, by='transID')
        #self.cSNP = CommonSNP()
        
        self.cSNP = VCF(dbsnp_path)
        self.cSNP.set_index(index_path=dbsnp_index_path)  # .csi 인덱스 자동 사용
        self.maf_dir      = maf_dir
        self.query_asm    = query_assembly
        self.ref_asm      = ref_assembly
        self.transid      = transid
        self.tile_length  = tile_length
        
        #test
        self.t0   = time.time()
        self.endtime = 0
        #
        
        self.transInfo    = self._get_transInfo()
        if transid != None:
            self.transName = self.transInfo['transName']
            self.strand = self.transInfo['strand']
            self.anti = '+' if self.transInfo['strand'] == '-' else '-'
            self.chrom = self.transInfo['chrom']
            self.txnSta = self.transInfo['txnSta']
            self.txnEnd   = self.transInfo['txnEnd']
        else:
            raise ValueError("transid is None")
        self.txn_tiles    = self._tile_region()
        self.txn_tile_seq = self._tile_txn_seq()


    # ───── 내부 헬퍼 ──────────────────────────────────────────
    def _get_transInfo(self):
        return [t for t in self.refFlat[self.transid]
                if len(t['chrom'].split('_')) == 1][0]

    def _tile_region(self):
        L     = self.tile_length
        return [locus(f"{self.chrom}:{i}-{i+L-1}{self.anti}")   # inclusive end = i+L-1
            for i in range(self.txnSta, self.txnEnd-L+1)]

    def _tile_txn_seq(self):
        L     = self.tile_length
        txn_seq = locus(f"{self.chrom}:{self.txnSta}-{self.txnEnd-1}{self.anti}").twoBitFrag().upper()
        result_array = [ txn_seq[i:i+L] for i in range(0, self.txnEnd-self.txnSta-L+1)]
        if self.anti == '-':
            result_array = result_array[::-1]
        return result_array

    def _chunks(self, seq, n):
        """seq 를 n 등분하여 순차적으로 yield"""
        k = max(1, len(seq)//n)
        for i in range(0, len(seq), k):
            yield seq[i:i+k]

    # ───── 퍼블릭 메서드 ─────────────────────────────────────
    def process_main(self, chunk_division=3, max_workers=3, wobble=2, to_df=True, gapmer_filtered=False, to_csv=False, output_path=None):
        print(f"#tiles={len(self.txn_tiles)}, assemblies={self.query_asm}, tile_length={self.tile_length}, wobble={wobble}")     
        #
        _all_results_locInfo = []
        for chunk_loc, chunk_seq in zip(self._chunks(self.txn_tiles, chunk_division), self._chunks(self.txn_tile_seq, chunk_division)):
            chunk_locInfo = [self.getlocInfo(tile_loc, tile_seq) for tile_loc, tile_seq in zip(chunk_loc, chunk_seq)] 
            _all_results_locInfo.extend(chunk_locInfo)
            #print("Elapsed:", self.endtime, "sec")
        all_results_locInfo = self._flatten_dict(_all_results_locInfo)

        #
        result_dict = {}
        for asm in self.query_asm:                      # 어셈블리별로 별도 풀
            print(f"Assembly: [{asm}] prcessing...")
            with concurrent.futures.ProcessPoolExecutor(
                    max_workers=max_workers,
                    initializer=_init_worker,
                    initargs=(self.maf_dir, self.ref_asm, asm),
            ) as ex:
                all_results_maf, all_editdist = [], []
                for chunk_loc, chunk_seq in zip(self._chunks(self.txn_tiles, chunk_division), self._chunks(self.txn_tile_seq, chunk_division)):
                    # region 문자열만 워커에 전달
                    #chunk_locInfo = list(ex.map(self.getlocInfo, chunk))
                    #chunk_locInfo = [self.getlocInfo(tile_loc, tile_seq) for tile_loc, tile_seq in zip(chunk_loc, chunk_seq)]
                    chunk_maf_results = list(ex.map(_query_region, chunk_loc))
                    dists   = [self._editdistance_safe(r) for r in chunk_maf_results]
                    #all_results_locInfo.extend(chunk_locInfo)
                    all_results_maf.extend(chunk_maf_results)
                    all_editdist.extend(dists)
                result_dict[asm] = {"maf_seq": all_results_maf, "coverage": all_editdist}
                #result_dict[asm] = {"maf_seq": all_results_maf, "coverage": all_editdist, "locInfo": all_results_locInfo}
                #print(asm, result_dict[asm]["maf_seq"][:3])
        #print("finished. first 3 results:", result_dict[self.query_asm[0]]["locInfo"][:1])
        #print("finished. first 3 results:", result_dict[self.query_asm[0]]["maf_seq"][:1])
        #print("finished. first 3 results:", result_dict[self.query_asm[0]]["coverage"][:1])
        sort_result_dict = self._sort_dict(result_dict)
        #(zip(result["mm39"]["hg38"], result["mm39"]["mm39"],  
        #print(sort_result_dict)
        for _assembly in sort_result_dict.keys():
            wobble_pair = []
            #for human, other in zip(sort_result_dict[_assembly]["hg38"], sort_result_dict[_assembly][_assembly]):
            #human, other = sort_result_dict[_assembly][f"maf_seq_{_assembly}"].split(":")
            for _maf_seq in sort_result_dict[_assembly][f"maf_seq_{_assembly}"]:
                _human, _other = _maf_seq.split(":")
                wobble_pair.append(check_wobble(r=_human, q=_other, wob=wobble, anti_strand=self.anti))
            #print(wobble_pair, "wobble")
            sort_result_dict[_assembly].update({f"wobble_{_assembly}": wobble_pair})
            all_results_locInfo.update(sort_result_dict[_assembly])
            #all_results_locInfo.update(sort_result_dict[_assembly].update({f"wobble_{_assembly}": wobble_pair}))
            #sort_result_dict[_assembly]["wobble"] = wobble_pair
            #sort_result_dict[_assembly].update(all_results_locInfo)
        #print(all_results_locInfo)
        all_results_locInfo_copy = all_results_locInfo.copy()
        result = {"default": all_results_locInfo_copy} 
        if output_path == None and to_csv:
            output_path = os.path.dirname(os.path.realpath(__file__))
        if to_csv:
            save_csv_pyarrow(data_dict=all_results_locInfo_copy,
                             path=os.path.join(output_path,f"{self.transName}_{self.transid}_{self.tile_length}mer_wobble_{wobble}.csv"),
                             toString=True)
        if gapmer_filtered:
            all_results_locInfo= gapmer(result=all_results_locInfo, middle_size=10, gapmer_coord='', target_assembly=self.query_asm)
            if to_csv:
                save_csv_pyarrow(data_dict=all_results_locInfo,
                                 path=os.path.join(output_path,f"{self.transName}_{self.transid}_{self.tile_length}mer_wobble_{wobble}_gapmer_filtered.csv"),
                                 toString=True)
            result["gapmer"] = all_results_locInfo
        if to_df:
            #result_df  = self.apply_df(sort_result_dict)
            #result_df = pd.DataFrame(all_results_locInfo)
            result["default"] = pd.DataFrame(all_results_locInfo_copy)
            if gapmer_filtered:
                result["gapmer"] = pd.DataFrame(all_results_locInfo)
            return result
        else:
            return result

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
            #flagL = [ e[3] for e in loc_tmp.regionType()]
            h_tmp = {self.transInfo['chrom']: [self.transInfo]}
            flagL = [ e[3] for e in getRegionType(h_tmp, loc_tmp)]
            flag = []
            
            for a in flagL:
                #print(a[0][0])
                if (('cds' in a[0][0]) or ('utr' in a[0][0])): #exon flag ['cds','utr']
                    flag.append('e') 
                elif 'int' in a[0][0] : #intron
                    flag.append('i')
            
            #regionT = '/'.join(map(str, loc.regionType())) #test
            regionT = '/'.join(map(str, getRegionType(h_tmp, loc))) #test
            homo, mono = RNAcofold2(sequence)
            
            flag = ':'.join(flag)
            #self.t0 = time.time()
            snp_data = containCommonSNP(loc, self.cSNP)
            #self.endtime = self.endtime+ (time.time() - self.t0)
            #'type','gene','transcriptID','locus','sequence','length','regionType','commonSNP','Gquad','CpG','GC_content','Homo_dimer','Monomer'
            a = {"Type": flag,
                    "Gene": self.transInfo['transName'],
                    "TranscriptID": self.transInfo['transID'],
                    "ASO_Locus": locStr, "ASO_Sequence": sequence, "Length": len(sequence), "RegionType": regionT,#test
                    "CommonSNP": snp_data,
                    "Gquad": containGquad2(sequence),
                    "CpG": countCpG(sequence),
                    "GC_Content": GCcontent(sequence),
                    "Homo_Dimer": homo, "Monomer": mono
                    }
            
            return a
        except Exception as e:
            print(traceback.format_exc())
            return e.args
        
    
    def _flatten_dict(self, data):
        """
        Flatten a list of dictionaries into a single dictionary.
        #support to this typo [{type: i:e}, {type: i:i}] -> {type: [i:e, i:i]}
        """
        try:
            if isinstance(data, dict):
                raise Exception
            #result = defaultdict(list); [result[k].append(v) for d in data for k, v in d.items()]
            # ❶ 결과를 담을 defaultdict 생성
            result = defaultdict(list)
            # 1) 모든 레코드에서 등장한 키 수집
            keys, seen = [], set()
            for rec in data:
                if rec is None:
                    continue
                for k in rec:
                    if k not in seen:
                        seen.add(k)
                        keys.append(k)

            # 2) 각 위치에 맞춰 값 채워 넣기
            result = defaultdict(list)
            for rec in data:
                for k in keys:
                    value = None if (rec is None) else rec.get(k, None)
                    result[k].append(value)
            return dict(result)
        except Exception as e:
            #print(traceback.format_exc())
            #print(e.args)
            return data
    
    def _sort_dict(self, result):
        remake_output_result = {}
        for _assembbly in result.keys():
            remake_output_assembly = {}
            for _key in result[_assembbly].keys():
                flattened_result = self._flatten_dict(result[_assembbly][_key]) 
                if _key == "coverage":
                    if not flattened_result or isinstance(flattened_result, list):
                        flattened_result = {f"coverage_{_assembbly}": result[_assembbly][_key]}
                    else:
                        flattened_result = {f"coverage_{_assembbly}": flattened_result} # using editdistance
                elif _key == "maf_seq":
                    #if isinstance(flattened_result, dict):
                    if not any(x is not None for x in flattened_result):
                        flattened_result = {"hg38": result[_assembbly][_key], _assembbly: result[_assembbly][_key]}
                    flattened_result={f"{_key}_{_assembbly}": (flattened_result["hg38"], flattened_result[_assembbly])}
                    _list1, _list2 = flattened_result[f"{_key}_{_assembbly}"]  # 튜플을 풀어서 리스트 2개로 나눔
                    assert len(_list1) == len(_list2), "두 리스트의 길이가 다릅니다."
                    flattened_result[f"{_key}_{_assembbly}"] = [f"{a}:{b}" for a, b in zip(_list1, _list2)]
                #print(_assembbly,flattened_result)
                remake_output_assembly.update(flattened_result)
            remake_output_result[_assembbly] = remake_output_assembly
        return remake_output_result
    
    def apply_df(self, result_dict):
        """
        Convert the result dictionary to a DataFrame.
        """
        try:
            result_df = {}
            for _assembly in result_dict.keys():
                df = pd.DataFrame(result_dict[_assembly])
                df["assembly"] = _assembly
                result_df[_assembly] = df
            return result_df
        except Exception as e:
            print(e.args)
            return e.args

def run_ASOdesign(transid="NM_002415",
                  refFlat_path="/Users/dowonkim/Dropbox/data/UCSC/hg38/refFlat/refFlat_200817.txt",
                  maf_dir='/Users/dowonkim/Dropbox/data/offtarget_test/maf',
                  dbsnp_path ="/Users/dowonkim/Dropbox/data/VCF/dbsnp.bcf",
                  dbsnp_index_path ="/Users/dowonkim/Dropbox/data/VCF/dbsnp.bcf.csi",  # .csi 인덱스 자동 사용
                  query_assembly=["mm39"],      # tuple/리스트 허용
                  ref_assembly="hg38",
                  k_min=17, k_max=17,
                  chunk_division=5, max_workers=1, wobble=0, to_df=False, gapmer_filtered=True, to_csv=True,
                  output_path=None):
    """
    Run the ASOdesign process with default parameters.
    """
    try:
        if k_min > k_max:
            raise ValueError("k_min should be less than or equal to k_max")
        if k_min < 17 or k_max < 17:
            raise ValueError("k_min and k_max should be greater than or equal to 17")
        if k_max > 25:
            raise ValueError("k_max should be less than or equal to 25")
        t0   = time.time()
        mp.set_start_method("spawn", force=True)   # macOS/Linux 안전
        result_list = []
        for tile_length in range(k_min, k_max+1):
            aso = ASOdesign(transid=transid,
                            refFlat_path=refFlat_path,
                            maf_dir=maf_dir,
                            dbsnp_path =dbsnp_path,
                            dbsnp_index_path=dbsnp_index_path,  # .csi 인덱스 자동 사용
                            query_assembly=query_assembly,      # tuple/리스트 허용
                            ref_assembly=ref_assembly,
                            tile_length=tile_length)
            
            result = aso.process_main(chunk_division=chunk_division, max_workers=max_workers, wobble=wobble,
                                        to_df=to_df,
                                        gapmer_filtered=gapmer_filtered, 
                                        to_csv=to_csv,
                                        output_path=output_path)
            result_list.append({"tile_length": tile_length, "result": result})
        print("Elapsed:", time.time() - t0, "sec")
        return result_list
    except Exception as e:
        print(traceback.format_exc())
        print(e.args)
        return e.args