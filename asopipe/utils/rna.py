
import os
import re
import traceback
from functools import lru_cache

import RNA
from diskcache import Cache
#from jklib.bioDB import CommonSNP

from cyvcf2 import VCF  

# ① 디스크 캐시: 10 GB 또는 항목 1 M개 선에서 LRU 자동 제거
CACHE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "./../../" ,".rnacofold_cache")
disk_cache = Cache(directory=str(CACHE_DIR), size_limit=10 * 1024 ** 3)
RNA.cvar.dangles = 2
RNA.cvar.noLonelyPairs = 1

def loadSNP(locStr, dbsnp_path=None):
    try:
        check_type_dbsnp = str(type(dbsnp_path)).lower()
        cSNP = dbsnp_path
        if "none" in check_type_dbsnp:
            dbsnp_path = "/Users/dowonkim/Dropbox/data/VCF/dbsnp.bcf"
            dbsnp_index_path = "/Users/dowonkim/Dropbox/data/VCF/dbsnp.bcf.csi"  
            cSNP = VCF(dbsnp_path)
            cSNP.set_index(index_path=dbsnp_index_path)
            #cSNP = VCF(dbsnp_path)          # .csi 인덱스 자동 사용
        elif 'vcf' in check_type_dbsnp:
            pass
        else:
            raise ValueError(f"Check your dpsnp_path argument.(Now: {dbsnp_path})")
        # “This format (chrom:chrSta-chrEnd) used to be accepted by dbSNP.
        locStr = locStr.rstrip('-').rstrip('+').replace("chr", "")
        resultL = [(v.CHROM, v.POS ,v.REF, v.ALT[0], [info for info in v.INFO])  for v in cSNP(locStr)] 
        return resultL
    except Exception as e:
        print(traceback.format_exc())
        return e.args # [] error default value is empty list. 
    

def containCommonSNP(loc, cSNP=None):
    try:
        locStr = loc.toString()
        resultL = loadSNP(locStr, cSNP)
        return resultL
        #return len(resultL)>0
    except:
        print(loc.toString())
        raise('Error')

def containGquad(sequence):
    return 'GGGG' in sequence or 'CCCC' in sequence

def containGquad2(sequence, loop=7):
    reg_seq = f"(G{3,})[ATCG]{1,loop}(G{3,})[ATCG]{1,loop}(G{3,})[ATCG]{1,loop}(G{3,})"
    G4_RE = re.compile(reg_seq)
    return bool(G4_RE.search(sequence))

def countCpG(sequence):
    return sequence.count('CG')

def GCcontent(sequence):
    comp = dict([(x, sequence.count(x)/len(sequence)) for x in 'ATCG'])
    return comp['G'] + comp['C']

def RNAcofold(sequence):
    seq = sequence + '&' + sequence
    tokL = RNA.co_pf_fold(seq)
    return tokL[4], tokL[2]

@lru_cache(maxsize=1024)
def _cofold_in_memory(dimer):
    """
    가장 먼저 메모리 LRU를 조회하고,
    miss이면 diskcache → 계산 순으로 진행.
    """
    # ── 1) 디스크 캐시 조회
    hit = disk_cache.get(dimer, default=None)
    if hit is not None:
        return hit  # (FB, FcAB) tuple

    # ── 2) 실제 계산
    tokL = RNA.co_pf_fold(dimer)
    result = (tokL[4], tokL[2])      # FB, FcAB

    # ── 3) 디스크 캐시 저장 후 반환
    disk_cache.set(dimer, result)
    return result

def RNAcofold2(sequence):
    """
    Parameters
    ----------
    sequence : str  ― 단량체 서열 (A/C/G/T/U)

    Returns
    -------
    FB, FcAB : float, float
        기존 구현과 완전 동일한 두 자유에너지 값
    """
    seq = sequence.upper().replace('T', 'U')
    dimer = f"{seq}&{seq}"
    homo, mono = _cofold_in_memory(dimer)
    return homo, mono