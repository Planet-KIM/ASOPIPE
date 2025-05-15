
import os
import re
from functools import lru_cache

import RNA
from diskcache import Cache
from jklib.bioDB import CommonSNP

# ① 디스크 캐시: 10 GB 또는 항목 1 M개 선에서 LRU 자동 제거
CACHE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "./../../" ,".rnacofold_cache")
disk_cache = Cache(directory=str(CACHE_DIR), size_limit=10 * 1024 ** 3)
RNA.cvar.dangles = 2
RNA.cvar.noLonelyPairs = 1

def containCommonSNP(loc, cSNP=None):
    try:
        if cSNP == None:
            cSNP = CommonSNP()
        return len(cSNP.query(loc))>0
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