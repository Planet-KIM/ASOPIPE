"""
MAF 파일의 인덱스를 생성하고, 특정 유전체 좌표에 해당하는 다중 정렬 정보를 조회하는 모듈입니다.
"""
import os
import traceback
import collections

from io import TextIOWrapper
from bx import interval_index_file
import bx.align.maf


def compl(seq):
    return seq.translate(seq.maketrans('ACGTacgt','TGCAtgca'))

# 현재는 original에서는 지원하지 않는데 해당이유는 모르겠음 어떤 값을 넣어도 wobble 2이 이하로만 실행 하지만 original에서는 아에 0취급
def check_wobble(r, q, anti_strand="-", wob=2):
    """
    Check whether two DNA/RNA sequences can pair with at most `wob`
    tolerated wobble mismatches.

    Parameters
    ----------
    r : str   # reference (human) sequence, 5'→3'
    q : str   # query    (mouse) sequence, 5'→3'
    wob : int # max tolerated wobble sites (default 2)

    Returns
    -------
    (bool, dict or None)
        True  + dict : 통과.  딕트 구조는 원본 그대로
        False + None : 탈락
    """
    try:
        if r == None or q == None:
            raise ValueError
        r, q = r.upper().rstrip('-'), q.upper().rstrip('-')
        if anti_strand == '+':
            # complement만 하는 이유를 아직 파악하지 못함 
            r,q = compl(r), compl(q)
        # 1) 갭·길이 불일치 즉시 탈락
        if '-' in r or '-' in q or len(r) != len(q):
            raise ValueError
            #return False
            #return False, None

        # 2) 불일치 위치 파악
        mismatch_idx = [i for i, (rb, qb) in enumerate(zip(r, q)) if rb != qb]
        if len(mismatch_idx) > wob:
            raise ValueError
            #return False
            #return False, None

        # 3) 네 가지 wobble 클래스(중복 허용)
        gu_humanC, gu_otherC, i_humanC, i_otherwise = [], [], [], []
        for idx in mismatch_idx:
            rb, qb = r[idx], q[idx]

            # GU wobble
            if (rb == 'C' and qb == 'T') or (rb == 'A' and qb == 'G'):
                gu_humanC.append(idx)
            if (rb == 'T' and qb == 'C') or (rb == 'G' and qb == 'A'):
                gu_otherC.append(idx)

            # I wobble
            if rb == 'C' and qb != 'G':
                i_humanC.append(idx)
            if rb not in 'GC' and qb != 'G':
                i_otherwise.append(idx)

        # 4) 허용되지 않은 불일치가 섞여 있나?
        # 위조건이 전부매치하여도 위의 네 가지 속하지 않을 수 있음
        if set(gu_humanC + gu_otherC + i_humanC + i_otherwise) != set(mismatch_idx):
            raise ValueError
            #return False
            #return False, None

        # 5) 원본과 완전히 같은 딕트 포맷 반환
        wobble_dict = {
            'GU_humanC':   gu_humanC,
            'GU_otherC':   gu_otherC,
            'I_humanC':    i_humanC,
            'I_otherwise': i_otherwise
        }
        if anti_strand == '-':
            # reverse complement
            for _key, _value in wobble_dict.items():
                if not _value:
                    continue
                # if len(r) == len(q) : True
                wobble_dict[_key] = sorted([len(r) - 1 - i for i in _value])
        return wobble_dict
        #return True, wobble_dict
    except ValueError:
        return False


def build_index(maf_file, species=None):
    """
    주어진 MAF 파일에 대해 인덱스 파일(maf_file.index)을 생성합니다.
    species가 지정되어 있으면 해당 종에 해당하는 컴포넌트만 인덱싱합니다.
    """
    index_file = maf_file + ".index"
    try:
        with open(maf_file, "rb") as maf_in:
            # 바이너리 파일을 ASCII 텍스트 스트림으로 변환
            with TextIOWrapper(maf_in, encoding="ascii") as maf_in_text:
                maf_reader = bx.align.maf.Reader(maf_in_text)
                indexes = interval_index_file.Indexes()
                # 파일 포인터의 위치를 기반으로 인덱스 생성
                for block in iter(lambda: next(maf_reader), None):
                    pos = maf_reader.file.tell()
                    for component in block.components:
                        if species is not None and component.src.split('.')[0] not in species:
                            continue
                        indexes.add(component.src,
                                      component.forward_strand_start,
                                      component.forward_strand_end,
                                      pos,
                                      max=component.src_size)
                # 생성된 인덱스를 파일에 기록
                with open(index_file, 'wb') as out:
                    indexes.write(out)
    except Exception as e:
        print(f"Error building index for {maf_file}: {e}")
        

class MultipleAlignmentReader:
    def __init__(self, **kwargs):
        """
        MultipleAlignmentReader를 초기화합니다.
        
        인수:
          - ref_assembly: 기준 유전체 (기본: 'hg38')
          - query_assembly: 비교 대상 유전체 (예: 'macFas5' 또는 다른 값)
          - maf_dir: MAF 파일이 저장된 디렉토리. 지정하지 않으면 REFERENCE에서 동적으로 가져옴.
        """
        options = {'ref_assembly': 'hg38', 'query_assembly': None, 'maf_dir': None}
        options.update(kwargs)
        
        self.ref_assembly = options['ref_assembly']
        self.query_assembly = options['query_assembly']
        self.maf_dir = options['maf_dir']
        
        if not self.query_assembly:
            raise ValueError("query_assembly must be specified.")
        
        if self.maf_dir is None:
            # eval 대신 getattr 사용: REFERENCE.Genome.UCSC.MultipleAlignment.<ref_assembly>
            try:
                self.maf_dir = getattr(REFERENCE.Genome.UCSC.MultipleAlignment, self.ref_assembly)
            except AttributeError:
                raise ValueError(f"Reference directory for {self.ref_assembly} not found in REFERENCE.")
        
        self.maf_file = os.path.join(self.maf_dir, f'{self.ref_assembly}.{self.query_assembly}.synNet.maf')
        self.index = collections.defaultdict(dict)
        
        if self.query_assembly.lower() != 'macfas5':
            self.load_mafnidx(self.maf_file)
        else:
            # macFas5의 경우 염색체별로 인덱스를 개별적으로 로드
            self.index['macFas5'] = collections.defaultdict(dict)
    
    def load_mafnidx(self, maf_file=None):
        """
        주어진 MAF 파일의 인덱스를 로드합니다.
        인덱스 파일이 없으면 build_index 함수를 호출하여 생성한 후 로드합니다.
        """
        if not maf_file:
            maf_file = self.maf_file
        if not os.path.exists(maf_file + '.index'):
            build_index(maf_file, species=[self.ref_assembly, self.query_assembly])
        self.index[self.query_assembly] = bx.align.maf.Indexed(maf_file, maf_file + '.index')
    
    def _get_index(self, chrom):
        """
        주어진 염색체(chrom)에 해당하는 인덱스 객체를 반환합니다.
        query_assembly가 'macfas5'인 경우, 염색체별로 인덱스를 동적으로 로드합니다.
        """
        if self.query_assembly.lower() == 'macfas5':
            if not self.index[self.query_assembly].get(chrom):
                maf_file = os.path.join(macfas5_maf_dir, f'{chrom}.maf')
                if not os.path.exists(maf_file + '.index'):
                    print('No index file found for', maf_file + '.index')
                    build_index(maf_file, species=[self.ref_assembly, self.query_assembly])
                self.index[self.query_assembly][chrom] = bx.align.maf.Indexed(maf_file, maf_file + '.index')
            return self.index[self.query_assembly][chrom]
        else:
            return self.index[self.query_assembly]
    
    def query(self, loc, verbose=True):
        """
        지정한 유전체 영역(loc)에 대해 다중 정렬 정보를 조회합니다.
        loc 객체는 아래와 같은 속성을 가져야 합니다.
          - chrom: 염색체 이름
          - chrSta: 시작 좌표 (1-base, UCSC 좌표 체계)
          - chrEnd: 종료 좌표
          - strand: (선택사항) 염색체의 방향 정보
        결과:
          - 해당 영역의 정렬된 서열 정보를 포함하는 딕셔너리 (없는 경우 None)
        """
        idx = self._get_index(loc.chrom)
        region_name = f"{self.ref_assembly}.{loc.chrom}"
        s, e = loc.chrSta - 1, loc.chrEnd  # 0-base 좌표로 보정
        
        for alignment in idx.get(region_name, s, e):
            try:
                #print(region_name, s, e)
                region_alignment = alignment.slice_by_component(component_index=region_name, start=s, end=e)
            except ValueError:
                #print(traceback.format_exc())
                continue
            
            seqs_by_org = {}
            for component in region_alignment.components:
                if component.src.startswith((self.ref_assembly, self.query_assembly)):
                    seqs_by_org[component.src] = component.text
            if not verbose:
                seqs_by_org = {k.split('.')[0]: v for k, v in seqs_by_org.items()}
            
            if len(seqs_by_org) < 2:
                # 충분한 서열 정보가 없는 경우, 한 베이스씩 조회하는 fallback
                return self.query_one_by_one(loc, verbose=verbose)
            return seqs_by_org
        
        # 조회 결과가 없는 경우 fallback을 시도
        return self.query_one_by_one(loc, verbose=verbose)
    
    def query_one_by_one(self, loc, verbose=True):
        """
        fallback 메서드로, 영역 내의 각 베이스를 개별적으로 조회하여 정렬된 서열을 얻습니다.
        결과적으로 각 위치의 결과를 병합하여 전체 서열을 구성합니다.
        """
        idx = self._get_index(loc.chrom)
        region_name = f"{self.ref_assembly}.{loc.chrom}"
        sta, end = loc.chrSta - 1, loc.chrEnd
        
        results = []
        for pos in range(sta, end):
            for alignment in idx.get(region_name, pos, pos + 1):
                try:
                    region_alignment = alignment.slice_by_component(region_name, pos, pos + 1)
                except ValueError:
                    continue
                seqs_by_org = {}
                for component in region_alignment.components:
                    if component.src.startswith((self.ref_assembly, self.query_assembly)):
                        seqs_by_org[component.src] = component.text
                if not verbose:
                    seqs_by_org = {k.split('.')[0]: v for k, v in seqs_by_org.items()}
                
                if len(seqs_by_org) < 2:
                    results.append(None)
                    break
                results.append(seqs_by_org)
                break
        #print('results:', results)
        if not results:
            return None
        if len(results) == 1:
            return results[0]
        
        merged = results[0]
        if not merged:
            return None
        
        for d in results[1:]:
            if not d:
                return None
            for key in merged.keys():
                try:
                    merged[key] += d[key]
                except Exception:
                    return None
        return merged
