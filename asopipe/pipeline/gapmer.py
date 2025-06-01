import traceback


def filter_gapmer(result, species_prefix="gapmer_filtered",
                       keep_rule=None):
    """
    ▸ species_prefix 로 시작하는 키들만 '종'으로 간주하여 필터링한다.
      (== gapmer_filtered_* 컬럼끼리 비교)
    ▸ 기본 keep_rule = 한 종이라도 True → keep

    매개변수
    ----------
    gapmer_json : dict
        입력 딕셔너리
    species_prefix : str
        종 키를 식별하는 접두어 (기본: "gapmer_filtered")
    keep_rule : callable 또는 None
        각 행의 [bool, bool, ...] 리스트를 받아 True/False 를 반환.
        • None → any(vals)  (한 개라도 True면 keep)
        • 예) lambda vals: all(vals)   # 모든 종이 True여야 keep
    """
    # 1) 종(key) 목록: prefix 로 시작하는 키
    species_keys = [k for k in result if k.startswith(species_prefix)]
    if not species_keys:
        raise ValueError(f"'{species_prefix}' 로 시작하는 키가 없습니다.")

    # 2) 길이 일치 여부 확인
    n = len(result[species_keys[0]])
    if not all(isinstance(v, list) and len(v) == n for v in result.values()):
        raise ValueError("모든 리스트의 길이가 일치해야 합니다.")

    # 3) keep_rule 기본값: 한 개라도 True → keep
    if keep_rule is None:
        keep_rule = lambda vals: any(vals)

    # 4) keep_mask 계산
    keep_mask = [
        bool(keep_rule([result[k][i] for k in species_keys]))
        for i in range(n)
    ]

    # 5) 마스킹 적용 후 새 dict 반환
    return {
        k: [v for v, keep in zip(lst, keep_mask) if keep]
        for k, lst in result.items()
    }

def _getWingCoord(result, middle_size=10, gapmer_coord=''):
    try:
        coordL = [] 
        for sequence, length in zip(result["ASO_Sequence"], result["Length"]):
            if length < middle_size:
                raise ValueError(f"Sequence Length is {length}")
            wing_coord = None
            gap_seq = ''
            if gapmer_coord!='':
                wing_coord = tuple(map(int, gapmer_coord.split('_')))
                gap_seq = sequence[wing_coord[0]:-wing_coord[2]]
                coordL.append( [(wing_coord, gap_seq), None] )
            else:
                # sequence length가 2의 배수 
                if length % 2 == 0:
                    # if length == 20
                    # (20-10) //2 = 5 -> (5, 10, 5) seq[5:-5] result: (5,10,5),가운데 10개 seq 
                    wing_size = (length-middle_size)//2
                    wing_coord= (wing_size,middle_size,wing_size)
                    gap_seq   = sequence[wing_size: -wing_size]
                    coordL.append( [(wing_coord, gap_seq), None] )
                else:
                    s_size = (length-middle_size)//2
                    l_size = s_size + 1 

                    prime = []
                    #5<3
                    wing_coord = (s_size, middle_size, l_size)
                    gap_seq = sequence[s_size:-l_size]
                    prime.append( (wing_coord, gap_seq) )
                    
                    #5>3
                    wing_coord = (l_size, middle_size, s_size)
                    gap_seq = sequence[l_size:-s_size]
                    prime.append( (wing_coord, gap_seq) )
                    coordL.append( prime ) 
        return coordL
                    
    except Exception as e:
        print(traceback.format_exc())
        return e.args
        


def flaten_list(t):
    # [[a,b], [c,d]] -> [a,b,c,d]
    return [e for st in t for e in st]

def gapmer(result, middle_size=10, gapmer_coord='',target_assembly=["mm39"]):
    print(f"#gapmer filter start: target_assembly: {target_assembly}")
    wingcoordL = _getWingCoord(result=result, middle_size=middle_size, gapmer_coord=gapmer_coord)

    # test conservation
    for assembly in target_assembly:
        #print(f"assembly: {assembly}")
        gapmer_filtered = []
        for _coverage, _wobble, coords in zip(result[f"coverage_{assembly}"], result[f"wobble_{assembly}"], wingcoordL):
            # test gap CpG
            gapL= [True for _item in [item[1].upper() for item in coords if item!=None] if "CG" in _item ]
            gap_condition = True if True in gapL else False
            wingL= [item[0] for item in coords if item!=None]
            
            if gap_condition:
                gapmer_filtered.append(False)
                continue
            if str(_coverage) == "0"  or str(_coverage) == "None": # coverage값이 없을때만 가능 
                if not _wobble or not _wobble["GU_humanC"]:
                    gapmer_filtered.append(False)
                    continue   
                # GU_humanC에 값이 있다는 전제조건
                # no wobble in gap site
                #elif any([wing_coord[0]<=i<wing_coord[0]+wing_coord[1] for i in _wobble['GU_humanC']]): 
                if any(
                    wing_coord[0] <= i < wing_coord[0] + wing_coord[1]
                    for i in _wobble['GU_humanC']
                    for wing_coord in wingL
                    ):
                    gapmer_filtered.append(False)
                    continue
                # GU_humanC외 다른 것을 합쳐서 중복을 제거하고 길이를 GU_humanC의 길이만큼 빼고 측정했을 때 0보다 크다면 contiune 
                # 그렇다면 두개가 동일한 길이면 0이니까 통과할 수 있음 
                # no additional wobble site
                if len(set(flaten_list([v for k,v in _wobble.items() if k != 'GU_humanC'])) - set(_wobble['GU_humanC'])) > 0: 
                    gapmer_filtered.append(False)
                    continue    
                gapmer_filtered.append(True)
            else:
                # distnace가 있다는 것이기에 true로 처리
                gapmer_filtered.append(True)
                continue
        result[f"gapmer_filtered_{assembly}"] = gapmer_filtered
        result[f'gapmer_coords_{assembly}'] = [ ":".join(['_'.join(map(str, wing_coord)) for wing_coord in wingL])] * len(gapmer_filtered)
    result = filter_gapmer(result, species_prefix="gapmer_filtered", keep_rule=None)
    return result

     