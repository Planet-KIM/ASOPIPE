
def levenshtein_distance(s1, s2, cost_sub=1, cost_ins=1, cost_del=1):
    """
    두 서열 s1과 s2 사이의 Levenshtein distance(편집 거리)를 동적 계획법으로 계산합니다.
    이 함수는 substitution, insertion, deletion 각 비용을 반영합니다.
    
    Parameters:
      s1, s2 (str): 비교할 두 서열 (예, aligned sequences)
      cost_sub (int): 대체(cost for substitution)
      cost_ins (int): 삽입(cost for insertion)
      cost_del (int): 삭제(cost for deletion)
      
    Returns:
      int: s1을 s2로 변환하는 데 필요한 최소 편집 비용
    """
    m, n = len(s1), len(s2)
    # dp[i][j]: s1[0:i]를 s2[0:j]로 변환할 때 최소 비용
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    # 초기값 설정: 한 쪽이 빈 문자열일 때
    for i in range(m + 1):
        dp[i][0] = i * cost_del
    for j in range(n + 1):
        dp[0][j] = j * cost_ins
        
    # 동적 계획법을 통해 계산
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s1[i - 1] == s2[j - 1]:
                cost = 0
            else:
                cost = cost_sub  # 대체 비용
            dp[i][j] = min(
                dp[i - 1][j] + cost_del,       # 삭제 연산
                dp[i][j - 1] + cost_ins,       # 삽입 연산
                dp[i - 1][j - 1] + cost        # 대체 연산 (또는 일치하는 경우 0)
            )
    return dp[m][n]


def average_edit_distance(msa, cost_sub=1, cost_ins=1, cost_del=1):
    """
    다중 서열 정렬(MSA)에 포함된 모든 서열 쌍에 대해 편집 거리를 계산하고 그 평균을 반환합니다.
    
    Parameters:
      msa (list of str): 동일 길이(정렬된)의 서열 리스트
      cost_sub, cost_ins, cost_del: edit distance 계산 시 각 연산의 비용
      
    Returns:
      float: 모든 서열 쌍에 대한 평균 편집 거리
    """
    n = len(msa)
    if n < 2:
        return 0  # 한 개 이하의 서열이면 비교 불가
    
    total_distance = 0
    count = 0
    # 모든 서열 쌍을 순회하며 편집 거리를 계산
    for i in range(n):
        for j in range(i + 1, n):
            d = levenshtein_distance(msa[i], msa[j], cost_sub, cost_ins, cost_del)
            total_distance += d
            count += 1
    avg_distance = total_distance / count
    return int(avg_distance)