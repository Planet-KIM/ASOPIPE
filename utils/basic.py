import collections
import gzip


def generate_long_dna_sequence(seq_length):
    """
    Generate a single long DNA sequence of given seq_length using GPU.
    Each element is a random integer 0~3, corresponding to A, C, T, G.
    """
    # Generate a random integer array of shape (seq_length,)
    random_indices = cp.random.randint(0, 4, size=(seq_length,), dtype=cp.int8)
    # Mapping for DNA bases in the order A, C, T, G
    mapping = np.array(['A', 'C', 'T', 'G'])
    # Bring the indices from GPU to CPU
    indices_cpu = random_indices.get()
    # Map each integer to its corresponding DNA base and join into one string
    long_sequence = ''.join(mapping[indices_cpu])
    return long_sequence


def tile_sequence(seq, tile_length):
    """
    Tile the input sequence using a sliding window of length tile_length.
    Returns a list of substrings; each substring is a tile.
    For example, if seq = "AATGGG" and tile_length = 3, then the tiles are:
    ["AAT", "ATG", "TGG", "GGG"].
    """
    num_tiles = len(seq) - tile_length + 1
    return [seq[i:i + tile_length] for i in range(num_tiles)]


def loadBlatOutput(blatOutputPath,by='transID',blacklist=['NR_106988']):
    """
    Parameters
    ----------
    blatOutputPath: string
    by: string
        user select the target name
        default value is 'transID'
    blacklist: array
        delete the blacklist
        default value is ['NR_106988']

    Results
    -------
    h: collections.defaultdict
        target are 'txnSta', 'txnEnd'
    """

    h = collections.defaultdict(list)

    if blatOutputPath.endswith('.gz'):
        f = gzip.open(blatOutputPath)
    else:
        f = open(blatOutputPath)

    for line in f:

        if line[0] == '#':
            continue

        r = _processBlatLine(line)

        if r['transID'] in blacklist:
            continue

        h[r[by]].append(r)

    from operator import itemgetter, attrgetter

    for k,vL in list(h.items()):
        h[k] = sorted(vL,key=itemgetter('txnSta','txnEnd'))

    return h



def _processBlatLine(line):
    """
    user not use this method
    coordinates: 0-base [,)
    you need to convert coordinate system to use this in jklib functions

    Parameters
    ----------
    line: string

    Results
    -------
    h : dictionary

    """
    tokL = line.rstrip().split('\t')

    h = {}

    h['transName'] = tokL[0]
    h['transID'] = tokL[1]
    h['chrom'] = tokL[2]
    h['chrNum'] = tokL[2][3:]
    h['strand'] = tokL[3]
    h['txnSta'] = int(tokL[4])
    h['txnEnd'] = int(tokL[5])
    h['txnLen'] = h['txnEnd'] - h['txnSta']
    h['cdsSta'] = int(tokL[6])
    h['cdsEnd'] = int(tokL[7])
    h['exnList'] = list(map(lambda x,y: (int(x),int(y)), tokL[9].split(',')[:-1], tokL[10].split(',')[:-1]))
    h['exnLenList'] = [e-s for (s,e) in h['exnList']]
    h['exnLen'] = sum(h['exnLenList'])

    if len(tokL) > 12:
        h['geneName'] = tokL[12]

    h['cdsList'] = []
    frontL, backL = [],[]

    if h['cdsSta'] != h['cdsEnd']:

        for (s,e) in h['exnList']:

            frontL.append((min(s,h['cdsSta']),min(e,h['cdsSta'])))
            h['cdsList'].append((max(s,h['cdsSta']),min(e,h['cdsEnd'])))
            backL.append((max(s, h['cdsEnd']), max(e, h['cdsEnd'])))

        frontL = [x for x in frontL if x[0] < x[1]]
        h['cdsList'] = [x for x in h['cdsList'] if x[0] < x[1]]
        backL = [x for x in backL if x[0] < x[1]]

    if h['strand'] == '+':
        h['utr5'] = frontL
        h['utr3'] = backL
    elif h['strand'] == '-':
        h['utr5'] = backL
        h['utr3'] = frontL
    else:
        raise Exception

    h['utr5Len'] = sum([e-s for (s,e) in h['utr5']])
    h['utr3Len'] = sum([e-s for (s,e) in h['utr3']])
    h['cdsLen'] = sum([e-s for (s,e) in h['cdsList']])

    h['intron'] = []

    for i in range(len(h['exnList'])-1):
        h['intron'].append((h['exnList'][i][1],h['exnList'][i+1][0]))

    return h

"""
# 설정: 원본 시퀀스 길이 300,000, 타일 길이 17
LONG_SEQ_LENGTH = 300_000
TILE_LENGTH = 17

# GPU를 이용해 길이 300,000의 DNA 시퀀스를 생성합니다.
long_dna_sequence = generate_long_dna_sequence(LONG_SEQ_LENGTH)
print("Generated long DNA sequence (first 50 bases):")
print(long_dna_sequence[:50], len(long_dna_sequence))
# 생성된 시퀀스를 타일링: 슬라이딩 윈도우 방식으로 타일 추출 (stride=1)
dna_tiles = tile_sequence(long_dna_sequence, TILE_LENGTH)

# 예시: 첫 5개의 타일을 출력합니다.
print("\nFirst 5 tiles of length {}:".format(TILE_LENGTH))
for tile in dna_tiles[:5]:
    print(tile)

print("\nTotal number of tiles generated: {}".format(len(dna_tiles)))
"""

