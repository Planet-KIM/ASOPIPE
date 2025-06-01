import os
import sys
sys.path.append("../")  # 상위 디렉토리 추가
from asopipe.main import run_ASOdesign

# ────────────────────────────
# ──────────────────────────────────────────────────────────────
# 3) 실행 예시  ────────────────────────────────────────────────
if __name__ == "__main__":
    transid ="NM_133379"#TTN
    #transid="NM_002415" #MIF
    if not os.path.isdir(f"./{transid}"):
        os.makedirs(f"{transid}")
    result  = run_ASOdesign(transid=transid,
                     query_assembly=["mm39","mm10", "macFas5", "rn7", "rn6", "rn5"],
                     refFlat_path = "/home/commons/Reference/UCSC/hg38/refseq/refFlat_200817.txt",
                     maf_dir='/home/commons/Reference/UCSC/hg38/maf',
                     dbsnp_path ="/home/jkportal/data/VCF/dbsnp.bcf",
                     dbsnp_path ="/home/jkportal/data/VCF/dbsnp.bcf.csi",
                     ref_assembly="hg38",
                     k_min=19, k_max=20,
                     chunk_division=5, max_workers=1, wobble=0, to_df=False, gapmer_filtered=True,
                     to_csv=True, output_path=f"./{transid}")