import sys
sys.path.append("../")  # 상위 디렉토리 추가
from asopipe.main import ASOdesign
import multiprocessing as mp
import time
# ────────────────────────────
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