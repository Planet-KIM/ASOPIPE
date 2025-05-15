import sys
sys.path.append("..")  # 상위 디렉토리 추가
from mypkg.asopipe import ASOdesign, run_aso_design
import multiprocessing as mp
import time
# ──────────────────────────────────────────────────────────────
# 4) 실행 예시  ────────────────────────────────────────────────
if __name__ == "__main__":
    mp.set_start_method("spawn", force=True)   # macOS/Lsinux 안전

    t0   = time.time()
    aso  = ASOdesign(transid="NM_002415",
                     query_assembly=["mm39", "mm10"],
                     ref_assembly="hg38",
                     tile_length=17)

    # Celery worker에서 실행할 때
    tiles = aso.txn_tiles  # 타일 목록
    result = run_aso_design.apply_async(args=[tiles, aso.maf_dir, aso.ref_asm, aso.query_asm])
    print(result.get())  # Celery Task 결과 가져오기
    print("Elapsed:", time.time() - t0, "sec")  # 실행 시간