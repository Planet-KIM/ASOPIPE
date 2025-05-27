import csv
import pyarrow as pa
import pyarrow.csv as pacsv
import polars as pl

def save_csv_std(data_dict, path):
    # 헤더 순서: 사전 키 그대로 사용
    header = list(data_dict)
    n_rows = len(data_dict[header[0]])

    with open(path, "w", newline="", encoding="utf-8", buffering=1<<20) as f:
        w = csv.writer(f)
        w.writerow(header)                 # 헤더 한 줄
        for i in range(n_rows):            # 행-단위 스트리밍
            w.writerow([data_dict[k][i] for k in header])
            


def save_csv_pyarrow(data_dict, path, toString=False):
    if toString:
        # wobble 컬럼의 값들을 문자열로 변환
        for _key, _value in data_dict.items():
            if 'wobble' in _key:
                data_dict[_key] = [str(_value_item) for _value_item in _value ]
    # 1) Arrow Table 생성
    table = pa.Table.from_pydict(data_dict)

    # 2) write_csv – 멀티스레드 옵션은 use_threads
    pacsv.write_csv(
        table,
        path,
        write_options=pacsv.WriteOptions(include_header=True),
    )
    


def save_csv_polars(data_dict, path):
    pl.DataFrame(data_dict).write_csv(path)   # 내부에서 코어 수만큼 병렬 쓰기