"""
**Numba 지원** 문자열·숫자 겸용 고속 인덱스 라이브러리

Public API
----------
    build_index(a, engine='auto') -> Index
    find_indices(a, b, engine='auto') -> list[int] | numpy.ndarray[int]

Notes
-----
* `engine='auto'` : array 크기·데이터 종류·Numba 설치 여부로 적절 backend 선택
* Numba backend 는 **숫자(int/float)** 뿐만 아니라 **문자열(str)** 도 지원 (Numba 0.56+) 
  * 숫자 → `typed.Dict[int64,int64]`
  * 문자열 → `typed.Dict[unicode,int64]`
"""

import numpy as np

#  Optional Numba backend ----------------------------------------------------
try:
    from numba import njit, types  # type: ignore
    from numba.typed import Dict as NumbaDict  # type: ignore

    _HAS_NUMBA = True
except ModuleNotFoundError:  # pragma: no cover
    _HAS_NUMBA = False

__all__ = ["build_index", "find_indices", "Index"]

#  Core container

class Index:
    """Lookup object – call with a sequence to fetch indices (‑1 if absent)."""

    __slots__ = ("_map", "_engine")

    def __init__(self, mapping, engine):
        self._map = mapping
        self._engine = engine  # 'python' | 'numpy' | 'numba_num' | 'numba_str'

    # ---------------- public API ---------------- #
    def lookup(self, b):
        if self._engine == "numba_num":
            return _lookup_num(d=self._map, b_np=np.asarray(b, dtype=np.int64))
        if self._engine == "numba_str":
            return _lookup_str(d=self._map, b_list=list(b))
        return [self._map.get(x, -1) for x in b]

    __call__ = lookup  # alias
    def __getitem__(self, item):
        return self._map.get(item, -1)

    def __repr__(self):  # pragma: no cover
        return f"<Index engine={self._engine} size={len(self._map)}>"

# ---------------------------------------------------------------------------
#  Builder -------------------------------------------------------------------
# ---------------------------------------------------------------------------

def build_index(a, engine="auto"):
    """Pre‑compute an index structure for *a* (first occurrence per key).

    Parameters
    ----------
    a : sequence‑like (list/ndarray)
    engine : {"auto", "python", "numpy", "numba"}
    """
    if engine == "auto":
        engine = _pick_engine(a)

    # ---------- pure‑python ---------- #
    if engine == "python":
        return Index({v: i for i, v in enumerate(a)}, "python")

    # ---------- numpy backend ---------- #
    arr = np.asarray(a)
    if engine == "numpy":
        uniq, first_pos = _first_pos_numpy(arr)
        return Index(dict(zip(uniq.tolist(), first_pos.tolist())), "numpy")

    # ---------- numba backend ---------- #
    if engine == "numba":
        if not _HAS_NUMBA:
            raise RuntimeError("Numba not installed; choose another engine.")

        kind = arr.dtype.kind
        if kind in {"i", "u", "f"}:  # numeric
            mapping = _build_num(arr.astype(np.int64))
            return Index(mapping, "numba_num")

        # treat as string (unicode/object)
        py_list = [str(x) for x in a]  # ensure Python str
        mapping = _build_str(py_list)
        return Index(mapping, "numba_str")

    raise ValueError(
        f"Unknown engine {engine!r} – choose 'auto', 'python', 'numpy', or 'numba'.")


def find_indices(a, b, engine="auto"):
    """Convenience: build index then look up *b*."""
    return build_index(a, engine).lookup(b)

# ---------------------------------------------------------------------------
#  Numpy helper --------------------------------------------------------------
# ---------------------------------------------------------------------------

def _first_pos_numpy(arr):
    order = arr.argsort(kind="mergesort")  # stable sort to keep first position
    sorted_arr = arr[order]
    mask = np.empty(len(arr), dtype=bool)
    mask[0] = True
    mask[1:] = sorted_arr[1:] != sorted_arr[:-1]
    return sorted_arr[mask], order[mask]

# ---------------------------------------------------------------------------
#  Numba helpers -------------------------------------------------------------
# ---------------------------------------------------------------------------
if _HAS_NUMBA:

    # ---- numeric ---- #
    @njit(cache=True)
    def _build_num(a_np):
        d = NumbaDict.empty(key_type=types.int64, value_type=types.int64)
        for i in range(a_np.size):
            v = int(a_np[i])
            if v not in d:
                d[v] = i
        return d

    @njit(cache=True)
    def _lookup_num(d, b_np):
        out = np.empty(b_np.size, dtype=np.int64)
        for i in range(b_np.size):
            k = int(b_np[i])
            out[i] = d.get(k, -1)
        return out

    # ---- string ---- #
    @njit(cache=True)
    def _build_str(a_list):
        d = NumbaDict.empty(key_type=types.unicode_type, value_type=types.int64)
        for i in range(len(a_list)):
            v = a_list[i]
            if v not in d:
                d[v] = i
        return d

    @njit(cache=True)
    def _lookup_str(d, b_list):
        out = np.empty(len(b_list), dtype=np.int64)
        for i in range(len(b_list)):
            k = b_list[i]
            out[i] = d.get(k, -1)
        return out
else:  # No numba

    def _build_num(*_a, **_k):
        raise RuntimeError("Numba not available – cannot use engine='numba'.")

    _lookup_num = _build_str = _lookup_str = _build_num  # alias to raise

# ---------------------------------------------------------------------------
#  Engine heuristic ----------------------------------------------------------
# ---------------------------------------------------------------------------

def _pick_engine(a):
    n = len(a)
    if n < 1_000_000:
        return "python"
    if _HAS_NUMBA:
        # try numba even for strings (0.56+)
        return "numba"
    return "numpy"
