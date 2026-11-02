"""Official CEC 2022 functions (Suganthan's C code) via ctypes. One function per process at a time."""
import os, ctypes, numpy as np
_lib = ctypes.CDLL(os.path.join(os.environ.get("CEC_LIB_DIR", os.path.dirname(os.path.abspath(__file__))), "libcec22.so"))
_lib.cec22_test_func.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                 ctypes.c_int, ctypes.c_int, ctypes.c_int]
BIAS = {1: 300, 2: 400, 3: 600, 4: 800, 5: 900, 6: 1800, 7: 2000, 8: 2200, 9: 2300, 10: 2400, 11: 2600, 12: 2700}
def make(fid, dim):
    out = np.zeros(1)
    def f(x):
        x = np.ascontiguousarray(x, dtype=np.float64)
        _lib.cec22_test_func(x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                             out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), dim, 1, fid)
        return out[0]
    return f, BIAS[fid]
