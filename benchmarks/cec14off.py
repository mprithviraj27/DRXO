"""Official CEC 2014 functions (Suganthan's C code) via ctypes. Optimum of F_k is 100*k.
HELD-OUT SUITE: do not run any algorithm on it before the SO engine is frozen."""
import os, ctypes, numpy as np
_lib = ctypes.CDLL(os.path.join(os.environ.get("CEC_LIB_DIR", os.path.dirname(os.path.abspath(__file__))), "libcec14.so"))
_lib.cec14_test_func.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                 ctypes.c_int, ctypes.c_int, ctypes.c_int]
def make(fid, dim):
    out = np.zeros(1)
    def f(x):
        x = np.ascontiguousarray(x, dtype=np.float64)
        _lib.cec14_test_func(x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                             out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), dim, 1, fid)
        return out[0]
    return f, 100.0 * fid
