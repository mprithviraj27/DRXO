"""CEC 2020 real-world constrained suite (official problems, C++ port by V. Stanovov) with Deb's feasibility rules
encoded in a scalar key: key = f if the solution is feasible, otherwise C + v, where v is the mean constraint violation
of the competition rules (equality tolerance 1e-4) and C = 1e10 exceeds every feasible objective value of RC15-RC33.
Comparing keys therefore reproduces Deb's rules exactly; the optimizer itself is not modified."""
import os, ctypes, numpy as np
from numpy.ctypeslib import ndpointer
_lib = ctypes.CDLL(os.path.join(os.environ.get("CEC_LIB_DIR", os.path.dirname(os.path.abspath(__file__))), "cec2020rwc.so"))
D_ = [9,11,7,6,9,38,48,2,3,3,7,7,5,10,7,14,3,4,4,2,5,9,5,7,4,22,10,10,4,3,4,5,30,118,153,158,126,126,126,76,74,86,86,30,25,25,25,30,30,30,59,59,59,59,64,64,64]
GN = [0,0,14,1,2,0,0,2,1,3,4,9,3,10,11,15,4,4,5,3,8,10,8,7,7,86,3,9,1,8,0,6,30,0,0,0,0,0,0,0,0,0,0,91,24,24,24,29,29,29,14,14,14,14,0,0,0]
HN = [8,9,0,4,4,32,38,0,1,0,4,0,0,0,0,0,0,0,0,0,0,1,3,0,0,0,0,0,0,0,0,0,0,108,148,148,116,116,116,76,74,76,76,0,1,1,1,1,1,1,1,1,1,1,6,6,6]
FSTAR = {15: 2.9944244658E+03, 16: 3.2213000814E-02, 17: 1.2665232788E-02, 18: 5.8853327736E+03, 19: 1.6702177263E+00,
         20: 2.6389584338E+02, 21: 2.3524245790E-01, 22: 5.2576870748E-01, 23: 1.6069868725E+01, 24: 2.5287918415E+00,
         25: 1.6254428092E+03, 26: 3.5359231973E+01, 27: 5.2445076066E+02, 28: 1.4614135715E+04, 29: 2.9648954173E+06,
         30: 2.6138840583E+00, 31: 0.0, 32: -3.0665538672E+04, 33: 2.6393464970E+00}
C = 1e10
dbl = ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
_lib.cec20_func.argtypes = [dbl, ctypes.c_int, dbl, dbl, dbl]; _lib.cec20_func.restype = None
_lib.get_bounds.argtypes = [ctypes.c_int, dbl, dbl]; _lib.get_bounds.restype = None
def maxfes(D): return 100000 if D <= 10 else 200000 if D <= 30 else 400000 if D <= 50 else 800000 if D <= 150 else 1000000
def make(k):
    D, ng, nh = D_[k - 1], GN[k - 1], HN[k - 1]
    lb, ub = np.zeros(D), np.zeros(D); _lib.get_bounds(k, lb, ub)
    fv, gv, hv = np.zeros(1), np.zeros(max(ng, 1)), np.zeros(max(nh, 1))
    def raw(x):
        _lib.cec20_func(np.array(x, dtype=np.float64, copy=True), k, fv, gv, hv)   # the C code may round x in place
        G = np.maximum(gv[:ng], 0).sum(); H = np.where(np.abs(hv[:nh]) - 1e-4 > 0, np.abs(hv[:nh]), 0).sum()
        return fv[0], ((G + H) / (ng + nh) if ng + nh > 0 else 0.0)
    def key(x):
        f, v = raw(x)
        if not np.isfinite(f): return C * 10
        return f if v <= 0 else C + v
    return key, raw, lb, ub, D
def decode(key):
    return (key, 0.0) if key < C else (np.nan, key - C)
