"""Minimal example: run the frozen DRXO on one CEC 2014 function and on a simple sphere function.
Build the benchmark library first (see README), or skip the CEC part."""
import os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "drxo")); sys.path.insert(0, os.path.join(ROOT, "benchmarks"))
import numpy as np
from drxo import drxo_u_final as drxo

# 1) any Python objective: minimise f over the box [lb, ub]^dim
sphere = lambda x: float(np.sum(x ** 2))
x, fx, info = drxo(sphere, -100, 100, 10, max_fe=100000, seed=1)
print("sphere D=10:", fx)

# 2) CEC 2014 F17 at D = 30 with the held-out seed of run 0 (seed = 14000 + 100*f + r)
try:
    from cec14off import make
    f, bias = make(17, 30)
    x, fx, info = drxo(f, -100, 100, 30, max_fe=300000, seed=14000 + 100 * 17 + 0)
    print("CEC 2014 F17 D=30, run 0: error =", fx - bias)
except OSError as e:
    print("CEC 2014 library not built:", e)
