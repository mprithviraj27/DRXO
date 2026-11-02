"""Unified DRX engine - single-objective driver (development version).

All new switches off  ->  identical to the frozen DRXO v0.7 (verified by test_equiv.py).
DRX-off skeleton      ->  use_drive=False, use_nucleation=False, nuc_adapt=False, gb_nuc_site=False,
                          use_Z=False, use_jmak_phi=False (L-SHADE-type search with v0.7 memories).

New modules (each one switch, screened one at a time):
  temp='adaptive'   Thermal state s driven by the energy released in each pass (success rate),
                    with thermal inertia c_heat: s <- (1 - c_heat) s + c_heat * SR.
                    Temperature theta = th_min + (th_max - th_min) s;  Z = exp(Q (1/theta - 1/th_max)).
                    Couplings (each switchable):
                      mob='thermal'  boundary mobility M ~ N(mu(theta), 0.02), Arrhenius mean
                                     mu = mu_max exp(-Q_M (1/theta - 1))
                      front='thermal' growth-front fraction phi = phi0 exp(-beta_phi s)
                      Z also drives critical strain and nucleus size (as in v0.7, but adaptive).
  r1_mode='boltz'   thermal-fluctuation grain chosen with Boltzmann weight exp(-E / (kappa theta))
          'rank'    exponential rank weights exp(-3 rank / N) (as in L-SRTDE)
  tex_p > 0         texture-guided migration: with probability tex_p a grain exchanges boundary
                    segments in the texture frame (eigenvectors of the covariance of the
                    best half of the grains) instead of the sample frame.
  reheat=True       multi-pass processing: if the grain structure has collapsed (spread < reheat_tol
                    of the initial spread) or the best energy has not dropped for reheat_stall passes,
                    the material is reheated: every grain except the best is re-nucleated in the box,
                    strain and memories are reset.
  srx_frac > 0      static recrystallization: during the final srx_frac of the budget the deformation
                    stops and the best grain grows by a (mu/mu_w, lambda) evolution strategy whose
                    initial shape is the covariance of the surviving grains (CMA-ES update).
  cdrx=True         continuous/discontinuous dual mode: grains whose misorientation to their growth
                    front is below theta_c (low-angle boundary) rotate as subgrains - they sample
                    around the front with the texture covariance - instead of migrating.
"""
import numpy as np

# ----------------------------------------------------------------------------- helpers

def _cauchy_pos(rng, loc, scale=0.1):
    v = loc + scale * np.tan(np.pi * (rng.random(loc.shape) - 0.5))
    bad = v <= 0
    while bad.any():
        v[bad] = loc[bad] + scale * np.tan(np.pi * (rng.random(bad.sum()) - 0.5))
        bad = v <= 0
    return np.minimum(v, 1.0)


def _norm_trunc(rng, mu, sd, n):
    v = rng.normal(mu, sd, n)
    bad = (v <= 0) | (v > 1)
    while bad.any():
        v[bad] = rng.normal(mu if np.isscalar(mu) else mu[bad], sd, bad.sum()); bad = (v <= 0) | (v > 1)
    return v


def read_shockley(theta, theta_m=1.0):
    t = np.clip(theta / theta_m, 1e-12, 1.0)
    return t * (1.0 - np.log(t))


def texture_frame(X, F, frac=0.5):
    """Eigenvectors of the covariance of the best grains (crystallographic texture)."""
    k = max(X.shape[1] + 1, int(round(frac * len(X))))
    k = min(k, len(X))
    Xb = X[np.argsort(F)[:k]]
    C = np.cov(Xb, rowvar=False) + 1e-30 * np.eye(X.shape[1])
    w, B = np.linalg.eigh(C)
    return B, np.maximum(w, 0.0), C


def cma_polish(f, x0, sigma0, C0, lb, ub, budget, rng, fx0):
    """(mu/mu_w, lambda)-CMA-ES started at x0 with shape C0 (static recrystallization)."""
    n = len(x0); lam = 4 + int(3 * np.log(n)); mu = lam // 2
    w = np.log(mu + 0.5) - np.log(np.arange(1, mu + 1)); w /= w.sum(); mueff = 1 / (w ** 2).sum()
    cc = (4 + mueff / n) / (n + 4 + 2 * mueff / n); cs = (mueff + 2) / (n + mueff + 5)
    c1 = 2 / ((n + 1.3) ** 2 + mueff); cmu = min(1 - c1, 2 * (mueff - 2 + 1 / mueff) / ((n + 2) ** 2 + mueff))
    damps = 1 + 2 * max(0, np.sqrt((mueff - 1) / (n + 1)) - 1) + cs; chiN = np.sqrt(n) * (1 - 1 / (4 * n) + 1 / (21 * n * n))
    m = x0.copy(); sigma = sigma0
    tr = np.trace(C0) / n
    C = C0 / tr if tr > 0 else np.eye(n)
    C = (C + C.T) / 2 + 1e-12 * np.eye(n)
    pc = np.zeros(n); ps = np.zeros(n); used = 0; bx, bf = x0.copy(), fx0; g = 0
    while used + lam <= budget:
        try:
            Dw, B = np.linalg.eigh(C); Dw = np.sqrt(np.maximum(Dw, 1e-20))
        except np.linalg.LinAlgError:
            break
        Z = rng.standard_normal((lam, n)); Y = (Z * Dw) @ B.T
        Xs = m + sigma * Y
        Xs = np.where(Xs < lb, (lb + m) / 2, Xs); Xs = np.where(Xs > ub, (ub + m) / 2, Xs)
        Fs = np.array([f(x) for x in Xs]); used += lam
        o = np.argsort(Fs)
        if Fs[o[0]] < bf: bf, bx = Fs[o[0]], Xs[o[0]].copy()
        Yw = (Y[o[:mu]] * w[:, None]).sum(0)
        m = m + sigma * Yw
        Cinvsqrt = (B / Dw) @ B.T
        ps = (1 - cs) * ps + np.sqrt(cs * (2 - cs) * mueff) * (Cinvsqrt @ Yw)
        g += 1
        hs = np.linalg.norm(ps) / np.sqrt(1 - (1 - cs) ** (2 * g)) / chiN < 1.4 + 2 / (n + 1)
        pc = (1 - cc) * pc + hs * np.sqrt(cc * (2 - cc) * mueff) * Yw
        Ymu = Y[o[:mu]]
        C = (1 - c1 - cmu) * C + c1 * (np.outer(pc, pc) + (1 - hs) * cc * (2 - cc) * C) + cmu * (Ymu.T * w) @ Ymu
        C = (C + C.T) / 2
        sigma *= np.exp((cs / damps) * (np.linalg.norm(ps) / chiN - 1))
        if sigma * np.sqrt(Dw.max() ** 2) < 1e-14 or not np.isfinite(sigma): break
    return bx, bf, used


# ----------------------------------------------------------------------------- engine

def drx_so(f, lb, ub, dim, max_fe=30000, seed=None,
           r_N=18, N_min=4, H=6, r_arc=2.0,
           phi_max=0.2, phi_min=0.05, dscale=0.5,
           lam=0.5, th_max=1.0, th_min=0.1, Q=0.3, eps_c0=3.0, m=0.5,
           d0=0.2, n=0.5, k=5.0, p=2.0, fade_pow=2.0,
           use_archive=True, use_nucleation=True, use_drive=True, use_jmak_phi=True,
           nuc_adapt=True, nuc_life=5, nuc_alpha=0.1, nuc_floor=0.05,
           gb_nuc_site=True, theta_m=1.0, use_Z=True,
           # ---- new modules
           temp='schedule', c_heat=1.0, mob='memory', mu_max=0.65, Q_M=0.054, sd_M=0.02,
           front='jmak', phi0=0.7, beta_phi=7.0, sat='jmak',
           r1_mode='uniform', kappa=0.3,
           tex_p=0.0, tex_frac=0.5,
           reheat=False, reheat_tol=1e-6, reheat_stall=0, reheat_min_left=0.1,
           srx_frac=0.0,
           cdrx=False, theta_c=0.1, cdrx_scale=1.0,
           accept='parent', mb_init=0.5, cr_mode='shade', cr_sd=0.1, bound='mid',
           trace=None):
    rng = np.random.default_rng(seed)
    lb = np.full(dim, lb, float); ub = np.full(dim, ub, float); span = ub - lb
    N_init = max(int(r_N * dim), N_min + 2); N = N_init
    X = lb + rng.random((N, dim)) * span
    F = np.array([f(x) for x in X]); fe = N
    strain = np.zeros(N)
    MM = np.full(H, 0.5); MB = np.full(H, mb_init); km = 0
    archive = np.empty((0, dim))
    gF = F.min(); gbest = X[np.argmin(F)].copy()
    nuc_age = np.full(N, -1); nuc_ref = np.zeros(N); rho = 0.5
    s_th = 0.5                        # thermal state (adaptive temperature)
    spread0 = np.linalg.norm(X.std(axis=0)) + 1e-300
    last_imp_gen = 0; gen = 0
    srx_budget = int(srx_frac * max_fe)
    stop_fe = max_fe - srx_budget
    n_reheat = 0; fifo_ptr = 0; age = np.zeros(N)

    while fe < stop_fe:
        tau = fe / max_fe
        Xj = 1 - np.exp(-k * tau ** p)                       # JMAK fraction (time based)
        if temp == 'adaptive':
            theta = th_min + (th_max - th_min) * s_th
        else:
            theta = th_max - (th_max - th_min) * tau
        z = np.exp(Q * (1 / theta - 1 / th_max)) if use_Z else 1.0
        if front == 'thermal':
            phi = phi0 * np.exp(-beta_phi * s_th)
        else:
            phi = phi_max - (phi_max - phi_min) * Xj if use_jmak_phi else 0.11

        fmin, fmax = F.min(), F.max()
        E = (F - fmin) / (fmax - fmin + 1e-300)
        order = np.argsort(F); rank = np.empty(N, int); rank[order] = np.arange(N)
        nfront = max(2, int(round(phi * N)))

        mem = rng.integers(0, H, N)
        if mob == 'thermal':
            mu = mu_max * np.exp(-Q_M * (1 / max(theta, 1e-3) - 1))
            M = _norm_trunc(rng, mu, sd_M, N)
        else:
            M = _cauchy_pos(rng, MM[mem])
        B = np.where(MB[mem] < 0, 0.0, np.clip(rng.normal(MB[mem], cr_sd), 0, 1))

        g = order[rng.integers(0, nfront, N)]
        drive = dscale * (1 + np.maximum(E - E[g], 0.0)) if use_drive else np.full(N, 0.5)
        pool = np.vstack([X, archive]) if (use_archive and len(archive)) else X
        if r1_mode == 'uniform':
            r1 = rng.integers(0, N, N)
        else:
            wsel = np.exp(-3.0 * rank / N) if r1_mode == 'rank' else np.exp(-E / (kappa * theta))
            cw = np.cumsum(wsel); r1 = np.searchsorted(cw, rng.random(N) * cw[-1]); r1 = np.minimum(r1, N - 1)
        r1 = np.where(r1 == np.arange(N), (r1 + 1) % N, r1)
        r2 = rng.integers(0, len(pool), N)
        clash = (r2 == np.arange(N)) | (r2 == r1)
        r2[clash] = rng.integers(0, len(pool), clash.sum())

        V = X + (M * drive)[:, None] * (X[g] - X) + M[:, None] * (X[r1] - pool[r2])
        if bound == 'mid':
            V = np.where(V < lb, (lb + X) / 2, V)
            V = np.where(V > ub, (ub + X) / 2, V)
        else:
            out_ = (V < lb) | (V > ub); V = np.where(out_, lb + rng.random(V.shape) * span, V)
        mask = rng.random((N, dim)) < B[:, None]
        mask[np.arange(N), rng.integers(0, dim, N)] = True
        U = np.where(mask, V, X)
        need_frame = tex_p > 0 or cdrx
        if need_frame and N > 2:
            Bt, wt, Ct = texture_frame(X, F, tex_frac)
        if tex_p > 0 and N > 2:
            tsel = rng.random(N) < tex_p
            if tsel.any():
                Xr = X[tsel] @ Bt; Vr = V[tsel] @ Bt
                mr = rng.random(Xr.shape) < B[tsel][:, None]
                mr[np.arange(tsel.sum()), rng.integers(0, dim, tsel.sum())] = True
                Ut = np.where(mr, Vr, Xr) @ Bt.T
                Ut = np.where(Ut < lb, (lb + X[tsel]) / 2, Ut); Ut = np.where(Ut > ub, (ub + X[tsel]) / 2, Ut)
                U[tsel] = Ut
        if cdrx and N > 2:
            ref = np.sqrt((X.std(axis=0) ** 2).sum()) + 1e-300
            mis = np.linalg.norm(X - X[g], axis=1) / ref
            low = (mis < theta_c) & (g != np.arange(N))
            if low.any():
                L = Bt * np.sqrt(wt)
                Zs = rng.standard_normal((low.sum(), dim))
                Uc = X[g[low]] + cdrx_scale * M[low, None] * (Zs @ L.T)
                Uc = np.where(Uc < lb, (lb + X[low]) / 2, Uc); Uc = np.where(Uc > ub, (ub + X[low]) / 2, Uc)
                U[low] = Uc

        nev = min(N, stop_fe - fe)
        FU = np.full(N, np.inf)
        FU[:nev] = [f(u) for u in U[:nev]]; fe += nev

        age += 1
        imp = FU < F
        SR = imp.sum() / N
        s_th = (1 - c_heat) * s_th + c_heat * SR
        if imp.any():
            w = (F[imp] - FU[imp]); w = w / w.sum()
            sM, sB = M[imp], B[imp]
            MM[km] = np.sum(w * sM ** 2) / np.sum(w * sM)
            if cr_mode == 'smooth':   # bulge fraction actually realised, smoothed memory update
                sBa = mask[imp].mean(axis=1)
                MB[km] = 0.5 * ((np.sum(w * sBa ** 2) / np.sum(w * sBa) if np.sum(w * sBa) > 1e-8 else 1.0) + MB[km])
            else:
                MB[km] = -1.0 if np.sum(w * sB) == 0 else np.sum(w * sB ** 2) / np.sum(w * sB)
            km = (km + 1) % H
            if use_archive and accept == 'parent': archive = np.vstack([archive, X[imp]])
        strain = np.where(imp, strain * lam, strain + np.maximum(E, 0.1))
        acc = FU <= F
        if accept == 'parent' or not acc.any():
            X[acc], F[acc] = U[acc], FU[acc]
        else:
            # the new grain grows into another grain: the oldest slot (fifo) or the most strained grain
            na = int(acc.sum()); tidx = np.where(acc)[0]
            if accept == 'fifo':
                slots = (fifo_ptr + np.arange(na)) % N; fifo_ptr = (fifo_ptr + na) % N
            elif accept in ('age', 'age_keep', 'age_worse'):   # oldest grains carry the largest accumulated deformation
                ag = age + rng.random(N) * 0.5
                if accept == 'age_keep': ag[np.argmin(F)] = -np.inf           # the lowest-energy grain is stable
                slots = np.argsort(-ag, kind='stable')[:na]
                if accept == 'age_worse':      # a new grain cannot consume a grain of lower energy: it grows into its parent
                    slots = np.where(F[slots] >= FU[tidx], slots, tidx)
                    _, first = np.unique(slots, return_index=True); first = np.sort(first)
                    slots, tidx = slots[first], tidx[first]
            elif accept == 'age_tour':        # the new grain consumes the older of its parent and a neighbouring grain
                cand = rng.integers(0, N, na)
                slots = np.where(age[cand] > age[tidx], cand, tidx)
                _, first = np.unique(slots, return_index=True); first = np.sort(first)
                slots, tidx = slots[first], tidx[first]
            else:
                st = strain.copy()
                if accept == 'strain_keep': st[np.argmin(F)] = -np.inf
                slots = np.argsort(-st, kind='stable')[:na]
            Ua, FUa = U[tidx], FU[tidx]
            if use_archive: archive = np.vstack([archive, X[slots]])     # prior boundaries of consumed grains
            X[slots], F[slots] = Ua, FUa
            strain[slots] = 0.0
            nuc_age[slots] = -1; age[slots] = 0
        if F.min() < gF: gF = F.min(); gbest = X[np.argmin(F)].copy(); last_imp_gen = gen
        if nuc_adapt:
            tr = nuc_age >= 0
            grew = tr & (F < nuc_ref)
            nuc_age[tr] += 1
            dead = tr & ~grew & (nuc_age >= nuc_life)
            for _ in range(int(grew.sum())): rho = (1 - nuc_alpha) * rho + nuc_alpha
            for _ in range(int(dead.sum())): rho = (1 - nuc_alpha) * rho
            nuc_age[grew | dead] = -1

        if use_nucleation:
            eps_c = eps_c0 * z ** m
            p_nuc = (1 - Xj) ** fade_pow if sat == 'jmak' else 1.0
            if nuc_adapt: p_nuc *= min(1.0, 2 * max(rho, nuc_floor))
            d_nuc = d0 * z ** (-n) * 2 * X.std(axis=0)
            order = np.argsort(F); rank = np.empty(N, int); rank[order] = np.arange(N)
            for i in np.where(strain >= eps_c)[0]:
                if fe >= stop_fe: break
                if rank[i] == 0 or rng.random() >= p_nuc: continue
                if gb_nuc_site and rank[i] > 1:
                    cand = order[:rank[i]]
                    ref = np.sqrt((X.std(axis=0) ** 2).sum()) + 1e-300
                    gam = read_shockley(np.linalg.norm(X[cand] - X[i], axis=-1) / ref, theta_m) + 1e-12
                    j = cand[np.searchsorted(np.cumsum(gam), rng.random() * gam.sum())]
                else:
                    j = order[rng.integers(rank[i])]
                xn = np.clip(X[i] + rng.random(dim) * (X[j] - X[i])
                             + d_nuc * (2 * rng.random(dim) - 1), lb, ub)
                if nuc_adapt and nuc_age[i] < 0: nuc_age[i], nuc_ref[i] = 0, F[i]
                X[i], F[i] = xn, f(xn); fe += 1; strain[i] = 0.0; age[i] = 0
                if F[i] < gF: gF, gbest = F[i], X[i].copy(); last_imp_gen = gen

        N_new = int(round(N_init + (N_min - N_init) * fe / max_fe))
        if N_new < N:
            keep = np.argsort(F)[:N_new]
            if nuc_adapt:
                gone = np.setdiff1d(np.arange(N), keep)
                for _ in range(int((nuc_age[gone] >= 0).sum())): rho = (1 - nuc_alpha) * rho
            X, F, strain, age = X[keep], F[keep], strain[keep], age[keep]
            fifo_ptr = fifo_ptr % N_new
            nuc_age, nuc_ref = nuc_age[keep], nuc_ref[keep]; N = N_new
        A_max = int(round(r_arc * N))
        if len(archive) > A_max:
            archive = archive[rng.choice(len(archive), A_max, replace=False)]
        gen += 1

        # ---- multi-pass processing: reheat a collapsed or stalled structure
        if reheat and fe < (1 - reheat_min_left) * stop_fe and N > 4:
            spread = np.linalg.norm(X.std(axis=0)) / spread0
            stalled = reheat_stall > 0 and gen - last_imp_gen >= reheat_stall
            if spread < reheat_tol or stalled:
                b = np.argmin(F)
                Xn = lb + rng.random((N, dim)) * span; Xn[0] = X[b]
                nv = min(N - 1, stop_fe - fe)
                Fn = np.full(N, np.inf); Fn[0] = F[b]
                Fn[1:1 + nv] = [f(x) for x in Xn[1:1 + nv]]; fe += nv
                X, F = Xn, Fn
                strain = np.zeros(N); age = np.zeros(N); nuc_age = np.full(N, -1); nuc_ref = np.zeros(N)
                MM[:] = 0.5; MB[:] = mb_init; archive = np.empty((0, dim)); s_th = 0.5
                last_imp_gen = gen; n_reheat += 1
                if F.min() < gF: gF = F.min(); gbest = X[np.argmin(F)].copy()
        if trace is not None: trace.append((fe, gF, theta, phi, SR))

    # ---- static recrystallization (final polishing of the best grain)
    if srx_budget > 0 and fe < max_fe:
        C0 = np.cov(X, rowvar=False) if N > 2 else np.eye(dim)
        sig0 = max(np.sqrt(np.trace(C0) / dim), 1e-12)
        bx, bf, used = cma_polish(f, gbest.copy(), sig0, C0, lb, ub, max_fe - fe, rng, gF)
        fe += used
        if bf < gF: gF, gbest = bf, bx
    return gbest, gF, dict(fe=fe, reheats=n_reheat)


# ======================= FROZEN SINGLE-OBJECTIVE MODE OF THE UNIFIED ENGINE (DRXO-U) =======================
# Frozen 2026-10-24 before any run on the held-out CEC 2014 suite. Selected on CEC 2017 (D10, D30, D50) and
# CEC 2022 (D10) in five screening rounds (uni/SCREEN_LOG.md). Do not change.
# DRX mechanisms: stored-energy driving force; strain-triggered nucleation (Read-Shockley site, must-grow rule,
# critical strain eps_c0 = 6 scaled by Zener-Hollomon); adaptive temperature from released energy driving an
# Arrhenius boundary mobility; grain consumption (new grains consume the oldest grains; no prior-boundary archive).
DRXO_U = dict(temp='adaptive', c_heat=1.0, mob='thermal', accept='age', r_arc=0.0, eps_c0=6.0,
              cr_mode='smooth', mb_init=1.0, cr_sd=0.05)
# DRX-off skeleton: same DE engine and bulge memory, L-SHADE archive, parent replacement, SHADE mobility memory.
SKELETON = dict(use_drive=False, use_nucleation=False, nuc_adapt=False, gb_nuc_site=False, use_Z=False,
                use_jmak_phi=False, cr_mode='smooth', mb_init=1.0, cr_sd=0.05)
ABLATION_SWITCHES = {'use_drive', 'use_nucleation', 'temp', 'mob', 'accept', 'r_arc', 'use_Z'}


def drxo_u_final(f, lb, ub, dim, max_fe=30000, seed=None, **ablation):
    bad = set(ablation) - ABLATION_SWITCHES
    if bad: raise ValueError(f"not an ablation switch: {bad}")
    p = dict(DRXO_U); p.update(ablation)
    return drx_so(f, lb, ub, dim, max_fe=max_fe, seed=seed, **p)


def drxo_u_skeleton(f, lb, ub, dim, max_fe=30000, seed=None):
    return drx_so(f, lb, ub, dim, max_fe=max_fe, seed=seed, **SKELETON)
