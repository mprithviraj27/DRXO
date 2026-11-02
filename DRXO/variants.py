"""Named configurations of the unified engine (single-objective). Screening log is kept in SCREEN_LOG.md."""
OFF = dict(use_drive=False, use_nucleation=False, nuc_adapt=False, gb_nuc_site=False, use_Z=False, use_jmak_phi=False)
V = {
    "v07": {},                       # frozen DRXO v0.7 (reproduced exactly)
    "skel": dict(OFF),               # DRX-off skeleton of v0.7
}
# ---- round 1: one module at a time on top of v0.7
AT = dict(temp='adaptive', c_heat=1.0)
V.update({
    "at_full":  dict(AT, mob='thermal', front='thermal'),
    "at_mob":   dict(AT, mob='thermal'),
    "at_front": dict(AT, front='thermal'),
    "at_nuc":   dict(AT),                         # adaptive Z only drives nucleation
    "boltz":    dict(r1_mode='boltz'),
    "tex50":    dict(tex_p=0.5),
    "srx05":    dict(srx_frac=0.05),
    "cdrx10":   dict(cdrx=True, theta_c=0.1),
    "reheat":   dict(reheat=True, reheat_tol=1e-4),
})
# ---- round 2: grain consumption (new grains consume the oldest grains) on top of adaptive mobility
AM = dict(temp='adaptive', c_heat=1.0, mob='thermal')
V.update({
    "A1":      dict(AM, accept='age'),
    "A1fifo":  dict(AM, accept='fifo'),
    "A1tex":   dict(AM, accept='age', tex_p=0.5),
    "A1boltz": dict(AM, accept='age', r1_mode='boltz'),
    "A1front": dict(AM, accept='age', front='thermal'),
})
# ---- round 3: no prior-boundary archive (consumed grains stay as old grains instead), texture, nucleation threshold
B0 = dict(AM, accept='age', r_arc=0.0, tex_p=0.5)
V.update({
    "B0":       dict(B0),
    "B0eps6":   dict(B0, eps_c0=6.0),
    "B0d0eps":  dict(B0, eps_c0=6.0, d0=0.05),
    "B0nonuc":  dict(B0, use_nucleation=False),
    "B0notex":  dict(B0, tex_p=0.0),
})
# ---- round 4: candidate C1 (texture dropped: no measurable effect), and realised-bulge memory
C1 = dict(AM, accept='age', r_arc=0.0, eps_c0=6.0)
V.update({
    "C1":   dict(C1),
    "C1sm": dict(C1, cr_mode='smooth', mb_init=1.0, cr_sd=0.05),
})
# DRX-off skeleton for the unified engine: same DE engine (realised-bulge memory), L-SHADE archive,
# parent replacement, SHADE mobility memory; no driving force, nucleation, temperature or consumption.
V["skel2"] = dict(OFF, cr_mode='smooth', mb_init=1.0, cr_sd=0.05)
# ---- round 5: consumption only of higher-energy grains (physical driving-force condition)
V["C2"] = dict(C1, cr_mode='smooth', mb_init=1.0, cr_sd=0.05, accept='age_worse')
# ---- ablation of the frozen single-objective mode (DRXO-U = C1sm); one DRX mechanism removed at a time
V.update({
    "abl_nodrive": dict(V["C1sm"], use_drive=False),
    "abl_nonuc":   dict(V["C1sm"], use_nucleation=False),
    "abl_notemp":  dict(V["C1sm"], temp='schedule', mob='memory'),
    "abl_parent":  dict(V["C1sm"], accept='parent'),
    "abl_noZ":     dict(V["C1sm"], use_Z=False),
})
# ---- sensitivity of the frozen configuration (one parameter at a time)
V.update({
    "sens_eps3":   dict(V["C1sm"], eps_c0=3.0),
    "sens_eps9":   dict(V["C1sm"], eps_c0=9.0),
    "sens_mu055":  dict(V["C1sm"], mu_max=0.55),
    "sens_mu075":  dict(V["C1sm"], mu_max=0.75),
    "sens_heat03": dict(V["C1sm"], c_heat=0.3),
})
# ---- L-SRTDE-type DRX-free baseline (success-rate scale factor, cyclic newest-population overwrite, rank-based r1,
#      success-rate elite fraction 0.7 exp(-7 SR), smoothed CR memory, no archive, no driving force/strain/nucleation)
V["srskel"] = dict(OFF, temp='adaptive', c_heat=1.0, mob='lsrtde', accept='fifo', r_arc=0.0, front='thermal', phi0=0.7,
                   beta_phi=7.0, r1_mode='rank', cr_mode='smooth', mb_init=1.0, cr_sd=0.05)
# ---- factorial decomposition between SR-skel (L-SRTDE-type) and DRXO (C1sm), within one engine
SR = V["srskel"]; DX = V["C1sm"]
V.update({
    # add one DRXO element to SR-skel
    "F_sr_front": dict(SR, front='jmak', use_jmak_phi=True),          # JMAK growth-front schedule instead of 0.7 exp(-7 SR)
    "F_sr_r1":    dict(SR, r1_mode='uniform'),                         # uniform r1 instead of rank-based
    "F_sr_map":   dict(SR, mob='thermal'),                             # Arrhenius mobility instead of tanh map
    "F_sr_age":   dict(SR, accept='age'),                              # age-based consumption instead of cyclic overwrite
    "F_sr_drive": dict(SR, use_drive=True),                            # stored-energy driving force
    "F_sr_nuc":   dict(SR, use_nucleation=True, nuc_adapt=True, gb_nuc_site=True, use_Z=True, eps_c0=6.0),  # strain + nucleation
    # remove one element from DRXO (driving force and nucleation removals are abl_nodrive, abl_nonuc)
    "F_dx_front": dict(DX, front='thermal', phi0=0.7, beta_phi=7.0),
    "F_dx_r1":    dict(DX, r1_mode='rank'),
    "F_dx_map":   dict(DX, mob='lsrtde'),
    "F_dx_age":   dict(DX, accept='fifo'),
})
# ---- post-hoc diagnostics (after the held-out test; exploratory, not part of DRXO)
V["post_late20"] = dict(V["C1sm"], late_sw=0.2)     # parent replacement in the final 20 % of the budget
V["post_rankE"] = dict(V["C1sm"], energy='rank')    # rank-based stored energy
