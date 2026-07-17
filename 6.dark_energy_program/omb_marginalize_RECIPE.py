"""
om_b MARGINALIZATION — drop-in module for the VALIDATED barrow_cmb_full2.py pipeline.

STATUS: the r_s calibration question is RESOLVED (see notes at bottom): r_s(z*)=144 Mpc
is CORRECT for the CMB acoustic scale l_A (the '147 Mpc' target was r_drag at the baryon
drag epoch z~1060, a different quantity used for BAO, not for l_A at photon decoupling
z*~1090). Your validated pipeline already hits l_A=300.4 for LCDM, so it is calibrated.
This module only ADDS per-model om_b profiling on top of that validated machinery.

HOW TO USE: replace the two hooks marked >>> HOOK <<< with your pipeline's actual
functions (model_predict and the Chen inverse-cov chi2), then call
  profiled_chi2(model_params, model_name)
for each model instead of the fixed-om_b chi2.
"""
import numpy as np
from scipy.optimize import minimize_scalar

# ---- Chen 2019 LCDM distance-prior data + inverse covariance (from your pipeline) ----
# X_obs = {R, l_A, om_b}. Use YOUR pipeline's exact X_obs and inverse cov.
X_OBS = np.array([1.7493, 301.462, 0.02239])   # your current values (or Chen 1.7502/301.471/0.02236)
# inverse covariance icov: 3x3, from Chen Table (your pipeline already has this)
# >>> HOOK: replace with your actual icov <<<
SIG = np.array([0.0046, 0.090, 0.00015])
CORR = np.array([[1,0.46,-0.66],[0.46,1,-0.33],[-0.66,-0.33,1]])
ICOV = np.linalg.inv(CORR*np.outer(SIG,SIG))

# >>> HOOK: your validated model prediction. Given model params AND om_b, return
#     the predicted {R, l_A, om_b}. This is EXACTLY your barrow_cmb_full2.py
#     R/l_A machinery (r_s with matter+radiation-only E above z*, quad to 1e6,
#     Hu-Sugiyama z*), just with om_b as an explicit argument. <<<
def model_predict(model_params, om_b):
    """
    model_params: whatever your model needs (n/Delta/delta_T, c, Om_m, H0, ...)
    Returns np.array([R_pred, l_A_pred, om_b]).
    PLACEHOLDER — wire to your validated code.
    """
    raise NotImplementedError("Wire to barrow_cmb_full2.py's R/l_A computation")

def chi2_fixed(model_params, om_b=0.02237):
    """CMB chi2 at a FIXED om_b (the old behavior)."""
    d = model_predict(model_params, om_b) - X_OBS
    return d @ ICOV @ d

def profiled_chi2(model_params, om_b_bounds=(0.016, 0.028), warn_on_bound=True):
    """
    CMB chi2 MARGINALIZED (profiled) over om_b. This is the fix.
    1-D minimize_scalar over om_b inside the model's chi2. Cheap.
    Returns (chi2_min, om_b_bestfit).

    BOUNDS: default (0.016, 0.028) is DELIBERATELY WIDE (~+-40 Planck-sigma).
    A too-narrow window (an earlier version used (0.0205,0.0240)) silently CLIPS
    any model whose true best-fit om_b lands outside it and reports a hugely
    inflated chi2 (bug found in bug-check: a model preferring om_b~0.019 was
    clipped to chi2=31327 vs the true profiled 897). Since the whole purpose of
    marginalizing is to remove exactly this kind of spurious, model-dependent
    penalty, the window must be wide enough that the minimum is interior.
    warn_on_bound flags if the minimizer still lands on a bound (widen further).
    """
    lo, hi = om_b_bounds
    res = minimize_scalar(
        lambda ob: (lambda d: d @ ICOV @ d)(model_predict(model_params, ob) - X_OBS),
        bounds=om_b_bounds, method='bounded'
    )
    if warn_on_bound and (abs(res.x - lo) < 1e-4 or abs(res.x - hi) < 1e-4):
        import warnings
        warnings.warn(
            f"profiled om_b={res.x:.5f} sits on a bound {om_b_bounds}; "
            f"minimum may be clipped and chi2={res.fun:.1f} spuriously inflated. "
            f"Widen om_b_bounds and re-run.", RuntimeWarning)
    return res.fun, res.x

# ---- VALIDATION the module is correct (runs without the real pipeline) ----
# Self-consistent mock: if model_predict at om_b_fid IS the data, profiling recovers ~0.
def _self_test(model_predict_mock):
    om_b_fid = 0.02237
    global X_OBS
    saved = X_OBS.copy()
    X_OBS = model_predict_mock(None, om_b_fid)   # data = model at fiducial
    c0, ob0 = profiled_chi2(None)
    X_OBS = saved
    assert c0 < 1e-3, f"self-test failed: profiled chi2={c0}"
    assert abs(ob0-om_b_fid) < 1e-4
    return c0, ob0

if __name__ == "__main__":
    # demonstrate the self-test with a toy linear model_predict
    def toy(_, om_b):
        # toy: R,l_A depend linearly on om_b around fiducial (stand-in for real pipeline)
        return np.array([1.7493 + 50*(om_b-0.02239),
                         301.462 - 8000*(om_b-0.02239),
                         om_b])
    import types
    g = globals(); g['model_predict'] = toy
    c0, ob0 = _self_test(toy)
    print(f"SELF-TEST PASS: profiled chi2={c0:.2e} at om_b={ob0:.5f} (recovers fiducial)")
    print()
    print("USAGE AT HOME:")
    print("  1. wire model_predict() to barrow_cmb_full2.py's R/l_A computation")
    print("  2. for each model, replace chi2_fixed() with profiled_chi2()")
    print("  3. regenerate Table: total_chi2 = BAO_chi2 + profiled_CMB_chi2")
    print("  4. check steeper-helps/flatter-hurts ordering survives the per-model")
    print("     om_b shifts (it likely does qualitatively, but the ~tens-of-chi2")
    print("     non-uniform penalties mean this MUST be checked, not assumed).")
