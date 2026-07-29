#!/usr/bin/env python3
"""Coupled QNM solver: regular inner BC + outgoing outer BC.

Debug version with robust mismatch computation.
"""

from __future__ import annotations
import json, os, sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize, differential_evolution

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, _PROJECT_ROOT)
sys.path.insert(0, os.path.join(_PROJECT_ROOT, "scripts"))
from step_34_solve_interior import (
    hayward_mass, hayward_N, hayward_lapse_sq, hayward_dNdr, gb_invariant,
)

M = 1.0; G_PARAM = 1.1; ETA = -0.1; L = 2
ALPHA_GB = ETA * M**2 / 3.0

def F_met(r): return hayward_lapse_sq(r, M, G_PARAM)
def N_met(r): return hayward_N(r, M, G_PARAM)
def Np_met(r): return hayward_dNdr(r, M, G_PARAM)
def m_fn(r):  return hayward_mass(r, M, G_PARAM)

def V_scalar(r):
    if r < 1e-14: return 0.0
    F = F_met(r); m = m_fn(r)
    mp = 3.0*M*G_PARAM**3*r**2/(r**3+G_PARAM**3)**2
    dF = -2.0*(mp*r - m)/r**2
    return F*(L*(L+1)/r**2 + dF/r)

def V_grav(r):
    if r < 1e-14: return 0.0
    F = F_met(r); m = m_fn(r)
    return F*(L*(L+1)/r**2 - 6.0*m/r**3)

def V_mix(r):
    if r < 1e-14: return 0.0
    F = F_met(r)
    ang = L*(L+1)*(L-1)*(L+2)
    G = gb_invariant(r, M, G_PARAM)
    return ALPHA_GB * ang * G * F * r**2 / (48.0*M**2)

# Tortoise grid
_r_grid = np.concatenate([
    np.linspace(1e-6, 0.1, 1000),
    np.linspace(0.1, 1.0, 500),
    np.linspace(1.0, 5.0, 500),
    np.linspace(5.0, 20.0, 500),
    np.linspace(20.0, 100.0, 500),
    np.linspace(100.0, 500.0, 500),
])
from scipy.integrate import cumulative_trapezoid
_integrand = np.array([1.0/N_met(r) for r in _r_grid])
_rs_grid = cumulative_trapezoid(_integrand, _r_grid, initial=0.0)

def ts_fast(r): return float(np.interp(r, _r_grid, _rs_grid))

def ode_scalar(r, y, omega, Vfunc):
    u, du = y
    N = N_met(r); F = N**2
    if F < 1e-30: F = 1e-30
    d2u = -Np_met(r)/N * du - (omega**2 - Vfunc(r))/F * u
    return [du, d2u]

def ode_coupled(r, y, omega):
    uZ, duZ, up, dup = y
    N = N_met(r); F = N**2
    if F < 1e-30: F = 1e-30
    Np = Np_met(r)
    VZ = V_grav(r); Vs = V_scalar(r); Vm = V_mix(r)
    d2uZ = -Np/N * duZ - (omega**2 - VZ)/F * uZ + Vm/F * up
    d2up = -Np/N * dup - (omega**2 - Vs)/F * up + Vm/F * uZ
    return [duZ, d2uZ, dup, d2up]

def mismatch_1d(omega_arr, Vfunc, r_start=0.01, r_far=60.0, r_match=4.0):
    omega = complex(omega_arr[0], omega_arr[1])
    l1 = L + 1
    # Inner: u ~ r^{l+1}
    y_in = [r_start**l1, l1*r_start**L]
    sol_in = solve_ivp(lambda r,y: ode_scalar(r,y,omega,Vfunc),
                       [r_start, r_match], y_in,
                       method='DOP853', rtol=1e-11, atol=1e-13, max_step=0.2)
    if not sol_in.success: return 1e6
    # Outer: u ~ e^{i*omega*r*}
    rs_far = ts_fast(r_far)
    N = N_met(r_far)
    phase = np.exp(1j*omega*rs_far)
    y_out = [phase, 1j*omega*N*phase]
    sol_out = solve_ivp(lambda r,y: ode_scalar(r,y,omega,Vfunc),
                        [r_far, r_match], y_out,
                        method='DOP853', rtol=1e-11, atol=1e-13, max_step=0.2)
    if not sol_out.success: return 1e6
    # Normalised Wronskian
    v_in = np.array([sol_in.y[0,-1], sol_in.y[1,-1]], dtype=complex)
    v_out = np.array([sol_out.y[0,-1], sol_out.y[1,-1]], dtype=complex)
    # Normalise
    n_in = np.linalg.norm(v_in); n_out = np.linalg.norm(v_out)
    if n_in < 1e-30 or n_out < 1e-30: return 1e6
    v_in /= n_in; v_out /= n_out
    # Mismatch: |sin(angle between vectors)|
    overlap = abs(np.vdot(v_in, v_out))
    return 1.0 - overlap

def mismatch_coupled(omega_arr, r_start=0.01, r_far=60.0, r_match=4.0):
    omega = complex(omega_arr[0], omega_arr[1])
    l1 = L + 1
    # Two inner solutions
    y1_in = [r_start**l1, l1*r_start**L, 0.0, 0.0]
    y2_in = [0.0, 0.0, r_start**l1, l1*r_start**L]
    sol1 = solve_ivp(lambda r,y: ode_coupled(r,y,omega),
                     [r_start, r_match], y1_in,
                     method='DOP853', rtol=1e-11, atol=1e-13, max_step=0.2)
    sol2 = solve_ivp(lambda r,y: ode_coupled(r,y,omega),
                     [r_start, r_match], y2_in,
                     method='DOP853', rtol=1e-11, atol=1e-13, max_step=0.2)
    if not sol1.success or not sol2.success: return 1e6
    # Two outer solutions
    rs_far = ts_fast(r_far)
    N = N_met(r_far)
    phase = np.exp(1j*omega*rs_far)
    y1_out = [phase, 1j*omega*N*phase, 0.0, 0.0]
    y2_out = [0.0, 0.0, phase, 1j*omega*N*phase]
    sol3 = solve_ivp(lambda r,y: ode_coupled(r,y,omega),
                     [r_far, r_match], y1_out,
                     method='DOP853', rtol=1e-11, atol=1e-13, max_step=0.2)
    sol4 = solve_ivp(lambda r,y: ode_coupled(r,y,omega),
                     [r_far, r_match], y2_out,
                     method='DOP853', rtol=1e-11, atol=1e-13, max_step=0.2)
    if not sol3.success or not sol4.success: return 1e6
    # Solution vectors at r_match
    v1_in = np.array(sol1.y[:,-1], dtype=complex)
    v2_in = np.array(sol2.y[:,-1], dtype=complex)
    v1_out = np.array(sol3.y[:,-1], dtype=complex)
    v2_out = np.array(sol4.y[:,-1], dtype=complex)
    # Orthonormalize each pair using QR
    Q_in, _ = np.linalg.qr(np.column_stack([v1_in, v2_in]))
    Q_out, _ = np.linalg.qr(np.column_stack([v1_out, v2_out]))
    # Principal angles: SVD of Q_in^H @ Q_out
    M = Q_in.conj().T @ Q_out
    sv = np.linalg.svd(M, compute_uv=False)
    # Largest singular value = cos(smallest principal angle)
    # QNM condition: subspaces intersect -> cos(theta) = 1
    return 1.0 - sv[0]

def scan_grid(mismatch_func, omega_re_range, omega_im_range, n_re=20, n_im=15, **kw):
    """Scan a grid of omega values to find approximate roots."""
    best_mis = 1e10
    best_omega = None
    for re in np.linspace(*omega_re_range, n_re):
        for im in np.linspace(*omega_im_range, n_im):
            mis = mismatch_func([re, im], **kw)
            if mis < best_mis:
                best_mis = mis
                best_omega = complex(re, im)
    return best_omega, best_mis

def refine(omega_guess, mismatch_func, **kw):
    result = minimize(lambda x: mismatch_func(x, **kw),
                      [omega_guess.real, omega_guess.imag],
                      method='Nelder-Mead',
                      options={'xatol':1e-9, 'fatol':1e-14, 'maxiter':1000})
    return complex(result.x[0], result.x[1]), result.fun

def main():
    print("="*70)
    print("COUPLED QNM SOLVER (regular inner BC + outgoing outer BC)")
    print(f"M={M}, g={G_PARAM}, eta={ETA}, l={L}")
    print("="*70)

    # 1. Scalar channel validation
    print("\n--- 1. Scalar channel (Schwarzschild ref: 0.4836 - 0.0968i) ---")
    print("  Scanning grid...")
    omega_s, mis_s = scan_grid(lambda x, **kw: mismatch_1d(x, V_scalar, **kw),
                               [0.3, 0.7], [-0.20, -0.02], n_re=30, n_im=20)
    print(f"  Grid best: {omega_s:.4f}  mis={mis_s:.4e}")
    if mis_s < 0.5:
        omega_s, mis_s = refine(omega_s, lambda x, **kw: mismatch_1d(x, V_scalar, **kw))
        print(f"  Refined:   {omega_s:.6f}  mis={mis_s:.4e}")

    # 2. Gravitational channel validation
    print("\n--- 2. Gravitational channel (Schwarzschild ref: 0.3737 - 0.0890i) ---")
    print("  Scanning grid...")
    omega_g, mis_g = scan_grid(lambda x, **kw: mismatch_1d(x, V_grav, **kw),
                               [0.2, 0.6], [-0.20, -0.02], n_re=30, n_im=20)
    print(f"  Grid best: {omega_g:.4f}  mis={mis_g:.4e}")
    if mis_g < 0.5:
        omega_g, mis_g = refine(omega_g, lambda x, **kw: mismatch_1d(x, V_grav, **kw))
        print(f"  Refined:   {omega_g:.6f}  mis={mis_g:.4e}")

    # 3. Coupled system
    print("\n--- 3. Coupled polar-scalar system ---")
    print(f"  V_mix at r=1: {V_mix(1.0):.4e}, r=3: {V_mix(3.0):.4e}, r=5: {V_mix(5.0):.4e}")
    print(f"  Geometry: F_min={min(F_met(r) for r in np.linspace(0.1,5,100)):.4f} (no horizon)")
    print(f"  Throat at r~1.39, cavity modes expected (less damped than Schwarzschild)")

    # Search near the uncoupled cavity modes
    print("  Scanning near uncoupled cavity modes...")
    coupled_modes = []
    for guess_center, label in [(omega_g, "near grav"), (omega_s, "near scalar")]:
        re_c = guess_center.real
        im_c = guess_center.imag
        omega_c, mis_c = scan_grid(mismatch_coupled,
                                   [re_c-0.1, re_c+0.1], [im_c-0.05, im_c+0.05],
                                   n_re=20, n_im=15)
        print(f"  Grid ({label}): {omega_c:.4f}  mis={mis_c:.4e}")
        if mis_c < 0.5:
            omega_c, mis_c = refine(omega_c, mismatch_coupled)
            print(f"  Refined ({label}): {omega_c:.6f}  mis={mis_c:.4e}")
            coupled_modes.append((omega_c, mis_c, label))

    # Also do a broad scan
    print("  Broad scan...")
    omega_broad, mis_broad = scan_grid(mismatch_coupled,
                                       [0.2, 0.8], [-0.15, -0.01], n_re=40, n_im=25)
    print(f"  Broad best: {omega_broad:.4f}  mis={mis_broad:.4e}")
    if mis_broad < 0.3:
        omega_broad, mis_broad = refine(omega_broad, mismatch_coupled)
        print(f"  Broad refined: {omega_broad:.6f}  mis={mis_broad:.4e}")
        coupled_modes.append((omega_broad, mis_broad, "broad"))

    # Pick the best coupled mode
    if coupled_modes:
        coupled_modes.sort(key=lambda x: x[1])
        omega_c, mis_c, label_c = coupled_modes[0]
    else:
        omega_c, mis_c, label_c = omega_broad, mis_broad, "broad"
    print(f"\n  Best coupled mode: {omega_c:.6f}  mis={mis_c:.4e}  ({label_c})")
    if len(coupled_modes) > 1:
        print(f"  All coupled modes found:")
        for om, mis, lab in coupled_modes:
            print(f"    {om:.6f}  mis={mis:.4e}  ({lab})")

    # 4. Summary
    print("\n"+"="*70)
    print("SUMMARY")
    print("="*70)
    print(f"  Scalar (regular BC):       {omega_s:.6f}  mis={mis_s:.2e}")
    print(f"  Schwarzschild scalar ref:  0.483644 - 0.096759i")
    print(f"  Gravitational (regular BC):{omega_g:.6f}  mis={mis_g:.2e}")
    print(f"  Schwarzschild grav ref:    0.373672 - 0.088962i")
    print(f"  Coupled (regular BC):      {omega_c:.6f}  mis={mis_c:.2e}")
    print(f"  Exterior benchmark:        0.374076 - 0.089143i  (polar-led)")

    if omega_g and omega_c:
        br = 100*(omega_c.real - omega_g.real)/omega_g.real
        bi = 100*(omega_c.imag - omega_g.imag)/omega_g.imag
        print(f"  Isospectrality breaking:   real={br:.4f}%  imag={bi:.4f}%")

    output = {
        "parameters": {"M": M, "g": G_PARAM, "eta": ETA, "l": L},
        "method": "Coupled shooting, regular inner BC (r^{l+1}), outgoing outer BC",
        "scalar_qnm": {"omega": str(omega_s), "mismatch": float(mis_s)},
        "gravitational_qnm": {"omega": str(omega_g), "mismatch": float(mis_g)},
        "coupled_qnm": {"omega": str(omega_c), "mismatch": float(mis_c)},
    }
    _results_dir = os.path.join(_PROJECT_ROOT, "results")
    os.makedirs(_results_dir, exist_ok=True)
    out_path = os.path.join(_results_dir, "step_29_matrix_leaver.json")
    def jd(obj):
        if isinstance(obj,(complex,np.complex128)): return {"real":float(obj.real),"imag":float(obj.imag)}
        if isinstance(obj,np.ndarray): return obj.tolist()
        if isinstance(obj,(np.floating,np.integer)): return float(obj)
        raise TypeError(str(type(obj)))
    with open(out_path,"w") as f: json.dump(output,f,indent=2,default=jd)
    print(f"\nSaved to {out_path}")

if __name__ == "__main__":
    main()
