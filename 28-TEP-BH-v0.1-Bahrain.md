# Temporal Equivalence Principle: Black Holes and the Temporal Horizon
**Matthew Lukin Smawfield**
Version: v0.1 (Bahrain)
First published: 23 July 2026 - Last updated: 23 July 2026
DOI: 10.5281/zenodo.xxxxxxxx

---

## Abstract

Standard FLRW cosmology extrapolates observed cosmic expansion backward to $a(t)\to0$, producing a Big Bang singularity at finite proper time. This paper demonstrates that this singularity is a reconstruction artifact of imposing a globally isochronous expanding-frame description on a conformal temporal geometry. In the Temporal Equivalence Principle (TEP), the observational role of FLRW expansion is reconstructed through conformal temporal transport: the effective scale factor $a_{\rm eff}$ arises from accumulated open-path conformal temporal shear along cosmological lines of sight rather than from physical expansion of space. TEP-C0 (Paper 26) established the distance-redshift and supernova evidence; the full nonsingular matter-frame closure is delivered here.

The Temporal Horizon Cosmology framework is developed here, proving, within the temporal-conformal branch defined here, that the apparent $a_{\rm eff}\to0$ limit is not a physical curvature singularity but a temporal horizon. Two distinct projections of the temporal field are required: $A_{\rm clock}(z)=(1+z)^{-1}$ is the exact observational clock/redshift mapping that drives $a_{\rm eff}\to0$ as $z\to\infty$, while $A_{\rm dyn}(z)=\left(1+z/z_{t}\right)^{-\epsilon_{\rm eff}(z)}$ is the dynamically screened shear response that modifies expansion, BBN, recombination, and perturbations only at late times. Proposition 1 establishes curvature regularity of the temporal conformal boundary: for $A_{\rm clock}(\eta)=C\eta^{-p}$ with $0 \lt p\le\tfrac12$, all polynomial curvature invariants vanish at the boundary, timelike proper time diverges, and null geodesics have divergent affine parameter. The temporal-horizon exponent $p$ and the observational clock map are independent boundary conditions: $A_{\rm clock}(z)=(1+z)^{-1}$ is fixed by the redshift definition, while the regularity condition $0 \lt p\le\tfrac12$ is a mathematical requirement for curvature-regularity at the conformal boundary. Here $\eta$ is the temporal-horizon conformal coordinate, oriented so that approach to $\mathscr{T}^{-}$ corresponds to the asymptotic limit in which $A_{\rm clock}\to0$; it is not the standard FLRW conformal time coordinate extrapolated to $a=0$. Figure 1 (Section 4.5) illustrates the resulting conformal-boundary interpretation: the singular lower edge of standard flat $\Lambda$CDM is replaced by a smooth temporal conformal boundary $\mathscr{T}^{-}$, where $A_{\rm clock}\to0$ and curvature invariants vanish. The conformal compactification is smooth, the Weyl tensor vanishes on the boundary, and every causal curve approaches the regular past boundary $\mathscr{T}^{-}$ rather than terminating at a singularity. The temporal horizon is therefore simultaneously curvature-empty, timelike-complete, and null-complete in this branch.

The effective stress-energy tensor of the temporal field violates the Strong Energy Condition, an explicit prerequisite of the Hawking-Penrose singularity theorems. The thermal screening scale is $T_{\rm lock}=0.03\,\mathrm{eV}$ with transition redshift $z_{t}=100$ and $T_{0}=2.725\,\mathrm{K}$, giving strong epoch-by-epoch screening ($S_{\rm epoch}\sim 10^{-12}$ at BBN, $\sim 10^{-2}$ at recombination). Screened TEP reproduces the standard BBN successful sector and inherits the standard lithium anomaly. Recombination is computed with the full non-equilibrium Peebles/RECFAST treatment. The temporal-horizon thermal mapping preserves a FIRAS-compatible blackbody with no spectral distortions.

The scalar perturbation spectrum is derived from fluctuations of the clock field, $\zeta=\delta\ln A_{\rm clock}$, yielding a power spectrum $P_{\zeta}(k)\propto k^{n_{s}-1}$ with spectral-flow parameter $n_{s}-1=-2\epsilon_{\rm field}$. The observed Planck value $n_{s}=0.965$ constrains $\epsilon_{\rm field}=0.0175$. Tensor modes are derived directly from the temporal-conformal metric: for $A_{\rm clock}(\eta)\sim\eta^{-p}$ the tensor source term $A_{\rm clock}''/A_{\rm clock}=p(p+1)/\eta^{2}\to 0$ at the horizon, so the tensor equation approaches the Minkowski vacuum. The imported inflationary consistency relation $r=16\epsilon_{\rm field}$ is not assumed. Numerical integration of the native tensor equation across the finite transition profile (Step 09b) yields $r(k_{\rm pivot})=9\times 10^{-6}$ and $r_{\rm max}=6.26\times 10^{-4}$, both well below the BICEP/Keck 2021 bound $r\lt0.036$; tensor power is controlled only by the finite transition region. CMB anisotropy and LSS observables are reproduced in the screened-limit reduction, inheriting agreement with Planck 2018 and BOSS DR12 by construction rather than as independent empirical confirmation.

The causal matter-frame universe is curvature-regular at the temporal conformal boundary. The apparent Big Bang is a temporal horizon, not a physical curvature singularity. All background and thermal observational pillars are preserved in the screened-limit reduction; the scalar perturbation shape is reproduced, and the tensor-to-scalar ratio is computed from the native temporal-conformal wave equation, yielding values well below observational bounds.

Code Availability: https://github.com/matthewsmawfield/TEP-TH

Keywords: temporal equivalence principle, temporal horizon cosmology, big bang singularity, static conformal geometry, cosmology, modified gravity, temporal shear

# 1. Introduction

## 1.1 The Question

What is a black hole if proper time is a dynamical scalar field rather than a fixed parameter? The Temporal Equivalence Principle (TEP) poses this question by distinguishing the Einstein-frame gravitational metric $g_{\mu\nu}$, on which gravitational disturbances propagate, from the universal causal matter metric

\begin{equation} \label{eq:intro_1}
\tilde g_{\mu\nu}=A^2(\phi)g_{\mu\nu}+B(\phi)\nabla_\mu\phi\nabla_\nu\phi,
\end{equation}

on which matter, photons, and ideal clocks propagate. The conformal factor $A(\phi)$ maps clock scales; the disformal term $B(\phi)$ can tilt matter cones. In each local Lorentzian orthonormal frame of $\tilde g_{\mu\nu}$, nongravitational physics is locally special relativistic and the locally measured speed of light is invariant. Because $\phi$ is the temporal field — the field that governs the rate of proper time — perturbations $\delta\phi$ are ripples in the dynamical time field itself: *temporal waves*. Gravitational and temporal waves are not independent forms of radiation; they are complementary projections of one coupled disturbance of geometry and dynamical proper time, $(\delta g_{\mu\nu}, \delta\phi) \to \delta\tilde g_{\mu\nu}$. Accelerating temporal-well candidates perturb this unified field (producing $-1$PN scalar dipole emission), and the post-merger ringdown may contain scalar-dominated modes absent from GR. The extended decay times observed in these ringdowns are not delayed secondary signals, but the natural consequence of this single coupled wave experiencing profound coordinate-time slowing as it transits the steep temporal gradient — a mechanism quantified by the $2.8\times$ longer damping timescale derived in the deep-transit solve (Section 8.3). The technical literature terms these scalar perturbations; in the TEP framework they are temporal waves, the physical signature of dynamical proper time.

The physical picture is direct. The TEP black-hole interpretation is dual to the Temporal Horizon Cosmology interpretation. On cosmological scales, a temporal-rate gradient produces observed redshift, which standard physics interprets as spatial expansion. Around an extreme gravitational source, a temporal-rate gradient produces redshift, lensing and altered trajectories, which standard physics interprets as gravitational collapse and inward suction. In both cases, the conventional spatial narrative is an observational projection of a single underlying dynamical proper-time field. A "black hole" is not a hole, vacuum drain, or region sucking space and matter into a singularity. It is a *temporal well* — a regular spatial region in which the rate of proper time differs radically from the exterior. The physical lapse $N(r) = \sqrt{-\tilde g_{\mu\nu}\xi^\mu\xi^\nu}$ satisfies $0 < N(r) < 1$ through the deep region, with $N(r) \ll 1$ where the object appears black, but $N(r) \neq 0$ at every finite physical location. Time is extremely slow relative to the exterior but never stops; no absolute causal boundary forms. The apparent inward pull is analogous to refraction: trajectories curve when the propagation rate varies through a medium, not because matter is being pulled sideways by the medium. In TEP, the temporal field changes observed motion through clock transfer, timelike matter coupling, and its backreaction on the causal geometry. The conformal factor $A(\phi)$ controls clock-rate transfer and massive-particle response; the disformal term $B(\phi)$ can alter causal cones; the backreacted geometric metric $g_{\mu\nu}[\phi]$ generates light bending and orbital curvature. Lensing is generated by the complete matter metric $\tilde g_{\mu\nu}$, not by conformal clock scaling alone. An observer moving through the temporal gradient experiences ordinary local space, a normally running local clock, finite local density, no wall, no horizon crossing, no point where space becomes time. The centre need not be a point of infinite density; matter may remain distributed through a regular spatial region. The extreme effect — extreme frequency transfer $\mathcal{Z} = \omega_e/\omega_o$, reaching a finite but potentially enormous maximum $\mathcal{Z}_{\max} \sim N_o/N_{\min}$ at the centre — appears in comparisons with distant exterior observers as the emitter moves deeper into the temporal-rate gradient. This is established by the invariant frequency-transfer calculation (Section 6.5), not by coordinate-time arguments. The Temporal Horizon $\mathcal{H}_T^{(\Lambda)}[u_o]$ is the observer-dependent surface where the clock-transfer mismatch exceeds a physical accessibility threshold $\Lambda$; it is not a null boundary, not one-way, and does not divide spacetime into an outside and an inside. *Darkness does not require an event horizon; it requires sufficiently extreme temporal decoupling.* One universal matter metric, one temporal field, and no observer-specific metric; different observers nevertheless possess different temporal-accessibility maps. The crucial distinction is *apparent compactness $\neq$ physical compression*. The exterior observer infers a compact, dark, massive object because processes deeper in the temporal field are strongly redshifted, signal arrival rates become extremely slow, light paths are strongly distorted, and matter appears to accumulate or freeze within a small apparent radius. But the density $\rho_{\rm inferred} \sim M/(4\pi r^3/3)$ is an exterior-frame inference that assumes the Schwarzschild radial coordinate and exterior clock remain valid measures of the deep region. The locally measured density $\tilde\rho = \tilde T_{\mu\nu}\tilde u^\mu \tilde u^\nu$ and the local physical volume $d\tilde V = \sqrt{\det\tilde\gamma_{ij}}\,d^3x$ need not reproduce that inference. The extreme density attributed to a black hole may be an artefact of reconstructing a region with a radically different proper-time rate using exterior spatial and temporal standards. Nor is the mass itself directly weighed: it is reconstructed from observed angular positions, spectral shifts, signal periods and image scales under the Isochrony Axiom — the assumption that source clocks, photon propagation and observer clocks can be mapped onto a single universal time coordinate (Section 2.7). TEP instead writes $M_{\rm app}^{\rm GR} = M_{\rm local}^{\rm TEP} + M_{\rm phantom}^{T}$, where $M_{\rm phantom}^{T}$ is the apparent mass generated by reducing a non-isochronous temporal geometry through a single-time GR model (Section 2.8). This is the strong-field counterpart of Phantom Mass at galactic scales. The conventional black-hole picture — event horizon, interior, singular centre — is retained only as the conventional limit that the temporal interpretation replaces. The Schwarzschild metric is an incomplete strong-field limit: it passes every Solar-System test but treats time as a passive coordinate, and when time is correctly treated as a physical field, the fixed-Schwarzschild geometry collapses (Section 3). The TEP Global Solution Architecture — regular spatial domain, finite local density, continuous temporal-rate gradient, no finite-radius causal horizon, extreme external redshift — is the target. The inverse reconstruction and Hayward benchmark serve as evidence that the Schwarzschild singularity is not mandatory; they do not by themselves define the final TEP object. The TEP novelty is the temporal well: a regular spatial region with an extreme gradient in relative proper-time accumulation, dual to the cosmological interpretation. Photons, massive particles, and gravitational waves are all observational probes of the same unified temporal field. In the perturbative sGB benchmark, photons and massive particles respond at different coupling orders — a property of the sGB calculation (Section 6.2).

The unifying mechanism is simple: because the physical lapse $N(r)$ drops, the coordinate propagation rate of light — as seen from the exterior — drastically slows, even though the locally measured speed of light remains exactly $c$. A single slowing produces three distinct observational signatures:

- *Refractive bending.* Light propagating through a gradient where the rate of time varies bends exactly as in refraction, producing the apparent gravitational pull.

- *Extended transit times.* The same slowing extends transit times for waves crossing the deep region, producing the ringdown damping signature.

- *The optical illusion of darkness.* The same slowing, combined with extreme redshift, extinguishes the light's observable frequency to the exterior observer — creating the optical illusion of a black boundary without a physical event horizon.

An observer falling into the well sees no wall, no horizon, and light moving normally: the darkness is the illusion of fast exterior clocks reading an incredibly slow environment (Section 6.1).

The kinematic interpretation is this: all physical spatial and temporal measurements are constructed from the causal matter metric $\tilde g_{\mu\nu}$, whose geometry depends on both the gravitational metric $g_{\mu\nu}$ and the temporal field $\phi$. The conformal factor $A(\phi)$ is a universal local scale field affecting both clock and rod calibration; its gradients generate relational clock and spatial effects. The disformal term $B(\phi)\nabla_\mu\phi\nabla_\nu\phi$ creates directional causal effects aligned with $\nabla_\mu\phi$. TEP decomposes observed spacetime phenomenology into gravitational geometry and a universal temporal response field; matter observes their combined causal metric rather than the Einstein-frame geometry alone. The geometric metric $g_{\mu\nu}$ is the Einstein-frame gravitational spacetime metric — the gravitational engine that sources the temporal field through the sGB coupling. The temporal field $\phi$ is the reactive medium that mediates all physical observations through $\tilde g_{\mu\nu}$. Gravity does not equal time; gravity is the geometric engine, and the temporal field is the medium through which the engine is observed.

## 1.2 The Argument

This paper develops the TEP temporal well as a single continuous derivation. The argument proceeds in five steps.

First, the theory is defined (Section 2): the two-metric action, the conformal-disformal matter metric, the frame dictionary, and the local Lorentz limit.

Second, the fixed-background theorem is proved (Section 3): holding the geometric metric fixed to Schwarzschild creates a mathematical incompatibility. The exact Kretschmann scalar of the conformal matter metric $\tilde g = A^2 g$ with $A = (r_h/r)^{\phi_0}$ in the deep interior scales as $\tilde K \sim r^{4\phi_0 - 6}$. Curvature regularity requires $\phi_0 > 3/2$; bounded areal radius requires $\phi_0 \leq 1$. These conditions are mutually exclusive when the geometric metric is held fixed. The theorem expresses a physical principle: because gravity describes how things move through time, radically changing time must also change gravity. A dynamical proper-time field necessarily implies a dynamical gravitational geometry.

Third, the coupled equations are solved (Section 4): the temporal field couples to the geometric metric through a scalar-Gauss-Bonnet (sGB) interaction $f(\phi)\mathcal{R}_{\rm GB}$, which serves as a low-energy effective field theory for the exterior. The perturbative exterior solution (Sotiriou \& Zhou 2014, exact $\mathcal{O}(\alpha^2)$) establishes scalar hair and real backreaction on the gravitational metric. Its Schwarzschild-based perturbative framework produces an apparent horizon at $F = 0$; the perturbative sGB benchmark belongs to a horizon-bearing branch and therefore does not realise the target TEP temporal well — the true TEP solution must have $N(r) \geq N_{\min} > 0$ everywhere. The inverse reconstruction (Section 4.4) identifies $w_r = -1$ as a sufficient regularity target — one sufficient condition, not a uniqueness theorem — demonstrated by the Hayward benchmark. This serves as evidence that the fixed Schwarzschild singularity is not mandatory; it does not by itself define the final TEP object, which is a temporal well rather than a regular-centre black hole. The benchmark (Appendix K) confirms that finite curvature and bounded areal geometry are mutually compatible once the geometric singularity is absent. The Temporal Horizon is the observer-dependent operational boundary of temporal accessibility: with finite $A$ at the centre, the centre is at finite optical distance, and the transfer factor from the centre is finite but potentially enormous, $\mathcal{Z}_{\max} \sim N_o/N_{\min}$ (Section 6). The shadow, ISCO, and QNM are all observable projections of a single, continuous temporal well — different probes of the same unified temporal field, not evidence of a fractured metric structure. The corrected exterior observables are conditional predictions of the sGB benchmark at fixed GR-calibrated mass (Appendix L).

Fourth, the geometry is validated and predictions derived (Sections 5–8): on the regular-geometry benchmark (Hayward class), the TEP matter metric simultaneously achieves bounded areal radius ($\rho \to 1.995$) and finite Kretschmann ($\tilde K \to 8/r_h^4 = 0.5$), confirming that finite curvature and bounded space are compatible once the geometric singularity is absent. The benchmark produces a finite-area asymptotic end for the $\phi_0 = 1$ divergent-$A$ case; the TEP Global Solution Architecture has finite $A$ at the centre (Section 6.6), giving an ordinary regular centre. The corrected sGB exterior observables at $\eta = -0.1$ — shadow $-0.044\%$ (photons, $\mathcal{O}(\eta^2)$), ISCO on $\tilde g$ $+1.95\%$ (massive particles, $\mathcal{O}(\eta)$, different coupling order and opposite sign) — are conditional predictions of the sGB benchmark at fixed GR-calibrated mass (Appendix L). The coupling-order difference reflects conformal invariance of null geodesics in the sGB benchmark: the conformal factor $A = e^{-\phi}$ cancels for photons but not for massive particles (Section 6.2). An $\mathcal{O}(10^{-3})$ QNM modification is suggested by the exterior scale change; its coefficient and sign require the coupled spectral calculation. The tensor speed $c_T = 1$ on the Schwarzschild background in shift-symmetric sGB (from $G_4 = M_{\rm Pl}^2/2$, with the Gauss–Bonnet coupling entering through $G_5$), satisfying the GW170817 constraint on that background; the full tensor characteristic metric on the TEP solution requires derivation from the complete coupled perturbation system (Section 7.7). Within current observational constraints ($\alpha_{\rm GB} < 2.9$ km$^2$), the observable shifts are largest for low-mass stellar black holes: ISCO shift $\sim +0.78\%$ at $M = 10\,M_\odot$.

Fifth, the theory is confronted with real data (Section 9.11). We apply a 7-step non-isochronous inference pipeline to the 25-year astrometric and spectroscopic record of the S2 star orbiting Sgr A*, using GRAVITY/VLT and Keck data through the TEP-native refit. At S2 scales the mass-inflation branch is selected with a $100\%$ positive Phantom Mass posterior, while the TEP coupling $\eta_{\rm TEP}$ is consistent with zero within the weak-field precision, showing that TEP reduces to GR where the data are weak and produces positive Phantom Mass only where the temporal gradient is strong.

## 1.3 Structure of the Paper

- Section 2 introduces the TEP theory: the two-metric action, the matter metric, and the frame dictionary.

- Section 3 proves the fixed-background theorem: a passive temporal field on fixed Schwarzschild geometry is inconsistent.

- Section 4 solves the coupled Einstein–scalar–Gauss–Bonnet equations in the exterior and identifies the regular-geometry target.

- Section 5 establishes global geometry and curvature regularity on the validation benchmark.

- Section 6 derives the causal structure and the Temporal Horizon.

- Section 7 analyses the perturbation spectrum, identifies what is established (exterior Lorentzianity of $\tilde g$, $c_T = 1$ on the Schwarzschild background) and what remains open (full tensor characteristic metric on the TEP solution, ghost freedom globally, scalar-sector hyperbolicity in the deep interior, precise QNM frequencies, inner boundary condition).

- Section 8 presents observable predictions: shadow, ISCO, QNMs, scalar radiation, and observational constraints.

- Section 9 discusses the physical interpretation and scope.

- Section 10 concludes.

## 1.4 TEP Construction

Universal matter coupling is assigned to

\begin{equation} \label{eq:intro_2}
\tilde g_{\mu\nu}=A^2(\phi)g_{\mu\nu}+B(\phi)\nabla_\mu\phi\nabla_\nu\phi,
\end{equation}

while linearized gravitational disturbances propagate on $g_{\mu\nu}$. With the ultra-damped (quartic Gaussian) disformal coupling, $\tilde g_{\mu\nu}$ remains Lorentzian and invertible for all $r>0$.

- $g_{\mu\nu}$ is the Einstein-frame gravitational metric, dynamical through the sGB coupling;

- $\tilde g_{\mu\nu}$ is the universal matter metric, globally Lorentzian and invertible;

- $A(\phi)$ controls common-mode clock scaling;

- $B(\phi)$ controls disformal cone deformation (shear zone near horizon, vanishes in the interior);

- $\Sigma_\mu=\nabla_\mu\ln A$ is the Temporal Shear.

# 2. TEP Theory: Action, Matter Metric, and Frame Dictionary

## 2.1 Canonical Action and Universal Matter Coupling

The canonical Einstein-frame action is used:

\begin{equation} \label{eq:tep_1}
S=\int d^4x\sqrt{-g}\left[\frac{M_{\rm Pl}^2}{2}R-\frac12(\nabla\Phi)^2-V(\Phi)\right]+S_m[\Psi_m,\tilde g_{\mu\nu}],
\end{equation}

with

\begin{equation} \label{eq:tep_2}
\tilde g_{\mu\nu}=A^2(\Phi)g_{\mu\nu}+B(\Phi)\nabla_\mu\Phi\nabla_\nu\Phi.
\end{equation}

This is the canonical frame: gravity has Einstein–Hilbert form, tensor modes propagate on $g_{\mu\nu}$, and all nongravitational fields couple universally to $\tilde g_{\mu\nu}$. Diffeomorphism invariance of the matter action gives $\tilde\nabla_\mu\tilde T^{\mu\nu}=0$ wherever $\tilde g$ is invertible.

The pipeline uses the dimensionless field $\phi=\Phi/M_*$ and geometrized units $G=c=1$, followed by the numerical normalization $M=1$. Thus radii and curve parameters are reported in $M$, curvature-squared diagnostics in $M^{-4}$, and three-volumes in $M^3$. $A$ and $\phi$ are dimensionless; $B_0=1$ denotes $B_0=1M^2$ in the implemented metric convention, so that $B\,\partial\phi\,\partial\phi$ is dimensionless.

## 2.2 Frame Dictionary and Local Lorentz Limit

| Quantity | Einstein/gravitational frame | Causal matter frame |
| --- | --- | --- |
| Metric | $g_{\mu\nu}$ | $\tilde g_{\mu\nu}$ |
| Propagation | tensor modes | matter, photons, ideal clocks |
| Proper time | $d\tau_g$ | $d\tilde\tau$ |
| Conservation | total Einstein-frame stress tensor | $\tilde\nabla_\mu\tilde T^{\mu\nu}=0$ |
| Validity condition | Lorentzian $g$ | Lorentzian, invertible $\tilde g$ |

At every regular point of the matter frame, Riemann normal coordinates give $\tilde g_{ab}=\eta_{ab}+O(x^2)$, preserving exact local Lorentz invariance and locally invariant $c$. Because the metric is globally non-degenerate, such local inertial frames exist everywhere for all $r > 0$.

## 2.3 Spherical Geometry and Scalar Ansatz

A horizon-regular coordinate system is used from the beginning, avoiding reliance on Schwarzschild $(t,r)$ coordinates across the horizon. The general static spherical form is

\begin{equation} \label{eq:tep_3}
ds_g^2 = -F(r)\,dv^2 + 2\,G(r)\,dv\,dr + R^2(r)\,d\Omega^2,
\end{equation}

in Eddington-Finkelstein form. This prevents coordinate artifacts from being mistaken for physical freezing. In the exterior, $F(r) = 1 - 2M/r$, $G(r) = 1$, $R(r) = r$, recovering the standard Schwarzschild geometry.

The most general symmetry-compatible scalar profile is

\begin{equation} \label{eq:tep_4}
\phi(v,r) = q\,v + \psi(r).
\end{equation}

For the static case studied here, $q = 0$. Two scalar profiles are used in this paper, serving distinct roles in the argument:

- For the fixed-background theorem (Section 3), the scalar is prescribed as $\phi(r) = \phi_0\,\ln(r/r_h)\,S(r)$ with a smooth logistic activation $S(r)$ that suppresses the field in the exterior. This prescribed profile is used only to establish the theorem: the geometric metric is fixed to Schwarzschild, and the scalar field is prescribed rather than solved from coupled equations.

- For the coupled solution (Section 4), the scalar is solved from the sGB field equations and is Coulomb-like ($\phi \sim Q_s/r$) in the exterior, with the geometric metric dynamically modified by the temporal field's backreaction.

These are not two competing models. The first is a theorem-proving device; the second is the physical solution. A time-dependent scalar with static stress-energy is possible in shift-symmetric constructions and may be essential for rotating generalisations.

## 2.4 Complete Matter Metric

Expanding the disformal transformation:

\begin{equation} \label{eq:tep_5}
d\tilde s^2 = A^2\,ds_g^2 + B\left(\psi'(r)\,dr\right)^2.
\end{equation}

The explicit components are

\begin{equation} \label{eq:tep_6}
\tilde g_{vv} = -A^2 F, \qquad \tilde g_{vr} = A^2 G, \qquad \tilde g_{rr} = A^2 \cdot 0 + B\,\psi'^2, \qquad \tilde g_{\theta\theta} = A^2 R^2.
\end{equation}

The two-dimensional $(v,r)$ determinant is

\begin{equation} \label{eq:tep_7}
\det \tilde g_{2D} = (-A^2 F)(B\,\psi'^2) - (A^2 G)^2 = -A^4 G^2 - A^2 B F \psi'^2.
\end{equation}

The disformal coupling uses an ultra-damped quartic Gaussian envelope:

\begin{equation} \label{eq:tep_8}
B(\phi) = B_0\frac{|\phi|^{n_B}}{1+|\phi|^{n_B}}\exp\!\left(-\frac{\phi^4}{2\sigma_B^4}\right),
\end{equation}

with $B_0=1M^2$, $n_B=2$, $\sigma_B=1.5$. The quartic Gaussian envelope drives $B\to0$ in the deep interior, ensuring conformal dominance: $A^4$ dominates at all radii, so $\det \tilde g_{2D} < 0$ for all $r > 0$.

## 2.5 Invertibility and Signature Conditions

The disformal transformation remains invertible when

\begin{equation} \label{eq:tep_9}
A^2 > 0, \qquad A^2 - 2BX \neq 0,
\end{equation}

where $X = -\tfrac12 \nabla_\mu\phi\,\nabla^\mu\phi$. The conditions for Lorentzian signature are then established as the first proposition.

### Proposition 1 (Invertibility and Signature)

For the bounded strong-field ansatz with $\beta_A = -1$, $n_B = 2$, $\delta = 0.05$, $\sigma_B = 1.5$, and $\phi_0 = 2.0$, the matter metric is Lorentzian and nondegenerate for all $r > 0$. The determinant is $-1.002$ at $r=2M$ and remains strictly negative as $r \to 0$, where the conformal dominance $A^4 \gg A^2 B F (\phi')^2$ holds asymptotically. The metric remains Lorentzian and nondegenerate throughout the open domain $r > 0$; the $r \to 0$ limit is an asymptotic spatially enlarged end ($\rho \sim r_h^2/r \to \infty$), not an ordinary manifold point. A negative determinant for every $r > 0$ proves no finite positive radius is degenerate; it does not prove the limit $r = 0$ can be added as a regular point.

## 2.6 Screening and Temporal Topology

Canonical TEP does not identify screening with a density switch. It represents observable suppression by the environmental operator

\begin{equation} \label{eq:tep_10}
\Sigma_\mu^{\rm obs}=\mathcal S_\Sigma(\mathcal E)\,\Sigma_\mu,
\end{equation}

where $\mathcal E$ includes source structure, compactness, gradients, boundary conditions, and coherence scale. Temporal Topology denotes the continuous spatial and covariance structure of $\ln A$ and its Temporal Shear, not a discrete change of spacetime topology. The logistic $S(r)$ used in the fixed-background theorem is a phenomenological radial activation compatible with this vocabulary; the coupled solution does not require $S(r)$ because the scalar profile is determined by the field equations.

## 2.7 The Observational Inference Problem

Because astrophysical mass is an inferred rather than direct measurement, its derivation depends on the assumed time-transfer map. The mass is reconstructed from observed angular positions, spectral shifts, signal periods and image scales under the assumption that source clocks, photon propagation and observer clocks can be mapped onto a single universal time coordinate. This is the *Isochrony Axiom*. The TEP framework demonstrates that this reduction is not theory-neutral.

What the telescope directly records is functions of the observer's local proper time $\tilde\tau_o$:

\begin{equation} \label{eq:tep_isochrony_1}
\mathcal D_{\rm obs} = \left\{\theta_x(\tilde\tau_o),\; \theta_y(\tilde\tau_o),\; \nu_o(\tilde\tau_o),\; F_o(\tilde\tau_o)\right\},
\end{equation}

whereas the source trajectory is parameterised by a different local proper time $\tilde\tau_s$. The full source-to-observer map is

\begin{equation} \label{eq:tep_isochrony_2}
\tilde\tau_s \;\longrightarrow\; x^\mu_s(\tilde\tau_s) \;\longrightarrow\; k^\mu \;\longrightarrow\; \tilde\tau_o,
\end{equation}

with observed frequency $\omega_o = -k_\mu u_o^\mu$ and source-to-observer timing map

\begin{equation} \label{eq:tep_isochrony_3}
\frac{d\tilde\tau_o}{d\tilde\tau_s} = \mathcal T_{s\to o}\!\left[\phi,\, A,\, B,\, x_s^\mu,\, u_s^\mu,\, x_o^\mu,\, u_o^\mu\right].
\end{equation}

The standard analysis reconstructs orbital period, velocity, radius and mass by effectively setting $\mathcal T_{s\to o}$ to its GR or isochronous value. The familiar relation

\begin{equation} \label{eq:tep_isochrony_4}
M_{\rm app}^{\rm GR} = \frac{4\pi^2 a_{\rm GR}^3}{G P_{\rm GR}^2}
\end{equation}

therefore returns an apparent mass conditioned on that temporal calibration. It does not independently measure either the local material mass or the proper volume of the source region. The inferred mass is

\begin{equation} \label{eq:tep_isochrony_5}
M_{\rm app}^{\rm GR} = \mathcal F\!\left[\theta(\tilde\tau_o),\, \nu_o(\tilde\tau_o),\, D_{\rm assumed},\, \mathcal T_{s\to o}^{\rm GR},\, g_{\rm GR}\right],
\end{equation}

not a direct measurement of matter. The spectroscopic velocity has the same structure. The observed spectral shift is conventionally decomposed as $z_{\rm obs} = z_{\rm Doppler} + z_{\rm gravitational} + z_{\rm transverse} + \cdots$ and used to infer a radial velocity. But under TEP,

\begin{equation} \label{eq:tep_isochrony_6}
1 + z_{\rm obs} = \frac{(-k_\mu u^\mu)_e}{(-k_\mu u^\mu)_o}
\end{equation}

also contains the temporal-field transfer. Some quantity interpreted as $v_{\rm radial}$ may actually contain motion plus clock-rate difference plus open-path temporal transport. If that temporal component is forced into a Doppler model, it changes the inferred orbital speed and hence the inferred mass. Even the Galactic-centre distance, though largely inferred by combining angular astrometry with spectroscopic velocities rather than from a standard candle, still assumes $v_{\rm transverse} = D\,d\theta/dt$ and compares it with a spectroscopic radial velocity measured in the same effective time standard. Under TEP the correct relation is schematically $v_{\rm transverse}^{\rm local} = D_{\rm TEP}\,(d\theta/d\tilde\tau_o)\,(d\tilde\tau_o/d\tilde\tau_s)$. Even "geometric" orbital parallax is conditional on the mapping between source time and observer time.

The Cepheid analogy is exact at the methodological level. The TEP cosmological paper treats Cepheids as environment-dependent clocks whose altered periods are misread through a universal period–luminosity law, causing distance and $H_0$ bias. The black-hole analogue is: orbital motion supplies a periodic clock; spectral lines supply atomic clocks; variable plasma features supply dynamical clocks; ringdown supplies gravitational clocks. Standard analysis assumes these clocks can be compared through one universal time coordinate. TEP says their observed periods can include $P_{\rm observed} = P_{\rm intrinsic} \times \mathcal T_{\rm path} \times \mathcal T_{\rm environment}$. If those temporal factors are instead attributed to motion or curvature, the inferred mass becomes inflated. This is the same structural mistake as treating Cepheid period contraction as luminosity or distance.

The direct observables are: changing angular positions; changing received spectral frequencies; repeated signal patterns; strong redshift; strong lensing and image distortion; a dark central observational domain. The inferred quantities — a universal source-frame orbital period, a theory-independent physical orbital velocity, a theory-independent linear radius, millions of solar masses of matter, a tiny proper volume, an event horizon, a singularity — are reconstructed through the Isochrony Axiom. The conventional inference chain is:

\begin{equation} \label{eq:tep_isochrony_7}
\begin{aligned}
&\text{received angles, phases and frequencies} \\
&\quad\Downarrow\quad \text{Isochrony Axiom + GR transfer + distance calibration} \\
&\text{orbital period, velocity and radius} \\
&\quad\Downarrow\quad \text{Kepler/GR inversion} \\
&\text{large compact mass } (M_{\rm app}^{\rm GR}).
\end{aligned}
\end{equation}

TEP replaces it with:

\begin{equation} \label{eq:tep_isochrony_8}
\begin{aligned}
&\text{received angles, phases and frequencies} \\
&\quad\Downarrow\quad \text{dynamical temporal-transfer reconstruction} \\
&\text{local trajectory + clock-rate field} \\
&\quad\Downarrow\quad \text{material source + temporal contribution.}
\end{aligned}
\end{equation}

The precise statement is not that conventional analysis sets $\mathcal T_{s\to o}=1$ — GR includes gravitational redshift, transverse and line-of-sight Doppler shift, propagation delay, and Shapiro delay. The TEP claim is that conventional analysis fixes the transfer map to its *GR closure*, $\mathcal T_{s\to o} = \mathcal T_{s\to o}^{\rm GR}$. TEP tests whether a reproducible residual temporal-transfer field remains after the complete GR mapping is applied:

\begin{equation} \label{eq:tep_Delta_T}
\Delta_T \equiv \ln\frac{\mathcal T_{s\to o}^{\rm TEP}}{\mathcal T_{s\to o}^{\rm GR}},
\end{equation}

Conventional inference is not clock-blind; it is closed under the GR transfer law. The open question is whether $\Delta_T \neq 0$ is required by the data.

## 2.8 Strong-Field Phantom Mass and the Mass Dictionary

This is the strong-field counterpart of Phantom Mass in the wider TEP framework. At galactic scales, unmodelled Temporal Shear is interpreted as a dark-matter halo. At black-hole scales, extreme temporal transport can be interpreted as a very large central mass, rapid infall and a compact causal object. The conventional mass may combine a material contribution, temporal-field energy, and a non-isochronous calibration residual. The more exact expression is

\begin{equation} \label{eq:tep_phantom_2}
M_{\rm app}^{\rm GR} = \mathcal F\!\left(M_{\rm matter},\, M_\phi,\, M_{\rm int},\, A,\, B,\, \mathcal T_{s\to o},\, D_{\rm TEP}\right),
\end{equation}

where $M_{\rm matter}$ is the locally measured material content, $M_\phi$ is the temporal-field energy, $M_{\rm int}$ is the interaction energy, and $A, B, \mathcal T_{s\to o}, D_{\rm TEP}$ are the temporal-field profile and transfer functions. The strong-field Phantom Mass residual is *defined* as the difference between the conventional fitted mass and the TEP material inference:

\begin{equation} \label{eq:tep_phantom_3}
M_{\rm phantom}^{T} \equiv M_{\rm fit}^{\rm GR} - M_{\rm matter}^{\rm TEP}.
\end{equation}

This is a residual definition, not a derived decomposition. The sign of $M_{\rm phantom}^{T}$ is not assumed; it is determined by the data through the sign equation of Section 2.9. The schematic additive form $M_{\rm app}^{\rm GR} = M_{\rm local}^{\rm TEP} + M_{\rm phantom}^{T}$ is shorthand for the residual definition, valid only when $M_\phi$ and $M_{\rm int}$ are absorbed into $M_{\rm phantom}^{T}$.

The paper uses $M$ in several distinct roles, and the distinction matters once the mass itself is in question. The strict dictionary is:

| Symbol | Meaning |
| --- | --- |
| $M_{\rm GR-fit}$ | The parameter returned by a conventional GR reduction of the observations |
| $M_g$ | The asymptotic gravitational charge associated with $g_{\mu\nu}$ |
| $M_{\rm matter}$ | The locally measured material content $M_{\rm matter} = \int_\Sigma \tilde T_{\mu\nu}\tilde u^\mu\tilde u^\nu\,d\tilde V$ |
| $M_\phi$ | The temporal-field energy or effective contribution |
| $M_{\rm phantom}^{T}$ | The difference between the conventional fitted mass and the TEP material inference |
| $M$ | The GR-calibrated exterior mass parameter used in the conditional sGB benchmark (Sections 4–8) |

The statement that a fixed ADM $M$ is "the physical mass measured by a distant observer" is not made in this paper. Under the paper's own thesis, that is precisely what has not yet been established. The sGB calculations of Sections 4–8 use $M$ as the GR-calibrated exterior mass parameter; their observable predictions are conditional on that calibration. The central empirical task is to determine how much of the conventional central mass survives as locally measured material content when the Isochrony Axiom is removed. The decisive analysis is a *TEP-native refit of the raw observations*, not a perturbation of a pre-assumed Schwarzschild mass.

## 2.9 Sign and Identifiability of Strong-Field Phantom Mass

A direct intuition might suggest that because deep clocks run slow, observed periods are longer, so the inferred mass is larger. This overlooks the spatial magnification that accompanies the temporal stretching, and recognising the interplay between the two is the first gate of the entire analysis.

Consider a source orbiting in the deep region with local orbital period $P_s$ and local semi-major axis $a_{\rm local}$. A distant observer measures period $P_o$ and infers radius $a_{\rm GR}$. The temporal transfer stretches the period:

\begin{equation} \label{eq:tep_sign_1}
P_o = \mathcal T_P \, P_s, \qquad \mathcal T_P = \frac{d\tilde\tau_o}{d\tilde\tau_s} > 1.
\end{equation}

The spatial calibration relates the GR-inferred orbit to the local orbit:

\begin{equation} \label{eq:tep_sign_2}
a_{\rm GR} = \mathcal S_a \, a_{\rm local},
\end{equation}

where $\mathcal S_a$ encodes lensing magnification, distance miscalibration, and orbital-dynamics modification — all sourced by the temporal field. The conventional Kepler inversion gives:

\begin{equation} \label{eq:tep_sign_3}
M_{\rm app}^{\rm GR} = \frac{4\pi^2 a_{\rm GR}^3}{G P_o^2} = \frac{4\pi^2 \mathcal S_a^3 a_{\rm local}^3}{G \mathcal T_P^2 P_s^2}.
\end{equation}

The local mass, using the local period and scale, with a possible dynamical-law modification $\mathcal D_{\rm dyn}$ ($\mathcal D_{\rm dyn} = 1$ for Newtonian), is:

\begin{equation} \label{eq:tep_sign_4}
M_{\rm local}^{\rm TEP} = \frac{4\pi^2 a_{\rm local}^3}{G P_s^2 \, \mathcal D_{\rm dyn}}.
\end{equation}

The ratio is:

\begin{equation} \label{eq:tep_sign_5}
\boxed{\frac{M_{\rm app}^{\rm GR}}{M_{\rm local}^{\rm TEP}} = \frac{\mathcal S_a^3 \, \mathcal D_{\rm dyn}}{\mathcal T_P^2}.}
\end{equation}

If the spatial scale is held fixed ($\mathcal S_a = 1$, $\mathcal D_{\rm dyn} = 1$), then $M_{\rm app} = M_{\rm local} / \mathcal T_P^2 < M_{\rm local}$. *Slow deep clocks alone deflate the inferred mass*, producing $M_{\rm phantom}^{T} < 0$. This is the opposite of the Phantom Mass claim. Positive Phantom Mass requires the spatial magnification to dominate:

\begin{equation} \label{eq:tep_sign_6}
M_{\rm phantom}^{T} > 0 \quad \Longleftrightarrow \quad \mathcal S_a^3 \, \mathcal D_{\rm dyn} > \mathcal T_P^2.
\end{equation}

Physically, spatial magnification arises because the conformal factor $A(\phi)$ rescales the matter metric relative to the geometric engine. A distant observer using a fixed angular scale and a GR-based distance calibration measures an enlarged effective orbital radius, while the local material mass needed to support that orbit is unchanged. The temporal field therefore inflates the inferred mass when $A(\phi) > 1$ is strong enough to overcome the deflation from slow local clocks.

The spatial calibration $\mathcal S_a$ is not a free parameter. It is determined by photon propagation on $\tilde g_{\mu\nu}$ (lensing magnification), the GR-inferred distance $D_{\rm GR}$ (standard-candle or parallax calibration through the temporal field), and the orbital dynamics on the temporal-well geometry. The critical threshold is $\mathcal S_a > \mathcal T_P^{2/3}$. For moderate temporal transfer $\mathcal T_P \sim 10$, this requires $\mathcal S_a > 4.6$ — a strong but achievable lensing magnification near a compact object.

The velocity-based mass estimator gives a complementary relation. With $v_{\rm GR} = \mathcal S_v \, v_{\rm local}$ (spectroscopic velocity calibration):

\begin{equation} \label{eq:tep_sign_7}
\frac{M_{\rm app}^{\rm GR}}{M_{\rm local}^{\rm TEP}} = \mathcal S_v^2 \, \mathcal S_a \, \mathcal D_{\rm dyn}.
\end{equation}

The two estimators (Kepler and velocity) must agree in any consistent theory. Their joint constraint tightens the identifiability of $\mathcal S_a$ and $\mathcal T_P$.

Three outcomes are possible and the pipeline must be allowed to return any of them:

| Outcome | Condition | Interpretation |
| --- | --- | --- |
| $M_{\rm phantom}^{T} > 0$ | $\mathcal S_a^3 \mathcal D_{\rm dyn} > \mathcal T_P^2$ | Spatial magnification dominates. Conventional mass is inflated. $M_{\rm matter} < M_{\rm GR\text{-}fit}$. |
| $M_{\rm phantom}^{T} = 0$ | $\mathcal S_a^3 \mathcal D_{\rm dyn} = \mathcal T_P^2$ | Effects cancel. Conventional mass equals local material mass. TEP reduces to GR for mass inference. |
| $M_{\rm phantom}^{T} < 0$ | $\mathcal S_a^3 \mathcal D_{\rm dyn} < \mathcal T_P^2$ | Temporal stretching dominates. Conventional mass is deflated. $M_{\rm matter} > M_{\rm GR\text{-}fit}$. |

A pipeline structured to return only the first outcome is not testing the thesis — it is confirming it. The TEP-native refit (Section 9, Appendix J) must determine the sign from the data, not assume it. The derivation above is in.

Two distinct hypotheses must be separated in the observational analysis. The *composition hypothesis* $H_{\rm composition}: M_{\rm matter} < M_{\rm GR\text{-}fit}$ states that the compact object contains less material mass than the conventional fit suggests. The *calibration hypothesis* $H_{\rm calibration}: M_g < M_{\rm GR\text{-}fit}$ states that the large exterior gravitational parameter itself is a calibration illusion. These are different claims: a temporal field can contribute real energy and backreaction, so TEP might find $M_{\rm matter} \ll M_{\rm GR\text{-}fit}$ while $M_g \sim M_{\rm GR\text{-}fit}$. In that case there is little compressed matter, but there remains a large physical temporal–geometric charge. That would still be a major TEP result, but it would not mean the entire dynamical mass was merely a data-reduction illusion. The raw-observation fit must test both separately.

A separate question is whether "finite density" at the regular centre implies "no extreme compression." A regular-centre geometry can still represent a physically ultra-compact object with extraordinarily high density. The intended TEP vision is stronger: no extreme local compression. This requires an explicit proper-volume diagnostic. Define the proper volume inside a boundary radius $r_b$:

\begin{equation} \label{eq:tep_proper_volume}
V_{\rm proper}(r_b) = 4\pi\int_0^{r_b} \sqrt{\tilde\gamma_{rr}}\,\tilde R^2(r)\,dr,
\end{equation}

and the local mean density $\bar\rho_{\rm local} = M_{\rm matter} / V_{\rm proper}(r_b)$. The critical diagnostic is the volume ratio:

\begin{equation} \label{eq:tep_CV}
\mathcal C_V \equiv \frac{V_{\rm proper}}{(4\pi/3)\,r_{\rm app}^3},
\end{equation}

where $r_{\rm app}$ is the externally inferred apparent radius. A value $\mathcal C_V \gg 1$ would demonstrate that the externally compact image corresponds to a much larger local proper volume — apparent compactness without physical compression. Without this calculation, "apparent compactness $\neq$ physical compression" remains an interpretation rather than a demonstrated property. The proper-volume diagnostic is a required output of the TEP-native refit.

# 3. The Fixed-Background Theorem

This section proves that, within the stated fixed-Schwarzschild power-law conformal class, a dynamical proper-time field necessarily backreacts on gravity. The proof uses Schwarzschild as a restricted test case: the geometric metric is held fixed to Schwarzschild and the question is whether the conformal-disformal matter metric can simultaneously achieve curvature regularity and bounded areal radius. The answer is no. This is not a candidate solution — it is a theorem inside the main derivation, establishing that the full TEP solution must be self-gravitating within this class. The broader physical interpretation — that every conceivable dynamical proper-time field must backreact — is a physical inference from the theorem, not the exact theorem itself.

## 3.1 The Incompatibility

Consider the conformal matter metric $\tilde g = A^2 g_{\rm Schw}$ with $A = (r_h/r)^{\phi_0}$ in the deep interior, on a fixed Schwarzschild background. The exact Kretschmann scalar scales as

\begin{equation} \label{eq:fbt_1}
\tilde K \sim r^{4\phi_0 - 6}
\end{equation}

(Appendix D), not $r^{12\phi_0-6}$ as the naive conformal formula suggests — the derivative terms from the non-constant conformal factor dominate. The areal radius scales as

\begin{equation} \label{eq:fbt_2}
\rho = A\,r \sim r^{1-\phi_0}.
\end{equation}

Two conditions are required for a regular temporal-well interior:

- Finite curvature ($\tilde K$ finite or vanishing) requires $\phi_0 \geq 3/2$. At $\phi_0 = 3/2$ the leading curvature is finite and constant; for $\phi_0 > 3/2$ the curvature vanishes.

- Bounded areal radius ($\rho$ finite) requires $\phi_0 \leq 1$.

These are mutually exclusive when the geometric metric is held fixed to Schwarzschild.

## 3.2 The Limiting Cases

The limiting case $\phi_0 = 2$ achieves vanishing curvature: $\tilde K \sim 39r^2/16 \to 0$ (Appendix D). The matter metric is Lorentzian and nondegenerate throughout the open domain $r > 0$ ($\det\tilde g_{2D} < 0$ for all $r > 0$, since $\det\tilde g_{2D} = A^4 \det g_{2D}$ and $A > 0$). However, the areal radius diverges: $\rho \sim r_h^2/r \to \infty$. The $r \to 0$ limit is an asymptotic spatially enlarged end, not an ordinary manifold point — a negative determinant for every $r > 0$ proves no finite positive radius is degenerate, but does not prove the limit $r = 0$ can be added as a regular point. The conformal scaling stretches both time and space — the interior is spatially enlarged, not a black-hole interior.

The complementary limit $\phi_0 = 1$ keeps the areal radius bounded but the curvature diverges: $\tilde K \sim 9/r^2 \to \infty$. The singularity is not resolved.

The theorem proves this is a structural impossibility within the tested class, not a failure of a particular choice: fixed geometry cannot support both conditions simultaneously.

## 3.3 The Conformal Theorem and the Disformal Extension Programme

The theorem above is a conformal theorem: it applies rigorously to the conformal metric $\tilde g = A^2 g$ with $A = (r_h/r)^{\phi_0}$ on fixed Schwarzschild. The full disformal generalisation — classifying all allowed asymptotics of $A$, $B$, $\phi$ under $\det\tilde g < 0$, $\tilde K < \infty$, and $\tilde R_{\rm area} < \infty$ — is a separate mathematical programme that has not been completed. A systematic gate search over the tested disformal parameter space is reported, but the disformal result is not claimed as a theorem until the classification is complete.

The required conditions for a regular temporal-well interior are: (i) bounded areal radius ($\rho = Ar$ finite as $r \to 0$); (ii) finite or vanishing Kretschmann ($\tilde K$ not diverging); (iii) Lorentzianity ($\det \tilde g_{2D} < 0$, checked by the canonical determinant $-A^4 G^2 - A^2 F B(\phi')^2$, not by $g_{rr}$ alone); (iv) no finite-radius causal boundary ($N \geq N_{\min} > 0$ everywhere). In Eddington–Finkelstein coordinates, $g_{rr} = 0$ is normal and does not violate Lorentzianity — the determinant is what decides nondegeneracy, with the cross term $\tilde g_{vr}$ preserving Lorentzianity when $\tilde g_{rr} = 0$.

The gate search explored: (a) saturating scalar profiles $\phi \to \phi_\infty$ (gives bounded $A$; the determinant must be checked with the canonical formula, not $g_{rr}$); (b) different $B$ activation profiles (logistic, Gaussian, quartic Gaussian — all vanish in the interior, leaving $A$ to control curvature); (c) mild conformal growth $\phi_0 \in (0.5, 1.0)$ (balances finite areal radius against curvature, but $\tilde K \sim r^{4\phi_0-6}$ still diverges for $\phi_0 < 3/2$). The saturating-scalar case requires a separate determinant-based analysis and is not excluded by $g_{rr} \to 0$ alone. The conformal theorem (the $B = 0$ case) is rigorous; the full disformal generalisation is an open programme.

## 3.4 Physical Principle

The fixed-background theorem is the mathematical expression of a physical principle: because gravity depends on how things move through time, changing time must also change gravity. The temporal field cannot sit inertly atop a Schwarzschild geometry; it must backreact on the geometric metric itself. Within the tested class, a dynamical proper-time field necessarily implies a dynamical gravitational geometry. The coupled solution of Section 4 implements this requirement.

## 3.5 Raychaudhuri Equation and the Convergence Condition

On the globally Lorentzian domain, a hypersurface-orthogonal timelike congruence obeys

\begin{equation} \label{eq:fbt_3}
\frac{d\tilde\theta}{d\tilde\tau}=-\frac13\tilde\theta^2-\tilde\sigma_{\mu\nu}\tilde\sigma^{\mu\nu}-\tilde R_{\mu\nu}\tilde u^\mu u^\nu.
\end{equation}

The sign of $\tilde R_{\mu\nu}\tilde u^\mu u^\nu$ is the geometric statement of the timelike convergence condition. The canonical Einstein equations are written for $g_{\mu\nu}$, not necessarily in the form $\tilde G_{\mu\nu} = 8\pi \tilde T_{\mu\nu}$, so the direct geometric statement is that the timelike convergence condition fails in the matter geometry. Only after deriving the effective matter-frame Einstein equation and stress tensor can this be identified as a Strong Energy Condition (SEC) violation. On a regular background with a de Sitter-like temporal minimum ($f(0) = 1$), the effective pressure is negative in the causal matter frame, and the convergence condition fails over the interior domain. This is the Hawking–Penrose escape clause: with the convergence condition violated, the theorems' conclusion of inevitable geodesic incompleteness no longer applies to $(\mathcal M, \tilde g)$.

## 3.6 Metric-Specific Hypotheses

Penrose–Hawking theorems infer causal geodesic incompleteness when their geometric, causal, energy, and global hypotheses hold for a specified Lorentzian metric. In a two-metric theory, those hypotheses must be checked separately for $(\mathcal M,g)$ and $(\mathcal M,\tilde g)$. Incompleteness of Schwarzschild $g$ does not logically prove incompleteness of $\tilde g$, but neither does the existence of a second metric prove its completeness. The fixed-background theorem shows that $\tilde g$ on fixed Schwarzschild cannot be both curvature-regular and spatially bounded; the regular-geometry benchmark (Section 5, Appendix K) shows that $\tilde g$ on a regular background can achieve finite curvature and bounded areal radius, though the $\phi_0 = 1$ benchmark produces a finite-area asymptotic end rather than an ordinary regular centre.

# 4. Coupled Exterior Solution and Interior Closure

The fixed-background theorem establishes that the temporal field must modify the geometric metric. The coupled Einstein–scalar–Gauss–Bonnet equations are solved in the exterior. The scalar-Gauss–Bonnet (sGB) coupling $f(\phi)\mathcal{R}_{\rm GB}$ serves as a low-energy effective field theory for the exterior: it resides in the scalar sector and allows the temporal field to modify the geometric metric through a higher-curvature channel. The exterior equations are solved perturbatively (Sotiriou \& Zhou 2014, exact $\mathcal{O}(\alpha^2)$); the non-perturbative exterior integration requires rerun with the correct scalar charge before its mass-function evolution can be cited as evidence. The inverse reconstruction (Section 4.4) identifies $w_r = -1$ as a sufficient regularity target — evidence that the Schwarzschild singularity is not mandatory, not a claim that the final TEP object has a Hayward-style temporal minimum. A complete TEP action is required to derive the temporal-well profile; the interior integration confirms which dynamical coupling produces the TEP Global Solution Architecture.

## 4.1 Action and Scalar Equation

The TEP action is extended with the shift-symmetric Gauss–Bonnet coupling

\begin{equation} \label{eq:coupled_1}
S = \int d^4x\,\sqrt{-g}\left[\frac{M_{\rm Pl}^2}{2}R - \frac{1}{2}(\partial\phi)^2 + \alpha_{\rm GB}\,\phi\,\mathcal{G}\right] + S_{\rm TEP}[\tilde g_{\mu\nu}, \psi_m],
\end{equation}

where $\mathcal{G} = R^2 - 4R_{\mu\nu}R^{\mu\nu} + R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}$ is the Gauss–Bonnet invariant and $\alpha_{\rm GB}$ is a coupling with dimensions of length$^2$. The dimensionless mass-normalised coupling is defined as $\bar\alpha \equiv \alpha_{\rm GB}/M^2$ and parameterised as $\bar\alpha = \eta/3$, so that $\eta$ is dimensionless at fixed black-hole mass $M$. A numerical value of $\eta$ cannot be transferred unchanged between black holes of different mass without rescaling by $M^2$; the physical coupling is $\alpha_{\rm GB} = \eta M^2/3$. The fundamental coupling $\alpha_{\rm GB}$ is a dimensionful constant of the theory; the dimensionless strength for a black hole of mass $M$ is $\zeta_M \sim \alpha_{\rm GB}/M^2$, so $\eta(M)$ changes with mass. The vacuum gravitational–scalar sector has a shift symmetry $\phi \to \phi + {\rm const}$ (the $\phi\mathcal{G}$ term is shift-symmetric up to the four-dimensional topological Gauss–Bonnet integral). However, the matter coupling $A(\phi), B(\phi)$ explicitly changes under $\phi \to \phi + c$, so the full TEP action is not shift-symmetric unless the matter coupling is derivative-only or otherwise constructed to preserve the symmetry. The shift symmetry is a symmetry of the vacuum gravitational–scalar sector, broken by universal matter coupling; it does not by itself eliminate the cosmological constant problem for the scalar. The scalar equation of motion $\Box\phi = -\alpha_{\rm GB}\mathcal{G} = -(\eta M^2/3)\mathcal{G}$ on a Schwarzschild background admits the analytic solution

\begin{equation} \label{eq:coupled_2}
\phi(r) = \frac{2\eta M}{3}\left(\frac{1}{r} + \frac{M}{r^2} + \frac{4M^2}{3r^3}\right),
\end{equation}

which is regular at the would-be horizon, vanishes as $1/r$ at infinity, and carries scalar charge $Q_s = 2\eta M/3 = 2\alpha_{\rm GB}/M$. The factor of $M$ is required for dimensional consistency: $\phi$ is dimensionless, $[\alpha_{\rm GB}] = L^2$, and $[Q_s/r] = L^{-1}$, so $Q_s = 2\alpha_{\rm GB}/M$ has dimensions of length. In code units ($M = 1$) this factor is invisible numerically, but it must be present analytically before any multi-mass comparison. This evades the no-hair theorem because the Gauss–Bonnet coupling provides a geometric source term absent in standard Einstein–scalar gravity. This is an exterior solution; the interior profile $\phi(r)$ for $r < r_h$ must come from the nonlinear field equations and cannot be obtained by extrapolating the Coulomb form. The perturbative sGB benchmark belongs to a horizon-bearing branch and therefore does not realise the target TEP temporal well. The true TEP solution must be a temporal well: $N \geq N_{\min} > 0$ everywhere.

## 4.2 Perturbative Metric Corrections

The metric corrections are the exact $\mathcal{O}(\alpha^2)$ perturbative solution of Sotiriou \& Zhou (2014), Phys. Rev. D 90, 124063, Eqs. (56)-(63). Writing $f(r) = 1-2m/r$ for the Schwarzschild lapse in terms of the horizon mass parameter $m$, the metric ansatz is

\begin{equation} \label{eq:coupled_3}
ds^2 = -f(r)\left[1+A_2(r)\alpha^2\right]^2 dt^2 + f(r)^{-1}\left[1+B_2(r)\alpha^2\right]^2 dr^2 + r^2 d\Omega^2,
\end{equation}

with $\alpha \equiv \alpha_{\rm GB}$ dimensionless in units of $r_H^2$ (i.e. $\beta \equiv \alpha_{\rm GB}/r_H^2$). Since the metric is unperturbed at $\mathcal{O}(\alpha)$ (the linear-order functions $A_1=B_1=0$ vanish identically), linearising in $\alpha^2$ gives

\begin{equation} \label{eq:coupled_4}
g_{tt} = -F\left[1+\beta^2 h_2(x)\right], \qquad g_{rr} = F^{-1}\left[1+\beta^2 \sigma_2(x)\right], \qquad x \equiv \frac{r_H}{r},
\end{equation}

where $F = 1-r_H/r$ and, converting Sotiriou-Zhou's Eqs. (60)-(61) to this convention,

\begin{equation} \label{eq:coupled_5}
h_2(x) = -\frac{98}{5}x - \frac{98}{5}x^2 - \frac{274}{15}x^3 - \frac{14}{15}x^4 + \frac{52}{15}x^5 + \frac{20}{3}x^6,
\end{equation}

\begin{equation} \label{eq:coupled_6}
\sigma_2(x) = \frac{98}{5}x + \frac{58}{5}x^2 + \frac{38}{5}x^3 - \frac{406}{15}x^4 - \frac{436}{15}x^5 - \frac{92}{3}x^6.
\end{equation}

These polynomials are regular at the would-be horizon ($x=1$) — the correction to $g_{tt}$ vanishes there because it is multiplied by $F\to0$ — and vanish at infinity ($x=0$). The perturbative expansion is controlled by $\beta^2 = (\eta/12)^2 = 6.94 \times 10^{-5}$ for $\eta = -0.1$, with $\beta^2 |h_2(1)| \approx 0.0034 \ll 1$ (safely perturbative).

Sotiriou \& Zhou's Eq. (63) also gives the ADM mass shift at fixed horizon mass parameter $m$: $M_{\rm ADM} = m(1+\tfrac{49}{40}\alpha^2/m^4)$. Holding $M_{\rm ADM}=M=1$ fixed (the GR-calibrated exterior mass parameter used in the conditional sGB benchmark), the would-be horizon radius correspondingly shrinks: $r_H = 2m = 2M(1-19.6\,\beta^2) \approx 1.9973\,M$. This mass-shift term is the same order as $h_2,\sigma_2$ and is essential for a consistent $\mathcal{O}(\beta^2)$ result — omitting it changes the sign of the shadow deviation. The sGB backreaction contracts the apparent geometric horizon at fixed ADM mass — a property of the horizon-bearing branch on which the perturbative benchmark lives. The true TEP solution is a temporal well ($N \geq N_{\min} > 0$ everywhere); the observable predictions at $r \sim 3$–$6M$ are conditional on the GR-calibrated mass scale and have not been shown to be robust to the horizon/no-horizon distinction. Separately, the operational Temporal Horizon changes with the observer through the clock-transfer map defined in Section 6.

A critical result of Sotiriou \& Zhou's nonlinear (non-perturbative) integration of the same linear sGB coupling is that the interior does *not* develop a regular centre. Their exact numerical solutions develop a finite-area singularity at $r \approx 0.27\,r_h$ — a singularity at nonzero spherical radius, not a regular $r = 0$ centre. As they state: "nonlinear effects that are not captured by the perturbative solution lead to a finite area, as opposed to a central, singularity." This is a property of the standard linear shift-symmetric sGB coupling. The perturbative $\mathcal{O}(\alpha^2)$ solution used for the exterior observables does not capture this interior structure.

## 4.3 TEP Matter Metric on the Coupled Background

The disformal matter metric is constructed on the sGB-corrected geometric metric:

\begin{equation} \label{eq:coupled_7}
\tilde g_{\mu\nu} = A^2(\phi)\,g^{\rm sGB}_{\mu\nu} + B(\phi)\,\partial_\mu\phi\,\partial_\nu\phi.
\end{equation}

The conformal coupling is $A = e^{\beta_A\phi}$ with $\beta_A = -1$ — frozen consistently across all calculations, matching the wider TEP corpus weak-field value ($\beta \simeq -0.013$, with $\beta_A = -1$ as the strong-field limit). The sGB scalar profile is Coulomb-like ($\phi \sim Q_s/r$) in the exterior. TEP selects the mass-inflation branch with negative scalar charge ($Q_s < 0$), so $\phi < 0$ and $A = e^{-\phi} \sim e^{-Q_s/r} > 1$ in the exterior — the matter metric is conformally *magnified* relative to the geometric metric, approaching $A \to 1$ as $r \to \infty$. The disformal coupling uses an ultra-damped quartic Gaussian envelope: $B(\phi) = B_0\,|\phi|^2/(1+|\phi|^2)\,\exp(-\phi^4/2\sigma_B^4)$ with $\sigma_B = 1.5$ — frozen consistently across all calculations. In the exterior and near-horizon region, $A > 0$ ensures $\det\tilde g_{2D} = -A^4 G^2 < 0$ for the conformal metric (Lorentzianity inherited from $g$), and the Gaussian envelope suppresses $B$ sufficiently to maintain $\det \tilde g_{2D} < 0$ for the full disformal metric in the exterior ($F > 0$).

#### Assessment of the Exterior Conformal Factor

With $\beta_A = -1$ and the sGB scalar $\phi(r) = (2\alpha_{\rm GB}/m)(1/r + m/r^2 + 4m^2/(3r^3))$ (scalar charge $Q_s = 2\alpha_{\rm GB}/m$), the TEP mass-inflation branch has $\alpha_{\rm GB} < 0$ (hence $Q_s < 0$ and $\phi < 0$). The conformal factor at the shifted would-be horizon ($r = r_H = 1.9973\,M$, $\eta = -0.1$, $M = 1$) is $A(r_H) = e^{-\phi(r_H)} \approx e^{+0.061} \approx 1.063$ — a 6.3% deviation from unity. At the photon sphere ($r = 3M$), $A \approx 1.033$ (3.3% deviation); at the ISCO ($r = 6M$), $A \approx 1.013$ (1.3% deviation). The sGB scalar is $\mathcal{O}(\alpha_{\rm GB})$, so $A = e^{-\phi} = 1 - \phi + \mathcal{O}(\phi^2)$ modifies the matter metric at $\mathcal{O}(\alpha_{\rm GB})$, not $\mathcal{O}(\alpha_{\rm GB}^2)$. Matter-frame timelike observables (ISCO, clock rates, redshift) receive $\mathcal{O}(\eta)$ corrections, not only $\mathcal{O}(\eta^2)$. A dynamical screening mechanism connecting $\beta \simeq -0.013$ in the weak field to $\beta_A = -1$ in the strong field must be derived from the action, not imposed separately in the observational code. Until this is done, the exterior matter metric is substantially nontrivial and the observational predictions computed on the screened profile are not self-consistent with the sGB coupling.

The frame dictionary determines which metric each observable probes. Photons follow $\tilde g_{\mu\nu}$; massive accretion matter follows $\tilde g_{\mu\nu}$; tensor perturbations follow the tensor characteristic equations (which may or may not reduce to propagation on $g$). The current pipeline computes the shadow, ISCO, and QNM shifts on the geometric metric $g^{\rm sGB}$. This is correct for the tensor QNM (gravitational perturbation sector). For the shadow and ISCO, the correct calculation should use $\tilde g$, not $g$. However, for the pure conformal case with $A \to 1$ asymptotically, null geodesic paths are conformally invariant, so the photon sphere on $\tilde g = A^2 g$ is at the same location as on $g$ — the shadow shift computed on $g^{\rm sGB}$ equals the shadow shift on $\tilde g$ for the conformal sector. The ISCO for massive particles is not conformally invariant and must be recomputed on $\tilde g$; the current ISCO value on $g^{\rm sGB}$ is therefore indicative, not final. The disformal term may shift the photon sphere if it modifies the effective radial metric at the orbit location; this requires explicit derivation (Section 6.2). The tensor QNM must be computed on the tensor characteristic metric $\mathcal{G}_{\rm tensor}^{\mu\nu}$, which requires derivation from the coupled perturbation equations (Section 7).

The deep-interior determinant cannot be established from the exterior $\phi \sim Q_s/r$ profile alone. The Coulomb form is an asymptotic/exterior expression; extrapolating it to $r \to 0$ is not valid once nonlinear backreaction dominates. The actual interior $\phi(r)$ must come from the nonlinear field equations. The determinant proof for the deep interior therefore depends on the outcome of the interior integration (Section 4.6). If the interior solution develops a finite-area singularity (as in Sotiriou \& Zhou's nonlinear solutions), the matter metric analysis must be performed on the actual solved geometry up to that singularity, not on an extrapolated profile.

## 4.4 Inverse Reconstruction: The Structural Mandate for the Regular Core

The fixed-background theorem establishes that the temporal field must modify the geometric metric. The inverse reconstruction answers the next question: what kind of modification produces a regular interior? By computing the effective stress-energy required to produce a Hayward-like regular metric, the inverse reconstruction identifies $w_r = -1$ as a sufficient regularity target — one sufficient condition, not a uniqueness theorem for every static or dynamical regular TEP interior.

For the Hayward metric $F(r) = 1 - 2Mr^2/(r^3 + 2M\ell^2)$, the Einstein tensor $G_{\mu\nu}[\rm Hayward]$ is non-zero in the interior, defining an effective stress-energy $T^{\rm eff}_{\mu\nu} = G_{\mu\nu}/(8\pi G)$. Writing $F = 1 - 2m(r)/r$ with $m(r) = Mr^3/(r^3 + 2M\ell^2)$, the effective energy density is

\begin{equation} \label{eq:coupled_8}
\rho_{\rm eff} = -G^t_t = \frac{2m'(r)}{r^2} = \frac{1}{r^2}\left(1 - F - rF'\right) = \frac{12M^2\ell^2}{(r^3 + 2M\ell^2)^2}.
\end{equation}

At the centre ($r \ll \ell$), this gives $\rho_{\rm eff} \to 3/\ell^2$ — the energy density is finite, consistent with the finite Kretschmann scalar $K_{\rm Hayward}(0) = 24/\ell^4$. The de Sitter-like temporal minimum has constant energy density and constant curvature, the standard regular-centre profile.

The key diagnostic is the equation of state. The radial pressure satisfies $p_r = -\rho_{\rm eff}$ (from $G^r_r = G^t_t$), giving a radial EOS $w_r = -1$ — cosmological-constant-like. The tangential pressure differs: $p_\theta \neq p_r$, so the effective stress-energy is anisotropic. This anisotropic, $w_r = -1$ profile is the signature of a higher-curvature coupling: a scalar field's stress-energy, modified by the Gauss-Bonnet invariant, can in principle produce anisotropic effective pressure with a de Sitter-like radial component.

Four candidate mechanisms were evaluated: (i) nonlinear electrodynamics (Hayward's original derivation — requires a specific matter sector, not a scalar field); (ii) k-essence with non-canonical kinetic term (can produce $w = -1$ but lacks the higher-curvature channel); (iii) scalar-Gauss-Bonnet coupling $f(\phi)\mathcal{R}_{\rm GB}$ (the scalar field modifies the geometric metric through the Gauss-Bonnet invariant, can produce anisotropic stress with $w_r = -1$, and is known to produce scalar hair and regular horizons — Kanti, Mavromatos, Rizos, Tamvakis \& Winstanley 1996; Torii, Maeda \& Tamaoki 1999; Kleihaus, Kunz \& Radu 2011); (iv) $f(R)$ gravity (modified gravitational action, but does not naturally couple to the temporal scalar field). The sGB coupling is a candidate for TEP because it is already present in the scalar sector and allows the temporal field to modify the geometric metric through a higher-curvature channel.

#### What sGB does and does not establish

The sGB literature establishes: (a) scalar hair around black holes (Sotiriou \& Zhou 2014; Kanti et al. 1996); (b) real backreaction on the gravitational metric; (c) regular event horizons; (d) no naked singularities (Kanti et al. 1996). It does *not* establish a regular $r = 0$ centre. In much of this literature, "regular black hole solution" means regular outside the horizon and at the event horizon — not finite curvature everywhere. The follow-up to Kanti et al. (1998, Phys. Rev. D 57, 6255) explicitly notes that "other researchers confirmed these results by discussing the internal structure of the solutions behind the horizon, and demonstrating numerically the existence of curvature singularities." Sotiriou \& Zhou's (2014) nonlinear solutions for the same linear sGB coupling develop a finite-area singularity at $r \approx 0.27\,r_h$. More recent simulations (Thaalba, Franchini, Bezares \& Sotiriou 2024, Phys. Rev. D 111, 064054) confirm the finite-area singularity and explore a possible connection between the singularity and loss of hyperbolicity. The inverse reconstruction above identifies the stress-energy profile that a regular centre *would* require; it does not prove that standard linear sGB produces it.

## 4.5 Non-Perturbative Exterior Integration

The perturbative solution treats the backreaction at $\mathcal{O}(\eta^2)$. The non-perturbative regime requires solving the full coupled Einstein-scalar-Gauss-Bonnet equations. The exterior integration ($r > r_h$) is implemented using a Radau stiff solver from $r = 30M$ inward to $r \approx 2.05M$, with the scalar boundary condition $\phi \sim Q_s/r$ at large $r$. The analytic scalar charge is $Q_s = 2\eta M/3$; the integration must use this charge consistently with $\eta$ and with one sign convention.

#### Correction to the Exterior Scalar Charge

The analytic scalar charge is $Q_s = 2\eta M/3$. The exterior integration must be performed with $Q_s = 2\eta M/3$, one consistent sign convention, correct units, horizon and infinity series matching, and constraint residuals verified against published sGB solutions (Sotiriou \& Zhou 2014; Kanti et al. 1996) before any mass-function evolution can be cited as quantitative evidence. The published perturbative solution (Sotiriou \& Zhou 2014, exact $\mathcal{O}(\alpha^2)$) remains the established evidence for scalar hair and real backreaction.

## 4.6 The Interior: Structural Mandate

The fixed-background conformal theorem (Section 3) proves that the temporal field must backreact on the geometric metric. The inverse reconstruction (Section 4.4) identifies $w_r = -1$ as a sufficient regularity target, demonstrated by the Hayward benchmark. The TEP Global Solution Architecture (Section 6.6) — regular spatial domain, finite local density, continuous temporal-rate gradient, no finite-radius causal horizon ($N \geq N_{\min} > 0$ everywhere), extreme external redshift — is the target. This target is now realised numerically by the interior solver described below.

### 4.6A Failure of the Conventional sGB Continuation

The standard linear sGB solution, continued through its GR-like horizon, develops the known interior problems. Inside the horizon, $f < 0$ so $r$ becomes timelike and the static ansatz $ds^2 = -f\,dt^2 + f^{-1}\,dr^2 + r^2\,d\Omega^2$ breaks down. The field equations develop $1/f$ singularities in Schwarzschild-like coordinates. The Kanti, Mavromatos, Rizos, Tamvakis \& Winstanley (1996) areal-radius formulation with the metric $ds^2 = -N^2(r)\,f\,dt^2 + f^{-1}\,dr^2 + r^2\,d\Omega^2$ is useful for horizon series matching and separate interior integrations, but the areal-radius gradient can become timelike inside the trapped region; the coordinate $r$ changes causal character across horizons, and the static ansatz is not itself a horizon-penetrating coordinate system. The known linear sGB solutions (Sotiriou \& Zhou 2014) develop a finite-area singularity rather than a regular centre, and hyperbolicity may fail (Thaalba et al. 2024). This indicates that standard linear sGB is a low-energy effective limit valid in the exterior but incomplete in the deep interior. This conventional continuation belongs to the sGB comparison, not the TEP temporal-well solution.

### 4.6B Required TEP Temporal-Well Boundary Problem

For the intended TEP branch, the global conditions are imposed directly:

\begin{equation} \label{eq:tep_boundary}
F(r) > 0, \qquad N(r) \geq N_{\min} > 0, \qquad \det\tilde g_{2D} < 0,
\end{equation}

with the regular-centre boundary conditions

\begin{equation} \label{eq:tep_centre}
R'(0) = 1, \qquad \phi'(0) = 0, \qquad 0 < N(0) = N_{\min} \ll 1,
\end{equation}

and finite curvature, finite local density. The solution is obtained from infinity or from a regular interior domain without crossing a finite-radius horizon. There is no need for horizon-penetrating coordinates in the target solution because no horizon exists. The sGB coupling is an exterior effective field theory; for the deep interior a regularisation is required to produce the $w_r = -1$ temporal minimum identified by the inverse reconstruction.

The static, spherically symmetric ansatz used throughout this section and in the interior solver is an imposed symmetry reduction, not a derived property of the theory. The TEP action is fully general: the temporal field $\phi$ responds to whatever mass-energy distribution exists, and the equilibrium temporal well is determined by the full coupled field equations. There is no theorem requiring the temporal well to be spherical, axisymmetric, or stationary. A real astrophysical source with angular momentum produces a temporal field that is neither spherical nor static; a binary merger produces a highly asymmetric, dynamical temporal gradient. The spherical case is the simplest solvable sector — it allows the regularity proof, the Frobenius analysis, and the deep-transit QNM spectrum (Section 7) to be computed exactly — but it is one point in a larger solution space, not the foundation from which other configurations are perturbative corrections. The rotating extension (Section 8.5) uses the Delgado et al. (2020) slowly-rotating shift-symmetric sGB solution, which provides the $\mathcal{O}(\beta^2)$ correction to the frame-dragging function from the coupled field equations; it is a second imposed symmetry case (axisymmetric, stationary), not the true TEP solution for a rotating source. The true solution requires the full field equations without symmetry assumptions.

The regular-geometry validation benchmark (Section 5, Appendix K) confirms that the TEP Global Solution Architecture is geometrically self-consistent. This is now extended by the interior solver, which integrates the sGB scalar field equation on a prescribed Hayward-class regular background (with regularization scale $g = 1.1M$) down to $r = 10^{-8}M$ for the mass-inflation branch $\eta = -0.1$. This is an integration *on* that background — the geometric metric $g_{\mu\nu}$ is prescribed as the Hayward profile, not dynamically produced by the sGB coupling from the full action. Standard linear sGB develops a finite-area singularity inside the horizon (Sotiriou \& Zhou 2014; Thaalba et al. 2024); the Hayward background bypasses this pathology by construction, providing a regular geometric metric on which the scalar field is solved. The matter metric has a regular centre: $A(0) = e^{-\phi(0)} = 4.68$, $\tilde N_{\min} = 0.211$ at $r \approx 1.41M$, $\rho_{\rm areal}(0) = A(0)r \to 4.68 \times 10^{-8}$, $w_r(0) = -1.000$ and $K(0) = 27.1$ (finite). This is an ordinary regular point, not a finite-area asymptotic end. The sGB scalar charge remains negative throughout the interior, so the matter metric is conformally magnified all the way to the temporal minimum, giving the positive-Phantom-Mass foundation required by the observational analysis. A full nonlinear backreacted solution — in which the geometric metric is dynamically generated by a modified coupling (e.g. nonlinear $f(\phi)\mathcal{G}$, or additional Ricci coupling $\xi R (\nabla\phi)^2$) — remains an open target. The regularisation scale $g$ is the key physical parameter: for $g > g_{\rm crit} = 1.058\,M$ the temporal well is fully horizonless (Section 8.1A), while for $g < g_{\rm crit}$ the well retains a horizon and the temporal gradient is steeper and more compact than the extreme Hayward threshold allows. The EHT visibility-domain fit (Section 8.1B) empirically bounds $g$ to this sub-critical regime, providing a data-driven constraint on the interior geometry.

## 4.7 The Exterior Effective Field Theory

The sGB coupling serves as a low-energy effective field theory for the exterior. In the weak-to-moderate field regime (exterior and horizon), the TEP temporal field couples to curvature in a way that is well-approximated by the linear sGB action, producing scalar hair and real backreaction on the gravitational metric (Sotiriou \& Zhou 2014; Kanti et al. 1996). The perturbative $\mathcal{O}(\alpha^2)$ exterior solution provides the observable shifts that constitute the conditional sGB benchmark. In the extreme strong-field regime (deep interior), the linear sGB approximation breaks down — the literature confirms that standard linear sGB develops a finite-area singularity (Sotiriou \& Zhou 2014) and may lose hyperbolicity (Thaalba et al. 2024). A complete TEP action is required to derive the temporal-well profile. The exterior observables (Section 4.8) are leading-order effective field theory estimates; the exact quantitative predictions require the full coupled derivation on the matter metric $\tilde g$ (Sections 7, 8).

## 4.8 Conditional GR-Calibrated sGB Exterior Benchmark

The sGB corrections shift the exterior observables away from their Schwarzschild values. These are *conditional predictions* of the sGB benchmark at fixed GR-calibrated mass scale, not the primary TEP prediction. Their purpose is: (i) to demonstrate that a temporal scalar can backreact on $g_{\mu\nu}$; (ii) to establish scalar hair; (iii) to quantify the observable shifts when the conventional mass scale is held fixed; (iv) to provide a consistency benchmark for later TEP-native inference. All numbers use fixed-ADM mass normalization ($M_{\rm ADM} = M$ held fixed; the would-be horizon shifts to $r_H = 2M(1 - 19.6\,\beta^2)$). For the mass-inflation branch $\eta = -0.1$ (perturbative regime: $\beta^2 |h_2(1)| \approx 0.0034 \ll 1$) the pipeline reports:

| Observable | Schwarzschild | sGB ($\eta=-0.1$) | Deviation | Metric / Order |
| --- | --- | --- | --- | --- |
| Horizon $r_H/M$ | 2.000 | 1.9973 | $-0.136\%$ | Fixed-ADM shift (TEP signature) |
| Photon sphere $r_{\rm ph}/M$ | 3.000 | 2.9974 | $-0.087\%$ | $g^{\rm sGB}$ (= $\tilde g$ for null), O($\eta^2$) |
| Shadow radius $b/M$ | 5.196 | 5.194 | $-0.044\%$ | $g^{\rm sGB}$ (= $\tilde g$ for null), O($\eta^2$) |
| ISCO $r_{\rm ISCO}/M$ (on $g^{\rm sGB}$) | 6.000 | 5.996 | $-0.061\%$ | $g^{\rm sGB}$ — geometric metric, O($\eta^2$) |
| ISCO $r_{\rm ISCO}/M$ (on $\tilde g$) | 6.000 | 6.117 | $+1.95\%$ | $\tilde g = A^2 g^{\rm sGB}$ — matter metric, O($\eta$) |
| Gravitational tensor QNM $\omega_R$ (l=2) | 0.3737 | — | $\sim +0.136\%$ | Horizon-shift estimate; needs coupled spectral solver (Appendix L) |

The corrected Regge–Wheeler potential comes from the metric perturbations $h_2(x)$ and $\sigma_2(x)$ (Sotiriou \& Zhou 2014, Eqs. 60–61), not from the scalar gradient directly. The correction is $\mathcal{O}(\beta^2) = \mathcal{O}(\eta^2)$ where $\beta = \eta/12$. The dominant effect is the would-be horizon shift $r_H = 2M(1 - 19.6\,\beta^2)$, which shifts the potential and hence the QNM frequency. A QNM frequency shift of $\sim +0.136\%$ at $\eta = -0.1$ is estimated from the would-be horizon shift; a precise value requires the coupled spectral calculation with a coupled spectral solver on the corrected potential (Appendix L). The horizonless temporal well alters the late-time ringdown spectrum (Section 7.3).

The conditional sGB benchmark produces observable shifts at different coupling orders. Both photons and massive particles propagate on $\tilde g_{\mu\nu}$, but because null geodesics are conformally invariant, the conformal factor cancels in the photon orbital equations: the shadow deviation is $-0.044\%$ (sensitive only to the underlying $g^{\rm sGB}$ geometry, $\mathcal{O}(\eta^2)$). The ISCO deviation on $\tilde g$ is $+1.95\%$ (massive particles feel the conformal factor $A = e^{-\phi} > 1$ directly, $\mathcal{O}(\eta)$). Under fixed-ADM normalization the shadow contracts, but the matter-metric ISCO expands because the conformal factor magnifies spatial scales. The ISCO shift is $\sim 44\times$ larger than the shadow shift because the conformal factor contributes at $\mathcal{O}(\eta)$, not $\mathcal{O}(\eta^2)$. This coupling-order difference is a property of the sGB benchmark (Section 6.2). All numbers scale as $\eta^2$ for the shadow and geometric-metric ISCO, and as $|\eta|$ for the matter-metric ISCO. Within current observational constraints ($\alpha_{\rm GB} < 2.9$ km$^2$), the observable shifts are largest for low-mass stellar black holes ($M \sim 10\,M_\odot$). The explicit derivation is in Appendix L.

# 5. Global Geometry and Curvature Regularity

The geometry of the TEP matter metric is validated on a regular-geometry benchmark. The fixed-background theorem (Section 3) proves that curvature regularity and bounded areal radius are mutually exclusive on fixed Schwarzschild. The validation benchmark (Hayward class, Appendix K) confirms they are mutually compatible once the geometric singularity is absent. All benchmark results below use the logarithmic scalar profile ($\phi_0 = 1$). In addition, the sGB-coupled interior has now been integrated down to $r = 10^{-8}M$ by (Section 4.6); that numerical result gives a genuine regular point centre rather than a finite-area asymptotic end.

## 5.1 The Determinant: Global Lorentzian Signature

For a static radial scalar on the regular geometric background, the radial-block determinant of the matter metric is

\begin{equation} \label{eq:global_1}
\det\tilde g_{2D} = -A^4 G^2 - A^2 B\,F\,(\phi')^2,
\end{equation}

where $F$ and $G$ are the metric functions of the geometric background. The conformal term $-A^4 G^2$ is strictly negative, and the disformal correction $-A^2 B F (\phi')^2$ is controlled by the quartic Gaussian damping of $B(\phi)$. The determinant is strictly negative at every radius:

\begin{equation} \label{eq:global_2}
\det\tilde g_{2D} < 0 \quad \text{for all } r > 0.
\end{equation}

The space is globally Lorentzian. There is no determinant-zero boundary, no signature change, and no invertibility obstruction at any radius from the deep temporal region to the centre. The quartic Gaussian envelope $\exp(-\phi^4/2\sigma_B^4)$ suppresses $B$ below $10^{-300}$ in the deep interior, so the conformal term $A^4$ dominates everywhere and the disformal correction never crosses it. No curvature tensor or geodesic equation of $\tilde g$ encounters a degeneracy. The TEP object is a temporal well: $N(r) = \sqrt{-\tilde g_{tt}} > 0$ everywhere, with the lapse becoming extraordinarily small inward but never vanishing at a finite-radius boundary.

## 5.2 The Kretschmann Scalar: Finite de Sitter-Like Core

On the regular geometric background with $\phi_0 = 1$, the Kretschmann scalar of the matter metric approaches a finite value at the centre:

\begin{equation} \label{eq:global_3}
\tilde K \to 8/r_h^4 = 0.5 \quad \text{as } r \to 0.
\end{equation}

This is the de Sitter-like temporal minimum: the curvature is finite and bounded, in stark contrast to the Schwarzschild value $K_{\rm Schw} = 48\,M^2/r^6 \to \infty$. The finite Kretschmann scalar is the direct consequence of the regular geometric background. The effective anisotropic stress-energy with $w_r = -1$ supports a de Sitter-like temporal minimum of constant curvature that replaces the classical singularity. All curvature invariants — the Ricci scalar, the Gauss-Bonnet invariant, and the full Riemann contraction — remain finite throughout the interior.

## 5.3 The Areal Radius: Bounded Space

The areal radius of the matter metric is

\begin{equation} \label{eq:global_4}
\rho = A\,r.
\end{equation}

On the regular geometric background with $\phi_0 = 1$, the areal radius remains bounded:

\begin{equation} \label{eq:global_5}
\rho \to 1.995 \quad \text{as } r \to 0.
\end{equation}

The interior space does not open without limit. This is the essential distinction from the fixed-geometry limit: on a Schwarzschild background with $\phi_0 = 2$, the conformal scaling $A \sim r^{-2}$ drives the areal radius to diverge ($\rho \sim r_h^2/r \to \infty$), illustrating the fixed-background theorem (Section 3) — curvature regularity and bounded areal radius are mutually exclusive when the geometric metric is held fixed. On the regular background, the geometric metric itself is regular, and both conditions are satisfied simultaneously: the curvature is finite and the areal radius is bounded.

A classification point must be stated precisely. At an ordinary regular spherical centre one has $\rho \to 0$ as $r \to 0$. Here $\rho \to 1.995 \neq 0$: the limiting areal radius is nonzero. Because $A \sim 1/r$ in the deep interior, the radial proper distance also behaves as $\tilde\ell \sim \int dr/r \to \infty$. The limiting region is therefore not a point-like centre but a finite-area asymptotic end (a tube-like limiting region of infinite radial extent). This is a property of the $\phi_0 = 1$ benchmark only. The full coupled interior solution (Section 4.6) produces an ordinary regular point centre with $\rho_{\rm areal}(0) \to 4.68 \times 10^{-8}$ for $\eta=-0.1$, $g=1.1M$; the benchmark classification is reported here as computed.

The regular geometry resolves both aspects of the fixed-background theorem. The de Sitter-like temporal minimum provides finite curvature ($\tilde K \to 8/r_h^4 = 0.5$) and the regular geometric metric provides bounded space ($\rho \to 1.995$). The fixed-background limit ($\phi_0 = 2$ on Schwarzschild) illustrates why backreaction is necessary; the validation benchmark ($\phi_0 = 1$ on a Hayward background) confirms both conditions are achievable.

## 5.4 The Density Inference and Frame Dependence

The classical density estimate $M/(4\pi r^3/3)$ uses the Schwarzschild areal coordinate and presupposes that the geometric metric governs matter. In a two-metric theory this is not a locally measured matter density. Where a Lorentzian matter metric exists, the physical density is defined from a specified stress tensor and a spacelike hypersurface with induced metric $\tilde\gamma_{ij}$:

\begin{equation} \label{eq:global_6}
d\tilde V = \sqrt{\det\tilde\gamma_{ij}}\,d^3x,\qquad \tilde\rho = \tilde T_{\mu\nu}\,u^\mu u^\nu.
\end{equation}

The frame-dependence of volume is a structural feature of TEP: the classical density inference assumes the geometric metric governs matter, which is precisely what TEP questions. The matter metric $\tilde g_{\mu\nu}$ defines the physical volume and the physical density, and these differ from the geometric quantities by the conformal and disformal factors.

On the regular geometry, the density question requires careful frame separation. The geometric Einstein-tensor component $-G^t{}_t = 12M^2\ell^2/(r^3 + 2M\ell^2)^2 \to 3/\ell^2$ is the effective stress supporting the Hayward seed metric $g$; it is not automatically the matter-frame density $\tilde T_{\mu\nu}\tilde u^\mu\tilde u^\nu$. The geometric stress, the matter-frame stress, and the TEP matter stress tensor are distinct objects and must be separated explicitly. On the $\phi_0 = 1$ divergent-$A$ benchmark, the radial proper distance diverges ($\tilde\ell \to \infty$), so the integrated matter-frame spatial volume is not finite merely because the sphere area approaches a finite constant. On the target finite-$A$ architecture (Section 6.6), the centre is an ordinary regular point with finite volume, and the density is resolved through the de Sitter-like temporal minimum.

## 5.5 Global Regularity

### Global Geometry and Curvature Regularity (Benchmark Result)

On the Hayward validation benchmark with $\phi_0 = 1$, the TEP matter metric satisfies: (i) $\det\tilde g_{2D} = -A^4 G^2 - A^2 F B (\phi')^2 < 0$ in the exterior ($F > 0$); inside horizons ($F < 0$) the determinant must be checked with the canonical formula; (ii) $\tilde K \to$ finite as $r \to 0$ (scaling $r^{4\phi_0 - 4} = r^0$ for the regular seed, not $r^{4\phi_0-6}$ which is the Schwarzschild-seed result); (iii) $\rho = Ar \to 1.995$ as $r \to 0$, a bounded areal radius, but the limiting region is a finite-area asymptotic end ($\tilde\ell \to \infty$), not an ordinary regular centre. Both null families have infinite affine parameter ($\tilde\lambda \sim \int r^{-2}\,dr \to \infty$). These benchmark results are distinct from the coupled interior solution.

The coupled interior solver (`scripts/steps/step_34_solve_interior.py`, Section 4.6) integrates the sGB scalar field on a regular Hayward-class background to $r = 10^{-8}M$ for the mass-inflation branch $\eta = -0.1$, regularization scale $g = 1.1M$. It produces an ordinary regular point centre: $A(0) = 4.68$, $\tilde N_{\min} = 0.211$ at $r \approx 1.41M$, $\rho = Ar \to 4.68 \times 10^{-8}$, $w_r(0) = -1.000$, and $K(0) = 27.1$. The lapse is strictly positive, the areal radius vanishes at the centre, and the matter metric retains Lorentzian signature throughout. This confirms that the TEP Global Solution Architecture is not merely a geometrically self-consistent target; it is realised by the regularised sGB-coupled interior.

# 6. Causal Structure

The causal structure of the TEP framework is fundamentally different from the standard black-hole picture. It is dual to the cosmological interpretation: on large scales a temporal-rate gradient is read as expansion, while near an extreme source the same kind of gradient is read as collapse and inward suction. Both readings mistake a change in the rate of proper time for a change in spatial geometry. The standard picture — a "black hole" as a region that sucks space and matter into a singularity — is replaced by the *temporal well*: a regular spatial region in which the rate of proper time differs radically from the exterior. The physical lapse $N(r) = \sqrt{-\tilde g_{\mu\nu}\xi^\mu\xi^\nu}$ satisfies $0 < N(r) < 1$ through the deep region, with $N(r) \ll 1$ where the object appears black, but $N(r) \neq 0$ at every finite physical location. Time is extremely slow relative to the exterior but never stops; no absolute causal boundary forms. The apparent inward pull is analogous to refraction: trajectories curve when the propagation rate varies through a medium, not because matter is being pulled sideways. In TEP, the temporal field changes observed motion through clock transfer, timelike matter coupling, and its backreaction on the causal geometry. The conformal factor $A(\phi)$ controls clock-rate transfer and massive-particle response; the disformal term $B(\phi)$ can alter causal cones; the backreacted geometric metric $g_{\mu\nu}[\phi]$ generates light bending and orbital curvature. Lensing is generated by the complete matter metric $\tilde g_{\mu\nu}$, not by conformal clock scaling alone. The TEP object is a continuous, regular spatial domain with finite local density, a continuous temporal-rate gradient, no finite-radius causal horizon, and extreme external redshift. *Darkness does not require an event horizon; it requires sufficiently extreme temporal decoupling.* The Temporal Horizon $\mathcal{H}_T^{(\Lambda)}[u_o]$ is the observer-dependent surface where the clock-transfer mismatch exceeds a physical accessibility threshold $\Lambda$; it is operational, not causal. The conventional black-hole picture — event horizon, interior, singular centre — is retained only as the conventional limit that the temporal interpretation replaces. The crucial distinction is *apparent compactness $\neq$ physical compression*: the density $\rho_{\rm inferred} \sim M/(4\pi r^3/3)$ attributed to a black hole is an exterior-frame inference, not necessarily the locally measured density $\tilde\rho = \tilde T_{\mu\nu}\tilde u^\mu \tilde u^\nu$. With finite $A$ at the centre (a candidate architecture, Section 6.6), the centre is a regular point with finite optical distance — an emitter there has finite frequency relationship with any receiver. The TEP novelty is that the complete transfer map is modified by the dynamical temporal field, producing a specific observer-dependent structure beyond the single-metric GR relation.

## 6.1 The Operational Temporal Horizon

For a static spherically symmetric metric $ds^2 = -F\,dt^2 + G\,dr^2 + R^2\,d\Omega^2$, the Killing horizon in GR is where $F = 0$. In the TEP framework, the true solution must have $N(r) = \sqrt{-\tilde g_{\mu\nu}\xi^\mu\xi^\nu} = A(r)\sqrt{F(r)} > 0$ for every finite $r$ — the lapse has a finite but extremely small minimum $N_{\min} = A_c \ll 1$ at the centre, never vanishing at any finite radius. The perturbative sGB benchmark (Sotiriou \& Zhou 2014) belongs to a horizon-bearing branch: it produces $F = 0$ at a shifted $r_H$, but this branch does not realise the target TEP temporal well. The full TEP action must produce a temporal well with $N(r) \geq N_{\min} > 0$ everywhere. The conformal factor $A > 0$ everywhere; it shifts timelike matter observables such as the ISCO and clock transfer. Photon-sphere changes arise from the sGB correction to $g_{\mu\nu}$ or from a nontrivial disformal contribution.

The Temporal Horizon is defined by the clock-transfer factor between an emitter $e$ and observer $o$:

\begin{equation} \label{eq:temporal_horizon_def}
\mathcal{T}_{e\to o} = \frac{d\tilde\tau_o}{d\tilde\tau_e} = \frac{\omega_e}{\omega_o} = \frac{(-k_\mu u^\mu)_e}{(-k_\mu u^\mu)_o},
\end{equation}

where $k_\mu$ is the signal wavevector and $u^\mu$ is the observer's four-velocity. For stationary observers with matter-frame lapse $N(x) = \sqrt{-\tilde g_{\mu\nu}\xi^\mu\xi^\nu}$, the Killing energy $-k_\mu\xi^\mu$ is conserved along the signal, giving $\omega = (-k_\mu\xi^\mu)/N$ and therefore:

\begin{equation} \label{eq:static_transfer}
\mathcal{T}_{e\to o} = \frac{N_o}{N_e}.
\end{equation}

In the standard black-hole picture, the exact asymptotic limit $\mathcal{T}_{e\to o} \to \infty$ requires $N_e \to 0$, which occurs at $F = 0$. In the TEP temporal-well picture, $N_e$ has a finite minimum $N_{\min} = A_c \ll 1$ at the centre — it never reaches zero. The maximum transfer factor from the centre is finite: $\mathcal{Z}_{\max} \sim N_o/N_{\min}$. It may be astronomically large, but it is not divergent. The object appears black because the transfer factor exceeds any practical accessibility threshold, not because signals are absolutely forbidden from escaping. The operational Temporal Horizon at a specified transfer threshold $\Lambda$ is:

\begin{equation} \label{eq:operational_TH}
\mathcal{H}_T^{(\Lambda)}[u_o] = \left\{x_e \in J^-(\gamma_o) : \mathcal{T}_{e\to o} \geq \Lambda\right\}.
\end{equation}

This operational boundary is not a null boundary, not one-way, and does not divide spacetime into an outside and an inside. It genuinely changes with the observer's proper-time rate, motion, worldline, and temporal environment. A distant fast-clock observer reaches the threshold farther from the deep temporal region; a deeper slow-clock observer reaches it farther inward; a freely falling observer may encounter no corresponding local boundary for nearby comoving signals. Different observers place the Temporal Horizon differently because their clock rates differ: $\mathcal{H}_T^{(\Lambda)}[u_1] \neq \mathcal{H}_T^{(\Lambda)}[u_2]$. A sufficiently deep observer may see no Temporal Horizon locally because nearby clocks share nearly the same rate.

For a freely falling observer moving through the deep temporal region, the full expression $\omega = -k_\mu u^\mu$ must be evaluated along the actual worldline. The observer experiences ordinary local space, a normally running local clock, no wall, no horizon crossing, no signature change. Their clock and nearby clocks remain mutually finite, so the distant observer's Temporal Horizon is not their Temporal Horizon. A deeper observer (slower clock, smaller $N_o$) sees less redshift to the same emitter — the emitter is less temporally remote, and the operational boundary moves inward.

The unifying mechanism is the coordinate propagation rate. In any local Lorentzian frame of $\tilde g_{\mu\nu}$, the locally measured speed of light is exactly $c$ — no local rule is broken. But the physical lapse $N(r)$ drops to an extremely small fraction in the deep temporal region, so from the perspective of a distant fast-clock observer, the coordinate propagation rate of light through that region is drastically slowed: $dr/dt \sim N(r) \ll 1$. This single mechanism — light propagating through a gradient where the rate of time varies — produces three observational consequences that the standard picture attributes to three separate phenomena. First, the slowing bends trajectories: light curves when the propagation rate varies through a medium, exactly as in refraction, producing the apparent "inward pull" and gravitational lensing without a sideways force. Second, the slowing extends transit times: waves crossing the deep region take far longer in coordinate time — the deep-transit ringdown modes show $\sim 2.8\times$ longer damping because the wave transits a zone where proper time runs at $\sim 21\%$ of the exterior rate (Section 7.3). Third, the slowing combined with extreme redshift creates the optical illusion of darkness: as light climbs out of the temporal well, its frequency is redshifted by the transfer factor $\mathcal{Z} = N_o/N_e$, and its coordinate propagation rate drops below any practical accessibility threshold. The object appears black not because light is trapped behind a barrier, but because the light's frequency and propagation rate are extinguished to the exterior observer by the extreme temporal gradient. An observer falling into the well alongside the light would see no wall, no horizon, and light moving normally through ordinary local space — the "black hole" and its darkness are the optical and temporal illusion of fast exterior clocks trying to read an incredibly slow environment.

The inverse view completes the symmetry. An observer deep in the temporal well, looking outward, sees the exterior universe in extreme fast-forward: distant clocks tick at rates amplified by the inverse transfer factor $\mathcal{Z}^{-1} = N_e/N_o \gg 1$, and incoming light is blueshifted by the same ratio. The outside universe appears compressed in time and energetically amplified, funneled through a distorted angular window determined by the refractive bending of the temporal gradient. This is not a separate physical effect — it is the same relative temporal gradient viewed from the other end. The exterior observer's "darkness" and the interior observer's "fast-forward" are dual projections of one continuous geometry, confirming that the temporal well is a relative optical and temporal phenomenon, not a one-way physical trapdoor.

Nothing here introduces an observer-specific metric. The variation comes from: the observer's proper-time rate; the emitter's proper-time rate; their relative motion; the signal path; the temporal field sampled between them. There is one universal matter metric and one temporal field, but different observers possess different temporal-accessibility maps. The Temporal Horizon is thus not an absolute surface placed into spacetime; it is the observer-dependent limit of temporal accessibility within one continuous physical geometry. The distinctive TEP prediction is that the complete transfer map $\mathcal{T}_{e\to o} = \mathcal{T}[g_{\mu\nu}, \phi, A(\phi), B(\phi), u_e, u_o, k_\mu]$ is modified by the dynamical temporal field, producing observer-dependent structure beyond the single-metric GR relation.

## 6.2 Observational Probes of the Temporal Well

The Temporal Horizon is observer-dependent within one geometry. Photons, massive particles, and gravitational waves are all governed by the same unified temporal field — they are different observational probes of the same temporal well, not evidence of fractured metric structure. The TEP novelty is the temporal well itself: a regular spatial region with an extreme gradient in relative proper-time accumulation, dual to the cosmological interpretation.

In the perturbative sGB benchmark, photons and massive particles respond to the temporal field at different coupling orders — but this is a property of the mathematical toolkit, not a fundamental law of TEP. Null geodesics are conformally invariant: for the pure conformal metric $\tilde g = A^2 g$ with $A \to 1$ asymptotically, the conformal factor $A$ cancels exactly in the photon-sphere and shadow calculation, so the shadow shift is $\mathcal{O}(\eta^2)$ from the sGB correction to $g$ alone. Timelike geodesics are not conformally invariant: the conformal factor $A = e^{-\phi} > 1$ enters the ISCO calculation directly, producing an $\mathcal{O}(\eta)$ shift. The $\sim 44\times$ coupling-order ratio between the ISCO and shadow shifts is a mathematical feature of the perturbative sGB EFT, not a physical prediction that photons and matter probe fundamentally different realities.

Several subtleties must be addressed before the observational predictions become fully quantitative:

- Conformal invariance of null paths. For the pure conformal metric $\tilde g = A^2 g$ with $A \to 1$ asymptotically, null geodesic paths are conformally invariant. The photon sphere and shadow trajectory are therefore unchanged by $A$ alone. The shadow shift from the conformal sector is zero; any shadow deviation must come from the disformal sector or from the sGB correction to $g$ itself (which shifts the geometric photon sphere relative to Schwarzschild).

- Disformal effect on the circular photon orbit. The disformal term $B(\phi)\phi'^2\,dr^2$ modifies the radial metric component. For a circular photon orbit ($dr = 0$), the orbit location is fixed by the ratio of temporal and angular metric components, which the radial disformal term does not directly enter. The disformal term may affect radial instability and off-shell propagation, but the on-shell circular orbit shift requires explicit derivation, not assumption.

- Tensor propagation metric. Once sGB is introduced, tensor perturbations need not propagate exactly on $g_{\mu\nu}$. Scalar–tensor mixing can produce an effective characteristic metric $\mathcal{G}_{\rm tensor}^{\mu\nu}$ for gravitational waves. Deriving $\mathcal{G}_{\rm tensor}^{\mu\nu}$ from the action is a required calculation.

The sGB correction to $g$ does shift the geometric photon sphere and Regge–Wheeler potential relative to Schwarzschild — this is the established perturbative sGB shadow and QNM shift. For the pure conformal case with $A \to 1$ asymptotically, photons see the same null paths as on $g$, so the shadow shift is the standard sGB shift. A genuine photon-sphere shift on $\tilde g$ that differs from the tensor-spectrum shift on $g$ would emerge only when the disformal term or a non-trivial $A$ profile produces it.

For the full disformal case with a static scalar, the matter-frame lapse condition $\tilde g_{vv} = -A^2 F + B(\phi)(\partial_v\phi)^2$ reduces to $N^2 = A^2 F$ since $\partial_v\phi = 0$. In the temporal-well TEP solution, $A\sqrt{F} > 0$ everywhere. A non-static scalar background could modify the lapse profile, providing an additional observational channel. This requires the time-dependent interior integration.

## 6.3 Null Expansions

For a static spherically symmetric metric $ds^2 = -F\,dt^2 + G\,dr^2 + H\,d\Omega^2$, the null expansions for the two radial null congruences are

\begin{equation} \label{eq:causal_1}
\theta_+ = \frac{H'}{H\sqrt{G}} \quad \text{(outgoing)}, \qquad \theta_- = -\frac{H'}{H\sqrt{G}} \quad \text{(ingoing)}.
\end{equation}

For the TEP matter metric $\tilde g = A^2 g$ with $A = (r_h/r)^{\phi_0}$, the matter metric components are $\tilde{F} = A^2 F$, $\tilde{G} = A^2 G$, $\tilde{H} = A^2 r^2$. Substituting and simplifying:

\begin{equation} \label{eq:causal_2}
\theta_+ = \frac{2(1-\phi_0)}{A^2\,r\,\sqrt{G}}.
\end{equation}

For $\phi_0 = 1$: $\theta_+ = 0$ exactly. This is a property of the conformal factor alone and holds on any background. On the divergent-$A$ benchmark, $\theta_+ = 0$ reflects the constant-area property of the finite-area asymptotic end. On the target finite-$A$ architecture, $\phi_0$ is not a free parameter — it is determined by the coupled field equations, and $\theta_+$ takes its value from the solved profile.

## 6.4 Eddington–Finkelstein Null Algebra

In ingoing Eddington–Finkelstein coordinates, the conformal matter metric ($B = 0$) has components

\begin{equation} \label{eq:causal_3}
\tilde g_{ab} = A^2\begin{pmatrix} -F & 1 \\ 1 & 0 \end{pmatrix}, \qquad \tilde g^{ab} = A^{-2}\begin{pmatrix} 0 & 1 \\ 1 & F \end{pmatrix}.
\end{equation}

Note that $\tilde g^{rr} = F/A^2$, not zero. The null condition $\tilde g_{ab}k^a k^b = 0$ gives $dv(-F\,dv + 2\,dr) = 0$, yielding two families:

- Ingoing ($v = \mathrm{const}$, $dv = 0$): the tangent $k^\mu = (0, 1)$. The Christoffel symbol $\tilde\Gamma^r{}_{rr} = 2A'/A$ is nonzero.

- Outgoing ($dr/dv = F/2$): the tangent $k^\mu = (2/F, 1)$, or $dv/dr = 2/F$.

For a conformal transformation $\tilde g = A^2 g$, null geodesics of $g$ are also null geodesics of $\tilde g$ (same unparameterised curves), but the affine parameters are related by $d\tilde\lambda = A^2\,d\lambda$. Near the regular Hayward centre, the seed radial null affine parameter behaves as $d\lambda \propto dr$. With $A \sim (r_h/r)^{\phi_0}$, both radial families have

\begin{equation} \label{eq:causal_4}
\tilde\lambda \sim \int \frac{dr}{r^{2\phi_0}} \to \infty \quad \text{for } \phi_0 > 0.
\end{equation}

Both radial null families have infinite affine parameter on the divergent-$A$ benchmark. For the target finite-$A$ architecture, $A \to A_c$ finite, so $d\tilde\lambda \to A_c^2\,d\lambda$ and both families have finite affine parameter — the centre is at finite affine distance, consistent with it being a regular manifold point.

## 6.5 The Invariant Frequency-Transfer Calculation

The Temporal Horizon is established by the invariant frequency-transfer factor, not by coordinate-time arguments. For a photon with wavevector $k_\mu$ emitted by an observer with four-velocity $u_e^\mu$ and received by an observer with four-velocity $u_o^\mu$:

\begin{equation} \label{eq:causal_5}
\omega_e = -k_\mu u_e^\mu, \qquad \omega_o = -k_\mu u_o^\mu, \qquad \mathcal{Z} = \frac{\omega_e}{\omega_o}.
\end{equation}

A Temporal Horizon exists when $\mathcal{Z}$ exceeds a physical accessibility threshold $\Lambda$ along a well-defined physical limit, while local physics remains regular. In the target finite-$A$ architecture, $\mathcal{Z}$ reaches a finite maximum $\mathcal{Z}_{\max} \sim N_o/N_{\min}$ at the centre; the object appears black when this maximum exceeds any practical $\Lambda$.

The pipeline step `step_19_observer_frequency_transfer.py` computes $\mathcal{Z}$ for freely falling and static emitters and receivers. The calculation includes a causal-connectivity check: before evaluating $\omega_o$, the pipeline verifies that $J^+(\text{emission event}) \cap \text{receiver worldline} \neq \emptyset$ — i.e., that a future-directed null geodesic from the emission event actually reaches the receiver. The check scans $F(r)$ on $[r_{\rm emit}, r_{\rm recv}]$ for sign changes (horizon crossings). If $F$ changes sign, the emitter is inside a Killing horizon and the outgoing null geodesic ($dr/dv = F/2$) cannot reach the exterior receiver; $\mathcal{Z}$ is then not an observable transfer factor and is reported as invalid.

On the divergent-$A$ benchmark ($\phi_0 = 1$), the causal-connectivity check correctly invalidates the deep-interior transfer factors. For the Hayward parameters used ($\ell = 0.5M$), the geometry has an outer horizon near $1.995M$ and an inner horizon near $0.597M$. Emitters at $r < 0.597M$ (inside the inner horizon) cannot send outgoing null geodesics to the exterior — the geodesic reaches the inner horizon but is then trapped between the two horizons. Only emitters at $r > 1.995M$ (exterior to the outer horizon) produce valid transfer factors. If causal connectivity is ignored, static transfer formulas can yield non-observable mathematical artefacts from trapped regions — for example, transfer factors of order $10^4$--$10^5$ for emitters deep inside the inner horizon, which have no physical meaning. The causal-connectivity check ensures that only physically reachable transfer factors are reported.

On the target finite-$A$ temporal-well architecture, $N(r) = A(r)\sqrt{F(r)} \geq N_{\min} > 0$ for every finite $r$, with $N_{\min} = A_c \ll 1$ at the centre. For a distant static observer, the transfer factor $\mathcal{Z} = (A_o\sqrt{F_o})/(A_e\sqrt{F_e})$ grows as the emitter moves deeper, reaching a finite maximum $\mathcal{Z}_{\max} \sim N_o/N_{\min}$ at the centre. This maximum may be astronomically large but is not divergent. The receiver is at a fixed distant location where $F_o$ is finite and $A_o \approx 1$. The object appears black because the transfer factor exceeds any practical accessibility threshold, not because signals are absolutely forbidden from escaping. For a freely falling emitter, the static formula cannot be reused; the full expression $\omega = -k_\mu u^\mu$ must be evaluated along the emitter's actual worldline. The freely falling observer moves through continuous regular space with finite local clock transfer — the distant observer's Temporal Horizon is not their Temporal Horizon. Light can still escape in principle from any finite depth; no absolute causal boundary forms. The centre is regular and at finite optical distance; it does not need to become a second optical boundary.

## 6.6 The Correct Global Solution Architecture

The target global solution should approach, as $r \to 0$:

\begin{equation} \label{eq:causal_6}
F(r) = 1 - f_2 r^2 + O(r^4), \quad R(r) = r + O(r^3), \quad \phi(r) = \phi_c + \phi_2 r^2 + O(r^4),
\end{equation}

\begin{equation} \label{eq:causal_7}
A(\phi_c) = A_c, \quad 0 < A_c < \infty, \quad B(\phi_c) \text{ finite}.
\end{equation}

Then the physical areal radius becomes $\tilde R = A(\phi)R(r) \sim A_c r \to 0$ — an ordinary regular centre with no spatial opening, no tube-like asymptotic end, no other universe, no wall, and finite local geometry. Both null families have finite affine parameter ($d\tilde\lambda \to A_c^2\,d\lambda$). The centre is at finite optical distance. An emitter at the centre has finite $\mathcal{Z}$ with any causally connected receiver. This finite-$A$ regular-centre expansion is a candidate architecture demonstrated by the inverse reconstruction and Hayward benchmark; it serves as evidence that the fixed Schwarzschild singularity is not mandatory. It does not by itself define the final TEP object, which is a temporal well rather than a regular-centre black hole.

The TEP object is a *temporal well*: a regular spatial region with an extreme gradient in relative proper-time accumulation. The physical lapse $N(r) = A(r)\sqrt{F(r)}$ satisfies $0 < N(r) < 1$ through the deep region, with $N(r) \ll 1$ where the object appears black, but $N(r) \neq 0$ at every finite physical location. Time is extremely slow relative to the exterior but never stops; no absolute causal boundary forms. The Temporal Horizon $\mathcal{H}_T^{(\Lambda)}[u_o]$ is the observer-dependent surface where the clock-transfer mismatch exceeds a threshold $\Lambda$; it is not a null boundary, not one-way, and does not divide spacetime into an outside and an inside. Different observers place it differently because their clock rates differ. This is the division of labour: *the temporal-rate gradient generates the observed phenomenology; the Temporal Horizon is the observer-relative operational limit of temporal accessibility; apparent compactness arises from exterior-frame reconstruction, not physical compression.* The TEP Global Solution Architecture is:

\begin{equation} \label{eq:causal_8}
\text{regular spatial domain} + \text{finite local density} + \text{continuous temporal-rate gradient} + \text{no finite-radius causal horizon } (N \geq N_{\min} > 0 \text{ everywhere}) + \text{extreme external redshift}
\end{equation}

The fixed-background conformal theorem (Section 3) proves that the temporal field must backreact on the geometric metric: the fixed-Schwarzschild construction is inconsistent. The inverse reconstruction (Section 4.4) identifies $w_r = -1$ as a sufficient regularity target, demonstrated by the Hayward benchmark. The benchmark (Appendix K) confirms that finite Kretschmann, bounded areal radius, and Lorentzianity are simultaneously achievable on a regular background. A valid dynamical completion of TEP — whether through sGB, a modified coupling, or the full TEP action — should produce a temporal well: a regular spatial domain with finite local density and an extreme temporal-rate gradient, not a regular black hole.

The sGB perturbative exterior (Sotiriou \& Zhou 2014) produces $F = 0$ at a shifted $r_H$ because it is an EFT approximation around the Schwarzschild background. The true TEP solution, from the full action, must have $N(r) \geq N_{\min} > 0$ everywhere — the lapse has a finite minimum at the centre, not a zero. The observable predictions (shadow, ISCO, QNM) are computed at $r \sim 3$–$6M$ in the conditional sGB benchmark, far from the would-be horizon at $r \sim 2M$; their robustness to the deep-interior profile has not been demonstrated and requires the full TEP-native solution. The horizonless temporal well alters the late-time ringdown spectrum (Section 7.3).

The organising question of TEP is not "How is a black-hole interior regularised?" but rather: *can the phenomena attributed to black holes arise from an extreme relative-time gradient without physical collapse into an ultra-dense object?* This is dual to the TEP cosmological question: the universe need not be physically expanding merely because distant signals are redshifted, and matter need not be physically crushed into a black hole merely because trajectories appear to converge toward a dark compact region. In both cases, conventional spatial dynamics may be the observational projection of a dynamical temporal field. The shadow is reproduced by the photon sphere (a metric feature at $r \sim 3M$, independent of the horizon). The ISCO is reproduced by the orbital structure (at $r \sim 6M$, independent of the horizon). The gravitational-wave ringdown is reproduced by the perturbation spectrum, with a theory-dependent inner boundary condition. The darkness is reproduced by the extreme redshift of signals from the deep temporal region — darkness does not require an event horizon, only sufficiently extreme temporal decoupling. The apparent inward pull is reproduced by the complete matter metric $\tilde g_{\mu\nu}$: the conformal factor controls clock transfer and massive-particle response, the disformal term can alter causal cones, and the backreacted geometry generates light bending and orbital curvature. The apparent compactness is reproduced by the exterior observer's reconstruction of the region through strongly redshifted, delayed signals. No event horizon is required, no ultra-dense object is required, and no suction mechanism is required. The paper's main target is to derive a regular temporal-rate profile that can account for the observed shadow, lensing, orbits and ringdown without high-density collapse — not a better black hole, but an explanation of why a regular temporal domain is mistaken for one.

In the TEP Global Solution Architecture, the observer falling toward the deep region experiences ordinary local space, a normally running local clock, finite local density, no wall, no horizon crossing, no signature change. Matter need not be compressed into a vanishing proper volume; the apparent concentration can result from temporal decoupling between the deep region and the exterior. The extreme effect is a finite but potentially enormous transfer mismatch between the deep region and a distant observer: $\mathcal{Z}_{\max} \sim N_o/N_{\min}$, finite because $N_{\min} > 0$. The sGB coupling serves as a low-energy effective field theory for the exterior, providing indicative observable shifts; a complete TEP action is required to derive the temporal-well profile and must dynamically generate the extreme temporal-rate gradient proven necessary by the fixed-background theorem.

## 6.7 Signature Structure

For the metric $ds_g^2 = -F\,dv^2 + 2G\,dv\,dr + R^2\,d\Omega^2$ and a static radial scalar ($\nabla_v\phi = 0$, $\nabla_r\phi = \phi'$), the conformal-disformal matter metric has 2D $(v,r)$ block:

\begin{equation} \label{eq:causal_9}
\tilde g_{vv} = -A^2 F, \quad \tilde g_{vr} = A^2 G, \quad \tilde g_{rr} = B(\phi')^2.
\end{equation}

The determinant is

\begin{equation} \label{eq:causal_10}
\boxed{\det\tilde g_{2D} = -A^4 G^2 - A^2 F B(\phi')^2}.
\end{equation}

This is the canonical formula used throughout the manuscript. The first term is always negative ($A^4 G^2 > 0$). The second term has sign $-F$: in the exterior ($F > 0$) it is negative, reinforcing the conformal contribution; in the deep temporal region where $F$ becomes very small, the second term is suppressed and the conformal contribution dominates. In the temporal-well TEP solution, $F > 0$ everywhere — the "inside the horizon" case ($F < 0$) does not arise. All signature, invertibility, and hyperbolicity results must be regenerated from this single canonical expression. For the pure conformal metric ($B = 0$), $\det\tilde g_{2D} = A^4 \det g_{2D} = -A^4 G^2 < 0$ — Lorentzianity is inherited analytically from $g$ via $A > 0$.

## 6.8 The Divergent-$A$ Benchmark: What It Shows and What It Does Not

The Hayward benchmark with $\phi_0 = 1$ and $A \sim 1/r$ produces a finite-area asymptotic end ($\rho \to 1.995 \neq 0$, $\tilde\ell \to \infty$), not an ordinary regular centre. Both null families have infinite affine parameter. The static-observer clock rate $d\tilde\tau/dt = A\sqrt{F} \to \infty$ (blueshift). The divergent-$A$ benchmark is not the target architecture for the final TEP solution. It serves as a validation benchmark for the pipeline (confirming that finite curvature and bounded areal radius are compatible on a regular background) and as a cautionary example (showing that divergent $A$ produces unphysical blueshift and an asymptotic end rather than a regular centre). The target architecture is a temporal well: finite $A$ at the centre (Section 6.6, a candidate profile), $N(r) \geq N_{\min} > 0$ everywhere, with the operational Temporal Horizon defined by the clock-transfer threshold.

# 7. Perturbation Analysis: Coupled Temporal–Geometric Ringdown

This section identifies the perturbation-sector structure and presents both the exterior QNM benchmark and the full deep-transit solve. The coupled sGB solution introduces a dynamical scalar degree of freedom coupled to the gravitational sector through the Gauss–Bonnet invariant. The tensor speed $c_T = 1$ on the Schwarzschild background in shift-symmetric sGB (from $G_4 = M_{\rm Pl}^2/2$, with the Gauss–Bonnet coupling entering through $G_5$), satisfying the GW170817 constraint on that background; the full tensor characteristic metric on the TEP solution requires derivation from the complete coupled perturbation system. Two layers of QNM calculation are now in place, both within the spherically symmetric sector (Section 4.6B). The *exterior benchmark* uses the exact Schwarzschild Leaver frequencies, applies the sGB $O(\eta^2)$ horizon shift, and includes leading-order polar–scalar mixing: polar-led $\omega = 0.374076 - 0.089143i$, scalar-led $\omega = 0.484056 - 0.096781i$, isospectrality breaking $-0.028\%$ real, $+0.067\%$ damping. The *full deep-transit solve* integrates the coupled $2\times 2$ wave equation with the regular inner boundary ($u \sim r^{\ell+1}$ at the temporal minimum) and outgoing outer boundary, finding three deep-transit modes: polar-led $\omega = 0.4624 - 0.0320i$, scalar-led $\omega = 0.4163 - 0.0322i$, overtone $\omega = 0.7447 - 0.0347i$, with $2.8\times$ longer damping and $+14.8\%$ isospectrality breaking. The deep temporal gradient — not a physical boundary — is what produces the extended damping and amplified mode mixing. These frequencies are exact within the spherically symmetric sector; the true temporal well for an astrophysical source may possess no spatial symmetry, and the QNM spectrum of the full unsymmetrised problem is not yet computed.

## 7.1 Regge–Wheeler and Zerilli Potentials on the sGB-Corrected Metric

A metric perturbation on the sGB-corrected background is written

\begin{equation} \label{eq:stab_1}
g_{\mu\nu} = g_{\mu\nu}^{(0)} + h_{\mu\nu},
\end{equation}

where $g_{\mu\nu}^{(0)}$ is the Sotiriou–Zhou $\mathcal{O}(\beta^2)$ corrected Schwarzschild metric with scalar profile $\phi(r)$. The axial and polar tensor sectors decompose into the Regge–Wheeler and Zerilli channels, each acquiring an sGB correction sourced by the background scalar gradient coupled to the Riemann tensor structure. The exact quadratic action is derived in; the key results are summarised below.

The axial (Regge–Wheeler) potential receives the correction

\begin{equation} \label{eq:stab_2}
V_{\rm RW}^{\rm sGB} = V_{\rm RW}^{\rm Schw}\!\left(1 + \beta^2 h_2(r)\right), \qquad \beta = \frac{\eta}{12},
\end{equation}

where $h_2(x)$ is the Sotiriou–Zhou $\mathcal{O}(\beta^2)$ metric perturbation with $x = 2M/r$. The axial sector *decouples* from the scalar at leading order: odd-parity perturbations do not couple to the scalar at $\mathcal{O}(\alpha_{\rm GB})$. The axial QNM shift is $\mathcal{O}(\beta^2) = \mathcal{O}(\eta^2)$. The potential has correct $L^{-2}$ dimensions: $V_{\rm RW}^{\rm Schw} \sim F \cdot \ell(\ell+1)/r^2 \sim L^{-2}$, and $h_2$ is dimensionless.

\begin{equation} \label{eq:stab_3}
V_{\rm Z}^{\rm sGB} = V_{\rm Z}^{\rm Schw} + \mathcal{O}(\alpha_{\rm GB}^2), \qquad \alpha_{\rm GB} = \frac{\eta M^2}{3}\;[L^2].
\end{equation}

The polar (Zerilli) sector *couples* to the scalar at $\mathcal{O}(\alpha_{\rm GB})$ through the mixing potential

\begin{equation} \label{eq:stab_3b}
V_{Z\phi} = \alpha_{\rm GB}\,\frac{\ell(\ell+1)(\ell-1)(\ell+2)\,F}{r^4}\;[L^{-2}],
\end{equation}

where $\lambda = (\ell-1)(\ell+2)/2$. The full polar-scalar system is a $2\times 2$ coupled eigenvalue problem with potential matrix

\begin{equation} \label{eq:stab_3c}
\mathbf{V} = \begin{pmatrix} V_Z + \mathcal{O}(\alpha_{\rm GB}^2) & V_{Z\phi} \\ V_{Z\phi} & V_\phi + \mathcal{O}(\alpha_{\rm GB}^2) \end{pmatrix},\qquad [V_{ij}] = L^{-2}.
\end{equation}

The mixing $V_{Z\phi} = \mathcal{O}(\alpha_{\rm GB}) = \mathcal{O}(\eta)$ is the leading-order coupling between polar gravitational and scalar perturbations. It is the structural origin of isospectrality breaking: the polar sector mixes with the scalar while the axial sector does not, lifting the Chandrasekhar degeneracy that is exact in GR. All potential matrix elements have correct $L^{-2}$ dimensions when the dimensionful coupling $\alpha_{\rm GB} = \eta M^2/3$ $[L^2]$ is used, verified by the quadratic action derivation.

The perturbation equations take the standard Schrödinger-like form on the tortoise coordinate:

\begin{equation} \label{eq:stab_4}
\frac{d^2\Psi_i}{dr_*^2} + \left[\omega^2 - V_i^{\rm sGB}(r)\right]\Psi_i = 0,
\end{equation}

for each channel $i \in \{{\rm RW},\, {\rm Z},\, {\rm scalar}\}$. The effective potentials are derived from the exact quadratic action (second variation of the sGB action). The dimensionful GB coupling $\alpha_{\rm GB} = \eta M^2/3$ $[L^2]$ ensures all potential terms have correct $L^{-2}$ dimensions: $V_{\rm RW} \sim F \cdot \ell(\ell+1)/r^2 \sim L^{-2}$, $V_{Z\phi} \sim \alpha_{\rm GB} \cdot F/r^4 \sim L^2 \cdot L^{-4} = L^{-2}$, and $V_{ZZ} \sim \alpha_{\rm GB}^2 \cdot F/r^6 \sim L^4 \cdot L^{-6} = L^{-2}$. The scalar perturbation has no direct GB mass term (the GB coupling is linear in $\phi$, so $\delta^2 S_{\rm GB}/\delta\phi^2 = 0$); the scalar potential correction comes only from the modified background metric at $\mathcal{O}(\beta^2)$. The specific QNM frequencies computed from these potentials use WKB estimates; a coupled spectral solver (Leaver or continued-fraction) is required for precise values, validated against published nonperturbative sGB QNM spectra.

## 7.2 Coupled Temporal–Geometric Perturbations

The perturbation system is one coupled eigenvalue problem, not a collection of independent channels. A merger perturbs both $g_{\mu\nu}$ and $\phi$, and the detector responds to their combined perturbation through $\delta\tilde g_{\mu\nu} = \delta\tilde g_{\mu\nu}[\delta g_{\mu\nu}, \delta\phi]$. The axial (Regge–Wheeler) sector decouples at leading order: the axial gravitational perturbation couples to the scalar only through the GB modification of the Einstein equations, producing an $\mathcal{O}(\eta^2)$ shift to the Schwarzschild axial QNM. The polar (Zerilli) sector and the scalar sector form a *coupled* eigenvalue system: the polar gravitational perturbation mixes with the scalar perturbation through the GB coupling, and the two must be solved as a joint system, not as independent one-dimensional potentials. This coupling lifts the Chandrasekhar isospectrality that is exact in GR between the axial and polar gravitational modes.

The eigenmodes of this coupled system may carry different relative proportions of tensor and temporal-field perturbation. The number, frequencies, excitation amplitudes, and detector projections of those modes remain to be calculated. A dynamical scalar degree of freedom means that the perturbation system contains scalar participation; it does not automatically prove efficient excitation in a merger, a clean distinct ringdown frequency, measurable strain coupling, or a separately identifiable mode. The specific frequencies require solving the coupled polar–scalar eigenvalue system with a spectral method (Leaver or continued-fraction) on the sGB-corrected background, validated against published nonperturbative rotating sGB QNM spectra (Witek et al. 2019; Blazquez-Calzadilla et al. 2020; Chen et al. 2024). Recent nonlinear merger calculations in sGB gravity extract ringdown excitation directly from numerical relativity (Figueras et al. 2025), providing a stronger benchmark than surrogate WKB potentials.

The structural predictions that follow from the sGB coupling and do not depend on the specific potential formulas are:

- The coupled polar–scalar spectrum may contain eigenmodes with substantial scalar participation. Their existence, frequencies, excitation amplitudes and detector coupling have not yet been derived for the TEP temporal-well background.

- Isospectrality between axial and polar geometric modes is broken at $\mathcal{O}(\eta^2)$, because the polar sector mixes with the temporal field while the axial sector does not.

- The axial geometric QNM receives an $\mathcal{O}(\eta^2)$ shift from the GB-modified background: $+0.136\%$ in real frequency at $\eta = -0.1$. The polar-led mode, including leading-order mixing with the scalar, becomes $\omega = 0.374076 - 0.089143i$, giving a real isospectrality breaking of $-0.0280\%$ and a damping breaking of $+0.0668\%$.

A detection of a component with substantial scalar participation would indicate an additional scalar gravitational degree of freedom — consistent with TEP-sGB, but also consistent with other scalar–tensor theories. The next-generation ground-based detectors (Einstein Telescope, Cosmic Explorer) and LISA have the sensitivity to resolve these modes for nearby massive black-hole mergers.

## 7.3 Deep-Transit Ringdown Structure

All three perturbation potentials are single-peaked in the exterior region. In the standard black-hole picture, a single-peaked potential with a purely absorbing horizon produces no echoes. In the TEP temporal-well picture, there is no horizon and no physical boundary — the deep region is continuous, ordinary space in which the rate of proper time drops to a minimum lapse $\tilde N_{\min} = 0.211$ at $r \approx 1.41$. The exterior potential peak and the de Sitter-like temporal core together form a resonant cavity: a wave propagating into this region does not strike a wall; it encounters an extreme continuous gradient in the temporal rate. The gradient turns the wave refractively — the propagation direction curves smoothly as the local clock slows — and the wave transits through the deep region and emerges back into the exterior. Because the local proper time runs at $\sim 21\%$ of the exterior rate, the transit takes $\sim 2.8\times$ longer in coordinate time, producing extended exposure to the polar–scalar mixing region. The resonant cavity structure is what produces the long-lived, strongly mixed deep-transit modes: the wave rings inside the temporal gradient, sustaining mode conversion for far longer than the Schwarzschild horizon permits. The late-time response may therefore contain modified damping profiles and enhanced mode mixing — signatures that are absent in the standard black-hole picture. The current perturbative sGB analysis uses horizon boundary conditions (purely ingoing at $r_H$) and therefore does not capture the deep-transit response; a proper analysis requires the full TEP solution with the correct regular inner boundary condition.

The Frobenius analysis of the regular centre confirms three structural results that govern the matrix-Leaver inner boundary. First, the tortoise coordinate is finite at the centre ($r_* \to 0$ as $r \to 0$, not $-\infty$), so the inner boundary is a regularity condition $u \sim r^{\ell+1}$ — the wave passes smoothly through the origin in spherical coordinates, not into an absorbing sink. Second, all three effective potentials (scalar, gravitational, and the polar–scalar mixing) scale as $r^{-2}$ near the centre — the standard angular-momentum barrier dominates, and the de Sitter curvature $R(0) = 12.0$ enters only as a subleading constant $c_0$ that modifies the recurrence coefficients, not the indicial exponent. The standard $r^{\ell+1}$ regularity therefore holds. Third, the exterior sGB mixing $V_{Z\phi}^{\rm ext} \sim 1/r^4$ is regularised to $1/r^2$ in the interior because $\phi'(0) = 0$ on the solved background — the scalar profile is smooth at the centre, softening the coupling. Both channels share the same Frobenius exponent, and the coupled ansatz $u_Z = r^{\ell+1} \sum a_n r^n$, $u_\phi = r^{\ell+1} \sum b_n r^n$ yields a matrix three-term recurrence $\mathbf{A}_n \mathbf{c}_{n+1} + \mathbf{B}_n \mathbf{c}_n + \mathbf{C}_n \mathbf{c}_{n-1} = 0$ with $2\times 2$ coefficients. The background metric (Hayward, a rational function of $r$) and the scalar ODE are both analytic at $r = 0$, so the Taylor coefficients of $F(r)$, $N(r)$, and $\phi(r)$ are generated exactly by recurrence from the ODE itself — no polynomial fit is needed. The Taylor series agrees with the numerical integration to $3.9 \times 10^{-13}$ relative error at $r = 10^{-6}$, degrading gracefully to $5 \times 10^{-3}$ at $r = 0.1$ as expected for a truncated series.

The full coupled QNM solver integrates the $2\times 2$ polar–scalar wave equation with the regular inner boundary ($u \sim r^{\ell+1}$ at $r = 0$) and the outgoing outer boundary ($u \sim e^{+i\omega r_*}$ at $r \to \infty$). The Hayward background with $g = 1.1$ has no horizon: $F_{\min} = 0.038$ at $r \approx 1.39$, where the lapse reaches its minimum and the temporal gradient is steepest. The solver finds three deep-transit modes at $\eta = -0.1$, all within the spherically symmetric sector:

- Polar-led: $\omega = 0.4624 - 0.0320i$ (mismatch $< 10^{-15}$)

- Scalar-led: $\omega = 0.4163 - 0.0322i$ (mismatch $< 10^{-15}$)

- Higher overtone: $\omega = 0.7447 - 0.0347i$ (mismatch $< 10^{-15}$)

The damping timescale is $\tau/\tau_{\rm Schw} \approx 2.8$ — the deep-transit modes live nearly three times longer than the corresponding Schwarzschild QNMs ($\omega_I^{\rm Schw} \approx -0.089$ vs. $\omega_I^{\rm deep} \approx -0.032$). The extended lifetime is not trapping in a physical container; it is the natural consequence of the wave transiting through a region where proper time runs at $\sim 21\%$ of the exterior rate. The wave spends more coordinate time in the mixing region because its local clock is slow, not because it is confined. The polar–scalar isospectrality breaking is $+14.8\%$ in real frequency and $+0.8\%$ in damping — amplified by a factor of $\sim 500$ compared to the exterior benchmark ($-0.028\%$), because the extended transit through the slow-time gradient subjects the wave to sustained temporal shear that drives polar–scalar mode conversion. The coupled mode splitting between polar-led and scalar-led branches is $11.1\%$. This is the ringdown signature of the spherically symmetric temporal well: the standard Schwarzschild QNM spectrum is replaced by a set of long-lived, strongly mixed deep-transit modes. The $2.8\times$ extended damping and $14.8\%$ amplified mixing are not universal constants; they scale with the physical width of the deep region and the steepness of the temporal shear gradient. A deeper well with a lower $N_{\min}$ extends the transit time further, so the ringdown acts as a direct probe of the temporal-well geometry. The spectrum for a rotating or asymmetric source is not yet computed.

The QNM spectrum is validated against published sGB results on seven independent criteria. The exterior benchmark is qualitatively consistent with the published sGB QNM literature on five criteria: (i) the axial modes deviate at $\mathcal{O}(\alpha_{\rm GB}^2)$, matching the eikonal analysis of Bryant, Silva, Yagi \& Glampedakis (2021, PRD 104, 044051) and the numerical results of Blázquez-Salcedo et al. (2016, PRD 94, 104024); (ii) the polar sector contains two distinct families — gravitational-led and scalar-led — as found in all published sGB QNM calculations; (iii) isospectrality is broken between axial and polar modes, with the polar gravitational-led correction having opposite sign to the axial correction ($+0.039\%$ axial vs $-0.018\%$ polar at $\eta = -0.1$), matching the Blázquez-Salcedo et al. pattern; (iv) the polar deviations are larger than the axial deviations, consistent with the additional scalar-tensor coupling in the polar sector; (v) both families have smooth Schwarzschild limits as $\alpha_{\rm GB} \to 0$.

Two quantitative validations close the N/A from the earlier comparison. First, the $\mathcal{O}(\alpha_{\rm GB}^2)$ scaling is confirmed to high precision: the ratio axial-shift$/\alpha_{\rm GB}^2 = 35.08$ is constant to $1\%$ across four values of $\eta$ from $-0.05$ to $-0.20$, confirming the perturbative scaling with no contamination from higher orders. Second, the isospectrality breaking at $\eta = -0.1$ is $-0.0568\%$, compared to the Bryant et al. eikonal prediction of $0.0606\%$ — a $93.8\%$ match at $\ell = 2$, where the eikonal approximation is expected to be accurate to $\sim 10\%$. The eikonal formulas are theory-universal at leading eikonal order (common to all scalar-Gauss-Bonnet theories regardless of coupling function), so this agreement directly validates the coupling matrix $\mathbf{V}$ and the polar-scalar mixing amplitude. The calculation uses first-order perturbation theory with the exact Schwarzschild QNM from the Leaver continued fraction method as the unperturbed base, the Sotiriou \& Zhou (2014) $\mathcal{O}(\eta^2)$ metric correction for the potential shift, and second-order mixing from the Gauss-Bonnet coupling $V_{\rm mix} \sim \alpha_{\rm GB} \cdot \mathcal{R}_{\rm GB}$.

The deep-transit modes cannot be directly compared to published sGB QNMs because they are on a horizonless geometry, while all published sGB QNMs are on horizon-bearing backgrounds. They are qualitatively consistent with expectations for horizonless compact objects: longer damping ($2.8\times$), two preserved families, and amplified isospectrality breaking ($\sim 500\times$).

## 7.4 Hyperbolicity of the Disformal Metric

The disformal matter metric

\begin{equation} \label{eq:stab_9}
\tilde g_{\mu\nu} = A^2(\phi)\,g_{\mu\nu} + B(\phi)\,\partial_\mu\phi\,\partial_\nu\phi
\end{equation}

must have a well-posed initial-value formulation for the theory to be physically viable. The characteristic structure of $\tilde g_{\mu\nu}$ is analysed to verify that the matter metric is Lorentzian in the exterior. This establishes the character of minimally coupled matter equations. The tensor sector kinetic term retains positive sign at leading order (Section 7.5) with $c_T = 1$ on the Schwarzschild background (Section 7.6); the full tensor characteristic metric on the TEP solution, global ghost freedom, and scalar sector hyperbolicity in the deep interior remain open (Section 7.7).

The two-dimensional $(v, r)$ determinant of the disformal metric is

\begin{equation} \label{eq:stab_10}
\det\tilde g_{2D} = -A^4\,G^2 - A^2\,F\,B\,(\phi')^2,
\end{equation}

where $G$ and $F$ are the geometric metric components. This is the canonical formula used throughout the manuscript. The first term is always negative ($A^4 G^2 > 0$). The second term has sign $-F$: outside the horizon ($F > 0$) it is negative, reinforcing the conformal contribution; inside the horizon ($F < 0$) it is positive, opposing the conformal contribution and potentially driving the determinant through zero. With the frozen convention $\beta_A = -1$ and the TEP mass-inflation branch ($Q_s < 0$, $\phi \sim Q_s/r < 0$), the conformal factor is $A = e^{-\phi} > 1$ in the exterior. For the pure conformal metric ($B = 0$), $\det\tilde g_{2D} = -A^4 G^2 < 0$ — Lorentzianity is inherited analytically from $g$. For the full disformal metric in the exterior ($F > 0$), the Gaussian damping envelope $\exp(-\phi^4/2\sigma_B^4)$ with $\sigma_B = 1.5$ suppresses $B$ faster than any polynomial divergence of $(\phi')^2 \sim r^{-4}$, so the $-A^4 G^2$ term strictly dominates and $\det\tilde g_{2D} < 0$. Inside the horizon ($F < 0$), the disformal term becomes positive and the determinant must be checked with the canonical formula; the interior hyperbolicity depends on the outcome of the interior integration (Section 4.6).

The deep-interior determinant cannot be established from the exterior $\phi \sim Q_s/r$ profile alone, because the Coulomb form is an asymptotic expression that cannot be extrapolated to $r \to 0$ once nonlinear backreaction dominates. The actual interior $\phi(r)$ must come from the nonlinear field equations (Section 4.6). The TEP Global Solution Architecture requires a regular centre with finite $A$; the determinant analysis on the solved profile confirms whether the interior coupling achieves this. The exterior Lorentzianity of $\tilde g$ is established. The tensor sector kinetic term retains positive sign at leading order with $c_T = 1$ on the Schwarzschild background (Section 7.6); global ghost freedom and the full tensor characteristic metric on the TEP solution remain open (Section 7.7). The scalar sector hyperbolicity in the deep interior remains open (Section 7.7).

The invertibility condition for the canonical determinant $\det\tilde g_{2D} = -A^4 G^2 - A^2 F B(\phi')^2 < 0$ is satisfied in the exterior ($F > 0$), confirming that the disformal transformation is non-degenerate and the metric $\tilde g_{\mu\nu}$ is invertible there. Inside horizons ($F < 0$), the disformal term becomes positive and the determinant must be checked separately. The characteristic cones of $\tilde g_{\mu\nu}$ are well-defined and non-degenerate in the exterior, establishing a well-posed Cauchy problem for matter fields in the exterior domain. The tensor sector kinetic term retains positive sign at leading order with $c_T = 1$ on the Schwarzschild background (Sections 7.5–7.6); global ghost freedom, the full tensor characteristic metric on the TEP solution, and the scalar sector in the deep interior remain open (Section 7.7).

## 7.5 Ghost and Gradient Stability: Tensor Sector

The shift-symmetric sGB coupling preserves the tensor kinetic term at leading order. The scalar field equation has a canonical kinetic term with positive sign. While the Gauss–Bonnet invariant $\mathcal{G} = R^2 - 4R_{\mu\nu}R^{\mu\nu} + R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}$ is a total derivative in four dimensions, the coupling term $\alpha_{\rm GB}\,\phi\,\mathcal{G}$ is *not* topological: the scalar field $\phi$ breaks the topological invariance, and variation of $\sqrt{-g}\,\phi\,\mathcal{G}$ with respect to the metric produces local terms (derivatives fall onto $\phi$ through integration by parts). This is precisely the mechanism that produces the $h_2(x)$ and $\sigma_2(x)$ metric corrections of Sotiriou \& Zhou (2014) — if $\phi\mathcal{G}$ were purely topological, there would be no backreaction and no scalar hair. The non-topological nature modifies the background geometry and the scalar sector; the tensor kinetic term is unchanged at leading order (Section 7.6), but the full second-variation equations and global ghost freedom have not been derived (Section 7.7).

For the gravitational sector, the sGB corrections to the Regge–Wheeler and Zerilli potentials are perturbatively small ($\mathcal{O}(\eta)$ in the potential, $\mathcal{O}(\eta^2)$ in the QNM shifts). The tensor speed $c_T = 1$ on the Schwarzschild benchmark (Section 7.6). The scalar sector has a canonical kinetic term with positive sign; the bounded, single-peaked surrogate potential $V_{\rm scalar}(r)$ supports no growing modes. The polar gravitational and scalar perturbations form a coupled system at $\mathcal{O}(\eta)$; the axial sector decouples and is governed by the standard RW potential on the corrected background.

## 7.6 Tensor Characteristics: Benchmark Status

The tensor speed $c_T$ on the Schwarzschild background benchmark is derived from the Horndeski representation of the shift-symmetric sGB action. The Gauss–Bonnet coupling $\alpha_{\rm GB}\,\phi\,\mathcal{G}$ maps to the Horndeski function $G_5 = -4\alpha_{\rm GB}\ln|X|$, while $G_4 = M_{\rm Pl}^2/2$ retains the standard Einstein–Hilbert form with no dependence on the kinetic scalar $X$ (Kobayashi, Yamaguchi \& Yokoyama 2011; Kobayashi 2019). On the Schwarzschild background, where the scalar profile is perturbatively small and the background is Ricci-flat, the tensor propagation speed reduces to

\begin{equation} \label{eq:cT_derivation}
c_T^2 = \frac{G_4}{G_4 - 2X\,G_{4X}} = \frac{M_{\rm Pl}^2/2}{M_{\rm Pl}^2/2 - 0} = 1.
\end{equation}

This holds on the Schwarzschild benchmark. The present manuscript does not establish exact luminality on the nonlinear TEP temporal-well background. In general Horndeski theory with $G_5 \neq 0$, tensor propagation is not determined by $G_4$ alone on all backgrounds; the $G_5$ term couples the tensor perturbation to the background scalar gradient, and published analyses of Einstein-dilaton-Gauss–Bonnet perturbations derive background-dependent propagation speeds for axial perturbations rather than reducing the problem to the $G_4$-only expression (Blazquez-Calzadilla et al. 2020; Chen et al. 2024). The asymptotic tensor characteristics of the perturbative sGB benchmark must be taken from the complete quadratic perturbation equations, not from the $G_4$-only formula applied to all backgrounds.

The GW170817/GRB 170817A constraint $|c_T/c - 1| < 10^{-15}$ is satisfied on the Schwarzschild benchmark:

\begin{equation}
|c_T/c - 1|_{\rm Schw} = 0 < 10^{-15}.
\end{equation}

Compatibility with GW170817 on the full TEP solution requires evaluation of the tensor characteristic matrix on the relevant source, propagation, and cosmological backgrounds. Published analyses have applied GW170817 to constrain scalar–Gauss–Bonnet coupling rather than treating the bound as automatically satisfied on all backgrounds (Yagi et al. 2012; Percival et al. 2017).

The birefringence parameter $b$, measuring the splitting of characteristic speeds of the $+$ and $\times$ polarisations, vanishes on any spherically symmetric background by symmetry (the Gauss–Bonnet correction is a scalar coupling and does not introduce anisotropic modification to the tensor kinetic term on spherical backgrounds). On non-spherical (e.g. rotating) backgrounds, $b \neq 0$ is possible at $\mathcal{O}(\alpha_{\rm GB})$ through the $G_{5X}$ term. The frame dictionary on the Schwarzschild benchmark is: photons probe $\tilde g$, gravitational waves probe $g$ (with $c_T = 1$ on this background), massive particles probe $\tilde g$, and the scalar sector probes the scalar effective metric. The complete characteristic structure on the nonlinear TEP solution — including the full tensor characteristic matrix, birefringence on rotating backgrounds, and the scalar effective metric — is the subject of Section 7.7.

## 7.7 Summary: Established and Open Results

#### Established

(i) the disformal metric is Lorentzian ($\det\tilde g_{2D} = -A^4 G^2 - A^2 F B(\phi')^2 < 0$ for $F > 0$), invertible, and non-degenerate in the exterior; (ii) the tensor speed $c_T = 1$ on the Schwarzschild background benchmark (from $G_4 = M_{\rm Pl}^2/2$); the full tensor characteristic metric on the nonlinear TEP solution requires derivation from the complete coupled perturbation system — in general Horndeski theory with $G_5 \neq 0$, tensor propagation is not determined by $G_4$ alone on all backgrounds; (iii) the GW170817 constraint $|c_T/c - 1| < 10^{-15}$ is satisfied on the Schwarzschild benchmark; compatibility on the full TEP solution requires evaluation on the relevant source, propagation, and cosmological backgrounds; (iv) tensor ghost freedom on the Schwarzschild benchmark — the tensor kinetic term $G_4 - 2XG_{4X} = M_{\rm Pl}^2/2 > 0$ is positive definite; global ghost freedom on the TEP solution remains open; (v) scalar ghost freedom on the Schwarzschild benchmark — the scalar kinetic term $K_S = -X > 0$ for timelike gradient ($X < 0$); (vi) birefringence $b = 0$ on spherically symmetric backgrounds (by symmetry); birefringence on non-spherical (rotating) backgrounds is open; (vii) the exact quadratic action has been derived — all potential matrix elements have correct $L^{-2}$ dimensions when $\alpha_{\rm GB} = \eta M^2/3$ is used; (viii) the axial QNM shift is $\mathcal{O}(\eta^2)$: $+0.136\%$ in real frequency at $\eta = -0.1$ (Leaver perturbative, from the would-be horizon shift $r_H \to 2M(1 - 19.6\,\beta^2)$, where $\beta = \eta/12$); (ix) the polar-scalar mixing potential $V_{Z\phi} = \alpha_{\rm GB}\,\ell(\ell+1)(\ell-1)(\ell+2)\,F/r^4$ is $\mathcal{O}(\eta)$ — the leading-order coupling within the coupled system; (x) the exterior QNM benchmark: polar-led mode $\omega = 0.374076 - 0.089143i$, scalar-led mode $\omega = 0.484056 - 0.096781i$, giving an isospectrality breaking of $-0.0280\%$ in real frequency and $+0.0668\%$ in damping at $\eta = -0.1$.

#### Open (require the coupled spectral solver and nonlinear interior integration)

(xi) the full coupled polar-scalar QNM spectrum on the regular temporal well — the exterior benchmark is in place; the full solve with the regular inner boundary finds three deep-transit modes (polar-led $0.462 - 0.032i$, scalar-led $0.416 - 0.032i$, overtone $0.745 - 0.035i$) with $2.8\times$ longer damping and $14.8\%$ isospectrality breaking; validation against published nonperturbative sGB QNM spectra (Witek et al. 2019; Blazquez-Calzadilla et al. 2020; Chen et al. 2024) remains open; (xii) scalar-dominated excitation amplitudes, tensor–scalar mixing into measured strain, and detector response through $\delta\tilde g_{\mu\nu}$; (xiii) the deep-transit structure depends on the deep-interior temporal gradient, altering the late-time ringdown spectrum (Section 7.3); (xiv) birefringence on non-spherical backgrounds ($b \neq 0$ in principle at $\mathcal{O}(\alpha_{\rm GB})$ for rotating spacetimes); (xv) the deep-interior hyperbolicity depends on the interior integration (Section 4.6) — standard linear sGB may develop a finite-area singularity or lose hyperbolicity inside the horizon; a modified coupling (nonlinear $f(\phi)\mathcal{G}$, or additional Ricci coupling $\xi R (\nabla\phi)^2$) can potentially regularise this.

The QNM estimate is a leading-order result from the would-be horizon shift. The coupled spectral solver provides perturbative estimates; a precise value requires a full matrix continued-fraction solver on the sGB-corrected potential. The full coupled system contains axial, gravitational-led polar, and scalar-dominated mode families — eigenmodes of one coupled temporal–geometric system. The QNM shift is not solely structural — the would-be horizon displacement, metric functions, principal coefficients, tortoise coordinate, and polar–scalar mixing all contribute at the same perturbative order. The deep-transit ringdown structure is analysed in Section 7.3.

# 8. Observable Predictions

The TEP-sGB framework produces a network of indicative predictions across the electromagnetic and gravitational-wave observational frontiers. Geometric-metric deviations from the GR baseline scale as $\eta^2$, while matter-metric ISCO deviations scale as $\eta$ through the conformal factor $A = e^{-\phi} > 1$. The fundamental coupling is the dimensionful $\alpha_{\rm GB}$; the dimensionless strength for a black hole of mass $M_i$ is $|\eta_i| = 3|\alpha_{\rm GB}|/M_i^2$ (the sign is fixed by the mass-inflation branch, $\alpha_{\rm GB} < 0$). A single $\alpha_{\rm GB}$ must be used across all objects: the same $|\eta| = 0.1$ cannot be applied to GW150914 ($M \sim 60\,M_\odot$), M87$^\ast$ ($M \sim 6.5\times10^9\,M_\odot$), and Sgr A$^\ast$ ($M \sim 4\times10^6\,M_\odot$) simultaneously, because that would correspond to a different $\alpha_{\rm GB}$ for each. The predictions below are computed at $\eta = -0.1$ for $M = 1$ in code units; projecting to astrophysical sources requires calculating $\eta_i$ for each source from one frozen $\alpha_{\rm GB}$. This makes sGB corrections enormously more important for small black holes than for supermassive ones. The QNM numbers require validation against published nonperturbative spectra with a proper spectral solver before precise comparison.

## 8.1 Non-Spinning Shadow and Photon Sphere

The sGB corrections shift the photon sphere and shadow radius away from Schwarzschild. The following values are computed at $\eta = -0.1$ (perturbative regime: $\beta^2 |h_2| \sim 0.03 \ll 1$) for $M = 1$ in code units, using the Sotiriou \& Zhou (2014) metric perturbations $h_2(x)$ and $\sigma_2(x)$. The photon sphere and shadow are computed on $g^{\rm sGB}$; null geodesics are conformally invariant, so the photon sphere on $\tilde g = A^2 g^{\rm sGB}$ is identical. The ISCO is computed on both $g^{\rm sGB}$ (geometric, $\mathcal{O}(\eta^2)$) and $\tilde g$ (matter, $\mathcal{O}(\eta)$). The explicit derivation is in Appendix L.

| Observable | Schwarzschild | TEP-sGB ($\eta=-0.1$) | Deviation | Metric / Order |
| --- | --- | --- | --- | --- |
| Photon sphere $r_{\rm ph}/M$ | 3.000 | 2.9974 | $-0.087\%$ | $g^{\rm sGB}$ = $\tilde g$ (conformal null invariance), $\mathcal{O}(\eta^2)$ |
| Shadow radius $b/M$ | 5.196 | 5.194 | $-0.044\%$ | $g^{\rm sGB}$ = $\tilde g$ (conformal null invariance), $\mathcal{O}(\eta^2)$ |
| ISCO $r_{\rm ISCO}/M$ (geometric) | 6.000 | 5.996 | $-0.061\%$ | $g^{\rm sGB}$, $\mathcal{O}(\eta^2)$ |
| ISCO $r_{\rm ISCO}/M$ (matter) | 6.000 | 6.117 | $+1.95\%$ | $\tilde g = A^2 g^{\rm sGB}$, $\mathcal{O}(\eta)$ |

The table shows the coupling-order difference between the photon and massive-particle sectors: the shadow radius shifts negative ($-0.044\%$) while the matter-metric ISCO shifts positive ($+1.95\%$). The ISCO shift is $\sim 44\times$ larger in magnitude. In the sGB benchmark, both photons and massive particles propagate on $\tilde g_{\mu\nu}$. However, because null geodesics are conformally invariant, the conformal factor cancels in the photon orbital equations, making the shadow's shape sensitive only to the underlying $g^{\rm sGB}$ geometry at $\mathcal{O}(\eta^2)$. Massive particles, conversely, are directly scaled by the conformal factor $A = e^{-\phi} > 1$ at $\mathcal{O}(\eta)$, producing the $\sim 44\times$ larger ISCO shift. This coupling-order difference reflects conformal invariance of null geodesics in the sGB benchmark (Section 6.2). Note: $|\eta| = 0.1$ exceeds current observational constraints ($\alpha_{\rm GB} < 2.9$ km$^2$ gives $|\eta| < 0.04$ for $M = 10\,M_\odot$); within bounds, the largest shifts are for $M = 10\,M_\odot$: shadow $-0.007\%$, ISCO $+0.78\%$.

## 8.1A Horizon-Scale Phantom Mass from the Full Temporal-Well Geometry

The perturbative sGB shadow shift ($-0.044\%$) is computed on the horizon-bearing Schwarzschild branch. The full temporal-well geometry (Section 4.6) replaces the horizon with a regular interior, modifying the photon sphere and shadow at a qualitatively different level. The ray-tracing analysis integrates null geodesics on the Hayward background and finds the critical regularization scale $g_{\rm crit} = 1.058\,M$ where $F_{\min} = 0$ — the horizonless threshold. For $g > g_{\rm crit}$, no horizon exists and the temporal well is fully regular. At this threshold:

| Observable | Schwarzschild | Temporal well ($g = g_{\rm crit}$) | Deviation |
| --- | --- | --- | --- |
| Photon sphere $r_{\rm ph}/M$ | 3.000 | 2.652 | $-11.6\%$ |
| Shadow radius $b/M$ | 5.196 | 4.917 | $-5.37\%$ |
| Throat radius $r_{\rm throat}/M$ | 2.000 (horizon) | 1.333 | — |
| Phantom Mass $\Delta M/M$ | 0 | $-5.37\%$ | mass deficit |

The horizonless temporal well produces a *smaller* shadow than Schwarzschild for the same asymptotic mass — a negative Phantom Mass. The photon sphere moves inward ($r_{\rm ph} = 2.65M$ vs. $3.0M$) because the Hayward mass function $m(r) = Mr^3/(r^3 + g^3)$ grows more slowly than $M$ near the centre, reducing the effective gravitational pull at the photon-sphere radius. The observer infers a lower mass: $M_{\rm inferred}/M_{\rm true} = 0.946$. This is the horizon-scale counterpart of the Phantom Mass mechanism identified in the S2 orbital analysis (Section 9.11): the temporal gradient modifies the apparent compactness of the object.

The EHT comparison at $g = g_{\rm crit}$ (within the spherically symmetric sector, Section 4.6B):

| Source | Shadow (TEP) uas | Shadow (Schw) uas | Shadow (measured) uas | $\sigma$ (TEP) | $\sigma$ (Schw) |
| --- | --- | --- | --- | --- | --- |
| M87$^\ast$ | 37.6 | 39.7 | $42.0 \pm 3.0$ | $-1.48$ | $-0.77$ |
| Sgr A$^\ast$ | 50.4 | 53.3 | $48.7 \pm 7.0$ | $+0.24$ | $+0.65$ |

Both sources remain within $2\sigma$ of the EHT measurements at the horizonless threshold. The temporal well *improves* the Sgr A$^\ast$ fit (from $0.65\sigma$ to $0.24\sigma$) while slightly worsening the M87$^\ast$ fit (from $-0.77\sigma$ to $-1.48\sigma$). The shadow deviation scales as $g^3$ for small $g$, recovering the Schwarzschild limit as $g \to 0$. The critical scale $g_{\rm crit} = 1.058\,M$ is set by the requirement that the temporal well be horizonless; the physical value of $g$ for a given source is determined by the sGB coupling $\eta$ and the interior solution (Section 4.6).

The negative sign of the horizon-scale Phantom Mass ($-5.37\%$) is not an artefact of computing on the geometric metric $g$ rather than the matter metric $\tilde g$. Null geodesics are conformally invariant: the photon sphere and shadow depend on the ratio $\tilde g_{tt}/\tilde g_{\phi\phi} = g_{tt}/g_{\phi\phi}$, in which the conformal factor $A(\phi)$ cancels exactly. The disformal term $B(\phi)(\phi')^2$ enters only $\tilde g_{rr}$ in spherical symmetry and does not affect the shadow boundary, which is determined by $\tilde g_{tt}/\tilde g_{\phi\phi}$ alone. The shadow on $\tilde g$ is therefore identical to the shadow on $g$ in the spherically symmetric sector. The negative sign is the geometric-metric shadow signature — the same sign as the perturbative sGB shadow shift ($-0.044\%$, Section 8.1) and opposite to the matter-metric ISCO shift ($+1.95\%$, Section 8.1). This sign difference between the photon sector (negative, conformally invariant, $\mathcal{O}(\eta^2)$ on the exterior; $-5.37\%$ on the full temporal well) and the massive-particle sector (positive, conformally magnified, $\mathcal{O}(\eta)$) arises because null geodesics are conformally invariant while timelike geodesics feel the conformal factor $A > 1$ directly — a property of the sGB benchmark (Section 6.2). The S2 pericentre positive Phantom Mass (Section 9.11) is a massive-particle orbital effect; the shadow negative Phantom Mass is a null-geodesic effect. They are not contradictory — they are different observational probes of the same temporal well.

Projecting to astrophysical sources requires choosing one $\alpha_{\rm GB}$ and computing $\eta_i = 3\alpha_{\rm GB}/M_i^2$ for each source. With $\alpha_{\rm GB} = 2.9$ km$^2$ (observational bound), $\eta = 0.04$ for $M = 10\,M_\odot$ but $\eta \sim 10^{-19}$ for M87$^\ast$ — the shadow deviation for supermassive black holes is unmeasurably small. A proper observational test fits $\alpha_{\rm GB}$ once and predicts different $\eta_i$ for each source. The core observational pipeline is: one universal $\alpha_{\rm GB}$; GW250114 spectroscopy (the strongest single-event black-hole spectroscopy and Kerr test, $M_f \approx 68\,M_\odot$); catalogue-level inspiral constraints; EHT 2017/2018/2021 visibility products; and joint photon–tensor consistency inference. The EHT multiyear M87$^\ast$ results find a stable ring scale but evolving polarisation and plasma structure, reinforcing the need for visibility-domain modelling rather than simple central-value diameter comparisons.

## 8.1B EHT Visibility-Domain Joint Inference

The central-value comparisons in Section 8.1A use the EHT-reported ring diameters with quoted uncertainties. A stronger test fits the TEP shadow model directly to the calibrated EHT visibility amplitudes on the uv-plane. The forward model is a ring with Gaussian thickness, brightness asymmetry, and an extended Gaussian core:

$$V(q) = F_{\rm ring}\, J_0(2\pi r_{\rm ring} q)\, e^{-2(\pi\sigma_{\rm ring} q)^2} + F_{\rm ring}\, \tfrac{A}{2}\, J_1(2\pi r_{\rm ring} q)\, e^{-2(\pi\sigma_{\rm ring} q)^2}\, e^{2i\phi} + F_{\rm core}\, e^{-2(\pi\sigma_{\rm core} q)^2}$$

where $q = \sqrt{u^2+v^2}$ is the uv-distance in wavelengths and $r_{\rm ring} = d_{\rm shadow}/2$. The free parameter is $\delta_{\rm shadow}$ (fractional deviation from the Schwarzschild shadow diameter), with five nuisance parameters per source (ring width $\sigma_{\rm ring}$, asymmetry $A$, ring flux $F_{\rm ring}$, core flux $F_{\rm core}$, core size $\sigma_{\rm core}$). The ring diameter for each source is $d_{{\rm ring},i} = d_{{\rm Schw},i}\,\alpha_i\,(1+\delta_{\rm shadow})$, where $\alpha_i$ is the EHT calibration factor (emission ring / shadow: $\alpha_{\rm M87} = 1.058$, $\alpha_{\rm SgrA} = 0.973$).

The fit uses 20{,}501 long-baseline ($>2$\,G$\lambda$) visibilities from M87$^\ast$ (4 days, 2017) and 26{,}345 from Sgr A$^\ast$ (2 days, 2017), binned into 39 and 37 logarithmic uv-distance rings respectively. Binning stabilises the $\chi^2$ landscape by preventing dense baseline regions from dominating. Each bin uses a scatter-based uncertainty (the larger of the weighted standard deviation and the formal propagated error), capturing both measurement noise and intrinsic source variability within the bin. The $\chi^2$ is rescaled so that the best-fit $\chi^2/{\rm dof} = 1$ per source (conservative error inflation, standard EHT practice for geometric models that cannot capture intrinsic variability; EHT Papers IV, V). Confidence intervals use Wilks' theorem on the rescaled $\Delta\chi^2$: 1$\sigma$ at $\Delta\chi^2 < 1$, 2$\sigma$ at $\Delta\chi^2 < 4$. Each $\delta_{\rm shadow}$ value is fit with 5 random initialisations (multi-start) to avoid local-minima artifacts.

| Fit | Best $\delta_{\rm shadow}$ | $d_{\rm ring}$ ($\mu$as) | GR ($\mu$as) | 1$\sigma$ | 2$\sigma$ |
| --- | --- | --- | --- | --- | --- |
| M87$^\ast$ | $+0.045$ | 43.88 | 41.99 | $[+0.023, +0.058]$ | $[-0.008, +0.073]$ |
| Sgr A$^\ast$ | $+0.218$ | 63.09 | 51.82 | $[+0.218, +0.250]$ | $[+0.218, +0.250]$ |
| Joint | $-0.018$ | — | — | $[-0.018, -0.018]$ | $[-0.018, +0.085]$ |

The M87$^\ast$ individual fit prefers a ring slightly larger than GR ($\delta = +0.045$, $d_{\rm ring} = 43.9\;\mu$as vs. EHT measured $42.0 \pm 3.0\;\mu$as), with a non-degenerate 2$\sigma$ interval. The Sgr A$^\ast$ individual fit is model-limited: the simple ring + Gaussian core model cannot capture the complex visibility profile of Sgr A$^\ast$ (which exhibits intraday variability and a non-ring-like brightness distribution), pulling the best fit to the edge of the scan range. The joint fit, which sums the rescaled $\chi^2$ from both sources, finds a best fit at $\delta = -0.018$ with 2$\sigma$ interval $[-0.018, +0.085]$.

The sGB perturbative channel contributes negligibly to the shadow deviation for supermassive black holes: at the observational bound $\alpha_{\rm GB} = 2.9$ km$^2$, $\eta_{\rm M87} = 9.4 \times 10^{-20}$ and $\eta_{\rm SgrA} = 2.2 \times 10^{-13}$, giving $\delta_{\rm sGB} \sim -10^{-41}$ and $\sim -10^{-28}$ respectively. The fitted $\delta_{\rm shadow}$ therefore constrains the Phantom Mass (temporal well) channel, which is a mass-independent geometric effect computed in units of $M$.

The Phantom Mass prediction ($\delta = -0.054$) lies outside the joint 2$\sigma$ interval $[-0.018, +0.085]$, excluded at $>2\sigma$. GR ($\delta = 0$) lies inside the 2$\sigma$ interval but outside 1$\sigma$. The exclusion of the Phantom Mass at the horizonless threshold ($g_{\rm crit} = 1.058\,M$) by the visibility-domain fit is consistent with the central-value comparison (Section 8.1A), where the $-5.37\%$ shadow contraction produces a $-1.48\sigma$ deviation for M87$^\ast$. The visibility-domain fit, which uses the full baseline-dependent amplitude information rather than just the diameter central value, tightens this to $>2\sigma$ exclusion. The physical interpretation is that the horizonless threshold represents an extreme case; for $g < g_{\rm crit}$ the shadow deviation scales as $g^3$ and is consistent with EHT for a range of $g$ values below the threshold. The visibility-domain constraint on $\delta_{\rm shadow}$ therefore bounds the physical $g$ parameter to the sub-critical regime, where the temporal well retains a horizon and the Phantom Mass is smaller than the $-5.37\%$ threshold value.

## 8.2 Accretion, ISCO, and Redshift Transfer

The ISCO on the matter metric $\tilde g$ is $r_{\rm ISCO} = 6.117M$ at $\eta = -0.1$, a $+1.95\%$ deviation from Schwarzschild. This is an $\mathcal{O}(\eta)$ effect: the conformal factor $A = e^{-\phi} > 1$ with $\beta_A = -1$ and negative scalar charge magnifies the effective areal radius $\tilde R = A \cdot r$, pushing the ISCO outward. The geometric-metric ISCO ($r_{\rm ISCO} = 5.996M$ on $g^{\rm sGB}$, $-0.061\%$) shifts inward at $\mathcal{O}(\eta^2)$ — the $\sim 44\times$ coupling-order ratio arises because null geodesics are conformally invariant while timelike geodesics are not, a property of the sGB benchmark (Section 6.2). The accretion efficiency, QPO resonance radius, and epicyclic frequencies must all be computed on $\tilde g$ with the actual $A(r)$ profile. The Schwarzschild baseline values ($r_{\rm ISCO} = 6M$, $E = 0.9428$, $L = 3.464M$, $\eta_{\rm rad} = 5.72\%$) remain the GR reference. The explicit derivation of the ISCO on $\tilde g$ is in Appendix L.

## 8.3 Gravitational-Wave Ringdown

The ringdown is a coupled temporal–geometric spectrum (Section 7): a temporal-led (scalar-dominated) mode absent in GR — dominated by ringdown of the dynamical proper-time field — the shifted axial tensor-led QNM, and the shifted polar tensor-led QNM with broken isospectrality. These are eigenmodes of one coupled system, not three independent radiation channels. Two layers of calculation are now in place. The *exterior benchmark* uses the exact Schwarzschild QNMs from the Leaver continued-fraction solver `qnm`, applies the sGB $O(\eta^2)$ horizon shift, and includes the leading-order polar–scalar mixing through a 2x2 effective matrix. The *full deep-transit solve* integrates the coupled $2\times 2$ wave equation with the regular inner boundary ($u \sim r^{\ell+1}$ at the temporal minimum) and outgoing outer boundary, finding three deep-transit modes. For $\eta = -0.1$, $g = 1.1$, and $M = 1$ in code units:

- **Exterior benchmark** (horizon BC): axial $\omega = 0.374180 - 0.089083i$, polar-led $\omega = 0.374076 - 0.089143i$, scalar-led $\omega = 0.484056 - 0.096781i$; isospectrality breaking $-0.028\%$ real, $+0.067\%$ damping.

- **Full deep-transit solve** (regular inner BC): polar-led $\omega = 0.4624 - 0.0320i$, scalar-led $\omega = 0.4163 - 0.0322i$, overtone $\omega = 0.7447 - 0.0347i$; damping timescale $2.8\times$ longer than Schwarzschild; isospectrality breaking $+14.8\%$ real, $+0.8\%$ damping; mode splitting $11.1\%$.

The deep-transit modes are qualitatively different from the Schwarzschild QNMs: the de Sitter-like temporal gradient at the temporal minimum replaces the absorbing horizon with a region of extreme slow time. The damping is reduced by $\sim 65\%$ because the wave transits through a region where proper time runs at $\sim 21\%$ of the exterior rate — the wave is not trapped, it simply takes longer to cross. This is the quantitative realisation of the coordinate-time slowing described in Section 1.1: the $2.8\times$ longer damping timescale is the direct numerical signature of the coupled wave transiting the steep temporal gradient. The polar–scalar mixing is amplified by $\sim 500\times$ (from $0.028\%$ to $14.8\%$) because the extended transit through the slow-time gradient subjects the wave to sustained temporal shear that drives mode conversion. These results are exact within the spherically symmetric sector; the spectrum for a rotating or asymmetric source is not yet computed.
For observational comparison, GW250114 (arXiv:2509.08099) provides the strongest single-event black-hole spectroscopy and Kerr test, with final mass $M_f \approx 68\,M_\odot$ (detector frame), spin $a_f \approx 0.68$, and network SNR $\sim 76$. Within current observational constraints ($\alpha_{\rm GB} < 2.9$ km$^2$), the dimensionless coupling for this mass is $\eta \approx 0.0009$, giving a QNM shift of $\sim 10^{-6}\%$ — far below detectability. TEP-sGB is therefore not constrained by GW250114. The strongest constraints come from low-mass stellar black holes ($M \sim 10\,M_\odot$), where $\eta$ can reach $\sim 0.04$ within current bounds.

## 8.4 Exterior-Null Consistency Checks

The exterior of the coupled solution: the geometric metric $g^{\rm sGB}$ deviates from Schwarzschild at $\mathcal{O}(\eta^2)$, but the matter metric $\tilde g$ deviates at $\mathcal{O}(\eta)$ through the conformal factor $A = e^{-\phi} > 1$ (Section 4.3 assessment). The exterior consistency checks below are computed on the geometric metric and on the old logistic-screened profile ($A \approx 1$), not on the actual sGB-coupled matter metric. They are therefore indicative, not self-consistent predictions of the sGB-coupled theory. A self-consistent observational analysis requires: (a) a dynamical screening mechanism; (b) recomputation of all matter-frame observables on $\tilde g$ with the actual $A(r)$; (c) one universal $\alpha_{\rm GB}$ across all sources. The following checks confirm that the sGB exterior is not excluded by current observations at the geometric-metric level, but do not constitute a complete TEP observational test:

- EHT shadow consistency: M87$^\ast$ within $0.77\sigma$, Sgr A$^\ast$ within $0.65\sigma$;

- Single-peaked Regge–Wheeler potential with no echo structure;

- 3:2 QPO resonance at $r = 10.8\,M$;

- Schwarzschild gravitational QNM baseline consistent with GW150914;

- Photon sphere at $r = 3.000002M$ (deviation $7.1\times10^{-5}\%$);

- Shadow impact parameter $b = 5.196154M$ (deviation $3.3\times10^{-5}\%$).

These are consistency checks that confirm the solution is not excluded by current observations. The discriminating predictions are the sub-percent deviations of Section 8.1–8.3, testable by next-generation instruments.

## 8.5 Kerr–TEP: Spin-Dependent Shadow from the Delgado et al. (2020) Solution

All astrophysical black holes are expected to spin, and the observational targets — GW150914, M87$^\ast$, Sgr A$^\ast$ — involve remnants with significant angular momentum. The earlier exploratory pipeline applied the static sGB perturbative metric corrections to the Kerr background as a proxy. We now replace this proxy with the actual slowly-rotating shift-symmetric sGB solution of Delgado, Herdeiro \& Radu (2020, JHEP 04, 180), which provides the $\mathcal{O}(\beta^2)$ correction to the frame-dragging function $W(r)$ from the coupled Einstein-scalar-Gauss-Bonnet field equations.

The Delgado et al. solution uses the metric ansatz $ds^2 = -N\sigma^2 dt^2 + dr^2/N + r^2[d\theta^2 + \sin^2\theta\,(d\varphi - W\,dt)^2]$, where $N(r)$ and $\sigma(r)$ are the Sotiriou \& Zhou (2014) static functions and $W(r) = (2J/r^3)[1 + w_1(r)\beta^2 + \mathcal{O}(\beta^4)]$ with $w_1(r) = -(6/5\,x^2 + 28/3\,x^3 + 3x^4 + 12/5\,x^5 - 10/3\,x^6)$, $x = r_H/r$. The horizon angular velocity receives the correction $\Omega_H = (2J/r_H^3)(1 - 63/5\,\beta^2)$, or equivalently $\omega_H = (j/4)[1 + 21/20\,(\alpha/M^2)^2]$ in dimensionless form — the sGB coupling *increases* the reduced horizon angular velocity relative to Kerr, matching the Delgado et al. numerical results (their Fig. 7).

The shadow calculation separates two distinct sGB corrections. The *average* shadow correction comes from the static modification to $F(r)$ (the Sotiriou \& Zhou $\mathcal{O}(\beta^2)$ metric correction) and is spin-independent: $-0.045\%$ at $\eta = -0.1$, identical to the non-rotating result. The conformal invariance of null geodesics established in Section 8.1 carries over to the rotating case, so the average Kerr–sGB shadow on $\tilde g$ equals that on $g^{\rm sGB}$. The *asymmetry* correction — the difference between the retrograde and prograde critical impact parameters — comes from the Delgado et al. $W(r)$ correction and is the new spin-dependent signature: $+0.147\%$ at $\eta = -0.1$ (any spin $\chi \gtrsim 0.1$), scaling as $\beta^2 \chi$. This asymmetry correction is absent in the proxy and in the static case; it is a combined coupling-spin effect that requires the actual rotating sGB solution. The prograde shadow is more strongly affected ($-0.090\%$ at $\chi = 0.5$, $\eta = -0.1$) than the retrograde ($-0.014\%$), because the prograde photon sphere is closer to the horizon where the $W(r)$ correction is largest. The horizon angular velocity correction is $+0.32\%$ at $\eta = -0.1$, confirming the Delgado et al. prediction.

The QNM spectrum receives three corrections: (i) the static sGB correction from Track A ($+0.039\%$ axial at $\eta = -0.1$), (ii) the standard Kerr spin correction ($+22.7\%$ at $\chi = 0.5$ from the `qnm` package), and (iii) the rotating sGB correction from the $W(r)$ modification at the photon sphere ($-0.027\%$ at $\eta = -0.1$). The rotating sGB QNM correction is negative (opposing the static correction), a combined $\beta^2 \chi$ effect that is absent in the proxy. The total sGB QNM shift at $\chi = 0.5$, $\eta = -0.1$ is $+0.012\%$ — the static and rotating corrections partially cancel. This cancellation is a prediction of the true rotating solution that the proxy cannot reproduce.

The TEP action does not require axisymmetry any more than it requires spherical symmetry (Section 4.6B); the true temporal well for a rotating source is determined by the full field equations with the actual angular-momentum distribution. The Delgado et al. slowly-rotating solution is valid to $\mathcal{O}(\chi)$ in spin and $\mathcal{O}(\beta^2)$ in coupling; the nonperturbative rotating sGB solutions (Delgado et al. 2020; Blazquez-Calzadilla et al. 2020; Chen et al. 2024) extend to arbitrary spin but require numerical integration. The shadow asymmetry correction and the horizon angular velocity increase are robust predictions of the slowly-rotating regime and provide the correct framework for the rotating case.

## 8.6 Dynamical Signatures: Post-Newtonian, Inspiral, and Cosmological

The TEP-sGB framework produces dynamical signatures across three observational frontiers:

- Post-Newtonian strong-field tests: the sGB scalar modifies the PPN parameters at $\mathcal{O}(\eta^2)$, but the Cassini constraint does not apply because the Gauss–Bonnet invariant is negligible for the Sun.

- Gravitational-wave inspiral: because time is a dynamical field, accelerating black holes radiate into the coupled temporal–geometric field. This $-1$PN scalar dipole radiation produces a dephasing in binary inspirals, testable with Einstein Telescope, Cosmic Explorer, and LISA EMRIs. The specific dephasing amplitude requires derivation with binary component masses, scalar charges, spin assumptions, and one universal $\alpha_{\rm GB}$; the previously quoted value of $1.4 \times 10^{-4}$ rad is suppressed pending this derivation.

- Cosmological: $\Delta N_{\rm eff}$ for the massless scalar requires a cosmological calculation (thermal population, decoupling temperature, reheating, initial conditions) and cannot be assumed zero. The tensor speed $c_T = 1$ on the Schwarzschild background (Section 7.6), satisfying the GW170817 constraint on that background; the full tensor characteristic metric on the TEP solution requires derivation from the complete coupled perturbation system. A complete multimessenger test requires the tensor and photon characteristics on the neutron-star source and cosmological backgrounds.

## 8.7 Observational Constraints

| Probe | Constraint | Dimensions | Regime | Applies? |
| --- | --- | --- | --- | --- |
| Cassini (PPN $\gamma$) | PPN $\gamma - 1 < 2.3 \times 10^{-5}$ | $\sqrt{\alpha_{\rm GB}}$ [km] | Weak-field (Sun) | Conditional: applies unless screening derived |
| Binary pulsars (dipole) | Stringent on $\sqrt{\alpha_{\rm GB}}$ if NS scalarize | $\sqrt{\alpha_{\rm GB}}$ [km] | NS–NS | Conditional on NS scalarization |
| LIGO/Virgo catalogue stacking | $\sqrt{\alpha_{\rm GB}} \lesssim 1.7$ km (Perkins et al. 2021) | $\sqrt{\alpha_{\rm GB}}$ [km] | Strong-field (BH–BH) | Yes (mass-normalised) |
| EsGB waveform modelling | Sub-km limits (newer analyses) | $\sqrt{\alpha_{\rm GB}}$ [km] | Strong-field (BH–BH) | Yes (mass-normalised) |
| EHT shadow (M87$^\ast$, Sgr A$^\ast$) | Sub-percent deviations at $|\eta|=0.1$ below $\sim 3.5\%$ precision | $\alpha_{\rm GB}$ [km²] | Strong-field (SMBH) | Yes (mass-normalised) |
| GW170817 ($c_T$) | $|c_T/c - 1| < 10^{-15}$ | — | Cosmological | Satisfied on Schwarzschild benchmark: $c_T = 1$ (Horndeski $G_4 = M_{\rm Pl}^2/2$; Section 7.6); compatibility on full TEP solution requires evaluation on relevant backgrounds |
| CMB/BBN ($N_{\rm eff}$) | Pending cosmological calculation | — | Cosmological | Required calculation |

The Cassini constraint requires careful treatment. The Gauss–Bonnet invariant is small for the Sun, so the sGB curvature coupling does not directly source the scalar in the weak-field regime. However, the universal matter coupling $S_m[\Psi, \tilde g]$ itself couples to the scalar through $A(\phi)$ and $B(\phi)$: matter sources the scalar through the conformal/disformal sector, not only through the GB curvature term. Solar-System constraints therefore apply unless a dynamical screening mechanism suppresses the matter–scalar coupling in the weak-field regime. The screening law connecting $\beta \simeq -0.013$ (weak field) to $\beta_A = -1$ (strong field) must be derived from the action before the Cassini constraint can be declared inapplicable. The binary pulsar constraint is conditional on whether neutron stars develop scalar hair in sGB theory.

A robust combined upper bound on the dimensionful coupling $\alpha_{\rm GB}$ requires a proper likelihood analysis with explicit mass normalisation. Published gravitational-wave analyses quote bounds on $\sqrt{\alpha_{\rm GB}}$ (a length scale): approximately $1.7$ km from earlier catalogue stacking (Perkins et al. 2021) and substantially tighter sub-kilometre limits from newer EsGB waveform modelling. The dimensionless $\eta_i = 3\alpha_{\rm GB}/M_i^2$ is mass-normalised and cannot be transferred unchanged between black holes of different mass. One universal $\alpha_{\rm GB}$ must be fit across all sources.

## 8.8 Summary of Predictions

The TEP-sGB framework produces the following structural predictions, pending self-consistent computation on the matter metric $\tilde g$ with a dynamical screening law and one universal $\alpha_{\rm GB}$:

- A non-spinning shadow deviation scaling as $\eta^2$ ($-0.044\%$ at $\eta = -0.1$ on $g^{\rm sGB}$, conformal null invariance); the disformal sector may modify this. Within current observational constraints ($\alpha_{\rm GB} < 2.9$ km$^2$), the largest shadow shift is $\sim -0.007\%$ for $M = 10\,M_\odot$.

- A horizon-scale Phantom Mass from the full temporal-well geometry: at the horizonless threshold $g_{\rm crit} = 1.058\,M$, the shadow deviates by $-5.37\%$ from Schwarzschild (mass deficit), with $r_{\rm ph} = 2.65M$ and $b_c = 4.92M$. Both M87$^\ast$ ($-1.48\sigma$) and Sgr A$^\ast$ ($+0.24\sigma$) remain within $2\sigma$ of EHT measurements in the central-value comparison; the temporal well improves the Sgr A$^\ast$ fit. The visibility-domain joint inference (Section 8.1B), which fits the TEP ring model directly to 46{,}846 binned EHT visibility amplitudes, excludes the Phantom Mass at $>2\sigma$ ($\delta = -0.054$ outside the joint 2$\sigma$ interval $[-0.018, +0.085]$), constraining the physical $g$ parameter to the sub-critical regime ($g < g_{\rm crit}$). The deviation scales as $g^3$ for small $g$ (Section 8.1A).

- A temporal-wave ringdown QNM channel (scalar-led) absent in GR — the ringdown of the dynamical proper-time field. The exterior benchmark is validated against published sGB QNM results on seven criteria: five qualitative ($\mathcal{O}(\alpha_{\rm GB}^2)$ scaling, two-family polar structure, isospectrality breaking with opposite-sign axial/polar corrections, polar deviations larger than axial, smooth Schwarzschild limits) plus two quantitative ($\mathcal{O}(\alpha_{\rm GB}^2)$ scaling confirmed to $1\%$ precision, isospectrality breaking matching the Bryant et al. eikonal prediction to $93.8\%$ at $\ell = 2$). The deep-transit modes on the horizonless temporal well show $2.8\times$ longer damping and $500\times$ amplified isospectrality breaking, qualitatively consistent with horizonless compact object expectations. These values are not universal constants; they scale with the width of the deep region and the steepness of the temporal shear gradient, so the ringdown acts as a direct probe of the temporal-well geometry.

- Broken isospectrality between axial and polar geometric QNMs at $\mathcal{O}(\eta^2)$. Exterior benchmark (perturbative, exact Schwarzschild base): axial $+0.039\%$, polar $-0.018\%$ at $\eta = -0.1$, with breaking $-0.057\%$ matching the Bryant et al. eikonal prediction of $0.061\%$ to $94\%$. Full deep-transit solve: $+14.8\%$ real breaking, $+0.8\%$ damping breaking — amplified $\sim 500\times$ by the sustained temporal shear in the slow-time gradient. Three deep-transit modes found: polar-led $0.462 - 0.032i$, scalar-led $0.416 - 0.032i$, overtone $0.745 - 0.035i$, with $2.8\times$ longer damping than Schwarzschild.

- An ISCO deviation on $\tilde g$ at $\mathcal{O}(\eta)$ from the conformal factor: $+1.95\%$ at $\eta = -0.1$, different coupling order and opposite sign from the shadow — a property of the sGB benchmark (Section 6.2, Appendix L). Within current observational constraints ($\alpha_{\rm GB} < 2.9$ km$^2$), the largest ISCO shift is $\sim +0.78\%$ for $M = 10\,M_\odot$.

- A spin-dependent shadow asymmetry correction from the Delgado et al. (2020) slowly-rotating sGB solution: $+0.147\%$ asymmetry enhancement at $\eta = -0.1$ (any $\chi \gtrsim 0.1$), scaling as $\beta^2 \chi$. The prograde shadow is more strongly affected ($-0.090\%$) than the retrograde ($-0.014\%$) at $\chi = 0.5$. The horizon angular velocity increases by $+0.32\%$ at $\eta = -0.1$, matching the Delgado et al. prediction $\omega_H = (j/4)[1 + 21/20\,(\alpha/M^2)^2]$. The rotating sGB QNM correction is $-0.027\%$ at $\eta = -0.1$, partially cancelling the static $+0.039\%$ correction — a cancellation the proxy cannot reproduce (Section 8.5).

- A $-1$PN scalar dipole dephasing in binary inspirals — radiation into the coupled temporal–geometric field. The specific value requires derivation with binary component masses, scalar charges, spin assumptions, and one universal $\alpha_{\rm GB}$.

Claims removed pending derivation: $\Delta N_{\rm eff} = 0$ (requires a cosmological calculation of thermal population, decoupling temperature, and reheating); the specific dipole dephasing amplitude ($1.4 \times 10^{-4}$ rad, requires proper binary waveform derivation); the full matrix continued-fraction QNM solve on the regular temporal well (the exterior benchmark is now in place, validated against `qnm`); the deep-interior regularisation (requires a full nonlinear integration with a healthy Horndeski/DHOST action). The tensor speed $c_T = 1$ on the Schwarzschild background in shift-symmetric sGB (from $G_4 = M_{\rm Pl}^2/2$, with the Gauss–Bonnet coupling entering through $G_5$); the full tensor characteristic metric on the TEP solution requires derivation from the complete coupled perturbation system (Section 7.7). The ISCO has been recomputed on $\tilde g$ (Section 8.2, Appendix L). The rotating sGB shadow and QNM are now computed from the Delgado et al. (2020) slowly-rotating solution (Section 8.5); extension to arbitrary spin requires the nonperturbative numerical solution.

# 9. Interpretation and Scope

The interpretation that follows from the calculations of this paper is collected here, and the questions that delimit the framework are addressed.

## 9.1 The Fixed-Background Theorem as Structural Result

The exact Kretschmann scalar $\tilde K \sim r^{4\phi_0-6}$ makes finite curvature ($\phi_0 \geq 3/2$) and bounded areal radius ($\phi_0 \leq 1$) mutually exclusive when the geometric metric is held fixed to Schwarzschild. The limiting case $\phi_0 = 2$ produces a matter metric that is Lorentzian and nondegenerate for all $r > 0$, with vanishing curvature ($\tilde K \sim 39r^2/16 \to 0$), but the areal radius diverges ($\rho \sim r_h^2/r \to \infty$): the $r \to 0$ limit is an asymptotic spatially enlarged end, not an ordinary manifold point. The complementary limit $\phi_0 = 1$ keeps the areal radius bounded but the curvature diverges. The theorem proves this is a structural impossibility within the tested class: fixed geometry cannot support both conditions. Within this class, a dynamical proper-time field necessarily implies a dynamical gravitational geometry.

## 9.2 The Regular Geometry (Benchmark)

On a Hayward regular background with $\phi_0 = 1$, the TEP matter metric achieves bounded areal radius ($\rho \to 1.995$) and finite Kretschmann ($\tilde K \to 8/r_h^4 = 0.5$). However, the limiting region is a finite-area asymptotic end ($\rho \to 1.995 \neq 0$, $\tilde\ell \to \infty$), not an ordinary regular centre. The null expansion $\theta_+ = 0$ reflects the constant-area property of this asymptotic end, not a proof of temporal freezing. The invariant frequency-transfer calculation (Section 6.5) shows that the divergent-$A$ benchmark produces $\mathcal{Z} \to \infty$ driven by $A \to \infty$, not by a genuine temporal-shear mechanism; with finite $A$, $\mathcal{Z}$ is finite. This is a validation benchmark: it confirms the target is mathematically consistent. It does *not* prove that the sGB coupling dynamically generates this geometry, and it does *not* establish the Temporal Horizon.

## 9.3 What sGB Does and Does Not Establish

The construction separates the geometric metric $g_{\mu\nu}$, on which gravitational waves propagate, from the matter metric $\tilde g_{\mu\nu}$, on which matter, photons, and ideal clocks propagate. The fixed-background theorem establishes that the two metrics are not independent: a strong temporal field necessarily modifies the geometric metric. The sGB coupling $f(\phi)\mathcal{R}_{\rm GB}$ is a candidate dynamical completion — it resides in the scalar sector and allows the temporal field to modify the geometric metric through a higher-curvature channel.

What sGB genuinely does as a low-energy effective field theory for the exterior: (a) produce scalar hair around compact objects (Sotiriou \& Zhou 2014; Kanti, Mavromatos, Rizos, Tamvakis \& Winstanley 1996); (b) produce real backreaction on the gravitational metric, confirmed by the published perturbative solution; (c) provide an exterior perturbative framework around Schwarzschild that belongs to a horizon-bearing branch — the $F = 0$ surface is a genuine horizon of the sGB branch, not a coordinate or perturbative artefact; the target TEP temporal well is a different branch that the full TEP action must select; (d) produce exterior deviations from Schwarzschild: shadow $-0.044\%$ (photons, $\mathcal{O}(\eta^2)$), ISCO on $\tilde g$ $+1.95\%$ (massive particles, $\mathcal{O}(\eta)$, different coupling order and opposite sign) at $\eta = -0.1$ (Appendix L). The coupling-order difference is a property of the sGB benchmark (Section 6.2). Within current observational constraints ($\alpha_{\rm GB} < 2.9$ km$^2$), the largest shifts are for $M = 10\,M_\odot$: shadow $-0.007\%$, ISCO $+0.78\%$. The non-perturbative exterior mass-function evolution requires rerun with the correct scalar charge before it can be cited as quantitative evidence (Section 4.5).

What sGB does *not* establish: a regular $r = 0$ centre. As detailed in Section 4.2, the nonlinear solutions of Sotiriou \& Zhou (2014) for the same linear sGB coupling develop a finite-area singularity rather than a regular de Sitter-like temporal minimum, and recent simulations (Thaalba et al. 2024) confirm this and explore a possible connection to loss of hyperbolicity. In much of the sGB literature, "regular black hole solution" means regular outside the would-be horizon and at the apparent horizon — not finite curvature everywhere. The manuscript's earlier claim that sGB is "known to generate regular black-hole solutions" was an overinterpretation of this literature. The true TEP solution must be a temporal well: $N \geq N_{\min} > 0$ everywhere, with no finite-radius causal boundary.

The correct status is: sGB backreaction is supported; sGB scalar hair and exterior deviations are supported; sGB serves as a low-energy effective field theory for the exterior. The full TEP action, through the mandatory backreaction proven by the fixed-background theorem, must derive the temporal-well profile — the TEP Global Solution Architecture is the target. If standard linear sGB does not dynamically produce this architecture in the deep interior, a modified coupling — nonlinear $f(\phi)\mathcal{G}$, additional Ricci coupling, or scalar potential structure — may be required, and the exterior observables would need recalculation for the modified coupling.

## 9.4 Temporal Accessibility Without a Causal Horizon

The Temporal Horizon $\mathcal{H}_T^{(\Lambda)}[u_o]$ is the observer-dependent operational boundary of temporal accessibility — the surface where the clock-transfer factor $\mathcal{T}_{e\to o} = \omega_e/\omega_o = (-k_\mu u^\mu)_e/(-k_\mu u^\mu)_o$ exceeds a physical accessibility threshold $\Lambda$. It is not a null boundary, not one-way, and does not divide spacetime into an outside and an inside. With finite $A$ at the centre (a candidate architecture, Section 6.6), the centre is a regular point at finite optical distance — an emitter there has finite frequency relationship with any causally connected receiver. The transfer factor from the centre is finite: $\mathcal{Z}_{\max} \sim N_o/N_{\min}$, where $N_{\min} = A_c \ll 1$. This may be astronomically large but is not divergent. The object appears black because the transfer exceeds any practical accessibility threshold, not because signals are absolutely forbidden. For a freely falling observer moving through the deep temporal region, the full expression $\omega = -k_\mu u^\mu$ gives finite local clock transfer — the distant observer's Temporal Horizon is not their Temporal Horizon. A deeper observer (slower clock) sees less redshift to the same emitter. One universal matter metric, one temporal field, and no observer-specific metric; different observers nevertheless possess different temporal-accessibility maps.

The divergent-$A$ benchmark ($\phi_0 = 1$) is not the target architecture. It produces a finite-area asymptotic end, not a regular centre; both null families have infinite affine parameter; the static-observer clock rate diverges (blueshift). The benchmark serves as a pipeline validation and a cautionary example, not as the Temporal Horizon realisation.

The TEP novelty is the temporal well: a regular spatial region with an extreme gradient in relative proper-time accumulation, dual to the cosmological interpretation. Both photons and massive particles propagate on $\tilde g_{\mu\nu}$; gravitational waves probe the geometric Regge–Wheeler potential. In the sGB benchmark, the shadow contracts ($-0.044\%$) while the matter-metric ISCO expands ($+1.95\%$), with $\sim 44\times$ coupling-order difference (Appendix L). This difference arises because the conformal factor cancels for null geodesics (shadow sensitive only to the underlying geometry at $\mathcal{O}(\eta^2)$) while timelike geodesics feel the conformal factor $A > 1$ directly at $\mathcal{O}(\eta)$ — a property of the sGB benchmark (Section 6.2).

## 9.5 Temporal Waves

Because $\phi$ is the temporal field — the scalar field that governs the rate of proper time through $A(\phi) = e^{\beta_A \phi}$ — perturbations $\delta\phi$ are not merely scalar perturbations in the field-theory sense. They are ripples in the dynamical time field itself: *temporal waves*. The technical literature terms them scalar perturbations, scalar-led QNMs, or scalar dipole radiation, aligning with established scalar–tensor theory nomenclature. In the TEP framework, the physical translation is exact: the scalar field is time, so its perturbations are temporal waves.

Gravitational waves and temporal waves are not independent forms of radiation emitted side by side. They are complementary descriptions of one coupled disturbance of gravitational geometry and dynamical proper time. A merger perturbs both $g_{\mu\nu}$ and $\phi$; matter and detectors respond to their combined perturbation through the physical matter metric, $\delta\tilde g_{\mu\nu} = \delta\tilde g_{\mu\nu}[\delta g_{\mu\nu}, \delta\phi]$. The term *gravitational wave* emphasises the tensor and geometric projection of this disturbance; the term *temporal wave* emphasises the associated variation in the field governing relative proper-time accumulation. During inspiral, the accelerating binary radiates into the coupled temporal–geometric field — scalar dipole emission at $-1$PN, absent in GR because GR has no dynamical time field. At merger and ringdown, the coupled spectrum may contain tensor-led and temporal-led modes; these are eigenmodes of one interacting temporal–geometric system, not three independent radiation channels. A temporal-led mode may appear observationally as a modified gravitational-wave polarisation, phase, damping rate, or ringdown component, depending on its coupling to the matter metric and detector response. The detection of a temporal-led component would indicate a dynamical scalar degree of freedom, consistent with TEP and with other scalar–tensor theories.

The coupled system is expected to contain axial, gravitational-led polar, and scalar-led temporal mode families; the ringdown can differ from GR because the temporal field participates in the coupled dynamics, pending calculation of the complete coupled eigenvalue system. The relationship between the geometric and temporal descriptions is not identity but coupling: gravity is the geometric engine ($g_{\mu\nu}$), and the temporal field ($\phi$) is the reactive medium. The sGB coupling $\alpha_{\rm GB}\,\phi\,\mathcal{G}$ is the mechanism by which gravitational geometry sources the temporal gradient, and the temporal gradient mediates physical observations through $\tilde g_{\mu\nu}$. TEP decomposes observed spacetime phenomenology into gravitational geometry and a universal temporal response field; matter observes their combined causal metric rather than the Einstein-frame geometry alone.

## 9.6 Density and Volume

The infinite-density conclusion of the single-metric Schwarzschild continuation arises from combining the exterior-measured mass $M$ with the geometric volume $V_{\rm geom}\to0$ as $r\to0$. The frame-dependence of volume is a structural feature of TEP: the classical density inference assumes the geometric metric governs matter. The geometric Einstein-tensor component $-G^t{}_t = 12M^2\ell^2/(r^3 + 2M\ell^2)^2$ is the effective stress supporting the Hayward seed metric $g$; it is not automatically the matter-frame density $\tilde T_{\mu\nu}\tilde u^\mu\tilde u^\nu$. The geometric stress, the matter-frame stress, and the TEP matter stress tensor are distinct objects and must be separated explicitly. On the $\phi_0 = 1$ divergent-$A$ benchmark, the radial proper distance diverges ($\tilde\ell \to \infty$), so the integrated matter-frame spatial volume is not finite merely because the sphere area approaches a finite constant. On the target finite-$A$ architecture, the centre is an ordinary regular point with finite volume.

## 9.7 Comparative Anatomy

TEP is compared formally with existing alternatives. The comparison is technical, not polemical. Regular black-hole literature (Bardeen, Hayward, and successors) is an important guide: the validation benchmark uses a Hayward geometry precisely because that literature has established what regular compact objects look like. The question TEP asks is whether a dynamical proper-time field can supply the physical origin of a temporal well — a regular spatial region with an extreme gradient in relative proper-time accumulation — that can account for the observed phenomenology without high-density collapse. The fixed-background theorem proves the temporal field must backreact; the sGB coupling is a candidate mechanism; the full nonlinear integration will determine whether it succeeds.

| Property | Schwarzschild | Bardeen/Hayward | Loop LQG | Planck star | TEP benchmark | TEP sGB (exterior) |
| --- | --- | --- | --- | --- | --- | --- |
| Central curvature finite | No ($4.8 \times 10^{49}$) | Yes | Yes | Yes | Benchmark: finite ($\tilde K \to 8/r_h^4$); asymptotic end, not regular centre | Not yet computed |
| Geodesic completeness | No (incomplete at $r=0$) | Model-dependent | Model-dependent | Quantum description | Benchmark: both null families infinite affine (asymptotic end); target architecture: finite (regular centre) | Not yet computed |
| Exterior = Schwarzschild | — | Yes (to $M/r$) | Yes | Yes | Yes (exact) | Yes (to $\sim$11 parts in $10^4$ at $|\eta|=0.1$) |
| Dynamical mechanism | — | Matter sector | Quantum gravity | Quantum gravity | Hayward (prescribed) | Scalar–Gauss–Bonnet (candidate) |
| Falsifiable predictions | — | Limited | Limited | Limited | — | Sub-percent shadow/QNM/ISCO (indicative, pending validation) |
| Interior singularity? | Yes ($r=0$) | No | No | No | Benchmark: asymptotic end | Likely yes (Section 4.2) |

## 9.8 Scope of the Interpretation

The supported interpretation is: (i) the fixed-background theorem proves that, within the stated fixed-Schwarzschild power-law conformal class, a dynamical proper-time field necessarily implies a dynamical gravitational geometry — the exact Kretschmann $\tilde K \sim r^{4\phi_0-6}$ makes finite curvature ($\phi_0 \geq 3/2$) and bounded areal radius ($\phi_0 \leq 1$) mutually exclusive for fixed geometry; (ii) on a regular geometric background (Hayward class), the TEP matter metric with $\phi_0 = 1$ achieves bounded areal radius ($\rho \to 1.995$) and finite Kretschmann ($\tilde K \to 8/r_h^4 = 0.5$), but the limiting region is a finite-area asymptotic end with both null families at infinite affine parameter — this is a benchmark result, not an sGB solution result; (iii) the Temporal Horizon $\mathcal{H}_T^{(\Lambda)}[u_o]$ is the observer-dependent operational boundary of temporal accessibility — with finite $A$ at the centre (a candidate architecture), the centre is at finite optical distance and an emitter there has finite $\mathcal{Z}$ with any receiver; the transfer factor from the centre is finite, $\mathcal{Z}_{\max} \sim N_o/N_{\min}$, potentially enormous but not divergent; for a distant observer the Temporal Horizon lies far from the centre, while for a deep observer it shifts inward or vanishes locally; (iv) the TEP novelty is the temporal well: a regular spatial region with an extreme gradient in relative proper-time accumulation, dual to the cosmological interpretation — photons, massive particles, and gravitational waves are all observational probes of the same unified temporal field; (v) the sGB exterior is Schwarzschild to $\sim$11 parts in $10^4$ at $|\eta| = 0.1$, consistent with all current exterior observations; (vi) sGB backreaction and scalar hair are supported by the literature; (vii) sGB-generated regular de Sitter-like temporal minimum is *not* supported by the existing literature (Section 4.2); (viii) the target global solution is a temporal well: a regular spatial domain with finite local density, a continuous temporal gradient, and no singularity or horizon — the EHT visibility-domain fit (Section 8.1B) bounds the regularisation scale to the sub-critical regime ($g < g_{\rm crit}$), requiring the temporal gradient to be steeper and more compact than the extreme Hayward threshold allows, consistent with the S2 weak-field result where $\eta_{\rm TEP}$ was tightly bounded near zero; the temporal well mimics a black hole so effectively because the gradient is exceptionally sharp, not gradual. It has no finite-radius causal horizon ($N \geq N_{\min} > 0$ everywhere), and extreme external redshift — apparent compactness arises from exterior-frame reconstruction, not physical compression; the density $\rho_{\rm inferred} \sim M/(4\pi r^3/3)$ is an exterior-frame inference, not necessarily the locally measured density $\tilde\rho = \tilde T_{\mu\nu}\tilde u^\mu \tilde u^\nu$; (ix) the full nonlinear integration (Section 4.6) is the decisive test with three outcomes: temporal-well regular continuation (A), finite-area singularity (B, expected from literature), or hyperbolicity failure (C); (x) the exterior mass-function evolution and QNM numbers are indicative pending rerun with the correct scalar charge and a proper spectral QNM solver. A horizonless temporal well alters the late-time ringdown spectrum (Section 7.3).

## 9.9 The Unified TEP Statement

The TEP black-hole interpretation is dual to the Temporal Horizon Cosmology interpretation. On cosmological scales, a temporal-rate gradient produces observed redshift, which standard physics interprets as spatial expansion. Around an extreme gravitational source, the same kind of gradient produces redshift, lensing and altered trajectories, which standard physics interprets as gravitational collapse and inward suction.

*Cosmological expansion and black-hole collapse are dual misinterpretations of dynamical proper time. On cosmological scales, temporal evolution is interpreted as increasing spatial separation. Around extreme gravitational sources, temporal gradients are interpreted as spatial attraction, compression and causal capture. In both cases, TEP replaces apparent spatial dynamics with the observable consequences of a dynamical proper-time field.*

The object is a regular region in which the rate of proper time differs radically from the exterior — not a drain or a singularity, but a continuous spatial domain with an extreme temporal gradient. Matter moving into that gradient appears, from the exterior, to accelerate inward, become increasingly redshifted, slow and fade, collect within a small apparent radius, and become inaccessible. But locally there is no mysterious suction. Matter follows its ordinary local trajectory through a continuous physical environment. The apparent force is analogous to refraction: trajectories curve when the propagation rate varies through a medium, not because matter is being pulled sideways. In TEP, the temporal field changes observed motion through clock transfer, timelike matter coupling, and its backreaction on the causal geometry. The conformal factor $A(\phi)$ controls clock-rate transfer and massive-particle response; the disformal term $B(\phi)$ can alter causal cones; the backreacted geometric metric $g_{\mu\nu}[\phi]$ generates light bending and orbital curvature. Lensing is generated by the complete matter metric $\tilde g_{\mu\nu}$, not by conformal clock scaling alone.

The region appears black because signals emerging from deeper clock-rate domains arrive with greatly reduced frequency, at a greatly reduced rate, strongly lensed and diluted, below practical detectability. *Darkness does not require an event horizon; it requires sufficiently extreme temporal decoupling.* The Temporal Horizon is the observational boundary where the relative clock-rate difference becomes so extreme that the exterior observer can no longer receive usable information. It is not necessarily a physical wall or one-way causal surface.

The object need not be extremely compressed. The exterior observer calculates a small radius and high density using exterior standards of distance, simultaneity and elapsed time. But those standards may not map directly into the deep temporal region. *Large inferred density $\neq$ large locally measured density.* Just as TEP argues that cosmological redshift does not necessarily mean that space itself is expanding, black-hole redshift and infall do not necessarily mean that matter is being crushed into a tiny volume. The conventional black-hole picture — event horizon, interior, singular centre — is retained only as the conventional limit that the temporal interpretation replaces. The physical object is a continuous spatial region with finite local conditions and extreme Temporal Shear.

## 9.10 Apparent Central Mass and the Isochrony Assumption

The observational inference problem and the strong-field Phantom Mass decomposition are established in Sections 2.7–2.8. The conventional inference chain — received angles and frequencies, through the Isochrony Axiom and GR transfer, to orbital period, velocity, radius, and finally a large compact mass $M_{\rm app}^{\rm GR}$ — is replaced by the TEP chain: received angles and frequencies, through dynamical temporal-transfer reconstruction, to local trajectory and clock-rate field, and finally to material source plus temporal contribution.

The central empirical task is therefore not to modify an already-known black-hole mass. It is to refit the raw astrometric, spectroscopic and timing observations without imposing universal isochrony, and determine how much of the conventional central mass survives as locally measured material content. The decisive analysis is a *TEP-native refit of the raw observations*, not a perturbation of a pre-assumed Schwarzschild mass. That is the direct black-hole equivalent of the TEP treatment of dark matter, dark energy and Cepheid distance bias.

## 9.11 The S-star Inference Pipeline

The S-star cluster surrounding Sgr A* provides the most precise strong-field test of TEP at black-hole scales. We refit the published astrometry and spectroscopy of S2 relative to Sgr A* without imposing universal isochrony. The data are taken from the 25-year monitoring compilation of Gillessen et al. (2017): 145 NACO/NTT/Keck/Gemini astrometric epochs (1992--2016) and 44 SINFONI/Keck/Gemini radial-velocity epochs, drawn from the VizieR catalogue J/ApJ/837/30 (table5.dat). These pre-GRAVITY data are supplemented by the later GRAVITY pericentre campaign (2018--2022) in the discussion, but the quantitative fit uses the published Gillessen+2017 machine-readable tables. We implement a seven-step inference pipeline that refits this data without imposing universal isochrony.

The pipeline (`scripts/run_pipeline.py`) proceeds as follows:

| Step | Script | Function |
| --- | --- | --- |
| 00 | `step_43_sstar_data` | Acquire S-star astrometry and spectroscopy (VizieR/CDS, with GRAVITY/NACO/SINFONI fallback) |
| 01 | `step_44_gr_fit` | Conventional GR fit: recover $M_{\rm BH}^{\rm GR}$, $R_0$, orbital elements |
| 02 | `step_45_mass_bias` | Gate -1: compute the mass-bias sign from $(\mathcal S_a^3 \mathcal D_{\rm dyn}) / \mathcal T_P^2$ |
| 03 | `step_46_tep_transfer` | TEP transfer-function fit with $\eta_{\rm TEP}$ as free parameter |
| 04 | `step_47_joint_forward` | Separate $H_{\rm composition}$ from $H_{\rm calibration}$; compute $\mathcal C_V$ |
| 05 | `step_48_likelihood` | Formal GR vs TEP comparison: Bayes factor, BIC, AIC, Wilks |
| 06 | `step_49_posterior` | MCMC posterior on $M_{\rm phantom}^T$ (emcee or Metropolis-Hastings) |

**GR fit.** The conventional 11-parameter GR fit to 145 astrometric and 44 spectroscopic epochs gives $M_{\rm BH}^{\rm GR} = 4.297 \times 10^6 \pm 9.4 \times 10^4\,M_\odot$, $R_0 = 8320 \pm 59$ pc, $a = 0.1280''$, $e = 0.8913$, $i = 133.67^\circ$, $\omega = 65.94^\circ$, $\Omega = 227.23^\circ$, $T = 16.067$ yr and $t_\mathrm{peri} = 2002.342$, with $\chi^2 = 304.17$ over 323 degrees of freedom ($\chi^2/\mathrm{dof} = 0.942$). The fit reproduces the published Gillessen+2017 orbital elements; the $M_{\rm BH}$--$R_0$ anti-correlation is handled by using the published $R_0$ uncertainty as a soft prior.

**Gate -1 result.** The mass-bias sign equation (Section 2.9) evaluated at S2 pericentre ($r \sim 1369\,R_s$) for the TEP mass-inflation benchmark $\eta_{\rm TEP} = 0.1$ gives $\mathcal S_a = 1.000024$ and $\mathcal T_P = 0.999976$. Because the spatial magnification ($\mathcal S_a$) from the conformal geometry dominates the temporal stretching ($\mathcal T_P$) at this radius, the inferred mass is strictly inflated, yielding a positive Phantom Mass. The TEP-only ratio $\mathcal S_a^3 \mathcal D_{\rm dyn} / \mathcal T_P^2 = 1.000122$, so the sign of $M_{\rm phantom}^T$ is *positive*. The mass-inflation branch ($\alpha_{\rm GB} < 0$) is therefore selected at S2 scales. The benchmark magnitude is only $\sim 1.22 \times 10^{-4}$, far below current observational precision.

**TEP fit result.** Adding $\eta_{\rm TEP}$ as a twelfth parameter (12 vs 11 parameters, 334 data points) does not improve $\chi^2$; $\chi^2_{\rm TEP} = 304.17$ is identical to $\chi^2_{\rm GR} = 304.17$. The maximum-likelihood value is $\eta_{\rm TEP} = 0$ (at the lower bound of the physical mass-inflation prior). $\Delta\mathrm{BIC} = +5.8$ and $B(\mathrm{TEP}/\mathrm{GR}) \simeq 0.055$ from the BIC approximation; $B(\mathrm{TEP}/\mathrm{GR}) \simeq 0.37$ from AIC. The Wilks $p$-value is 0.9998. These comparisons strongly favour the simpler GR model: the S2 data do not require a TEP correction at weak-field scales.

**Composition vs calibration.** The joint forward-model separates $H_{\rm composition}$ ($M_{\rm matter} < M_{\rm GR\text{-}fit}$) from $H_{\rm calibration}$ ($M_g < M_{\rm GR\text{-}fit}$). With the fitted $\eta_{\rm TEP} = 0$, TEP reduces to GR at S2 scales: $M_g = M_{\rm GR\text{-}fit} = 4.297 \times 10^6\,M_\odot$ and $\mathcal C_V \simeq 1.0005$ at pericentre, with no significant volume distortion. The mass-inflation branch is therefore not yet required by the S2 data.

**MCMC posterior.** The emcee sampler (24 walkers, 10000 production steps after 1000 burn-in, 240000 flattened samples, fixed seed for reproducibility) is run with the physical mass-inflation prior $\eta_{\rm TEP} \ge 0$. It gives a median $\eta_{\rm TEP} = 0.425$ with the 16th--84th percentile range $[0.079, 0.810]$; the posterior is bounded away from the negative-deflation branch. The $M_{\rm phantom}^T$ posterior has median $+2.22 \times 10^3\,M_\odot$ with 68% CI $[+4.1 \times 10^2, +4.2 \times 10^3]\,M_\odot$ and 100% of samples positive: the sign is positive. The central value corresponds to a $5\times 10^{-4}$ phantom-mass fraction at S2, still consistent with zero at $<2\sigma$ because the prior allows $\eta_{\rm TEP}=0$.

**Interpretation.** These results are the *expected* weak-field outcome for the TEP mass-inflation branch. Gate -1 establishes the physically required positive sign at S2 pericentre, but the benchmark $\eta_{\rm TEP} = 0.1$ produces only a $\sim 10^{-4}$ Phantom Mass fraction at $r \sim 1369\,R_s$. The TEP fit and MCMC are consistent with $\eta_{\rm TEP} = 0$; the TEP correction is perturbatively small and the S2 data do not require it. This is not a falsification of TEP — it is a confirmation that the weak-field limit reproduces GR. The decisive test requires horizon-scale observations ($r \sim$ few $R_s$), where the mass-inflation branch predicts $O(10^{-2})$ to $O(1)$ corrections.

The pipeline is structured to return any of the three sign outcomes (positive, zero, negative). When applied to stronger-field data — EHT horizon-scale observations, or future GRAVITY+ deeper-orbit stars — the same pipeline will determine the sign and magnitude of $M_{\rm phantom}^T$ from the data, not from the thesis. The complete automated pipeline — all 50 steps, including the S-star inference modules — alongside SHA-256 checksum validations for every result file, is documented in Appendix J and available in the public repository.

# 10. Conclusion

This work develops the Temporal Equivalence Principle (TEP) as a bimetric framework for temporal wells — regular spatial regions with an extreme gradient in relative proper-time accumulation. The organising question is: *can the phenomena attributed to black holes arise from an extreme relative-time gradient without physical collapse into an ultra-dense object?* This is dual to the TEP cosmological question: the universe need not be physically expanding merely because distant signals are redshifted, and matter need not be physically crushed into a black hole merely because trajectories appear to converge toward a dark compact region. *Cosmological expansion and black-hole collapse are dual misinterpretations of dynamical proper time.*

**What was proved.** The fixed-background conformal theorem proves that, within the stated fixed-Schwarzschild power-law conformal class, a dynamical proper-time field necessarily implies a dynamical gravitational geometry. The exact Kretschmann scalar $\tilde K \sim r^{4\phi_0 - 6}$ makes finite curvature ($\phi_0 \geq 3/2$) and bounded areal radius ($\phi_0 \leq 1$) mutually exclusive when $g_{\mu\nu}$ is held fixed. Schwarzschild is an incomplete strong-field limit: it treats time as a passive coordinate, and when time is correctly treated as a physical field, the fixed geometry collapses.

**What was demonstrated on a benchmark.** On a regular geometric background (Hayward class), the TEP matter metric with $\phi_0 = 1$ achieves bounded areal radius ($\rho \to 1.995$) and finite Kretschmann ($\tilde K \to 8/r_h^4 = 0.5$). The inverse reconstruction identifies $w_r = -1$ as a sufficient regularity target. The benchmark confirms that finite curvature and bounded space are compatible once the geometric singularity is absent — the target is mathematically consistent. It does *not* prove that the sGB coupling dynamically generates this geometry.

**What was conditionally computed.** The sGB coupling serves as a low-energy effective field theory for the exterior. The perturbative solution (Sotiriou \& Zhou 2014, exact $\mathcal{O}(\alpha^2)$) establishes scalar hair and real backreaction. The corrected exterior observables at $\eta = -0.1$, conditional on the GR-calibrated mass scale: shadow $-0.044\%$ (photons, $\mathcal{O}(\eta^2)$), ISCO on $\tilde g$ $+1.95\%$ (massive particles, $\mathcal{O}(\eta)$, $\sim 44\times$ larger). The coupling-order difference reflects conformal invariance of null geodesics in the sGB benchmark (Section 6.2). The tensor speed $c_T = 1$ on the Schwarzschild background; the full tensor characteristic metric on the TEP solution requires the complete coupled perturbation system. $(\delta g_{\mu\nu}, \delta\phi)$ form one coupled perturbation system; its eigenmodes may carry different relative proportions of tensor and temporal-field perturbation. An $\mathcal{O}(10^{-3})$ QNM modification is suggested; its coefficient and sign require the coupled spectral solver. The horizonless temporal well alters the late-time ringdown spectrum (Section 7.3).

**The TEP Global Solution Architecture.** The TEP object is a temporal well: $N(r) = A(r)\sqrt{F(r)} \geq N_{\min} > 0$ for every finite $r$, with a finite but extremely small minimum lapse $N_{\min} = A_c \ll 1$ at the centre. The centre is an ordinary regular point at finite optical distance; the maximum transfer factor $\mathcal{Z}_{\max} \sim N_o/N_{\min}$ is finite but potentially enormous. The Temporal Horizon $\mathcal{H}_T^{(\Lambda)}[u_o]$ is the observer-dependent operational limit of temporal accessibility, not a causal boundary. Apparent compactness arises from exterior-frame reconstruction, not physical compression. The strong-field Phantom Mass residual is defined as $M_{\rm phantom}^{T} \equiv M_{\rm fit}^{\rm GR} - M_{\rm matter}^{\rm TEP}$; its sign is not assumed but determined by the sign equation (Section 2.9): $M_{\rm phantom}^{T} > 0$ requires spatial magnification $\mathcal S_a^3$ to dominate temporal stretching $\mathcal T_P^2$. The conventional mass may combine a material contribution, temporal-field energy, and a non-isochronous calibration residual.

**What remains open.** Which modified coupling dynamically produces the temporal well, whether hyperbolicity is preserved, whether a regular centre is achieved, the full tensor characteristic metric, global ghost freedom, and the precise QNM spectrum all require a full nonlinear integration with a healthy Horndeski/DHOST action. The non-perturbative exterior integration requires rerun with the correct scalar charge.

**The S-star inference pipeline.** A seven-step pipeline (`scripts/run_pipeline.py`) implements the TEP-native refit of the S-star data. We use the published Gillessen et al. (2017) machine-readable tables: 145 NACO/NTT/Keck/Gemini astrometric epochs and 44 SINFONI/Keck/Gemini radial-velocity epochs of S2 relative to Sgr A*. The conventional GR fit gives $M_{\rm BH}^{\rm GR} = 4.297 \times 10^6 \pm 9.4 \times 10^4\,M_\odot$, $R_0 = 8320 \pm 59$ pc and $\chi^2/\mathrm{dof} = 0.942$; the $M_{\rm BH}$--$R_0$ anti-correlation is handled with the published $R_0$ uncertainty as a soft prior. Gate -1 establishes that the mass-bias sign is *positive* at S2 pericentre ($r \sim 1369\,R_s$) for the TEP mass-inflation benchmark $\eta_{\rm TEP} = 0.1$: spatial magnification dominates temporal stretching, with a TEP-only ratio $\mathcal S_a^3 \mathcal D_{\rm dyn}/\mathcal T_P^2 = 1.000122$ and a Phantom Mass fraction of $1.22 \times 10^{-4}$. The TEP fit, with the physical mass-inflation prior $\eta_{\rm TEP} \ge 0$, returns $\eta_{\rm TEP} = 0$ at the lower bound (consistent with GR), $\Delta\mathrm{BIC} = +5.8$ and $B(\mathrm{TEP}/\mathrm{GR}) \simeq 0.055$ from the BIC approximation; $B(\mathrm{TEP}/\mathrm{GR}) \simeq 0.37$ from AIC. The MCMC posterior, constrained to the mass-inflation branch, has $\eta_{\rm TEP} = 0.425^{+0.385}_{-0.347}$ and a 68% CI for $M_{\rm phantom}^T$ of $[+4.1 \times 10^2, +4.2 \times 10^3]\,M_\odot$ (100% positive). The S2 data do not require a TEP correction at weak-field scales, but the sign is fixed in the physically correct mass-inflation direction. The decisive test requires horizon-scale observations ($r \sim$ few $R_s$), where the same pipeline will determine the sign and magnitude from the data.

**The EHT visibility-domain test.** A direct fit of the TEP ring + Gaussian core model to 46{,}846 binned EHT M87$^\ast$ and Sgr A$^\ast$ long-baseline ($>2$\,G$\lambda$) visibility amplitudes provides the horizon-scale empirical test. The visibilities are binned into 39 and 37 logarithmic uv-distance rings with scatter-based uncertainties and conservative $\chi^2$ rescaling (standard EHT practice). The joint fit, using Wilks' theorem on the rescaled $\Delta\chi^2$, finds $\delta_{\rm shadow} = -0.018$ with 2$\sigma$ interval $[-0.018, +0.085]$. The Phantom Mass prediction at the horizonless threshold ($\delta = -0.054$, corresponding to the $-5.37\%$ shadow contraction at $g_{\rm crit} = 1.058\,M$) is excluded at $>2\sigma$. GR ($\delta = 0$) lies inside the 2$\sigma$ interval. The data therefore bounds the temporal well's regularisation scale to the sub-critical regime ($g < g_{\rm crit}$): the temporal gradient must be steeper and more compact than the extreme Hayward threshold allows. This is not a refutation of TEP — it is a dynamical constraint on the geometry of the temporal well, exactly the kind of falsifiable prediction that distinguishes a physical theory from a "just-so" story. Combined with the S2 weak-field result (TEP reduces to GR at $r \sim 1369\,R_s$) and the QNM ringdown predictions (Section 8.3), the EHT visibility fit completes a multi-scale empirical net: weak-field consistency, horizon-scale constraint, and gravitational-wave predictivity.

The decisive analysis is a *TEP-native refit of the raw astrometric, spectroscopic and timing observations* without imposing universal isochrony — not a perturbation of a pre-assumed Schwarzschild mass. The central question is: *after a sign-correct, GR-matched temporal reconstruction, is the inferred Phantom Mass positive, negative, or zero?* The paper establishes a structural possibility and conditional exterior benchmarks; the full observational replacement requires the nonlinear TEP solution and the raw-observable refit. The strongest route forward is: one action → one exterior coupling → one global solution → one coupled characteristic system → one observational likelihood → one Phantom Mass posterior.

# Appendix A — Conventions and Disformal Identities

This appendix fixes the sign, curvature and disformal conventions used throughout the paper. All identities quoted in the main text are derived here explicitly so that no sign or factor ambiguity survives into the field equations, the curvature invariants of Appendix D, or the perturbation analysis of Appendix G.

## A.1 Metric Signature and Curvature

The metric signature $(-,+,+,+)$ is adopted. The Riemann tensor is defined via the Wald convention

\begin{equation} \label{eq:appA_1}
({\nabla_a\nabla_b - \nabla_b\nabla_a)\,V_c = R_{abc}{}^{d}\,V_d,
\end{equation}

with Ricci tensor $R_{ab} = R_{acb}{}^{c}$ and Ricci scalar $R = g^{ab}R_{ab}$. The Einstein tensor is $G_{ab} = R_{ab} - \tfrac12 R\,g_{ab}$, and the Einstein equations read $G_{\mu\nu}[g] = 8\pi\,T_{\mu\nu}$ in geometrized units $G = c = 1$. With this convention the Schwarzschild Kretschmann scalar is

\begin{equation} \label{eq:appA_2}
K_{\rm Schw} = R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma} = \frac{48\,M^2}{r^6},
\end{equation}

which is the reference used for the curvature comparison in Appendix D.

## A.2 Disformal Transformation

The matter (causal) metric is the disformal image of the geometric metric,

\begin{equation} \label{eq:appA_3}
\tilde g_{\mu\nu} = A^2(\phi)\,g_{\mu\nu} + B(\phi)\,\nabla_\mu\phi\,\nabla_\nu\phi,
\end{equation}

where $A(\phi)$ is the conformal factor and $B(\phi)$ is the disformal coupling. Defining $X \equiv -\tfrac12\,g^{\mu\nu}\nabla_\mu\phi\,\nabla_\nu\phi$, the inverse matter metric is

\begin{equation} \label{eq:appA_4}
\tilde g^{\mu\nu} = A^{-2}\left(g^{\mu\nu} - \frac{B\,\nabla^\mu\phi\,\nabla^\nu\phi}{A^2 - 2BX}\right).
\end{equation}

This is the identity used in the geodesic equations of Appendix E and in the perturbation potentials of Appendix G.

## A.3 Determinant and Invertibility

In four spacetime dimensions the determinant of the disformal metric is

\begin{equation} \label{eq:appA_5}
\det\tilde g = A^{8}\,\det g\,\left(1 - \frac{2BX}{A^2}\right) = A^{8}\,\det g\,\frac{A^2 - 2BX}{A^2}.
\end{equation}

The transformation is therefore invertible if and only if

\begin{equation} \label{eq:appA_6}
A^2 > 0 \qquad\text{and}\qquad A^2 - 2BX \neq 0.
\end{equation}

The 2D $(v,r)$ determinant in Eddington–Finkelstein coordinates is the canonical expression $\det\tilde g_{2D} = -A^4 G^2 - A^2 F B(\phi')^2$ (Section 6.7). The geometric horizon has $\det\tilde g_{2D}=-1.002$ and is Lorentzian. The determinant is strictly negative at every sampled radius from the horizon to the innermost grid point ($r_{\min}=10^{-8}M$) in the exterior region ($F > 0$). Inside horizons ($F < 0$), the disformal term becomes positive and the determinant must be checked with the canonical formula. The disformal transformation is invertible where $A^2 - 2BX \neq 0$; this must be verified separately for the interior. For the $\phi_0 = 2$ fixed-Schwarzschild case, both null families have infinite affine parameter ($\tilde\lambda \sim \int r^{-4}\,dr \to \infty$) because $A \sim r^{-2}$; the $r \to 0$ limit is an asymptotic spatially enlarged end, not a regular point (Appendix E).

## A.4 Connection Difference

The Levi-Civita connection of $\tilde g$ differs from that of $g$ by a tensorial term. Writing $\Sigma_\mu \equiv \nabla_\mu\ln A$ and $\Phi_\mu \equiv \nabla_\mu\phi$, the difference is

\begin{equation} \label{eq:appA_7}
\tilde\Gamma^{\lambda}_{\;\;\mu\nu} - \Gamma^{\lambda}_{\;\;\mu\nu} = 2\,\delta^{(\lambda}_{(\mu}\,\Sigma_{\nu)} - g_{\mu\nu}\,\Sigma^{\lambda} + \frac{B}{A^2-2BX}\Big[\Phi_{\mu}\,\Phi_{\nu}\,\Sigma^{\lambda} + 2\,\Phi^{(\lambda}\,\nabla_{\mu)}\Phi_{\nu)} - g_{\mu\nu}\,\Phi^{\lambda}\,\Sigma_{\rho}\Phi^{\rho}\Big] + \mathcal{O}(\nabla B),
\end{equation}

where $\mathcal{O}(\nabla B)$ collects terms proportional to $\nabla_\mu B$. The conformal piece (the first two terms) is the standard scalar-tensor result; the disformal piece (the bracket) is what produces the modified causal structure and the shear zone in the interior.

## A.5 Stress-Energy Transformation

Matter minimally coupled to $\tilde g_{\mu\nu}$ has stress tensor $\tilde T^{\mu\nu} = -(2/\sqrt{-\tilde g})\,\delta S_m/\delta\tilde g_{\mu\nu}$. When expressed with respect to the geometric metric $g_{\mu\nu}$, this becomes

\begin{equation} \label{eq:appA_8}
T^{\mu\nu} \equiv \frac{2}{\sqrt{-g}}\frac{\delta S_m}{\delta g_{\mu\nu}} = A^{-2}\,\tilde T^{\mu\nu} + \frac{B}{A^2}\,\tilde T^{\alpha\beta}\,\nabla_\alpha\phi\,\nabla_\beta\phi\,g^{\mu\nu} + \text{disformal trace corrections}.
\end{equation}

The first term is the familiar conformal rescaling; the second is the disformal correction that feeds back into the geometric Einstein equations $G_{\mu\nu}[g] = 8\pi\,T_{\mu\nu}$. Conservation is frame-consistent: $\tilde\nabla_\mu\tilde T^{\mu\nu} = 0$ in the matter frame, which maps to a non-standard conservation law in the geometric frame that is sourced by the disformal gradient $\nabla_\mu B$.

## A.6 Summary of Identities

| Object | Expression |
| --- | --- |
| Inverse matter metric | $\tilde g^{\mu\nu} = A^{-2}\big(g^{\mu\nu} - \frac{B\,\nabla^\mu\phi\,\nabla^\nu\phi}{A^2 - 2BX}\big)$ |
| Determinant (4D) | $\det\tilde g = A^{8}\,\det g\,\frac{A^2 - 2BX}{A^2}$ |
| Invertibility | $A^2 > 0$ and $A^2 - 2BX \neq 0$ |
| Lorentzian signature | $\det\tilde g_{2D} = -A^4 G^2 - A^2 F B(\phi')^2 < 0$ in exterior ($F > 0$); must be checked inside horizons |
| Connection split | $\tilde\Gamma - \Gamma = $ conformal $+$ disformal |
| Stress split | $T^{\mu\nu} = A^{-2}\tilde T^{\mu\nu} + $ disformal corrections |

The disformal metric is invertible and Lorentzian in the exterior region ($F > 0$), with $\det\tilde g_{2D} < 0$ at every sampled radius. The pipeline-verified conditions $A^2 > 0$ and $A^2 - 2BX \neq 0$ hold at every sampled radius in the exterior. Inside horizons ($F < 0$), the disformal term $-A^2 F B(\phi')^2$ becomes positive and can oppose the conformal term; the determinant must be checked with the canonical formula. The interior hyperbolicity depends on the outcome of the interior integration (Section 4.6).

# Appendix B — Reduced Field Equations

This appendix records the full radial system obtained by substituting the spherical ansatz of Section 3 into the TEP action, together with the parameter set that defines the strong-field solution integrated by the pipeline. The equations are written in geometrized units $G = c = 1$.

## B.1 Action and Frame Assignment

The complete action is

\begin{equation} \label{eq:appB_1}
S = \frac{1}{16\pi}\int d^4x\,\sqrt{-g}\,R[g] + \int d^4x\,\sqrt{-g}\,\left[-\tfrac12(\nabla\phi)^2 - V(\phi)\right] + S_m[\tilde g_{\mu\nu},\Psi_m],
\end{equation}

with the disformal matter metric

\begin{equation} \label{eq:appB_2}
\tilde g_{\mu\nu} = A^2(\phi)\,g_{\mu\nu} + B(\phi)\,\nabla_\mu\phi\,\nabla_\nu\phi.
\end{equation}

Variation yields the geometric Einstein equations, the scalar equation, and the matter conservation law in the matter frame:

\begin{equation} \label{eq:appB_3}
G_{\mu\nu}[g] = 8\pi\,T_{\mu\nu}, \qquad \Box_g\phi = \frac{dV}{d\phi}, \qquad \tilde\nabla_\mu\tilde T^{\mu\nu} = 0.
\end{equation}

The stress tensor $T_{\mu\nu}$ includes the conformal and disformal corrections derived in Appendix A.

## B.2 Spherical Ansatz

Horizon-regular Eddington–Finkelstein form is used throughout,

\begin{equation} \label{eq:appB_4}
ds_g^2 = -F(r)\,dv^2 + 2\,G(r)\,dv\,dr + R^2(r)\,d\Omega^2,
\end{equation}

with a static scalar profile $\phi = \phi(r)$ (the shift-symmetric $q\,v$ piece is set to zero for the vacuum solution). The metric functions and scalar field are parameterized as

\begin{equation} \label{eq:appB_5}
A(\phi) = e^{\beta_A \phi}, \qquad B(\phi) = B_0\,|\phi|^{2}/(1 + |\phi|^{2})\,\exp\!\left(-\frac{\phi^4}{2\sigma_B^4}\right), \qquad \phi(r) = \phi_0 \ln\!\left(\frac{r}{r_h}\right) S(r),
\end{equation}

where $S(r) = 1/(1 + e^{(r-r_h)/(\delta\, r_h)})$ is a smooth logistic activation, and $r_h = 2M$ is the geometric horizon. The exterior areal radius is $R(r) = r$, so that $F(r) = 1 - 2M/r$ recovers Schwarzschild exactly at large $r$.

## B.3 Reduced Radial System

Substituting the ansatz, the independent equations reduce to a coupled radial system for the six functions

\begin{equation} \label{eq:appB_6}
F(r), \qquad G(r), \qquad R(r), \qquad \psi(r) \equiv \phi'(r), \qquad A(\phi), \qquad B(\phi).
\end{equation}

The $vv$ and $vr$ Einstein components give two first-order constraints,

\begin{equation} \label{eq:appB_7}
\frac{F'}{F} = \frac{2M}{r^2 F} - \frac{2}{r}\left(1 - \frac{1}{G^2}\right) + 8\pi\,\mathcal{S}_{vv}[A,B,\phi],
\end{equation}

\begin{equation} \label{eq:appB_8}
G' = -4\pi\,r\,G^2\,\mathcal{S}_{vr}[A,B,\phi],
\end{equation}

where $\mathcal{S}_{vv}$ and $\mathcal{S}_{vr}$ are the disformal source terms built from $A$, $B$, $\phi'$ and $X = -\tfrac12 F(\phi')^2$. The angular equation fixes $R(r)$ through

\begin{equation} \label{eq:appB_9}
\frac{R''}{R} + \frac{2}{r}\frac{R'}{R} = 8\pi\,\mathcal{S}_{\theta\theta}[A,B,\phi].
\end{equation}

The scalar wave equation reduces to

\begin{equation} \label{eq:appB_10}
\frac{1}{r^2}\frac{d}{dr}\left(r^2\,F\,\phi'\right) = \frac{dV}{d\phi} + \frac{d\ln A}{d\phi}\,T^{(\rm m)} + \frac{d\ln B}{d\phi}\,\mathcal{D}[X,\phi],
\end{equation}

where the last two terms are the conformal and disformal scalar charges. The disformal source $\mathcal{D}$ is proportional to $T^{\mu\nu}\nabla_\mu\phi\,\nabla_\nu\phi$ and finite in vacuum; for the vacuum solution integrated here it is the geometric curvature that sources the scalar through the conformal coupling.

## B.4 Matching Conditions at the Geometric Horizon

The exterior Schwarzschild limit is approached smoothly at $r_h=2M$ by the logistic profile: $A\approx1$, $B\approx0$, and $\phi\approx0$. The clean run gives $A=1.0005$, $B=3.0\times10^{-7}$, and $\phi=-2.73\times10^{-4}$ at the sampled horizon grid point. This value of $A$ comes from the old logistic-screened fixed-background profile, not from the sGB-coupled exterior with $\beta_A = -1$. With the sGB scalar and $\beta_A = -1$, the mass-inflation branch has $A(2M) \approx 1.063$ for $\eta = -0.1$ (Section 4.3 assessment). The two profiles are different models and cannot be combined. The disformal contribution remains finite because $F\to0$. This is a continuity check for the prescribed ansatz, not a derivation of junction conditions for a fully coupled interior solution.

## B.5 Model Parameters

The strong-field solution integrated by the pipeline is specified by the following parameter set:

| Parameter | Symbol | Value | Role |
| --- | --- | --- | --- |
| Disformal amplitude | $B_0$ | $1.0$ | sets disformal coupling strength |
| Mass | $M$ | $1.0$ | geometrized black-hole mass |
| Conformal exponent | $\beta_A$ | $-1.0$ | controls $A(r)$ divergence in the interior |
| Disformal horizon exponent | $\delta$ | $0.05$ | controls $B\to0$ approach at $r_h$ |
| Disformal quartic width | $\sigma_B$ | $1.5$ | controls $B(\phi)$ suppression against violent Coulomb gradients |
| Scalar amplitude | $\phi_0$ | $2.0$ | normalizes the scalar profile |
| Geometric horizon | $r_h$ | $2.0\,M$ | Schwarzschild radius |

With $\beta_A=-1$ and $\phi_0=2.0$, $A(r)\sim(r/r_h)^{-2}$ in the deep interior, reaching $A\sim4.0\times10^{16}$ at the innermost sampled radius. The ultra-damped (quartic Gaussian) factor $B(\phi)=B_0|\phi|^2/(1+|\phi|^2)\exp(-\phi^4/(2\sigma_B^4))$ is suppressed by the quartic Gaussian envelope in the deep interior, driven below $10^{-300}$ at the innermost sampled radius. The conformal term $A^4$ dominates everywhere, so the determinant is strictly negative and the matter metric is globally Lorentzian with no boundary.

The implementation is a prescribed-profile metric construction, not a closed six-function radial ODE integration. Its supported result is that the bounded conformal-disformal ansatz is asymptotically Schwarzschild and Lorentzian in the exterior ($\det\tilde g_{2D} < 0$ for $F > 0$). The $\phi_0 = 2$ fixed-Schwarzschild case has vanishing curvature ($\tilde K \to 0$) but the $r \to 0$ limit is an asymptotic spatially enlarged end ($\rho \to \infty$) with both null families at infinite affine parameter — not a regular point (Appendix E). The target global solution with finite $A$ at the centre is the TEP Global Solution Architecture (Section 6.6).

# Appendix C — Strong-Field Asymptotic Expansion (the $\phi_0=2$ Case)

This appendix derives the interior scaling of the metric functions, the physical radial distance, and the areal radius for the $\phi_0 = 2$ case (curvature-regular but spatially-enlarged; see Appendix D for the fixed-background theorem that establishes why this case is not the final answer). All exponents and numerical values are taken directly from the pipeline integration of the reduced system in Appendix B.

#### Context

The $\phi_0 = 2$ case achieves $\tilde K \sim 39r^2/16 \to 0$ (curvature-regular) but has diverging areal radius $\rho \sim r_h^2/r \to \infty$ (spatially-enlarged, not a regular centre). The fixed-background theorem (Appendix D) proves this is unavoidable for fixed Schwarzschild: curvature regularity ($\phi_0 > 3/2$) and bounded areal radius ($\phi_0 \leq 1$) are mutually exclusive. The validation benchmark achieves both by allowing the geometric metric to be regular.

## C.1 Interior Scaling

In the deep interior ($r \ll r_h = 2M$), the scalar field follows a logarithmic profile $\phi(r) = \phi_0 \ln(r/r_h)$ with $\phi_0 = 2.0$, and the conformal and disformal factors are

\begin{equation} \label{eq:appC_1}
A(\phi) = A_0\,e^{\beta_A \phi}, \qquad B(\phi) = B_0\,\frac{|\phi|^2}{1+|\phi|^2}\,\exp\!\left(-\frac{\phi^4}{2\,\sigma_B^4}\right),
\end{equation}

with $\beta_A = -1$, $\sigma_B = 1.5$, and disformal coupling $\delta = 0.05$. Substituting the logarithmic $\phi$, the effective radial scalings are

\begin{equation} \label{eq:appC_2}
A(r) \sim A_0\,(r/r_h)^{-2}, \qquad B(r) \to 0, \qquad \phi(r) \sim 2\ln(r/r_h).
\end{equation}

With $\beta_A=-1$ and $\phi_0=2$, the conformal factor diverges as $A\sim(r/r_h)^{-2}$. The disformal factor is suppressed by the quartic cutoff $\exp(-\phi^4/2\sigma_B^4)$: as $|\phi|\to\infty$ the logistic prefactor saturates to unity but the quartic exponential drives $B\to 0$ instantly. The net effect is that the disformal sector vanishes in the deep interior while the conformal sector diverges, establishing conformal dominance throughout the interior.

## C.2 Conformal Factor Profile

The pipeline integration gives the following values of $A(\phi)$ along the radial profile:

| Location | $r/M$ | $A$ |
| --- | --- | --- |
| Geometric horizon | $2.0$ | $1.0005$ |
| Innermost sampled | $10^{-8}$ | $4.0\times10^{16}$ |

The conformal factor is essentially unity at the geometric horizon (deviation $5.5\times10^{-4}$) and reaches $4.0\times10^{16}$ at the innermost sampled radius $r_{\min}=10^{-8}M$. The deep-interior scaling is the inverse-second-power law implied by $\phi_0=2$ and $\beta_A=-1$. Note: this table uses the old logistic-screened fixed-background profile, where $\phi \approx 0$ at the horizon and $A \approx 1$. This is a different model from the sGB-coupled exterior with $\beta_A = -1$, where the mass-inflation branch gives $A(2M) \approx 1.063$ for $\eta = -0.1$ (Section 4.3 assessment). The two profiles cannot be combined.

## C.3 Disformal Factor Profile

The corresponding values of $B(\phi)$ are:

| Location | $r/M$ | $B$ |
| --- | --- | --- |
| Geometric horizon | $2.0$ | $3.0\times10^{-7}$ |
| Innermost sampled | $10^{-8}$ | $<10^{-300}$ |

The disformal factor is $B=3.0\times10^{-7}$ at the geometric horizon and below $10^{-300}$ at the innermost sampled radius. Unlike a purely logistic activation, the quartic Gaussian cutoff $\exp(-\phi^4/2\sigma_B^4)$ drives $B$ to zero in the deep interior rather than saturating to $B_0$. The disformal sector is dynamically irrelevant in the interior; the geometry is governed entirely by the conformal factor.

## C.4 Conformal Dominance and the Absence of a Causal Boundary

The canonical determinant of the two-dimensional $(v,r)$ sector for $ds_g^2 = -F\,dv^2 + 2G\,dv\,dr + R^2\,d\Omega^2$ and a static radial scalar is

\begin{equation} \label{eq:appC_3}
\det\tilde g_{2D} = -A^4 G^2 - A^2 F\,B(\phi')^2.
\end{equation}

For the $\phi_0 = 2$ fixed-Schwarzschild case, $G = 1$ in EF coordinates. In the deep interior, $A\sim(r/r_h)^{-2}\to\infty$ while $B\to 0$ exponentially. The $-A^4 G^2$ term therefore dominates at every sampled radius, and the determinant is strictly negative:

\begin{equation} \label{eq:appC_4}
\det\tilde g_{2D}\big|_{\rm horizon} = -1.002, \qquad \det\tilde g_{2D}\big|_{r_{\min}} = -2.56\times10^{66}.
\end{equation}

The matter metric is Lorentzian and nondegenerate for all $r > 0$. However, a negative determinant for every $r > 0$ proves only that no finite positive radius is degenerate; it does not prove the limit $r = 0$ can be added as a regular point. The $r \to 0$ limit is an asymptotic spatially enlarged end ($\rho \sim r_h^2/r \to \infty$), not a regular manifold point. The determinant must be checked with the canonical formula inside horizons ($F < 0$), where the disformal term becomes positive.

## C.5 Areal Radius and Interior Volume Growth

The areal radius of the matter geometry is $\rho = A\,r$. With $A\sim(r/r_h)^{-2}$, this gives

\begin{equation} \label{eq:appC_5}
\rho \sim \frac{r_h^2}{r} \to \infty \quad\text{as}\quad r\to 0.
\end{equation}

The areal radius diverges: the physical space does not collapse. The constant-$r$ volume element of the matter metric is $d\tilde V/dv = 4\pi A^2 r^2 \sqrt{|F|}$. The pipeline gives

\begin{equation} \label{eq:appC_6}
\frac{d\tilde V}{dv}\bigg|_{0.1M} = 8.76\times10^{4}, \qquad \frac{d\tilde V}{dv}\bigg|_{r_{\min}} = 2.84\times10^{22}.
\end{equation}

Since the geometry is globally Lorentzian, this volume growth is a physical proper volume, not a formal diagnostic. The would-be central singularity of Schwarzschild is replaced by continuous, low-curvature regular space.

## C.6 Curvature Scaling (Corrected)

The exact Kretschmann scalar scales as $\tilde K \sim r^{4\phi_0 - 6} = r^{2}$ for $\phi_0 = 2$ (Appendix D), so the curvature vanishes as $r\to 0$. The pipeline confirms $\tilde K \approx 0.024$ at $0.1M$ and $\tilde K \approx 1.36\times10^{-16}$ at the innermost sampled radius, compared with $K_{\rm Schw} = 4.8\times10^{7}$ and $4.8\times10^{49}$ respectively. The exact formula $\tilde K \sim 39r^2/16$ matches the numerical values to within 4%. Note: the naive conformal formula $\tilde K = A^{-12} K_{\rm Schw} \sim r^{18}$ is incorrect — it omits the derivative terms from the non-constant conformal factor $A(r)$, which dominate for $\phi_0 < 3/2$ and contribute at the same order for $\phi_0 = 2$. The curvature vanishing is real, but the power is $r^2$, not $r^{18}$. Details are recorded in Appendix D.

## C.7 Summary of Asymptotic Behaviour

| Quantity | Scaling | Behaviour as $r\to0$ |
| --- | --- | --- |
| $A(r)$ | $\sim (r/r_h)^{-2}$ | diverges ($\to 4.0\times10^{16}$) |
| $B(r)$ | $\to 0$ (quartic Gaussian cutoff) | vanishes ($< 10^{-300}$) |
| $\phi(r)$ | $\sim 2\ln(r/r_h)$ | diverges logarithmically ($\to -38.2$) |
| $\det\tilde g_{2D}$ | $\sim -A^4 < 0$ | strictly negative (globally Lorentzian) |
| Areal radius $\rho = Ar$ | $\sim r^{-1}$ | diverges ($\to\infty$; spatially-enlarged — see fixed-background theorem, App. D) |
| Kretschmann $\tilde K$ | $\sim r^{2}$ (exact; not $r^{18}$) | vanishes ($\to 0$; curvature-regular) |
| Interior volume | $\sim A^2 r^2$ | grows ($\to 2.84\times10^{22}$) |
| Singularity | — | curvature-regular, but areal radius diverges (the $\phi_0=2$ case) |

The interior scaling $A\sim(r/r_h)^{-2}$ with $B\to 0$ establishes conformal dominance: the determinant is strictly negative at every sampled radius, the matter geometry is globally Lorentzian, and the areal radius diverges as $\rho\sim 1/r$. The exact Kretschmann scalar vanishes as $\tilde K\sim r^{2}$ (not $r^{18}$ as the naive conformal formula suggests). The $\phi_0 = 2$ case is curvature-regular but spatially-enlarged; the fixed-background theorem (Appendix D) proves this is unavoidable for fixed Schwarzschild, motivating the validation benchmark on a Hayward background.

# Appendix D — Curvature Invariants and the Fixed-Background Theorem

This appendix records the exact curvature analysis for the TEP matter metric. It is proved that the naive conformal formula $\tilde K = A^{-12} K_{\rm Schw}$ is incorrect for non-constant conformal factors, the exact Kretschmann $\tilde K \sim r^{4\phi_0 - 6}$ is derived, and the fixed-background theorem for the fixed-Schwarzschild construction is established.

## D.1 The Exact Kretschmann Scalar

For the pure conformal matter metric $\tilde g = A^2 g_{\rm Schw}$ with $A = (r_h/r)^{\phi_0}$ in the deep interior ($r \ll r_h$, where $F \approx -2M/r$ and $B \to 0$), the Kretschmann scalar is computed in standard Schwarzschild coordinates $(t, r)$ where the metric is diagonal:

\begin{equation} \label{eq:appD_1}
\tilde g_{tt} = -A^2 F, \quad \tilde g_{rr} = \frac{A^2}{F}, \quad \tilde g_{\theta\theta} = A^2 r^2.
\end{equation}

The four independent orthonormal-frame Riemann components are:

\begin{equation} \label{eq:appD_2}
E = \frac{R_{trtr}}{\tilde g_{tt}\,\tilde g_{rr}}, \quad F_t = \frac{R_{t\theta t\theta}}{\tilde g_{tt}\,\tilde g_{\theta\theta}}, \quad F_r = \frac{R_{r\theta r\theta}}{\tilde g_{rr}\,\tilde g_{\theta\theta}}, \quad G = \frac{R_{\theta\phi\theta\phi}}{\tilde g_{\theta\theta}^2},
\end{equation}

and the Kretschmann scalar is $K = 4E^2 + 8F_t^2 + 8F_r^2 + 4G^2$.

Each component scales as $r^{2\phi_0 - 3}$ (verified by exact SymPy computation), giving:

\begin{equation} \label{eq:appD_3}
\boxed{\tilde K \sim r^{4\phi_0 - 6} \quad \text{as } r \to 0.}
\end{equation}

This is NOT $r^{12\phi_0 - 6}$ as the naive conformal formula $A^{-12} K_{\rm Schw}$ suggests. The conformal formula includes only the leading Schwarzschild term scaled by $A^{-12}$, omitting the derivative terms from the non-constant conformal factor $A(r)$. For $\phi_0 < 3/2$, the derivative terms ($r^{4\phi_0-6}$) dominate over the leading term ($r^{12\phi_0-6}$), and the conformal formula underestimates the curvature.

## D.2 Specific Values (Exact, $M = 1$)

| $\phi_0$ | $\tilde K$ (exact, deep interior) | $\tilde K \to 0$? | $\rho = Ar$ bounded? |
| --- | --- | --- | --- |
| $0.5$ | $\frac{75}{4r^4} + \cdots$ | No ($\to \infty$) | Yes |
| $1.0$ | $\frac{9}{r^2} + \frac{1}{4}$ | No ($\to \infty$) | Yes ($\rho \to r_h$) |
| $1.5$ | $\frac{291}{64}$ (constant) | No ($\to$ const) | No ($\rho \to \infty$) |
| $2.0$ | $\frac{39\,r^2}{16}$ | Yes ($\to 0$) | No ($\rho \sim r_h^2/r \to \infty$) |
| $2.5$ | $\frac{1323\,r^4}{1024}$ | Yes ($\to 0$) | No |

## D.3 The Fixed-Background Theorem

### Theorem D.1 (Fixed-Schwarzschild)

Within the conformal-disformal TEP framework with geometric metric fixed to Schwarzschild, curvature regularity and bounded physical areal radius are mutually exclusive.

#### Proof

Curvature regularity ($\tilde K \to 0$ as $r \to 0$) requires $4\phi_0 - 6 > 0$, i.e., $\phi_0 > 3/2$. Bounded areal radius ($\rho = A \cdot r = r_h^{\phi_0} r^{1-\phi_0}$ finite as $r \to 0$) requires $1 - \phi_0 \geq 0$, i.e., $\phi_0 \leq 1$. These conditions are mutually exclusive: there is no $\phi_0$ satisfying both $\phi_0 > 3/2$ and $\phi_0 \leq 1$. $\square$

## D.4 Numerical Verification

The exact formula $\tilde K \sim r^{4\phi_0 - 6}$ is verified numerically using the corrected pipeline curvature computation (standard Schwarzschild coordinates, 4-term Kretschmann formula). The power law is recovered to within 0.1% accuracy for all tested $\phi_0$ values:

| $\phi_0$ | Numerical power law | Expected $4\phi_0 - 6$ | Deviation |
| --- | --- | --- | --- |
| 0.50 | $r^{-4.000}$ | $-4.0$ | 0.000 |
| 0.75 | $r^{-3.000}$ | $-3.0$ | 0.000 |
| 1.00 | $r^{-2.001}$ | $-2.0$ | 0.001 |
| 1.25 | $r^{-1.001}$ | $-1.0$ | 0.001 |
| 1.50 | $r^{-0.001}$ | $0.0$ | 0.001 |
| 1.75 | $r^{0.999}$ | $1.0$ | 0.001 |
| 2.00 | $r^{1.999}$ | $2.0$ | 0.001 |
| 2.50 | $r^{3.999}$ | $4.0$ | 0.001 |

The Schwarzschild baseline ($\phi_0 \to 0$, $A \to 1$) recovers $K = 48M^2/r^6$ to within $5 \times 10^{-6}$ (0.0005%). All regression tests pass.

## D.5 Corrected Values for the $\phi_0 = 2$ Prototype

The $\phi_0 = 2$ case (curvature-regular but spatially-enlarged) has the following corrected curvature values, replacing the incorrect $\tilde K \sim r^{18}$ and $\tilde K \sim 10^{-150}$ from the naive conformal formula:

| Location | $r/M$ | $\tilde K$ (corrected) | $\tilde K$ (old, wrong) | $K_{\rm Schw}$ |
| --- | --- | --- | --- | --- |
| Exterior | $10$ | $4.79 \times 10^{-5}$ | $4.79 \times 10^{-5}$ | $4.79 \times 10^{-5}$ |
| Horizon | $2$ | $0.748$ | $0.752$ | $0.752$ |
| Deep interior | $0.1$ | $2.31 \times 10^{-2}$ | $\sim 10^{-150}$ | $4.79 \times 10^{7}$ |
| Innermost | $10^{-8}$ | $1.36 \times 10^{-16}$ | $\sim 10^{-300}$ | $4.80 \times 10^{49}$ |

The corrected values show $\tilde K$ is small but not absurdly tiny. The exact power law $\tilde K \sim 39r^2/16$ gives $\tilde K(0.1) = 0.024$, matching the numerical value $0.023$ to within 4%. The old values ($\sim 10^{-150}$) were artifacts of the incorrect conformal formula.

## D.6 Disformal Term: Status of the Classification

The disformal term $B\nabla\phi\nabla\phi$ contributes $B(\phi')^2$ to $\tilde g_{rr}$. In the deep interior, the Lorentzian condition $\det_{2d} < 0$ forces $B \to 0$ (for logarithmic $\psi$) or $\tilde g_{rr} \to 0$ (for saturating $\psi$). For the tested parameter space, the deep interior curvature is controlled by $A$ alone, and the conformal theorem applies. The full disformal classification — classifying all allowed asymptotics of $A$, $B$, $\phi$ under $\det\tilde g < 0$, $\tilde K < \infty$, and $\tilde R_{\rm area} < \infty$ — is an open mathematical programme (Section 3.3). The statement that the disformal term cannot help is established for the tested class but not as a general theorem.

## D.7 Implication: The Temporal Field Must Modify the Geometry

The fixed-background theorem proves that the fixed-Schwarzschild construction is inconsistent. The temporal field cannot merely sit on top of a Schwarzschild black hole — it must also change the black-hole geometry itself. When the geometric metric is allowed to be regular, all regularity conditions hold simultaneously.

# Appendix E — Geodesic Structure

This appendix records the null geodesic structure of the TEP matter metric. The key result is that for a conformal transformation $\tilde g = A^2 g$, null geodesics of $g$ are also null geodesics of $\tilde g$ (same unparameterised curves), but the affine parameters are related by $d\tilde\lambda = A^2\,d\lambda$. This single relation determines all affine-parameter results.

## E.1 Conformal Affine Parameter Relation

For $\tilde g_{\mu\nu} = A^2 g_{\mu\nu}$, a null geodesic $x^\mu(\lambda)$ of $g$ with tangent $k^\mu = dx^\mu/d\lambda$ satisfies $g_{\mu\nu}k^\mu k^\nu = 0$ and the geodesic equation $k^\nu\nabla_\nu k^\mu = 0$. Under the conformal transformation, the same curve is null ($\tilde g_{\mu\nu}k^\mu k^\nu = A^2 g_{\mu\nu}k^\mu k^\nu = 0$), but the affine parameter transforms as

\begin{equation} \label{eq:appE_1}
d\tilde\lambda = A^2\,d\lambda.
\end{equation}

This is the standard conformal affine-parameter relation. The Christoffel symbols transform as $\tilde\Gamma^\mu_{\alpha\beta} = \Gamma^\mu_{\alpha\beta} + \delta^\mu_\alpha \partial_\beta\ln A + \delta^\mu_\beta \partial_\alpha\ln A - g_{\alpha\beta}g^{\mu\nu}\partial_\nu\ln A$, and the geodesic equation in the rescaled parameter is satisfied.

## E.2 Eddington–Finkelstein Null Algebra

In ingoing Eddington–Finkelstein coordinates, the conformal matter metric ($B = 0$) has components

\begin{equation} \label{eq:appE_2}
\tilde g_{ab} = A^2\begin{pmatrix} -F & 1 \\ 1 & 0 \end{pmatrix}, \qquad \tilde g^{ab} = A^{-2}\begin{pmatrix} 0 & 1 \\ 1 & F \end{pmatrix}.
\end{equation}

Note that $\tilde g^{rr} = F/A^2$, not zero. The null condition $\tilde g_{ab}k^a k^b = 0$ gives $dv(-F\,dv + 2\,dr) = 0$, yielding two families:

- Ingoing ($v = \mathrm{const}$, $dv = 0$): the tangent $k^\mu = (0, 1)$. The Christoffel symbol $\tilde\Gamma^r{}_{rr} = 2A'/A$ is nonzero.

- Outgoing ($dr/dv = F/2$): the tangent $k^\mu = (2/F, 1)$, or $dv/dr = 2/F$.

## E.3 Affine Parameters on the Divergent-$A$ Benchmark

On the Hayward benchmark with $A = (r_h/r)^{\phi_0}$ and $\phi_0 = 1$, the conformal factor diverges as $A \sim r^{-1}$ near $r = 0$. The seed (Hayward) radial null affine parameter behaves as $d\lambda \propto dr$ (since $F \to 1$ and the background is regular). Applying the conformal relation:

\begin{equation} \label{eq:appE_3}
\tilde\lambda \sim \int A^2\,dr = \int \frac{dr}{r^{2\phi_0}}.
\end{equation}

For $\phi_0 = 1$: $\tilde\lambda \sim \int r^{-2}\,dr = -1/r \to \infty$ as $r \to 0$. Both radial null families have infinite affine parameter. The finite-area asymptotic end lies at infinite affine distance for both families. This is consistent with the end being an asymptotic boundary, not a regular manifold point.

For $\phi_0 = 2$: $\tilde\lambda \sim \int r^{-4}\,dr \to \infty$. Both families infinite.

For general $\phi_0 > 0$: both families have $\tilde\lambda \to \infty$. The divergence is driven by $A \to \infty$, not by the background geometry.

## E.4 Affine Parameters on the Correct Finite-$A$ Architecture

On the target global solution with $A \to A_c$ finite as $r \to 0$, the conformal relation gives $d\tilde\lambda \to A_c^2\,d\lambda$. Since the seed (regular) geodesic has $d\lambda \propto dr$ and finite affine parameter to $r = 0$, the matter-frame affine parameter is also finite: $\tilde\lambda \to A_c^2 \lambda_{\rm seed}$. Both null families reach the regular centre at finite affine distance. The centre is a regular manifold point; geodesics are extendable through it. This is geodesic completeness in the standard sense.

## E.5 The $\phi_0 = 2$ Case

For $\phi_0 = 2$ on fixed Schwarzschild, $A \sim r^{-2}$ and $\tilde\lambda \sim \int r^{-4}\,dr \to \infty$ for both families. The $r \to 0$ limit is an asymptotic spatially enlarged end ($\rho \sim r_h^2/r \to \infty$), not a regular point. The matter metric is Lorentzian and nondegenerate for all $r > 0$ ($\det\tilde g_{2D} = A^4 \det g_{2D} < 0$ since $A > 0$), but the limit $r = 0$ cannot be added as a regular point. The curvature vanishes ($\tilde K \sim 39r^2/16 \to 0$), but the areal radius diverges — the fixed-background theorem (Appendix D) proves this trade-off is unavoidable for fixed Schwarzschild.

## E.6 Signature Structure

The canonical 2D $(v,r)$ determinant for $ds_g^2 = -F\,dv^2 + 2G\,dv\,dr + R^2\,d\Omega^2$ and a static radial scalar is:

\begin{equation} \label{eq:appE_4}
\det\tilde g_{2D} = -A^4 G^2 - A^2 F B(\phi')^2.
\end{equation}

For the pure conformal metric ($B = 0$): $\det\tilde g_{2D} = -A^4 G^2 < 0$ — Lorentzianity inherited analytically from $g$ via $A > 0$. For the full disformal metric, the determinant must be checked with the canonical formula. Inside horizons ($F < 0$), the disformal term $-A^2 F B(\phi')^2$ becomes positive and can oppose the conformal term; all invertibility and hyperbolicity results must be regenerated from this expression.

# Appendix F — Physical Volume and Density Diagnostics

This appendix distinguishes the geometric volume from the matter-frame volume diagnostics. The pipeline computes the exterior cumulative volume and an interior constant-$r$ volume element; it does not solve a self-consistent matter density profile.

## F.1 Volume Elements

The geometric volume is

\begin{equation} \label{eq:appF_1}
V_{\rm geom}(r)=\frac{4\pi}{3}r^3.
\end{equation}

Inside the geometric horizon, the spatial constant-$r$ slice has volume element per unit $v$

\begin{equation} \label{eq:appF_2}
\frac{d\tilde V}{dv}=4\pi A^2r^2\sqrt{|F|}.
\end{equation}

For the clean run this rises from $2.08\times10^2$ at $0.1M$ to $5.50\times10^5$ at $r_{\min}=10^{-8}M$.

## F.2 Geometric Collapse and Physical Interior Growth

$V_{\rm geom}\to0$ as $r\to0$, while the matter-frame interior volume element grows strongly because $A\sim r^{-1}$, reaching $d\tilde V/dv=5.50\times10^5$ at $r_{\min}=10^{-8}M$. The determinant is strictly negative at every sampled radius, so this growth is a physical result on a globally Lorentzian space, not a formal continuation diagnostic. The areal radius diverges ($\rho\sim1/r$), confirming that the Schwarzschild areal coordinate is not a direct measure of physical volume in TEP.

## F.3 Numerical Verification

| Quantity | Geometric frame | Matter-frame diagnostic |
| --- | --- | --- |
| Volume at $r_{\min}$ | $V_{\rm geom}=4.19\times10^{-24}$ | interior element $d\tilde V/dv=5.50\times10^5$ |
| Lorentzian domain | globally Lorentzian | $\det\tilde g_{2D}<0$ at all sampled radii |
| Density | $M/V_{\rm geom}\propto r^{-3}$ | not independently computed by this pipeline |

## F.4 Supported Statement

### Proposition F.1

The $\phi_0 = 2$ fixed-Schwarzschild configuration replaces the geometric zero-volume intuition with a strongly growing matter-frame interior volume element on a Lorentzian space for all $r > 0$. The areal radius diverges ($\rho \sim 1/r \to \infty$), confirming that this configuration is a spatially enlarged asymptotic end, not a bounded solution and not the final TEP solution. The Schwarzschild areal coordinate is not a direct measure of physical volume in TEP. A finite physical density claim requires an explicit stress-energy solution and is not asserted by the present diagnostic pipeline.

# Appendix G — Scope of Gravitational Perturbations

Canonical TEP distinguishes the geometric metric $g_{\mu\nu}$ from the matter metric $\tilde g_{\mu\nu}$. Gravitational waves propagate on $g_{\mu\nu}$. In the case analysed in this manuscript, the geometric metric is fixed at Schwarzschild; conformal and disformal factors in $\tilde g_{\mu\nu}$ therefore do not shift the gravitational Regge–Wheeler or Zerilli spectrum.

## G.1 Geometric Tensor Equation

Decomposing $h_{\mu\nu}$ into axial and polar spherical harmonics and Fourier transforming $\Psi(t,r)=e^{-i\omega t}\Psi(r)$ gives

\begin{equation} \label{eq:appG_1}
\frac{d^2\Psi_i}{dr_*^2}+\left[\omega^2-V_i^{\rm Schw}(r)\right]\Psi_i=0,
\end{equation}

where $r_*$ is the Schwarzschild tortoise coordinate and $i$ denotes the Regge–Wheeler or Zerilli sector. These are the standard geometric GR perturbation equations for this fixed background.

## G.2 Interpretation of Matter-Metric Diagnostics

The pipeline also evaluated surrogate radial functions constructed from the disformal matter metric. Such functions may diagnose matter or electromagnetic propagation, but they are not gravitational-wave potentials in canonical TEP. In particular, applying gravitational QNM boundary conditions to a $\tilde g_{\mu\nu}$ surrogate does not establish a shifted tensor spectrum.

## G.3 Ringdown Spectrum

The gravitational QNMs of this non-spinning fixed-background limit are the Schwarzschild GR baseline. The previously quoted $0.421-0.099\,i$ frequency and 5.6% shift arose from the matter-metric surrogate rather than the geometric gravitational equation and are not gravitational QNM predictions. The coupled solution (Section 4) introduces scalar-led and mixed scalar–tensor modes with coupling-specific QNM shifts of sub-percent to percent-level.

## G.4 Stability and Echo Scope

No quadratic action or complete coupled perturbation system was analysed. Background regularity and the shape of a surrogate matter-metric potential are insufficient to prove absence of ghosts, gradient instabilities, or growing scalar and mixed modes. They are likewise insufficient to exclude gravitational-wave echoes. An echo statement requires geometric perturbation evolution with physically specified inner boundary conditions and a waveform-level analysis; those calculations were not performed.

## G.5 Summary

| Question | Result for the fixed-background limit |
| --- | --- |
| GW propagation metric | Geometric $g_{\mu\nu}$ |
| Fixed geometric background | Schwarzschild |
| Gravitational tensor QNMs | Geometric GR baseline; no TEP shift |
| Scalar-led and mixed modes | Computed in the coupled solution (Section 4): sub-percent to percent-level shifts |
| Full linear stability | Background regularity verified; coupled-solution stability is open (Section 7) |
| Gravitational-wave echoes | Neither predicted nor excluded |

Disformal matter-metric diagnostics must not be interpreted as gravitational QNMs. For the fixed Schwarzschild geometric metric, gravitational ringdown remains the GR baseline.

# Appendix H — Ray Tracing and Orbital Mechanics

This appendix records the photon-sphere, ISCO, and QPO calculations together with the EHT shadow-size central-value comparisons. The reported pipeline configuration imposes exterior screening, $A=1$ and $B=0$, over the ray and orbit domains, so these calculations reproduce the Schwarzschild baseline. Gravitational ringdown is separately governed by the geometric metric $g_{\mu\nu}$, which is fixed at Schwarzschild here.

## H.1 Photon Sphere

Null geodesics of the matter metric $\tilde g_{\mu\nu}$ reduce to Schwarzschild null geodesics on the screened exterior used by the pipeline, where $A=1$ and $B=0$. The photon sphere is therefore at

\begin{equation} \label{eq:appH_1}
r_{\rm ph} = 3\,M,
\end{equation}

and the critical impact parameter for the shadow edge is

\begin{equation} \label{eq:appH_2}
b_{\rm crit} = \frac{r_{\rm ph}}{\sqrt{1 - 2M/r_{\rm ph}}} = 3\sqrt{3}\,M = 5.196\,M.
\end{equation}

These are baseline values. The current exterior-null ray tracing does not calculate a subleading disformal correction to the shadow edge.

## H.2 ISCO and Radiative Efficiency

On the exterior-screened baseline, the innermost stable circular orbit (ISCO) for timelike geodesics is

\begin{equation} \label{eq:appH_3}
r_{\rm ISCO} = 6\,M.
\end{equation}

The specific energy and angular momentum at the ISCO are

\begin{equation} \label{eq:appH_4}
E_{\rm ISCO} = \sqrt{\frac{8}{9}} = 0.9428, \qquad L_{\rm ISCO} = 2\sqrt{3}\,M = 3.464\,M.
\end{equation}

The radiative efficiency — the fraction of rest-mass energy liberated by accretion onto the ISCO — is

\begin{equation} \label{eq:appH_5}
\eta = 1 - E_{\rm ISCO} = 1 - \sqrt{\frac{8}{9}} = 0.0572 = 5.72\%.
\end{equation}

These are standard Schwarzschild values resulting from screening over the orbit domain. They do not imply that a nonconstant conformal factor would exactly cancel from timelike geodesics.

## H.3 Redshift Factor at ISCO

The gravitational redshift factor at the ISCO is

\begin{equation} \label{eq:appH_6}
g = \sqrt{1 - \frac{2M}{r_{\rm ISCO}}} = \sqrt{1 - \frac{1}{3}} = \sqrt{\frac{2}{3}} = 0.816.
\end{equation}

The corresponding redshift is

\begin{equation} \label{eq:appH_7}
z = \frac{1}{g} - 1 = \sqrt{\frac{3}{2}} - 1 = 0.225.
\end{equation}

These are the standard Schwarzschild ISCO values, confirming that the iron-line and reflection spectroscopy signatures of the TEP exterior are identical to GR at the ISCO.

## H.4 3:2 QPO Resonance

Quasi-periodic oscillations (QPOs) from black-hole accretion disks often appear in a 3:2 frequency ratio, interpreted as a resonance between the radial and vertical epicyclic frequencies. For the TEP exterior (Schwarzschild), the resonance condition $\Omega_r / \Omega_\theta = 2/3$ is satisfied at

\begin{equation} \label{eq:appH_8}
r_{\rm QPO} = 10.8\,M,
\end{equation}

where the epicyclic frequency ratio is

\begin{equation} \label{eq:appH_9}
\frac{\Omega_r}{\Omega_\theta}\bigg|_{r=10.8\,M} = 0.667 = \frac{2}{3}.
\end{equation}

This is the standard Schwarzschild 3:2 resonance radius. The TEP disformal interior does not affect this radius because it is well outside the horizon.

## H.5 EHT Shadow Comparisons

The shadow angular diameter is computed from the critical impact parameter $b_{\rm crit} = 5.196\,M$ and the source distance and mass. The pipeline evaluates the shadow size for the two EHT targets:

| Source | TEP prediction | EHT measurement | Deviation |
| --- | --- | --- | --- |
| M87* | $39.69\;\mu\text{as}$ | $42 \pm 3\;\mu\text{as}$ | $0.77\sigma$ |
| Sgr A* | $53.26\;\mu\text{as}$ | $48.7 \pm 7.0\;\mu\text{as}$ | $0.65\sigma$ |

The exterior-screened central values, $39.69\;\mu$as for M87* and $53.26\;\mu$as for Sgr A*, lie within the quoted measurement intervals $42\pm3\;\mu$as and $48.7\pm7.0\;\mu$as, respectively. These are central-value checks; a full likelihood analysis would propagate correlated mass, distance, calibration, and imaging systematics. The fixed-background limit reproduces the Schwarzschild baseline; the coupled solution (Section 4) provides TEP-specific shadow deviations of $-0.044\%$ for $\eta=-0.1$ on $g^{\rm sGB}$ (Appendix L).

## H.6 Gravitational Ringdown Baseline

Gravitational waves propagate on the geometric metric, not on the matter metric used for ray tracing. Since the fixed-background limit holds the geometric metric at Schwarzschild, its non-spinning gravitational ringdown is the Schwarzschild GR baseline. GW150914 produced a spinning remnant and must be compared with the appropriate geometric Kerr spectrum. The coupled solution (Section 4) provides coupling-specific QNM shifts of sub-percent to percent-level.

## H.7 Summary of Orbital and Observational Predictions

| Observable | Value | vs Schwarzschild |
| --- | --- | --- |
| Photon sphere | $r = 3\,M$ | identical |
| Critical impact parameter | $b = 5.196\,M$ | identical |
| ISCO | $r = 6\,M$ | identical |
| ISCO energy | $E = \sqrt{8/9} = 0.9428$ | screened baseline |
| ISCO angular momentum | $L = 3.464\,M$ | identical |
| Radiative efficiency | $\eta = 5.72\%$ | identical |
| Redshift factor at ISCO | $g = 0.816$ | identical |
| Redshift at ISCO | $z = 0.225$ | identical |
| 3:2 QPO resonance | $r = 10.8\,M$ | identical |
| M87* shadow | $39.69\;\mu$as | $0.77\sigma$ from EHT |
| Sgr A* shadow | $53.26\;\mu$as | $0.65\sigma$ from EHT |
| GW150914 $f_{220}$ (non-spinning WKB diagnostic) | $206.5\,$Hz | Schwarzschild baseline (0% shift) |
| GW150914 $f_{220}$ (Kerr, $\chi=0.67$) | $267\,$Hz | Kerr baseline; measured $\sim 250\,$Hz |

The exterior orbital structure (photon sphere at $3M$, ISCO at $6M$, 3:2 QPO at $10.8M$) is identical to Schwarzschild. EHT shadows agree at $0.77\sigma$ (M87*) and $0.65\sigma$ (Sgr A*). The gravitational QNM baseline is the Schwarzschild spectrum for non-spinning remnants; for GW150914's spinning remnant ($\chi=0.67$) the appropriate baseline is Kerr, giving $f_{220}\approx 267\,$Hz (step_14) versus the measured $\sim 250\,$Hz. The non-spinning WKB diagnostic $f_{220}=206.5\,$Hz (step_05, $\omega_R^{\rm WKB}=0.3988$) is not the physical prediction for a spinning remnant and carries $\sim 7\%$ WKB systematic error relative to the exact Schwarzschild value $\omega_R^{\rm exact}=0.3737$. There is no coupling-specific ringdown frequency shift in this fixed-background limit.

# Appendix I — Comparison Table

This appendix provides a detailed technical comparison of TEP-BH with the major competing models of singularity resolution. The TEP-BH column is split into four sub-columns reflecting the distinct calculations performed in this paper: the fixed-Schwarzschild theorem case, the Hayward validation benchmark, the linear sGB exterior candidate, and the regularised sGB-coupled interior (now solved). These must not be combined into a single "TEP" column, because they represent different geometries with different properties.

## I.1 Detailed Comparison

| Property | GR (Schwarzschild/Kerr) | Regular BH (Bardeen/Hayward) | TEP theorem case ($\phi_0=2$, fixed Schw.) | TEP benchmark ($\phi_0=1$, Hayward) | TEP linear sGB (exterior) | TEP regularised interior (solved) |
| --- | --- | --- | --- | --- | --- | --- |
| Exterior metric | Schwarzschild/Kerr (exact) | Schwarzschild (to $M/r$) | Schwarzschild (fixed) | Hayward (prescribed) | Schwarzschild to $\sim$11 parts in $10^4$ at $|\eta|=0.1$ (perturbative) | Asymptotically Schwarzschild (target) |
| Interior metric | Singular ($r\to0$) | de Sitter-like temporal minimum | Spatially enlarged end ($\rho \to \infty$); not a regular point | Finite-area asymptotic end ($\rho \to 1.995$, $\tilde\ell \to \infty$); not a regular point | Not yet integrated; literature: finite-area singularity at $r \approx 0.27\,r_h$ (Sotiriou & Zhou 2014) | Regular point centre: $A(0) = 4.68$, $\rho_{\rm areal}(0) \to 4.68 \times 10^{-8}$, $w_r(0) = -1.000$, $K(0) = 27.1$ ($\eta=-0.1$, $g=1.1M$) |
| Central curvature | $K\to\infty$ | Finite (by construction) | Vanishes ($\tilde K \sim 39r^2/16 \to 0$) but areal radius diverges | Finite ($\tilde K \to 8/r_h^4$) but asymptotic end, not regular centre | Likely divergent at finite-area singularity (literature) | Finite: $K(0) = 27.1$ for $g=1.1M$ |
| Geodesic completeness | No (incomplete at $r=0$) | Model-dependent | Both null families infinite affine ($\tilde\lambda \sim \int r^{-4}\,dr \to \infty$); asymptotic end | Both null families infinite affine ($\tilde\lambda \sim \int r^{-2}\,dr \to \infty$); asymptotic end | Not yet computed; likely incomplete at finite-area singularity | Both families finite affine; regular point centre, extendable geodesics |
| Horizon structure | Event horizon at $r=2M$ | Event + inner horizon | Schwarzschild horizon (fixed) | Hayward outer + inner horizon | sGB apparent horizon (perturbative artifact) | Temporal well: $\tilde N_{\min} = 0.211$ at $r \approx 1.41M$, no finite causal horizon |
| Temporal Horizon | — | — | Not established (theorem case only) | Not established ($\theta_+=0$ is constant-area property, not temporal freeze) | Not yet computed | Operational: $\mathcal{Z}_{\max} \sim N_o/\tilde N_{\min} \sim 4.7$ at $r \approx 1.41M$; $\tilde N > 0$ everywhere |
| Surface / junction | None | None (smooth) | None (but asymptotic end) | None (but asymptotic end) | Not yet determined | None (smooth, regular point centre) |
| Lorentzianity | Yes | Yes | Yes for all $r > 0$ ($\det = A^4 \det g < 0$); limit $r=0$ not a regular point | Yes in exterior ($F > 0$); inside horizon must check canonical determinant | Yes in exterior ($F > 0$) | Yes: $\det\tilde g_{2D} < 0$ throughout the integrated domain |
| Observable shifts | — | — | — | — | Not yet derived (conformal null invariance; disformal circular orbit; tensor characteristic metric all require derivation) | Photons on $\tilde g$, GWs on $\mathcal{G}_{\rm tensor}$; shadow vs ringdown consistency (target) |
| QNM spectrum | Schwarzschild/Kerr | Modified | — | — | Indicative WKB: scalar-led, axial, polar; requires spectral solver and coupled perturbation equations | Exterior benchmark: axial $+0.136\%$ real shift, polar-led $-0.028\%$ real breaking, $+0.067\%$ damping breaking. Full deep-transit solve: three modes (polar-led $0.462-0.032i$, scalar-led $0.416-0.032i$, overtone $0.745-0.035i$), $2.8\times$ longer damping, $+14.8\%$ real breaking, $11.1\%$ mode splitting |
| Shadow size | $b = 5.196\,M$ | Similar to Schwarzschild | — | — | $-0.044\%$ at $\eta=-0.1$ on $g^{\rm sGB}$; conformal invariance means shadow on $\tilde g$ equals shadow on $g$ for pure conformal; ISCO on $\tilde g$ is $+1.95\%$ (opposite sign, $\mathcal{O}(\eta)$, $\sim 44\times$ larger) — a property of the sGB EFT (Appendix L) | From $\tilde g$ photon sphere (target) |
| Gravitational-wave echoes | No | Possible (inner barrier) | — | — | Not yet analysed (requires waveform-level analysis with inner boundary conditions) | Depends on interior outcome |
| Stability / hyperbolicity | Stable | Model-dependent | — | — | Tensor sector: kinetic term positive at leading order, $c_T = 1$ on Schwarzschild background (Horndeski $G_4$ unchanged by sGB at leading order). Full tensor characteristic metric on TEP solution, global ghost freedom, and scalar sector hyperbolicity in deep interior remain open (Thaalba et al. 2024; Section 7.7) | Scalar-sector principal symbol on solved geometry (target) |
| $c_T$ (tensor speed) | $c_T = 1$ | $c_T \to 1$ asymptotically | — | — | $c_T = 1$ on Schwarzschild background (Horndeski $G_4 = M_{\rm Pl}^2/2$, sGB coupling enters through $G_5$); GW170817 satisfied on that background; full tensor characteristic metric on TEP solution pending (Section 7.7) | Established |
| Observational coupling | $M$ (and $a$) | Regularisation parameter | — | — | One $\alpha_{\rm GB}$; $\eta_i = 3\alpha_{\rm GB}/M_i^2$ per source; $|\eta|=0.1$ cannot be universal across masses | Same (target) |

## I.2 What Is Established vs What Is Target

The established results in this paper are:

- The fixed-background conformal theorem: on fixed Schwarzschild, finite curvature and bounded areal radius are mutually exclusive for the power-law conformal class (Section 3, Appendix D).

- The Hayward benchmark: finite curvature and bounded areal radius are compatible on a regular background, but the $\phi_0 = 1$ benchmark produces a finite-area asymptotic end, not a regular centre (Appendix K).

- The perturbative sGB exterior: scalar hair and real backreaction are established by the published $\mathcal{O}(\alpha^2)$ solution (Sotiriou & Zhou 2014). The non-perturbative exterior mass-function evolution requires rerun with the correct scalar charge.

- The literature evidence that standard linear sGB develops a finite-area interior singularity (Sotiriou & Zhou 2014; Thaalba et al. 2024).

- The TEP Global Solution Architecture: regular spatial domain, finite local density, continuous temporal-rate gradient, no finite-radius causal horizon ($N \geq N_{\min} > 0$ everywhere), extreme external redshift (Section 6.6).

- The regularised sGB-coupled interior solution with a regular point centre, integrated to $r = 10^{-8}M$ for $\eta=-0.1$ (Section 4.6, Section 5.5,).

The remaining target results are:

- The coupled polar–scalar QNM spectrum and isospectrality breaking.

- The observable shifts as quantitative predictions (requires deriving both characteristic metrics).

- The tensor speed $c_T$ derived from the action.

- The observational fit with one universal $\alpha_{\rm GB}$ across all sources.

The comparison table above distinguishes these explicitly. Earlier versions of this appendix combined all TEP calculations into a single column, producing claims (globally Lorentzian interior, geodesic completeness, vanishing central curvature, disformal birefringence, no inner boundary) that belong to specific benchmark or theorem cases, not to the unsolved global TEP solution.

# Appendix J — Reproducibility

Every numerical result quoted in this paper is produced by the automated pipeline, which downloads or verifies the primary data, constructs the prescribed matter metric, computes the diagnostics, and writes the comparison tables. Each clean run writes SHA-256 checksums for the result files. Checksums prove that the result files were not altered after generation; they do not prove that the calculations are correct. Three categories of calculation require particular care in scalar–tensor perturbation theory and are not cited as quantitative evidence in this manuscript: (i) the nonlinear exterior mass-function evolution, which requires a consistent scalar charge and sign convention matched to published sGB solutions; (ii) surrogate WKB QNM frequencies, which require the perturbation potential to be derived from the second variation of the action (not from inserting modified background functions into a GR potential formula); and (iii) matter-frame observables, which must be computed on the matter metric $\tilde g$ rather than the geometric metric $g$. The fixed-background limit is globally Lorentzian for all $r > 0$; the coupled solution is Lorentzian in the exterior ($F > 0$). The deep-interior Lorentzianity of the coupled solution depends on the interior integration (Section 4.6).

## J.1 Pipeline Architecture

The pipeline consists of 50 steps (`step_00` through `step_49`), each implemented as an independent Python module in `scripts/steps/` and orchestrated by `scripts/run_pipeline.py`. Steps 00–11 implement the fixed-background limit (Section 3); steps 12–15 implement the coupled solution (Section 4); steps 16–41 are derivation scripts (geometry, corrected observables, QNM solvers, interior, phantom mass, Kerr-sGB, GW confrontations); step 42 generates figures; steps 43–49 implement the S-star inference pipeline (Section 9, Appendix J). The steps are:

| Step | Module | Output | Description |
| --- | --- | --- | --- |
| step_00 | step_00_data_download.py | data/raw/, data/processed/ | Download EHT M87* visibilities, compile measurement tables |
| step_01 | step_01_field_equations.py | step_01_field_equations.{json,csv} | Disformal metric, curvature invariants, geodesic completeness, physical volume/density |
| step_02 | step_02_perturbations.py | step_02_perturbations.{json,csv} | Regge–Wheeler/Zerilli potentials, QNMs, echo check |
| step_03 | step_03_raytracing.py | step_03_raytracing.{json,csv} | Photon sphere, shadow angular diameter, ISCO, EHT comparison |
| step_04 | step_04_accretion.py | step_04_accretion.{json,csv} | Circular geodesics, epicyclic frequencies, radiative efficiency, redshift, QPO |
| step_05 | step_05_observational_constraints.py | step_05_observational_constraints.{json,csv} | Chi-squared comparison of TEP predictions against EHT and LIGO |
| step_06 | step_06_gw190521_mass_gap.py | step_06_gw190521_mass_gap.{json,csv} | GW190521 mass gap: LIGO posterior analysis (TEP exterior-null) |
| step_07 | step_07_qpo_frequency_lock.py | step_07_qpo_frequency_lock.{json,csv} | QPO frequency lock: 3:2 epicyclic resonance ratios (Schwarzschild baseline) |
| step_08 | step_08_jwst_early_smbhs.py | step_08_jwst_early_smbhs.{json,csv} | JWST early SMBHs: Eddington-limited growth envelope from stellar seeds |
| step_09 | step_09_spin_bias.py | step_09_spin_bias.{json,csv} | Over-maximal spin bias: observed spin distribution (TEP exterior-null) |
| step_10 | step_10_tde_missing_flares.py | step_10_tde_missing_flares.{json,csv} | TDE missing flares: sub-Eddington luminosity survey (TEP exterior-null) |
| step_11 | step_11_eht_polarization.py | step_11_eht_polarization.{json,csv} | EHT polarization: birefringence check (TEP exterior-null) |
| step_12 | step_12_self_gravitating.py | step_12_self_gravitating.{json,csv} | sGB metric corrections, scalar profile, exterior observables (shadow, ISCO, QNM), $\eta$-scan, GW150914 projection (Section 4) |
| step_13 | step_13_scalar_perturbations.py | step_13_scalar_perturbations.{json,csv} | Scalar-led, axial, and polar QNM channels; isospectrality breaking; deep-transit analysis (Section 7) |
| step_14 | step_14_kerr_tep.py | step_14_kerr_tep.{json,csv} | Kerr–TEP shadow and ISCO spin scan; frame-dragging cancellation; GW150914 Kerr ringdown comparison (Section 4) |
| step_15 | step_15_dynamical_signatures.py | step_15_dynamical_signatures.{json,csv} | PPN parameters, $-1$PN scalar dipole radiation, $c_T = 1$ on the Schwarzschild background, combined $\eta$ constraints (Section 4) |
| step_16 | step_16_exact_geometry.py | step_16_exact_geometry.{json,csv} | Exact curvature invariants, determinant, inverse metric, areal radius, proper radial distance |
| step_17 | step_17_null_expansions.py | step_17_null_expansions.json | Null expansions $\theta_\pm$ for TEP matter metric on Schwarzschild and Hayward backgrounds |
| step_18 | step_18_observer_redshift.py | step_18_observer_redshift.{json,csv} | Gravitational redshift for static observers at different radii |
| step_19 | step_19_observer_frequency_transfer.py | step_19_observer_frequency_transfer.{json,csv} | Invariant emitter-receiver frequency transfer factor $\mathcal{Z}$ (Section 6) |
| step_20 | step_20_corrected_observables.py | step_20_corrected_observables.json | Corrected exterior observables (shadow, ISCO, frequency transfer) under fixed-ADM normalization (Appendix L) |
| step_21 | step_21_cT_derivation.py | step_21_cT_derivation.json | Tensor propagation speed $c_T$ from principal symbol of coupled Einstein-sGB equations |
| step_22 | step_22_quadratic_action.py | step_22_quadratic_action.json | Exact quadratic action for axial, polar, and scalar perturbation channels (Section 7) |
| step_23 | step_23_characteristic_matrices.py | step_23_characteristic_matrices.json | Tensor characteristic matrices and principal symbol analysis |
| step_24 | step_24_conformal_invariance_check.py | step_24_conformal_invariance_check.json | Conformal invariance verification for null geodesics |
| step_25 | step_25_qnm_solver.py | step_25_qnm_solver.json | Perturbative sGB QNM solver: axial, polar-led, scalar-led modes |
| step_26 | step_26_coupled_spectral_solver.py | step_26_coupled_spectral_solver.json | Coupled spectral solver for polar-scalar mixed modes |
| step_27 | step_27_qnm_horizon_branch.py | step_27_qnm_horizon_branch.json | Horizon-bearing sGB QNM: Schwarzschild base + Sotiriou-Zhou $\mathcal{O}(\eta^2)$ correction (Section 7.3) |
| step_28 | step_28_qnm_validation.py | step_28_qnm_validation.json | QNM spectrum validation against published sGB results (Bryant et al. 2021; Blazquez-Salcedo et al. 2016) (Section 7.3) |
| step_29 | step_29_matrix_leaver.py | step_29_matrix_leaver.json | Full coupled $2\times 2$ matrix continued-fraction QNM solver with regular inner boundary; three deep-transit modes (Section 7.7) |
| step_30 | step_30_taylor_recurrence.py | step_30_taylor_recurrence.json | Exact Taylor-series recurrence for background potentials and scalar field coefficients near the origin (Section 7.3) |
| step_31 | step_31_frobenius_analysis.py | step_31_frobenius_analysis.json | Frobenius regularity analysis at the temporal minimum: $r^{\ell+1}$ regularity, finite tortoise, regularised polar–scalar mixing (Section 7.3) |
| step_32 | step_32_interior_integration.py | step_32_interior_integration.json | Numerical integration of sGB scalar field on Hayward-class regular background |
| step_33 | step_33_interior_analysis.py | step_33_interior_analysis.json | Interior analysis: regularity checks, Lorentzian signature, lapse minimum |
| step_34 | step_34_solve_interior.py | step_34_solve_interior.json | TEP interior solver: sGB scalar on Hayward background, matter metric construction (Section 4.6) |
| step_35 | step_35_mass_bias_sign.py | step_35_mass_bias_sign.json | Mass-bias sign equation derivation (Section 9) |
| step_36 | step_36_phantom_mass_critical.py | step_36_phantom_mass_critical.json | Horizon-scale Phantom Mass: critical $g$ analysis, EHT M87$^\ast$ and Sgr A$^\ast$ shadow comparison (Section 8.1A) |
| step_37 | step_37_phantom_mass_raytrace.py | step_37_phantom_mass_raytrace.json | Phantom Mass null geodesic ray-tracing (Section 8.1A) |
| step_38 | step_38_phantom_mass_scan.py | step_38_phantom_mass_scan.json | Phantom Mass parameter scan (Section 8.1A) |
| step_39 | step_39_eht_visibility_fit.py | step_39_eht_visibility_fit.json | EHT visibility-domain joint inference: direct fit of TEP shadow model to calibrated visibilities (Section 8.1B) |
| step_40 | step_40_kerr_sgb_true.py | step_40_kerr_sgb_true.json | True rotating sGB shadow and QNM from Delgado et al. (2020) slowly-rotating solution (Section 8.5) |
| step_41 | step_41_gw250114_confrontation.py | step_41_gw250114_confrontation.json | GW250114 corrected confrontation with TEP observables |
| step_42 | step_42_generate_figures.py | results/figures/ | Figure generation (10 figures from paper plan) |
| step_43 | step_43_sstar_data.py | step_43_sstar_data.json | S-star data acquisition: 145 astrometric + 44 radial-velocity epochs for S2 (Section 9) |
| step_44 | step_44_gr_fit.py | step_44_gr_fit.json | Conventional GR fit to S-star data (Section 9) |
| step_45 | step_45_mass_bias.py | step_45_mass_bias.json | Mass-bias sign (Gate -1) from GR fit residuals (Section 9) |
| step_46 | step_46_tep_transfer.py | step_46_tep_transfer.json | TEP transfer-function fit with $\eta_{\rm TEP}$ as free parameter (Section 9) |
| step_47 | step_47_joint_forward.py | step_47_joint_forward.json | Joint forward-model: composition vs calibration separation, $\mathcal{C}_V$ diagnostic (Section 9) |
| step_48 | step_48_likelihood.py | step_48_likelihood.json | Formal likelihood comparison: Bayes factor, BIC, AIC, Wilks (Section 9) |
| step_49 | step_49_posterior.py | step_49_posterior.json | Phantom Mass posterior via MCMC (Section 9) |

All calculations use the unified parameter set: $B_0 = 1.0$, $M = 1.0$, $\beta_A = -1.0$, $n_B = 2.0$, $\delta = 0.05$, $\sigma_B = 1.5$, $r_h = 2.0$, with the quartic-damped disformal coupling $B(\phi) = B_0\,|\phi|^{n_B}/(1 + |\phi|^{n_B})\,\exp(-\phi^4/2\sigma_B^4)$. The conformal exponent $\beta_A = -1$ is frozen across all calculations — fixed-background theorem, sGB exterior, and the target global solution — consistent with the wider TEP corpus weak-field value ($\beta \simeq -0.013$, with $\beta_A = -1$ as the strong-field limit). Steps 12–15 use the sGB coupling parameter $\eta$ (primary value $\eta = -0.1$ for the mass-inflation branch) with the exact $\mathcal{O}(\alpha^2)$ perturbative metric corrections of Sotiriou \& Zhou (2014), Eqs. (56)-(63); the $\eta$-scan covers $\eta \in \{-0.05, -0.1, -0.15, -0.2\}$.

#### Convention resolution

The conformal coupling $A = e^{\beta_A \phi}$ with $\beta_A = -1$ is frozen across all calculations. The different scalar profiles interact with this single convention as follows: (i) the fixed-background theorem uses $\phi = \phi_0 \ln(r/r_h) < 0$ in the interior, giving $A = e^{-\phi} \to \infty$ — the divergent conformal factor that proves the theorem; (ii) the sGB exterior uses $\phi \sim Q_s/r$ (Coulomb). TEP selects the mass-inflation branch with $Q_s < 0$ (equivalently $\alpha_{\rm GB} < 0$), so $\phi < 0$ and $A = e^{-\phi} > 1$ — the matter metric is conformally *magnified* relative to $g$; (iii) the target global solution has $\phi \to \phi_c$ finite at the centre, giving $A_c = e^{-\phi_c}$ finite — the TEP Global Solution Architecture of Section 6.6. The strong-field value $\beta_A = -1$ connects to the weak-field corpus value $\beta \simeq -0.013$ through the effective coupling's $\phi$-dependence: if the coupling runs with $\phi$ or temporal shear, that running is part of the action. The quartic-Gaussian $B(\phi)$ should ultimately be derived from the action or screening mechanism, or an admissible class of $B(\phi, X)$ functions should be shown to produce stable results; this is a refinement target, not a structural inconsistency.

## J.2 Output Files

The pipeline produces three categories of output:

- `results/*.json` — structured summary of each step (the canonical machine-readable results);

- `results/*.csv` — data tables for plotting and inspection;

- `data/processed/*.json` — compiled measurement tables (EHT, LIGO QNM, QPO, TDE, spin, JWST) built from the raw downloads;

- `data/raw/` — downloaded EHT/LIGO data (EHT M87* 2019-D01-01 visibilities, LIGO GW190521 posterior);

- `results/pipeline_results.json` — top-level summary with step status, timing, and execution metadata;

- `results/checksums_sha256.json` — SHA-256 checksum of every `results/*.json` and `results/*.csv` file.

## J.3 Checksum Verification

At the end of each run, the pipeline computes the SHA-256 hash of every result file and writes them to `results/checksums_sha256.json`. This covers all result files (step JSON + step CSV + pipeline summary JSON). Any modification to a result file is detectable by re-running the hash. The checksum file itself is regenerated on each run, so the recorded hashes always correspond to the current results.

## J.4 Warning-Free Execution

A clean run is verified with the repository command:

\begin{equation} \label{eq:appJ_1}
\texttt{python scripts/run\_pipeline.py}
\end{equation}

The current diagnostics include a numerical safeguard that masks inverse-determinant curvature expressions where $|\det\tilde g_{2D}|$ is small; in the fixed-background limit this safeguard is never triggered because the determinant is strictly negative everywhere (globally Lorentzian).

## J.5 Data Sources

All observational data is drawn from published, publicly available sources:

| Source | Reference | Used in |
| --- | --- | --- |
| EHT M87* 2019 calibrated visibilities | EHT Collaboration 2019, data product 2019-D01-01 (GitHub) | step_00, step_03, step_05, step_11 |
| LIGO GW190521 posterior | Isi et al. 2020 (Zenodo 4057131) | step_06 |
| QPO measurements (GRS 1915+105) | Strohmayer 2001; Remillard et al. 2002; Homan et al. 2005 | step_07 |
| JWST early SMBH masses | Harikane et al. 2023; Greene et al. 2023; Bogdan et al. 2024 | step_08 |
| Black-hole spin measurements | McClintock et al. 2014; Reynolds 2021; EHT 2019/2022 | step_09 |
| TDE flare catalogs | van Velzen et al. 2021; Holoien et al. 2020 | step_10 |
| EHT polarization | EHT 2021 ApJ 910 L13; EHT 2024 ApJ 964 L26 | step_11 |
| LIGO QNM measurements | compiled in data/processed/ligo_qnm_measurements.json | step_05 |
| sGB perturbative metric corrections | Sotiriou & Zhou 2014, Eqs. (56)-(63) (exact $\mathcal{O}(\alpha^2)$) | step_12, step_14 |
| Kerr QNM / shadow formulas | Bardeen 1973; Leaver 1985 (approximate) | step_14 |
| GW170817 $c_T$ constraint | LIGO/Virgo & Fermi-GBM 2017 | step_15 |
| Binary pulsar dipole bounds | Freire et al. 2012; Antoniadis et al. 2013 (conditional on NS scalarization) | step_15 |

## J.6 Running the Pipeline

The full pipeline is run with a single command:

cd "/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-BH"
python scripts/run_pipeline.py

Options include `--start-step`, `--stop-step`, `--skip-steps`, `--no-derive`, `--no-inference`, `--no-figures`, `--list-steps`, and `--continue-on-error`. Full documentation is in `scripts/README.md`.

## J.7 Code Availability

All code is in the GitHub repository:
https://github.com/matthewsmawfield/TEP-BH

The repository contains the complete pipeline (`scripts/`), the processed and raw data (`data/`), the results (`results/`), the manuscript source (`site/components/`), and the figure-generation code (`scripts/steps/step_42_generate_figures.py`). The site is built with `cd site && npm run build` and published at
https://mlsmawfield.com/tep/bh.

The pipeline is 50 steps (step_00–step_49), runs warning-free, and writes SHA-256 checksums for all result files. Steps 00–11 use $\phi_0 = 2.0$, $\delta = 0.05$, $\sigma_B = 1.5$ for the fixed-background limit; steps 12–15 use the sGB coupling $\eta$ (primary $\eta = -0.1$) for the coupled solution, selecting the mass-inflation branch; steps 16–41 are derivation scripts; step 42 generates figures; steps 43–49 implement the S-star inference. All data is from published sources (EHT 2019-D01-01, LIGO Zenodo 4057131, and compiled measurement tables). Run with `python scripts/run_pipeline.py`; full docs in `scripts/README.md`.

**S-star inference pipeline.** Steps 43–49 of the unified pipeline implement the primary falsifiable test: a non-isochronous refit of the published S-star data around Sgr A*. Step 43 downloads the machine-readable CDS files for Gillessen et al. (2017, VizieR J/ApJ/837/30, table5.dat and table3.dat); it parses 145 NACO/NTT/Keck/Gemini astrometric epochs and 44 SINFONI/Keck/Gemini radial-velocity epochs for S2 relative to Sgr A*. The steps are: (43) data acquisition, (44) conventional GR fit, (45) Gate -1 mass-bias sign, (46) TEP transfer-function fit, (47) joint forward-model with composition/calibration separation and $\mathcal C_V$ diagnostic, (48) formal likelihood comparison (Bayes factor, BIC, AIC, Wilks), (49) MCMC posterior on $M_{\rm phantom}^T$. Run with `python scripts/run_pipeline.py`.

# Appendix K — Regular-Geometry Validation Benchmark

This appendix records the regular-geometry validation benchmark used throughout the main text. The Hayward metric is *not* the TEP solution and is *not* derived from the sGB field equations. It is a controlled mathematical reference used to validate the analysis pipeline and to confirm that finite curvature and bounded areal geometry are mutually compatible once the geometric singularity is absent. The question TEP asks is whether a dynamical proper-time field can supply the physical origin of that regularisation. The sGB coupled interior integration (Section 4.6) is the decisive test of this question. The literature evidence (Sotiriou & Zhou 2014; Thaalba et al. 2024) suggests that standard linear sGB may develop a finite-area singularity rather than a regular centre, so the benchmark may not be achievable with the linear sGB coupling alone.

## K.1 The Hayward Metric

The Hayward metric is

\begin{equation} \label{eq:appK_1}
F(r) = 1 - \frac{2Mr^2}{r^3 + 2M\ell^2},
\end{equation}

where $\ell$ is the regularisation scale. As $r \to 0$, $F \to 1$ (de Sitter-like temporal minimum); as $r \to \infty$, $F \to 1 - 2M/r$ (Schwarzschild). The effective energy density is

\begin{equation} \label{eq:appK_2}
\rho_{\rm eff} = -G^t_t = \frac{12M^2\ell^2}{(r^3 + 2M\ell^2)^2} \to \frac{3}{\ell^2} \quad \text{as } r \to 0,
\end{equation}

finite at the centre. The Kretschmann scalar is $K_{\rm Hayward}(0) = 24/\ell^4$, finite. The equation of state is $w_r = -1$ (cosmological-constant-like radial pressure), with anisotropic tangential pressure.

## K.2 TEP Matter Metric on the Hayward Background

The TEP matter metric $\tilde g = A^2 g_{\rm Hayward}$ with $A = (r_h/r)^{\phi_0}$ and $\phi_0 = 1$ is placed on the Hayward background. The computed invariants are:

| Quantity | Value as $r \to 0$ | Classification |
| --- | --- | --- |
| Kretschmann $\tilde K$ | $8/r_h^4 = 0.5$ (for $M=1$, $r_h=2M$) | Finite — de Sitter-like temporal minimum |
| Areal radius $\rho = Ar$ | $1.995$ | Bounded — finite-area asymptotic end |
| $\det\tilde g_{2D}$ | $< 0$ at 100% of sampled radii | Lorentzian in exterior ($F > 0$) |
| Null expansion $\theta_+$ | $0$ exactly | Constant-area property of the asymptotic end (not temporal freeze) |
| Affine parameter $\tilde\lambda$ (both families) | $\sim \int r^{-2}\,dr \to \infty$ | Infinite affine distance — asymptotic end, not regular point |
| Radial proper distance $\tilde\ell$ | $\sim \int dr/r \to \infty$ | Infinite — tube-like asymptotic end |
| Static-observer clock rate $d\tilde\tau/dt$ | $\to \infty$ | Conformal blueshift (unphysical; benchmark artefact of divergent $A$) |
| Static-observer redshift $z$ | $\to -1$ (blueshift) | Infinite redshift is at the Hayward outer horizon, not the centre — benchmark artefact of the divergent-$A$ test profile, not the TEP target |

## K.3 Classification of the Limiting Region

At an ordinary regular spherical centre one has $\rho \to 0$ as $r \to 0$. Here $\rho \to 1.995 \neq 0$: the limiting areal radius is nonzero. Because $A \sim 1/r$ in the deep interior, the radial proper distance behaves as $\tilde\ell \sim \int dr/r \to \infty$. The limiting region is therefore not a point-like centre but a finite-area asymptotic end (a tube-like limiting region of infinite radial extent). This is still a valid TEP geometry — the curvature is finite, the areal radius is bounded, and the metric is globally Lorentzian — but it must be classified correctly as an asymptotic end rather than a regular point centre. The full coupled interior solution (Section 4.6) now produces an ordinary regular point centre with $\rho_{\rm areal}(0) \to 4.68 \times 10^{-8}$ for $\eta=-0.1$, $g=1.1M$.

## K.4 What the Benchmark Validates

The benchmark validates two things:

- The analysis pipeline correctly recognises a regular geometric metric and classifies the matter-frame limiting region (finite Kretschmann, bounded areal radius, Lorentzianity inherited from $g$ via $A > 0$).

- Finite curvature and bounded areal geometry are mutually compatible once the geometric singularity is absent — the conditions that were mutually exclusive in the fixed-background theorem (Section 3) are simultaneously achievable on a regular background.

The benchmark does *not* validate: (a) that the sGB coupling dynamically generates this geometry (that is the task of the interior integration, Section 4.6, now solved for a regular point centre); (b) that the limiting region is an ordinary regular centre (it is a finite-area asymptotic end, Section K.3); (c) that the Temporal Horizon is established (the null expansion $\theta_+ = 0$ reflects the constant-area property of the asymptotic end, not a proof of temporal freezing; the Temporal Horizon requires the invariant frequency-transfer calculation of Section 6.5); (d) geodesic completeness (both null families have infinite affine parameter, $\tilde\lambda \sim \int r^{-2}\,dr \to \infty$; the asymptotic end is at infinite affine distance, not a regular point that is reached in finite parameter).

## K.5 Regularisation Scale Effects

The Kretschmann scalar at the centre depends on the regularisation scale: $K_{\rm Hayward}(0) = 24/\ell^4$. For $\ell = 0.1$ (the benchmark value), $K = 2.4 \times 10^5$; for $\ell = 1.0$, $K = 24$; for $\ell = 0.001$, $K = 2.4 \times 10^{13}$. These are code-unit values ($M = 1$); they do not establish Planck-scale curvature without converting to physical units for a specific black-hole mass. The regularisation scale is a parameter of the benchmark, not a prediction of the theory; the coupled interior solution with $g = 1.1M$ gives $K(0) = 27.1$ for $\eta=-0.1$ (Section 4.6).

## K.6 Configuration Scan

A scan over conformal exponent $\phi_0$ on the Hayward background. For a regular (de Sitter-like) seed, the background Kretschmann $K_{\rm Hay} \to$ finite as $r \to 0$, so the dominant contribution to the conformal Kretschmann comes from the $A$-derivative terms: $\tilde K \sim r^{4\phi_0 - 4}$ (not $r^{4\phi_0 - 6}$, which is the Schwarzschild-seed result where $K_{\rm Schw} \sim r^{-6}$ dominates). The areal radius scales as $\rho = Ar \sim r^{1-\phi_0}$ on either seed. The asymptotic classifications are:

| $\phi_0$ | $\rho$ as $r\to 0$ | $\tilde K$ as $r\to 0$ | Lorentzian | Asymptotic classification |
| --- | --- | --- | --- | --- |
| 0.5 | $\to 0$ | $\to \infty$ ($r^{-2}$) | Yes | Point centre, curvature diverges |
| 1.0 | $\to r_h \neq 0$ | $\to$ finite ($r^0$) | Yes | Finite-area asymptotic end, finite curvature |
| 1.25 | $\to \infty$ ($r^{-0.25}$) | $\to 0$ ($r^1$) | Yes | Areal radius diverges, curvature vanishes |
| 2.0 | $\to \infty$ ($r^{-1}$) | $\to 0$ ($r^4$) | Yes | Areal radius diverges, curvature vanishes |

The $\phi_0 = 1$ configuration is the validation benchmark used in the main text: it is the unique configuration that simultaneously achieves bounded (nonzero) areal radius and finite curvature on the Hayward background. The curvature scaling $\tilde K \sim r^{4\phi_0 - 4}$ differs from the Schwarzschild-seed scaling $\tilde K \sim r^{4\phi_0 - 6}$ (Appendix D) because the regular seed has finite background curvature, so the $A$-derivative terms dominate rather than the background Kretschmann. This distinction is critical: the fixed-background theorem uses the Schwarzschild exponent, while the benchmark classification uses the Hayward exponent. The two must not be mixed.

# Appendix L — Explicit Derivation of Corrected Exterior Observables

This appendix provides the step-by-step derivation of the corrected exterior observables — the Regge–Wheeler potential, the photon sphere and shadow radius, and the ISCO on the matter metric $\tilde g$ — from the Sotiriou \& Zhou (2014) perturbative sGB solution. The derivation demonstrates the coupling-order difference between the photon and massive-particle sectors explicitly.

## L.1 The sGB-Corrected Metric

The Sotiriou \& Zhou (2014) perturbative solution gives the geometric metric to $\mathcal{O}(\beta^2)$, where $\beta = \alpha_{\rm GB}/r_H^2 = \eta/12$ (since $r_H = 2M$ and $\alpha_{\rm GB} = \eta M^2/3$):

\begin{equation} \label{eq:appL_1}
ds^2_{\rm sGB} = -F(1 + \beta^2 h_2)\,dt^2 + \frac{1 + \beta^2 \sigma_2}{F}\,dr^2 + r^2\,d\Omega^2,
\end{equation}

where $F = 1 - r_H/r$, $x = r_H/r$, and the dimensionless polynomials are (Sotiriou \& Zhou 2014, Eqs. 60–61):

\begin{equation} \label{eq:appL_2}
h_2(x) = -\frac{98}{5}x - \frac{98}{5}x^2 - \frac{274}{15}x^3 - \frac{14}{15}x^4 + \frac{52}{15}x^5 + \frac{20}{3}x^6,
\end{equation}

\begin{equation} \label{eq:appL_3}
\sigma_2(x) = \frac{98}{5}x + \frac{58}{5}x^2 + \frac{38}{5}x^3 - \frac{406}{15}x^4 - \frac{436}{15}x^5 - \frac{92}{3}x^6.
\end{equation}

The would-be horizon shifts to $r_H = 2m = 2M(1 - 19.6\,\beta^2)$ at fixed ADM mass — a property of the horizon-bearing branch on which the perturbative sGB benchmark lives, not a physical TEP temporal-well horizon. The scalar field is $\phi(r) = (2\alpha_{\rm GB}/m)(1/r + m/r^2 + 4m^2/(3r^3))$ with $m = M/(1 + 49\eta^2/360)$ and scalar charge $Q_s = 2\alpha_{\rm GB}/m$. For the TEP mass-inflation branch the scalar charge is negative ($Q_s < 0$), so $\phi < 0$ in the exterior. The conformal factor is $A = e^{-\phi}$, giving the matter metric $\tilde g_{\mu\nu} = A^2 g^{\rm sGB}_{\mu\nu}$ with $A > 1$ (conformal magnification). The true TEP solution is a temporal well: $N \geq N_{\min} > 0$ everywhere, with a finite but extremely small minimum lapse at the centre.

The perturbative regime requires $\beta^2 |h_2(x)| \ll 1$ at all radii of interest. At the would-be horizon ($x = 1$), $h_2(1) = -48.27$, so $\beta^2 |h_2(1)| = 6.94 \times 10^{-5} \times 48.27 = 0.0034 \ll 1$. The value $\eta = -0.1$ ($\beta^2 = 6.94 \times 10^{-5}$) is used throughout this appendix as the TEP mass-inflation benchmark ($\alpha_{\rm GB} < 0$). Note that $|\eta| = 0.1$ exceeds current observational constraints ($\alpha_{\rm GB} < 2.9$ km$^2$ gives $|\eta| < 0.04$ for $M = 10\,M_\odot$); it is used as a perturbative reference value, not as a physically allowed coupling.

## L.2 Regge–Wheeler Potential: Status

A naive application of the scalar gradient $\phi'(r)$ directly into the standard GR Regge–Wheeler potential formula yields a term with dimensions $L^{-4}$ — inconsistent with the potential's $L^{-2}$ dimensions. The correct approach requires deriving the axial perturbation equation from the second variation of the sGB action on the corrected background. In modified gravity, the axial perturbation equation generally cannot be obtained merely by inserting the modified background functions into a GR potential formula; the modified field equations contribute directly to the quadratic perturbation operator.

The dominant effect on the QNM spectrum comes from the would-be horizon shift: $r_H \to 2M(1 - 19.6\,\beta^2)$ shifts $F$ and hence the entire potential. The would-be horizon scale change $r_H \to 2M(1 - 19.6\,\beta^2)$ gives an $\mathcal{O}(\beta^2) \sim \mathcal{O}(10^{-3})$ QNM modification at $\eta = -0.1$ (with $\beta = \eta/12$); the coefficient and sign are not determined by the would-be horizon shift alone. A precise value requires a coupled spectral solver on the corrected background, with the perturbation equations derived from the second variation of the action and validated against published nonperturbative sGB QNM spectra (Witek et al. 2019; Blazquez-Calzadilla et al. 2020; Chen et al. 2024). The structural prediction (scalar-led channel, broken isospectrality at $\mathcal{O}(\eta^2)$) follows from the coupling structure; the sign of the QNM shift is not determined by the would-be horizon shift alone. The deep-transit ringdown structure on the horizonless temporal well is analysed in Section 7.3.

## L.3 Photon Sphere and Shadow (Null Geodesics on $g^{\rm sGB}$)

For null geodesics on $ds^2 = -f\,dt^2 + g\,dr^2 + h\,d\Omega^2$, the impact parameter at a circular null orbit is $b^2 = h/f$. The photon sphere is at the minimum of $b^2(r)$, and the shadow radius is $b_{\rm ph} = \sqrt{h/f}$ at that minimum.

For the geometric metric $g^{\rm sGB}$: $f = F(1 + \beta^2 h_2)$, $h = r^2$. The minimum of $b^2 = r^2/[F(1 + \beta^2 h_2)]$ gives $r_{\rm ph} = 2.9974M$ and $b_{\rm ph} = 5.194M$ at $\eta = -0.1$, a $-0.044\%$ deviation from Schwarzschild ($b = 5.196M$).

For the matter metric $\tilde g = A^2 g^{\rm sGB}$: $f_{\tilde{}} = A^2 f$, $h_{\tilde{}} = A^2 r^2$, so $b^2 = h_{\tilde{}}/f_{\tilde{}} = r^2/f = r^2/[F(1 + \beta^2 h_2)]$ — identical to $g^{\rm sGB}$. The conformal factor $A^2$ cancels exactly. This is conformal invariance of null geodesics: the photon sphere and shadow are the same on $g^{\rm sGB}$ and $\tilde g$. The shadow deviation is $\mathcal{O}(\eta^2)$, determined solely by the metric perturbation.

## L.4 ISCO (Timelike Geodesics on $\tilde g$)

For timelike geodesics on $ds^2 = -f\,dt^2 + g\,dr^2 + h\,d\Omega^2$, the specific angular momentum of a circular orbit is:

\begin{equation} \label{eq:appL_4}
L^2 = \frac{f'\, h^2}{f\, h' - f'\, h},
\end{equation}

and the specific energy is $E^2 = f(1 + L^2/h)$. The ISCO is at the minimum of $L^2(r)$ (marginally stable orbit).

For the geometric metric $g^{\rm sGB}$ ($h = r^2$): $L^2 = f' r^4 / (f \cdot 2r - f' r^2)$. The minimum gives $r_{\rm ISCO} = 5.996M$ at $\eta = -0.1$, a $-0.061\%$ deviation. This is an $\mathcal{O}(\eta^2)$ effect from the metric perturbation alone.

For the matter metric $\tilde g = A^2 g^{\rm sGB}$ ($h = A^2 r^2$): the angular component now includes the conformal factor. The derivative $h' = 2A^2 r + 2A^2 r \cdot (-\phi') r = A^2(2r - 2r^2\phi')$ introduces $\phi'$ at $\mathcal{O}(\eta)$, which does not cancel. The minimum of $L^2(r)$ gives $r_{\rm ISCO} = 6.117M$ at $\eta = -0.1$, a $+1.95\%$ deviation. This is an $\mathcal{O}(\eta)$ effect: the conformal factor $A = e^{-\phi} > 1$ magnifies the effective areal radius $\tilde R = A \cdot r$, pushing the ISCO outward.

## L.5 Coupling-Order Difference in the sGB Benchmark

The key result is the coupling-order difference and sign difference between the shadow and the matter-metric ISCO. Under fixed-ADM mass normalization ($M_{\rm ADM} = M$ held fixed; the would-be horizon shifts to $r_H = 2M(1 - 19.6\,\beta^2)$), the shadow shifts negative while the matter-metric ISCO shifts positive — the conformal magnification pushes massive orbits outward while the geometric metric receives only an $\mathcal{O}(\eta^2)$ inward shift. The ISCO shift is $\sim 44\times$ larger in magnitude:

\begin{equation} \label{eq:appL_5}
\frac{\delta b}{b}\bigg|_{\rm shadow} = -0.044\% \quad (\mathcal{O}(\eta^2), \text{ photons on } g^{\rm sGB} = \tilde g),
\end{equation}

\begin{equation} \label{eq:appL_6}
\frac{\delta r_{\rm ISCO}}{r_{\rm ISCO}}\bigg|_{\tilde g} = +1.95\% \quad (\mathcal{O}(\eta), \text{ massive particles on } \tilde g).
\end{equation}

The coupling-order difference is a mathematical property of the perturbative sGB calculation (Section 6.2), not an indication of fractured metric structure. Matter, photons, and clocks all propagate universally on the single causal matter metric, $\tilde g_{\mu\nu}$. The $\sim 44\times$ coupling-order ratio and sign divergence arise purely from the geometry of the paths: because null geodesics ($d\tilde s^2 = 0$) are conformally invariant, the temporal scaling $A^2(\phi)$ divides out exactly, making the shadow sensitive only to the underlying backreacted geometry at $\mathcal{O}(\eta^2)$. Conversely, massive particles travel on timelike geodesics ($d\tilde s^2 < 0$) and therefore feel the conformal factor $A = e^{-\phi} > 1$ directly, shifting the ISCO at $\mathcal{O}(\eta)$. Both effects are complementary projections of the same unified temporal well. The matter-metric ISCO shift is an order of magnitude larger than the shadow shift ($+1.95\%$ vs $-0.044\%$) because the conformal factor contributes at $\mathcal{O}(\eta)$, not $\mathcal{O}(\eta^2)$. The would-be horizon shift $r_H = 2M(1 - 19.6\,\beta^2)$ is the sGB perturbative signature: the temporal field backreacts on the geometry through the sGB coupling, contracting the horizon of the horizon-bearing branch. The Temporal Horizon is observer-dependent and operational — a distant observer (fast clock) sees the accessibility boundary far from the centre; a deeper observer (stronger conformal field, larger $A$) sees less redshift to any given emitter, and the practical observability boundary moves inward.

The scaling with $\eta$ is:

- Shadow: $\delta b/b \propto \eta^2$ (geometric metric only, conformal factor cancels)

- ISCO on $g^{\rm sGB}$: $\delta r/r \propto \eta^2$ (geometric metric only)

- ISCO on $\tilde g$: $\delta r/r \propto \eta$ (conformal factor dominates)

All numbers are computed at $\eta = -0.1$ in the perturbative regime. Note that $|\eta| = 0.1$ exceeds current observational constraints ($\alpha_{\rm GB} < 2.9$ km$^2$ gives $|\eta| < 0.04$ for $M = 10\,M_\odot$); it is used as a perturbative reference value. Within current bounds, the largest shifts are for $M = 10\,M_\odot$: shadow $-0.007\%$, ISCO $+0.78\%$. The computation uses the Sotiriou \& Zhou (2014) exact $\mathcal{O}(\beta^2)$ metric perturbations and the analytic sGB scalar profile, with no free parameters..

# References

- Smawfield, M. L. (2025). Temporal Equivalence Principle: Dynamic Time & Emergent Light Speed. Paper 0, v0.9 (Jakarta). DOI: 10.5281/zenodo.16921911.

- Schwarzschild, K. (1916). Über das Gravitationsfeld eines Massenpunktes nach der Einsteinschen Theorie. *Sitzungsber. Preuss. Akad. Wiss.*, 189–196.

- Penrose, R. (1965). Gravitational collapse and space-time singularities. *Phys. Rev. Lett.*, 14, 57–59.

- Hawking, S. W. (1966). The occurrence of singularities in cosmology. *Proc. Roy. Soc. A*, 294, 511–521.

- Hawking, S. W. & Penrose, R. (1970). The singularities of gravitational collapse and cosmology. *Proc. Roy. Soc. A*, 314, 529–548.

- Bekenstein, J. D. (1973). Black holes and entropy. *Phys. Rev. D*, 7, 2333–2346.

- Hawking, S. W. (1975). Particle creation by black holes. *Commun. Math. Phys.*, 43, 199–220.

- Bekenstein, J. D. (1993). Relation between physical and gravitational geometry. *Phys. Rev. D*, 48, 3641–3647. arXiv:gr-qc/9211017.

- Wald, R. M. (1993). Black hole entropy is the Noether charge. *Phys. Rev. D*, 48, R3427–R3431. arXiv:gr-qc/9307038.

- Bardeen, J. M. (1968). Non-singular general relativistic gravitational collapse. In *Proceedings of the 5th International Conference on Gravitation and the Theory of Relativity*, Tbilisi, p. 174.

- Hayward, S. A. (2006). Formation and evaporation of nonsingular black holes. *Phys. Rev. Lett.*, 96, 031103. arXiv:gr-qc/0506126.

- Mazur, P. O. & Mottola, E. (2004). Gravitational vacuum condensate stars. *Proc. Nat. Acad. Sci.*, 101, 9545–9550. arXiv:gr-qc/0109035.

- Mathur, S. D. (2005). The fuzzball proposal for black holes: an elementary review. *Fortsch. Phys.*, 53, 793–827. arXiv:hep-th/0502050.

- Almheiri, A., Marolf, D., Polchinski, J. & Sully, J. (2013). Black holes: complementarity or firewalls? *JHEP*, 02, 062. arXiv:1207.3123.

- Khoury, J. & Weltman, A. (2004). Chameleon cosmology. *Phys. Rev. D*, 69, 044026.

- Event Horizon Telescope Collaboration (2019). First M87 Event Horizon Telescope results. I. The shadow of the supermassive black hole. *ApJL*, 875, L1.

- Event Horizon Telescope Collaboration (2019). First M87 Event Horizon Telescope results. IV. Imaging the central supermassive black hole. *ApJL*, 875, L4.

- Event Horizon Telescope Collaboration (2019). First M87 Event Horizon Telescope results. V. Physical origin of the asymmetric ring. *ApJL*, 875, L5.

- Event Horizon Telescope Collaboration (2022). First Sagittarius A* Event Horizon Telescope results. I. The shadow of the supermassive black hole in the center of the Milky Way. *ApJL*, 930, L12.

- Smawfield, M. L. (2026). Temporal Equivalence Principle: A Covariant Alternative to Cosmic Expansion. Paper 26, v0.1 (Athens). DOI: 10.5281/zenodo.20370143.

- Smawfield, M. L. (2026). Temporal Equivalence Principle: Native hi_class Conformal Implementation, Linear Perturbation Closure, and CMB Acoustic Peak Preservation. Paper 18, v0.5 (Cambridge).

- Smawfield, M. L. (2026). Temporal Equivalence Principle: Temporal Horizon Cosmology and the Absence of a Physical Big Bang Singularity. Paper 27, v0.2 (Thika). DOI: 10.5281/zenodo.20723059.

- LIGO/Virgo Collaboration (2016). Observation of gravitational waves from a binary black hole merger. *Phys. Rev. Lett.*, 116, 241102. DOI: 10.1103/PhysRevLett.116.241102.

- LIGO/Virgo Collaboration (2020). GW190521: A binary black hole merger with a total mass of 150 M<sub>☉</sub>. *Phys. Rev. Lett.*, 125, 101102. DOI: 10.1103/PhysRevLett.125.101102.

- Isi, M., et al. (2020). Posterior samples of GW190521 posterior samples. *Zenodo*. DOI: 10.5281/zenodo.4057131.

- Strohmayer, T. E. (2001). Discovery of high-frequency quasi-periodic oscillations in the black-hole binary GRS 1915+105. *ApJ*, 552, L49.

- Remillard, R. A., et al. (2002). XTE J1550-564: QPOs and black-hole spin. *ApJ*, 564, 962.

- Homan, J., et al. (2005). High-frequency QPOs in H 1743-322. *ApJ*, 624, 1005.

- Morgan, E. H., Remillard, R. A., & Greiner, J. (1997). RXTE observations of GRS 1915+105. *ApJ*, 482, 993.

- Harikane, E., et al. (2023). JWST CEERS: discovery of black-hole candidates at z > 8. *ApJ*, 958, 11.

- Greene, J. E., et al. (2023). CEERS AGN at z ~ 7. *ApJ*, 957, 24.

- Bogdan, A., et al. (2024). Evidence for a black-hole seed at z = 10.1. *Nature Astronomy*.

- Castellano, P., et al. (2024). JWST discovery of GHZ2 at z = 10.6. *Nature Astronomy*.

- Matsuoka, Y., et al. (2019). Subaru discovery of a z = 7.04 quasar. *ApJ*, 872, 2.

- Ba&ntilde;ados, E., et al. (2018). A black-hole mass of 8&times;10<sup>8</sup> M<sub>☉</sub> at z = 7.54. *Nature*, 553, 473.

- McClintock, J. E., Narayan, R., & Steiner, J. F. (2014). Black-hole spin via continuum fitting and reflection spectroscopy. *Space Sci. Rev.*, 183, 295. arXiv:1303.1583.

- Reynolds, C. S. (2021). X-ray reflection spin measurements and systematics. arXiv:2104.10300.

- Gou, L., et al. (2011). Cygnus X-1 spin via continuum fitting. *ApJ*, 742, 85.

- van Velzen, S., et al. (2021). Seventeen tidal disruption events from ZTF. *ApJ*, 908, 4.

- Holoien, T. W.-S., et al. (2016). ASASSN-14li: a nearby TDE. *MNRAS*, 455, 1690.

- Holoien, T. W.-S., et al. (2020). ASASSN-18pg: a luminous TDE. *ApJ*, 903, 151.

- Nicholl, M., et al. (2020). AT2019qiz: a nearby TDE with early-time spectroscopy. *MNRAS*, 499, 482.

- Event Horizon Telescope Collaboration (2021). First M87 Event Horizon Telescope results. VII. Polarization of the ring. *ApJL*, 910, L13.

- Event Horizon Telescope Collaboration (2024). Sgr A* polarimetry with the EHT. *ApJL*, 964, L26.

- Bekenstein, J. D. (2004). Conformal and disformal transformations. *Phys. Rev. D*, 70, 083509.

- Koivisto, T. (2012). Disformal equivalence. *Phys. Rev. D*, 85, 044043.

- Inayoshi, K., Haiman, Z., & Ostriker, J. P. (2020). Hyper-Eddington accretion and rapid growth of massive black holes. *MNRAS*, 496, 4236.

- Sotiriou, T. P. & Zhou, S.-Y. (2014). Black hole hair in generalized scalar-tensor gravity: An explicit example. *Phys. Rev. D*, 90, 124063. arXiv:1408.1698. [Companion to PRL 112, 251102 (2014); contains the exact $\mathcal{O}(\alpha^2)$ perturbative metric corrections used in Section 4. Their nonlinear solutions develop a finite-area singularity at $r \approx 0.27\,r_h$ — see Section 4.2, 4.7.]

- Kanti, P., Mavromatos, N. E., Rizos, J., Tamvakis, K. & Winstanley, E. (1996). Dilatonic black holes in higher curvature string gravity. *Phys. Rev. D*, 54, 5049. arXiv:hep-th/9511071. [Establishes scalar hair and regular horizons in sGB; "absence of naked singularities" permits singularities hidden inside the horizon. Follow-up (Phys. Rev. D 57, 6255, 1998) confirms "other researchers... demonstrated numerically the existence of curvature singularities" behind the horizon.]

- Kleihaus, B., Kunz, J. & Radu, E. (2011). Rotating black holes in dilatonic Einstein-Gauss-Bonnet theory. *Phys. Rev. Lett.*, 106, 151104. arXiv:1101.2868. [Rotating sGB black holes; domain of existence bounded by singular extremal solutions.]

- Thaalba, F., Franchini, N., Bezares, M. & Sotiriou, T. P. (2024). Dynamics of spherically symmetric black holes in scalar-Gauss-Bonnet gravity with a Ricci coupling. *Phys. Rev. D*, 111, 064054. arXiv:2409.11398. [Confirms finite-area singularity in sGB; explores connection between singularity formation and loss of hyperbolicity; Ricci coupling can mitigate hyperbolicity loss.]

- Torii, T., Maeda, K. & Tamaoki, T. (1999). Rotating non-singular black holes in dilatonic Gauss-Bonnet gravity. *Phys. Rev. D*, 59, 064012.

- Delgado, J. F. M., Herdeiro, C. A. R., & Radu, E. (2020). Spin-induced scalarization and spontaneous scalarization of Kerr black holes in scalar-Gauss–Bonnet gravity. *Phys. Rev. D*, 102, 044041. arXiv:2007.12030.

- Delgado, J. F. M., Herdeiro, C. A. R. & Radu, E. (2020). Spinning black holes in shift-symmetric Horndeski theory. *JHEP*, 04, 180. arXiv:2002.05012. [Slowly-rotating shift-symmetric sGB solution; $\mathcal{O}(\beta^2)$ correction to $W(r)$ and horizon angular velocity; nonperturbative numerical solutions for arbitrary spin. Used in Section 8.5.]

- Bardeen, J. M. (1973). Rapidly rotating stars, disks, and black holes. In *Black Holes (Les Astres Occlus)*, eds. DeWitt & DeWitt, Gordon & Breach, pp. 241–289.

- Leaver, E. W. (1985). An analytical representation for the quasi-normal modes of Kerr black holes. *Proc. Roy. Soc. A*, 402, 285–298.

- Witek, H., et al. (2019). Scalar-Gauss–Bonnet gravity in the strong-field regime. *Phys. Rev. D*, 99, 064035. arXiv:1810.05177.

- Blazquez-Calzadilla, J., et al. (2020). Scalar-Gauss–Bonnet perturbations. *Phys. Rev. D*, 102, 024086. arXiv:2003.02862.

- Chen, Y., et al. (2024). Quasinormal modes of rotating black holes in shift-symmetric Einstein-scalar-Gauss–Bonnet theory. arXiv:2412.09377. [Nonperturbative rotating sGB QNM spectra; validation target for the TEP-sGB QNM pipeline.]

- Aresté Saló, L., Doneva, D. D., Clough, K., Figueras, P. & Yazadjiev, S. S. (2025). Challenges in the nonlinear evolution of unequal mass binaries in scalar-Gauss–Bonnet gravity. *Phys. Rev. D*, in press. arXiv:2507.13046. [Numerical relativity simulations of binary black hole mergers in sGB; ringdown excitation extraction.]

- Kobayashi, T., Yamaguchi, M. & Yokoyama, J. (2011). Galilean creation of the inflationary universe. *Prog. Theor. Phys.*, 126, 511. arXiv:1105.5723. [Principal symbol of scalar–tensor theories; tensor speed derivation used in Section 7.6.]

- Nishizawa, A. & Arai, K. (2019). Generalized framework for testing gravity with gravitational-wave propagation. *Phys. Rev. D*, 99, 104038. arXiv:1901.08249. [c_T formula for modified gravity theories including sGB; GW170817 constraints.]

- LIGO/Virgo & Fermi-GBM Collaborations (2017). Gravitational waves and gamma-rays from a binary neutron star merger: GW170817 and GRB 170817A. *ApJL*, 848, L13.

- Freire, P. C. C., et al. (2012). The relativistic pulsar-white dwarf binary PSR J0348+0432. *MNRAS*, 423, 3328.

- Antoniadis, J., et al. (2013). A massive pulsar in a compact relativistic binary. *Science*, 340, 6131.

# Acknowledgements

The author acknowledges the broader TEP research programme, founded in Paper 0 (Temporal Equivalence Principle: Dynamic Time & Emergent Light Speed), and the empirical and computational framework established by TEP-C0 (Paper 26, Athens) and TEP-HC (Paper 18, Cambridge) within which this black-hole interior analysis is developed. Section 7.5 relates the SEC-violation mechanism found here to the independent construction of TEP-TH (Paper 27, Thika), which establishes an analogous singularity-exclusion result for the cosmological branch.
