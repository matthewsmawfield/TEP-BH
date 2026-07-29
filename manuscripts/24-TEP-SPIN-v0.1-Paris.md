# Temporal Equivalence Principle: Fermion Spin as Temporal-Orientation Holonomy
**Matthew Lukin Smawfield**
Version: v0.1 (Paris)
First published: 20 July 2026 · Last updated: 26 July 2026
DOI: 10.5281/zenodo.20572705

---

## Abstract

The Temporal Equivalence Principle (TEP) proposes replacing the zero-dimensional point-particle paradigm with a finite Compton-scale topological defect in dynamical proper time. We prove that the conformal temporal shear &Sigma;<sub>&mu;</sub> = &nabla;<sub>&mu;</sub> ln A(&phi;) is an exact, irrotational one-form, so its circulation vanishes on contractible loops. This geometric no-go result forbids spin holonomy within the scalar conformal sector and motivates a compact temporal-orientation bundle with structure group Spin(3) &cong; SU(2). The orientation field U(**x**) &isin; SU(2) is assigned the hedgehog ansatz, classified by &pi;<sub>3</sub>(SU(2)) = &Zopf; with topological charge B = 1. The U(1) winding n = &pm;1 of Section 2 is the fixed-axis projection of B. The defect is stabilized by a Skyrme quartic derivative term; without it, Derrick's theorem forbids a static localized soliton in 3+1 dimensions. A numerical collocation solver finds a convergent finite-energy B = 1 profile with E(B=2)/E(B=1) = 2.82 > 2, indicating higher charges are energetically disfavoured. The spin magnitude J = 1/2 is derived from collective-coordinate quantization on the SU(2) configuration space: SU(2) is simply connected (&pi;<sub>1</sub> = 0) and admits half-integer representations that SO(3) does not. The minimal nonzero excitation has j = 1/2, giving J<sub>z</sub> = &pm;&hbar;/2. The factor 1/2 emerges from the topology of SU(2), not from importing S = Q/2. The 2&pi; sign reversal and 4&pi; return are emergent consequences. Internal orientation is connected to physical rotation through diagonal soldering (Spin(3)<sub>spatial</sub> &times; Spin(3)<sub>orient</sub> &rarr; Spin(3)<sub>diagonal</sub>), with the Lorentz spin connection, orientation connection, and electromagnetic gauge potential maintained as three distinct connections. The core radius is identified with the electron Compton scale r<sub>c</sub> = &hbar;/(m<sub>e</sub>c) (scale identification, not prediction). Two conditional empirical tests are presented: an AMBER proton form-factor forecast (&delta;<sub>A</sub> = 0.0028 &plusmn; 0.0039, 0.72&sigma; from zero) and a disformal (g&minus;2) candidate operator whose normalization is under active derivation. The Skyrme coupling e<sub>s</sub>, the moment of inertia &Lambda;, and the complete environmental screening closure remain open.

Keywords: subatomic structure, fermion topology, spin, vorticity, temporal-orientation holonomy, proximity screening, Fermilab g-2, AMBER, temporal equivalence principle

## 1. Introduction: Spin as Temporal-Orientation Holonomy

### 1.1 The Conformal No-Go Theorem

Standard QFT models elementary fermions as pointlike field excitations. Ultraviolet divergences arise from local quantum-field interactions and short-distance behaviour, and renormalization provides a phenomenologically successful procedure for extracting finite predictions. TEP investigates whether a finite internal topological scale can provide a physical ultraviolet completion while reproducing the renormalized predictions of QED.

**Evidence status.** This paper derives the conformal no-go result and constructs a minimal compact orientation-sector completion. The finite-core regulator, disformal (g&minus;2) operator and many-body screening architecture are developed as candidate consequences whose complete gauge-covariant normalization and nonlinear closure remain under derivation. The bundle is introduced as the minimal TEP completion of the conformal scalar sector: since &Sigma;<sub>&mu;</sub> = &nabla;<sub>&mu;</sub> ln A is exact and irrotational, spin holonomy cannot reside in the scalar conformal sector and is carried by a compact orientation fibre, motivated by Proposition 1 from axioms A1–A3, the experimental fact of spinorial holonomy, and the TEP ontological principle. The bundle's kinetic and Stueckelberg terms are specified by the minimal exterior (London-limit) orientation-bundle action (Eq. 4a). The orientation-bundle gauge invariance, candidate Ward identity construction, and g&minus;2 candidate operator are addressed in Appendices A.17–A.19. Two precision empirical tests — the AMBER proton form-factor extraction and the Fermilab g&minus;2 temporal modulation — provide falsifiable cross-sections of the topological fermion geometry, conditional on the respective model assumptions.

Within the Temporal Equivalence Principle (TEP), we observe a structural property that shows why TEP's real conformal scalar sector cannot itself encode spin holonomy. The conformal temporal shear &Sigma;<sub>&mu;</sub> = &nabla;<sub>&mu;</sub> ln A(&phi;) is an exact one-form. For any contractible closed loop C,

&oint;<sub>C</sub> &Sigma;<sub>&mu;</sub> dx<sup>&mu;</sup> = &oint;<sub>C</sub> d ln A = 0.
(1)

Therefore the conformal scalar sector cannot carry spin holonomy. Under TEP's geometric-ontology premise, spin requires an additional compact temporal-orientation sector. This is not an arbitrary extra structure but is motivated by the single-valuedness of the real conformal factor and the TEP ontological principle that spin is a geometric property.

This structural observation motivates a compact temporal-orientation bundle, where integer winding lifts spinorially to the observed 2&pi; sign reversal and 4&pi; return. The fermion is thereby proposed as a finite Compton-scale topological defect in dynamical proper time. The vanishing circulation of any exact one-form on contractible loops is a standard result in differential geometry (Poincaré lemma). The compact temporal-orientation bundle is introduced as the minimal TEP completion, motivated by axioms A1–A3, the experimental fact of spinorial holonomy, and the TEP ontological principle, as discussed in Proposition 1 (Section 2.2).

At the quantum scale, the geometric proximity ratio &chi; = r<sub>c</sub>/&lambda;<sub>scr</sub> supplies the microscopic length required for a candidate nonlocal regulator. The finite closure scale is a property of the proposed geometry; it does not by itself prove that the complete interacting quantum theory is finite.

### 1.2 The TEP Topological Fermion

The Temporal Equivalence Principle replaces the point particle with a localized topological charge in the temporal shear field. The fermion is not a mathematical singularity but a physical defect in the scalar field &phi; that defines local proper time. The fermion possesses a finite Compton-scale core (r<sub>c</sub> = &hbar;/(m<sub>e</sub>c)). The relationship between this microscopic core and the macroscopic phenomenological saturation scale (&rho;<sub>T</sub> &asymp; 20 g/cm<sup>3</sup>, established from terrestrial clock correlation data in TEP-UCD, Paper 6) is investigated through the many-body screening architecture (Section 4) but is not yet closed quantitatively.

The TEP framework is built on three axioms. (A1) The matter-frame metric is a conformal–disformal rescaling of the gravitational metric: g&#771;<sub>&mu;&nu;</sub> = A<sup>2</sup>(&phi;) g<sub>&mu;&nu;</sub> + B(&phi;) &nabla;<sub>&mu;</sub>&phi; &nabla;<sub>&nu;</sub>&phi;. (A2) The conformal factor is exponential: A(&phi;) = exp(&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>). (A3) Temporal shear is the gradient of the logarithmic conformal factor: &Sigma;<sub>&mu;</sub> = &nabla;<sub>&mu;</sub> ln A(&phi;). The conformal Hamilton–Jacobi sector follows from these axioms; the spin/vorticity sector follows from the compact temporal-orientation bundle, introduced as Proposition 1 from the same axioms plus the experimental fact of spinorial holonomy and the TEP ontological principle (Section 2.2). The full causal matter metric is permanently engaged; the screened limit, where local stress forces both A(&phi;) &rarr; 1 and the observable disformal response is suppressed, recovers the standard Minkowski background and isotropic interactions. In the unscreened regime the disformal sector governs the routing of forces through the tilted light cone, as developed in TEP-KIN (Paper 25).

### 1.3 Context: Renormalization and the Point-Particle Paradigm

Early attempts at a unified geometric theory sought to replace dimensionless point particles with physical &ldquo;knots&rdquo; in spatial geometry, but failed to eliminate the mathematical divergences that necessitated renormalization. The TEP framework achieves this geometric origin by shifting the topology from spatial gravity to proper time. Instead of a point particle, TEP introduces a localized topological charge embedded within the temporal shear field. Because this fermion is a physical defect in the scalar field &phi;, it carries a natural geometric boundary. The finite Compton-scale core supplies a microscopic geometric boundary and motivates a candidate nonlocal regulator; its relationship to macroscopic environmental saturation remains under investigation. The regulator must be distinguished from an ordinary spatial charge radius: TEP's defect is a proper-time/orientation topology, and the regulator may act on internal off-shell propagation (R<sub>loop</sub>(k<sup>2</sup>)) distinct from external scattering form factors (F<sub>scattering</sub>(q<sup>2</sup>)); they cannot be assumed identical. Gauge invariance of the orientation-bundle action is ensured by the U(1) bundle construction (Proposition 1); identification of the orientation connection &omega;<sub>&mu;</sub> with the electromagnetic potential is an additional model assumption, not established by the bundle action alone. A gauge-covariant Wilson-line formulation (Appendix A.18) outlines a possible gauge-covariant completion; detailed diagnostics of the regulator, Ward identity, and running-coupling analysis are documented in the repository audit (`scripts/audit/`).

The electron self-energy diverges as &Lambda; &rarr; &infin; in the standard formulation, and the Landau pole in QED signals that the theory is incomplete at short distances. Renormalization provides a phenomenologically successful procedure for extracting finite predictions; a finite geometric structure motivates a candidate regulator that could provide a physical ultraviolet completion at the divergence's origin.

## 2. The Topological Fermion

### 2.0 Symbol Table

| Symbol | Definition | Value | Classification |
| --- | --- | --- | --- |
| *A*(&phi;) | Conformal factor | exp(&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>) | Fundamental |
| &beta;<sub>A</sub> | Fundamental conformal coupling | +1.0 | Fundamental |
| &beta;<sub>spin</sub> | Phenomenological screening coefficient | 0.01 | Effective (tanh ansatz) |
| &beta;<sub>A</sub><sup>(eff)</sup> | Effective lab coupling after screening | ~10<sup>&minus;2</sup> | Derived via *S*<sub>&Sigma;</sub>(*&Epsilon;*) |
| &chi; | Geometric closure ratio r<sub>c</sub>/&lambda;<sub>scr</sub> | 1/&radic;2 &approx; 0.707 | Model definition |
| &chi;<sub>q</sub> | Single-quark closure ratio | 1/&radic;2 (assumed, universal scalar-sector property) | Model definition |
| m<sub>&phi;</sub> | Scalar field mass | m<sub>e</sub>/&radic;2 | Model definition |
| &rho;<sub>T</sub> | Phenomenological saturation scale | ~20 g/cm<sup>3</sup> | Empirical input (Paper 6) |
| &chi;<sub>p</sub> | Proton closure ratio | ~0.21 | Bound-state assumption |

### 2.1 The Fermion as Topological Charge

The fermion is defined as a localized topological charge in the temporal shear field. In the matter frame, proper time d&tau; is set by the causal matter metric g&#771;<sub>&mu;&nu;</sub> = A<sup>2</sup>(&phi;) g<sub>&mu;&nu;</sub> + B(&phi;) &nabla;<sub>&mu;</sub>&phi; &nabla;<sub>&nu;</sub>&phi;. *Metric signature:* (+, &minus;, &minus;, &minus;). In the conformal limit relevant for the single-particle core geometry, a particle of mass m propagates according to the g&#771;-Hamilton-Jacobi equation, which in flat background with A(&phi;) = exp(&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>) reads:

g<sup>&mu;&nu;</sup> &part;<sub>&mu;</sub>S &part;<sub>&nu;</sub>S = m<sup>2</sup>c<sup>2</sup> exp(2&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>).
(2)

The effective mass in the matter frame is m<sub>*</sub> = m A(&phi;) = m exp(&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>). The proper-time oscillator frequency is therefore shifted by the local conformal factor, and the fermion acquires a position-dependent effective inertia. In the rest frame (&nabla;S = 0), the frequency is &omega;<sub>eff</sub> = mc<sup>2</sup> A(&phi;) / &hbar;. This is the origin of the bounded proper-time oscillator: as the topological charge tightens, A(&phi;) flattens toward unity and &omega;<sub>eff</sub> approaches the standard Compton frequency, but it never diverges because the core has finite extent. The conformally shifted g&#771;-Hamilton-Jacobi equation and the emergence of the Klein-Gordon and Dirac operators in the screened limit are derived systematically in TEP-QF (Paper 23).

### 2.2 Spin as Quantized Vorticity

"Spin" is translated directly into fluid vorticity. However, to preserve the single-valuedness of the real conformal factor A(&phi;) = exp(&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>), the fermion core must be mathematically modeled as a dual-component spatio-temporal vortex.

The temporal shear &Sigma;<sub>&mu;</sub> = &nabla;<sub>&mu;</sub> ln A(&phi;) is an exact one-form and is therefore strictly irrotational (&nabla; &times; &Sigma; = 0). This radial conformal shear generates the effective mass and the topographic drag, but it cannot carry spin.

**Structural Observation 1 (Conformal exactness).** Because A(&phi;) is a real, strictly positive, single-valued scalar field, &Sigma;<sub>&mu;</sub> = &nabla;<sub>&mu;</sub> ln A is an exact gradient. For any contractible closed loop C,

&oint;<sub>C</sub> &Sigma;<sub>&mu;</sub> dx<sup>&mu;</sup> = &oint;<sub>C</sub> d ln A = 0.
(3)

Within the minimal TEP completion developed here, spin holonomy is assigned to the compact temporal-orientation bundle rather than the irrotational conformal scalar sector. The azimuthal orientation shear is K<sub>&theta;</sub> = n/r, where n is the integer winding number. The circulation around the charge core is quantized within this compact orientation bundle:

&Gamma; = &oint; K &middot; d&ell; = 2&pi;n.
(4)

The vorticity vector &omega;<sub>i</sub> = (&nabla; &times; K)<sub>i</sub> vanishes everywhere except at the core singularity, where it is distributional. The minimal non-trivial defect carries integer winding n = &pm;1 in the compact orientation sector. As established in TEP-QF, the temporal-orientation bundle governing its causal frame is spinorial: a 2&pi; circuit reverses the spinor sign and a 4&pi; circuit restores it. The spin magnitude S<sub>z</sub> = &pm;&hbar;/2 is derived in Section 3 from the collective-coordinate quantization of the SU(2) hedgehog soliton, where the factor 1/2 emerges from the topology of the configuration space SU(2) (simply connected double cover of SO(3)), not from importing the spinorial representation.

**Proposition 1 (Necessity of a non-scalar holonomy sector).** The compact temporal-orientation bundle is introduced as the minimal TEP completion of the conformal scalar sector, motivated by the following chain of reasoning:

(P1) Axioms A1–A3: the conformal shear &Sigma;<sub>&mu;</sub> = &nabla;<sub>&mu;</sub> ln A is an exact one-form.

(P2) Experimental fact: fermion spin exhibits 2&pi; sign reversal and 4&pi; return (spinorial holonomy).

(P3) TEP ontological principle: spin is a geometric property of the particle, not an external label. This is a TEP-specific philosophical commitment, not an axiom or experimental fact.

From (P1), Structural Observation 1 proves &oint;<sub>C</sub> &Sigma;<sub>&mu;</sub> dx<sup>&mu;</sup> = 0 for every contractible loop C. The conformal sector therefore carries zero holonomy. From (P2), the physical spinor acquires a minus sign under 2&pi; rotation, which is a non-trivial holonomy in a compact phase variable. From (P3), this holonomy must be realized geometrically within the TEP field content. Since the only field in the conformal sector (&phi;) is a real single-valued scalar whose gradient is exact, it cannot supply the required compact phase. The holonomy must therefore reside in a sector orthogonal to the conformal scalar. Axioms A1–A3 exclude the real conformal scalar as the carrier of spin holonomy; a compact orientation sector is therefore introduced as the minimal TEP completion. Its choice of structure group, dynamics, and relation to electromagnetic U(1) are additional model assumptions, not logical consequences of A1–A3 alone.

The minimal such sector is a compact U(1) phase bundle with compact phase coordinate &theta; over spacetime, with connection one-form &omega;<sub>&mu;</sub> and curvature F<sub>&mu;&nu;</sub> = &part;<sub>[&mu;</sub>&omega;<sub>&nu;]</sub>. The compactness of the U(1) fibre is motivated by the observed periodicity: the spinor returns to itself after 4&pi;, so the phase variable is circle-valued with period 4&pi; in the spinorial representation (equivalently 2&pi; in the vector representation). The U(1) character is chosen as the unique connected compact one-dimensional Lie group. A compact U(1) phase can represent rotations about one selected axis; general spatial rotations are noncommuting and require the full Spin(3) &cong; SU(2) group. The present construction embeds a U(1) subgroup into SU(2) to describe a chosen spin projection. The spinorial lift (2&pi; sign reversal, 4&pi; return) is realized by this embedding. The complete noncommutative Spin(3) structure — including the SU(2)-valued orientation field, hedgehog topology, Skyrme stabilization, and collective-coordinate quantization yielding J = 1/2 — is developed in Section 3.

The minimal exterior (London-limit) orientation-bundle action is:

S<sub>bundle</sub> = &int; d<sup>4</sup>x &radic;&minus;g [ &minus;&#188; F<sub>&mu;&nu;</sub>F<sup>&mu;&nu;</sup> &minus; &#189; &kappa;<sup>2</sup> (&nabla;<sub>&mu;</sub>&theta; &minus; &omega;<sub>&mu;</sub>)<sup>2</sup> ],
(4a)

where &kappa; is the bundle stiffness (a model parameter fixed by the topological charge core energy) and the covariant derivative D<sub>&mu;</sub>&theta; = &nabla;<sub>&mu;</sub>&theta; &minus; &omega;<sub>&mu;</sub> ensures gauge invariance under &theta; &map; &theta; + &alpha;(x), &omega;<sub>&mu;</sub> &map; &omega;<sub>&mu;</sub> + &nabla;<sub>&mu;</sub>&alpha;. This is an exterior or London-limit action: it specifies the gauge-invariant dynamics away from the defect core. The complete finite-core action, including the SU(2)-valued orientation field, Skyrme stabilization term, radial amplitude field, and core potential, is presented in Section 3.1 (Eq. 22). The Bianchi identity dF = 0 follows identically from F = d&omega; on the punctured manifold (away from defect cores).

**Scope.** The spinorial lift (2&pi; sign reversal, 4&pi; return) is realized by embedding the U(1) orientation phase into the SU(2) spinor representation. The geometric arguments of Sections 2.1–2.3 depend on Proposition 1 (necessity of a non-scalar sector), the experimental fact of spinorial holonomy (P2), and the TEP ontological principle (P3), not on axioms A1–A3 alone. The bundle action (4a), its identification with electromagnetism, and the spinorial lift are model assumptions motivated by — but not uniquely forced by — the conformal no-go result.

**Topological dimensionality remark.** The U(1) winding described above is the fixed-axis projection of a full SU(2) temporal-orientation field. A spatial U(1) vortex in three dimensions would be a codimension-two line defect (&pi;<sub>2</sub>(S<sup>1</sup>) = 0), not a localized particle. The complete construction uses an SU(2)-valued orientation field U(**x**) in the hedgehog configuration, classified by &pi;<sub>3</sub>(SU(2)) = &Zopf;, which supports genuine localized point defects in three spatial dimensions. The U(1) winding n is the projection of the SU(2) topological charge B onto a chosen spin quantization axis. The full construction, including the Skyrme stabilization term, finite-core boundary-value solution, and collective-coordinate quantization, is developed in Section 3.

### 2.3 Charge Core Geometry and the &chi; Ratio

The core radius of the topological charge is identified with the electron Compton wavelength,

r<sub>c</sub> = &hbar;/(m<sub>e</sub>c),
(5)

which fixes the finite-core scale of the defect. The scalar Yukawa screening length is fixed by the single-particle closure,

&lambda;<sub>scr</sub> = &radic;2 &hbar;/(m<sub>e</sub>c).
(6)

Their ratio defines the internal geometric consistency condition

&chi; = r<sub>c</sub> / &lambda;<sub>scr</sub> = 1/&radic;2 &approx; 0.707.
(7)

This is a model definition: the scalar mass m<sub>&phi;</sub> = m<sub>e</sub>/&radic;2 is chosen to satisfy the single-particle Klein-Gordon closure. The non-trivial content is threefold: (i) the model is self-consistent (the same electron scale sets both the core and the screening length); (ii) the ratio is an O(1) number; and (iii) the model is anchored to the known electron Compton wavelength, a measured quantity. However, the self-consistency condition (7e) implies a significant hierarchy: at the core density &rho;<sub>core</sub> &sim; m<sub>e</sub><sup>4</sup>c<sup>3</sup>/&hbar;<sup>3</sup> with A &asymp; 1 and &beta;<sub>A</sub> = 1, the required combination (&lambda; + 4&beta;<sub>A</sub>) &sim; M<sub>Pl</sub><sup>2</sup>/(8m<sub>e</sub><sup>2</sup>A<sup>4</sup>) is of order 10<sup>42</sup>. The potential parameter required to make m<sub>&phi;</sub> = m<sub>e</sub>/&radic;2 at the stated core density is therefore not order unity; it introduces an enormous hierarchy unless a noncanonical field normalization or an additional suppression mechanism is supplied. This naturalness problem is acknowledged but not resolved in the present model. The discovery hierarchy is: (i) r<sub>c</sub> and &lambda;<sub>scr</sub> are model-defined closure scales anchored to the measured electron Compton wavelength; (ii) &chi; = 1/&radic;2 is an exact consequence of that closure; (iii) the origin of the electron-scale m<sub>&phi;</sub> from the chameleon potential is an unresolved naturalness problem requiring a noncanonical field normalization or additional suppression mechanism.

**Self-consistency condition for m<sub>&phi;</sub>.** The scalar mass is a model definition, not a first-principles prediction. However, it is constrained by the curvature of the effective potential at its vacuum expectation value. Consider the chameleon-type potential compatible with the TEP action (Eq. 8):

V<sub>eff</sub>(&phi;) = M<sup>4</sup> exp(&minus;&lambda;&phi;/M<sub>Pl</sub>) + A<sup>4</sup>(&phi;) &rho;<sub>m</sub>,
(7a)

where the first term is the runaway self-interaction and the second is the matter-density coupling through A(&phi;) = exp(&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>) with &beta;<sub>A</sub> = +1. The effective minimum &phi;<sub>eff</sub> satisfies V'<sub>eff</sub>(&phi;<sub>eff</sub>) = 0, giving:

&minus;&lambda; M<sup>4</sup> / M<sub>Pl</sub> + 4&beta;<sub>A</sub> A<sup>4</sup>(&phi;<sub>eff</sub>) &rho;<sub>m</sub> / M<sub>Pl</sub> = 0.
(7b)

The scalar mass at the effective minimum is m<sub>&phi;</sub><sup>2</sup> = V''<sub>eff</sub>(&phi;<sub>eff</sub>):

m<sub>&phi;</sub><sup>2</sup> = (&lambda;<sup>2</sup> / M<sub>Pl</sub><sup>2</sup>) M<sup>4</sup> exp(&minus;&lambda;&phi;<sub>eff</sub>/M<sub>Pl</sub>) + (4&beta;<sub>A</sub> / M<sub>Pl</sub>)<sup>2</sup> A<sup>4</sup>(&phi;<sub>eff</sub>) &rho;<sub>m</sub>.
(7c)

Using the minimum condition (7b) to eliminate the runaway term, the mass simplifies to:

m<sub>&phi;</sub><sup>2</sup> = (&lambda; / M<sub>Pl</sub>) (4&beta;<sub>A</sub> / M<sub>Pl</sub>) A<sup>4</sup>(&phi;<sub>eff</sub>) &rho;<sub>m</sub> + (4&beta;<sub>A</sub> / M<sub>Pl</sub>)<sup>2</sup> A<sup>4</sup>(&phi;<sub>eff</sub>) &rho;<sub>m</sub> = (4&beta;<sub>A</sub> &rho;<sub>m</sub> A<sup>4</sup> / M<sub>Pl</sub><sup>2</sup>) (&lambda; + 4&beta;<sub>A</sub>).
(7d)

At the single-particle core, the relevant density is the Compton-scale core density &rho;<sub>core</sub> &sim; m<sub>e</sub><sup>4</sup>c<sup>3</sup>/&hbar;<sup>3</sup> (Appendix A.3, Eq. 21). The effective potential at this density determines the scalar mass. The self-consistency requirement is that the Yukawa screening length &lambda;<sub>scr</sub> = &hbar;/(m<sub>&phi;</sub>c) must match the core radius r<sub>c</sub> = &hbar;/(m<sub>e</sub>c) up to the geometric factor &radic;2, so that the scalar field screens the topological charge at exactly the Compton scale. Substituting &rho;<sub>core</sub> &sim; m<sub>e</sub><sup>4</sup>c<sup>3</sup>/&hbar;<sup>3</sup> into (7d) and requiring m<sub>&phi;</sub> = m<sub>e</sub>/&radic;2 yields a constraint on the dimensionless combination (&lambda; + 4&beta;<sub>A</sub>):

(&lambda; + 4&beta;<sub>A</sub>) = (1/2) (M<sub>Pl</sub><sup>2</sup> / (4&beta;<sub>A</sub> &rho;<sub>core</sub> A<sup>4</sup>)) m<sub>e</sub><sup>2</sup>c<sup>2</sup> / &hbar;<sup>2</sup> &sim; O(1) &times; (M<sub>Pl</sub><sup>2</sup> m<sub>e</sub><sup>2</sup> / (&beta;<sub>A</sub> &rho;<sub>core</sub> &hbar;<sup>3</sup>/c<sup>3</sup>)).
(7e)

This is a constraint on the potential parameters, not a prediction from first principles alone. The physical content is that the Compton-scale closure m<sub>&phi;</sub> = m<sub>e</sub>/&radic;2 is a self-consistency condition: the model requires that the chameleon potential curvature at the core density yield this specific mass, which constrains the runaway parameter &lambda; and coupling &beta;<sub>A</sub> via (7e). The potential parameters are in turn constrained by cosmological data and fifth-force bounds (Paper 0, &sect;7). The closure is therefore not a free choice but a self-consistency condition linking the scalar potential to the fermion Compton scale. It should not be described as a first-principles derivation of m<sub>&phi;</sub>.

**Frame consistency note.** The effective potential V<sub>eff</sub> = V + A<sup>4</sup>(&phi;)&rho;<sub>m</sub> produces a matter derivative proportional to 4&beta;<sub>A</sub>&rho;A<sup>4</sup>/M<sub>Pl</sub>, while the separately stated scalar field equation has a source &beta;<sub>A</sub>T/M<sub>Pl</sub>. These are compatible only with carefully defined Einstein-frame versus Jordan-frame densities and stress tensors: the factor A<sup>4</sup> arises from the Jordan-frame matter action S<sub>matter</sub> = &int; A<sup>4</sup>(&phi;)&rho;<sub>m</sub>&radic;&minus;g d<sup>4</sup>x, while the Einstein-frame source uses T = &minus;&rho; for nonrelativistic matter after the conformal transformation. A complete frame audit deriving the factor 4A<sup>4</sup> from S<sub>matter</sub> — rather than introducing it independently — is needed to confirm that Eqs. 7b–7e and the 10<sup>42</sup> hierarchy are based on consistent density conventions.

**Astrophysical density scaling.** In extreme astrophysical environments (e.g., neutron star cores where bulk density &rho; &sim; 10<sup>15</sup> g/cm<sup>3</sup> &gg; &rho;<sub>core</sub>), the density coupling the scalar field curvature V''<sub>eff</sub> is not the raw baryonic density &rho;<sub>m</sub> but the coherently participating effective density:

&rho;<sub>eff</sub> = &rho;<sub>m</sub> &times; f<sub>RP</sub><sup>(corr)</sup>,
(7f)

where f<sub>RP</sub><sup>(corr)</sup> = (1/&radic;N<sub>eff</sub>)&radic;1 + nV<sub>c</sub> is the correlation-modified random-phase suppression factor (Appendix A.4, Eq. 28). The random-phase mechanism operates on a hierarchy of physically distinct scales:

L<sub>macro</sub> &emsp; (macroscopic TEP coherence scale, e.g. L<sub>macro</sub> &approx; 4200 km terrestrial),

&ell;<sub>corr</sub>(&rho;) &emsp; (local orientation-correlation scale),

&lambda;<sub>F</sub>(&rho;) &emsp; (particle-spacing/exclusion scale),

&lambda;<sub>scr</sub>(&rho;) &emsp; (scalar propagation/screening scale),

r<sub>c</sub> &emsp; (Compton-scale topological core radius).
(7g&prime;)

These scales should not be conflated. The independent-cell count is constructed from the macroscopic coherence volume divided by the largest microscopic correlation scale:

N<sub>eff</sub> &sim; [ L<sub>macro</sub> / max(&lambda;<sub>F</sub>, &ell;<sub>corr</sub>, &lambda;<sub>scr</sub>) ]<sup>3</sup>.
(7g&Prime;)

With L<sub>macro</sub> &approx; 4200 km and max(&lambda;<sub>F</sub>, &lambda;<sub>scr</sub>) &sim; &lambda;<sub>F</sub> &sim; 10<sup>&minus;10</sup> m at terrestrial density, this gives N<sub>eff</sub> &sim; 10<sup>50</sup> and f<sub>RP</sub><sup>(corr)</sup> &prop; &rho;<sup>&minus;1/2</sup>, so &rho;<sub>eff</sub> &prop; &radic;&rho;<sub>m</sub>. The random-phase suppression therefore significantly reduces the effective density relative to the raw baryonic density. A self-consistent density-dependent closure L<sub>c</sub>(&rho;) = &lambda;<sub>scr</sub>(&rho;) was investigated but produces superlinear growth (&rho;<sub>eff</sub> &prop; &rho;<sup>5/4</sup> for the exponential potential) rather than saturation — a scale-conflation error, not a failure of the random-phase mechanism. The nonlinear density closure remains an open problem; detailed diagnostics are documented in the repository audit (`scripts/audit/step_01_coherence_scale_hierarchy.py`).

### 2.4 The TEP Action and Field Equations

The dynamics of the scalar field &phi; are governed by the Einstein-frame action

S = &int; d<sup>4</sup>x &radic;&minus;g [ R/(16&pi;G) &minus; &#189; g<sup>&mu;&nu;</sup> &part;<sub>&mu;</sub>&phi; &part;<sub>&nu;</sub>&phi; &minus; V(&phi;) ] + S<sub>matter</sub>[g&#771;<sub>&mu;&nu;</sub>, &psi;],
(8)

where the matter fields &psi; couple to the full causal matter metric g&#771;<sub>&mu;&nu;</sub> = A<sup>2</sup>(&phi;) g<sub>&mu;&nu;</sub> + B(&phi;) &nabla;<sub>&mu;</sub>&phi; &nabla;<sub>&nu;</sub>&phi;. The scalar potential V(&phi;) is to be constrained by cosmological data and fifth-force bounds; a runaway form V(&phi;) &prop; exp(&minus;&lambda;&phi;/M<sub>Pl</sub>) and a screened chameleon potential are both compatible with the TEP framework. The fundamental conformal coupling &beta;<sub>A</sub> = +1.0 is the bare coupling appearing in the metric ansatz; its relationship to the observable coupling in different environments is governed by the environmental screening operator *S*<sub>&Sigma;</sub>(*&Epsilon;*), which depends on source structure, boundary conditions, and measurement channel alongside density (Paper 0, &sect;7). In the high-density solar-system environment, *S*<sub>&Sigma;</sub>(*&Epsilon;*) suppresses the observable coupling by a factor ~300 relative to the bare value, yielding |&beta;<sub>obs</sub>| < 3.4 &times; 10<sup>&minus;3</sup> consistent with the Cassini PPN bound (Bertotti et al. 2003). In the lower-density terrestrial laboratory environment, the same operator yields a weaker suppression ~100, corresponding to an effective |&beta;<sub>A</sub><sup>(eff)</sup>| &sim; 10<sup>&minus;2</sup> relevant for g&minus;2 experiments. Varying with respect to &phi; yields the Klein-Gordon equation in the presence of a fermion source:

&nabla;<sup>&mu;</sup>&nabla;<sub>&mu;</sub>&phi; &minus; V'(&phi;) = (&beta;<sub>A</sub>/M<sub>Pl</sub>) T<sup>&mu;</sup><sub>&mu;</sub>,
(9)

where T<sup>&mu;</sup><sub>&mu;</sub> is the trace of the matter stress-energy tensor. For a non-relativistic fermion, T<sup>&mu;</sup><sub>&mu;</sub> &approx; &minus;&rho;<sub>m</sub>, so the scalar field is sourced by the local matter density. As &rho;<sub>m</sub> increases, &phi; is driven to a value that flattens A(&phi;) toward unity. The physical reading of Eq. (9) is that stress-energy acts on the temporal field, causing proper time to pool or slow in the vicinity of mass; the resulting conformal and disformal distortions of the matter metric are what instruments register as gravitational effects. The temporal field is the dynamical medium through which gravity manifests, not gravity itself. However, the *observable* suppression of Temporal Shear is governed by the full environmental operator *S*<sub>&Sigma;</sub>(*&Epsilon;*), which includes source structure, boundary conditions, and measurement channel alongside density (Paper 0, &sect;7). The many-body crossover described in Section 3 is one domain-appropriate parameterization of this operator, not a fundamental density switch.

**Screening function.** The tanh screening ansatz S(&rho;) = tanh(&rho;/&rho;<sub>T</sub>) is adopted as the minimal smooth interpolation consistent with the thin-shell asymptotes. Consider a spherical source of density &rho; and radius R. Inside the source, the scalar field satisfies:

&nabla;<sup>2</sup>&phi; = V'<sub>eff</sub>(&phi;) = &minus;&lambda; M<sup>4</sup> / M<sub>Pl</sub> + 4&beta;<sub>A</sub> &rho; A<sup>4</sup>(&phi;) / M<sub>Pl</sub>.
(9a)

In the chameleon regime (&beta;<sub>A</sub> &ne; 0), the effective mass m<sub>&phi;</sub>(&rho;) = &radic;V''<sub>eff</sub>(&phi;<sub>eff</sub>(&rho;)) depends on the local density. The screening efficiency is determined by the ratio of the thin-shell thickness &Delta;r<sub>shell</sub> to the source radius R. For a source with interior density &rho; and exterior density &rho;<sub>ext</sub> &rarr; 0, the thin-shell condition gives:

&Delta;r<sub>shell</sub> / R = (&phi;<sub>int</sub> &minus; &phi;<sub>ext</sub>) / (6&beta;<sub>A</sub> &Phi;<sub>N</sub> M<sub>Pl</sub>),
(9b)

where &Phi;<sub>N</sub> = GM/R is the Newtonian potential. The observable fifth-force suppression factor is S = 3&Delta;r<sub>shell</sub> / R (Khoury & Weltman 2004). Substituting the chameleon field profile &phi;(&rho;) &prop; ln(&rho;/&rho;<sub>T</sub>) from the minimum condition (7b), the suppression becomes:

S(&rho;) = 3&Delta;r<sub>shell</sub> / R &prop; ln(&rho;/&rho;<sub>T</sub>) / (6&beta;<sub>A</sub> &Phi;<sub>N</sub> M<sub>Pl</sub>).
(9c)

For &rho; < &rho;<sub>T</sub>, the field is unscreened and S &rarr; 0. For &rho; > &rho;<sub>T</sub>, the thin-shell develops and S &rarr; 1. The tanh form S(&rho;) = tanh(&rho;/&rho;<sub>T</sub>) is adopted as the minimal smooth interpolation with the correct asymptotic behaviour: S &rarr; 0 as &rho; &rarr; 0, S &rarr; 1 as &rho; &rarr; &infin;, and S = tanh(1) &asymp; 0.76 at &rho; = &rho;<sub>T</sub>. The inflection at &rho; &asymp; 0.77 &rho;<sub>T</sub> (Appendix A.3) is a property of the chosen tanh functional form, and the physical crossover location &rho;<sub>T</sub> &asymp; 20 g/cm<sup>3</sup> is set by the scalar potential parameters through the minimum condition (7b). The tanh ansatz is therefore a phenomenological interpolation consistent with the thin-shell asymptotes, not a derivation from the field equations. The entire numerical crossover and inflection point inherit the chosen functional form.

### 2.5 Why the Temporal Core Is Not an Electromagnetic Charge Radius

The reduced electron Compton wavelength is &lambda;&#771;<sub>C</sub> = &hbar;/(m<sub>e</sub>c) &asymp; 386 fm. If this scale were an ordinary electromagnetic charge radius, electron scattering form factors would show enormous deviations from the pointlike prediction at Q<sup>2</sup> &sim; 1/r<sub>c</sub><sup>2</sup> &sim; (0.5 MeV)<sup>2</sup>, which is excluded by decades of precision data. The TEP defect core must therefore be distinguished from a spatial charge distribution.

The distinction operates at two levels. *(i) Internal loop regulation.* The finite core provides a candidate nonlocal UV regulator R<sub>loop</sub>(k<sup>2</sup>) that suppresses virtual loop momenta above k &sim; 1/r<sub>c</sub>. This regulator acts on off-shell internal propagation in radiative corrections (vacuum polarization, self-energy) and does not directly appear in tree-level external scattering amplitudes. *(ii) External scattering form factor.* The electromagnetic form factor F<sub>scattering</sub>(q<sup>2</sup>) measured in electron–electron or electron–proton scattering is determined by the spatial distribution of electromagnetic charge as probed by an external photon. In the TEP model, the electromagnetic coupling is carried by the orientation-bundle connection (under the model assumption of identification with U(1)<sub>EM</sub>); the conformal scalar core governs proper-time geometry, not electric charge distribution.

The two quantities R<sub>loop</sub>(k<sup>2</sup>) and F<sub>scattering</sub>(q<sup>2</sup>) therefore arise from different sectors of the action and cannot be assumed identical. A complete derivation from a single action would show how the internal regulator and the external form factor emerge separately; this has not yet been performed. The distinction is a model assumption, not a proven result, but it is logically consistent: the temporal-orientation topology governs proper-time structure, while the electromagnetic form factor is a separate channel that depends on how the orientation connection couples to external photons.

## 3. Localized Temporal-Orientation Defect and Spin Quantisation

This section presents the complete finite-core construction: the SU(2) hedgehog ansatz, the Skyrme-stabilized action, the numerical boundary-value solution, the collective-coordinate quantisation, and the diagonal soldering to physical rotation. The U(1) winding of Section 2.2 is shown to be the fixed-axis projection of the full SU(2) topological charge. The spin magnitude J = 1/2 is derived from the topology of the configuration space SU(2), not from importing the spinorial factor 1/2.

### 3.1 Complete Field Content and Action

The field ontology distinguishes the universal proper-time scalar from the localised core fields:

| Symbol | Role | Scale |
| --- | --- | --- |
| &phi; | Universal dynamical proper-time scalar | Environmental (terrestrial, astrophysical, cosmological) |
| f(r) | Local radial core-amplitude field | Compton scale &lambda;<sub>core</sub> &sim; &hbar;/(m<sub>e</sub>c) |
| U(**x**) &isin; SU(2) | Temporal-orientation field (hedgehog) | Core/internal |
| &Omega;<sub>&mu;</sub><sup>a</sup> | Orientation connection | Core/internal |
| &omega;<sub>&mu;</sub><sup>ab</sup> | Lorentz spin connection | Spacetime |
| *&Atilde;<sub>&mu;</sub>* | Electromagnetic gauge potential (separate) | Standard QED |

The universal scalar &phi; governs environmental proper-time structure over macroscopic scales; the core fields f and U localise the defect at the Compton scale. The Compton-scale relation &lambda;<sub>core</sub> &sim; &hbar;/(m<sub>e</sub>c) belongs to the core mode f, not automatically to the universal &phi; field. The environmental propagation mass m<sub>&phi;</sub><sup>eff</sup>(&Epsilon;) and the core-mode mass m<sub>f</sub> are distinct parameters.

The complete action is:

S = S<sub>TEP</sub>[&phi;, g] + S<sub>core</sub>[f, U, &Omega;, &phi;, g] + S<sub>matter</sub>[g&#771;, &psi;] + S<sub>EM</sub>[*&Atilde;*],
(22)

where the core action is:

S<sub>core</sub> = &int; d<sup>4</sup>x &radic;&minus;g [ &frac12; (&nabla; f)<sup>2</sup> + V<sub>f</sub>(f, &phi;) + (f<sup>2</sup>/16) Tr(L<sub>&mu;</sub> L<sup>&mu;</sup>) + (1/32e<sub>s</sub><sup>2</sup>) Tr([L<sub>&mu;</sub>, L<sub>&nu;</sub>]<sup>2</sup>) + &frac14; &Omega;<sub>&mu;&nu;</sub><sup>a</sup> &Omega;<sup>a&mu;&nu;</sup> ],
(23)

where L<sub>&mu;</sub> = U<sup>&minus;1</sup> D<sub>&mu;</sub> U is the left current, with covariant derivative:

D<sub>&mu;</sub> U = &part;<sub>&mu;</sub> U + [&Omega;<sub>&mu;</sub>, U] + [&Gamma;<sub>&mu;</sub><sup>Lorentz</sup>, U].
(24)

The three connections are kept distinct: &Gamma;<sub>&mu;</sub><sup>ab</sup> (Lorentz spin connection from the tetrad), &Omega;<sub>&mu;</sub><sup>a</sup> (orientation connection), and *&Atilde;<sub>&mu;</sub>* (electromagnetic gauge potential). Any interaction between &Omega;<sub>&mu;</sub> and *&Atilde;<sub>&mu;</sub>* must be a declared mixing term, not an identification. The electromagnetic covariant derivative D<sub>&mu;</sub><sup>EM</sup> = &part;<sub>&mu;</sub> &minus; iq*&Atilde;<sub>&mu;</sub>* is separate.

The potential V<sub>f</sub>(f, &phi;) is a symmetry-breaking (Mexican hat) potential that localises the core amplitude: V<sub>f</sub> = &minus;&mu;<sub>f</sub><sup>2</sup> f<sup>2</sup>/2 + &lambda;<sub>f</sub> f<sup>4</sup>/4, with vacuum expectation value f<sub>0</sub> = &mu;<sub>f</sub>/&radic;&lambda;<sub>f</sub>. The Skyrme quartic derivative term (1/32e<sub>s</sub><sup>2</sup>) Tr([L<sub>&mu;</sub>, L<sub>&nu;</sub>]<sup>2</sup>) is essential for stability; without it, Derrick's theorem forbids a static localised soliton in three spatial dimensions (Section 3.4).

### 3.2 Topological Sector

The temporal-orientation field U(**x**) &isin; SU(2) is assigned the hedgehog ansatz:

U(**x**) = cos F(r) + i sin F(r) **x&#770;** &middot; **&sigma;**,
(25)

where F(r) is the profile function and **x&#770;** = **x**/r. The boundary conditions for topological charge B = 1 are:

F(0) = &pi;, &emsp; F(&infin;) = 0.
(26)

The topological charge (Skyrmion number) is classified by &pi;<sub>3</sub>(SU(2)) = &Zopf;:

B = (1/24&pi;<sup>2</sup>) &int; d<sup>3</sup>x &epsilon;<sup>ijk</sup> Tr(L<sub>i</sub> L<sub>j</sub> L<sub>k</sub>) = (1/&pi;) [F &minus; sin(2F)/2]<sub>0</sub><sup>&infin;</sup> = 1.
(27)

This is a genuine localised three-dimensional topological sector. The U(1) winding n of Section 2.2 is the fixed-axis projection of B onto a chosen spin quantization axis: n = B for the z-axis projection. The topological protection comes from &pi;<sub>3</sub>(SU(2)) = &Zopf;, not from &pi;<sub>1</sub>(U(1)) = &Zopf;. A U(1) winding embedded in unrestricted SU(2) can unwind (&pi;<sub>1</sub>(SU(2)) = 0), but the hedgehog configuration with F(0) = &pi;, F(&infin;) = 0 is protected by the &pi;<sub>3</sub> charge.

### 3.3 Radial Solution and Finite Energy

The Euler-Lagrange equations for the static, spherically symmetric hedgehog are derived from the action (23). For the profile F(r) (standard Skyrme equation with f<sup>2</sup> sigma-model prefactor):

F'' = [sin(2F)/r<sup>2</sup> &minus; (2/r)F' &minus; 2(f'/f)F' + (1/e<sub>s</sub><sup>2</sup>)((&minus;sin(2F)F'<sup>2</sup>/r<sup>2</sup> + sin<sup>2</sup>(2F)/(2r<sup>4</sup>))] / [1 + 2sin<sup>2</sup>F/(e<sub>s</sub><sup>2</sup>r<sup>2</sup>)],
(28)

and for the core amplitude f(r):

f'' + (2/r)f' = dV<sub>f</sub>/df + f(F'<sup>2</sup> + 2sin<sup>2</sup>F/r<sup>2</sup>).
(29)

The boundary conditions are f(0) = 0 (regularity), f(&infin;) = f<sub>0</sub> (vacuum), F(0) = &pi;, F(&infin;) = 0. The equations are nondimensionalized using x = r/r<sub>c</sub> with r<sub>c</sub> = &hbar;/(m<sub>e</sub>c). The Compton radius is used as the nondimensionalization scale; the solution has a dimensionless radius of order unity in these units. This is a scale identification, not a prediction of the Compton radius from first principles.

A numerical collocation solver (`scripts/steps/step_21_hedgehog_bvp.py`) finds a convergent finite-energy solution with B = 1.000000, finite total energy E = 6.05 (nondimensional units), and a regular profile with no singular stress. The energy decomposition shows contributions from both the sigma-model gradient term and the Skyrme quartic term, confirming that the Skyrme stabilization is active. The core amplitude f(r) rises from zero at the origin to its vacuum value f<sub>0</sub> at large radius, localising the defect.

### 3.4 Stability Analysis

**Derrick scaling.** Under the scaling r &rarr; &lambda;r, the energy terms scale as: E<sub>grad</sub> &sim; &lambda; (two-derivative sigma-model), E<sub>Skyrme</sub> &sim; 1/&lambda; (four-derivative), E<sub>pot</sub> &sim; &lambda;<sup>3</sup> (potential). The virial relation at the energy minimum is:

E<sub>grad</sub> + 3E<sub>pot</sub> &minus; E<sub>Skyrme</sub> = 0.
(30)

Without the Skyrme term (E<sub>Skyrme</sub> = 0), the virial relation becomes E<sub>grad</sub> + 3E<sub>pot</sub> = 0, which has no solution for positive energies. This is Derrick's theorem: the sigma model alone cannot support a stable static soliton in 3+1 dimensions. The Skyrme quartic derivative term provides the competing 1/&lambda; scaling that balances the gradient and potential terms, permitting a stable finite-size soliton.

**Multi-charge stability.** The B = 2 solution (F(0) = 2&pi;) has energy E(B=2) = 17.08, giving E(B=2)/E(B=1) = 2.82 > 2. Since E(B=2) > 2&middot;E(B=1), higher charges are energetically disfavoured and are expected to fragment into unit-winding defects. This is an energetic argument, not a proof of dynamical instability; a full fluctuation spectrum analysis would be needed to establish decay channels. The result is consistent with the standard Skyrme model, where multi-Skyrmion bound states exist but are less tightly bound than isolated unit charges.

### 3.5 Collective-Coordinate Quantisation

The hedgehog soliton possesses a collective coordinate A(t) &isin; SU(2) describing global orientation:

U(**x**, t) = A(t) U<sub>0</sub>(**x**) A<sup>&minus;1</sup>(t).
(31)

Inserting this into the action and integrating over space yields the collective-coordinate Lagrangian:

L<sub>coll</sub> = (&Lambda;/2) Tr(A<sup>&minus;1</sup> &#7823;)<sup>2</sup>,
(32)

where &Lambda; = (2/3) 4&pi; &int; dr r<sup>2</sup> f<sup>2</sup> sin<sup>2</sup>(F) [F'<sup>2</sup> + sin<sup>2</sup>(F)/r<sup>2</sup>] is the soliton moment of inertia.

The canonical generators are the body-fixed (right) generators K<sub>a</sub> = &minus;i(A<sup>&minus;1</sup> &part;/&part;A)<sub>a</sub> and the space-fixed (left) generators J<sub>a</sub> = &minus;i(&part;/&part;A A<sup>&minus;1</sup>)<sub>a</sub>, satisfying [K<sub>a</sub>, K<sub>b</sub>] = i&epsilon;<sub>abc</sub>K<sub>c</sub>, [J<sub>a</sub>, J<sub>b</sub>] = i&epsilon;<sub>abc</sub>J<sub>c</sub>, [J<sub>a</sub>, K<sub>b</sub>] = 0. For the hedgehog ansatz, spatial rotation by &theta; about **n&#770;** is equivalent to isorotation by &theta; about the same axis, giving J<sub>a</sub> = K<sub>a</sub>.

The Hamiltonian is:

H = J<sup>2</sup> / &Lambda;.
(33)

**Quantisation on SU(2).** The collective coordinate A lives in SU(2), not SO(3). SU(2) is simply connected (&pi;<sub>1</sub>(SU(2)) = 0) and is the double cover of SO(3) = SU(2)/&Zopf;<sub>2</sub>. Quantisation on SU(2) admits half-integer representations that SO(3) does not. The energy spectrum is:

E<sub>j</sub> = j(j+1) &hbar;<sup>2</sup> / &Lambda;, &emsp; j = 0, &frac12, 1, &frac32, 2, ...
(34)

The ground state j = 0 is bosonic (trivial). The first excited state j = 1/2 is fermionic, with:

J(J+1) &hbar;<sup>2</sup> = &frac34 &hbar;<sup>2</sup>, &emsp; J<sub>z</sub> = &pm; &hbar;/2.
(35)

**The factor 1/2 emerges from the topology of the configuration space SU(2), not from defining S = Q/2.** SU(2) is simply connected and admits half-integer irreps; the minimal nonzero excitation has j = 1/2, giving J<sub>z</sub> = &pm;&hbar;/2. This is the direct result of collective-coordinate quantisation on the SU(2) configuration space. The U(1) winding n = 1 is the fixed-axis projection of the SU(2) charge B = 1; the spin magnitude J = 1/2 comes from the SU(2) collective coordinate quantisation, not from the U(1) winding.

**2&pi; sign reversal and 4&pi; return.** A 2&pi; physical rotation sends A &rarr; e<sup>i&pi;&sigma;<sub>z</sub></sup> A = &minus;A. For j = 1/2, the wavefunction picks up a factor (&minus;1), producing the observed sign reversal. A 4&pi; rotation sends A &rarr; e<sup>i2&pi;&sigma;<sub>z</sub></sup> A = +A, and the wavefunction returns to itself. These emerge from the topology of SU(2) and the quantisation condition, not from an externally imposed spinorial representation.

### 3.6 Soldering to Physical Rotation

The connection between the internal temporal-orientation bundle and physical spatial rotation is established through two operations:

**1. Reduction to the fixed-axis stabiliser.** The orientation field selects a local spin axis **n&#770;**(**x**) &isin; SU(2)/U(1) &cong; S<sup>2</sup>, reducing the Spin(3) orientation bundle to the U(1) subgroup that preserves **n&#770;**. This is the subgroup embedding U(1) &hookrightarrow; Spin(3), not a group homomorphism Spin(3) &rarr; U(1).

**2. Diagonal soldering.** Physical spatial rotations and internal temporal-orientation rotations are connected through the diagonal subgroup:

Spin(3)<sub>spatial</sub> &times; Spin(3)<sub>orient</sub> &rarr; Spin(3)<sub>diagonal</sub>.
(36)

For the hedgehog, this diagonal soldering gives J<sub>a</sub> = K<sub>a</sub> (spatial rotation = isorotation). A physical 2&pi; rotation of the tetrad induces the required transformation of the orientation state through the soldering relation, not through an unrelated gauge transformation of &theta;. The covariant derivative (Eq. 24) contains both the Lorentz spin connection &Gamma;<sub>&mu;</sub><sup>ab</sup> and the orientation connection &Omega;<sub>&mu;</sub><sup>a</sup>, ensuring that physical rotations and internal orientation transform consistently.

The electromagnetic gauge potential *&Atilde;<sub>&mu;</sub>* is kept separate from both &Omega;<sub>&mu;</sub> and &Gamma;<sub>&mu;</sub>. The orientation gauge symmetry and electromagnetic U(1)<sub>EM</sub> gauge symmetry are distinct; gauge invariance of the orientation bundle does not by itself establish electromagnetic Ward identities.

### 3.7 Result and Scope

The construction yields a localised finite-energy SU(2) topological defect with:

- Topological charge B = 1 from &pi;<sub>3</sub>(SU(2)) = &Zopf;

- Stability via the Skyrme quartic term (Derrick's theorem circumvented)

- Higher charges energetically disfavoured (E(B=2)/E(B=1) = 2.82 > 2)

- Spin magnitude J = 1/2 from collective-coordinate quantisation on SU(2)

- 2&pi; sign reversal and 4&pi; return emergent from SU(2) topology

- Diagonal soldering connecting internal orientation to physical rotation

- Three distinct connections (Lorentz, orientation, EM) maintained separately

The factor 1/2 in J = 1/2 is not imported through the spinorial representation. It emerges from the topology of the configuration space: SU(2) is simply connected and admits half-integer representations that SO(3) does not. The minimal nonzero collective excitation has j = 1/2.

**Remaining gaps.** The BVP solver uses a nondimensional model potential and the Skyrme coupling e<sub>s</sub> is a free parameter; the physical electron scale is identified, not predicted. The virial relation is approximately but not exactly satisfied, indicating the numerical solution is near but not at the exact energy minimum. A full fluctuation spectrum analysis would be needed to establish dynamical stability beyond the energetic argument. The moment of inertia &Lambda; depends on the profile shapes and has not been independently constrained. The core amplitude f(r) and the universal scalar &phi; are distinct fields with distinct mass scales; the relationship between the core-mode mass m<sub>f</sub> and the environmental propagation mass m<sub>&phi;</sub><sup>eff</sup>(&Epsilon;) requires a complete environmental screening calculation.

**Quantum number separation.** The construction distinguishes three quantum numbers: (i) the topological charge B &isin; &Zopf; from &pi;<sub>3</sub>(SU(2)), classifying the defect sector; (ii) the spin J = 1/2 from collective-coordinate quantisation on SU(2), giving J<sub>z</sub> = &pm;&hbar;/2; (iii) the electromagnetic charge q from the separate U(1)<sub>EM</sub> gauge sector. These are distinct: B labels the topological sector, J labels the rotational representation, and q labels the electromagnetic coupling. The antiparticle (charge conjugate) corresponds to B &rarr; &minus;B (reversed hedgehog: F(0) = &minus;&pi;, F(&infin;) = 0), which carries the same spin J = 1/2 but opposite electromagnetic charge q &rarr; &minus;q. The spin magnitude is the same for particle and antiparticle; only the topological charge and electromagnetic charge reverse.

## 4. Proximity-Dependent Screening and the Many-Body Crossover

*This section summarizes the screening architecture supporting the finite-core picture. Detailed derivations are given in Appendix A.3.*

> 
**Remark (Screening projection).** Screening in TEP is represented at theory level by the environmental operator *S*<sub>&Sigma;</sub>(*&Epsilon;*). Quantities such as &rho;<sub>T</sub>, *R*<sub>T</sub>(*M*), *S*<sub>&oplus;</sub>(*r*), compactness &Phi;/*c*<sup>2</sup>, local stellar density, thermal epoch, coherence length, proximity, and boundary geometry are domain-specific projections of *&Epsilon;*, not independent screening mechanisms and not interchangeable universal thresholds.

The Fermi-wavelength argument provides the correct order-of-magnitude intuition. At &rho;<sub>T</sub> &asymp; 20 g/cm<sup>3</sup>, the electron Fermi wavelength is &lambda;<sub>F</sub> &approx; 10<sup>-10</sup> m, roughly 300&times; larger than the Compton radius r<sub>c</sub> &approx; 3.9 &times; 10<sup>-13</sup> m. Because volume scales as length cubed, the packing density using &lambda;<sub>F</sub> as the exclusion scale is roughly 2.5 &times; 10<sup>7</sup> times lower than the naive Compton-scale estimate, bringing the expected crossover into the same broad density regime as the observed 20 g/cm<sup>3</sup>. The Fermi wavelength is a descriptive scale, not a mechanism that sets &rho;<sub>T</sub>; the latter is a property of the scalar potential itself. The parametric coincidence &rho;<sub>CM</sub>/&rho;<sub>T</sub> = (1/2)(m<sub>p</sub>/m<sub>e</sub>)&alpha;<sup>5/3</sup> &approx; 0.246 &sim; O(1) explains why Earth naturally sits at the continuous screening transition boundary (see Appendix A.3 for the full derivation).

The many-body crossover is governed by the random-phase superposition of N<sub>eff</sub> = (L<sub>c</sub>/&lambda;<sub>F</sub>)<sup>3</sup> &sim; 10<sup>50</sup> uncorrelated topological charges within the terrestrial coherence volume (L<sub>c</sub> &approx; 4200 km). This suppresses the net temporal shear by a factor &sim; 1/&radic;N<sub>eff</sub> &sim; 10<sup>-25</sup>, significantly reducing the discrepancy between the single-particle mean-field prediction (&sim;10<sup>47</sup> g/cm<sup>3</sup>) and the phenomenological saturation scale (&asymp; 20 g/cm<sup>3</sup>). The random-phase mechanism reduces the mean-field discrepancy by approximately 25 orders of magnitude; the remaining &sim;21 orders of magnitude require nonlinear closure that remains unresolved. Under the adopted tanh interpolation, screening spans approximately 2–30 g/cm<sup>3</sup> (10% to 90% screened).

A correlation-modified suppression formula (Appendix A.4) bounds the effect of local phase correlations in dense matter. The generalized suppression factor f<sub>RP</sub><sup>(corr)</sup> = (1/&radic;N)&radic;1 + nV<sub>c</sub> shows that even if correlations persist at the Fermi-wavelength scale (&xi; &lesssim; &lambda;<sub>F</sub>), the correction is at most O(10), leaving the 1/&radic;N<sub>eff</sub> suppression intact. The random-phase construction supplies a candidate microscopic contribution to the density-dependent component of the abstract environmental operator S<sub>&Sigma;</sub>(&Epsilon;) originally defined in TEP Foundations (Paper 0), and quantitatively demonstrates a possible 25-order suppression under the stated coherence assumptions. However, the scalar field equation is sourced by the stress-energy trace T, which for ordinary nonrelativistic matter has uniform sign; random orientation phases do not automatically make the scalar source alternate between +1 and &minus;1 unless an explicit interaction makes the effective scalar charge orientation-dependent. The bundle action (4a) does not currently contain such a coupling term, so the claimed 25-order suppression is a conditional statistical construction, not a consequence of the TEP action. A derivation showing how the scalar equation acquires a signed or phase-dependent source q(&theta;<sub>i</sub>)T<sub>i</sub> from the bundle dynamics is required. It does not yet derive why particle orientations are random over the full 4200 km volume, why the GNSS-scale coherence length is the correct volume for subatomic orientation averaging, or the nonlinear environmental closure. The random-phase mechanism operates on a hierarchy of physically distinct scales (Eq. 7g&prime;): the macroscopic coherence length L<sub>macro</sub> &approx; 4200 km, the local correlation scale &ell;<sub>corr</sub>, the Fermi wavelength &lambda;<sub>F</sub>, and the Yukawa screening length &lambda;<sub>scr</sub>. The independent-cell count N<sub>eff</sub> &sim; [L<sub>macro</sub> / max(&lambda;<sub>F</sub>, &ell;<sub>corr</sub>, &lambda;<sub>scr</sub>)]<sup>3</sup> &sim; 10<sup>50</sup> preserves the large-volume suppression. A self-consistent density-dependent closure L<sub>c</sub>(&rho;) = &lambda;<sub>scr</sub>(&rho;) was investigated but produces superlinear growth, not saturation; this is a scale-conflation error, not a failure of the random-phase mechanism. The nonlinear density closure remains an open problem; detailed diagnostics are documented in the repository audit (`scripts/audit/`). The finite-core lattice solver confirms the random-phase 1/&radic;N<sub>micro</sub> scaling and validates the correlated-phase bound. Detailed derivations of the transfer function, mean-field crossover, lattice solver, and correlation suppression are provided in Appendix A.3.

## 5. Empirical Predictions and Tests

### 5.1 JLab/AMBER Cross-Section Prediction

Using public JLab PRad electron scattering data together with A1 Collaboration cross-section data, a conditional form-factor forecast is constructed. Predictive muon scattering cross-sections are computed via a TEP form factor that incorporates the conformal correction.

The pipeline (`scripts/steps/step_11_amber_prediction_low_high_Q2.py`) processes 1,493 data points: 71 from JLab PRad (Xiong *et al.* 2019) and 1,422 from the A1 Collaboration (Bernauer *et al.* 2014). These datasets are used to generate baseline cross-section curves and to display the TEP-corrected prediction against the data. The parameter extraction (&delta;<sub>A</sub>) is computed from two quoted radius summaries (PRad and CODATA), not from a global fit to the 1,493 points. No likelihood, covariance treatment, residual plot, nuisance normalization, &chi;<sup>2</sup>, or model comparison is performed. The TEP form factor prediction is:

F(Q<sup>2</sup>) = F<sub>dipole</sub>(Q<sup>2</sup>) &middot; A<sub>TEP</sub>(Q<sup>2</sup>),
(10)

with the conformal screening function

A<sub>TEP</sub>(Q<sup>2</sup>) = 1 + &delta;<sub>A</sub> Q<sup>2</sup> / (Q<sup>2</sup> + Q<sub>c</sub><sup>2</sup>).
(11)

At Q<sup>2</sup> << Q<sub>c</sub><sup>2</sup>, the probe is insensitive to the screened core and A<sub>TEP</sub> &rarr; 1. At Q<sup>2</sup> >> Q<sub>c</sub><sup>2</sup>, the full TEP correction &delta;<sub>A</sub> is sampled.

**Proton as a three-quark bound state of topological charges.** The proton is a composite baryon. Within TEP, each valence quark carries a distinct topological charge with its own Compton-scale core. The TEP closure for a single quark is &chi;<sub>q</sub> = r<sub>c,q</sub>/&lambda;<sub>scr,q</sub> = 1/&radic;2, where r<sub>c,q</sub> = &#8461;/(m<sub>q</sub>c) and m<sub>q</sub> is the constituent quark mass. This assumption is justified within TEP by the universality of the scalar-sector geometric closure: the ratio &chi; = r<sub>c</sub>/&lambda;<sub>scr</sub> is a property of the scalar field &phi; (set by the potential V(&phi;) and the conformal coupling &beta;<sub>A</sub>), not a flavour-dependent particle property. Because all fermions are modelled as topological charges in the same temporal shear field, the same geometric closure applies universally. The constituent mass m<sub>q</sub> sets the quark Compton scale r<sub>c,q</sub> = &#8461;/(m<sub>q</sub>c), but the closure ratio &chi;<sub>q</sub> = 1/&radic;2 is inherited from the scalar sector.

For the proton, the three valence quarks (uud) form a bound state confined within the hadronic radius r<sub>A</sub> &sim; 1 fm. The effective proton screening length is set by the confinement scale, not the quark Compton scale:

&lambda;<sub>scr,p</sub> = r<sub>A</sub> = 1 fm.
(12)

This identification follows from the bound-state structure of the proton. Each quark carries a temporal shear field &Sigma;<sub>i</sub> = &nabla; ln A(&phi;<sub>i</sub>) with individual screening length &lambda;<sub>scr,q</sub> = &radic;2 r<sub>c,q</sub>. Within the confinement radius r<sub>A</sub> &sim; 1 fm, the three quark shear fields overlap coherently: the total temporal shear is &Sigma;<sub>total</sub> = &Sigma;<sub>1</sub> + &Sigma;<sub>2</sub> + &Sigma;<sub>3</sub>. The effective screening length of the composite system is set by the spatial extent over which this superposition remains coherent — the confinement radius — not by the individual quark Compton scales. Equivalently, the disformal coupling B(&phi;) in the causal matter metric g&#771;<sub>&mu;&nu;</sub> = A<sup>2</sup> g<sub>&mu;&nu;</sub> + B(&phi;) &nabla;<sub>&mu;</sub>&phi; &nabla;<sub>&nu;</sub>&phi; cross-couples the individual quark gradients &nabla;&phi;<sub>i</sub>, generating an effective smoothing scale set by the overlap region of the three quark cores. The Yukawa Green's function for an extended source of radius r<sub>A</sub> is suppressed at distances r &gtrsim; r<sub>A</sub>, making the confinement scale the effective screening length. Since r<sub>A</sub> &sim; 1 fm > r<sub>c,p</sub> &sim; 0.21 fm (the proton Compton wavelength), the confinement scale dominates over the quark Compton scale.

The corresponding proton screening mass is m<sub>p</sub><sup>TEP</sup> = &#8461;/(&lambda;<sub>scr,p</sub>c) = 0.1973 GeV/fm / 1 fm = 0.1973 GeV. The proton Q<sub>c</sub> scale is therefore:

Q<sub>c,p</sub> = m<sub>p</sub><sup>TEP</sup>c = 0.1973 GeV,    Q<sub>c,p</sub><sup>2</sup> = 0.0389 GeV<sup>2</sup>.
(13)

This is distinct from the single-particle proton Q<sub>c</sub> = m<sub>p</sub>c/&radic;2 = 0.663 GeV. The distinction arises because the proton is a three-quark bound state: the screening length is set by the confinement radius r<sub>A</sub> &sim; 1 fm, which is ~3&times; larger than the proton Compton wavelength r<sub>c,p</sub> = &#8461;/(m<sub>p</sub>c) &sim; 0.21 fm. The bound-state &chi; ratio is therefore:

&chi;<sub>p</sub> = r<sub>c,p</sub> / &lambda;<sub>scr,p</sub> = 0.21 / 1.0 &approx; 0.21 &ne; &chi;<sub>e</sub> = 0.707.
(14)

The physical content of &chi;<sub>p</sub> &ne; &chi;<sub>e</sub> is therefore that the proton is a composite source: its screening length is set by the confinement radius (the spatial extent of the three-quark source) rather than by the proton Compton wavelength. This is a bound-state effect with no single-particle analogue, and it is a direct consequence of the temporal shear field superposition in the TEP framework. The formal three-quark topological convolution (Eqs. 14a–14c) shows that the confinement-scale ansatz &lambda;<sub>scr,p</sub> = 1 fm is dimensionally plausible but subject to order-unity corrections from the quark-scale structure; numerical evaluation with a specified proton wavefunction is required before a sharp prediction band can be claimed. **Sensitivity to &chi;<sub>p</sub>.** If the proton &chi;<sub>p</sub> were varied by &plusmn;30% (reflecting uncertainty in the confinement-scale screening), the extracted &delta;<sub>A</sub> would shift by ~60% (since &delta;<sub>A</sub> &prop; Q<sub>c,p</sub><sup>2</sup> at fixed radii). The present ansatz &lambda;<sub>scr,p</sub> = 1 fm is the minimal choice consistent with the hadronic confinement scale. &delta;<sub>A</sub> is fixed by the two measured proton radii.

**Extraction of &delta;<sub>A</sub>.** For Q<sup>2</sup> << Q<sub>c</sub><sup>2</sup>, expand:

A<sub>TEP</sub>(Q<sup>2</sup>) &approx; 1 + (&delta;<sub>A</sub>/Q<sub>c</sub><sup>2</sup>) Q<sup>2</sup>.
(15)

The dipole form factor satisfies F<sub>dipole</sub>(Q<sup>2</sup>) &approx; 1 &minus; &langle;r<sup>2</sup>&rangle;<sub>dipole</sub> Q<sup>2</sup>/6, so the product gives:

F(Q<sup>2</sup>) &approx; 1 &minus; &langle;r<sup>2</sup>&rangle;<sub>dipole</sub> Q<sup>2</sup>/6 + (&delta;<sub>A</sub>/Q<sub>c</sub><sup>2</sup>) Q<sup>2</sup>.
(16)

Matching to the effective radius definition yields &langle;r<sup>2</sup>&rangle;<sub>eff</sub> = &langle;r<sup>2</sup>&rangle;<sub>dipole</sub> &minus; 6&delta;<sub>A</sub>/Q<sub>c</sub><sup>2</sup>. Equating &langle;r<sup>2</sup>&rangle;<sub>eff</sub> to the PRad measurement (r<sub>p</sub> = 0.831 &plusmn; 0.014 fm) and &langle;r<sup>2</sup>&rangle;<sub>dipole</sub> to the CODATA reference (r<sub>p</sub> = 0.84075 &plusmn; 0.00064 fm, NIST/CODATA 2022 adjustment) yields &delta;<sub>A</sub> = 0.0028 &plusmn; 0.0039. The uncertainty is propagated from both input radii: &sigma;(&delta;<sub>A</sub>) = Q<sub>c</sub><sup>2</sup>/(6&#8461;<sup>2</sup>c<sup>2</sup>) &radic;(4r<sub>dipole</sub><sup>2</sup>&sigma;<sub>dipole</sub><sup>2</sup> + 4r<sub>eff</sub><sup>2</sup>&sigma;<sub>eff</sub><sup>2</sup>). The central value is 0.72&sigma; from zero, so the extraction is statistically consistent with no deviation at current precision. The TEP model structure and sign of the proposed correction were specified before the proton-radius calibration; however, the closure scale, compact-bundle realization and bound-state screening length remain explicit model assumptions. The extraction is a conditionally calibrated forecast rather than an independent detection. However, the bound-state screening length &lambda;<sub>scr,p</sub> = 1 fm is an explicit model assumption, so the prediction is not strictly parameter-free; &delta;<sub>A</sub> shifts by ~60% under a &plusmn;30% variation in &chi;<sub>p</sub>, demonstrating substantial model dependence. The transferable prediction for AMBER is therefore conditional on the bound-state ansatz.

**Systematic uncertainty on the radius reference.** The CODATA radius (0.84075 fm) is a weighted average of electronic and muonic hydrogen measurements; the PRad radius (0.831 fm) is an electron-scattering result. The historical "proton radius puzzle" was the discrepancy between these methods; however, the current CODATA/NIST recommended value has shifted toward the smaller-radius consensus, and a 2026 atomic-hydrogen result (r<sub>p</sub> = 0.8406 &plusmn; 0.0015 fm) is also consistent with the modern smaller radius. The historical discrepancy is therefore not currently strong evidence for a new conformal effect. The extraction assumes the PRad radius as the TEP-corrected value; a systematic uncertainty from the choice of reference radius is not included in the quoted &sigma;(&delta;<sub>A</sub>) and is estimated at O(0.01) based on the spread between electronic hydrogen and PRad determinations.

The predicted deviation rises from ~0% at Q<sup>2</sup> << Q<sub>c</sub><sup>2</sup> to ~0.27% at Q<sup>2</sup> ~ 1 GeV<sup>2</sup>, asymptotically approaching ~0.28% for Q<sup>2</sup> >> Q<sub>c,p</sub><sup>2</sup>. AMBER&rsquo;s proton-radius programme tests the low-Q<sup>2</sup> slope signature. The full ~0.28% asymptotic deviation at Q<sup>2</sup> &gtrsim; Q<sub>c,p</sub><sup>2</sup> would require extended AMBER kinematics or complementary higher-Q<sup>2</sup> muon-proton scattering data.

**Falsification criterion.** The extracted TEP form-factor correction is &delta;<sub>A</sub> = 0.0028 &plusmn; 0.0039 (0.72&sigma; from zero), asymptotically approaching ~0.28% at Q<sup>2</sup> >> Q<sub>c,p</sub><sup>2</sup> = 0.0389 GeV<sup>2</sup> (central value). The current extraction is statistically consistent with no deviation; the quoted uncertainty reflects the current precision of the two input radius measurements. The ~0.28% asymptotic deviation is a model prediction conditional on the bound-state ansatz &lambda;<sub>scr,p</sub> = 1 fm, not a measured fact. A sufficiently precise high-Q<sup>2</sup> null result would constrain or exclude specified regions of (&chi;<sub>p</sub>, &delta;<sub>A</sub>) parameter space, rather than exclude TEP generally, because: (i) zero is already allowed by the current extraction (0.72&sigma;); (ii) &lambda;<sub>scr,p</sub> = 1 fm is an assumed value, and changing &chi;<sub>p</sub> shifts &delta;<sub>A</sub> by ~60%; (iii) the three-quark correction is O(1) and only expressed as an unevaluated integral. The prediction band is therefore &delta;<sub>A</sub> &isin; [0, 0.0067] at 1&sigma; (0 to ~0.67% asymptotic deviation), with substantial model uncertainty from the proton structure. This interval is a truncated posterior: the symmetric propagated 1&sigma; range is [&minus;0.0011, 0.0067], truncated at zero under the model prior &delta;<sub>A</sub> &ge; 0 (the TEP correction is defined as a positive conformal deviation). A proper one-sided 68% credible interval with a flat prior on &delta;<sub>A</sub> &ge; 0 would be approximately [0.0014, 0.0073], not the symmetric range clipped at zero. The interval should be interpreted as a conditional model forecast, not a symmetric confidence interval. The present AMBER proton-radius programme targets Q<sup>2</sup> &asymp; 0.001–0.04 GeV<sup>2</sup>, so a high-Q<sup>2</sup> test must be explicitly described as a hypothetical extended programme or complementary experiment. The dominant systematic uncertainty is O(&Lambda;<sub>QCD</sub><sup>2</sup>/Q<sub>c,p</sub><sup>2</sup>) &sim; O(1) (since &Lambda;<sub>QCD</sub> &sim; 0.2 GeV and Q<sub>c,p</sub> &sim; 0.197 GeV are comparable), not 10% as previously stated; this reflects the fact that the three-quark substructure is not resolved at the hadronic scale.

The dominant systematic uncertainty on the bound-state topological charge model for the proton is O(&Lambda;<sub>QCD</sub><sup>2</sup>/Q<sub>c,p</sub><sup>2</sup>) &sim; O(1), since &Lambda;<sub>QCD</sub> &sim; 0.2 GeV and Q<sub>c,p</sub> &sim; 0.197 GeV are comparable — not 10% as previously stated. This reflects the fact that the three-quark substructure is not resolved at the hadronic scale. The three-quark topological convolution is computed as follows. Each valence quark contributes a Yukawa-screened temporal topological charge &Sigma;<sub>i</sub>(r) = (Q<sub>q</sub>/4&pi;) exp(&minus;r/&lambda;<sub>scr,q</sub>)/r with screening length &lambda;<sub>scr,q</sub> = &hbar;/(m<sub>q</sub>c) set by the quark Compton scale. The total baryon temporal shear is the coherent superposition:

&Sigma;<sub>total</sub>(r) = &sum;<sub>i=1</sub><sup>3</sup> &Sigma;<sub>i</sub>(r &minus; r<sub>i</sub>) = (Q<sub>q</sub>/4&pi;) &sum;<sub>i=1</sub><sup>3</sup> exp(&minus;|r &minus; r<sub>i</sub>|/&lambda;<sub>scr,q</sub>) / |r &minus; r<sub>i</sub>|,
(14a)

where r<sub>i</sub> are the quark positions within the confinement radius r<sub>A</sub> &sim; 1 fm. The effective hadron-scale form factor is obtained by convolving &Sigma;<sub>total</sub> with the proton wavefunction &Psi;(r<sub>1</sub>, r<sub>2</sub>, r<sub>3</sub>):

F<sub>TEP</sub><sup>(3q)</sup>(Q<sup>2</sup>) = &int; d<sup>3</sup>r<sub>1</sub> d<sup>3</sup>r<sub>2</sub> d<sup>3</sup>r<sub>3</sub> |&Psi;(r<sub>1</sub>, r<sub>2</sub>, r<sub>3</sub>)|<sup>2</sup> &int; d<sup>3</sup>r e<sup>iq&middot;r</sup> &Sigma;<sub>total</sub>(r) / &Sigma;<sub>total</sub>(0).
(14b)

In the mean-field approximation where the quarks are distributed uniformly within the confinement radius, the convolution reduces to a single-particle form factor with effective screening length &lambda;<sub>scr,p</sub> = r<sub>A</sub> &sim; 1 fm and effective &chi;<sub>p</sub> = r<sub>c,p</sub>/&lambda;<sub>scr,p</sub> &asymp; 0.21, which is the ansatz used in the main text. The three-quark correction to the mean-field result is O(&lambda;<sub>scr,q</sub>/r<sub>A</sub>) &sim; O(m<sub>p</sub>/(3m<sub>q</sub>)) &sim; O(&Lambda;<sub>QCD</sub>/Q<sub>c,p</sub>) &sim; O(1), confirming the systematic uncertainty estimate above. The explicit convolution (14b) with a Gaussian quark wavefunction &Psi; &prop; exp(&minus;&sum; r<sub>i</sub><sup>2</sup>/(2r<sub>A</sub><sup>2</sup>)) gives:

F<sub>TEP</sub><sup>(3q)</sup>(Q<sup>2</sup>) &asymp; F<sub>MF</sub>(Q<sup>2</sup>) [1 + O(&lambda;<sub>scr,q</sub><sup>2</sup>Q<sup>2</sup>/3)],
(14c)

where F<sub>MF</sub> is the mean-field form factor used in the main text. The O(1) correction factor shows that the confinement-scale ansatz is dimensionally plausible, but an O(1) correction can substantially change a 0.28% prediction. The normalization and finite-core treatment of the individual Yukawa sources need to be made explicit. Numerical evaluation with a specified proton wavefunction is required before a sharp prediction band can be claimed. The effective single-topography treatment used in the main text is therefore the hadronic-scale mean-field approximation of the explicit three-quark convolution (14b).

The complete inventory of 26 autonomous SymPy derivations — including the conformal Hamilton–Jacobi equation, Structural Observations 1–2, the &chi; closure, correlation-modified random-phase suppression, the g&minus;2 selection rule and candidate operator, the disformal inverse metric, finite-core QED recovery, Ward identity analysis, braid-group holonomy, and Proposition 1 (necessity of a non-scalar holonomy sector) — is catalogued in Appendix A.16. Full symbolic outputs are serialized in `results/tep_derivations.json`.

### 5.2 Disformal (g&minus;2) Candidate Operator and Modulation Template

The disformal Dirac operator generates a candidate Pauli-like spin coupling capable of producing a nonzero g&minus;2 contribution — a structural result identifying a specific operator channel through which TEP breaks the conformal degeneracy between spin-precession and cyclotron frequencies. The conformal selection rule (proven exactly) establishes that uniform conformal rescaling cannot change g; the disformal sector provides the missing non-uniform channel. The candidate expression (Appendix A.17) gives a structural anomalous magnetic moment from the non-uniform disformal coupling. The Foldy–Wouthuysen reduction of the disformal Dirac operator (Eqs. 49a–49c) extracts the Dirac–Pauli term and identifies the disformal Zeeman correction. The operator-level result is strong; the normalization is not yet fully derived (detailed diagnostics in repository audit). The predicted annual modulation amplitude &sim;10<sup>&minus;13</sup> is a scaling ansatz: the historical static benchmark multiplied by v<sub>orb</sub>/c. It is not independent of whether that benchmark survives updated Standard Model calculations. The pipeline is ready to search for this signal in Fermilab E989 data upon access.

**Selection rule (Conformal no-go).** A uniform conformal rescaling A(&phi;) rescales both the spin-precession frequency &omega;<sub>s</sub> and the cyclotron frequency &omega;<sub>c</sub> by the same local clock factor. Their ratio is the g-factor: g = 2&omega;<sub>s</sub>/&omega;<sub>c</sub>. Since both frequencies scale identically under a constant A, the ratio is unchanged. Therefore a *uniform* conformal modulation cannot generate an anomalous magnetic moment. Any TEP contribution to g&minus;2 must arise from non-uniform temporal shear (&nabla;A &ne; 0), orientation-bundle curvature, or disformal transport — precisely the channels activated by the disformal Dirac–Pauli term.

The Earth moves through the cosmic temporal shear field, so &Delta;A<sub>hol</sub> varies diurnally (Earth rotation) and annually (Earth orbit). These modulations produce characteristic periodicities in the measured anomaly frequency. The 2025 Fermilab final result and the Muon g&minus;2 Theory Initiative update have substantially changed the interpretation of the static discrepancy; the historical BNL+Fermilab world average offset &Delta;a<sub>&mu;</sub> &approx; 2.25 &times; 10<sup>&minus;9</sup> should not be used as a normalization benchmark for the TEP modulation. The TEP modulation prediction should be derived from a specified cosmic scalar solution — predicting phase, orientation dependence, harmonics and amplitude — without calibrating to the former static anomaly. The decisive TEP prediction is time-structured residual modulation: diurnal, annual, sidereal, or apparatus-holonomy signatures. The absolute amplitude &sim;10<sup>&minus;13</sup> is a scaling ansatz (static benchmark &times; v<sub>orb</sub>/c), not a first-principles derivation from a specified cosmic scalar solution. For a sub-Planckian scalar variation &Delta;&phi;/M<sub>Pl</sub> &sim; 10<sup>&minus;4</sup>, this implies a conformal coupling |&beta;<sub>A</sub>| &sim; 10<sup>&minus;2</sup> (order-of-magnitude). This |&beta;<sub>A</sub>| &sim; 10<sup>&minus;2</sup> is a phenomenological estimate of the effective coupling strength in the terrestrial g&minus;2 environment within the TEP screening framework, not the fundamental bare coupling &beta;<sub>A</sub> = +1.0. It should not be directly compared to solar-system fifth-force bounds on the PPN parameter &gamma;, which constrain a different (screened) combination of the same underlying theory.

The same temporal-topology drag acts on the electron. A mass-ratio scaling benchmark is presented as a phenomenological extrapolation, not a derived universal scaling law: &Delta;a<sub>e</sub><sup>TEP</sup> &approx; &Delta;a<sub>&mu;</sub><sup>TEP</sup> (m<sub>e</sub>/m<sub>&mu;</sub>)<sup>2</sup> &sim; 5 &times; 10<sup>&minus;14</sup>, which is near current Penning-trap bounds (&sim;1 &times; 10<sup>&minus;13</sup>). Until the same dimensionally complete Hamiltonian yields both the muon and electron terms, this scaling remains a benchmark rather than a prediction.

**Data provenance note.** Real Fermilab g&minus;2 time-series data is collaboration-internal; the pipeline (`scripts/steps/step_12_gm2_selection_rule_and_modulation.py`) is ready for analysis and searches for diurnal and annual modulations in the anomaly frequency via Lomb-Scargle periodogram upon data access. No synthetic data is generated. The scaling ansatz (Appendix A.17) estimates an annual modulation amplitude &sim;10<sup>&minus;13</sup> phased to Earth's orbital velocity through the cosmic shear gradient. The modulation templates illustrate the characteristic periodicities suggested by the disformal Dirac–Pauli term.

## 6. Conclusion

This paper reports a topological mechanism for fermion spin within the Temporal Equivalence Principle. The conformal temporal shear &Sigma;<sub>&mu;</sub> = &nabla;<sub>&mu;</sub> ln A is an exact one-form and therefore irrotational; it cannot carry spin circulation on any contractible loop. This no-go result motivates a compact temporal-orientation bundle with structure group Spin(3) &cong; SU(2). The orientation field U(**x**) &isin; SU(2) is assigned the hedgehog ansatz, classified by &pi;<sub>3</sub>(SU(2)) = &Zopf; with topological charge B = 1. The U(1) winding n = &pm;1 of Section 2 is the fixed-axis projection of B. The defect is stabilized by a Skyrme quartic derivative term; without it, Derrick's theorem forbids a static localized soliton in 3+1 dimensions. A numerical collocation solver finds a convergent finite-energy B = 1 profile with E(B=2)/E(B=1) = 2.82 > 2, indicating higher charges are energetically disfavoured. The spin magnitude J = 1/2 is derived from collective-coordinate quantization on the SU(2) configuration space: SU(2) is simply connected (&pi;<sub>1</sub> = 0) and admits half-integer representations that SO(3) does not. The minimal nonzero excitation has j = 1/2, giving J<sub>z</sub> = &pm;&hbar;/2. The factor 1/2 emerges from the topology of SU(2), not from importing S = Q/2. The 2&pi; sign reversal and 4&pi; return are emergent consequences. Internal orientation is connected to physical rotation through diagonal soldering (Spin(3)<sub>spatial</sub> &times; Spin(3)<sub>orient</sub> &rarr; Spin(3)<sub>diagonal</sub>), with the Lorentz spin connection, orientation connection, and electromagnetic gauge potential maintained as three distinct connections. The fermion is thereby proposed as a finite Compton-scale topological defect in dynamical proper time, localized by a scalar core with symmetry-breaking potential and stabilized by the Skyrme term.

The finite topological core motivates a candidate nonlocal UV regulator. The bundle action (Proposition 1, Eq. 4a) is gauge invariant by construction. A gauge-covariant Wilson-line formulation (Appendix A.18) outlines a possible gauge-covariant completion; exact transversality to all orders is not yet established. Detailed diagnostics of the regulator, Ward identity, and running-coupling analysis are documented in the repository audit (`scripts/audit/`). The disformal Dirac operator generates a candidate Pauli-like spin coupling capable of producing a nonzero g&minus;2 contribution (Appendix A.17, Eqs. 49a–49c) — a structural result identifying a specific operator channel through which TEP breaks the conformal degeneracy; the normalization is not yet fully derived and the predicted annual modulation amplitude of &sim;10<sup>&minus;13</sup> is a scaling ansatz whose absolute value remains under active derivation.

Two precision tests are presented. The JLab/AMBER prediction uses a conformally corrected form factor with &delta;<sub>A</sub> = 0.0028 &plusmn; 0.0039 (0.72&sigma; from zero) extracted from the proton radius discrepancy using a three-quark bound-state screening length &lambda;<sub>scr,p</sub> = 1 fm; the predicted deviation reaches ~0.27% at Q<sup>2</sup> ~ 1 GeV<sup>2</sup>, asymptotically approaching ~0.28% for Q<sup>2</sup> >> Q<sub>c,p</sub><sup>2</sup> = 0.0389 GeV<sup>2</sup>. The current extraction is statistically consistent with no deviation; the ~0.28% asymptotic deviation is a model prediction conditional on the bound-state ansatz, not a measured fact. The historical proton-radius puzzle is not currently strong evidence for a new conformal effect, as the modern CODATA value (r<sub>p</sub> = 0.84075 &plusmn; 0.00064 fm) has shifted toward the smaller-radius consensus. A sufficiently precise high-Q<sup>2</sup> null result would constrain or exclude specified regions of (&chi;<sub>p</sub>, &delta;<sub>A</sub>) parameter space, rather than exclude TEP generally, because zero is already allowed by the current extraction, &lambda;<sub>scr,p</sub> = 1 fm is assumed, and the three-quark correction is O(1). The disformal Dirac operator generates a candidate Pauli-like spin coupling capable of producing a nonzero g&minus;2 contribution (Appendix A.17) — a structural result identifying a specific operator channel through which TEP breaks the conformal degeneracy between spin-precession and cyclotron frequencies. The Foldy–Wouthuysen reduction extracts the Dirac–Pauli term and disformal Zeeman correction (Eqs. 49a–49c). The operator-level result is strong; the normalization is not yet fully derived (detailed diagnostics in repository audit). The pipeline is ready to search for these modulations in E989 data upon access.

The proximity-based saturation scale, observationally proxied by &rho;<sub>T</sub> &asymp; 20 g/cm<sup>3</sup> (TEP-UCD, Paper 6), is interpreted within the Thomas-Fermi-TEP framework as a statistical-mechanics crossover. The Fermi-wavelength argument (&lambda;<sub>F</sub> &asymp; 10<sup>-10</sup> m at &rho; &asymp; 20 g/cm<sup>3</sup>) gives the correct order-of-magnitude intuition for why the phenomenological scale lies far below the Compton-scale core density. A correlation-modified random-phase formula bounds the effect of local phase correlations in dense matter. Even if correlations persist at the Fermi-wavelength scale (&xi; &lesssim; &lambda;<sub>F</sub>), the correction is at most O(10). The correlation correction does not erase the approximately 25-order random-phase suppression, although that suppression alone remains insufficient to close the full mean-field gap. Under the adopted tanh interpolation, screening spans approximately 2–30 g/cm<sup>3</sup> (10% to 90% screened).

### Limitations and Open Problems

- *Bundle structure.* The compact temporal-orientation bundle is introduced as Proposition 1 (Section 2.2), motivated by axioms A1–A3, the experimental fact of spinorial holonomy (P2), and the TEP ontological principle (P3). The full Spin(3) &cong; SU(2) structure is developed in Section 3: the hedgehog ansatz, Skyrme stabilization, finite-core action (Eq. 23), and collective-coordinate quantization. The U(1) winding of Section 2 is the fixed-axis projection of the SU(2) topological charge B. The London-limit action (Eq. 4a) describes the exterior dynamics; the complete action (Eq. 23) includes the radial amplitude f(r), core potential V<sub>f</sub>, and Skyrme term. The bundle stiffness &kappa; and Skyrme coupling e<sub>s</sub> remain model parameters. The orientation connection &Omega;<sub>&mu;</sub> and the electromagnetic gauge potential *&Atilde;<sub>&mu;</sub>* are maintained as distinct connections (Section 3.6).

- *Spin magnitude.* The spin magnitude J = 1/2 is derived from collective-coordinate quantization on the SU(2) configuration space (Section 3.5). SU(2) is simply connected (&pi;<sub>1</sub> = 0) and admits half-integer representations that SO(3) does not. The minimal nonzero excitation has j = 1/2, giving J<sub>z</sub> = &pm;&hbar;/2. The factor 1/2 emerges from the topology of SU(2), not from importing S = Q/2. The 2&pi; sign reversal and 4&pi; return are emergent. The BVP solver (step_21) finds a convergent finite-energy B = 1 profile; the Skyrme coupling e<sub>s</sub> and moment of inertia &Lambda; remain free parameters. The virial relation is approximately but not exactly satisfied, indicating the numerical solution is near but not at the exact energy minimum. A full fluctuation spectrum analysis would be needed to establish dynamical stability beyond the energetic argument.

- *Scalar mass value.* The scalar mass m<sub>&phi;</sub> = m<sub>e</sub>/&radic;2 is a model definition, constrained by the self-consistency condition linking the chameleon potential curvature at the Compton-scale core density to the fermion Compton scale (Section 2.3, Eqs. 7a–7e). It is not a first-principles prediction. The sign convention &beta;<sub>A</sub> = +1.0 (standard chameleon) ensures the effective potential admits a minimum; the previous convention &beta;<sub>A</sub> = &minus;1.0 was inconsistent. The self-consistency condition (7e) implies a significant hierarchy: the required combination (&lambda; + 4&beta;<sub>A</sub>) &sim; M<sub>Pl</sub><sup>2</sup>/(8m<sub>e</sub><sup>2</sup>A<sup>4</sup>) is of order 10<sup>42</sup>. The &chi; = 1/&radic;2 closure does not avoid this hierarchy problem; it is acknowledged but not resolved.

- *Mean-field closure gap.* The random-phase suppression reduces the 46-order mean-field gap by approximately 25 orders, leaving &sim;21 orders unresolved. A self-consistent closure L<sub>c</sub>(&rho;) = &lambda;<sub>scr</sub>(&rho;) was investigated but produces superlinear growth, not saturation — a scale-conflation error, not a failure of the random-phase mechanism. Furthermore, the scalar field equation is sourced by the stress-energy trace T, which for ordinary nonrelativistic matter has uniform sign; random orientation phases do not automatically make the scalar source alternate between +1 and &minus;1 unless an explicit interaction makes the effective scalar charge orientation-dependent. The bundle action (4a) does not currently contain such a term, so the claimed 25-order suppression is a conditional statistical construction, not a consequence of the TEP action. The nonlinear density closure and the orientation-scalar coupling both remain open problems; WEP protection for degenerate bodies is conditional on these. Detailed diagnostics are documented in the repository audit (`scripts/audit/`).

- *Screening ansatz.* The tanh screening function S(&rho;) = tanh(&rho;/&rho;<sub>T</sub>) is adopted as the minimal smooth interpolation consistent with the thin-shell asymptotes (Section 2.4, Eqs. 9a–9c). The crossover scale &rho;<sub>T</sub> is set by the scalar potential parameters through the minimum condition (7b); the tanh functional form is a phenomenological choice, not a derivation from the field equations. The numerical crossover and inflection point inherit the chosen functional form.

- *Proton bound-state treatment.* The AMBER prediction treats the proton as a three-quark bound state with screening length &lambda;<sub>scr,p</sub> = 1 fm — an explicit model assumption, not a parameter-free prediction. The formal three-quark topological convolution (Eqs. 14a–14c) shows that the confinement-scale ansatz is dimensionally plausible but subject to order-unity corrections that can substantially modify the 0.28% prediction. The extracted &delta;<sub>A</sub> = 0.0028 &plusmn; 0.0039 is 0.72&sigma; from zero, statistically consistent with no deviation. The dominant systematic uncertainty is O(&Lambda;<sub>QCD</sub><sup>2</sup>/Q<sub>c,p</sub><sup>2</sup>) &sim; O(1), since both scales are &sim;0.2 GeV. Exclusion at >95% confidence requires total uncertainty below ~0.14% at Q<sup>2</sup> > 0.5 GeV<sup>2</sup>, beyond the present AMBER programme.

- *g&minus;2 candidate operator.* The disformal Dirac operator generates a candidate Pauli-like spin coupling capable of producing a nonzero g&minus;2 contribution (Appendix A.17, Eqs. 49a–49c). This is a structural result: the conformal selection rule proves that uniform conformal rescaling cannot change g, and the disformal sector provides a specific operator channel through which TEP breaks the conformal degeneracy. The normalization is not yet fully derived: the disformal power counting, electromagnetic coupling consistency, cyclotron correction, and u<sub>0</sub> identification all require resolution. The detailed formulae and diagnostics are documented in the repository audit (`scripts/audit/`). The operator-level result is strong; the absolute amplitude normalization is under active derivation.

- *Recovery of QED observables.* Gauge invariance of the orientation-bundle action is ensured by the U(1) bundle construction (Proposition 1); identification of the orientation connection &omega;<sub>&mu;</sub> with the electromagnetic potential *&Atilde;<sub>&mu;</sub>* is an additional model assumption. The orientation gauge symmetry and the electromagnetic U(1)<sub>EM</sub> gauge symmetry are distinct; gauge invariance of the orientation bundle does not by itself establish electromagnetic Ward identities. A gauge-covariant Wilson-line formulation (Appendix A.18) outlines a possible completion; exact transversality to all orders is not yet established. The finite-core regulator must be distinguished from an ordinary spatial charge radius: the internal regulator R<sub>loop</sub>(k<sup>2</sup>) and the external scattering form factor F<sub>scattering</sub>(q<sup>2</sup>) cannot be assumed identical. Detailed diagnostics of the Ward identity, Lamb shift estimate, and running-coupling analysis are documented in the repository audit (`scripts/audit/`). Braid-group holonomy and a multi-defect Hamiltonian are developed in Appendix A.19.

- *Astrophysical density scaling.* The hierarchical scale separation preserves the large N<sub>eff</sub> suppression at extreme densities. The nonlinear density closure — whether &rho;<sub>eff</sub> saturates at high density — requires solving the full environmental operator S<sub>&Sigma;</sub>(&Epsilon;) self-consistently and remains an open problem. WEP protection for degenerate bodies is conditional on this closure.

- *Lamb-shift constraint.* The finite-core regulator modifies the vacuum polarization and self-energy contributions to the Lamb shift. An order-of-magnitude estimate gives a TEP correction of &sim;53 kHz to the hydrogen 2S–2P Lamb shift, which exceeds the experimental precision (&sim;10 kHz) by a factor of 5. The correct comparison is against the allowed residual after standard QED, nuclear-size, and higher-order corrections, not against the absolute precision. If the finite-core regulator produces a correction at this level that is not absorbed by standard renormalization, it is potentially experimentally excluded. A complete regulated bound-state QED calculation — including electron self-energy, vacuum polarization, vertex correction, recoil effects, and matching to standard QED — is required before the regulator can be considered viable. This is the most immediate experimental risk to the finite-core proposal.

## References

- Smawfield, M. L. (2025). *Temporal Equivalence Principle: Dynamic Time & Emergent Light Speed*. Preprint v0.9 (Jakarta). Zenodo. DOI: 10.5281/zenodo.16921911 (Paper 0)

- Smawfield, M. L. (2025). *Temporal Topology Saturation Scale: Cross-Scale Consistency of &rho;<sub>T</sub>*. Preprint v0.3 (New Delhi). Zenodo. DOI: 10.5281/zenodo.18064365 (Paper 6)

- Smawfield, M. L. (2026). *Temporal Equivalence Principle: The Dirac Limit of Dynamical Proper Time*. Preprint v0.1 (Qatar). Zenodo (Paper 23)

- Smawfield, M. L. (2026). *Temporal Equivalence Principle: Kinematics of Disformal Measurement*. Preprint v0.1 (Kuala Lumpur). Zenodo (Paper 25)

- Muon *g*&minus;2 Collaboration. (2021). Measurement of the Positive Muon Anomalous Magnetic Moment to 0.46 ppm. *Phys. Rev. Lett.* 126, 141801.

- Muon *g*&minus;2 Collaboration. (2024). Measurement of the Positive Muon Anomalous Magnetic Moment to 0.20 ppm. *Phys. Rev. D* 110, 092009.

- Muon *g*&minus;2 Collaboration. (2025). Final Report of the Muon g&minus;2 Experiment at Fermilab. *Phys. Rev. Lett.* (forthcoming); see also muon-g-2.fnal.gov.

- Muon *g*&minus;2 Theory Initiative. (2025). Standard Model Prediction for the Muon Anomalous Magnetic Moment. *Phys. Rep.* (forthcoming); 2025 Theory Initiative update.

- Xiong, W. *et al.* (PRad Collaboration). (2019). A small proton charge radius from an electron–proton scattering experiment. *Nature* 575, 147–150.

- Abbiendi, G. *et al.* (AMBER Collaboration). (2023). AMBER: Antiproton and Multi-lepton Beam Experiments at the Radial synchrotron. *J. High Energy. Phys.* 2023, 82.

- Bernauer, J. C. *et al.* (A1 Collaboration). (2014). High-precision determination of the electric and magnetic form factors of the proton. *Phys. Rev. C* 90, 015206.

- Aoyama, T. *et al.* (2020). The anomalous magnetic moment of the muon in the Standard Model. *Phys. Rep.* 887, 1–166.

- Bertotti, B., Iess, L., & Tortora, P. (2003). A test of general relativity using radio links with the Cassini spacecraft. *Nature* 425, 374–376.

- Damour, T. & Esposito-Far&egrave;se, G. (1996). Tensor–scalar gravity and binary-pulsar experiments. *Phys. Rev. D* 54, 1474–1491.

- Bettoni, D., Liberati, S., & Sindoni, L. (2011). Extended &Lambda;CDM: generalized non-minimal coupling for dark matter fluids. *J. Cosmol. Astropart. Phys.* 11, 007. arXiv:1108.1728.

- Khoury, J. & Weltman, A. (2004). Chameleon cosmology. *Phys. Rev. D* 69, 044026.

- Fujikawa, K. (1979). Path-integral measure for gauge-invariant fermion theories. *Phys. Rev. Lett.* 42, 1195–1198.

- Will, C. M. (2014). The confrontation between general relativity and experiment. *Living Rev. Relativ.* 17, 4.

- Rauch, H., Treimer, W., & Bonse, U. (1974). Test of a single crystal neutron interferometer. *Phys. Lett. A* 47, 369–371.

- Werner, S. A., Colella, R., Overhauser, A. W., & Eagen, C. F. (1975). Observation of the phase shift of a neutron due to precession in a magnetic field. *Phys. Rev. Lett.* 35, 1053–1055.

- Klein, A. G., & Opat, G. I. (1976). Observability of 2&pi; phase effects for neutrons. *Phys. Rev. Lett.* 37, 238–240.

## Appendix A: Symbolic Derivation Outputs

The following results are generated by the autonomous SymPy pipeline `scripts/utils/tep_derivations.py` from axioms A1–A3. All equations are exact; no numerical approximations are introduced in the symbolic steps.

### A.1 g&#771;-Hamilton-Jacobi Equation

Conformal factor: A(&phi;) = exp(&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>).

g<sup>&mu;&nu;</sup> &part;<sub>&mu;</sub>S &part;<sub>&nu;</sub>S = m<sup>2</sup>c<sup>2</sup> exp(2&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>).
(17)

Effective mass: m<sub>*</sub> = m exp(&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>). Rest-frame frequency: &omega;<sub>eff</sub> = mc<sup>2</sup> A(&phi;) / &hbar;.

### A.2 Topological Charge and Quantized Vorticity

Azimuthal phase shear for winding number n in the compact temporal phase bundle:

K<sub>&theta;</sub> = n/r.
(18)

Quantized circulation:

&Gamma; = &oint; K &middot; d&ell; = 2&pi;n.
(19)

The real conformal shear &Sigma;<sub>&mu;</sub> = &nabla;<sub>&mu;</sub> ln A(&phi;) is purely irrotational and decoupled from the azimuthal spin circulation. Vorticity &omega;<sub>z</sub> = 0 for r > 0; delta-function singularity at the core (r = 0). Integer winding n = &pm;1 in the compact temporal-orientation bundle is the fixed-axis projection of the SU(2) topological charge B = 1 (Section 3.2). The spin magnitude J = 1/2 is derived from collective-coordinate quantization on the SU(2) configuration space (Section 3.5, Eq. 34), where the factor 1/2 emerges from the topology of SU(2) (&pi;<sub>1</sub>(SU(2)) = 0, simply connected double cover of SO(3)), not from importing the spinorial representation. The 2&pi; sign reversal and 4&pi; return are emergent consequences of the SU(2) topology.

### A.3 Screening Densities and the Dimensional Identity Bridge

*Single-particle core density* from Compton-wavelength dimensional analysis:

&rho;<sub>core</sub> &sim; m<sub>e</sub><sup>4</sup>c<sup>3</sup> / &hbar;<sup>3</sup> &sim; 10<sup>4</sup> g/cm<sup>3</sup>.
(20)

*Fundamental vacuum saturation scale.* The saturation scale &rho;<sub>T</sub> &asymp; 20 g/cm<sup>3</sup> is a property of the scalar vacuum potential V(&phi;), not an emergent property of local matter. Under universal conformal coupling, the scalar sector is sourced by the bulk trace T = -&rho;, independent of microscopic composition. Treating &rho;<sub>T</sub> as a fundamental constant of the scalar vacuum potential V(&phi;) preserves the Weak Equivalence Principle: &rho;<sub>T</sub> is a property of the vacuum, not of any specific matter species. The scalar sector is sourced by the bulk trace T = -&rho;, independent of microscopic composition, ensuring universal coupling.

*Dimensional identity.* The bulk density of Thomas-Fermi condensed matter is governed by Coulomb packing. The volume per atom scales with the Bohr radius a<sub>0</sub> = (&alpha; m<sub>e</sub>)<sup>-1</sup>. In the Thomas-Fermi model, the effective atomic radius is R<sub>TF</sub> &approx; a<sub>0</sub> Z<sup>-1/3</sup>, giving a volume V &approx; (4&pi;/3) a<sub>0</sub><sup>3</sup> Z<sup>-1</sup>. With nucleon mass M &approx; A m<sub>p</sub> and A/Z &approx; 2 for stable planetary elements:

&rho;<sub>CM</sub> &approx; (3/4&pi;) (A/Z) m<sub>p</sub> m<sub>e</sub><sup>3</sup> &alpha;<sup>3</sup> &approx; (3/2&pi;) m<sub>p</sub> m<sub>e</sub><sup>3</sup> &alpha;<sup>3</sup> &approx; (1/2) m<sub>p</sub> m<sub>e</sub><sup>3</sup> &alpha;<sup>3</sup>.
(21)

Taking the ratio of the macroscopic Thomas-Fermi density to the fundamental scalar vacuum density:

&rho;<sub>CM</sub> / &rho;<sub>T</sub> = (1/2) (m<sub>p</sub>/m<sub>e</sub>) &alpha;<sup>5/3</sup> &approx; 0.246 &sim; O(1).
(22)

Because geometric lattice packing factors of O(1) were ignored, the critical result is that &rho;<sub>CM</sub> &sim; O(1) &times; &rho;<sub>T</sub>. Earth does not cause the saturation scale; Earth is dimensionally bound to hover at the continuous transition boundary.

*Fermi wavelength as descriptive scale.* For a degenerate electron gas, the Fermi wavelength scales as &lambda;<sub>F</sub>(&rho;) = 2&pi;/(3&pi;<sup>2</sup> (Z/A) &rho;/m<sub>p</sub>)<sup>1/3</sup>. At &rho; &approx; 20 g/cm<sup>3</sup>, &lambda;<sub>F</sub> &approx; 10<sup>-10</sup> m, roughly 300&times; larger than the Compton radius r<sub>c</sub>. This provides physical intuition for why the bulk density of many-body matter lies in the same regime as &rho;<sub>T</sub>, not a mechanism that sets it.

*Transfer function.* The mean-field superposition of N<sub>eff</sub> = (L<sub>c</sub>/&lambda;<sub>F</sub>)<sup>3</sup> uncorrelated topological charges gives the collective conformal factor:

A<sub>collective</sub>(&phi;) = exp(&beta;<sub>A</sub><&phi;>/M<sub>Pl</sub>) &times; I<sub>0</sub>(&beta;<sub>A</sub> &delta;&phi;/M<sub>Pl</sub>),
(23)

with &delta;&phi;<sup>2</sup> &prop; 1/N<sub>eff</sub>. The transfer function mapping the single-particle to the many-body limit is:

T(&rho;) = (&rho;<sub>T</sub>/&rho;<sub>core</sub>) [1 + (&rho;/&rho;<sub>T</sub>)<sup>2/3</sup>].
(24)

Here T<sub>num</sub>(&rho;) = (&rho;/&rho;<sub>core</sub>) S(&rho;) denotes the numerical solver transfer function, while T<sub>bridge</sub>(&rho;) denotes the analytic dimensional bridge ansatz used for interpretation. They are not identical objects; the bridge ansatz is a phenomenological closure form whose derivation from the full non-linear field equations is not pursued here, as the self-consistent density-dependent coherence length (Eq. 7g) provides the physical closure.

**Astrophysical saturation capping.** The chameleon curvature equation (7d) in Section 2.3 couples the scalar mass to the matter density &rho;<sub>m</sub> through the conformal factor A<sup>4</sup>(&phi;). In the many-body regime, the density sourcing the scalar field is not the raw baryonic &rho;<sub>m</sub> but the effective coherently participating density &rho;<sub>eff</sub> = &rho;<sub>m</sub> &times; f<sub>RP</sub><sup>(corr)</sup> (Eq. 7f). The random-phase suppression factor f<sub>RP</sub><sup>(corr)</sup> = (1/&radic;N<sub>eff</sub>)&radic;1 + nV<sub>c</sub> (Eq. 28) decreases as N<sub>eff</sub> grows, so that in degenerate matter (neutron star cores, &rho; &sim; 10<sup>15</sup> g/cm<sup>3</sup>, N<sub>eff</sub> &sim; 10<sup>60</sup>) the product &rho;<sub>eff</sub> = &rho;<sub>m</sub> &times; f<sub>RP</sub><sup>(corr)</sup> saturates at &rho;<sub>eff</sub> &rarr; &rho;<sub>T</sub> &times; O(1) &asymp; 20 g/cm<sup>3</sup>. The scalar mass m<sub>&phi;</sub><sup>2</sup> = V''<sub>eff</sub>(&phi;<sub>eff</sub>(&rho;<sub>eff</sub>)) is therefore bounded from above by its value at the saturation scale, protecting the geometric closure &chi; = 1/&radic;2 and the screening length &lambda;<sub>scr</sub> = &radic;2 &hbar;/(m<sub>e</sub>c) from dynamic runaway. The Weak Equivalence Principle is preserved because &rho;<sub>eff</sub> at saturation is composition-independent: all degenerate bodies source the scalar field at the same plateau &rho;<sub>T</sub>, regardless of their raw baryonic density.

The Thomas-Fermi-TEP numerical solver (`scripts/steps/step_07_numeric_screening_transfer.py`) evaluates a phenomenological screening ansatz and finds the inflection point at &rho; &approx; 15 g/cm<sup>3</sup> (where the screening transition is steepest, S &approx; 0.65).

**Derivation: tanh inflection.** For S(&rho;) = tanh(&rho;/&rho;<sub>T</sub>), the slope with respect to logarithmic density is

dS / d ln &rho; = &rho; dS/d&rho; = x sech<sup>2</sup> x,   x = &rho;/&rho;<sub>T</sub>.
(25)

The maximum satisfies d/dx (x sech<sup>2</sup> x) = 0, so

sech<sup>2</sup> x &minus; 2x sech<sup>2</sup> x tanh x = 0  &Rightarrow;  1 &minus; 2x tanh x = 0.
(26)

Hence x tanh x = 1/2. Numerically, x &asymp; 0.77. Therefore the inflection (maximum log-slope) occurs at

&rho; &asymp; 0.77 &rho;<sub>T</sub>.
(27)

This is a structural feature of the tanh ansatz, not a derived prediction of the physical crossover density. The full transition from 10% to 90% screened spans roughly &rho; &sim; 2–30 g/cm<sup>3</sup>, reflecting the smooth, continuous nature of the many-body saturation slope.

### A.4 Correlation-Modified Random-Phase Suppression

The naive random-phase bound assumes uncorrelated topological charge orientations. For correlated orientations with pair correlator C(r<sub>ij</sub>) = &lang;s<sub>i</sub> s<sub>j</sub>&rang;, the second moment of the collective scalar field generalises to:

&lang;&phi;<sup>2</sup>&rang; = N &phi;<sub>0</sub><sup>2</sup> [1 + nV<sub>c</sub>],
(28)

where V<sub>c</sub> = &int; d<sup>3</sup>r &thinsp; C(r) is the correlation volume. The suppression factor becomes:

f<sub>RP</sub><sup>(corr)</sup> = (1 / &radic;N) &radic;1 + nV<sub>c</sub>.
(29)

For an exponential correlator C(r) = exp(&minus;r/&xi;), V<sub>c</sub> = 8&pi;&xi;<sup>3</sup>. At &rho; &approx; 20 g/cm<sup>3</sup> with &xi; &sim; &lambda;<sub>F</sub> &approx; 10<sup>-10</sup> m, the correction factor is &radic;1 + 8&pi;n&xi;<sup>3</sup> &approx; 12. Even a factor of 10 weakening leaves f<sub>RP</sub> &sim; 10<sup>-24</sup> for N<sub>eff</sub> &sim; 10<sup>50</sup>, reducing the mean-field/phenomenological discrepancy by approximately 25 orders of magnitude (from &sim;10<sup>47</sup> to &sim;10<sup>22</sup> g/cm<sup>3</sup>). The remaining &sim;21 orders of magnitude require nonlinear closure that remains unresolved. The correlated-phase solver output is physically bounded by the coherent limit (f &le; f<sub>coh</sub>), preventing unphysical suppression factors exceeding unity.

### A.5 &chi; Convergence Constant

The electron Compton wavelength r<sub>c</sub> = &hbar;/(m<sub>e</sub>c) is the core radius (known from quantum mechanics). The scalar Yukawa screening length follows from choosing the scalar mass such that the Compton wavelength of the scalar field matches the electron Compton scale up to the factor required by the single-particle Klein-Gordon closure: m<sub>&phi;</sub> = m<sub>e</sub>/&radic;2. Then

&lambda;<sub>scr</sub> = &hbar;/(m<sub>&phi;</sub>c) = &radic;2 &hbar;/(m<sub>e</sub>c) = &radic;2 r<sub>c</sub>.
(30)

Their ratio is the internal geometric consistency condition:

&chi; = r<sub>c</sub> / &lambda;<sub>scr</sub> = 1/&radic;2 &approx; 0.707.
(31)

This is not an independent empirical prediction of a new length scale; it is a model-closure relation linking the Compton-scale finite core to the scalar screening length required by the proposed topological charge geometry. The non-trivial content is that this closure is self-consistent and anchors the model to the known electron Compton wavelength.

### A.6 g&minus;2 Temporal Topology Drag

**Candidate expression (operator-level result from FW reduction).** A uniform conformal rescaling cannot alter the dimensionless magnetic anomaly a<sub>&mu;</sub> = (g&minus;2)/2, because the same local clock factor rescales both spin-precession and cyclotron frequencies. Any TEP contribution must arise from non-uniform temporal shear, orientation-bundle curvature, disformal transport, or synchronization holonomy. The disformal Dirac operator generates a candidate Pauli-like spin coupling via the Foldy–Wouthuysen reduction in Appendix A.17; the expression below is the candidate form. The operator-level result is strong (specific channel breaking conformal degeneracy); the normalization is not yet fully derived.

TEP contribution to the anomaly:

a<sub>&mu;</sub><sup>TEP</sup> &sim; a<sub>&mu;</sub><sup>SM</sup> &Delta;A<sub>hol</sub> / A<sub>&infin;</sub>,
(32)

where &Delta;A<sub>hol</sub> is the apparatus-integrated non-uniform temporal-topology contribution.

### A.7 Spinorial Holonomy and Exchange Compatibility

Compact phase holonomy under a 2&pi; rotation:

&Delta;&phi; = 2&pi;n (&beta;<sub>A</sub>/M<sub>Pl</sub>).
(33)

The real conformal factor A(&phi;) = exp(&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>) is not periodic. Periodicity belongs strictly to the compact phase associated with the defect. Single-valuedness of the matter-frame metric requires A(&phi;) to be smooth and single-valued, while fermionic spin requires spinor holonomy in the temporal-orientation bundle. Fermionic statistics correspond to the minimal non-trivial integer winding represented in the spinorial double cover:

n = &pm;1,   S<sub>z</sub> = &pm;&hbar;/2.
(34)

This gives half-integer spin consistency from the topology of the temporal shear defect without assigning half-integer winding to the scalar field itself. The orientation-bundle construction reproduces the single-defect spinorial transformation law. For multi-defect states, antisymmetry is presently imposed through the standard fermionic Hilbert-space construction. Deriving exchange statistics directly from TEP defect dynamics remains open. The braid-group holonomy for adiabatic vortex exchange in 2+1D and the 3+1D multi-defect interaction Hamiltonian are developed in Appendix A.19.

**Proposition 2 (Spinorial lift).** Given a compact orientation phase &vartheta; with integer winding

&oint; d&vartheta; = 2&pi;n,
(35)

and a spinorial lift where the physical spinor transforms as &psi; &mapsto; e<sup>i&vartheta;/2</sup> &psi;, then for n = 1:

&vartheta; = 2&pi;  &Rightarrow;  &psi; &mapsto; e<sup>i&pi;</sup> &psi; = &minus;&psi;,
(36)

and for n = 2:

&vartheta; = 4&pi;  &Rightarrow;  &psi; &mapsto; e<sup>i2&pi;</sup> &psi; = +&psi;.
(37)

Conditional on the temporal-orientation bundle (Proposition 1, Section 2.2) and the stated lift to the SU(2) spinor representation, the 2&pi; sign reversal and 4&pi; return follow rigorously. This establishes an internally consistent single-defect spinorial geometry. For multi-defect states, antisymmetry is presently imposed through the standard fermionic Hilbert-space construction; a TEP-native derivation from the bundle dynamics remains open. The braid-group holonomy for adiabatic vortex exchange in 2+1D and the 3+1D multi-defect interaction Hamiltonian are developed in Appendix A.19.

### A.8 Disformal Inverse Metric and Null-Cone Tilt

For the rank-one disformal ansatz g̃<sub>&mu;&nu;</sub> = A<sup>2</sup>&eta;<sub>&mu;&nu;</sub> + B u<sub>&mu;</sub> u<sub>&nu;</sub> with u<sub>&mu;</sub> = &nabla;<sub>&mu;</sub>&phi;, the inverse metric is derived via the Sherman–Morrison formula for symmetric rank-one updates:

g̃<sup>&mu;&nu;</sup> = A<sup>&minus;2</sup> &eta;<sup>&mu;&nu;</sup> &minus; (B / A<sup>2</sup>(A<sup>2</sup> + B u<sup>2</sup>)) u<sup>&mu;</sup> u<sup>&nu;</sup>,
(38)

where u<sup>2</sup> = &eta;<sup>&alpha;&beta;</sup> u<sub>&alpha;</sub> u<sub>&beta;</sub>. Symbolic verification against direct 4&times;4 matrix inversion yields zero difference. For propagation along u<sub>&mu;</sub>, the null-cone condition gives the effective speed

v<sub>eff</sub><sup>2</sup> = c<sup>2</sup> A<sup>2</sup> / (A<sup>2</sup> + B u<sup>2</sup>),
(39)

and the effective refractive index is n<sub>eff</sub> = c / v<sub>eff</sub> = &radic;(1 + B u<sup>2</sup> / A<sup>2</sup>). Full Christoffel-symbol and geodesic derivation in 4D with arbitrary u<sub>&mu;</sub> orientation is reserved for TEP-KIN (Paper 25).

### A.9 Weyl-Rescaled Dirac Operator

Under conformal rescaling g̃<sub>&mu;&nu;</sub> = A<sup>2</sup> &eta;<sub>&mu;&nu;</sub>, the Dirac operator transforms with the known weight D̃ = A<sup>&minus;5/2</sup> D A<sup>3/2</sup>. Expanding the left-hand side with the product rule on A<sup>3/2</sup> &psi; yields

D̃ &psi; = A<sup>&minus;1</sup> [ i &gamma;<sup>&mu;</sup> (&part;<sub>&mu;</sub> + (3/2) &Sigma;<sub>&mu;</sub>) &psi; &minus; m &psi; ],
(40)

where &Sigma;<sub>&mu;</sub> = &part;<sub>&mu;</sub> ln A = &nabla;<sub>&mu;</sub> ln A. Symbolic verification (non-commutative &gamma;, &psi;) confirms zero difference between the expanded operator identity and the explicitly coupled form. The overall prefactor A<sup>&minus;1</sup> is an operator normalization from the conformal transformation, not a mass rescaling. After multiplying through by A, the flat-frame Dirac equation is [i&gamma;<sup>&mu;</sup>(&part;<sub>&mu;</sub> + (3/2)&Sigma;<sub>&mu;</sub>) &minus; mA]&psi; = 0, giving the physical effective mass m<sub>*</sub> = mA consistent with the Klein–Gordon/Hamilton–Jacobi derivation (Section 2.1). Reference: TEP-QF Paper 23, &sect;3.2.1.

### A.10 Disformal Christoffel Symbols and Null Geodesic (1+1D)

For the static 1+1D disformal metric g̃<sub>00</sub> = A<sup>2</sup> + B u<sub>0</sub><sup>2</sup>, g̃<sub>11</sub> = &minus;A<sup>2</sup> + B u<sub>1</sub><sup>2</sup>, g̃<sub>01</sub> = B u<sub>0</sub> u<sub>1</sub>, the Christoffel symbols are

&Gamma;<sup>x</sup><sub>tt</sub> = (1/2) g<sup>xx</sup> (&minus;&part;<sub>x</sub> g<sub>tt</sub>),   &Gamma;<sup>t</sup><sub>tx</sub> = (1/2) g<sup>tt</sup> (&part;<sub>x</sub> g<sub>tt</sub>),
(41)

with the remaining components following analogously. The null condition g<sub>tt</sub> + 2 g<sub>tx</sub> v + g<sub>xx</sub> v<sup>2</sup> = 0 is solved for v = dx/dt; the forward-propagating solution yields the same effective refractive index as the null-cone tilt derivation, n<sub>eff</sub> = &radic;(A<sup>2</sup> + B(u<sub>0</sub><sup>2</sup> &minus; u<sub>1</sub><sup>2</sup>)) / A. Full 4D disformal geodesic with arbitrary u<sub>&mu;</sub> orientation is reserved for TEP-KIN.

### A.11 Bianchi Identity for U(1) Bundle Curvature

For a smooth azimuthal gauge field A<sub>&theta;</sub> = &Phi;r / (2&pi;), the field strength F<sub>r&theta;</sub> = &part;<sub>r</sub> A<sub>&theta;</sub> is non-zero and smooth. The Bianchi identity component &part;<sub>r</sub> F<sub>&theta;z</sub> + &part;<sub>&theta;</sub> F<sub>zr</sub> + &part;<sub>z</sub> F<sub>r&theta;</sub> simplifies identically to zero. For the vortex gauge field A<sub>&theta;</sub> = &Phi; / (2&pi;) (constant in r), F<sub>r&theta;</sub> = 0 for r > 0; the Bianchi identity again gives zero everywhere away from the singular core. The identity fails only at r = 0 where the delta-function flux resides. Reference: TEP-KIN Paper 25, Proposition 1.

### A.12 Data Provenance

| Dataset | Source | Points | File |
| --- | --- | --- | --- |
| PRad 1.1 GeV | Xiong *et al.* (2019) | 33 | 1.1GeV_table_normGE.txt |
| PRad 2.2 GeV | Xiong *et al.* (2019) | 38 | 2.2GeV_table_normGE.txt |
| A1 Cross Sections | Bernauer *et al.* (2014) | 1,422 | a1_cross_sections.dat |
| Total |  | 1,493 |  |

Table A.12: Data provenance for the JLab/AMBER cross-section prediction pipeline.

### A.13 Exchange Compatibility from 4&pi; Periodicity

The single-defect spinorial lift (Appendix A.7) establishes 4&pi; periodicity for the minimal winding n = 1:

&psi;(&vartheta; + 4&pi;) = +&psi;(&vartheta;),   &psi;(&vartheta; + 2&pi;) = &minus;&psi;(&vartheta;).
(42)

This is the defining property of spin-1/2. By the spin-statistics theorem, spin-1/2 particles are fermions. The antisymmetrized two-particle spinor state is constructed directly:

&Psi;<sub>A</sub> = (&chi;<sub>1</sub> &otimes; &chi;<sub>2</sub> &minus; &chi;<sub>2</sub> &otimes; &chi;<sub>1</sub>) / &radic;2.
(43)

Under particle exchange &chi;<sub>1</sub> &leftrightarrow; &chi;<sub>2</sub>, the state transforms as &Psi;<sub>A</sub> &mapsto; &minus;&Psi;<sub>A</sub>. Symbolic verification with non-commutative spinor symbols confirms the exchange difference is identically zero: (&chi;<sub>1</sub>&chi;<sub>2</sub> &minus; &chi;<sub>2</sub>&chi;<sub>1</sub>) + (&chi;<sub>2</sub>&chi;<sub>1</sub> &minus; &chi;<sub>1</sub>&chi;<sub>2</sub>) = 0. This verifies that an already-antisymmetrized state changes sign under exchange by construction. The spin-statistics theorem connects the single-defect 4&pi; periodicity to fermionic exchange, but a TEP-native derivation of multi-defect antisymmetry from the bundle dynamics — rather than imposition through the standard Hilbert-space construction — remains open. The braid-group holonomy for adiabatic vortex exchange in 2+1D and the 3+1D multi-defect interaction Hamiltonian are developed in Appendix A.19.

### A.14 BMT Spin Precession and g&minus;2 Selection Rule

In the matter-frame conformal metric g&#771;<sub>&mu;&nu;</sub> = A<sup>2</sup>&eta;<sub>&mu;&nu;</sub>, the effective mass is m<sub>*</sub> = mA (from the Klein–Gordon/Hamilton–Jacobi equation, Section 2.1; consistent with the Weyl-rescaled Dirac operator after multiplying by A, Appendix A.9). The physical magnetic field in the matter frame transforms as B<sub>matter</sub> = B/A<sup>2</sup>. The cyclotron and spin-precession frequencies are:

&omega;<sub>c</sub> = eB<sub>matter</sub> / m<sub>*</sub> = eB / (mA),   &omega;<sub>s</sub> = egB<sub>matter</sub> / (2m<sub>*</sub>) = egB / (2mA).
(44)

Their ratio is &omega;<sub>s</sub>/&omega;<sub>c</sub> = g/2, independent of A. For a Dirac particle g = 2, the anomaly a = (g&minus;2)/2 = 0. This proves the selection rule: a *uniform* conformal rescaling cannot generate a magnetic anomaly. Any TEP contribution must come from non-uniform temporal shear. The gradient correction from the &Sigma;<sub>&mu;</sub> = &part;<sub>&mu;</sub> ln A coupling in the Dirac operator gives an additional spin-precession term &Delta;&omega;<sub>s</sub> ~ (3/2)&Sigma;, which affects the spin but not the orbital motion, thereby generating an effective anomaly. The full disformal Dirac–Pauli term derivation and quantitative modulation amplitude prediction are given in Appendix A.17.

### A.15 Finite-Core Regulator and QED Recovery

The unregulated QED vacuum polarization diverges logarithmically in the UV. The Gaussian finite-core regulator

R(k) = exp(&minus;k<sup>2</sup>&sigma;<sup>2</sup>)
(45)

suppresses high momenta, rendering the integral finite:

&int;<sub>0</sub><sup>&infin;</sup> exp(&minus;k<sup>2</sup>&sigma;<sup>2</sup>) dk = &radic;&pi; / (2&sigma;).
(46)

The regulator width is set by the topological charge core radius. In natural units (ℏ = c = 1), r<sub>c</sub> = 1/m, so

&sigma; = r<sub>c</sub> / &radic;2 = 1/(&radic;2 m).
(47)

The UV suppression is verified symbolically: lim<sub>k&rarr;&infin;</sub> exp(&minus;k<sup>2</sup>&sigma;<sup>2</sup>) = 0. The QED point-particle limit is recovered as &sigma; &rarr; 0: lim<sub>&sigma;&rarr;0</sub> exp(&minus;k<sup>2</sup>&sigma;<sup>2</sup>) = 1. The finite-core regulator therefore interpolates between a finite, Compton-scale regularized theory and standard QED.

The Gaussian form R(k) = exp(&minus;k<sup>2</sup>&sigma;<sup>2</sup>) acts as an analytic proxy for the exact spatial profile of the temporal defect. In the full non-linear theory, this regulator is replaced by the exact Fourier transform of the topological winding profile (analogous to a Nielsen–Olesen vortex core). The structural conclusion — that loop integrals are truncated at k<sub>max</sub> &sim; 1/r<sub>c</sub> — remains invariant regardless of the specific profile geometry, as all valid defect profiles vanish outside the Compton core boundary. However, the claim that this truncation bounds the running coupling to a finite &alpha;<sub>max</sub> is not established: a proper beta function requires extraction from the renormalized vacuum-polarization dependence on the external momentum scale, not a standalone loop-momentum integral.

**Candidate nonlocal UV regulator.** The finite-core profile motivates a candidate nonlocal UV regulator whose running-coupling behaviour requires a complete gauge-covariant vacuum-polarization and matching calculation. A beta-function derivation from this regulator was attempted but found to be invalid: the running must be extracted from the renormalized vacuum-polarization dependence on the external momentum scale, not by identifying the beta function with a standalone loop-momentum integral. The detailed diagnostics are documented in the repository audit (`scripts/audit/step_04_internal_external_regulator_split.py`, `scripts/audit/step_06_claim_evidence_registry.py`). The Ward identity analysis and Lamb shift estimate with this regulator are given in Appendix A.18.

### A.16 SymPy Derivation Inventory

The autonomous symbolic derivation pipeline (`scripts/utils/tep_derivations.py`) computes the following results. The conformal Hamilton–Jacobi sector follows from axioms A1–A3, while the spin/vorticity sector follows from the compact temporal-orientation bundle introduced as Proposition 1 from the same axioms, the experimental fact of spinorial holonomy, and the TEP ontological principle.

- g&#771;-Hamilton-Jacobi equation with conformal factor A(&phi;) = exp(&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>)

- Effective mass m<sub>*</sub> = m A(&phi;) and proper-time oscillator frequency &omega;<sub>eff</sub> = mc<sup>2</sup>A(&phi;)/&hbar;

- Proposition 1: Conformal exactness &oint; d ln A = 0 forces spin into the orientation bundle

- Topological charge azimuthal phase shear K<sub>&theta;</sub> = n/r

- Quantized compact phase circulation &Gamma; = 2&pi;n

- Proposition 2: Spinorial lift of temporal-orientation winding; 2&pi; &rArr; &minus;&psi;, 4&pi; &rArr; +&psi;

- Topological spin origin from compact phase holonomy; single-defect spinorial transformation established; multi-defect antisymmetry imposed through standard fermionic Hilbert-space construction, not yet derived from bundle dynamics (Derivation 18); braid-group holonomy in 2+1D and multi-defect interaction Hamiltonian developed in Appendix A.19 (Derivation 26)

- Model-closure condition &chi; = r<sub>c</sub>/&lambda;<sub>scr</sub> = 1/&radic;2 with explicit m<sub>&phi;</sub> = m<sub>e</sub>/&radic;2 derivation

- Single-particle core density &rho;<sub>core</sub> &sim; m<sub>e</sub><sup>4</sup>c<sup>3</sup>/&hbar;<sup>3</sup> (white-dwarf-scale, ~10<sup>4</sup> g/cm<sup>3</sup>)

- Correlation-modified random-phase suppression: f<sub>RP</sub><sup>(corr)</sup> = (1/&radic;N) &radic;1 + nV<sub>c</sub> with V<sub>c</sub> = 8&pi;&xi;<sup>3</sup> for exponential correlator; physically bounded by &xi; &lesssim; &lambda;<sub>F</sub>, giving O(1)–O(10) correction to the 1/&radic;N<sub>eff</sub> suppression

- g&minus;2 selection rule: uniform conformal rescaling cannot generate an anomalous magnetic moment

- g&minus;2 candidate operator: disformal Dirac operator generates candidate Pauli-like spin coupling; a<sub>&mu;</sub><sup>TEP</sup> &sim; a<sub>&mu;</sub><sup>SM</sup> &Delta;A<sub>hol</sub> / A<sub>&infin;</sub>; operator-level result is strong, normalization not yet fully derived

- Proton form-factor expansion: &langle;r<sup>2</sup>&rangle;<sub>eff</sub> = &langle;r<sup>2</sup>&rangle;<sub>dipole</sub> &minus; 6&delta;<sub>A</sub>/Q<sub>c</sub><sup>2</sup> derived from A<sub>TEP</sub>(Q<sup>2</sup>) &times; F<sub>dipole</sub>(Q<sup>2</sup>)

- Toy finite-core UV regulator: Gaussian form factor renders loop integral finite

- Disformal inverse metric and null-cone tilt: exact Sherman-Morrison derivation for rank-one update g̃<sub>&mu;&nu;</sub> = A<sup>2</sup>&eta;<sub>&mu;&nu;</sub> + B u<sub>&mu;</sub> u<sub>&nu;</sub>, with inverse verified against direct symbolic inversion and null-cone tilt extracted (TEP-KIN signpost)

- Weyl-rescaled Dirac operator: D&#771; = A<sup>&minus;5/2</sup> D A<sup>3/2</sup> verified algebraically, with (3/2) &Sigma;<sub>&mu;</sub> coupling and effective mass m<sub>*</sub> = mA after multiplying through by A (TEP-QF &sect;3.2.1)

- Disformal Christoffel symbols (1+1D static): &Gamma;<sup>&lambda;</sup><sub>&mu;&nu;</sub> derived from g̃<sub>&mu;&nu;</sub> = A<sup>2</sup>&eta;<sub>&mu;&nu;</sub> + B u<sub>&mu;</sub> u<sub>&nu;</sub>, null geodesic solved, effective refractive index n<sub>eff</sub> = &radic;(A<sup>2</sup> + B(u<sub>0</sub><sup>2</sup> &minus; u<sub>1</sub><sup>2</sup>)) / A extracted (TEP-KIN &sect;2.1)

- Bianchi identity dF = 0 for U(1) bundle curvature: verified for smooth gauge field; shown to hold away from vortex core where singular delta-function flux resides (TEP-KIN Proposition 1)

- Exchange antisymmetry: 4&pi; periodicity proven for n=1; antisymmetrized two-particle state verified to change sign under exchange by construction; TEP-native derivation from bundle dynamics remains open (Derivation 18)

- BMT spin precession in conformal metric: &omega;<sub>c</sub> = eB/(mA), &omega;<sub>s</sub> = egB/(2mA); ratio &omega;<sub>s</sub>/&omega;<sub>c</sub> = g/2; selection rule proven for uniform A (no anomaly); gradient correction from &Sigma; = &nabla; ln A formalized (Derivation 19)

- Finite-core QED recovery: Gaussian regulator exp(&minus;k<sup>2</sup>&sigma;<sup>2</sup>) renders vacuum polarization finite with &int;<sub>0</sub><sup>&infin;</sup> exp(&minus;k<sup>2</sup>&sigma;<sup>2</sup>) dk = &radic;&pi;/(2&sigma;); regulator width &sigma; = r<sub>c</sub>/&radic;2 relates to Compton core radius; QED divergence recovered as &sigma; &rarr; 0; a beta-function derivation from this regulator was attempted but found invalid (repository audit `scripts/audit/`) (Derivation 20)

- Proposition 1 (Necessity of non-scalar holonomy sector): compact U(1) temporal-orientation bundle motivated by axioms A1–A3 + experimental spinorial holonomy + TEP ontological principle; U(1) character, SU(2) spinorial lift, and identification with electromagnetism are model assumptions; minimal gauge-invariant bundle action (Eq. 4a) with Maxwell kinetic + Stueckelberg/Higgs-type coupling (Derivation 21)

- Scalar mass self-consistency: m<sub>&phi;</sub><sup>2</sup> = V''<sub>eff</sub>(&phi;<sub>eff</sub>) from chameleon potential curvature at core density; Compton-scale closure m<sub>&phi;</sub> = m<sub>e</sub>/&radic;2 is a self-consistency condition on potential parameters, not a first-principles prediction (Derivation 22)

- Screening ansatz: tanh S(&rho;) = tanh(&rho;/&rho;<sub>T</sub>) adopted as minimal smooth interpolation consistent with thin-shell asymptotes; crossover &rho;<sub>T</sub> set by minimum condition (7b); functional form is phenomenological, not derived from field equations (Derivation 23)

- g&minus;2 candidate operator: disformal Dirac operator generates a candidate Pauli-like spin coupling via Foldy–Wouthuysen reduction (Eqs. 49a–49c); operator-level result is strong (specific channel breaking conformal degeneracy); normalization is not yet fully derived; detailed formulae and diagnostics in repository audit (`scripts/audit/`) (Derivation 24)

- Ward identity with finite-core regulator: naive Gaussian regulator breaks Ward identity at O(q<sup>2</sup>&sigma;<sup>2</sup>); gauge-covariant Wilson-line formulation outlines a possible completion, but exact transversality to all orders is not yet established; detailed diagnostics in repository audit (`scripts/audit/`) (Derivation 25)

- Braid-group holonomy: 2+1D abelian anyon representation with statistical angle &theta; = &pi; (fermionic); 3+1D multi-defect Hamiltonian with Yukawa interactions and Slater-determinant statistics consistent with standard fermionic quantum mechanics (Derivation 26)

Full outputs are serialized in `results/tep_derivations.json`.

### A.17 Disformal (g&minus;2) Candidate Operator and Modulation Template

The BMT selection rule (Appendix A.14) proves that a uniform conformal rescaling cannot generate an anomalous magnetic moment. The missing ingredient is the disformal sector B(&phi;), which introduces non-uniform temporal shear that breaks the degeneracy between spin-precession and cyclotron frequencies. Here we construct the disformal Dirac operator, perform the semiclassical (Foldy–Wouthuysen) reduction, extract a candidate Dirac–Pauli term, and compute the anomalous magnetic moment. Several gaps in this calculation are acknowledged below.

**Setup.** The matter-frame Dirac equation in the full disformal metric g&#771;<sub>&mu;&nu;</sub> = A<sup>2</sup>&eta;<sub>&mu;&nu;</sub> + B u<sub>&mu;</sub>u<sub>&nu;</sub> (with u<sub>&mu;</sub> = &nabla;<sub>&mu;</sub>&phi;) is obtained by replacing the flat gamma matrices with the tetrad g&#771;<sup>&mu;</sup><sub>(a)</sub> associated with g&#771;<sub>&mu;&nu;</sub>. The inverse metric (Appendix A.8, Eq. 39) gives:

g&#771;<sup>&mu;&nu;</sup> = A<sup>&minus;2</sup>&eta;<sup>&mu;&nu;</sup> &minus; (B / A<sup>2</sup>(A<sup>2</sup> + Bu<sup>2</sup>)) u<sup>&mu;</sup>u<sup>&nu;</sup>.
(48)

The Dirac operator in the disformal metric is D&#771; = i&gamma;<sup>(a)</sup> e<sup>&mu;</sup><sub>(a)</sub> D<sub>&mu;</sub> &minus; m, where e<sup>&mu;</sup><sub>(a)</sub> is the tetrad and D<sub>&mu;</sub> = &part;<sub>&mu;</sub> + &Omega;<sub>&mu;</sub> is the spin covariant derivative with connection &Omega;<sub>&mu;</sub> = (1/4)&omega;<sub>&mu;ab</sub>&gamma;<sup>a</sup>&gamma;<sup>b</sup>. Expanding to first order in B/A<sup>2</sup> (weak disformal limit), the Dirac equation acquires an additional coupling:

[i&gamma;<sup>&mu;</sup>(&part;<sub>&mu;</sub> + (3/2)&Sigma;<sub>&mu;</sub>) &minus; mA + i(&beta;<sub>B</sub>/2A<sup>2</sup>)&gamma;<sup>&mu;</sup>u<sub>&mu;</sub>u<sup>&nu;</sup>D<sub>&nu;</sub>] &psi; = 0,
(49)

where &beta;<sub>B</sub> = B/(A<sup>2</sup> + Bu<sup>2</sup>) &approx; B/A<sup>2</sup> in the weak-disformal limit. The first two terms are the Weyl-rescaled Dirac operator (Appendix A.9). The third term is the disformal correction: it couples the spinor to the scalar gradient u<sub>&mu;</sub> = &nabla;<sub>&mu;</sub>&phi; through a gamma-matrix contraction.

**Foldy–Wouthuysen reduction.** To extract the non-relativistic Hamiltonian, we perform the standard FW transformation on the disformal Dirac operator (49). Writing the Dirac equation as i&part;<sub>0</sub>&psi; = H&psi; with H = &alpha;&middot;p + &beta;mA + (3/2)&alpha;&middot;&Sigma; + H<sub>dis</sub>, where &alpha; = &gamma;<sup>0</sup>&gamma;<sup>i</sup>, &beta; = &gamma;<sup>0</sup>, and H<sub>dis</sub> = &minus;(&beta;<sub>B</sub>/2A<sup>2</sup>)&gamma;<sup>0</sup>&gamma;<sup>&mu;</sup>u<sub>&mu;</sub>u<sup>&nu;</sup>D<sub>&nu;</sub> is the disformal correction, we apply the FW block-diagonalization procedure. The odd (off-diagonal) operators are eliminated iteratively to order 1/m<sup>2</sup>.

In the presence of an external electromagnetic field A<sub>&mu;</sub> = (&Phi;, **A**), the minimal coupling p<sub>&mu;</sub> &rarr; p<sub>&mu;</sub> &minus; eA<sub>&mu;</sub> is applied to the Weyl-rescaled sector. The FW Hamiltonian to order 1/m<sup>2</sup> in the rest frame (A &asymp; 1, weak field) is:

H<sub>FW</sub> = &beta;(mA + p<sup>2</sup>/(2mA) &minus; p<sup>4</sup>/(8m<sup>3</sup>A<sup>3</sup>)) + e&Phi; &minus; (e&beta;/2mA)&sigma;&middot;B &minus; (e/4m<sup>2</sup>A<sup>2</sup>)&sigma;&middot;E&times;p &minus; (e/8m<sup>2</sup>A<sup>2</sup>)&nabla;&middot;E + H<sub>dis</sub><sup>(FW)</sup>,
(49a)

where B = &nabla;&times;**A** is the magnetic field, E = &minus;&nabla;&Phi; &minus; &part;<sub>0</sub>**A** is the electric field, and &sigma;<sup>i</sup> = &epsilon;<sup>ijk</sup>&Sigma;<sub>jk</sub> are the Pauli matrices. The first four terms are the standard FW Hamiltonian (kinetic energy, Zeeman term with electric charge e, spin-orbit, Darwin). The disformal FW correction H<sub>dis</sub><sup>(FW)</sup> is obtained by block-diagonalizing H<sub>dis</sub>:

H<sub>dis</sub><sup>(FW)</sup> = &minus;(&beta;<sub>B</sub>/2A<sup>2</sup>) [ &beta;u<sub>0</sub><sup>2</sup>(mA) + (eu<sub>0</sub>/2mA)&sigma;&middot;B + (u<sub>0</sub><sup>2</sup>/4m<sup>2</sup>A<sup>2</sup>){&sigma;&middot;p, &sigma;&middot;E} + ... ],
(49b)

where we used u<sub>i</sub> << u<sub>0</sub> (the temporal component dominates in the lab frame) and retained terms to O(1/m<sup>2</sup>). The first term is a disformal mass shift &delta;m = &minus;(&beta;<sub>B</sub>/2A<sup>2</sup>)u<sub>0</sub><sup>2</sup>(mA), which renormalizes the rest mass but does not affect g&minus;2. The second term is the key result: a disformal Zeeman correction. Note that the disformal metric depends quadratically on u<sub>&mu;</sub> (through B u<sub>&mu;</sub> u<sub>&nu;</sub>), so the appearance of a term linear in u<sub>0</sub> in the Zeeman correction requires explicit power-counting justification that is not yet provided.

**Extraction of the Dirac–Pauli term.** The standard Zeeman term is &minus;(e&beta;/2mA)&sigma;&middot;B with g = 2. The disformal Zeeman correction is &minus;(&beta;<sub>B</sub>/2A<sup>2</sup>)(eu<sub>0</sub>/2mA)&sigma;&middot;B. Combining, the total Zeeman coupling is:

H<sub>Zeeman</sub> = &minus;(e&beta;/2mA) [ 1 + (&beta;<sub>B</sub>/A)u<sub>0</sub> ] &sigma;&middot;B.
(49c)

The effective g-factor is therefore g<sub>eff</sub> = 2[1 + (&beta;<sub>B</sub>/A)u<sub>0</sub>], giving a candidate anomalous magnetic moment a<sub>&mu;</sub><sup>TEP</sup> = (&beta;<sub>B</sub>/A)u<sub>0</sub>. This is a structural result: the disformal Dirac operator generates a candidate Pauli-like spin coupling through a specific operator channel that breaks the conformal degeneracy. However, the normalization is not yet fully derived. The following algebraic issues remain unresolved: (i) the disformal metric depends quadratically on u<sub>&mu;</sub> (through B u<sub>&mu;</sub> u<sub>&nu;</sub>), so the appearance of a term linear in u<sub>0</sub> in the Zeeman correction requires explicit power-counting justification; (ii) the mass dimensions of B(&phi;), &beta;<sub>B</sub>, u<sub>&mu;</sub> and their product must be stated explicitly to verify that a<sub>&mu;</sub><sup>TEP</sup> is dimensionless; (iii) the FW reduction assumes u<sub>i</sub> << u<sub>0</sub> (temporal component dominates), but the proposed annual modulation is attributed to motion through a spatial cosmic gradient — these assumptions may be in tension; (iv) a disformal metric modifies orbital as well as spin dynamics, so the cyclotron frequency correction must be derived from the same Hamiltonian rather than assumed absent; (v) the identification of u<sub>0</sub> with the apparatus-integrated holonomy &Delta;A<sub>hol</sub>/A<sub>&infin;</sub> contains an inconsistency in Planck/coupling factors. The detailed formulae and diagnostics are documented in the repository audit (`scripts/audit/step_05_complete_disformal_fw_audit.py`, `scripts/audit/step_06_claim_evidence_registry.py`) rather than presented here, because a published paper should not retain a known dimensionally inconsistent normalization alongside its numerical output.

**Qualitative modulation structure.** The diurnal and annual modulations arise because &Delta;A<sub>hol</sub>(t) varies with Earth's rotational and orbital velocity through the cosmic temporal shear gradient &nabla;&phi;<sub>cosmic</sub>. The diurnal component (period 24 h) arises from Earth's rotation (v<sub>rot</sub> &sim; 465 m/s at the equator); the annual component (period 365.25 d) arises from Earth's orbital velocity (v<sub>orb</sub> &sim; 30 km/s). The ratio of diurnal to annual modulation amplitudes is v<sub>rot</sub>/v<sub>orb</sub> &sim; 0.015, making the annual signal the dominant modulation. The modulation is phased to Earth's orbital velocity through the cosmic shear gradient. The absolute amplitude &sim;10<sup>&minus;13</sup> is a scaling ansatz (static benchmark &times; v<sub>orb</sub>/c), not a first-principles derivation from a specified cosmic scalar solution. The pipeline is ready to search for these modulations in Fermilab E989 data upon access.

### A.18 QED Recovery: Gauge Invariance and Ward Identities

The finite-core regulator (Appendix A.15) renders loop integrals finite. Here we examine gauge invariance and the Ward identity under the finite-core modification. The naive Gaussian regulator breaks translational invariance and the Ward identity at O(q<sup>2</sup>r<sub>c</sub><sup>2</sup>). We outline a gauge-covariant formulation that is a candidate for restoring transversality, though exact all-orders proof is not established.

**Gauge invariance from the U(1) bundle.** The bundle action (4a) is invariant under the U(1) gauge transformation &theta; &map; &theta; + &alpha;(x), &omega;<sub>&mu;</sub> &map; &omega;<sub>&mu;</sub> + &part;<sub>&mu;</sub>&alpha;, where &omega;<sub>&mu;</sub> is the orientation-bundle connection. The matter coupling is through the orientation covariant derivative D<sub>&mu;</sub><sup>(orient)</sub> = &part;<sub>&mu;</sub> &minus; iq&omega;<sub>&mu;</sub>, which transforms as D<sub>&mu;</sub><sup>(orient)</sup>&psi; &map; e<sup>iq&alpha;</sup> D<sub>&mu;</sub><sup>(orient)</sup>&psi;. The fermion Lagrangian L = &psi;&#771;(i&gamma;<sup>&mu;</sup>D<sub>&mu;</sub><sup>(orient)</sup> &minus; m)&psi; is therefore gauge invariant by construction. The orientation gauge current J<sup>&mu;</sup><sub>orient</sub> = q&psi;&#771;&gamma;<sup>&mu;</sup>&psi; is conserved: &part;<sub>&mu;</sub>J<sup>&mu;</sup><sub>orient</sub> = 0, by Noether's theorem. *Notation.* The orientation connection &omega;<sub>&mu;</sub> is distinct from the electromagnetic gauge potential *&Atilde;<sub>&mu;</sub>*. Identification of &omega;<sub>&mu;</sub> with *&Atilde;<sub>&mu;</sub>* — or an explicit mixing term such as &epsilon;&Omega;<sub>&mu;&nu;</sub>F<sub>EM</sub><sup>&mu;&nu;</sup> — is an additional model assumption not established by the bundle action alone. Gauge invariance of the orientation bundle does not by itself establish electromagnetic Ward identities; the latter require a separate U(1)<sub>EM</sub> gauge symmetry with potential *&Atilde;<sub>&mu;</sub>*.

**Ward identity with finite-core regulator.** The Ward identity states that the longitudinal part of the photon propagator is unrenormalized: q<sub>&mu;</sub>&Pi;<sup>&mu;&nu;</sup>(q) = 0, where &Pi;<sup>&mu;&nu;</sup> is the vacuum polarization tensor. In standard QED this follows from gauge invariance. With the finite-core regulator R(k) = exp(&minus;k<sup>2</sup>&sigma;<sup>2</sup>), the one-loop vacuum polarization becomes:

&Pi;<sup>&mu;&nu;</sup>(q) = &minus;e<sup>2</sup> &int; d<sup>4</sup>k/(2&pi;)<sup>4</sup> R(k) Tr[&gamma;<sup>&mu;</sup> (k&#771; + m)<sup>&minus;1</sup> &gamma;<sup>&nu;</sup> ((k&#771; + q&#771;) + m)<sup>&minus;1</sup>].
(56)

Contracting with q<sub>&mu;</sub>:

q<sub>&mu;</sub>&Pi;<sup>&mu;&nu;</sup>(q) = &minus;e<sup>2</sup> &int; d<sup>4</sup>k/(2&pi;)<sup>4</sup> R(k) Tr[&gamma;<sup>&nu;</sup> ((k&#771; + q&#771;) + m)<sup>&minus;1</sup> &minus; &gamma;<sup>&nu;</sup> (k&#771; + m)<sup>&minus;1</sup>],
(57)

where we used the identity q&#771; = (k&#771; + q&#771;) &minus; k&#771; and the cyclic property of the trace. The integrand is a total difference: f(k + q) &minus; f(k) with f(k) = R(k) Tr[&gamma;<sup>&nu;</sup>(k&#771; + m)<sup>&minus;1</sup>]. In standard QED (R = 1), this vanishes by translational invariance of the loop integral. With the finite-core regulator R(k) = exp(&minus;k<sup>2</sup>&sigma;<sup>2</sup>), the regulator breaks translational invariance: R(k + q) &ne; R(k). However, the Ward identity is preserved because the regulator depends on the virtual momentum k, not on the external momentum q. The regulated difference is:

q<sub>&mu;</sub>&Pi;<sup>&mu;&nu;</sup>(q) = &minus;e<sup>2</sup> &int; d<sup>4</sup>k/(2&pi;)<sup>4</sup> [R(k) &minus; R(k &minus; q)] Tr[&gamma;<sup>&nu;</sup>(k&#771; + m)<sup>&minus;1</sup>].
(58)

For the Gaussian regulator, R(k) &minus; R(k &minus; q) = exp(&minus;k<sup>2</sup>&sigma;<sup>2</sup>) [1 &minus; exp(2k&middot;q&sigma;<sup>2</sup> &minus; q<sup>2</sup>&sigma;<sup>2</sup>)]. This is not identically zero: the naive Gaussian regulator breaks the Ward identity at O(q<sup>2</sup>&sigma;<sup>2</sup>). The breaking vanishes as &sigma; &rarr; 0 (point-particle limit). At experimentally accessible momenta q << 1/r<sub>c</sub> &sim; m<sub>e</sub>c, the violation is suppressed, but it is not zero.

**Gauge-covariant formulation.** The O(q<sup>2</sup>r<sub>c</sub><sup>2</sup>) breaking arises because the naive regulator R(k) depends on the virtual momentum k rather than the covariant derivative. The gauge-covariant replacement is R(&minus;D<sup>2</sup>), where D<sub>&mu;</sub> = &part;<sub>&mu;</sub> &minus; iq*&Atilde;<sub>&mu;</sub>* is the electromagnetic gauge-covariant derivative. In position space, the covariantly smeared current is:

J<sup>&mu;</sup>(x) = &int; d<sup>4</sup>y   F(x &minus; y)   U(x, y)   j<sup>&mu;</sup><sub>pt</sub>(y)   U<sup>&dagger;</sup>(x, y),
(58a)

where F(x &minus; y) is the normalized form factor (&int; d<sup>4</sup>x F(x) = 1) whose Fourier transform is R(k) = exp(&minus;k<sup>2</sup>&sigma;<sup>2</sup>) with &sigma; = r<sub>c</sub>/&radic;2, and U(x, y) = exp(iq&int;<sub>y</sub><sup>x</sup> *&Atilde;<sub>&mu;</sub>*dz<sup>&mu;</sup>) is the Wilson line ensuring electromagnetic gauge covariance. Under *&Atilde;<sub>&mu;</sub>* &map; *&Atilde;<sub>&mu;</sub>* + &part;<sub>&mu;</sub>&alpha;, the Wilson line transforms as U(x, y) &map; e<sup>iq&alpha;(x)</sup> U(x, y) e<sup>&minus;iq&alpha;(y)</sup>, so J<sup>&mu;</sup>(x) &map; e<sup>iq&alpha;(x)</sup> J<sup>&mu;</sup>(x) e<sup>&minus;iq&alpha;(x)</sup>, which is gauge covariant. The vacuum polarization computed with this covariantly smeared current is:

&Pi;<sup>&mu;&nu;</sup><sub>cov</sub>(q) = &Pi;<sup>&mu;&nu;</sup><sub>naive</sub>(q) + &Pi;<sup>&mu;&nu;</sup><sub>seagull</sub>(q),
(58c)

where &Pi;<sup>&mu;&nu;</sup><sub>naive</sub> is the one-loop integral (56) with R(k), and &Pi;<sup>&mu;&nu;</sup><sub>seagull</sub> is the seagull (contact) contribution from the Wilson-line expansion. Expanding U(x, y) to first order in the gauge field, the Wilson line generates a two-photon vertex (seagull term) at the same spacetime point as the form-factor insertion:

&Pi;<sup>&mu;&nu;</sup><sub>seagull</sub>(q) = &minus;e<sup>2</sup> &int; d<sup>4</sup>k/(2&pi;)<sup>4</sup> [R(k) &minus; R(k &minus; q)] Tr[&gamma;<sup>&mu;</sup>(k&#771; + m)<sup>&minus;1</sup>&gamma;<sup>&nu;</sup>(k&#771; + m)<sup>&minus;1</sup>] + &delta;<sup>&mu;&nu;</sup> e<sup>2</sup> &int; d<sup>4</sup>k/(2&pi;)<sup>4</sup> R(k) Tr[(k&#771; + m)<sup>&minus;1</sup>].
(58d)

The first term is asserted to cancel the non-transverse part of Eq. 58, while the second term provides the transverse completion. The candidate Ward-Takahashi identity is:

q<sub>&mu;</sub> &Pi;<sup>&mu;&nu;</sup><sub>cov</sub>(q) = 0,
(58b)

This holds to all orders in q<sup>2</sup>&sigma;<sup>2</sup> because the Wilson-line construction ensures gauge covariance at each vertex. However, this result is asserted rather than demonstrated from a complete nonlocal action: (i) the Wilson-line path is not specified; (ii) the full tower of multi-photon vertices is not derived; (iii) the displayed seagull term is asserted to cancel the violation but is not shown to follow from the same action; (iv) the propagator, vertex, and regulator definitions are not shown to follow from one common action; (v) the orientation U(1) is still not identified consistently with electromagnetic U(1). The Wilson-line construction therefore outlines a possible gauge-covariant completion, not an established exact all-orders Ward–Takahashi identity. The modified fermion–photon vertex in the covariant formulation is:

&Gamma;<sup>&mu;</sup><sub>cov</sub>(p, p') = &gamma;<sup>&mu;</sup> R((p + p')/2) + (q&sigma;<sup>2</sup>) (p + p')<sup>&mu;</sup> R'((p + p')/2) + O(q<sup>2</sup>&sigma;<sup>4</sup>),
(58e)

where R'(k) = dR/dk<sup>2</sup>. The first term is the standard vertex with the form-factor evaluated at the average momentum. The second term is the Wilson-line correction that is a candidate for restoring the Ward–Takahashi identity q<sub>&mu;</sub>&Gamma;<sup>&mu;</sup><sub>cov</sub> = S<sup>&minus;1</sup>(p') &minus; S<sup>&minus;1</sup>(p), where S(p) = R(p)(p&#771; + m)<sup>&minus;1</sup> is the regulated propagator. At probe scales q &sim; 1/r<sub>c</sub>, the probe wavelength resolves the interior of the temporal shear defect. The generalized identity (58b) is a candidate result; the apparent Ward-identity "leakage" in the naive formulation (Eq. 58) is the finite-resolution signature of resolving a fluid vortex core rather than a Dirac point, and the seagull term (58d) is asserted to cancel it exactly, pending a complete nonlocal action derivation.

**Lamb shift estimate.** The finite-core regulator modifies the vacuum polarization and self-energy contributions to the Lamb shift. An order-of-magnitude estimate gives a TEP correction of &sim;53 kHz, which may already be experimentally constrained. A regulated bound-state QED calculation is required before this can be used as a consistency check. Detailed diagnostics are documented in the repository audit (`scripts/audit/`).

### A.19 Braid-Group Holonomy and Multi-Defect Hamiltonian

The single-defect spinorial lift (Appendix A.7) and the antisymmetrized two-defect construction (Appendix A.13) establish single-defect spinorial transformation and a consistent two-defect antisymmetric state. Here we develop the braid-group holonomy for adiabatic vortex exchange in 2+1D and the multi-defect interaction Hamiltonian in 3+1D.

**2+1D braid-group holonomy.** In 2+1 dimensions, the configuration space of N identical topological charges has fundamental group equal to the braid group B<sub>N</sub>, not the symmetric group S<sub>N</sub>. An adiabatic exchange of two vortices in the temporal orientation bundle corresponds to a braid &sigma;<sub>i</sub> &isin; B<sub>N</sub>. The holonomy of the U(1) connection around the braid path gives a phase:

U(&sigma;<sub>i</sub>) = exp(i&pi;n<sub>i</sub>) = (&minus;1)<sup>n<sub>i</sub></sup>,
(59)

where n<sub>i</sub> = &pm;1 is the winding number of the exchanged vortex. For fermions (n = 1), U(&sigma;<sub>i</sub>) = &minus;1, reproducing the exchange antisymmetry. The braid-group representation is:

&rho;: B<sub>N</sub> &rarr; U(1),   &rho;(&sigma;<sub>i</sub>) = &minus;1,   &rho;(&sigma;<sub>i</sub>&sigma;<sub>j</sub>) = &rho;(&sigma;<sub>j</sub>&sigma;<sub>i</sub>) for |i &minus; j| &ge; 2,   &rho;(&sigma;<sub>i</sub>&sigma;<sub>i+1</sub>&sigma;<sub>i</sub>) = &rho;(&sigma;<sub>i+1</sub>&sigma;<sub>i</sub>&sigma;<sub>i+1</sub>).
(60)

The last relation (Yang–Baxter equation) is satisfied trivially because &rho; is a one-dimensional representation. In 2+1D, this is the abelian anyon representation with statistical angle &theta; = &pi; (fermionic). The orientation bundle therefore admits the fermionic one-dimensional braid representation in 2+1 dimensions.

**3+1D multi-defect Hamiltonian.** In 3+1 dimensions, the braid group reduces to the symmetric group S<sub>N</sub> because braids can be unwound. The multi-defect Hamiltonian for N topological charges at positions x<sub>i</sub> is constructed from the single-defect Lagrangian by summing the scalar field contributions and including the disformal cross-coupling:

H<sub>N</sub> = &sum;<sub>i=1</sub><sup>N</sup> [ &radic;(p<sub>i</sub><sup>2</sup>c<sup>2</sup> + m<sub>i</sub><sup>2</sup>c<sup>4</sup>A<sup>2</sup>(&phi;(x<sub>i</sub>))) ] + &sum;<sub>i<j</sub> V<sub>int</sub>(|x<sub>i</sub> &minus; x<sub>j</sub>|),
(61)

where the interaction potential V<sub>int</sub> arises from the overlap of the scalar fields &phi;<sub>i</sub> and &phi;<sub>j</sub> through the disformal coupling:

V<sub>int</sub>(r<sub>ij</sub>) = &minus;(B&beta;<sub>A</sub><sup>2</sup> / M<sub>Pl</sub><sup>2</sup>) (&phi;<sub>0</sub><sup>2</sup> / 4&pi;) exp(&minus;m<sub>&phi;</sub> r<sub>ij</sub>) / r<sub>ij</sub>,
(62)

with &phi;<sub>0</sub> the single-defect scalar amplitude and r<sub>ij</sub> = |x<sub>i</sub> &minus; x<sub>j</sub>|. This is a Yukawa interaction mediated by the scalar field, with range &lambda;<sub>scr</sub> = &hbar;/(m<sub>&phi;</sub>c) = &radic;2 r<sub>c</sub>. This Yukawa interaction V<sub>int</sub> represents the static potential limit of the dynamic kinematic routing derived in TEP-KIN (Paper 25). The disformal light-cone tilt B(&phi;), which governs the geometry of measurement and interaction routing between distinct causal frames, mechanically manifests as this spatial potential when integrated over the non-perturbative overlap of the topological cores. Note that V<sub>int</sub> represents the pure scalar-disformal interaction native to the TEP geometry. Because the U(1) bundle action (Eq. 4a) preserves the standard Maxwell kinetic term, this Hamiltonian supplements, rather than replaces, the standard electromagnetic Coulomb interactions between charged defects. The full interaction Hamiltonian is therefore H<sub>int</sub> = H<sub>QED</sub> + V<sub>int</sub>. The exchange antisymmetry is enforced by standard fermionic quantum mechanics: the multi-defect wavefunction is the Slater determinant of single-defect spinors, which changes sign under any pair exchange by construction (Appendix A.13). The braid-group holonomy in 3+1D reduces to the symmetric group representation, and the multi-defect Hamiltonian (62) with Yukawa interactions and Slater-determinant statistics is the effective multi-defect Hamiltonian consistent with standard fermionic statistics.


## 6. Data Availability & Reproducibility


This work follows open-science practices. All theoretical derivations and numerical results
are fully reproducible using the documented code.



### Repository and Code


GitHub Repository: github.com/matthewsmawfield/TEP-SPIN

Zenodo DOI (v0.1 release): 10.5281/zenodo.20572705



The repository contains the analytical derivations and numerical verification scripts
for the TEP spin-coupling framework, screening model, and empirical constraints.
The Zenodo DOI provides an immutable, version-locked snapshot of the exact analytical
pipelines used for this draft.



### Repository Structure


TEP-SPIN/
├── data/
│   ├── gm2/                  # Muon g-2 summary data
│   └── jlab_prad/            # JLab PRad + A1 Collaboration data
├── scripts/
│   ├── steps/                # Analysis pipeline steps
│   └── utils/                # Shared utilities
├── core/                     # TEP shared constants and parameters
├── results/                  # Pipeline outputs and figures
├── site/
│   ├── components/           # Manuscript HTML sections
│   └── public/               # Built site assets
├── requirements.txt
├── CITATION.bib
└── README.md



### Software Environment


Key packages: NumPy, SciPy, SymPy, Matplotlib.
The scripts have been tested on Python 3.10+.



### License


All code and manuscripts are released under CC-BY-4.0.