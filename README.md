# Linearized Bicycle Model for Lateral Stability Analysis

Implementation of the benchmark linearized bicycle model from Meijaard, Papadopoulos, Ruina and Schwab (2007).

---

## Background

The stability of a single-track vehicle — bicycle or motorcycle — depends on the coupling between lean angle and steer angle through inertia, gyroscopic, gravitational, and centrifugal mechanisms. Analysing this coupling in its simplest form requires a linearized model valid for small perturbations about straight-ahead upright riding at constant speed.

Meijaard et al. (2007) derived and verified such a model using five independent methods (Newtonian, Lagrangian, two computer algebra systems, and linearization of the full nonlinear equations), resolving decades of conflicting results in the literature. They published explicit benchmark parameter values and benchmark eigenvalues that serve as a precise verification standard. This implementation reproduces those benchmarks exactly and extends the analysis with eigenvalue sweeps, stability boundary detection, root locus plots, and time-domain simulation.

---

## Equations of Motion

The linearized lateral dynamics are governed by:

```
M·q̈ + v·C₁·q̇ + (g·K₀ + v²·K₂)·q = f
```

where the generalised coordinates are:

```
q = [φ, δ]ᵀ
```

- `φ` — lean angle (roll of rear frame from vertical), rad  
- `δ` — steer angle (front assembly relative to rear frame), rad  
- `v` — forward speed, m/s (parameter, not a state)  
- `f` — external generalised force vector [lean torque; steer torque]

The four 2×2 coefficient matrices are:

| Matrix | Depends on | Physical origin | Symmetric |
|--------|-----------|-----------------|-----------|
| **M** | Parameters only | Total system inertia | Yes |
| **v·C₁** | Speed v | Gyroscopic torques + constraint forces | No |
| **g·K₀** | Gravity g | Gravitational potential energy | Yes |
| **v²·K₂** | Speed² | Centrifugal + trail effect | No |

Converting to first-order state-space form with `x = [q; q̇]ᵀ`:

```
ẋ = A·x

A = [  0₂ₓ₂      |      I₂ₓ₂      ]
    [ −M⁻¹·K_eff | −M⁻¹·(v·C₁)  ]

K_eff = g·K₀ + v²·K₂
```

The four eigenvalues of the 4×4 matrix **A** determine stability at each speed.

---

## The Three Stability Modes

The eigenvalue analysis reveals three physically distinct modes:

**Castering mode** — always a large negative real eigenvalue. The front wheel snaps back to centre when deflected. Stabilised by trail at all speeds.

**Weave mode** — a complex conjugate pair. Oscillatory coupling between lean, steer, and yaw. Unstable at low speed, becomes stable above the weave speed `v_weave`. Frequency increases approximately linearly with speed.

**Capsize mode** — a small real eigenvalue. Slow non-oscillatory lean divergence. Unstable at low speed, stabilised at intermediate speed by gyroscopic action, becomes mildly unstable again above the capsize speed `v_capsize`.

The **self-stable corridor** is the speed range where all four eigenvalues have negative real parts:

```
v_weave  <  v  <  v_capsize
```

For the Meijaard benchmark parameters: approximately **4.3 m/s to 6.0 m/s**.

---

## Benchmark Parameters

All parameters are taken directly from Table 1 of Meijaard et al. (2007). The vehicle represents a bicycle with a heavy rigid rider.

### Geometry

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Wheelbase | w | 1.02 | m |
| Trail | c | 0.08 | m |
| Steer axis tilt (from vertical) | λ | π/10 = 18° | rad |

### Rear wheel (R)

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Radius | r_R | 0.30 | m |
| Mass | m_R | 2.0 | kg |
| Spin inertia | I_Rxx | 0.0603 | kg·m² |
| Roll inertia | I_Ryy | 0.12 | kg·m² |

### Rear frame + rider (B)

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| CoM x-position | x_B | 0.30 | m |
| CoM z-position | z_B | −0.90 | m |
| Mass | m_B | 85.0 | kg |
| Roll inertia | I_Bxx | 9.2 | kg·m² |
| Product of inertia | I_Bxz | 2.4 | kg·m² |

### Front fork + handlebar (H)

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| CoM x-position | x_H | 0.90 | m |
| CoM z-position | z_H | −0.70 | m |
| Mass | m_H | 4.0 | kg |
| Roll inertia | I_Hxx | 0.05892 | kg·m² |
| Product of inertia | I_Hxz | −0.00756 | kg·m² |

### Front wheel (F)

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Radius | r_F | 0.35 | m |
| Mass | m_F | 3.0 | kg |
| Spin inertia | I_Fxx | 0.1405 | kg·m² |
| Roll inertia | I_Fyy | 0.28 | kg·m² |

---

## Benchmark Matrices

Reproduced from Table 2 of Meijaard et al. (2007):

```
M  = [ 80.81722,          2.31941332208709 ]
     [  2.31941332208709,  0.29784188199686 ]

C₁ = [  0,               33.86641391492494 ]
     [ -0.85035641456978,  1.68540397397560 ]

K₀ = [-80.95,            -2.59951685249872 ]
     [ -2.59951685249872, -0.80329488458618 ]

K₂ = [  0,               76.59734589573222 ]
     [  0,                2.65431523794604  ]
```

---


## Code Structure

The implementation is a single self-contained MATLAB script `motorcycle_stability_meijaard.m` with the following sections:

```
Section 1  — Benchmark parameters (Table 1)
Section 2  — Derived quantities: CoM, front assembly inertia,
             gyroscopic coefficients, steer-axis projections
Section 3  — System matrices M, C₁, K₀, K₂ (Table 2 values)
Section 4  — Eigenvalue sweep: v = 0 to 10 m/s
Section 5  — Stability boundary detection (weave speed, capsize speed)
Section 6  — Four-panel stability figure
Section 7  — Time-domain simulation at v = 4.5 m/s (ODE45)
Section 8  — Printed eigenvalue summary table
Local fn   — compute_eigs(M, C1, K0, K2, v, g)
```

---

## Output Figures

**Figure 1 — Linearized Bicycle Stability (2×2 subplot)**

- Top-left: Real parts of all four eigenvalues vs speed. Stable corridor shaded green.
- Top-right: Imaginary parts (oscillation frequencies) vs speed.
- Bottom-left: Root locus in the complex plane, coloured by speed.
- Bottom-right: Stability envelope — maximum real part vs speed, red/green shading.

**Figure 2 — Time-Domain Simulation**

- Top: Lean angle φ and steer angle δ vs time (initial condition: 3° lean).
- Bottom: Angular rates φ̇ and δ̇ vs time.
- Simulated at v = 4.5 m/s (inside the stable corridor) — perturbations decay to zero.

---

## Console Output

The script prints three blocks to the command window:

```
=== Benchmark Matrices (Meijaard et al. 2007) ===
M  = ...   C1 = ...   K0 = ...   K2 = ...

=== Stability Boundaries ===
Weave-stable onset  : v_weave  = 4.30 m/s
Capsize instability : v_capsize= 6.00 m/s

=== Speed-Dependent Eigenvalue Summary ===
v [m/s]   eig1           eig2           eig3           eig4
------------------------------------------------------------------------
0.00      ...
...
10.00     ...
```

---

## References

Meijaard, J. P., Papadopoulos, J. M., Ruina, A. and Schwab, A. L. (2007). Linearized dynamics equations for the balance and steer of a bicycle: a benchmark and review. *Proceedings of the Royal Society A*, 463(2084), 1955–1982. https://doi.org/10.1098/rspa.2007.1857

---


