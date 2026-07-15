# Classical & Quantum Double Pendulum (Fortran)

## Files

| File | What it is |
|---|---|
| `classical_double_pendulum.f90` | Classical double pendulum, Lagrangian mechanics + RK4 |
| `quantum_double_pendulum.f90` | Quantum double pendulum, full 2D-grid time-dependent Schrödinger equation |
| `plot_results.py` | Reads the `.dat` output of both programs and makes plots |
| `sample_*.dat`, `sample_*.png` | Reference output/plots (see note below) |

## Building and running

```bash
gfortran -O2 -o classical classical_double_pendulum.f90
./classical
# -> writes classical_trajectory.dat, prints a summary table

gfortran -O2 -o quantum quantum_double_pendulum.f90
./quantum
# -> writes quantum_observables.dat and quantum_density_snapshots.dat, prints a summary table

pip install numpy matplotlib   # if you don't have them
python3 plot_results.py
# -> writes overview.png and quantum_snapshots.png
```

**Important honest disclosure:** the sandbox I wrote this in doesn't have `gfortran`
installed and has no network access to install it, so I could not compile these
`.f90` files myself. What I *did* do instead:

- Independently re-implemented both algorithms in Python (vectorized numpy) and
  validated them numerically — see "Validation" below.
- Ran `numpy.f2py`'s Fortran parser over both files to catch gross structural
  errors (unmatched blocks, bad subroutine boundaries, etc. — it parsed cleanly).
- Manually cross-checked every line of the Fortran against the validated Python,
  term by term.

I'm confident in the physics and the algorithm, but since I never ran an actual
Fortran compiler on these files, there's a small chance of a stray syntax slip
(a missing comma, a mismatched parenthesis) that only `gfortran` itself would
catch. Compile it first thing and let me know if it doesn't build — I'll fix
it immediately.

The `sample_*` files were generated from the validated Python implementation
(exactly the same algorithm, translated line-for-line into the Fortran), **not**
from running the `.f90` files. They're there so you can see what the real output
looks like and test `plot_results.py` before you've even compiled anything.

## The physics

### Classical (`classical_double_pendulum.f90`)

Standard double pendulum: two point masses `m1`, `m2` on massless rigid rods
`L1`, `L2`, angles measured from vertical. The equations of motion are derived
from the 2×2 mass (inertia) matrix `M(theta)` and solved in closed form at each
RK4 stage — no linear-algebra library needed. Default parameters: `m1=m2=1 kg`,
`L1=L2=1 m`, `g=9.81 m/s^2`, released from rest at `theta1=120°, theta2=-10°`
(a well-known chaotic regime). Energy is conserved to ~1e-9 J over the whole
20-second run at `dt=0.0005 s` — that's your check that the RK4 step size is
small enough.

### Quantum (`quantum_double_pendulum.f90`)

This is the more unusual piece, so it's worth explaining the modeling choice.

The classical kinetic energy has a **position-dependent** mass matrix
(`M12 = m2 L1 L2 cos(theta1-theta2)` depends on the angles), so this isn't a
free particle on a torus — it's a particle with an anisotropic,
configuration-dependent effective mass. Quantizing `p -> -i*hbar d/dtheta` in
a Hamiltonian where the mass matrix doesn't commute with the momentum operators
requires picking an **operator ordering**. This program uses the standard
symmetric ("BenDaniel-Duke"-type) choice

```
T_hat psi = -(hbar^2/2) * sum_ab  d/dtheta_a ( G_ab(theta) d/dtheta_b psi ),   G = M^-1
```

which is Hermitian by construction, and discretizes it with a *conservative*
finite-difference scheme on a periodic 80×80 `(theta1,theta2)` grid. The
cross-derivative (`G12`) term needed a specific symmetrized discretization to
stay exactly Hermitian on the grid — I found this by trial and error during
development (a naive central-difference cross term broke Hermiticity by ~1%,
which was enough to make the norm blow up over a few hundred steps; the fix
was verified to give **exact** (0.0 to machine precision) Hermiticity by
building the operator as an explicit matrix and checking it against its own
conjugate transpose).

Time-stepping is explicit RK4 applied directly to
`d(psi)/dt = -(i/hbar) H psi`. RK4 isn't exactly unitary, but at the chosen
step size the norm drifts by less than 1e-8 and the energy expectation by
less than 1e-7 over the full run — both are reported in the program's summary
table so you can see the numerical quality for yourself.

**Units:** a real pendulum's quantum effects are far too small to see (its de
Broglie wavelength is ~1e-34 m), so this uses natural units
`hbar = m1 = m2 = L1 = L2 = g = 1`. It's a toy model of "a double pendulum at
quantum scale," not a quantitative prediction for a real one-metre pendulum —
treat the two programs' time axes as different unit systems (the comparison
plot in `plot_results.py` is deliberately labeled this way).

The initial state is a Gaussian wavepacket centered at the *same* angles as
the classical run (`120°, -10°`), released with zero mean momentum, so you can
watch where the quantum expectation values `<theta1>(t)`, `<theta2>(t)`
initially track the classical trajectory and then diverge from it as the
nonlinear potential spreads and deforms the wavepacket — a nice concrete look
at where the classical/quantum correspondence breaks down.

## Validation performed

- Classical RK4: verified energy conservation to ~1e-9 J over 40,000 steps.
- Quantum Hamiltonian: verified exact (machine-precision) Hermiticity by
  building the discrete operator as a dense matrix and comparing it to its
  conjugate transpose, for the precomputed-array structure used in the
  Fortran.
- Quantum time evolution: verified norm and energy expectation conservation
  to <1e-8 / <1e-7 respectively over the full 10,000-step run, at three
  different step sizes (confirms it's genuinely converged, not accidentally
  cancelling errors).
- Both `.f90` files parse cleanly with `numpy.f2py`'s Fortran front end.

## Tweaking parameters

Both files have a clearly marked block near the top of the `program` (masses,
lengths, initial angles, `dt`, run length, grid size) — edit those and
recompile. Grid size and `dt` are linked: halving the grid spacing roughly
quadruples the fastest frequency in the Hamiltonian, so you'll want a smaller
`dt` too if you make the quantum grid finer.
