# Data and Code for *Period Homeostasis Near Hopf Bifurcation*

This repository contains the Python code and data supporting the manuscript *Period Homeostasis Near Hopf Bifurcation* by
Steve Manns, Janet Best, and Martin Golubitsky.

The simulations use the dimensionless Kim–Forger single-negative-feedback
(SNF) model of the mammalian circadian clock. For a fixed value of the bifurcation parameter 
`K_d`, the code determines the interval of values of the input parameter `A` in which
a periodic solution exists, integrates the model, estimates the period
from successive maxima of `M(t)`, and locates critical points of the period
input-output function. These critical points are infinitesimal period
homeostasis (IPH) points.

## Installation

Python 3.10 or later is required. From the repository root, create an isolated
environment and install the package:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

## Reproducing the Figures

The notebooks use paths relative to the `notebooks/` directory. Launch Jupyter
from that directory and run the relevant notebook from top to bottom:

```bash
cd notebooks
jupyter lab
```

| Manuscript figure | Notebook | Output |
| --- | --- | --- |
| Figure 1 | `figure1_normal_form.ipynb` | `figure1a.png`, `figure1b.png` |
| Figure 2 | `figure2_parameter_space.ipynb` | `figure2_parameter_space.png` |
| Figure 3 | `figure3_hh_graphs.ipynb` | `figure3a.png`, `figure3b.png` |
| Figure 4 | `figure4_period_input_output.ipynb` | `figure4a.png`–`figure4d.png` |

Generated figures are written to `outputs/`.

Figure 2 has two reproducibility paths:

- **Plot from included data:** skip the grid-simulation and branch-assignment
  cells, then run the cells that load `data/iph_curves.csv` and make the plot.
- **Regenerate the data:** run the entire notebook. This performs the full IPH
  scan, writes `data/iph_points_all.csv` and `data/iph_curves.csv`, and then
  makes the plot. The scan took approximately 30–45 minutes on the machine used
  to prepare the manuscript; runtime will depend on hardware.

Figure 3 reads the included `data/iph_curves.csv`. Figure 4 recomputes its four
fixed-`K_d` period curves when run and takes approximately 30 seconds per panel
on the machine used to prepare the manuscript.

## Numerical Method and Data

The model is integrated with `scipy.integrate.solve_ivp` using LSODA. After
discarding the transient, maxima of `M(t)` are detected as downward zero
crossings of `dM/dt`; the period is the time between successive maxima. Period
derivatives with respect to `A` are estimated by finite differences, and IPH
points are interpolated where those derivatives change sign.

The full IPH scan uses 54 uniformly spaced `K_d` values from `1e-6` through
`2.66e-4`, separated by `5e-6`, and samples 200 values of `A` within each
periodic interval. Exact solver, interval, scan, and branch-assignment settings
are recorded in `data/simulation_settings.yml`. This YAML file is a
human-readable record; however, the executable defaults are defined in
`src/kim_forger_snf/period.py`, and the scan is defined in the Figure 2
notebook.

| File | Contents |
| --- | --- |
| `data/iph_points_all.csv` | All IPH points detected by the fixed-`K_d` scan before branch assignment |
| `data/iph_curves.csv` | The same IPH points with nearest-neighbor branch labels used in Figures 2 and 3 |
| `data/simulation_settings.yml` | Numerical settings used to generate the included data and Figure 4 slices |

All data in `data/` are generated computationally by this repository; no
experimental or third-party datasets are required.

## Repository Structure

```text
src/kim_forger_snf/
  model.py       Dimensionless Kim–Forger SNF model
  hopf.py        Analytical Hopf curve and periodic-interval utilities
  period.py      Numerical integration and period estimation
  iph.py         IPH detection from period curves

notebooks/       Figure-reproduction workflows
data/            Included derived data and simulation settings
outputs/         Generated manuscript figures
```

## Citation and Archival Availability

Please cite the associated manuscript:

> Steve Manns, Janet Best, and Martin Golubitsky, *Period Homeostasis Near Hopf
> Bifurcation*, manuscript in preparation.

A version-specific archival DOI for this repository will be added here after
the public release is deposited.
