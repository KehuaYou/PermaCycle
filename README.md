# PermaCycle

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20043352.svg)]([https://doi.org/XXXXX](https://doi.org/10.5281/zenodo.20043352))

**PermaCycle** is a one-dimensional MATLAB-based numerical model for
simulating the evolution of subsea permafrost dynamics and methane hydrate
stability zone along Arctic continental shelves over glacial--interglacial
cycles.

------------------------------------------------------------------------

## Overview

PermaCycle simulates coupled **hydrological, thermal, and geochemical
processes** controlling permafrost formation and degradation over long
geological timescales (up to 400 kyr). The model resolves vertical
profiles from the ground surface or seafloor to depths of up to 1,000 m.

------------------------------------------------------------------------

## Key Features

-   Coupled **fluid flow, heat transport, and salt transport**
-   Simulation of **ice formation and melting**
-   Prediction of:
    -   Temperature
    -   Pore pressure
    -   Salinity
    -   Ice saturation
-   Methane hydrate stability zone
-   Time-dependent boundary conditions driven by sea-level change
-   Fully implemented in **MATLAB**

------------------------------------------------------------------------

## Model Description

The model includes:

-   Heat transport via **conduction and advection**
-   Salt transport via **diffusion and advection**
-   Phase change between water and ice using a **freezing curve**
-   Dynamic sediment properties:
    -   Permeability changes with ice saturation
    -   Bulk thermal conductivity variations due to ice formation and melting
-   Methane hydrate stability derived from:
    -   Temperature
    -   Pore pressure
    -   Salinity

A detailed description of the mathematical model and numerical methods can be found in:

You, K. (2024). Biodegradation of Ancient Organic Carbon Fuels Seabed Methane Emission at the Arctic Continental Shelves, Global Biogeochemical Cycles, 38(2), e2023GB007999. https://doi.org/10.1029/2023GB007999

------------------------------------------------------------------------

## System requirements

- **Software:** MATLAB **R2025b**. No additional MATLAB toolboxes are required
  (base MATLAB only; the model uses only built-in functions such as `load` and
  `readmatrix`).
- **Versions tested on:** MATLAB R2025b on **Windows 11**. The code uses no
  operating-system-specific calls and is expected to run on any platform
  supported by MATLAB (Windows, Linux, macOS). <!-- Kehua: add Linux/macOS here if you have tested them -->
- **Non-standard hardware:** none required. PermaCycle runs on a standard
  desktop or laptop. (The pan-Arctic study ran many independent single-column
  simulations in parallel on a compute cluster, but no special hardware is
  needed to run the model.)

------------------------------------------------------------------------

## Installation

```bash
git clone https://github.com/KehuaYou/PermaCycle.git
```

Then, in MATLAB, change into the `src` directory and add the code to the path:

```matlab
cd path/to/PermaCycle/src
addpath(genpath(pwd))
```

**Typical install time:** under 5 minutes (download plus adding to the path; no
compilation required).

------------------------------------------------------------------------

## Instructions for use

1. Edit the inputs in `src/Initialization.m`.
Key model inputs are defined in initialization.m and include:

-   Current water depth of the modeled location (water_depth_simul)
-   Geothermal heat flux (qt)
-   Water freezing curve
	Slope (Ax)
	Residual unfrozen water saturation (Swr_freeze)
-   Thermal conductivity of solid grains (lambda_s)
-   E-folding depth for porosity (e_fold)

2. Provide the site's surface-temperature history (`surface_T.csv`) and relative
sea-level curve (`sealevel.xlsx`), then run `Main_loop` from the `src/`
directory. Results are saved as `.mat` files at the specified save interval.


------------------------------------------------------------------------
## Demo

A worked single-site example is provided in `example/` (present-day water
depth 19 m; 50 vertical cells; geothermal heat flux 0.07 W m⁻²). Run from the
`src/` directory so the input files `surface_T.csv`, `sealevel.xlsx`, and
`methane_hydrate_phase_boundary.mat` are found.

### Quick demo (~1 minute)

To verify the model runs, simulate a short interval. Set the simulation end time
to 5 kyr by editing the stop condition near the top of the main loop in
`src/Main_loop.m` (change `>401*(86400*365*1e3)` to `>5*(86400*365*1e3)`), then run:

```matlab
Main_loop
```

**Expected output:** predicted depth profiles of temperature, pore pressure,
salinity, and ice saturation for the first 5 kyr, saved as `.mat` files.
**Expected run time:** approximately **1 minute** on a normal desktop computer.

### Full example (~1.5 hours, reproduces the paper figures)

With the default settings, `Main_loop` runs the full 400-kyr simulation.

**Expected output:** the full evolution, matching the reference figures in
`example/Outputs/` (`pressure_temperature_salinity_icesaturation.fig` and
`ice_methane_hydrate_stability_zone.fig`). Run `plot_results` to regenerate the
figures from the saved output.
**Expected run time:** approximately **86 minutes** for the single-column
400-kyr run on a normal desktop computer.

------------------------------------------------------------------------

## Outputs

The model outputs:

-   Temperature (°C)
-   Pore pressure (Pa)
-   Salinity (wt.%)
-   Ice saturation (-)

------------------------------------------------------------------------

## Citation

> You, K. (2026). *PermaCycle v1.0.1: A Coupled Hydro–Thermal–Geochemical Model for Permafrost Dynamics* \[Software\]. Zenodo.
> https://doi.org/10.5281/zenodo.20043352

------------------------------------------------------------------------

## License

MIT License

------------------------------------------------------------------------

## Contact

Kehua You\
University of Texas at Austin\
khyouml@gmail.com
