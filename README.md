# PermaCycle

[![DOI](https://zenodo.org/badge/DOI/XXXXX.svg)](https://doi.org/XXXXX)

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

## Requirements

-   MATLAB (R20XX or later recommended)

------------------------------------------------------------------------

## Installation

``` bash
git clone https://github.com/KehuaYou/PermaCycle.git
```

Open MATLAB and run:

``` matlab
cd('path_to_PermaCycle')
addpath(genpath(pwd))
```

------------------------------------------------------------------------

## Usage

Run the main simulation:

``` matlab
main.m
```


------------------------------------------------------------------------

## Inputs

Key model inputs are defined in initialization.m and include:

-   Current water depth of the modeled location (water_depth_simul)
-   Geothermal heat flux (qt)
-   Water freezing curve
	Slope (Ax)
	Residual unfrozen water saturation (Swr_freeze)
-   Thermal conductivity of solid grains (lambda_s)
-   E-folding depth for porosity (e_fold)
------------------------------------------------------------------------

## Outputs

The model outputs:

-   Temperature (°C)
-   Pore pressure (Pa)
-   Salinity (wt.%)
-   Ice saturation (-)

------------------------------------------------------------------------

## Citation

> You, K. (2026). *PermaCycle v1.0.0* \[Software\]. Zenodo.
> https://doi.org/XXXXX

------------------------------------------------------------------------

## License

MIT License

------------------------------------------------------------------------

## Contact

Kehua You\
University of Texas at Austin\
khyouml@gmail.com
