# Stochastic SEICR FMD Model

This repository contains the MATLAB implementation for the Stochastic SEICR FMD Model under Periodic Environmental Forcing.

## File Structure
* `seicrParams.m`: Defines the biological parameters and noise scenarios.
* `seicrDrift.m`: Defines the deterministic drift vector \(f(x,t)\).
* `simulateSEICR.m`: Core simulation engine (supports Euler-Maruyama, Milstein, and NSFD methods).
* `trajectory.m`: Plots a single long-term trajectory.
* `highInfectedNoise.m`: Runs multiple simulations under high noise and calculates extinction probability.
* `run_stochasticThreshold.m`: Calculates the period-averaged noise-corrected threshold (\(\mathcal{R}_s\)).

## How to Cite

If you use this code in your research, please cite it as:

Palivamwe M Ndeevelo. (2026). *221092/SEICR-FMD-Model: v1.0.0 - Initial Release of the SEICR-FMD Model* (Version V01) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23889327
