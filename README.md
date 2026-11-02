# Stochastic SEICR FMD Model

This repository contains the MATLAB implementation for the Stochastic SEICR FMD Model under Periodic Environmental Forcing.

## File Structure
* `seicrParams.m`: Defines the biological parameters and noise scenarios.
* `seicrDrift.m`: Defines the deterministic drift vector \(f(x,t)\).
* `simulateSEICR.m`: Core simulation engine (supports Euler-Maruyama, Milstein, and NSFD methods).
* `trajectory.m`: Plots a single long-term trajectory.
* `highInfectedNoise.m`: Runs multiple simulations under high noise and calculates extinction probability.
* `run_stochasticThreshold.m`: Calculates the period-averaged noise-corrected threshold (\(\mathcal{R}_s\)).

