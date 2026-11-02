# Stochastic SEICR FMD Model

This repository contains the MATLAB implementation for the Stochastic SEICR FMD Model under Periodic Environmental Forcing.

## File Structure
* `seicrParams.m`: Defines the biological parameters and noise scenarios.
* `seicrDrift.m`: Defines the deterministic drift vector \(f(x,t)\).
* `simulateSEICR.m`: Core simulation engine (supports Euler-Maruyama, Milstein, and NSFD methods).
* `trajectory.m`: Plots a single long-term trajectory.
* `highInfectedNoise.m`: Runs multiple simulations under high noise and calculates extinction probability.
* `run_stochasticThreshold.m`: Calculates the period-averaged noise-corrected threshold (\(\mathcal{R}_s\)).

## How to Run
1. Open MATLAB and navigate to the folder containing these files.
2. To run the main trajectory plot, type:
   ```matlab
   trajectory
   run_stochasticThreshold
   highInfectedNoise
   
#### Step 4: Upload to GitHub
1. On your GitHub repository page, click **Add file** -> **Upload files**.
2. Drag and drop all the `.m` files you created in Step 2.
3. Scroll down and click **Commit changes**.

*(Alternatively, if you know basic Git, you can run `git init`, `git add .`, `git commit -m "Initial commit"`, and `git push` from your terminal).*

### Important Note on the Code:
I noticed a typo in the provided text in the `simulateSEICR.m` section. The text has `\end{listing}` instead of `\end{lstlisting}`. Also, when you copy the MATLAB code into `.m` files, **do not include the `%saved as filename.m` lines**; those are just comments in your LaTeX document. Start the actual `.m` files from the `function ...` line.
