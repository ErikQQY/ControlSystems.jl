# ControlSystems.jl

[![Build Status](https://github.com/JuliaControl/ControlSystems.jl/workflows/CI/badge.svg)](https://github.com/JuliaControl/ControlSystems.jl/actions?query=workflow%3ACI)
[![Documentation Status](https://github.com/JuliaControl/ControlSystems.jl/workflows/Docs/badge.svg)](https://github.com/JuliaControl/ControlSystems.jl/actions?query=workflow%3ADocs)

[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/C/ControlSystems.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![codecov](https://codecov.io/gh/JuliaControl/ControlSystems.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaControl/ControlSystems.jl)

[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliacontrol.github.io/ControlSystems.jl/latest)

A control systems design toolbox for Julia.

## Installation

To install, in the Julia REPL:

```julia
using Pkg; Pkg.add("ControlSystems")
```

## News

### 2021-01
- *Breaking*: Support for julia versions older than 1.3 is dropped
- *Breaking*: `c2d(::StateSpace)` now returns only the system, not the `x0map`. See `c2d_x0map` for the old functionality.
- System order can be specified in `baltrunc`.
- New discretization methods in `c2d`. We now support `:zoh,:foh,:fwdeuler,:tustin`
- Symbolic computation utilities in [SymbolicControlSystems.jl](https://github.com/JuliaControl/SymbolicControlSystems.jl)

More details under [releases](https://github.com/JuliaControl/ControlSystems.jl/releases).

### 2021-11
- Time-domain simuations now return a result structure (non breaking)
- *Breaking*: `lsimplot, stepplot, impulseplot` have been replaced by `plot(lsim())` etc.
- *Breaking*: `pole, tzero` has been renamed to their plural form, `poles, tzeros`.
- *Breaking*: `c2d` now no longer returns the `x0map` for statespace systems, see function `c2d_x0map` for the old behavior.
- *Breaking*: The array layout of time and frequency responses has been transposed, i.e., in `y,t,x,u = lsim(sys, ...)`, the output arrays `y,x,u` are now of shape `size(y) == (sys.ny, T)`.
- New functions `observer_controller, observer_predictor, placePI`.
- *Breaking*: The type `LQG` has been removed, see [RobustAndOptimalControl.jl](https://github.com/JuliaControl/RobustAndOptimalControl.jl/blob/master/src/lqg.jl) for its replacement.
- *Breaking*: `balreal` and `baltrunc` return an additional value, the applied similarity transform.
- A large set of bug fixes
- For a full list of changes, [see here](https://github.com/JuliaControl/ControlSystems.jl/pull/565/commits).


## Documentation

All functions have docstrings, which can be viewed from the REPL, using for example `?tf `.

A documentation website is available at [http://juliacontrol.github.io/ControlSystems.jl/latest/](http://juliacontrol.github.io/ControlSystems.jl/latest/).

Some of the available commands are:
##### Constructing systems
ss, tf, zpk
##### Analysis
poles, tzeros, norm, hinfnorm, linfnorm, ctrb, obsv, gangoffour, margin, markovparam, damp, dampreport, zpkdata, dcgain, covar, gram, sigma, sisomargin
##### Synthesis
care, dare, dlyap, lqr, dlqr, place, leadlink, laglink, leadlinkat, rstd, rstc, dab, balreal, baltrunc
###### PID design
pid, stabregionPID, loopshapingPI, pidplots
##### Time and Frequency response
step, impulse, lsim, freqresp, evalfr, bode, nyquist
##### Plotting
lsimplot, stepplot, impulseplot, bodeplot, nyquistplot, sigmaplot, marginplot, gangoffourplot, pidplots, pzmap, nicholsplot, pidplots, rlocus, leadlinkcurve
##### Other
minreal, sminreal, c2d
## Usage

This toolbox works similar to that of other major computer-aided control
systems design (CACSD) toolboxes. Systems can be created in either a [transfer
function](http://en.wikipedia.org/wiki/Transfer_function) or a [state
space](http://en.wikipedia.org/wiki/State-space_representation) representation.
These systems can then be combined into larger architectures, simulated in both
time and frequency domain, and analyzed for stability/performance properties.

### Example

Here we create a simple position controller for an electric motor with an
inertial load.

```julia
using ControlSystems

# Motor parameters
J = 2.0
b = 0.04
K = 1.0
R = 0.08
L = 1e-4

# Create the model transfer function
s = tf("s")
P = K/(s*((J*s + b)*(L*s + R) + K^2))
# This generates the system
# TransferFunction:
#                1.0
# ---------------------------------
# 0.0002s^3 + 0.160004s^2 + 1.0032s
#
#Continuous-time transfer function model

# Create an array of closed loop systems for different values of Kp
CLs = TransferFunction[kp*P/(1 + kp*P) for kp = [1, 5, 15]];

# Plot the step response of the controllers
# Any keyword arguments supported in Plots.jl can be supplied
plot(step.(CLs, 5), label=["Kp = 1" "Kp = 5" "Kp = 15"])
```

![StepResponse](/example/step_response.png)

### Additional examples
See the examples folder and [ControlExamples.jl](https://github.com/JuliaControl/ControlExamples.jl/)
