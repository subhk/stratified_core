## Stratified core 

This repository contains the Julia scripts used to generate the reduced-model figures for the manuscript:

**A stratified core lid as a spectral selector for trapped MAC modes**  

The code implements the reduced trapped-MAC-mode model described in the paper. It computes the trapped lid spectrum, modal damping, reduced envelope dynamics, single-mode analytic comparison, and detuning response used in the figures.

## Contents

- `paper_figures_common.jl`  
  Shared model definitions and numerical utilities. This file defines the baseline parameters, trapped-mode frequencies, damping rates, coupling coefficients, RK4 envelope solver, diagnostics, and exact single-mode solution.

- `fig_schematic.jl`  
  Generates the schematic geometry figure.

- `fig2_spectrum.jl`  
  Generates the trapped-mode spectrum and mode-resolvability figure.

- `fig3_transient.jl`  
  Generates the resonant transient response figure.

- `fig4_detuning.jl`  
  Generates the detuning-response figure.

## Requirements

The scripts were written for Julia and use CairoMakie for plotting.

Tested with:

- Julia 1.10
- CairoMakie
- LaTeXStrings

Install the plotting dependencies from the Julia REPL with:

```julia
using Pkg
Pkg.add(["CairoMakie", "LaTeXStrings"])
```

## Reproducing the figures

Run the scripts from this directory:

```bash
julia fig_schematic.jl
julia fig2_spectrum.jl
julia fig3_transient.jl
julia fig4_detuning.jl
```

## Model summary

The shared file `paper_figures_common.jl` implements the following reduced quantities:

```julia
kz(n, p)      = (n + 0.5) * pi / p.Hs
omegaA(n, p) = p.vA * kz(n, p)
omegaI(n, p) = 2 * p.Omega * kz(n, p) / k_total(n, p)
omegaN(n, p) = p.Ns * kh(p) / k_total(n, p)
omegaM(n, p) = abs(omegaA(n, p)) * sqrt(omegaA(n, p)^2 + omegaN(n, p)^2) / abs(omegaI(n, p))
```

The damping model used in the figures is

```julia
gamma_mode(n, p) = p.eta * kz(n, p)^2
```

The envelope equations are integrated by `rk4_envelope`, and the exact two-mode comparison is evaluated by `exact_single_mode`.

##  Parameters

The default parameter set is defined in `ModelParams`:

```julia
Omega    = 7.3e-5      # s^-1
vA       = 1.0e-2      # m s^-1
Ns       = 7.3e-5      # s^-1
Hs       = 100.0e3     # m
eta      = 1.0         # m^2 s^-1
lambda_h = 2.0e6       # m
```

The illustrative coupling profile is defined in `CouplingParams`. It is not a full interior-lid boundary-value solution; it is a reduced overlap model used to demonstrate the resonance mechanism.



