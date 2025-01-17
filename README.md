# GalerkinEtena
Discontinuous Galerkin code in written in Julia based on the Matlab code from [nodal-dg](https://github.com/tcew/nodal-dg)

![Build Status](https://github.com/aitzkora/GalerkinEtena.jl/actions/workflows/main.yml/badge.svg)
[![codecov](https://codecov.io/gh/aitzkora/GalerkinEtena.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/aitzkora/GalerkinEtena.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://aitzkora.github.io/GalerkinEtena.jl/dev)

# Current Status
- 1D problems : advection and Maxwell's equation
- 2D : implementation not yet achieved

# to use it

```bash
git clone https://github.com/aitzkora/GalerkinEtena.jl galerkin
julia
julia> using Pkg
julia> Pkg.activate("galerkin")
  Activating project at `~/galerkin`
julia> using GalerkinEtena
julia> include("galerkin/src/main_advec1D.jl")
9×10 Matrix{Float64}:
 -7.18383e-11  0.198669  0.389418  0.564642  0.717356  0.841471  0.932039  0.98545   0.999574  0.973848
  0.010024     0.208484  0.398632  0.572887  0.724304  0.846845  0.935625  0.987104  0.999231  0.971521
  0.0322758    0.230198  0.418943  0.590987  0.739469  0.858471  0.943249  0.990422  0.99811   0.966007
  0.0636452    0.260643  0.44725   0.616026  0.760244  0.874153  0.953212  0.994269  0.995689  0.957413
  0.0998334    0.29552   0.479426  0.644218  0.783327  0.891207  0.963558  0.997495  0.991665  0.9463
  0.13589      0.330008  0.510969  0.67156   0.805377  0.907087  0.972634  0.999405  0.986334  0.93394
  0.166933     0.359488  0.53771   0.694496  0.823594  0.919858  0.97945   0.999995  0.980673  0.922255
  0.188835     0.380166  0.556341  0.710336  0.836013  0.92836   0.983696  0.999816  0.976076  0.913423
  0.198669     0.389418  0.564642  0.717356  0.841471  0.932039  0.98545   0.999574  0.973848  0.909297
```
