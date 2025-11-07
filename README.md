# TSC-OPF Julia
This code solves the TSC-OPF using Julia Language.

To run the code and change the input parameters, use the file "main.jl"

The code stores Input Data in DataFrames. Some of the Output Data is also stored in DataFrames.

To install the required Julia packages, run:

```julia
using Pkg

Pkg.add("LinearAlgebra")
Pkg.add("SparseArrays")
Pkg.add("DataFrames")
Pkg.add("Printf")
Pkg.add("CSV")
Pkg.add("DataStructures")
Pkg.add("JuMP")
Pkg.add("Ipopt")
Pkg.add("Plots")
Pkg.add("LaTeXStrings")
Pkg.add("Measures")
Pkg.add("Trapz")
