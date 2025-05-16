
# Distributed computing

```@meta
CurrentModule = GridapEmbedded.Distributed
```

We support distributed computing through [GridapDistributed.jl](https://github.com/gridap/GridapDistributed.jl). As per usual, we design our libraries so that the high-level API is unchanged when using distributed computing. This means that for most users, the changes to your driver will be minimal.

The following features are currently supported:

- Level-Set Cutters
- STL Cutters

The folowing features are not yet supported:

- Aggregated FEM
- Moment-Fitted Quadratures

```@autodocs
Modules = [Distributed,]
Order   = [:type, :constant, :macro, :function]
```
