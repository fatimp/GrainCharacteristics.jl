# GrainCharacteristics.jl
[![CI](https://github.com/fatimp/GrainCharacteristics.jl/actions/workflows/test.yml/badge.svg)](https://github.com/fatimp/GrainCharacteristics.jl/actions/workflows/test.yml)

This library contains few functions used to characterize grains in grained
porous material, namely `equivalent_radius`, `sphericity`, `convexity` and
`elongation`. These functions can be applied to two- and three-dimensional
arrays with element type `Bool` where `true` means that the element belongs to
the grain and `false` means otherwise.
