using NPZ
using XUnit
import GrainCharacteristics as G

@testset "Test 3D objects" begin include("test-3d.jl") end
@testset "Test 2D objects" begin include("test-2d.jl") end
