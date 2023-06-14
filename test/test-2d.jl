expected_roundness = Dict{String, Float64}(
    "2d/tria.npy"   => sqrt(π*sqrt(3))/3,
    "2d/square.npy" => sqrt(π)/2
)

function ellipse(size, r, k)
    res = falses(size, size)
    for idx in CartesianIndices(res)
        x = idx[1] - size/2
        y = idx[2] - size/2
        
        x1 =  cos(π/4)*x + sin(π/4)*y
        y1 = -sin(π/4)*x + cos(π/4)*y
        
        if ((k*x1)^2 + y1^2 < r^2)
            res[idx] = true
        end
    end
    return res
end

@testcase "Test roundness" begin
    Base.GC.gc()
    for (name, expected) in expected_roundness
        @test isapprox(expected, name |> npzread |> G.sphericity; rtol = 0.04)
    end
end

@testcase "Test convexity" begin
    Base.GC.gc()
    for (name, _) in expected_roundness
        @test isapprox(1, name |> npzread |> G.convexity; rtol = 0.04)
    end
end

@testcase "Test elongation" begin
    Base.GC.gc()
    for k in [1, 2, 3]
        @test isapprox(1/k, ellipse(200, 80, k) |> G.elongation; rtol = 0.04)
    end

    Base.GC.gc()
    for k in [1, 1/2, 1/3]
        @test isapprox(k, ellipse(200, 30, 1/k) |> G.elongation; rtol = 0.04)
    end
end
