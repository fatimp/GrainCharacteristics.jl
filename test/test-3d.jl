expected_sphericity = Dict{String, Float64}(
    "3d/tetrahedron.npy"  => 0.671,
    "3d/cube.npy"         => 0.806,
    "3d/octahedron.npy"   => 0.846,
    "3d/dodecahedron.npy" => 0.910
)

function ellipsoid(size, r, k)
    res = falses(size, size, size)
    for idx in CartesianIndices(res)
        x = idx[1] - size/2
        y = idx[2] - size/2
        z = idx[3] - size/2
        
        x1 =  cos(π/4)*x + sin(π/4)*y
        y1 = -sin(π/4)*x + cos(π/4)*y
        z1 = z
        
        if ((k*x1)^2 + y1^2 + z1^2 < r^2)
            res[idx] = true
        end
    end
    return res
end

@testcase "Test sphericity" begin
    Base.GC.gc()
    for (name, expected) in expected_sphericity
        @test isapprox(expected, name |> npzread |> G.sphericity; rtol = 0.04)
    end
end

@testcase "Test convexity" begin
    Base.GC.gc()
    for (name, _) in expected_sphericity
        @test isapprox(1, name |> npzread |> G.convexity; rtol = 0.04)
    end

    Base.GC.gc()
    box = falses(200, 200, 200)
    box[51:150, 51:150, 51:150] .= true
    box[76:125, 76:125, 76:125] .= false
    @test isapprox(G.convexity(box), 0.875; rtol = 0.04)

    Base.GC.gc()
    box = falses(200, 200, 200)
    box[51:150, 51:150, 51:150] .= true
    box[56:145, 56:145, 56:145] .= false
    box[61:140, 61:140,  1:155] .= false
    @test isapprox(G.convexity(box), 0.207; rtol = 0.04)
end

@testcase "Test elongation" begin
    Base.GC.gc()
    for k in [1, 2, 3]
        @test isapprox(1/k, ellipsoid(200, 80, k) |> G.elongation; rtol = 0.04)
    end

    Base.GC.gc()
    for k in [1, 1/2, 1/3]
        @test isapprox(k, ellipsoid(200, 30, 1/k) |> G.elongation; rtol = 0.04)
    end
end
