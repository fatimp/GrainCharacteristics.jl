#
# A "grain" in this code is an array of booleans A for which a set X
# such that ∀x: A[x] = true ↔ x ∈ X is a (connected) set of points
# which represents an object we want to analyze.
#

function array_of_points(grain :: AbstractArray{Bool})
    points = findall(grain)
    return Array([idx[k] for k in 1:length(points[1]), idx in points]')
end

"""
    equivalent_radius(grain)

Return radius of a ball which has the same volume as `grain`.
"""
function equivalent_radius(grain :: AbstractArray{Bool})
    dim = ndims(grain)
    volume = sum(grain)
    return (volume * SF.gamma(dim/2 + 1))^(1/dim) / sqrt(π)
end

"""
    convexity(grain)

Return convexity of the grain, i.e. its volume devided by the volume
of its convex hull.
"""
function convexity(grain :: AbstractArray{Bool})
    hull = grain |> array_of_points |> QH.chull
    return sum(grain) / hull.volume
end

# Sphericity function gives quite big relative error because only
# exterior voxels contribute to the surface and neither there is a way
# we can reconstruct the original surface from edgy voxel data, nor
# can we unambiguously select voxels which belong to the surface. Here
# we use `extract_edges` function from our CorrelationFunctions.jl
# library to extract the surface which gives about 10% error for both
# round and low-poly objects. Another way is to use polygonisation
# algorithms like marching cubes (look at Meshing.jl). In this way we
# get nearly 0% error for low-poly objects, but too high error for
# round objects.
"""
    sphericty(grain; kernel = ConvKernel(5)

Return sphericity of the grain. Sphericity is the surface of a sphere
which has the same volume as the grain divided by the surface of the
grain.

`kernel` is a filter (either `ConvKernel` or `ErosionKernel`) which is
used to extract surface of the grain from an input image. For more
information look at `CorrelationFunctions.jl` library.

This function is also defined for other dimensions.
"""
function sphericity(grain :: AbstractArray{Bool}; kernel = U.ConvKernel(5))
    edge = U.extract_edges(grain, kernel, U.Periodic())
    surf = sum(edge)
    vol  = sum(grain)
    n    = ndims(grain)
    α    = (n-1) / n

    sphere_surf = vol^α * 2*sqrt(π)/SF.gamma(n/2)^(1/n) * (n/2)^α
    return sphere_surf / surf
end

function col_means(array :: AbstractMatrix)
    map(axes(array, 2)) do ax
        S.mean(array[:, ax])
    end
end

function covariance_matrix(grain :: AbstractArray{Bool})
    points = array_of_points(grain)
    m = col_means(points)
    dim = ndims(grain)
    mat = mapreduce(+, axes(points, 1)) do idx
        coord = points[idx, :] - m
        map(Iterators.product(1:dim, 1:dim)) do (i, j)
            coord[i]*coord[j]
        end
    end

    # Construct symmetric matrix to indicate that all eigenvalues are
    # real.
    sym = LA.Symmetric(mat)
    @assert sym ≈ mat
    return sym
end

"""
    elongation(grain)

Return elongation of the grain. A possible value of elongation is in
the range [0, 1].
"""
function elongation(grain :: AbstractArray{Bool})
    m = covariance_matrix(grain)
    eig = LA.eigvals(m)
    @assert all(x -> x > 0, eig)
    return sqrt(eig[begin]/eig[end])
end
