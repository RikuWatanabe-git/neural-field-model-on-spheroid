using DelimitedFiles

# Import the spherical_mesh module
include("spheroid_mesh.jl")

py"""
# Parameters for the surface; add points times : N, flattening : eps, semi-major axis : R1 R, semi-minor axis : R2
N, R, eps = 5, 1.0, 0.01
R1, R2 = R, R * (1 - eps)

from geographiclib.geodesic import Geodesic
geod = Geodesic(R1, eps)
def distance(n, points_3d, R1, R2):
    theta = np.arccos(points_3d[:, 2] / R2) / np.pi * 180 - 90
    phi = np.arctan2(points_3d[:, 1], points_3d[:, 0]) / np.pi * 180

    distance_matrix = np.zeros((n, n))
    for i in range(n):
        if i%100 == 0:
            print("i = ", i)
        for j in range(i+1, n):
            distance_matrix[i, j] = geod.Inverse(theta[i], phi[i], theta[j], phi[j])['s12']
    distance_matrix += distance_matrix.T

    return distance_matrix
"""
distance = py"distance"

N, R, ϵ = py"N", py"R", py"eps"
R1, R2 = py"R1", py"R2"

filename_dist = "R=$(R1)_eps=$(round(ϵ,digits=2)).txt"
if isfile(filename_dist)
    println("Distance matrix already exists")
    exit()
end
println("filename = $(filename_dist)")

# Get the points and simplices
n, simplices, points_3d = spherical_mesh_generator(N, R, R1, R2)
println("Number of points: ", n)

# Main calculation
println("Calculating distance matrix...")
dist = distance(n, points_3d, R1, R2)

# Save the distance matrix
open(filename_dist, "w") do io
    writedlm(io, dist)
end
println("Completed distance calculation")
