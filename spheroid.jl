using DelimitedFiles
using Base.Threads
using LinearAlgebra
using PyPlot
using Colors

# Import the spherical_mesh module
include("spheroid_mesh.jl")

# Parameters for the kernel function
A, B, C, D = 0.14, 0.9, 1.2, 0.45

# Parameters for the model
θT = 1.0

# Parameters for plotting
vmin, vmax = -1.8, 2.7

# Parameters for the surface; add points times : N, radius : R, flattening : ϵ, semi-major axis : R1, semi-minor axis : R2
N, R, ϵ = 5, 1.0, 0.01
R1, R2 = R, R * (1 - ϵ)

# Time parameters
dt = 0.01
Time, STEP = Int(25000), Int(100)   # simulation time , steps b/w output

# Initial condition : The equation in Example 3.1 is used
function u_0(θ)
    return 2π * (A - C) * (1 - cos(θT)) + π * (B - 3D) * sin(θT)^2 * cos(θ) + 2π * C * (sin(θT)^2 * cos(θT) * cos(θ)^2 + (cos(θT)^3 - 3cos(θT) + 2) / 3) + π * D * (2 * (1 - cos(θT)^4) * cos(θ)^3 + 3 * sin(θT)^4 * cos(θ) * sin(θ)^2)
end

# Threshold : The equation (3.4) is used
u_T = u_0(θT)

# Kernel function
function K(d)
    return A+B*cos(d)+C*cos(2*d)+D*cos(3*d)
end

# Main calculation function of the scheme
function update(U, K_matrix, Areas)
    # For a quick calculation, we only update the points that satisfies U(θ)>u_T.
    tgts = []
    for i in 1:n
        if U[i] > u_T
            push!(tgts, i)
        end
    end

    # Update the points
    Threads.@threads for i in 1:n
        temp = 0.0
        for j in tgts
            temp += K_matrix[i, j] * Areas[j]
        end
        
        U_next[i] = U[i] + dt * (temp - U[i])
    end
end

# Get the points and simplices
n, simplices, points_3d = spherical_mesh_generator(N, R, R1, R2)
println("Number of points: ", n)

# Load the distance matrix
filename_dist = "R=$(R1)_eps=$(round(ϵ,digits=2)).txt"
if isfile(filename_dist)
    println("Loading distance matrix from $(filename_dist)...")
    dist = readdlm(filename_dist)
    println("Completed loading distance matrix")
else
    println("The file name $(filename_dist) does not exist.")
    exit()
end

# Calculate the areas of the triangles
println("Calculating areas...")
Areas = Area(n, simplices, points_3d)
println("Completed area calculation")

# Calculate the kernel matrix
K_matrix = zeros(n, n) 
println("Calculating kernel matrix...")
Threads.@threads for i in 1:n
    for j in 1:n
        K_matrix[i, j] = K(dist[i, j])
    end
end
println("Completed kernel calculation")

# Arrays to store time evolution
U_next = zeros(n)
U_time = zeros(Int(Time / STEP)+1, n)

# Initial condition
U = zeros(n)
for i in 1:n
    U[i] = u_0(dist[29,i])
end

PyPlot.ion()
fig = figure()
ax = fig.add_subplot(111, projection="3d")

function scatter_plot(U, points_3d, vmin, vmax)
    ax.cla()
    x, y, z = points_3d[:,1], points_3d[:,2], points_3d[:,3]
    sc = ax.scatter(x, y, z, 
                    c=U, cmap="viridis", vmin=vmin, vmax=vmax, 
                    edgecolors="none", s=20)
    ax.view_init(elev=19, azim=-10)
    ax.set_box_aspect([1, 1, 1]) 
    display(fig)
end

# Time evolution function
function time_evolution(U, K_matrix, Areas, U_time, STEP)
    count = 1
    for tt in 0:Time

        if tt % STEP == 0
            U_time[count, :] = U[:]
            scatter_plot(U, points_3d, vmin, vmax)
            println("printed time = $(tt),$(maximum(U)) $(minimum(U)) $(argmax(U)) $(argmin(U))")
            count += 1
        end
        
        update(U, K_matrix, Areas)
        U .= U_next
    end
end

# Main loop
println("Start time evolution")
time_evolution(U, K_matrix, Areas, U_time, STEP)
println("End time evolution")

# Save the data
filename_U = "U_eps=$(round(ϵ,digits=2)).txt"
open(filename_U, "w") do io
    writedlm(io, U_time)
end
filename_xyz = "xyz_eps=$(round(ϵ,digits=2)).txt"
open(filename_xyz, "w") do io
    writedlm(io, points_3d)
end