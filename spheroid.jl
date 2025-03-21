using DelimitedFiles
using PyCall
using Base.Threads
using LinearAlgebra

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

# Define plot function via PyCall
py"""
import numpy as np
import matplotlib.pyplot as plt
fig3d = plt.figure()
ax3d = fig3d.add_subplot(111, projection='3d')

def scatter_plot(U_time,tt,count,points_3d,R1,R2,dt,vmin,vmax):
    ax3d.cla()
    lim=max(R1,R2)
    ax3d.scatter3D(points_3d[:,0], points_3d[:,1], points_3d[:,2],
                c=U_time[count-1],s=20, vmin=vmin, vmax=vmax) 
    ax3d.set_xlim(-lim,lim), ax3d.set_ylim(-lim,lim), ax3d.set_zlim(-lim,lim)
    ax3d.set_box_aspect((1,1,1))
    ax3d.view_init(elev=19, azim=-10)
    ax3d.set_title(f'time: {tt*dt:.2f}:, max index: {np.argmax(U_time[count-1]):.2f}, max: {np.max(U_time[count-1]):.2f}, min: {np.min(U_time[count-1]):.2f}')
    plt.pause(0.001)
"""
scatter_plot=py"scatter_plot"

# Time evolution function
function time_evolution(U, K_matrix, Areas, U_time, STEP)
    count = 1
    for tt in 0:Time

        if tt % STEP == 0
            U_time[count, :] = U[:]
            scatter_plot(U_time,tt,count,points_3d,R1,R2,dt,vmin,vmax)
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
