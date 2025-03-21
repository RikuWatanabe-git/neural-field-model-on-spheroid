using PyCall

py"""
import numpy as np
import scipy.spatial as ss

def spherical_mesh_generator(N, R, R1,R2):
    G = ( 1+np.sqrt(5) ) / 2
    r = np.sqrt( 1+G**2 )
    icosahedron = R / r * np.array([ [1, G, 0], [-1, -G, 0], [-1, G, 0], [1, -G, 0],
                                     [0, 1, G], [0, -1, -G], [0, -1, G], [0, 1, -G],
                                     [G, 0, 1], [-G, 0, -1], [G, 0, -1], [-G, 0, 1] ])

    points_3d = np.copy(icosahedron)
    n = len(points_3d)

    CH = ss.ConvexHull(points_3d)
    simplices1 = CH.simplices

    # add points
    for ii in range(N):
        new_points = []
        print( 'add points:', ii+1, 'times' )
        for i in range( len(simplices1[:,0]) ):
            v0, v1, v2 = simplices1[i,0], simplices1[i,1], simplices1[i,2]

            new_points.append([ (points_3d[v0,0] + points_3d[v1,0])/2,
                                (points_3d[v0,1] + points_3d[v1,1])/2,
                                (points_3d[v0,2] + points_3d[v1,2])/2 ])

            new_points.append([ (points_3d[v2,0] + points_3d[v1,0])/2,
                                (points_3d[v2,1] + points_3d[v1,1])/2,
                                (points_3d[v2,2] + points_3d[v1,2])/2 ])

            new_points.append([ (points_3d[v0,0] + points_3d[v2,0])/2,
                                (points_3d[v0,1] + points_3d[v2,1])/2,
                                (points_3d[v0,2] + points_3d[v2,2])/2 ])

        new_points = np.array( list(map(list, set(map(tuple, new_points)))) )

        nn = len(new_points)

        x, y, z = np.zeros(nn), np.zeros(nn), np.zeros(nn)
        THETA, PHI = np.zeros(nn), np.zeros(nn)

        THETA[:] = np.arccos(new_points[:,2] / np.linalg.norm(new_points[:,:], axis=1, ord=2))
        PHI[:]= np.arctan2(new_points[:,1], new_points[:,0]) + np.pi

        x[:] = R * np.cos(PHI[:]) * np.sin(THETA[:])
        y[:] = R * np.sin(PHI[:]) * np.sin(THETA[:])
        z[:] = R * np.cos(THETA[:])

        new_points = np.array( [ [x[k], y[k], z[k]] for k in range(nn) ] )

        points_3d = np.vstack( [points_3d, new_points] )

        CH = ss.ConvexHull(points_3d)
        simplices1  = CH.simplices

    n = len(points_3d)

    # rotation to avoid the singularity at the poles (z = +-R) in the 2d projection
    theta_x = 1e-4
    Rx = np.array([ [1, 0, 0],
                    [ 0, np.cos(theta_x), np.sin(theta_x) ],
                    [ 0,-np.sin(theta_x), np.cos(theta_x) ] ])

    for i in range(n):
        points_3d[i,:] = np.dot(Rx[:,:], points_3d[i,:])

    # 2d projection
    X, Y = np.zeros(n), np.zeros(n)
    X[:] = 2 * R * points_3d[:,0] / (R - points_3d[:,2])
    Y[:] = 2 * R * points_3d[:,1] / (R - points_3d[:,2])

    points = np.array( [[X[k], Y[k]] for k in range(n)] )

    tri    = ss.Delaunay(points, furthest_site = False)
    tri2   = ss.Delaunay(points, furthest_site = True)

    simplices = np.vstack([tri.simplices, tri2.simplices])
    simplices += 1

    x, y, z = np.zeros(n), np.zeros(n), np.zeros(n)
    THETA, PHI = np.zeros(n), np.zeros(n)

    THETA[:] = np.arccos(points_3d[:,2] / np.linalg.norm(points_3d[:,:], axis=1, ord=2))
    PHI[:]   = np.arctan2(points_3d[:,1], points_3d[:,0]) + np.pi

    x[:] =  R1 * np.cos(PHI[:]) * np.sin(THETA[:])
    y[:] =  R1 * np.sin(PHI[:]) * np.sin(THETA[:])
    z[:] =  R2 * np.cos(THETA[:])

    points_3d[:] = np.array([[x[k], y[k], z[k]] for k in range(n)])[:]

    return n, simplices, points_3d, THETA, PHI
"""
spherical_mesh_generator = py"spherical_mesh_generator"

function Area(n, simplices, points_3d)
    area = zeros(Float64, n)

    count = 0

    Threads.@threads for i in 1:n

        index = findall(simplices .== i)

        area_seg = 0.0
        m = 1

        for ind in index
            simplex_index = ind[1]
            simplex_pos = ind[2]

            # sort the vertices
            ordered_simplex = circshift(simplices[simplex_index, :], -simplex_pos)

            # get the vertices
            v0, v1, v2 = ordered_simplex[1], ordered_simplex[2], ordered_simplex[3]

            # calculate the length of each edge
            length0 = norm(points_3d[v1, :] - points_3d[v2, :])
            length1 = norm(points_3d[v2, :] - points_3d[v0, :])
            length2 = norm(points_3d[v0, :] - points_3d[v1, :])

            # calculate the angle of each edge
            angle1 = acos(dot(points_3d[v0, :] - points_3d[v1, :], points_3d[v2, :] - points_3d[v1, :]) / (length2 * length0))
            angle2 = acos(dot(points_3d[v0, :] - points_3d[v2, :], points_3d[v1, :] - points_3d[v2, :]) / (length1 * length0))

            # calculate the area of the triangle
            area_seg += (length1^2) / tan(angle1) + (length2^2) / tan(angle2)

            m += 1
            count += 1
        end

        # sum the area of the triangle
        area[i] = area_seg / 8.0
    end

    return area
end