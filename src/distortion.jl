# OpenCV Original License

#
#
#  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
#
#  By downloading, copying, installing or using the software you agree to this license.
#  If you do not agree to this license, do not download, install,
#  copy or use the software.
#
#
#                           License Agreement
#                For Open Source Computer Vision Library
#
# Copyright (C) 2000-2008, Intel Corporation, all rights reserved.
# Copyright (C) 2009, Willow Garage Inc., all rights reserved.
# Third party copyrights are property of their respective owners.
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
#
#   * Redistribution's of source code must retain the above copyright notice,
#     this list of conditions and the following disclaimer.
#
#   * Redistribution's in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#
#   * The name of the copyright holders may not be used to endorse or promote products
#     derived from this software without specific prior written permission.
#
# This software is provided by the copyright holders and contributors "as is" and
# any express or implied warranties, including, but not limited to, the implied
# warranties of merchantability and fitness for a particular purpose are disclaimed.
# In no event shall the Intel Corporation or contributors be liable for any direct,
# indirect, incidental, special, exemplary, or consequential damages
# (including, but not limited to, procurement of substitute goods or services;
# loss of use, data, or profits; or business interruption) however caused
# and on any theory of liability, whether in contract, strict liability,
# or tort (including negligence or otherwise) arising in any way out of
# the use of this software, even if advised of the possibility of such damage.
#
#M

using Images

"""

creategridcoordinates(xnum::Integer, ynum::Integer)::Matrix{<:Integer}

Creates a vector of grid 0-indexed coordinates up to xnum in the first dimension and up to ynum in the second dimension.

# Examples
```julia-repl
julia> creategridcoordinates(3,5)
```
"""
function creategridcoordinates(xnum::Integer, ynum::Integer)
    n = xnum * ynum
    res = zeros(Int32, n, 2) 
    for i = 0:(n - 1)
        res[i + 1, :] = [i ÷ ynum, i % ynum]
    end
    return res
end

"""

creategridcoordinates(xlength::Real, ylength::Real, xnum::Integer, ynum::Integer)

Creates a vector of grid coordinates which span the rowsize in the first dimension and the colsize in the second dimension.

# Arguments
- `xlength::Real`: Length to span in first dimension 
- `ylength::Real`: Length to span in second dimension 
- `xnum::Integer`: Number of points to use when spanning first dimension
- `ynum::Integer`: Number of points to use when spanning second dimension

# Examples
```julia-repl
julia> creategridcoordinates(640, 480, 10, 10)
```
"""
function creategridcoordinates(xlength::Real, ylength::Real, xnum::Integer, ynum::Integer)
    grid = creategridcoordinates(xnum, ynum)
    grid = grid .* [(xlength / (xnum - 1)) (ylength / (ynum - 1))]
    return grid
end


function distortpoints(src::AbstractMatrix{<:Real}; cameramatrix::AbstractMatrix{<:Real}, distcoeffs::AbstractVector{<:Real}, newcameramatrix::AbstractMatrix{<:Real})
    fx1 = newcameramatrix[1,1]
    fy1 = newcameramatrix[2,2] 
    cx1 = newcameramatrix[1,3]
    cy1 = newcameramatrix[2,3]
    ifx1 = 1 / fx1
    ify1 = 1 / fy1
    
    
    fx0 = cameramatrix[1,1]
    fy0 = cameramatrix[2,2] 
    cx0 = cameramatrix[1,3]
    cy0 = cameramatrix[2,3] 
    
    k1 = distcoeffs[1]
    k2 = distcoeffs[2]
    k3 = distcoeffs[5]
    p1 = distcoeffs[3]
    p2 = distcoeffs[4]
    
    # undo linear distortion 
    dst = (src .- [cx1 cy1]) .* [ifx1 ify1]
    
    # perform nonlinear distortion
    for i = 1:size(src)[1]
        x = dst[i, 1]
        y = dst[i, 2]
        
        x2 = x * x
        y2 = y * y
        
        r2 = x2 + y2
        r4 = r2 * r2
        r6 = r4 * r2
        
        rad = 1 + k1 * r2 + k2 * r4 + k3 * r6
        
        xp = x * rad + 2p1 * x * y + p2 * (r2 + 2x2)
        yp = y * rad + 2p2 * x * y + p1 * (r2 + 2y2)
        
        dst[i, :] = [xp yp]
    end
    
    # distort linearly
    dst = dst .* [fx0 fy0] .+ [cx0 cy0]
    return dst   
end

"""


    undistortpoints(src::AbstractMatrix{<:Real}; cameramatrix::AbstractMatrix{<:Real}, distcoeffs::Union{AbstractVector{<:Real}, Nothing}, newcameramatrix::Union{AbstractMatrix{<:Real}, Nothing})

Computes the ideal point coordinates from the observed point coordinates.

# Arguments
- `src::AbstractMatrix{<:Real}`: Observerd point coordinates, (N,2) 
- `cameramatrix::AbstractMatrix{<:Real}`: Camera matrix, (3,3)
- `distcoeffs::Union{AbstractVector{<:Real}, Nothing}`: Input vector of distortion coefficients, (5,) (k1, k2, p1, p2, k3)
- `newcameramatrix::Union{AbstractMatrix{<:Real}, Nothing}`: New camera matrix, (3,3). If the matrix is empty, the identity new camera matrix is used.


# Examples
```julia-repl
julia> undistortpoints([10 23; 100 100.3; 23.1 20.5], cameramatrix=[517.3 0 318.6; 0 516.5 244.3; 0 0 1], distcoeffs=[0.2624, -0.9531, -0.0054, 0.0026, 1.1633])
```
"""
function undistortpoints(src::AbstractMatrix{<:Real}; cameramatrix::AbstractMatrix{<:Real}, distcoeffs::Union{AbstractVector{<:Real},Nothing}, newcameramatrix::Union{AbstractMatrix{<:Real},Nothing})::AbstractMatrix{<:Real}
    # sanitize inputs
    if size(src)[2] != 2
        @error "src should have 2 columns"
        return
    end
    
    if size(cameramatrix) != (3, 3)
        @error "cameramatrix should be of size (3,3)"
        return
    end
    
    if newcameramatrix != nothing && size(newcameramatrix) != (3, 3)
        @error "newcameramatrix should be of size (3,3)"
        return
    end
    
    if distcoeffs != nothing && size(distcoeffs) != (5,)
        @error "distcoeffs must be of size (5,)"
        return
    end
    

    
    fx = cameramatrix[1,1]
    fy = cameramatrix[2,2] 
    cx = cameramatrix[1,3]
    cy = cameramatrix[2,3]

    ifx = 1 / fx
    ify = 1 / fy
    
    # undo linear distortion 
    dst = (src .- [cx cy]) .* [ifx ify]
    
    # undo nonlinear distortion numerically
    num_points = size(src)[1]
    if distcoeffs != nothing
    
        k1 = distcoeffs[1]
        k2 = distcoeffs[2]
        k3 = distcoeffs[5]
        p1 = distcoeffs[3]
        p2 = distcoeffs[4]
        
        for i = 1:num_points
            x0 = dst[i, 1]
            y0 = dst[i, 2]
            x = x0
            y = y0
            for j = 1:5
                r2 = x * x + y * y
                r4 = r2 * r2
                r6 = r4 * r2
                icdist = 1 / (1 + k1 * r2 + k2 * r4 + k3 * r6)
                deltaX = 2 * p1 * x * y + p2 * (r2 + 2 * x * x)
                deltaY = p1 * (r2 + 2 * y * y) + 2 * p2 * x * y
                x = (x0 - deltaX) * icdist
                y = (y0 - deltaY) * icdist
            end
            dst[i, :] = [x y]
        end
    end
    
    # linear distortion according to newcameramatrix
    if newcameramatrix != nothing
        fx = newcameramatrix[1,1]
        fy = newcameramatrix[2,2] 
        cx = newcameramatrix[1,3]
        cy = newcameramatrix[2,3]
        
        dst = dst .* [fx fy] .+ [cx cy] 
    end
    return dst
    
end

struct Rectangle 
    x::Real
    y::Real
    w::Real
    h::Real
end


"""


    getrectangles(cameramatrix::AbstractMatrix{<:Real}, distcoeffs::Union{AbstractVector{<:Real}, Nothing}, newcameramatrixUnion::{AbstractVector{<:Real}, Nothing}, imgsize::Tuple{Integer, Integer})

Compute the inscribed and circumscribed rectangle of a set of points in an image after undistorting them.

Returns outerrectangle, innerrectangle

# Arguments
- `cameramatrix::AbstractMatrix{<:Real}`: Camera matrix, (3,3)
- `distcoeffs::Union{AbstractVector{<:Real}, Nothing}`: Input vector of distortion coefficients, (5,) (k1, k2, p1, p2, k3)
- `newcameramatrix::Union{AbstractMatrix{<:Real}, Nothing}`: New camera matrix, (3,3). If the matrix is empty, the identity new camera matrix is used.
- `imgsize::Tuple{Integer, Integer}`: size of image represented as a tuple

# Examples
```julia-repl
julia> getrectangles(cameramatrix=[517.3 0 318.6; 0 516.5 244.3; 0 0 1], distcoeffs=[0.2624, -0.9531, -0.0054, 0.0026, 1.1633], imgsize=(640, 480))
```
"""
function getrectangles(;cameramatrix::AbstractMatrix{<:Real}, distcoeffs::Union{AbstractVector{<:Real},Nothing} = nothing, newcameramatrix::Union{AbstractMatrix{<:Real},Nothing} = nothing, imgsize::Tuple{Integer,Integer})::Tuple{Rectangle,Rectangle}
    N = 9
    pts = creategridcoordinates(imgsize[1], imgsize[2], N, N)
    pts = undistortpoints(pts, cameramatrix = cameramatrix, distcoeffs = distcoeffs, newcameramatrix = newcameramatrix)
    
    
    iX0 = -Inf; iX1 = Inf; iY0 = -Inf; iY1 = Inf
    oX0 = Inf; oX1 = -Inf; oY0 = Inf; oY1 = -Inf
    
    i = 1
    for x = 1:N
        for y = 1:N
            px = pts[i,1]
            py = pts[i,2]
            
            oX0 = min(oX0, px)
            oX1 = max(oX1, px)
            oY0 = min(oY0, py)
            oY1 = max(oY1, py)
            
            if x == 1
                iX0 = max(iX0, px)
            end
            if x == N
                iX1 = min(iX1, px)
            end
            if y == 1
                iY0 = max(iY0, py)
            end
            if y == N 
                iY1 = min(iY1, py)
            end
            i += 1
        end
    end
    return Rectangle(oX0, oY0, oX1 - oX0, oY1 - oY0), Rectangle(iX0, iY0, iX1 - iX0, iY1 - iY0)
    
end

"""


    getoptimal_newcameramatrix(;cameramatrix::AbstractMatrix{<:Real}, distcoeffs::Union{AbstractVector{<:Real}, Nothing} = nothing, imgsize::Tuple{Integer, Integer}, alpha::Real, newimgsize::Tuple{Integer, Integer})
 
Compute the optimal camera matrix such that after removing nonlinear distortion the image is remapped nicely.

# Arguments
- `cameramatrix::AbstractMatrix{<:Real}`: Camera matrix, (3,3)
- `distcoeffs::Union{AbstractVector{<:Real}, Nothing}`: Input vector of distortion coefficients, (5,) (k1, k2, p1, p2, k3)
- `newcameramatrix::Union{AbstractMatrix{<:Real}, Nothing}`: New camera matrix, (3,3). If the matrix is empty, the identity new camera matrix is used.
- `imgsize::Tuple{Integer, Integer}`: size of image represented as a tuple
- `alpha::Real`: interpolation factor between "remap without missing regions" to "remap fitting all data"
- `newimgsize::Tuple{Integer, Integer}`: desired size of image after remapping

"""
function get_optimalnewcameramatrix(;cameramatrix::AbstractMatrix{<:Real}, distcoeffs::Union{AbstractVector{<:Real},Nothing} = nothing, imgsize::Tuple{Integer,Integer}, alpha::Real, newimgsize::Tuple{Integer,Integer})
    
    # Get inscribed and circumscribed rectangles in normalized
    outer, inner = getrectangles(cameramatrix = cameramatrix, distcoeffs = distcoeffs, imgsize = imgsize)
    
    # Projection mapping inner rectangle to viewport
    fx0 = newimgsize[1] / inner.w
    fy0 = newimgsize[2] / inner.h
    cx0 = -fx0 * inner.x
    cy0 = -fy0 * inner.y
    
    # Projection mapping outer rectangle to viewport
    fx1 = newimgsize[1] / outer.w
    fy1 = newimgsize[2] / outer.h
    cx1 = -fx1 * outer.x
    cy1 = -fy1 * outer.y
    
    # Interpolate between the two optimal projections
    M = zeros(3, 3) # optimal camera matrix
    M[1,1] = fx0 * (1 - alpha) + fx1 * alpha
    M[2,2] = fy0 * (1 - alpha) + fy1 * alpha
    M[1,3] = cx0 * (1 - alpha) + cx1 * alpha
    M[2,3] = cy0 * (1 - alpha) + cy1 * alpha
    M[3,3] = 1
    return M
end

"""


    getoptimal_newcameramatrix(;cameramatrix::AbstractMatrix{<:Real}, distcoeffs::Union{AbstractVector{<:Real}, Nothing} = nothing, imgsize::Tuple{Integer, Integer}, alpha::Real, newimgsize::Tuple{Integer, Integer})
 
Computes the undistortion and rectification transformation map. Returns a tuple of maps in the form of matrices.

# Returns
- mapx: (u,v) -> x
- mapy: (u,v) -> y

# Arguments
- `cameramatrix::AbstractMatrix{<:Real}`: Camera matrix, (3,3)
- `distcoeffs::Union{AbstractVector{<:Real}, Nothing}`: Input vector of distortion coefficients, (5,) (k1, k2, p1, p2, k3)
- `newcameramatrix::Union{AbstractMatrix{<:Real}, Nothing}`: New camera matrix, (3,3). If the matrix is empty, the identity new camera matrix is used.
- `imgsize::Tuple{Integer, Integer}`: size of image represented as a tuple

"""
function init_undistortrectifymap(;cameramatrix::AbstractMatrix{<:Real}, distcoeffs::Union{AbstractVector{<:Real},Nothing}, newcameramatrix::AbstractMatrix{<:Real}, imgsize::Tuple{Integer,Integer})
    fx0 = cameramatrix[1,1]
    fy0 = cameramatrix[2,2] 
    cx0 = cameramatrix[1,3]
    cy0 = cameramatrix[2,3]
    ifx0 = 1 / fx0
    ify0 = 1 / fy0
    
    
    fx1 = newcameramatrix[1,1]
    fy1 = newcameramatrix[2,2] 
    cx1 = newcameramatrix[1,3]
    cy1 = newcameramatrix[2,3]
    ifx1 = 1 / fx1
    ify1 = 1 / fy1
    
    k1 = distcoeffs[1]
    k2 = distcoeffs[2]
    k3 = distcoeffs[5]
    p1 = distcoeffs[3]
    p2 = distcoeffs[4]
    
    xmap = zeros(Real, imgsize[2], imgsize[1])
    ymap = zeros(Real, imgsize[2], imgsize[1])
    
    
    for i = 0:(imgsize[1] - 1)
        x = (i - cx1) * ifx1
        x2 = x * x
        for j = 0:(imgsize[2] - 1)
            y = (j - cy1) * ify1
            y2 = y * y
            
            r2 = x2 + y2
            r4 = r2 * r2
            r6 = r2 * r4
            rad = (1 + k1 * r2 + k2 * r4 + k3 * r6)
            
            xp = x * rad + 2p1 * x * y + p2 * (r2 + 2x2)
            yp = y * rad + 2p2 * x * y + p1 * (r2 + 2y2)
            
            xp = xp * fx0 + cx0
            yp = yp * fy0 + cy0
            
            xmap[j + 1, i + 1] = xp
            ymap[j + 1, i + 1]  = yp
        end
    end
    return xmap, ymap
end

"""


remap(src::AbstractArray{T,2}, xmap::AbstractArray{<:Real,2}, ymap::AbstractArray{<:Real,2}) where T
 
Applies a generic geometrical transformation to an image.

# Arguments
- `src::AbstractMatrix{T,2}`: source image
- `xmap::AbstractArray{<:Real,2}: map taking (u,v) in destination to x in src
- `yxmap::AbstractArray{<:Real,2}: map taking (u,v) in destination to y in src

"""
function remap(src::AbstractArray{T,2}, xmap::AbstractArray{<:Real,2}, ymap::AbstractArray{<:Real,2}) where T
    dst = zeros(T, size(xmap))
    for j = 1:size(xmap)[1]
        for i = 1:size(xmap)[2]
            dst[j,i] = bilinear_interpolation(src, ymap[j, i], xmap[j, i])
        end
    end
    dst
end