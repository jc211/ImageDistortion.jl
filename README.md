# Image Distortion

Utility functions to distort and undistort an image or a set of points. The code is for the most part ported over from OpenCV to native Julia.

## Installation

```julia
Pkg.add("https://github.com/jc211/ImageDistortion.jl)"
using ImageDistortion
```


## Disortion Models

OpenCV implements the following [distortion model](https://docs.opencv.org/2.4/modules/calib3d/doc/camera_calibration_and_3d_reconstruction.html)

![Distortion Equations](https://github.com/jc211/ImageDistortion.jl/blob/master/img/distortionmodel.PNG "OpenCV distortion formulas")

This implementation currently only supports the k1, k2, k3, p1, p2 parameters. They are represented by distcoeffs where a sample distcoeffs is in the form

```julia
distcoeffs = [k1, k2, p1, p2, k3]
```

## Overview

![Overview](https://github.com/jc211/ImageDistortion.jl/blob/master/img/overview.PNG "Overview")

In general, this package has functions to move around the blocks illustrated above. 

## Usage

### Distorting Points

![](https://github.com/jc211/ImageDistortion.jl/blob/master/img/distortpoints.PNG "Distortion")

The function distortpoints takes a matrix of shape (n,2) representing points (e.g. [0 0; 10 15; 32 10]) and firstly undoes the linear distortion provided by newcameramatrix and subsequently applies nonlinear distortion by distcoeffs, followed by a linear distortion by cameramatrix

```julia
points = [10 23; 100 100.3; 23.1 20.5]
distortedpoints = distortpoints(points, cameramatrix=[517.3 0 318.6; 0 516.5 244.3; 0 0 1], distcoeffs=[0.2624, -0.9531, -0.0054, 0.0026, 1.1633], cameramatrix=[517.3 0 318.6; 0 516.5 244.3; 0 0 1])
```

### Undistorting Points

![](https://github.com/jc211/ImageDistortion.jl/blob/master/img/undistortpoints.PNG "Undistort Points")

The function undistortpoints takes a matrix of shape (n,2) representing points (e.g. [0 0; 10 15; 32 10]) and undoes the linear distortion given by cameramatrix, undoes the nonlinear distortion by distcoeffs, and optionally reapplies a linear distortion by newcameramatrix

```julia
points = [10 23; 100 100.3; 23.1 20.5]
distortedpoints = undistortpoints(points, cameramatrix=[517.3 0 318.6; 0 516.5 244.3; 0 0 1], distcoeffs=[0.2624, -0.9531, -0.0054, 0.0026, 1.1633])
```

### Get Optimal Camera Matrix

![](https://github.com/jc211/ImageDistortion.jl/blob/master/img/optimalcameramatrix.PNG "Optimal Camera Matrix")

The function get_optimalnewcameramatrix outputs a new camera matrix that transforms an undistorted image (in normalized coordinates) to one that contains all the information in the original source image if the parameter alpha = 0, or to one that has no missing information (i.e. no black parts in the image) if alpha = 1. A value of alpha in between does some mixture of the two extremes.

```julia
K_original = [517.3 0 318.6; 0 516.5 244.3; 0 0 1] # camera calibration of raw image
d_orignal = [0.2624, -0.9531, -0.0054, 0.0026, 1.1633]; # distortion coefficients of raw image

newcameramatrix = get_optimalnewcameramatrix(cameramatrix=K_original, distcoeffs=d_orignal, imgsize=(640,480), alpha=0, newimgsize=(640,480))
```

### Precompute Distortion Maps

![](https://github.com/jc211/ImageDistortion.jl/blob/master/img/undistortrectifymap.PNG "Optimal Camera Matrix")

The function init_undistortrectifymap precomputes the map shown above.

```julia
K_original = [517.3 0 318.6; 0 516.5 244.3; 0 0 1] # camera calibration of raw image
d_orignal = [0.2624, -0.9531, -0.0054, 0.0026, 1.1633]; # distortion coefficients of raw image

newcameramatrix = get_optimalnewcameramatrix(cameramatrix=K_original, distcoeffs=d_orignal, imgsize=(640,480), alpha=0, newimgsize=(640,480))
xmap, ymap = init_undistortrectifymap(cameramatrix=K_original, distcoeffs=d_orignal, newcameramatrix=newcameramatrix, imgsize=(640,480))

```

### Remap

The function remap uses bilinear interpolation to generate the domain of the maps which are generally computed by init_undistortrectifymap

```julia
K_original = [517.3 0 318.6; 0 516.5 244.3; 0 0 1] # camera calibration of raw image
d_orignal = [0.2624, -0.9531, -0.0054, 0.0026, 1.1633]; # distortion coefficients of raw image

newcameramatrix = get_optimalnewcameramatrix(cameramatrix=K_original, distcoeffs=d_orignal, imgsize=(640,480), alpha=0, newimgsize=(640,480))
xmap, ymap = init_undistortrectifymap(cameramatrix=K_original, distcoeffs=d_orignal, newcameramatrix=newcameramatrix, imgsize=(640,480))
newimg = remap(src, xmap, ymap)
```