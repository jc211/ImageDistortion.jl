module ImageDistortion

using('distortion.jl')

export
    distortpoints,
    undistortpoints,
    getrectangles,
    get_optimalnewcameramatrix,
    init_undistortrectifymap,
    remap

"""
Algorithms:
    - Distortion: `distortpoints`, `undistortpoints`, `get_optimalnewcameramatrix`, `init_undistortrectifymap`, `remap`
"""
end