'''
`tide` module contains the classes and functions to -

Submodules:

- `tide.atlas`
- `tide.interpolate`
- `tide.utilities`

Functionalities:

- Read tidal atlas dataset
- Interpolate in a grid, or list of points (x, y)
- TODO: Reconstruct a timeseries from a set of constituents and their amplitude and phase
    - In a time-series
    - In a grid
- Compute the nodal factor for a given time

# Generation of a bctides
1. Read hgrid.gr3 file
2. Get the boundary points from hgrid.gr3
3. Get a list of potential waves, amplitude, phase
    1. Correct for nodal before writing
4. Get a list of tidal waves, amplitude, phase
    1. Correct for nodal before writing
5. Read timeseries for a boundary (repeat for each flow boundary)
6. Write bctides, for a given time, and runtime
7. Write flux.th for a given time, and runtime

# Accessing tidal atlas
- interp_xy() -> dict(const:[amp, pha])
- interp_ds() -> dict(const:da)

'''
