# OpenFoam Fast Flow Approximation
Solver (`invEikonalFoam`) approximates the inverse distance field **G** by solving the inverse Eikonal equation introduced by Fares and Schröder (2002).
The distance field (/ fill time) can be obtained via postprocessing using a Paraview filter ([Git repository](https://git.scc.kit.edu/fk9002/paraview_eikonal.git)).

# How does it work
- Create an input deck with a CAE preprocessing tool (e.g. Hypermesh) or use `blockMesh` for mesh generation / boundary definition (like in *tutorials/basic/invEikonalFoam/simplePlate*)
- Define the element-wise thickness **H** in *0/H*.
Use `setWallThickness` (as in *tutorials/basic/invEikonalFoam/DEmiL_PIWZ_2mm_quadElements*) or `setFields` (as in *tutorials/basic/invEikonalFoam/nonUniformPlate*) to define the local wall thickness **H**.
In case of a part with uniform thickness, define it directly in *0/H*.
- Define the boundary conditions of **G** in *0/G*.
Boundary @inlet(s): **G_0=2/l_ref**, whereas **l_ref** corresponds to the maximum part dimension. @top and bottom sides: if inlet / gate in plane `empty`, otherwise `zeroGradient`. @other: `zeroGradient`.
Also check the boundary conditions in *constant/poymesh/boundary*.
- Declare the smoothing factor **sigma** and reference height **h** in *constant/eikonalProperties*. **Tau** corresponds to a dummy variable for the pseudo time term.
- Define simulation controls in *system/controlDict*, numerical approximation schemes in *system/fvSchemes* and solver controls in *system/fvSolution* of your simulation job path.
- Navigate to the job directory via `cd `
- Activate the OpenFoam environment via `ofrtm6`
- Run the solver via `invEikonalFoam`
- Convert the results to .vtk via `foamToVTK`
- Postprocess (to obtain the distance field and fill time) via the [paraview_eikonal](https://git.scc.kit.edu/fk9002/paraview_eikonal.git) module for Paraview (requires Version 5.8 or higher).

# References
Fares and Schroeder. A differential equation for approximate wall distance.
International Journal for Numerical Methods in Fluids. 2002. 39(8):743-764.
