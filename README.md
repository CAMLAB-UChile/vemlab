# VEMLab: a MATLAB library for the virtual element method

This repository contains the code for an open source MATLAB library for the virtual element method.
<h2>Features</h2>
<ul><li> Two-dimensional linear elastostatics (plane strain and plane stress) and two-dimensional Poisson problem.</li>
    <li> Solution methods: Linear VEM (polygonal elements), FEM (3-node triangles, 4-node quadrilateral).</li>
    <li> Boundary conditions: Dirichlet, Neumann on boundary edges; can be a constant or a function.</li>  
    <li> Meshers: PolyMesher, distmesh2d, quad4mesh; PolyMesher is customized for rectangular domain, wrench domain and plate with a hole domain; distmesh2d and quad4mesh are customized for rectangular domain only. Domains can be extended for any of the meshers, but it requires adjustments to some interface functions (see the instructions that are available in functions create_polygonal_mesh.m, create_quadrilateral_mesh.m and create_triangular_mesh.m in folder “mesher”).</li>  
    <li> Meshes need to be generated separately and saved to folder “test/mesh_files.”</li>
    <li> Meshes must be generated with the functions “create_” located in the folder “mesher.” Then, the files containing the generated meshes will be automatically saved to folder “test/mesh_files” for their use.</li> 
    <li> Solutions can be plotted to MATLAB figures, text files, <a href="https://www.gidhome.com/">GiD</a> files and <a href="https://www.vtk.org/">VTK</a> files.</li>  
</ul>
<h2>About plotting capabilities</h2>
<a>Currently, the MATLAB plotting of stresses and strains (linear elastostatic problem) or fluxes and gradients (Poisson problem) is very limited and won't recognize holes that might come with the domain geometries. It will work ok with rectangular geometries without holes, but the mesh will need to be very refined to obtain a colormap on the whole geometry ... just try a coarse mesh and then a refined mesh to see what is being said.
If quality colormap plots are required for stresses, strains, fluxes and gradients, it is highly recommended writing the GiD output files and visualizing them in "GiD the pre and postprocessor" (www.gidhome.com).<a/>
<a>On the other hand, the MATLAB plotting of primary field variables (like displacements) presents no limitations and produces quality colormaps.<a/>
<h2>Author</h2>
<a href="https://github.com/aaortizb">Alejandro Ortiz-Bernardin</a>, Assistant Professor, Department of Mechanical Engineering, Universidad de Chile.
<h2>Running VEMLab</h2>
<a>VEMLab is a library. You need to create a main .m file and place it inside the folder “test.” The main file has the typical structure of a FEM simulation. Simply follow the test problems (they are given with detailed comments) that are provided inside the folder “test” to write your own .m files or modify the ones provided. Alternatively, you can read the manual that is available in the folder 'doc.'</a>
<h2>License</h2>
<a>This project is licensed under the GPL3 License. This program is free software; it can be redistributed or modified under the terms of the GNU General Public License 3 as published by the Free Software Foundation.<a/>
