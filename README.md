# VEMLab: a MATLAB library for the virtual element method

This repository contains the code for an open source MATLAB library for the virtual element method.
<h2>Features</h2>
<ul><li> Two-dimensional linear elastostatics (plane strain and plane stress).</li>
    <li> Solution methods: VEM (polygonal elements), FEM (3-node triangles, 4-node quadrilateral).</li>
    <li> Boundary conditions: Dirichlet, Neumann on boundary edges; can be a constant or a function.</li>  
    <li> Meshers: PolyMesher, distmesh2d, quad4mesh; customized for rectangular domains only (requires adjustments for other domain types).</li>  
    <li> Meshes need to be generated separately and stored inside folder 'mesh_files' located in  the folder 'test.'</li>
    <li> Meshes must be generated with the functions 'create_' located in the folder 'mesher.'</li> 
    <li> Solutions can be plotted to MATLAB figures, text files and <a href="https://www.gidhome.com/">GiD</a> files.</li>  
</ul>
<h2>Author</h2><a href="https://github.com/aaortizb">Alejandro Ortiz-Bernardin</a>
<h2>Running VEMLab</h2><a>VEMLab is a library. You need to create a main .m file and place it inside the 'test' folder. The main file has the typical structure of a FEM simulation. Simply follow the test problems that are provided inside the folder 'test' to wrote your own .m files or modify the ones provided.</a>
<h2>License</h2><a>This project is licensed under the GPL3 License. This program is free software; it can be redistributed or modified under the terms of the GNU General Public License 3 as published by the Free Software Foundation.<a/>
