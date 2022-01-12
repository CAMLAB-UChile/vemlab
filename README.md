# VEMLAB: a MATLAB library for the virtual element method

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
<h2>First steps</h2>
<a>Once the software is downloaded there are two fonts that need to be installed so that the output figures are plotted correctly. If these fonts are not present, MATLAB will use default fonts. These fonts are Segoe UI Semibold and Good Times RG. For Windows systems, these fonts are included in folder “vemlab/utilities/fonts”. The instructions to install these fonts are also provided in the same folder. For completeness, the instruc-tions are also provided here.

Segoe UI Semibold: This font must be installed directly in your operating system. In Windows system this is done in Settings-->Fonts or just type "Font settings" in Window's search utility.

Good Times RG: This font must be added to the Java Runtime Environment (JRE) that ships with MATLAB (if that is the Java version being used by the program). Below are steps to install the fonts:

1.	Copy the TTF fonts to "<matlabroot>\sys\java\jre\win64\jre\lib\fonts", e.g., "C:\Program Files\MATLAB\R2019a\sys\java\jre\win64\jre\lib\fonts". Note: You might need to be an administrator to copy files to this folder.

2.	Restart MATLAB.

3.	You should now be able to see the fonts in MATLAB's font preference panel. Note: If using the default JVM shipped with MATLAB, you must do this for each MATLAB version installed.<a/>

<h2>About plotting capabilities</h2>
<a>Plots can be directly obtained through MATLAB’s figures by setting the following parame-ters in function plot_and_output_options.m that is located in folder “config”:

  create_matlab_contour_plots='yes'; % (for numerical solution)
  create_matlab_exact_contour_plots='yes'; % (for exact solution)
  
The figures can also be printed and saved to PDF files (the output PDFs are saved to folder “test/output_files/matlab_figures/”) by setting

  print_figures='yes'; % (for numerical solution)
  print_exact_figures='yes'; % (for exact solution)

The figures can also be saved to .fig files (the output .fig files are saved to folder “test/output_files/matlab_figures/”) by setting

  save_matlab_figures='yes'; % (for numerical solution) 
  save_exact_matlab_figures='yes'; % (for exact solution)   

The resolution of the MATLAB’s figures when plotting to PDF files or saving to .fig files is controlled by

  print_figures_resolution=600;

where it has been set to 600 dpi. This value can be changed as required by the user.
Many other customizations of the output figures can be made by setting the appropriate parameters in plot_and_output_options.m (see Section 8 of this primer for details). Sev-eral example output MATLAB’s figures are provided in Section 10 of this primer.

As an alternative, plots can be visualized externally in the GiD postprocessor. This is an independent process and can be performed even if plots were set to be displayed on MATLAB’s figures. When the following is set in plot_and_output_options.m:

  write_solutions_to_GiD_file='yes';

results are written to GiD files and saved to folder “test/output_files/GiD”. GiD can be downloaded from its webpage: https://www.gidhome.com/.<a/>
<h2>Author</h2>
<a href="https://github.com/aaortizb">Alejandro Ortiz-Bernardin</a>, Associate Professor, Department of Mechanical Engineering, Universidad de Chile.
<h2>Running VEMLab</h2>
<a>VEMLAB is a library. You need to create a main .m file and place it inside the folder “test.” The main file has the typical structure of a FEM simulation. Simply follow the test problems (they are given with detailed comments) that are provided inside the folder “test” to write your own .m files or modify the ones provided. Alternatively, you can read the manual that is available in the folder 'doc.'</a>
<h2>License</h2>
<a>This project is licensed under the GPL3 License. This program is free software; it can be redistributed or modified under the terms of the GNU General Public License 3 as published by the Free Software Foundation.<a/>
