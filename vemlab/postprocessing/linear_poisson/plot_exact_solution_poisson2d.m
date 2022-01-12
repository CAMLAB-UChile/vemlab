%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                     VEMLab
%           Source code  : http://camlab.cl/research/software/vemlab/
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:              plot_exact_solution_poisson2d 
%
% Created by : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
% Updated by :
% Updated by :
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% Plot scalar, flux and gradient field exact solutions.
%
% Usage
% =====
% [flux,grad] = plot_exact_solution_poisson2d(domainMesh,exact_solution_handle,...
%                                             matProps,config)
%
% Input
% =====
% domainMesh : structure containing mesh data (coords,connect,etc.)
% exact_solution_handle : handle to exact solution functions
% matProps : structure containing the material properties
% config : structure storing VEMLab configuration options and behavior
%
% Output
% ======
%
%-------------------------------------------------------------------------------
% References 
% ==========
%
%-------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Jan. 04, 2021: plotting functionalities improved
% Oct. 20, 2018: add option to switch off all matlab contour plots (by A. Ortiz-Bernardin)
% Apr. 19, 2018: improve the plotting of axis and fonts (by A. Ortiz-Bernardin)
% Mar. 17, 2018: first realease (by A. Ortiz-Bernardin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_exact_solution_poisson2d(domainMesh,exact_solution_handle,...
                                       matProps,config)
  if strcmp(config.create_matlab_exact_contour_plots,'yes')
    [exact_nodal_solution,exact_nodal_flux,exact_nodal_grad]=...
      exact_solutions_poisson2d(exact_solution_handle,domainMesh.coords,matProps);
    patch_plot_exact_nodal_solutions_poisson2d(domainMesh,exact_nodal_solution,...
                                     exact_nodal_flux,exact_nodal_grad,config);
  end
    
end

function patch_plot_exact_nodal_solutions_poisson2d(domainMesh,solution,flux,grad,config)  

  nodes=domainMesh.coords; 
  polygons=domainMesh.connect;

  fprintf('\n'); 
  fprintf('Plotting exact u field solution...\n');  
  
  if strcmp(config.poisson2d_plot_scalar_field.u,'yes')
    mytitleclb='$u$'; 
    mytitlefigfile='u'; 
    colorbar_limits=config.poisson2d_plot_scalar_field.clim.u;      
    patch_plot_figure(config,'on_nodes',polygons,nodes,solution,mytitleclb,mytitlefigfile,colorbar_limits,0.40,'exact');
  end
    
  if strcmp(config.poisson2d_plot_flux_and_gradient,'yes')
    fprintf('Plotting exact flux and gradient solutions...\n');   

    if strcmp(config.poisson2d_plot_flux.qx,'yes')
      mytitleclb='$q_{1}$';
      mytitlefigfile='q1'; 
      colorbar_limits=config.poisson2d_plot_flux.clim.qx;       
      patch_plot_figure(config,'on_nodes',polygons,nodes,flux.x,mytitleclb,mytitlefigfile,colorbar_limits,0.58,'exact');
    end

    if strcmp(config.poisson2d_plot_flux.qy,'yes')  
      mytitleclb='$q_{2}$';
      mytitlefigfile='q2'; 
      colorbar_limits=config.poisson2d_plot_flux.clim.qy;       
      patch_plot_figure(config,'on_nodes',polygons,nodes,flux.y,mytitleclb,mytitlefigfile,colorbar_limits,0.58,'exact');
    end 

    if strcmp(config.poisson2d_plot_flux.qnorm,'yes') 
      fluxNorm=sqrt((flux.x).*(flux.x)+(flux.y).*(flux.y)); 
      mytitleclb='$\scriptstyle{\|}\displaystyle{\mathbf{q}}\scriptstyle{\|}$';
      mytitlefigfile='q'; 
      colorbar_limits=config.poisson2d_plot_flux.clim.qnorm;       
      patch_plot_figure(config,'on_nodes',polygons,nodes,fluxNorm,mytitleclb,mytitlefigfile,colorbar_limits,0.80,'exact');
    end       

    if strcmp(config.poisson2d_plot_grad.dx,'yes')  
      mytitleclb='$\scriptstyle{\partial} \displaystyle{u} \scriptstyle{/ \partial} \displaystyle{x_1}$';
      mytitlefigfile='dudx';
      colorbar_limits=config.poisson2d_plot_grad.clim.dx;       
      patch_plot_figure(config,'on_nodes',polygons,nodes,grad.x,mytitleclb,mytitlefigfile,colorbar_limits,1.82,'exact');
    end  

    if strcmp(config.poisson2d_plot_grad.dy,'yes') 
      mytitleclb='$\scriptstyle{\partial} \displaystyle{u} \scriptstyle{/ \partial} \displaystyle{x_2}$';
      mytitlefigfile='dudy'; 
      colorbar_limits=config.poisson2d_plot_grad.clim.dy;       
      patch_plot_figure(config,'on_nodes',polygons,nodes,grad.y,mytitleclb,mytitlefigfile,colorbar_limits,1.82,'exact');
    end   

    if strcmp(config.poisson2d_plot_grad.dnorm,'yes')
      gradNorm=sqrt((grad.x).*(grad.x)+(grad.y).*(grad.y));
      mytitleclb='$\scriptstyle{\|}\nabla \displaystyle{u}\scriptstyle{\|}$';
      mytitlefigfile='gradu'; 
      colorbar_limits=config.poisson2d_plot_grad.clim.dnorm;       
      patch_plot_figure(config,'on_nodes',polygons,nodes,gradNorm,mytitleclb,mytitlefigfile,colorbar_limits,1.17,'exact');      
    end 

  end

end



