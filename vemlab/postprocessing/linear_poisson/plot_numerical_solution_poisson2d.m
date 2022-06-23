%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                       VEMLab
%           Source code  : http://camlab.cl/research/software/vemlab/
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:              plot_numerical_solution_poisson2d 
%
% Created by : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
% Updated by :
% Updated by :
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% Plot scalar, flux and gradient field numerical solutions.
%
% Usage
% =====
% [triangles_per_polygon,flux,grad] = plot_numerical_solution_poisson2d(domainMesh,solution,...
%                                                                       matProps,config)
%
% Input
% =====
% domainMesh : structure containing mesh data (coords,connect,etc.)
% solution : nodal scalar field solution
% matProps : structure containing the material properties
% config : structure storing VEMLab configuration options and behavior
%
% Output
% ======
% triangles_per_polygon : array containing the number of triangles that 
%                         subdivide each polygon (VEM case only)
% flux : struct. (flux.x and flux.y) storing numerical flux at Gauss points
% grad : struct. (grad.x and grad.y) storing numerical gradient at Gauss points
%
%-------------------------------------------------------------------------------
% References 
% ==========
%
%-------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Jan. 04, 2021: plotting functionalities improved
% Feb. 1, 2020: add a return array variable called triangles_per_polygon, which
%               is used to fix an error in the plotting of VEM fluxes and gradients 
%               into a text file stage (by A. Ortiz-Bernardin)
% Jan. 31, 2020: add a check on matProps to figure out if comes on an
%                element-by-element fashion or not (VEM and FEM cases).
% Oct. 20, 2018: add option to switch off all matlab contour plots (by A. Ortiz-Bernardin)
% Apr. 19, 2018: improve the plotting of axis and fonts (by A. Ortiz-Bernardin)
% Mar. 17, 2018: first realease (by A. Ortiz-Bernardin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [flux,grad] = ...
               plot_numerical_solution_poisson2d(domainMesh,solution,matProps,config)
             
  [flux,grad]=...
    plot_numerical_solution_vem_fem_poisson2d(solution,domainMesh,matProps,config);

end

function [flux,grad]=plot_numerical_solution_vem_fem_poisson2d(solution,domainMesh,matProps,config)

  if strcmp(config.vemlab_method,'VEM2D')
    if strcmp(config.poisson2d_plot_flux_and_gradient,'yes')
      [flux,grad]=...
        vem_flux_and_gradient_poisson2d(solution,domainMesh,matProps,config);
    else
      flux=[]; grad=[];
    end      
    if strcmp(config.create_matlab_contour_plots,'yes')
      patch_plot_VEM2D_FEM2DT3_poisson2d(domainMesh,solution,flux,grad,config); 
    end      
  elseif strcmp(config.vemlab_method,'FEM2DT3')
    if strcmp(config.poisson2d_plot_flux_and_gradient,'yes')    
      [flux,grad]=...
        fem2DT3_flux_and_gradient_poisson2d(solution,domainMesh,matProps,config);
    else
      flux=[]; grad=[];
    end         
    if strcmp(config.create_matlab_contour_plots,'yes')
      patch_plot_VEM2D_FEM2DT3_poisson2d(domainMesh,solution,flux,grad,config);
    end   
  elseif strcmp(config.vemlab_method,'FEM2DQ4')   
    if strcmp(config.poisson2d_plot_flux_and_gradient,'yes')       
      [flux,grad,gp_list,h_min,xmin,xmax,ymin,ymax]=...
        fem2DQ4_flux_and_gradient_poisson2d(solution,domainMesh,matProps,config); 
    else
      flux=[]; grad=[];
    end        
    if strcmp(config.create_matlab_contour_plots,'yes')
      patch_plot_FEM2DQ4_poisson2d(domainMesh,solution,flux,grad,gp_list,...
                                   h_min,xmin,xmax,ymin,ymax,config)                                     
    end      
  else
    throw_error('Error in plot_numerical_solution_poisson2d.m --> plot_numerical_solution_vem_fem_poisson2d: vemlab_method\n');
  end

end

function [flux,grad]=vem_flux_and_gradient_poisson2d(solution,domainMesh,matProps,config) 
     
  fprintf('Postprocessing %s flux and gradient...\n',...
          config.vemlab_method); 
  % mesh data
  coords=domainMesh.coords;
  connect=domainMesh.connect;    
  num_elem=length(connect);  
  % figure out whether material's data is given in an element-by-element fashion or not
  size_k=length(matProps.k);
  if size_k==num_elem
    k_is_unique = false;
  elseif size_k==1
    k_is_unique = true;
  end  
  % loop over elements  
  for i=1:num_elem
    nodes=connect{i};
    elem_coord=zeros(length(nodes),2);
    % element dofs global indices
    dofs=nodes;
    % element VEM nodal solution
    uh_elem_column=solution(dofs);
    % elemen nodal coordinates
    for h=1:length(nodes)
      elem_coord(h,1)=coords(nodes(h),1);
      elem_coord(h,2)=coords(nodes(h),2);
    end     
    % VEM matrix
    verts=elem_coord;    
    Wc=vem_Wc_poisson2d(verts);
    % VEM gradient
    pic_grad_uh=Wc'*uh_elem_column;    
    grad.dx(i)=pic_grad_uh(1);
    grad.dy(i)=pic_grad_uh(2);   
    % conductivity for isotropic material 
    if k_is_unique    
      k=matProps.k; % conductivity is the same for all the elements
    else
      k=matProps.k{i,1}; % conductivity is particular for the current element
    end
    size_k_i=length(k);      
    % VEM flux 
    if size_k_i == 1
      flux.qx(i)=-k*pic_grad_uh(1);
      flux.qy(i)=-k*pic_grad_uh(2);
    else
      throw_error('In plot_numerical_solution_poisson2d.m --> plot_flux_and_gradient --> vem_flux_and_gradient_poisson2d : either conductivity was not defined or multiple conductivities assigned to the element... not implemented for this condition');
    end     
  end    
end

function [flux,grad]=fem2DT3_flux_and_gradient_poisson2d(solution,domainMesh,matProps,config) 
     
  fprintf('Postprocessing %s flux and gradient...\n',...
          config.vemlab_method); 
  % mesh data
  coords=domainMesh.coords;
  connect=domainMesh.connect;    
  num_elem=length(connect);  
  % figure out whether material's data is given in an element-by-element fashion or not
  size_k=length(matProps.k);
  if size_k==num_elem
    k_is_unique = false;
  elseif size_k==1
    k_is_unique = true;
  end  
  % loop over elements  
  for i=1:num_elem
    nodes=connect{i};
    elem_coord=coords(nodes,:);
    % element dofs global indices
    dofs=nodes;
    % element FEM nodal solution
    uh_elem_column=solution(dofs); 
    % area of the element
    area=triangle_area(elem_coord);
    % deformation matrix
    xx=elem_coord(:,1); yy=elem_coord(:,2);  
    B=(1/(2*area))*[yy(2)-yy(3),yy(3)-yy(1),yy(1)-yy(2);...
                    xx(3)-xx(2),xx(1)-xx(3),xx(2)-xx(1)];                     
    % FEM gradient
    grad_uh=B*uh_elem_column;
    grad.dx(i)=grad_uh(1);
    grad.dy(i)=grad_uh(2);  
    % conductivity for isotropic material 
    if k_is_unique    
      k=matProps.k; % conductivity is the same for all the elements
    else
      k=matProps.k{i,1}; % conductivity is particular for the current element
    end
    size_k_i=length(k);      
    % FEM flux 
    if size_k_i == 1
      flux.qx(i)=-k*grad_uh(1);
      flux.qy(i)=-k*grad_uh(2); 
    else
      throw_error('In plot_numerical_solution_poisson2d.m --> plot_flux_and_gradient --> fem2DT3_flux_and_gradient_poisson2d : either conductivity was not defined or multiple conductivities assigned to the element... not implemented for this condition');
    end   
  end    
end

function [flux,grad,gp_list,h_min,xmin,xmax,ymin,ymax] = ...
           fem2DQ4_flux_and_gradient_poisson2d(solution,domainMesh,matProps,config)
  
  fprintf('Postprocessing %s flux and gradient at Gauss points...\n',...
          config.vemlab_method);
  % The bounding box must contain the domain. If the domain is a rectangle, the
  % bounding box coincides with the domain
  xmin=domainMesh.BdBox(1); ymin=domainMesh.BdBox(3);
  xmax=domainMesh.BdBox(2); ymax=domainMesh.BdBox(4);
  dx=abs(xmax-xmin); dy=abs(ymax-ymin);
  h_min=max(dx,dy);    
  % mesh data
  coords=domainMesh.coords;
  connect=domainMesh.connect;    
  num_elem=length(connect); 
  % figure out whether material's data is given in an element-by-element fashion or not
  size_k=length(matProps.k);
  if size_k==num_elem
    k_is_unique = false;
  elseif size_k==1
    k_is_unique = true;
  end
  
  % loop over elements
  gp=0;    
  for i=1:num_elem
    nodes=connect{i};
    elem_coord=coords(nodes,:);
    % element dofs global indices
    dofs=nodes;
    % element FEM nodal solution
    uh_elem_column=solution(dofs); 
    % compute flux and gradient at Gauss points
    int_order=2;
    xi=gauss_points_1d(int_order);
    eta=gauss_points_1d(int_order);
    for gpxi=1:length(xi)
      dN1deta=-(1-xi(gpxi))/4;
      dN2deta=-(1+xi(gpxi))/4;
      dN3deta=(1+xi(gpxi))/4;
      dN4deta=(1-xi(gpxi))/4;      
      for gpeta=1:length(eta)
        gp=gp+1;
        dN1dxi=-(1-eta(gpeta))/4;
        dN2dxi=(1-eta(gpeta))/4;
        dN3dxi=(1+eta(gpeta))/4;
        dN4dxi=-(1+eta(gpeta))/4;
        dNdxi=[dN1dxi;dN2dxi;dN3dxi;dN4dxi];
        dNdeta=[dN1deta;dN2deta;dN3deta;dN4deta];
        xcoord=elem_coord(:,1); ycoord=elem_coord(:,2);
        dxdxi=dot(dNdxi,xcoord);
        dydxi=dot(dNdxi,ycoord);
        dxdeta=dot(dNdeta,xcoord);
        dydeta=dot(dNdeta,ycoord);    
        detJ=dxdxi*dydeta-dydxi*dxdeta;    
        dN1dx=(dydeta*dN1dxi-dydxi*dN1deta)/detJ;
        dN1dy=(dxdxi*dN1deta-dxdeta*dN1dxi)/detJ;
        dN2dx=(dydeta*dN2dxi-dydxi*dN2deta)/detJ;
        dN2dy=(dxdxi*dN2deta-dxdeta*dN2dxi)/detJ;    
        dN3dx=(dydeta*dN3dxi-dydxi*dN3deta)/detJ;
        dN3dy=(dxdxi*dN3deta-dxdeta*dN3dxi)/detJ;  
        dN4dx=(dydeta*dN4dxi-dydxi*dN4deta)/detJ;
        dN4dy=(dxdxi*dN4deta-dxdeta*dN4dxi)/detJ; 
        % deformation matrix
        B=[dN1dx,dN2dx,dN3dx,dN4dx;...
           dN1dy,dN2dy,dN3dy,dN4dy];           
        % FEM strains
        grad_uh=B*uh_elem_column; 
        grad.dx(gp)=grad_uh(1);
        grad.dy(gp)=grad_uh(2);
        % conductivity for isotropic material 
        if k_is_unique    
          k=matProps.k; % conductivity is the same for all the elements
        else
          k=matProps.k{i,1}; % conductivity is particular for the current element
        end
        size_k_i=length(k);      
        % FEM flux 
        if size_k_i == 1
          flux.qx(gp)=-k*grad_uh(1);
          flux.qy(gp)=-k*grad_uh(2);
        else
          throw_error('In plot_numerical_solution_poisson2d.m --> plot_flux_and_gradient --> fem_flux_and_gradient_poisson2d : either conductivity was not defined or multiple conductivities assigned to the element... not implemented for this condition');
        end  
        % assign flux and gradient to the element Gauss points             
        N1=(1-xi(gpxi))*(1-eta(gpeta))/4;
        N2=(1+xi(gpxi))*(1-eta(gpeta))/4;
        N3=(1+xi(gpxi))*(1+eta(gpeta))/4;
        N4=(1-xi(gpxi))*(1+eta(gpeta))/4;
        Nv=[N1; N2; N3; N4];
        xcoord=dot(Nv,elem_coord(:,1));
        ycoord=dot(Nv,elem_coord(:,2));
        xy=[xcoord,ycoord]; % int. point in the form [ x(xi,eta) y(xi,eta) ]     
        gp_list.x(gp)=xy(1);
        gp_list.y(gp)=xy(2);          
      end
    end
    % h_min
    %h_size=max_edge_size(elem_coord);
    h_size=polyarea(elem_coord(:,1),elem_coord(:,2));        
    if h_size<h_min
      h_min=h_size;
    end        
  end
  h_min=sqrt(h_min); % square root of the minumum element's area among the elements of the mesh  
end

function patch_plot_VEM2D_FEM2DT3_poisson2d(domainMesh,solution,flux,grad,config)  

  solutionType=config.vemlab_method;
  nodes=domainMesh.coords; 
  polygons=domainMesh.connect;

  fprintf('\n'); 
  fprintf('Plotting %s u field solution...\n',solutionType);  
  
  if strcmp(config.poisson2d_plot_scalar_field.u,'yes')
    mytitleclb='$u_h$';   
    mytitlefigfile='uh'; 
    colorbar_limits=config.poisson2d_plot_scalar_field.clim.u;     
    patch_plot_figure(config,'on_nodes',polygons,nodes,solution,mytitleclb,mytitlefigfile,colorbar_limits,0.65,'numerical');
  end
    
  if ~isempty(flux)
    fprintf('Plotting %s flux and gradient solutions...\n',solutionType);   

    if strcmp(config.poisson2d_plot_flux.qx,'yes')      
      mytitleclb='$q_{1,h}$';
      mytitlefigfile='q1h'; 
      colorbar_limits=config.poisson2d_plot_flux.clim.qx;         
      patch_plot_figure(config,'on_faces',polygons,nodes,flux.qx',mytitleclb,mytitlefigfile,colorbar_limits,0.96,'numerical');
    end

    if strcmp(config.poisson2d_plot_flux.qy,'yes')  
      mytitleclb='$q_{2,h}$';
      mytitlefigfile='q2h';
      colorbar_limits=config.poisson2d_plot_flux.clim.qy;       
      patch_plot_figure(config,'on_faces',polygons,nodes,flux.qy',mytitleclb,mytitlefigfile,colorbar_limits,0.96,'numerical');
    end 

    if strcmp(config.poisson2d_plot_flux.qnorm,'yes') 
      fluxNorm=sqrt((flux.qx).*(flux.qx)+(flux.qy).*(flux.qy)); 
      mytitleclb='$\scriptstyle{\|}\displaystyle{\mathbf{q}_h}\scriptstyle{\|}$';
      mytitlefigfile='qh'; 
      colorbar_limits=config.poisson2d_plot_flux.clim.qnorm;       
      patch_plot_figure(config,'on_faces',polygons,nodes,fluxNorm',mytitleclb,mytitlefigfile,colorbar_limits,1.06,'numerical');
    end       

    if strcmp(config.poisson2d_plot_grad.dx,'yes')  
%       mytitle='$\frac{\partial u_h}{\partial x_1}$';
      mytitleclb='$\scriptstyle{\partial} \displaystyle{u_h} \scriptstyle{/ \partial} \displaystyle{x_1}$';    
      mytitlefigfile='dudxh'; 
      colorbar_limits=config.poisson2d_plot_grad.clim.dx;        
      patch_plot_figure(config,'on_faces',polygons,nodes,grad.dx',mytitleclb,mytitlefigfile,colorbar_limits,2.14,'numerical');
    end  

    if strcmp(config.poisson2d_plot_grad.dy,'yes') 
      mytitleclb='$\scriptstyle{\partial} \displaystyle{u_h} \scriptstyle{/ \partial} \displaystyle{x_2}$';
      mytitlefigfile='dudyh'; 
      colorbar_limits=config.poisson2d_plot_grad.clim.dy;       
      patch_plot_figure(config,'on_faces',polygons,nodes,grad.dy',mytitleclb,mytitlefigfile,colorbar_limits,2.14,'numerical');
    end   

    if strcmp(config.poisson2d_plot_grad.dnorm,'yes')
      gradNorm=sqrt((grad.dx).*(grad.dx)+(grad.dy).*(grad.dy));
      mytitleclb='$\scriptstyle{\|}\nabla \displaystyle{u_h}\scriptstyle{\|}$';
      mytitlefigfile='graduh';   
      colorbar_limits=config.poisson2d_plot_grad.clim.dnorm;     
      patch_plot_figure(config,'on_faces',polygons,nodes,gradNorm',mytitleclb,mytitlefigfile,colorbar_limits,1.47,'numerical');
    end 

  end

end

function patch_plot_FEM2DQ4_poisson2d(domainMesh,solution,flux,grad,gp_list,...
                                       h_min,xmin,xmax,ymin,ymax,config) 
                                         
  solutionType=config.vemlab_method;
  nodes=domainMesh.coords;
  polygons=domainMesh.connect;

  fprintf('\n'); 
  fprintf('Plotting %s u field solution...\n',solutionType);  

  if strcmp(config.poisson2d_plot_scalar_field.u,'yes')
    mytitleclb='$u_h$';  
    mytitlefigfile='uh';   
    colorbar_limits=config.poisson2d_plot_scalar_field.clim.u;     
    patch_plot_figure(config,'on_nodes',polygons,nodes,solution,mytitleclb,mytitlefigfile,colorbar_limits,0.65,'numerical');
  end
         
  if ~isempty(flux)  
    mt=0.25; %0.1;
    dx=mt*h_min;
    dy=dx;  
    xs=xmin:dx:xmax;
    ys=ymin:dy:ymax;         
    [xq,yq]=meshgrid(xs,ys); 

    fprintf('Plotting %s flux and gradient solutions...\n',solutionType);
    
    if strcmp(config.poisson2d_plot_flux.qx,'yes')
      mytitleclb='$q_{1,h}$';
      mytitlefigfile='q1h';  
      colorbar_limits=config.poisson2d_plot_flux.clim.qx;      
      mesh_plot_figure_fem2dQ4(config,polygons,nodes,gp_list,xq,yq,flux.qx',mytitleclb,mytitlefigfile,colorbar_limits,0.96,'numerical');
    end
    
    if strcmp(config.poisson2d_plot_flux.qy,'yes') 
      mytitleclb='$q_{2,h}$';
      mytitlefigfile='q2h';
      colorbar_limits=config.poisson2d_plot_flux.clim.qy;      
      mesh_plot_figure_fem2dQ4(config,polygons,nodes,gp_list,xq,yq,flux.qy',mytitleclb,mytitlefigfile,colorbar_limits,0.96,'numerical');
    end 
    
    if strcmp(config.poisson2d_plot_flux.qnorm,'yes')  
      fluxNorm=sqrt((flux.qx).*(flux.qx)+(flux.qy).*(flux.qy)); 
      mytitleclb='$\scriptstyle{\|}\displaystyle{\mathbf{q}_h}\scriptstyle{\|}$';
      mytitlefigfile='qh'; 
      colorbar_limits=config.poisson2d_plot_flux.clim.qnorm;      
      mesh_plot_figure_fem2dQ4(config,polygons,nodes,gp_list,xq,yq,fluxNorm',mytitleclb,mytitlefigfile,colorbar_limits,1.06,'numerical');
    end       
    
    if strcmp(config.poisson2d_plot_grad.dx,'yes')  
      mytitleclb='$\scriptstyle{\partial} \displaystyle{u_h} \scriptstyle{/ \partial} \displaystyle{x_1}$';
      mytitlefigfile='dudxh';  
      colorbar_limits=config.poisson2d_plot_grad.clim.dx;      
      mesh_plot_figure_fem2dQ4(config,polygons,nodes,gp_list,xq,yq,grad.dx',mytitleclb,mytitlefigfile,colorbar_limits,2.14,'numerical');
    end  
    
    if strcmp(config.poisson2d_plot_grad.dy,'yes') 
      mytitleclb='$\scriptstyle{\partial} \displaystyle{u_h} \scriptstyle{/ \partial} \displaystyle{x_2}$';
      mytitlefigfile='dudyh'; 
      colorbar_limits=config.poisson2d_plot_grad.clim.dy;       
      mesh_plot_figure_fem2dQ4(config,polygons,nodes,gp_list,xq,yq,grad.dy',mytitleclb,mytitlefigfile,colorbar_limits,2.14,'numerical');
    end   
    
    if strcmp(config.poisson2d_plot_grad.dnorm,'yes')  
      gradNorm=sqrt((grad.dx).*(grad.dx)+(grad.dy).*(grad.dy));
      mytitleclb='$\scriptstyle{\|}\nabla \displaystyle{u_h}\scriptstyle{\|}$';
      mytitlefigfile='graduh';
      colorbar_limits=config.poisson2d_plot_grad.clim.dnorm;       
      mesh_plot_figure_fem2dQ4(config,polygons,nodes,gp_list,xq,yq,gradNorm',mytitleclb,mytitlefigfile,colorbar_limits,1.47,'numerical');    
    end 
    
  end

end
