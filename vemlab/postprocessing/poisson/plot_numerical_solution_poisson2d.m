%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                     VEMLab
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
% [flux,grad] = plot_numerical_solution_poisson2d(domainMesh,solution,...
%                                                 matProps,config)
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
% Oct. 20, 2018: add option to switch off all matlab contour plots (by A. Ortiz-Bernardin)
% Apr. 19, 2018: improve the plotting of axis and fonts (by A. Ortiz-Bernardin)
% Mar. 17, 2018: first realease (by A. Ortiz-Bernardin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [flux,grad] = plot_numerical_solution_poisson2d(domainMesh,solution,...
                                                         matProps,config)
  % plot u field                                       
  plot_u_field(domainMesh,solution,config);
  % plot flux and gradient
  [flux,grad]=plot_flux_and_gradient(solution,domainMesh,matProps,config);
end


function plot_u_field(domainMesh,solution,config)
  if strcmp(config.create_matlab_contour_plots,'yes')
    solutionType=config.vemlab_method;
    plotMeshOverResults=config.plot_mesh_over_results;  

    fprintf('\n');
    fprintf('Plotting %s solution...\n',solutionType);   

    nodes=domainMesh.coords; 
    polygons=domainMesh.connect;
    titleSolution='$u^h$';  

    if strcmp(config.poisson2d_plot_scalar_field.u,'yes')
      figure; 
      title(titleSolution,'FontWeight','bold','FontSize',20,'FontName',...
            'Times New Roman','Interpreter','latex');
      maxNumVertices = max(cellfun(@numel,polygons));
      padFunc = @(vertList) [vertList' NaN(1,maxNumVertices-numel(vertList))];
      elements = cellfun(padFunc,polygons,'UniformOutput',false);
      elements = vertcat(elements{:});
      data = [nodes,solution];
      if strcmp(plotMeshOverResults,'yes')
        patch('Faces',elements,'Vertices',data,...
              'FaceColor','interp','CData',solution);  
      else
        patch('Faces',elements,'Vertices',data,...
              'FaceColor','interp','EdgeColor','interp','CData',solution);
      end  
    %   axis('square') 
      axis equal;  
      dx = 0; dy = 0; dz = 0;
      xlim([min(nodes(:, 1)) - dx, max(nodes(:, 1)) + dx])
      ylim([min(nodes(:, 2)) - dy, max(nodes(:, 2)) + dy]) 
      if min(solution)~=max(solution)
        zlim([min(solution) - dz, max(solution) + dz])
      end
      xlabel('$x$','FontWeight','bold','FontSize',20,'FontName',...
             'Times New Roman','Interpreter','latex');
      ylabel('$y$','FontWeight','bold','FontSize',20,'FontName',...
             'Times New Roman','Interpreter','latex');
      zlabel(titleSolution,'FontWeight','bold','FontSize',20,'FontName',...
             'Times New Roman','Interpreter','latex');
      colorbar('FontName','Times New Roman','FontSize',14,'FontWeight','bold');
      colormap jet
      %set(gcf,'Renderer','painters')    
      set(gcf,'InvertHardcopy','off','Color',[1 1 1])
      set(gca,'FontName','Times New Roman','FontSize',14,'FontWeight','bold');
    end
  end
end

function [flux,grad]=plot_flux_and_gradient(solution,domainMesh,matProps,config)
  if strcmp(config.poisson2d_plot_flux_and_gradient,'yes')
    if strcmp(config.vemlab_method,'VEM2D')
      [flux,grad,gp_list,h_min,xmin,xmax,ymin,ymax]=...
        vem_flux_and_gradient_poisson2d(solution,domainMesh,matProps,config);
    elseif strcmp(config.vemlab_method,'FEM2DT3')||strcmp(config.vemlab_method,'FEM2DQ4')
      [flux,grad,gp_list,h_min,xmin,xmax,ymin,ymax]=...
        fem_flux_and_gradient_poisson2d(solution,domainMesh,matProps,config);     
    else
      throw_error('Error in plot_numerical_solution_poisson2d.m --> plot_flux_and_gradient: vemlab_method\n');
    end
    if strcmp(config.create_matlab_contour_plots,'yes')
      plot_flux_and_gradient_poisson2d(domainMesh,flux,grad,gp_list,h_min,...
                                       xmin,xmax,ymin,ymax,config);
    end
  else
    flux=[]; grad=[];
  end
end

function [flux,grad,gp_list,h_min,xmin,xmax,ymin,ymax] = ...
       vem_flux_and_gradient_poisson2d(solution,domainMesh,matProps,config) 
     
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
  % loop over elements
  gp=0;   
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
    grad_uh.x(i)=pic_grad_uh(1);
    grad_uh.y(i)=pic_grad_uh(2);     
    % VEM flux
    k=matProps.k;  % isotropic material   
    flux_h.x(i)=-k*pic_grad_uh(1);
    flux_h.y(i)=-k*pic_grad_uh(2);    
    % triangulate and assign the constant flux and gradient to the subtriangles
    newconnect=triangulate_polygon(domainMesh,i);
    for tr_i=1:size(newconnect,1)
      gp=gp+1;   
      % assign constant flux and gradient to the subtriangle Gauss point (1-pt rule)
      subtriangle_coords=domainMesh.coords(newconnect(tr_i,:),:);
      [xy,~]=gauss_points_T3(1,subtriangle_coords); % int. point in the form [x y]       
      gp_list.x(gp)=xy(1);
      gp_list.y(gp)=xy(2);        
      % exact solution
      grad.dx(gp)=grad_uh.x(i);
      grad.dy(gp)=grad_uh.y(i);  
      flux.qx(gp)=flux_h.x(i);
      flux.qy(gp)=flux_h.y(i);          
      % minimum size of the subtriangular mesh
      nodes_indices_T3=newconnect(tr_i,1:3);
      verts=domainMesh.coords(nodes_indices_T3,:);
      h_size=max_edge_size(verts);
      if h_size<h_min
        h_min=h_size;
      end 
    end      
  end    
end

function [flux,grad,gp_list,h_min,xmin,xmax,ymin,ymax] = ...
           fem_flux_and_gradient_poisson2d(solution,domainMesh,matProps,config)
  
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
  if strcmp(config.vemlab_method,'FEM2DT3')
    % loop over elements
    gp=0;
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
      % FEM flux
      k=matProps.k;  % isotropic material   
      flux.qx(i)=-k*grad_uh(1);
      flux.qy(i)=-k*grad_uh(2);  
      % assign constant flux and gradient to the element Gauss point (1-pt rule)
      gp=gp+1;
      [xy,~]=gauss_points_T3(1,elem_coord); % int. point in the form [x y]         
      gp_list.x(gp)=xy(1);
      gp_list.y(gp)=xy(2);  
      % h_min
      h_size=max_edge_size(elem_coord);
      if h_size<h_min
        h_min=h_size;
      end
    end
  elseif strcmp(config.vemlab_method,'FEM2DQ4')
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
          % VEM flux
          k=matProps.k;  % isotropic material   
          flux.qx(gp)=-k*grad_uh(1);
          flux.qy(gp)=-k*grad_uh(2);    
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
      h_size=max_edge_size(elem_coord);
      if h_size<h_min
        h_min=h_size;
      end        
    end
  else
    throw_error('Error in plot_numerical_solution_poisson2d.m --> plot_flux_and_gradient --> fem_flux_and_gradient_poisson2d: vemlab_method');
  end
end

function plot_flux_and_gradient_poisson2d(domainMesh,flux,grad,gp_list,...
                                          h_min,xmin,xmax,ymin,ymax,config)  
  nodes=domainMesh.coords;  
  [xq,yq]=meshgrid(xmin:0.05*h_min:xmax,ymin:0.05*h_min:ymax);  
  fluxNorm=sqrt((flux.qx).*(flux.qx)+(flux.qy).*(flux.qy)); 
  gradNorm=sqrt((grad.dx).*(grad.dx)+(grad.dy).*(grad.dy));   
  
  if strcmp(config.poisson2d_plot_flux.qx,'yes')
    vq=griddata(gp_list.x(:),gp_list.y(:),flux.qx(:),xq,yq);
    figure
    mesh(xq,yq,vq,'FaceColor','interp','EdgeColor','interp');
    mytitle='$q_x^h$';
    title(mytitle,'FontWeight','bold','FontSize',20,'FontName',...
          'Times New Roman','Interpreter','latex');  
    view(2);
    grid off;
    axis equal;  
    dx = 0; dy = 0; dz = 0;
    xlim([min(nodes(:, 1)) - dx, max(nodes(:, 1)) + dx])
    ylim([min(nodes(:, 2)) - dy, max(nodes(:, 2)) + dy])  
    if min(flux.qx)~=max(flux.qx)
      zlim([min(flux.qx) - dz, max(flux.qx) + dz])
    end
    xlabel('$x$','FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    ylabel('$y$','FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    zlabel(mytitle,'FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    colorbar('FontName','Times New Roman','FontSize',14,'FontWeight','bold');
    colormap jet
    %set(gcf,'Renderer','painters')    
    set(gcf,'InvertHardcopy','off','Color',[1 1 1])
    set(gca,'FontName','Times New Roman','FontSize',14,'FontWeight','bold');
  end

  if strcmp(config.poisson2d_plot_flux.qy,'yes')  
    vq=griddata(gp_list.x(:),gp_list.y(:),flux.qy(:),xq,yq);
    figure
    mesh(xq,yq,vq,'FaceColor','interp','EdgeColor','interp');
    mytitle='$q_y^h$';
    title(mytitle,'FontWeight','bold','FontSize',20,'FontName',...
          'Times New Roman','Interpreter','latex');   
    view(2);
    grid off;
    axis equal;  
    dx = 0; dy = 0; dz = 0;
    xlim([min(nodes(:, 1)) - dx, max(nodes(:, 1)) + dx])
    ylim([min(nodes(:, 2)) - dy, max(nodes(:, 2)) + dy])  
    if min(flux.qy)~=max(flux.qy)
      zlim([min(flux.qy) - dz, max(flux.qy) + dz])
    end
    xlabel('$x$','FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    ylabel('$y$','FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    zlabel(mytitle,'FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    colorbar('FontName','Times New Roman','FontSize',14,'FontWeight','bold');
    colormap jet
    %set(gcf,'Renderer','painters')    
    set(gcf,'InvertHardcopy','off','Color',[1 1 1])
    set(gca,'FontName','Times New Roman','FontSize',14,'FontWeight','bold');
  end

  if strcmp(config.poisson2d_plot_flux.qnorm,'yes')  
    vq=griddata(gp_list.x(:),gp_list.y(:),fluxNorm,xq,yq);
    figure
    mesh(xq,yq,vq,'FaceColor','interp','EdgeColor','interp');
    mytitle='$||q^h||$';
    title(mytitle,'FontWeight','bold','FontSize',20,'FontName',...
          'Times New Roman','Interpreter','latex');   
    view(2);
    grid off;
    axis equal;  
    dx = 0; dy = 0; dz = 0;
    xlim([min(nodes(:, 1)) - dx, max(nodes(:, 1)) + dx])
    ylim([min(nodes(:, 2)) - dy, max(nodes(:, 2)) + dy])  
    if min(fluxNorm)~=max(fluxNorm)
      zlim([min(fluxNorm) - dz, max(fluxNorm) + dz])
    end
    xlabel('$x$','FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    ylabel('$y$','FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    zlabel(mytitle,'FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    colorbar('FontName','Times New Roman','FontSize',14,'FontWeight','bold');
    colormap jet
    %set(gcf,'Renderer','painters')    
    set(gcf,'InvertHardcopy','off','Color',[1 1 1])
    set(gca,'FontName','Times New Roman','FontSize',14,'FontWeight','bold');
  end

  if strcmp(config.poisson2d_plot_grad.dx,'yes')  
    vq=griddata(gp_list.x(:),gp_list.y(:),grad.dx(:),xq,yq);
    figure
    mesh(xq,yq,vq,'FaceColor','interp','EdgeColor','interp');
    mytitle='$\frac{\partial u^h}{\partial x}$';
    title(mytitle,'FontWeight','bold','FontSize',20,'FontName',...
          'Times New Roman','Interpreter','latex');   
    view(2);
    grid off;
    axis equal;  
    dx = 0; dy = 0; dz = 0;
    xlim([min(nodes(:, 1)) - dx, max(nodes(:, 1)) + dx])
    ylim([min(nodes(:, 2)) - dy, max(nodes(:, 2)) + dy])  
    if min(grad.dx)~=max(grad.dx)
      zlim([min(grad.dx) - dz, max(grad.dx) + dz])
    end
    xlabel('$x$','FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    ylabel('$y$','FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    zlabel(mytitle,'FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    colorbar('FontName','Times New Roman','FontSize',14,'FontWeight','bold');
    colormap jet
    %set(gcf,'Renderer','painters')    
    set(gcf,'InvertHardcopy','off','Color',[1 1 1])
    set(gca,'FontName','Times New Roman','FontSize',14,'FontWeight','bold'); 
  end

  if strcmp(config.poisson2d_plot_grad.dy,'yes')  
    vq=griddata(gp_list.x(:),gp_list.y(:),grad.dy(:),xq,yq);
    figure
    mesh(xq,yq,vq,'FaceColor','interp','EdgeColor','interp');
    mytitle='$\frac{\partial u^h}{\partial y}$';
    title(mytitle,'FontWeight','bold','FontSize',20,'FontName',...
          'Times New Roman','Interpreter','latex');    
    view(2);
    grid off;
    axis equal;  
    dx = 0; dy = 0; dz = 0;
    xlim([min(nodes(:, 1)) - dx, max(nodes(:, 1)) + dx])
    ylim([min(nodes(:, 2)) - dy, max(nodes(:, 2)) + dy])  
    if min(grad.dy)~=max(grad.dy)
      zlim([min(grad.dy) - dz, max(grad.dy) + dz])
    end
    xlabel('$x$','FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    ylabel('$y$','FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    zlabel(mytitle,'FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    colorbar('FontName','Times New Roman','FontSize',14,'FontWeight','bold');
    colormap jet
    %set(gcf,'Renderer','painters')    
    set(gcf,'InvertHardcopy','off','Color',[1 1 1])
    set(gca,'FontName','Times New Roman','FontSize',14,'FontWeight','bold');
  end

  if strcmp(config.poisson2d_plot_grad.dnorm,'yes')    
    vq=griddata(gp_list.x(:),gp_list.y(:),gradNorm,xq,yq);
    figure
    mesh(xq,yq,vq,'FaceColor','interp','EdgeColor','interp');
    mytitle='$||\nabla u^h||$';
    title(mytitle,'FontWeight','bold','FontSize',20,'FontName',...
          'Times New Roman','Interpreter','latex');   
    view(2);
    grid off;
    axis equal;  
    dx = 0; dy = 0; dz = 0;
    xlim([min(nodes(:, 1)) - dx, max(nodes(:, 1)) + dx])
    ylim([min(nodes(:, 2)) - dy, max(nodes(:, 2)) + dy])  
    if min(gradNorm)~=max(gradNorm)
      zlim([min(gradNorm) - dz, max(gradNorm) + dz])
    end
    xlabel('$x$','FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    ylabel('$y$','FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    zlabel(mytitle,'FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    colorbar('FontName','Times New Roman','FontSize',14,'FontWeight','bold');
    colormap jet
    %set(gcf,'Renderer','painters')    
    set(gcf,'InvertHardcopy','off','Color',[1 1 1])
    set(gca,'FontName','Times New Roman','FontSize',14,'FontWeight','bold');
  end
end

