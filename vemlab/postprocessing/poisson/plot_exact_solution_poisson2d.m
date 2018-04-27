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
% Mar. 17, 2018: first realease (by A. Ortiz-Bernardin)
% Apr. 19, 2018: improve the plotting of axis and fonts
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function plot_exact_solution_poisson2d(domainMesh,exact_solution_handle,...
                                       matProps,config)
  % plot u field 
  [exact_nodal_solution,~,~]=exact_solutions_poisson2d(exact_solution_handle,...
                                                       domainMesh.coords,...
                                                       matProps);
  plot_u_field(domainMesh,exact_nodal_solution,config);
  % plot flux and gradient
  plot_flux_and_gradient(exact_solution_handle,domainMesh,matProps,config);
end

function plot_u_field(domainMesh,solution,config)
  solutionType='exact';
  plotMeshOverResults=config.plot_mesh_over_results;
  
  fprintf('\n');
  fprintf('Plotting %s solution...\n',solutionType); 
  
  nodes=domainMesh.coords; 
  polygons=domainMesh.connect;
  titleSolution='$u$';  
  
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

function plot_flux_and_gradient(exact_solution_handle,domainMesh,matProps,config)
  [flux,grad,gp_list,h_min,xmin,xmax,ymin,ymax]=...
              exact_flux_and_gradient_poisson2d(exact_solution_handle,...
                                                domainMesh,matProps,config);
  plot_flux_and_gradient_poisson2d(domainMesh,flux,grad,gp_list,h_min,...
                                   xmin,xmax,ymin,ymax,config);
end

function [flux,grad,gp_list,h_min,xmin,xmax,ymin,ymax] = ...
                   exact_flux_and_gradient_poisson2d(exact_solution_handle,...
                                                     domainMesh,matProps,...
                                                     config) 
  fprintf('\n');
  fprintf('Postprocessing exact flux and gradient at Gauss points...\n');
  
  % The bounding box must contain the domain. If the domain is a rectangle, the
  % bounding box coincides with the domain
  xmin=domainMesh.BdBox(1); ymin=domainMesh.BdBox(3);
  xmax=domainMesh.BdBox(2); ymax=domainMesh.BdBox(4);
  dx=abs(xmax-xmin); dy=abs(ymax-ymin);
  h_min=max(dx,dy);
  
  gp=0;   
  num_elem=length(domainMesh.connect);  
  if strcmp(config.vemlab_method,'VEM2D')||strcmp(config.vemlab_method,'FEM2DT3')
    % loop over elements
    for i=1:num_elem
      if strcmp(config.vemlab_method,'VEM2D')
        connect=triangulate_polygon(domainMesh,i);
      elseif strcmp(config.vemlab_method,'FEM2DT3')
        node_indices=domainMesh.connect(i,:);
        connect=(node_indices{1}(:))';
      else
        throw_error('Error in plot_exact_solution_poisson2d.m --> plot_flux_and_gradient --> exact_flux_and_gradient_poisson2d: vemlab_method\n');
      end
      for tr_i=1:size(connect,1)
        gp=gp+1;   
        % assign constant flux and gradient to the subtriangle Gauss point (1-pt rule)
        subtriangle_coords=domainMesh.coords(connect(tr_i,:),:);
        [xy,~]=gauss_points_T3(1,subtriangle_coords); % int. point in the form [x y]       
        gp_list.x(gp)=xy(1);
        gp_list.y(gp)=xy(2);        
        % exact solution
        [~,flux_exact,grad_exact]=...
                  exact_solutions_poisson2d(exact_solution_handle,xy,matProps);
        grad.dx(gp)=grad_exact.x;
        grad.dy(gp)=grad_exact.y;  
        flux.qx(gp)=flux_exact.x;
        flux.qy(gp)=flux_exact.y;          
        % minimum size of the subtriangular mesh
        nodes_indices_T3=connect(tr_i,:);
        verts=domainMesh.coords(nodes_indices_T3,:);
        h_size=max_edge_size(verts);
        if h_size<h_min
          h_min=h_size;
        end 
      end
    end
  elseif strcmp(config.vemlab_method,'FEM2DQ4')
    % loop over elements
    gp=0;    
    coords=domainMesh.coords; 
    connect=domainMesh.connect;    
    for i=1:num_elem
      nodes=connect{i};
      elem_coord=coords(nodes,:);
      % compute flux and gradient at Gauss points
      int_order=2;
      xi=gauss_points_1d(int_order);
      eta=gauss_points_1d(int_order);
      for gpxi=1:length(xi)     
        for gpeta=1:length(eta)
          gp=gp+1;
          % Gauss point coordinates
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
          % exact flux and gradient
          [~,flux_exact,grad_exact]=...
                  exact_solutions_poisson2d(exact_solution_handle,xy,matProps);           
          grad.dx(gp)=grad_exact.x;
          grad.dy(gp)=grad_exact.y;  
          flux.qx(gp)=flux_exact.x;
          flux.qy(gp)=flux_exact.y;             
        end
      end
    end
  else
    throw_error('Error in plot_exact_solution_poisson2d.m --> plot_flux_and_gradient --> exact_flux_and_gradient_poisson2d: vemlab_method\n');
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
    mytitle='$q_x$';
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
    mytitle='$q_y$';
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
    mytitle='$||q||$';
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
    mytitle='$\frac{\partial u}{\partial x}$';
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
    mytitle='$\frac{\partial u}{\partial y}$';
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
    mytitle='$||\nabla u||$';
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

