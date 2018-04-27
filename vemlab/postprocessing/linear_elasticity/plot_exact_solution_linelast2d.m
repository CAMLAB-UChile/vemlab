%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                     VEMLab
%           Source code  : http://camlab.cl/research/software/vemlab/
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:              plot_exact_solution_linelast2d 
%
% Created by : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
% Updated by :
% Updated by :
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% Plot displacement, stress and strain field exact solutions.
%
% Usage
% =====
% [flux,grad] = plot_exact_solution_linelast2d(domainMesh,exact_solution_handle,...
%                                              matProps,config)
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
% stress : struct. storing numerical stress tensor at Gauss points
% strain : struct. storing numerical strain tensor at Gauss points
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

function plot_exact_solution_linelast2d(domainMesh,exact_solution_handle,...
                                        matProps,config)
  % plot u field       
  [exact_nodal_solution,~,~]=exact_solutions_linelast2d(exact_solution_handle,...
                                                        domainMesh.coords);  
  plot_u_field(domainMesh,exact_nodal_solution,config);
  % plot flux and gradient
  plot_stress_and_strain(exact_solution_handle,domainMesh,matProps,config);
end

function plot_u_field(domainMesh,solution,config)
  solutionType=config.vemlab_method;
  plotMeshOverResults=config.plot_mesh_over_results; 
  
  fprintf('\n'); 
  fprintf('Plotting %s solution...\n',solutionType);                               
  nodes=domainMesh.coords;
  nodesNumber=size(nodes,1);  
  polygons=domainMesh.connect;
  
  range=1:nodesNumber;
  displacementsX=solution(2*range-1);
  displacementsY=solution(2*range);  
  displacementsNorm=sqrt(displacementsX.*displacementsX+displacementsY.*displacementsY);
  
  titleSolutionX='$u_x$'; 
  titleSolutionY='$u_y$'; 
  titleSolutionNorm='$||u||$';    

  if strcmp(config.linelast2d_plot_displacement.unorm,'yes')
    figure;
    title(titleSolutionNorm,'FontWeight','bold','FontSize',20,'FontName',...
          'Times New Roman','Interpreter','latex');
    maxNumVertices = max(cellfun(@numel,polygons));
    padFunc = @(vertList) [vertList' NaN(1,maxNumVertices-numel(vertList))];
    elements = cellfun(padFunc,polygons,'UniformOutput',false);
    elements = vertcat(elements{:});
    data = [nodes,displacementsNorm];
    if strcmp(plotMeshOverResults,'yes')
      patch('Faces',elements,'Vertices',data,...
            'FaceColor','interp','CData',displacementsNorm);  
    else
      patch('Faces',elements,'Vertices',data,...
            'FaceColor','interp','EdgeColor','interp','CData',displacementsNorm);
    end  
  %   axis('square') 
    axis equal;  
    dx = 0; dy = 0; dz = 0;
    xlim([min(nodes(:, 1)) - dx, max(nodes(:, 1)) + dx])
    ylim([min(nodes(:, 2)) - dy, max(nodes(:, 2)) + dy])  
    if min(displacementsNorm)~=max(displacementsNorm)
      zlim([min(displacementsNorm) - dz, max(displacementsNorm) + dz])
    end
    xlabel('$x$','FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    ylabel('$y$','FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    zlabel(titleSolutionNorm,'FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    colorbar('FontName','Times New Roman','FontSize',14,'FontWeight','bold');
    colormap jet
    %set(gcf,'Renderer','painters')    
    set(gcf,'InvertHardcopy','off','Color',[1 1 1])
    set(gca,'FontName','Times New Roman','FontSize',14,'FontWeight','bold');
  end
  
  if strcmp(config.linelast2d_plot_displacement.ux,'yes')  
    figure; 
    title(titleSolutionX,'FontWeight','bold','FontSize',20,'FontName',...
          'Times New Roman','Interpreter','latex');
    maxNumVertices = max(cellfun(@numel,polygons));
    padFunc = @(vertList) [vertList' NaN(1,maxNumVertices-numel(vertList))];
    elements = cellfun(padFunc,polygons,'UniformOutput',false);
    elements = vertcat(elements{:});
    data = [nodes,displacementsX];
    if strcmp(plotMeshOverResults,'yes')
      patch('Faces',elements,'Vertices',data,...
            'FaceColor','interp','CData',displacementsX);  
    else
      patch('Faces',elements,'Vertices',data,...
            'FaceColor','interp','EdgeColor','interp','CData',displacementsX);
    end      
  %   axis('square') 
    axis equal;  
    dx = 0; dy = 0; dz = 0;
    xlim([min(nodes(:,1)) - dx, max(nodes(:,1)) + dx])
    ylim([min(nodes(:,2)) - dy, max(nodes(:,2)) + dy])
    if min(displacementsX)~=max(displacementsX)
      zlim([min(displacementsX) - dz, max(displacementsX) + dz])
    end
    xlabel('$x$','FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    ylabel('$y$','FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    zlabel(titleSolutionX,'FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    colorbar('FontName','Times New Roman','FontSize',14,'FontWeight','bold');
    colormap jet
    %set(gcf,'Renderer','painters')    
    set(gcf,'InvertHardcopy','off','Color',[1 1 1])
    set(gca,'FontName','Times New Roman','FontSize',14,'FontWeight','bold');
  end
  
  if strcmp(config.linelast2d_plot_displacement.uy,'yes')  
    figure; 
    title(titleSolutionY,'FontWeight','bold','FontSize',20,'FontName',...
          'Times New Roman','Interpreter','latex');
    maxNumVertices = max(cellfun(@numel,polygons));
    padFunc = @(vertList) [vertList' NaN(1,maxNumVertices-numel(vertList))];
    elements = cellfun(padFunc,polygons,'UniformOutput',false);
    elements = vertcat(elements{:});
    data = [nodes,displacementsY];
    if strcmp(plotMeshOverResults,'yes')
      patch('Faces',elements,'Vertices',data,...
            'FaceColor','interp','CData',displacementsY);  
    else
      patch('Faces',elements,'Vertices',data,...
            'FaceColor','interp','EdgeColor','interp','CData',displacementsY);
    end   
  %   axis('square') 
    axis equal;  
    dx = 0; dy = 0; dz = 0;
    xlim([min(nodes(:,1)) - dx, max(nodes(:,1)) + dx])
    ylim([min(nodes(:,2)) - dy, max(nodes(:,2)) + dy])
    if min(displacementsY)~=max(displacementsY)
      zlim([min(displacementsY) - dz, max(displacementsY) + dz])
    end
    xlabel('$x$','FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    ylabel('$y$','FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    zlabel(titleSolutionY,'FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    colorbar('FontName','Times New Roman','FontSize',14,'FontWeight','bold');
    colormap jet
    %set(gcf,'Renderer','painters')    
    set(gcf,'InvertHardcopy','off','Color',[1 1 1])
    set(gca,'FontName','Times New Roman','FontSize',14,'FontWeight','bold'); 
  end
end

function plot_stress_and_strain(exact_solution_handle,domainMesh,matProps,config)
  if strcmp(config.linelast2d_plot_stress_and_strain,'yes')
    [stress,strain,gp_list,h_min,xmin,xmax,ymin,ymax]=...
                  exact_stress_and_strain_linelast2d(exact_solution_handle,...
                                                     domainMesh,matProps,config);
    plot_stress_and_strain_linelast2d(domainMesh,stress,strain,gp_list,h_min,...
                                      xmin,xmax,ymin,ymax,config);
  end
end

function [stress,strain,gp_list,h_min,xmin,xmax,ymin,ymax] =...
                    exact_stress_and_strain_linelast2d(exact_solution_handle,...
                                                       domainMesh,matProps,...
                                                       config) 
  fprintf('\n');
  fprintf('Postprocessing exact stress and strain at Gauss points...\n');
  
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
        throw_error('Error in plot_exact_solution_linelast2d.m --> plot_stress_and_strain --> exact_stress_and_strain_linelast2d: vemlab_method\n');
      end
      for tr_i=1:size(connect,1)
        gp=gp+1;   
        % assign constant stress and strain to the subtriangle Gauss point (1-pt rule)
        subtriangle_coords=domainMesh.coords(connect(tr_i,:),:);
        [xy,~]=gauss_points_T3(1,subtriangle_coords); % int. point in the form [x y]       
        gp_list.x(gp)=xy(1);
        gp_list.y(gp)=xy(2);        
        % exact solution
        [~,dudx_exact,dudy_exact]=exact_solutions_linelast2d(exact_solution_handle,xy);
        % exact strain
        du_exact_xgp=[dudx_exact;dudy_exact]; % [duxdx; duydx; duxdy; duydy]
        str_xgp=[du_exact_xgp(1);du_exact_xgp(4);...
                 du_exact_xgp(3)+du_exact_xgp(2)];  % [e11,e22,2*e12]   
        if strcmp(matProps.plane_state,'plane_stress')
          eten_xgp=[str_xgp(1),str_xgp(3)/2,0;... % for tensor representation divide 
                    str_xgp(3)/2,str_xgp(2),0;...   % strain by 2: (2*e12)/2
                    0,0,-matProps.nu*(str_xgp(1)+str_xgp(2))];  
        elseif strcmp(matProps.plane_state,'plane_strain')
          eten_xgp=[str_xgp(1),str_xgp(3)/2,0;...
                    str_xgp(3)/2,str_xgp(2),0;...
                    0,0,0];      
        else
          throw_error('Error in plot_exact_solution_linelast2d.m --> plot_stress_and_strain --> exact_stress_and_strain_linelast2d: plane_state\n');
        end
        pstrain_xgp=sort(eig(eten_xgp),'descend'); % principal strain: pstrain1 > pstrain2 > pstrain3      
        strain.e11(gp)=eten_xgp(1,1);
        strain.e12(gp)=eten_xgp(1,2);   
        strain.e22(gp)=eten_xgp(2,2); 
        strain.e33(gp)=eten_xgp(3,3);  
        strain.e1(gp)=pstrain_xgp(1);
        strain.e2(gp)=pstrain_xgp(2);
        strain.e3(gp)=pstrain_xgp(3);        
        % exact stress
        D=matProps.D;    % This D is defined as per VEM for strain = [e11 e22 e12].
        D(3,3)=D(3,3)/4; % To avoid troubles, come back to the standard FEM definition of D 
                         % to obtain the stress by multiplying with [e11 e22 2*e12]
        svec_xgp=D*[eten_xgp(1,1);eten_xgp(2,2);2*eten_xgp(1,2)];
        if strcmp(matProps.plane_state,'plane_stress')
          sten_xgp=[svec_xgp(1),svec_xgp(3),0;svec_xgp(3),svec_xgp(2),0;0,0,0];
        elseif strcmp(matProps.plane_state,'plane_strain')
          sten_xgp=[svec_xgp(1),svec_xgp(3),0;svec_xgp(3),svec_xgp(2),0;0,0,matProps.nu*(svec_xgp(1)+svec_xgp(2))];      
        else
          throw_error('Error in plot_numerical_solution_linelast2d.m --> plot_stress_and_strain --> fem_stress_and_strain_linelast2d: plane_state\n');      
        end
        pstress_h = sort(eig(sten_xgp),'descend');  % principal stresses: pstress1 > pstress2 > pstress3  
        sxx_xgp=sten_xgp(1,1); syy_xgp=sten_xgp(2,2); szz_xgp=sten_xgp(3,3);
        sxy_xgp=sten_xgp(1,2); sxz_xgp=0.0; syz_xgp=0.0;
        VM_xgp=(1/sqrt(2))*sqrt((sxx_xgp-syy_xgp)^2+(syy_xgp-szz_xgp)^2+(szz_xgp-sxx_xgp)^2+6*(sxy_xgp^2+syz_xgp^2+sxz_xgp^2));    
        stress.s11(gp)=sten_xgp(1,1);
        stress.s12(gp)=sten_xgp(1,2);   
        stress.s22(gp)=sten_xgp(2,2); 
        stress.s33(gp)=sten_xgp(3,3);   
        stress.s1(gp)=pstress_h(1);
        stress.s2(gp)=pstress_h(2);
        stress.s3(gp)=pstress_h(3);  
        stress.vm(gp)=VM_xgp;
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
          [~,dudx_exact,dudy_exact]=exact_solutions_linelast2d(exact_solution_handle,xy);           
          % exact strain
          du_exact_xgp=[dudx_exact;dudy_exact]; % [duxdx; duydx; duxdy; duydy]
          str_xgp=[du_exact_xgp(1);du_exact_xgp(4);...
                   du_exact_xgp(3)+du_exact_xgp(2)];  % [e11,e22,2*e12]   
          if strcmp(matProps.plane_state,'plane_stress')
            eten_xgp=[str_xgp(1),str_xgp(3)/2,0;... % for tensor representation divide 
                      str_xgp(3)/2,str_xgp(2),0;...   % strain by 2: (2*e12)/2
                      0,0,-matProps.nu*(str_xgp(1)+str_xgp(2))];  
          elseif strcmp(matProps.plane_state,'plane_strain')
            eten_xgp=[str_xgp(1),str_xgp(3)/2,0;...
                      str_xgp(3)/2,str_xgp(2),0;...
                      0,0,0];      
          else
            throw_error('Error in plot_exact_solution_linelast2d.m --> plot_stress_and_strain --> exact_stress_and_strain_linelast2d: plane_state\n');
          end
          pstrain_xgp=sort(eig(eten_xgp),'descend'); % principal strain: pstrain1 > pstrain2 > pstrain3      
          strain.e11(gp)=eten_xgp(1,1);
          strain.e12(gp)=eten_xgp(1,2);   
          strain.e22(gp)=eten_xgp(2,2); 
          strain.e33(gp)=eten_xgp(3,3);  
          strain.e1(gp)=pstrain_xgp(1);
          strain.e2(gp)=pstrain_xgp(2);
          strain.e3(gp)=pstrain_xgp(3);        
          % exact stress
          D=matProps.D;    % This D is defined as per VEM for strain = [e11 e22 e12].
          D(3,3)=D(3,3)/4; % To avoid troubles, come back to the standard FEM definition of D 
                           % to obtain the stress by multiplying with [e11 e22 2*e12]
          svec_xgp=D*[eten_xgp(1,1);eten_xgp(2,2);2*eten_xgp(1,2)];
          if strcmp(matProps.plane_state,'plane_stress')
            sten_xgp=[svec_xgp(1),svec_xgp(3),0;svec_xgp(3),svec_xgp(2),0;0,0,0];
          elseif strcmp(matProps.plane_state,'plane_strain')
            sten_xgp=[svec_xgp(1),svec_xgp(3),0;svec_xgp(3),svec_xgp(2),0;0,0,matProps.nu*(svec_xgp(1)+svec_xgp(2))];      
          else
            throw_error('Error in plot_numerical_solution_linelast2d.m --> plot_stress_and_strain --> fem_stress_and_strain_linelast2d: plane_state\n');      
          end
          pstress_h = sort(eig(sten_xgp),'descend');  % principal stresses: pstress1 > pstress2 > pstress3  
          sxx_xgp=sten_xgp(1,1); syy_xgp=sten_xgp(2,2); szz_xgp=sten_xgp(3,3);
          sxy_xgp=sten_xgp(1,2); sxz_xgp=0.0; syz_xgp=0.0;
          VM_xgp=(1/sqrt(2))*sqrt((sxx_xgp-syy_xgp)^2+(syy_xgp-szz_xgp)^2+(szz_xgp-sxx_xgp)^2+6*(sxy_xgp^2+syz_xgp^2+sxz_xgp^2));    
          stress.s11(gp)=sten_xgp(1,1);
          stress.s12(gp)=sten_xgp(1,2);   
          stress.s22(gp)=sten_xgp(2,2); 
          stress.s33(gp)=sten_xgp(3,3);   
          stress.s1(gp)=pstress_h(1);
          stress.s2(gp)=pstress_h(2);
          stress.s3(gp)=pstress_h(3);  
          stress.vm(gp)=VM_xgp;
        end
      end
      % h_min
      h_size=max_edge_size(elem_coord);
      if h_size<h_min
        h_min=h_size;
      end        
    end
  else
    throw_error('Error in plot_exact_solution_linelast2d.m --> plot_stress_and_strain --> exact_stress_and_strain_linelast2d: vemlab_method\n');
  end
end

function plot_stress_and_strain_linelast2d(domainMesh,stress,strain,gp_list,...
                                           h_min,xmin,xmax,ymin,ymax,config)  
  nodes=domainMesh.coords;  
  [xq,yq]=meshgrid(xmin:0.05*h_min:xmax,ymin:0.05*h_min:ymax);  
  % plot numerical stress and strain
  if strcmp(config.linelast2d_plot_stress.s11,'yes')
    vq=griddata(gp_list.x(:),gp_list.y(:),stress.s11(:),xq,yq);
    figure
    mesh(xq,yq,vq,'FaceColor','interp','EdgeColor','interp');
    mytitle='$\sigma_x$';
    title(mytitle,'FontWeight','bold','FontSize',20,'FontName',...
          'Times New Roman','Interpreter','latex');  
    view(2);
    grid off;
    axis equal;  
    dx = 0; dy = 0; dz = 0;
    xlim([min(nodes(:, 1)) - dx, max(nodes(:, 1)) + dx])
    ylim([min(nodes(:, 2)) - dy, max(nodes(:, 2)) + dy])  
    if min(stress.s11)~=max(stress.s11)
      zlim([min(stress.s11) - dz, max(stress.s11) + dz])
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
  if strcmp(config.linelast2d_plot_stress.s12,'yes')  
    vq=griddata(gp_list.x(:),gp_list.y(:),stress.s12(:),xq,yq);
    figure
    mesh(xq,yq,vq,'FaceColor','interp','EdgeColor','interp');
    mytitle='$\sigma_{xy}$';
    title(mytitle,'FontWeight','bold','FontSize',20,'FontName',...
          'Times New Roman','Interpreter','latex');  
    view(2);
    grid off;
    axis equal;  
    dx = 0; dy = 0; dz = 0;
    xlim([min(nodes(:, 1)) - dx, max(nodes(:, 1)) + dx])
    ylim([min(nodes(:, 2)) - dy, max(nodes(:, 2)) + dy])  
    if min(stress.s12)~=max(stress.s12)
      zlim([min(stress.s12) - dz, max(stress.s12) + dz])
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
  if strcmp(config.linelast2d_plot_stress.s22,'yes')   
    vq=griddata(gp_list.x(:),gp_list.y(:),stress.s22(:),xq,yq);
    figure
    mesh(xq,yq,vq,'FaceColor','interp','EdgeColor','interp');
    mytitle='$\sigma_y$';
    title(mytitle,'FontWeight','bold','FontSize',20,'FontName',...
          'Times New Roman','Interpreter','latex');  
    view(2);
    grid off;
    axis equal;  
    dx = 0; dy = 0; dz = 0;
    xlim([min(nodes(:, 1)) - dx, max(nodes(:, 1)) + dx])
    ylim([min(nodes(:, 2)) - dy, max(nodes(:, 2)) + dy])  
    if min(stress.s22)~=max(stress.s22)
      zlim([min(stress.s22) - dz, max(stress.s22) + dz])
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
  if strcmp(config.linelast2d_plot_stress.s33,'yes')   
    vq=griddata(gp_list.x(:),gp_list.y(:),stress.s33(:),xq,yq);
    figure
    mesh(xq,yq,vq,'FaceColor','interp','EdgeColor','interp');
    mytitle='$\sigma_z$';
    title(mytitle,'FontWeight','bold','FontSize',20,'FontName',...
          'Times New Roman','Interpreter','latex');  
    view(2);
    grid off;
    axis equal;  
    dx = 0; dy = 0; dz = 0;
    xlim([min(nodes(:, 1)) - dx, max(nodes(:, 1)) + dx])
    ylim([min(nodes(:, 2)) - dy, max(nodes(:, 2)) + dy]) 
    if min(stress.s33)~=max(stress.s33)
      zlim([min(stress.s33) - dz, max(stress.s33) + dz])
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
  if strcmp(config.linelast2d_plot_stress.s1,'yes')    
    vq=griddata(gp_list.x(:),gp_list.y(:),stress.s1(:),xq,yq);
    figure
    mesh(xq,yq,vq,'FaceColor','interp','EdgeColor','interp');
    mytitle='$\sigma_1$';
    title(mytitle,'FontWeight','bold','FontSize',20,'FontName',...
          'Times New Roman','Interpreter','latex');  
    view(2);
    grid off;
    axis equal;  
    dx = 0; dy = 0; dz = 0;
    xlim([min(nodes(:, 1)) - dx, max(nodes(:, 1)) + dx])
    ylim([min(nodes(:, 2)) - dy, max(nodes(:, 2)) + dy])  
    if min(stress.s1)~=max(stress.s1)
      zlim([min(stress.s1) - dz, max(stress.s1) + dz])
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
  if strcmp(config.linelast2d_plot_stress.s2,'yes')     
    vq=griddata(gp_list.x(:),gp_list.y(:),stress.s2(:),xq,yq);
    figure
    mesh(xq,yq,vq,'FaceColor','interp','EdgeColor','interp');
    mytitle='$\sigma_2$';
    title(mytitle,'FontWeight','bold','FontSize',20,'FontName',...
          'Times New Roman','Interpreter','latex');  
    view(2);
    grid off;
    axis equal;  
    dx = 0; dy = 0; dz = 0;
    xlim([min(nodes(:, 1)) - dx, max(nodes(:, 1)) + dx])
    ylim([min(nodes(:, 2)) - dy, max(nodes(:, 2)) + dy])  
    if min(stress.s2)~=max(stress.s2)
      zlim([min(stress.s2) - dz, max(stress.s2) + dz])
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
  if strcmp(config.linelast2d_plot_stress.s3,'yes')     
    vq=griddata(gp_list.x(:),gp_list.y(:),stress.s3(:),xq,yq);
    figure
    mesh(xq,yq,vq,'FaceColor','interp','EdgeColor','interp');
    mytitle='$\sigma_3$';
    title(mytitle,'FontWeight','bold','FontSize',20,'FontName',...
          'Times New Roman','Interpreter','latex');  
    view(2);
    grid off;
    axis equal;  
    dx = 0; dy = 0; dz = 0;
    xlim([min(nodes(:, 1)) - dx, max(nodes(:, 1)) + dx])
    ylim([min(nodes(:, 2)) - dy, max(nodes(:, 2)) + dy])  
    if min(stress.s3)~=max(stress.s3)
      zlim([min(stress.s3) - dz, max(stress.s3) + dz])
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
  if strcmp(config.linelast2d_plot_stress.vm,'yes')     
    vq=griddata(gp_list.x(:),gp_list.y(:),stress.vm(:),xq,yq);
    figure
    mesh(xq,yq,vq,'FaceColor','interp','EdgeColor','interp');
    mytitle='$\sigma_{\mathrm{vm}}$';
    title(mytitle,'FontWeight','bold','FontSize',20,'FontName',...
          'Times New Roman','Interpreter','latex');  
    view(2);
    grid off;
    axis equal;  
    dx = 0; dy = 0; dz = 0;
    xlim([min(nodes(:, 1)) - dx, max(nodes(:, 1)) + dx])
    ylim([min(nodes(:, 2)) - dy, max(nodes(:, 2)) + dy])  
    if min(stress.vm)~=max(stress.vm)
      zlim([min(stress.vm) - dz, max(stress.vm) + dz])
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
  if strcmp(config.linelast2d_plot_strain.e11,'yes')     
    vq=griddata(gp_list.x(:),gp_list.y(:),strain.e11(:),xq,yq);
    figure
    mesh(xq,yq,vq,'FaceColor','interp','EdgeColor','interp');
    mytitle='$\varepsilon_x$';
    title(mytitle,'FontWeight','bold','FontSize',20,'FontName',...
          'Times New Roman','Interpreter','latex');  
    view(2);
    grid off;
    axis equal;  
    dx = 0; dy = 0; dz = 0;
    xlim([min(nodes(:, 1)) - dx, max(nodes(:, 1)) + dx])
    ylim([min(nodes(:, 2)) - dy, max(nodes(:, 2)) + dy])  
    if min(strain.e11)~=max(strain.e11)
      zlim([min(strain.e11) - dz, max(strain.e11) + dz])
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
  if strcmp(config.linelast2d_plot_strain.e12,'yes')    
    vq=griddata(gp_list.x(:),gp_list.y(:),strain.e12(:),xq,yq);
    figure
    mesh(xq,yq,vq,'FaceColor','interp','EdgeColor','interp');
    mytitle='$\varepsilon_{xy}$';
    title(mytitle,'FontWeight','bold','FontSize',20,'FontName',...
          'Times New Roman','Interpreter','latex');  
    view(2);
    grid off;
    axis equal;  
    dx = 0; dy = 0; dz = 0;
    xlim([min(nodes(:, 1)) - dx, max(nodes(:, 1)) + dx])
    ylim([min(nodes(:, 2)) - dy, max(nodes(:, 2)) + dy])  
    if min(strain.e12)~=max(strain.e12)
      zlim([min(strain.e12) - dz, max(strain.e12) + dz])
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
  if strcmp(config.linelast2d_plot_strain.e22,'yes')    
    vq=griddata(gp_list.x(:),gp_list.y(:),strain.e22(:),xq,yq);
    figure
    mesh(xq,yq,vq,'FaceColor','interp','EdgeColor','interp');
    mytitle='$\varepsilon_y$';
    title(mytitle,'FontWeight','bold','FontSize',20,'FontName',...
          'Times New Roman','Interpreter','latex');  
    view(2);
    grid off;
    axis equal;  
    dx = 0; dy = 0; dz = 0;
    xlim([min(nodes(:, 1)) - dx, max(nodes(:, 1)) + dx])
    ylim([min(nodes(:, 2)) - dy, max(nodes(:, 2)) + dy])  
    if min(strain.e22)~=max(strain.e22)
      zlim([min(strain.e22) - dz, max(strain.e22) + dz])
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
  if strcmp(config.linelast2d_plot_strain.e33,'yes')    
    vq=griddata(gp_list.x(:),gp_list.y(:),strain.e33(:),xq,yq);
    figure
    mesh(xq,yq,vq,'FaceColor','interp','EdgeColor','interp');
    mytitle='$\varepsilon_z$';
    title(mytitle,'FontWeight','bold','FontSize',20,'FontName',...
          'Times New Roman','Interpreter','latex');  
    view(2);
    grid off;
    axis equal;  
    dx = 0; dy = 0; dz = 0;
    xlim([min(nodes(:, 1)) - dx, max(nodes(:, 1)) + dx])
    ylim([min(nodes(:, 2)) - dy, max(nodes(:, 2)) + dy])  
    if min(strain.e33)~=max(strain.e33)
      zlim([min(strain.e33) - dz, max(strain.e33) + dz])
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
  if strcmp(config.linelast2d_plot_strain.e1,'yes')    
    vq=griddata(gp_list.x(:),gp_list.y(:),strain.e1(:),xq,yq);
    figure
    mesh(xq,yq,vq,'FaceColor','interp','EdgeColor','interp');
    mytitle='$\varepsilon_1$';
    title(mytitle,'FontWeight','bold','FontSize',20,'FontName',...
          'Times New Roman','Interpreter','latex');  
    view(2);
    grid off;
    axis equal;  
    dx = 0; dy = 0; dz = 0;
    xlim([min(nodes(:, 1)) - dx, max(nodes(:, 1)) + dx])
    ylim([min(nodes(:, 2)) - dy, max(nodes(:, 2)) + dy])  
    if min(strain.e1)~=max(strain.e1)
      zlim([min(strain.e1) - dz, max(strain.e1) + dz])
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
  if strcmp(config.linelast2d_plot_strain.e2,'yes')    
    vq=griddata(gp_list.x(:),gp_list.y(:),strain.e2(:),xq,yq);
    figure
    mesh(xq,yq,vq,'FaceColor','interp','EdgeColor','interp');
    mytitle='$\varepsilon_2$';
    title(mytitle,'FontWeight','bold','FontSize',20,'FontName',...
          'Times New Roman','Interpreter','latex');  
    view(2);
    grid off;
    axis equal;  
    dx = 0; dy = 0; dz = 0;
    xlim([min(nodes(:, 1)) - dx, max(nodes(:, 1)) + dx])
    ylim([min(nodes(:, 2)) - dy, max(nodes(:, 2)) + dy])  
    if min(strain.e2)~=max(strain.e2)
      zlim([min(strain.e2) - dz, max(strain.e2) + dz])
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
  if strcmp(config.linelast2d_plot_strain.e3,'yes')    
    vq=griddata(gp_list.x(:),gp_list.y(:),strain.e3(:),xq,yq);
    figure
    mesh(xq,yq,vq,'FaceColor','interp','EdgeColor','interp');
    mytitle='$\varepsilon_3$';
    title(mytitle,'FontWeight','bold','FontSize',20,'FontName',...
          'Times New Roman','Interpreter','latex');  
    view(2);
    grid off;
    axis equal;  
    dx = 0; dy = 0; dz = 0;
    xlim([min(nodes(:, 1)) - dx, max(nodes(:, 1)) + dx])
    ylim([min(nodes(:, 2)) - dy, max(nodes(:, 2)) + dy])  
    if min(strain.e3)~=max(strain.e3)
      zlim([min(strain.e3) - dz, max(strain.e3) + dz])
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

