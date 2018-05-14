%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                     VEMLab
%           Source code  : http://camlab.cl/research/software/vemlab/
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:              plot_numerical_solution_linelast2d 
%
% Created by : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
% Updated by :
% Updated by :
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% Plot displacement, stress and strain field numerical solutions.
%
% Usage
% =====
% [flux,grad] = plot_numerical_solution_linelast2d(domainMesh,solution,...
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


function [stress,strain] = plot_numerical_solution_linelast2d(domainMesh,...
                                                              solution,...
                                                              matProps,config)
  % plot u field                                       
  plot_u_field(domainMesh,solution,config);
  % plot flux and gradient
  [stress,strain]=plot_stress_and_strain(solution,domainMesh,matProps,config);
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
  
  titleSolutionX='$u_x^h$'; 
  titleSolutionY='$u_y^h$'; 
  titleSolutionNorm='$||u^h||$';   

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

function [stress,strain]=plot_stress_and_strain(solution,domainMesh,matProps,config)
  if strcmp(config.linelast2d_plot_stress_and_strain,'yes')
    if strcmp(config.vemlab_method,'VEM2D')
      [stress,strain,gp_list,h_min,xmin,xmax,ymin,ymax]=...
        vem_stress_and_strain_linelast2d(solution,domainMesh,matProps,config);
    elseif strcmp(config.vemlab_method,'FEM2DT3')||strcmp(config.vemlab_method,'FEM2DQ4')
      [stress,strain,gp_list,h_min,xmin,xmax,ymin,ymax]=...
        fem_stress_and_strain_linelast2d(solution,domainMesh,matProps,config);     
    else
      throw_error('Error in plot_numerical_solution_linelast2d.m --> plot_stress_and_strain: vemlab_method\n');
    end
    plot_stress_and_strain_linelast2d(domainMesh,stress,strain,gp_list,h_min,...
                                      xmin,xmax,ymin,ymax,config);
  else % return empty stress and strain
    stress=[]; strain=[];
  end
end

function [stress,strain,gp_list,h_min,xmin,xmax,ymin,ymax] = ...
       vem_stress_and_strain_linelast2d(solution,domainMesh,matProps,config) 
     
  fprintf('Postprocessing %s stress and strain at Gauss points...\n',...
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
    dofs=zeros(2*length(nodes),1);
    for node_i=1:length(nodes)
      dofs(2*node_i-1)=2*nodes(node_i)-1;
      dofs(2*node_i)=2*nodes(node_i);
    end
    % element VEM nodal solution
    uh_elem_column=solution(dofs);
    % elemen nodal coordinates
    for h=1:length(nodes)
      elem_coord(h,1)=coords(nodes(h),1);
      elem_coord(h,2)=coords(nodes(h),2);
    end     
    % VEM matrix
    verts=elem_coord;    
    Wc=vem_Wc_linelast2d(verts);
    % VEM strains
    pic_grad_uh=Wc'*uh_elem_column;
    if strcmp(matProps.plane_state,'plane_stress')
      eten_h=[pic_grad_uh(1),pic_grad_uh(3),0;...
            pic_grad_uh(3),pic_grad_uh(2),0;...
            0,0,-matProps.nu*(pic_grad_uh(1)+pic_grad_uh(2))];  
    elseif strcmp(matProps.plane_state,'plane_strain')
      eten_h=[pic_grad_uh(1),pic_grad_uh(3),0;...
            pic_grad_uh(3),pic_grad_uh(2),0;...
            0,0,0];      
    else
      throw_error('Error in plot_numerical_solution_linelast2d.m --> plot_stress_and_strain --> vem_stress_and_strain_linelast2d: plane_state\n');
    end
    pstrain_h=sort(eig(eten_h),'descend'); % principal strain: pstrain1 > pstrain2 > pstrain3      
    strain_h.e11(i)=eten_h(1,1);
    strain_h.e12(i)=eten_h(1,2);   
    strain_h.e22(i)=eten_h(2,2); 
    strain_h.e33(i)=eten_h(3,3);  
    strain_h.e1(i)=pstrain_h(1);
    strain_h.e2(i)=pstrain_h(2);
    strain_h.e3(i)=pstrain_h(3);    
    % VEM stresses
    D=matProps.D;    % This D is defined as per VEM for strain = [e11 e22 e12].
    D(3,3)=D(3,3)/4; % To avoid troubles, come back to the standard FEM definition of D 
                     % to obtain the stress by multiplying with [e11 e22 2*e12]    
    svec_h=D*[eten_h(1,1);eten_h(2,2);2*eten_h(1,2)];
    if strcmp(matProps.plane_state,'plane_stress')
      sten_h=[svec_h(1),svec_h(3),0;svec_h(3),svec_h(2),0;0,0,0];
    elseif strcmp(matProps.plane_state,'plane_strain')
      sten_h=[svec_h(1),svec_h(3),0;svec_h(3),svec_h(2),0;0,0,matProps.nu*(svec_h(1)+svec_h(2))];      
    else
      throw_error('Error in plot_numerical_solution_linelast2d.m --> plot_stress_and_strain --> vem_stress_and_strain_linelast2d: plane_state\n');      
    end
    pstress_h = sort(eig(sten_h),'descend');  % principal stresses: pstress1 > pstress2 > pstress3  
    sxx_h=sten_h(1,1); syy_h=sten_h(2,2); szz_h=sten_h(3,3);
    sxy_h=sten_h(1,2); sxz_h=0.0; syz_h=0.0;
    VM_h=(1/sqrt(2))*sqrt((sxx_h-syy_h)^2+(syy_h-szz_h)^2+(szz_h-sxx_h)^2+6*(sxy_h^2+syz_h^2+sxz_h^2));    
    stress_h.s11(i)=sten_h(1,1);
    stress_h.s12(i)=sten_h(1,2);   
    stress_h.s22(i)=sten_h(2,2); 
    stress_h.s33(i)=sten_h(3,3);   
    stress_h.s1(i)=pstress_h(1);
    stress_h.s2(i)=pstress_h(2);
    stress_h.s3(i)=pstress_h(3);  
    stress_h.vm(i)=VM_h; 
    % triangulate and assign the constant stress and strain to the subtriangles
    newconnect=triangulate_polygon(domainMesh,i);
    for tr_i=1:size(newconnect,1)
      gp=gp+1;   
      % assign constant stress and strain to the subtriangle Gauss point (1-pt rule)
      subtriangle_coords=domainMesh.coords(newconnect(tr_i,:),:);
      [xy,~]=gauss_points_T3(1,subtriangle_coords); % int. point in the form [x y]       
      gp_list.x(gp)=xy(1);
      gp_list.y(gp)=xy(2);        
      % exact solution
      strain.e11(gp)=strain_h.e11(i);
      strain.e12(gp)=strain_h.e12(i);
      strain.e22(gp)=strain_h.e22(i); 
      strain.e33(gp)=strain_h.e33(i);   
      strain.e1(gp)=strain_h.e1(i);   
      strain.e2(gp)=strain_h.e2(i); 
      strain.e3(gp)=strain_h.e3(i);      
      stress.s11(gp)=stress_h.s11(i);
      stress.s12(gp)=stress_h.s12(i);     
      stress.s22(gp)=stress_h.s22(i); 
      stress.s33(gp)=stress_h.s33(i);  
      stress.s1(gp)=stress_h.s1(i);  
      stress.s2(gp)=stress_h.s2(i);    
      stress.s3(gp)=stress_h.s3(i);
      stress.vm(gp)=stress_h.vm(i);             
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

function [stress,strain,gp_list,h_min,xmin,xmax,ymin,ymax] = ...
           fem_stress_and_strain_linelast2d(solution,domainMesh,matProps,config)
  
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
      dofs=zeros(2*length(nodes),1);
      for node_i=1:length(nodes)
        dofs(2*node_i-1)=2*nodes(node_i)-1;
        dofs(2*node_i)=2*nodes(node_i);
      end
      % element FEM nodal solution
      uh_elem_column=solution(dofs); 
      % area of the element
      area=triangle_area(elem_coord);
      % deformation matrix
      x=elem_coord(:,1); y=elem_coord(:,2);  
      % B is defined for [e11,e22,2*e12]
      B=(1/(2*area))*[y(2)-y(3),0,y(3)-y(1),0,y(1)-y(2),0;...
                      0,x(3)-x(2),0,x(1)-x(3),0,x(2)-x(1);...
                      x(3)-x(2),y(2)-y(3),x(1)-x(3),y(3)-y(1),x(2)-x(1),y(1)-y(2)];      
      % FEM strains
      str_h=B*uh_elem_column; % In FEM, B is defined for [e11,e22,2*e12]
      if strcmp(matProps.plane_state,'plane_stress')
        eten_h=[str_h(1),str_h(3)/2,0;... % for tensor representation divide 
              str_h(3)/2,str_h(2),0;... % strain coming from B by 2: (2*e12)/2
              0,0,-matProps.nu*(str_h(1)+str_h(2))];  
      elseif strcmp(matProps.plane_state,'plane_strain')
        eten_h=[str_h(1),str_h(3)/2,0;...
              str_h(3)/2,str_h(2),0;...
              0,0,0];      
      else
        throw_error('Error in plot_numerical_solution_linelast2d.m --> plot_stress_and_strain --> fem_stress_and_strain_linelast2d: plane_state\n');
      end
      pstrain_h=sort(eig(eten_h),'descend'); % principal strain: pstrain1 > pstrain2 > pstrain3      
      strain_h.e11(i)=eten_h(1,1);
      strain_h.e12(i)=eten_h(1,2);   
      strain_h.e22(i)=eten_h(2,2); 
      strain_h.e33(i)=eten_h(3,3);  
      strain_h.e1(i)=pstrain_h(1);
      strain_h.e2(i)=pstrain_h(2);
      strain_h.e3(i)=pstrain_h(3);    
      % FEM stresses
      D=matProps.D;    % This D is defined as per VEM for strain = [e11 e22 e12].
      D(3,3)=D(3,3)/4; % To avoid troubles, come back to the standard FEM definition of D 
                       % to obtain the stress by multiplying with [e11 e22 2*e12]
      svec_h=D*[eten_h(1,1);eten_h(2,2);2*eten_h(1,2)];
      if strcmp(matProps.plane_state,'plane_stress')
        sten_h=[svec_h(1),svec_h(3),0;svec_h(3),svec_h(2),0;0,0,0];
      elseif strcmp(matProps.plane_state,'plane_strain')
        sten_h=[svec_h(1),svec_h(3),0;svec_h(3),svec_h(2),0;0,0,matProps.nu*(svec_h(1)+svec_h(2))];      
      else
        throw_error('Error in plot_numerical_solution_linelast2d.m --> plot_stress_and_strain --> fem_stress_and_strain_linelast2d: plane_state\n');      
      end
      pstress_h = sort(eig(sten_h),'descend');  % principal stresses: pstress1 > pstress2 > pstress3  
      sxx_h=sten_h(1,1); syy_h=sten_h(2,2); szz_h=sten_h(3,3);
      sxy_h=sten_h(1,2); sxz_h=0.0; syz_h=0.0;
      VM_h=(1/sqrt(2))*sqrt((sxx_h-syy_h)^2+(syy_h-szz_h)^2+(szz_h-sxx_h)^2+6*(sxy_h^2+syz_h^2+sxz_h^2));    
      stress_h.s11(i)=sten_h(1,1);
      stress_h.s12(i)=sten_h(1,2);   
      stress_h.s22(i)=sten_h(2,2); 
      stress_h.s33(i)=sten_h(3,3);   
      stress_h.s1(i)=pstress_h(1);
      stress_h.s2(i)=pstress_h(2);
      stress_h.s3(i)=pstress_h(3);  
      stress_h.vm(i)=VM_h;  
      % triangulate and assign the constant stress and strain to the subtriangles
      newconnect=triangulate_polygon(domainMesh,i);
      for tr_i=1:size(newconnect,1)
        gp=gp+1;   
        % assign constant stress and strain to the subtriangle Gauss point (1-pt rule)
        subtriangle_coords=domainMesh.coords(newconnect(tr_i,:),:);
        [xy,~]=gauss_points_T3(1,subtriangle_coords); % int. point in the form [x y]       
        gp_list.x(gp)=xy(1);
        gp_list.y(gp)=xy(2);        
        % exact solution
        strain.e11(gp)=strain_h.e11(i);
        strain.e12(gp)=strain_h.e12(i);
        strain.e22(gp)=strain_h.e22(i); 
        strain.e33(gp)=strain_h.e33(i);   
        strain.e1(gp)=strain_h.e1(i);   
        strain.e2(gp)=strain_h.e2(i); 
        strain.e3(gp)=strain_h.e3(i);      
        stress.s11(gp)=stress_h.s11(i);
        stress.s12(gp)=stress_h.s12(i);     
        stress.s22(gp)=stress_h.s22(i); 
        stress.s33(gp)=stress_h.s33(i);  
        stress.s1(gp)=stress_h.s1(i);  
        stress.s2(gp)=stress_h.s2(i);    
        stress.s3(gp)=stress_h.s3(i);
        stress.vm(gp)=stress_h.vm(i);             
        % minimum size of the subtriangular mesh
        nodes_indices_T3=newconnect(tr_i,1:3);
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
    for i=1:num_elem
      nodes=connect{i};
      elem_coord=coords(nodes,:);
      % element dofs global indices
      dofs=zeros(2*length(nodes),1);
      for node_i=1:length(nodes)
        dofs(2*node_i-1)=2*nodes(node_i)-1;
        dofs(2*node_i)=2*nodes(node_i);
      end
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
          B=[dN1dx,0,dN2dx,0,dN3dx,0,dN4dx,0;...
             0,dN1dy,0,dN2dy,0,dN3dy,0,dN4dy;...
             dN1dy,dN1dx,dN2dy,dN2dx,dN3dy,dN3dx,dN4dy,dN4dx];
          % FEM strains
          str_h=B*uh_elem_column; % In FEM, B is defined for [e11,e22,2*e12]
          if strcmp(matProps.plane_state,'plane_stress')
            eten_h=[str_h(1),str_h(3)/2,0;... % for tensor representation divide 
                  str_h(3)/2,str_h(2),0;... % strain coming from B by 2: (2*e12)/2
                  0,0,-matProps.nu*(str_h(1)+str_h(2))];  
          elseif strcmp(matProps.plane_state,'plane_strain')
            eten_h=[str_h(1),str_h(3)/2,0;...
                  str_h(3)/2,str_h(2),0;...
                  0,0,0];      
          else
            throw_error('Error in plot_numerical_solution_linelast2d.m --> plot_stress_and_strain --> fem_stress_and_strain_linelast2d: plane_state');
          end    
          pstrain_h=sort(eig(eten_h),'descend'); % principal strain: pstrain1 > pstrain2 > pstrain3      
          strain.e11(gp)=eten_h(1,1);
          strain.e12(gp)=eten_h(1,2);   
          strain.e22(gp)=eten_h(2,2); 
          strain.e33(gp)=eten_h(3,3);  
          strain.e1(gp)=pstrain_h(1);
          strain.e2(gp)=pstrain_h(2);
          strain.e3(gp)=pstrain_h(3);    
          % FEM stresses
          D=matProps.D;    % This D is defined as per VEM for strain = [e11 e22 e12].
          D(3,3)=D(3,3)/4; % To avoid troubles, come back to the standard FEM definition of D 
                           % to obtain the stress by multiplying with [e11 e22 2*e12]
          svec_h=D*[eten_h(1,1);eten_h(2,2);2*eten_h(1,2)];
          if strcmp(matProps.plane_state,'plane_stress')
            sten_h=[svec_h(1),svec_h(3),0;svec_h(3),svec_h(2),0;0,0,0];
          elseif strcmp(matProps.plane_state,'plane_strain')
            sten_h=[svec_h(1),svec_h(3),0;svec_h(3),svec_h(2),0;0,0,matProps.nu*(svec_h(1)+svec_h(2))];      
          else
            throw_error('Error in plot_numerical_solution_linelast2d.m --> plot_stress_and_strain --> fem_stress_and_strain_linelast2d: plane_state');      
          end
          pstress_h = sort(eig(sten_h),'descend');  % principal stresses: pstress1 > pstress2 > pstress3  
          sxx_h=sten_h(1,1); syy_h=sten_h(2,2); szz_h=sten_h(3,3);
          sxy_h=sten_h(1,2); sxz_h=0.0; syz_h=0.0;
          VM_h=(1/sqrt(2))*sqrt((sxx_h-syy_h)^2+(syy_h-szz_h)^2+(szz_h-sxx_h)^2+6*(sxy_h^2+syz_h^2+sxz_h^2));    
          stress.s11(gp)=sten_h(1,1);
          stress.s12(gp)=sten_h(1,2);   
          stress.s22(gp)=sten_h(2,2); 
          stress.s33(gp)=sten_h(3,3);   
          stress.s1(gp)=pstress_h(1);
          stress.s2(gp)=pstress_h(2);
          stress.s3(gp)=pstress_h(3);  
          stress.vm(gp)=VM_h;   
          % assign stress and strain to the element Gauss points             
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
    throw_error('Error in plot_numerical_solution_linelast2d.m --> plot_stress_and_strain --> fem_stress_and_strain_linelast2d: vemlab_method');
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
    mytitle='$\sigma_x^h$';
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
    mytitle='$\sigma_{xy}^h$';
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
    mytitle='$\sigma_{y}^h$';
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
    mytitle='$\sigma_z^h$';
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
    mytitle='$\sigma_1^h$';
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
    mytitle='$\sigma_2^h$';
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
    mytitle='$\sigma_3^h$';
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
    mytitle='$\sigma_{\mathrm{vm}}^h$';
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
    mytitle='$\varepsilon_x^h$';
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
    mytitle='$\varepsilon_{xy}^h$';
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
    mytitle='$\varepsilon_y^h$';
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
    mytitle='$\varepsilon_z^h$';
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
    mytitle='$\varepsilon_1^h$';
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
    mytitle='$\varepsilon_2^h$';
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
    mytitle='$\varepsilon_3^h$';
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

