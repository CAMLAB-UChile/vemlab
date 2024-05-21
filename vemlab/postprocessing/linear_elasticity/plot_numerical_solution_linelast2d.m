%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                     VEMLAB
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
% [triangles_per_polygon,stress,strain] = plot_numerical_solution_linelast2d(domainMesh,solution,...
%                                                                            matProps,config)
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
% Jan. 04, 2021: plotting functionalities improved
% Feb. 1, 2020: add a return array variable called triangles_per_polygon, which
%               is used to fix an error in the plotting of VEM stresses and strains 
%               into a text file stage (by A. Ortiz-Bernardin)
% Oct. 20, 2018: add option to switch off all matlab contour plots (by A. Ortiz-Bernardin)
% Sept. 16, 2018: add option to plot deformed domain for linelast2d (by A. Ortiz-Bernardin)
% Apr. 19, 2018: improve the plotting of axis and fonts (by A. Ortiz-Bernardin)
% Mar. 17, 2018: first realease (by A. Ortiz-Bernardin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [stress,strain] = ...
                   plot_numerical_solution_linelast2d(domainMesh,solution,matProps,config)
  [stress,strain]=plot_numerical_solution_vem_fem_linelast2d(solution,domainMesh,matProps,config);
end

function [stress,strain]=plot_numerical_solution_vem_fem_linelast2d(solution,domainMesh,matProps,config)
  if strcmp(config.vemlab_method,'VEM2D')
    if strcmp(config.linelast2d_plot_stress_and_strain,'yes')
      [stress,strain]=...
        vem_stress_and_strain_linelast2d(solution,domainMesh,matProps,config);
    else
      stress=[]; strain=[];
    end
    if strcmp(config.create_matlab_contour_plots,'yes')
      patch_plot_VEM2D_NIVEM2D_FEM2DT3_linelast2d(domainMesh,solution,stress,strain,config);
    end      
  elseif strcmp(config.vemlab_method,'FEM2DT3')
    if strcmp(config.linelast2d_plot_stress_and_strain,'yes')
      [stress,strain]=...
        fem2DT3_stress_and_strain_linelast2d(solution,domainMesh,matProps,config);  
    else
      stress=[]; strain=[];
    end
    if strcmp(config.create_matlab_contour_plots,'yes')
      patch_plot_VEM2D_NIVEM2D_FEM2DT3_linelast2d(domainMesh,solution,stress,strain,config);
    end       
  elseif strcmp(config.vemlab_method,'FEM2DQ4') 
    if strcmp(config.linelast2d_plot_stress_and_strain,'yes')
      [stress,strain,gp_list,h_min,xmin,xmax,ymin,ymax]=...
        fem2DQ4_stress_and_strain_linelast2d(solution,domainMesh,matProps,config);
    else
      stress=[]; strain=[];
    end
    if strcmp(config.create_matlab_contour_plots,'yes')
      patch_plot_FEM2DQ4_linelast2d(domainMesh,solution,stress,strain,gp_list,h_min,...
                                    xmin,xmax,ymin,ymax,config);
    end
  else
    throw_error('Error in plot_numerical_solution_linelast2d.m --> plot_numerical_solution_vem_fem_linelast2d: vemlab_method\n');
  end
end

function [stress,strain] = ...
       vem_stress_and_strain_linelast2d(solution,domainMesh,matProps,config) 
     
  fprintf('Postprocessing %s stresses and strains...\n',config.vemlab_method); 
  
  mu=matProps.Ey/(2*(1+matProps.nu));
  kappa=matProps.Ey/(3*(1-2*matProps.nu));
  lambda=kappa-2*mu/3;  
  
  % mesh data
  coords=domainMesh.coords;
  connect=domainMesh.connect;    
  num_elem=length(connect); 
  % loop over elements  
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
    elseif strcmp(matProps.plane_state,'plane_strain')||strcmp(matProps.plane_state,'plane_strain_stokes')
      eten_h=[pic_grad_uh(1),pic_grad_uh(3),0;...
            pic_grad_uh(3),pic_grad_uh(2),0;...
            0,0,0];      
    else
      throw_error('Error in plot_numerical_solution_linelast2d.m --> plot_stress_and_strain --> vem_stress_and_strain_linelast2d: plane_state\n');
    end
    pstrain_h=sort(eig(eten_h),'descend'); % principal strain: pstrain1 > pstrain2 > pstrain3      
    strain.e11(i)=eten_h(1,1);
    strain.e12(i)=eten_h(1,2);   
    strain.e22(i)=eten_h(2,2); 
    strain.e33(i)=eten_h(3,3);  
    strain.e1(i)=pstrain_h(1);
    strain.e2(i)=pstrain_h(2);
    strain.e3(i)=pstrain_h(3);    
    % VEM stresses
    D=matProps.D;    % This D is defined as per VEM for strain = [e11 e22 e12].
    D(3,3)=D(3,3)/4; % To avoid troubles, come back to the standard FEM definition of D 
                     % to obtain the stress by multiplying with [e11 e22 2*e12]    
    svec_h=D*[eten_h(1,1);eten_h(2,2);2*eten_h(1,2)];
    if strcmp(matProps.plane_state,'plane_stress')
      sten_h=[svec_h(1),svec_h(3),0;svec_h(3),svec_h(2),0;0,0,0];
    elseif strcmp(matProps.plane_state,'plane_strain')||strcmp(matProps.plane_state,'plane_strain_stokes')
      sten_h=[svec_h(1),svec_h(3),0;svec_h(3),svec_h(2),0;0,0,matProps.nu*(svec_h(1)+svec_h(2))];      
    else
      throw_error('Error in plot_numerical_solution_linelast2d.m --> plot_stress_and_strain --> vem_stress_and_strain_linelast2d: plane_state\n');      
    end
    pstress_h = sort(eig(sten_h),'descend');  % principal stresses: pstress1 > pstress2 > pstress3  
    sxx_h=sten_h(1,1); syy_h=sten_h(2,2); szz_h=sten_h(3,3);
    sxy_h=sten_h(1,2); sxz_h=0.0; syz_h=0.0;
    VM_h=(1/sqrt(2))*sqrt((sxx_h-syy_h)^2+(syy_h-szz_h)^2+(szz_h-sxx_h)^2+6*(sxy_h^2+syz_h^2+sxz_h^2));
    % p_h=-lambda*(eten_h(1,1)+eten_h(2,2));
    % p_h=lambda*(eten_h(1,1)+eten_h(2,2));
    p_h=-1/3*(sten_h(1,1)+sten_h(2,2)+sten_h(3,3));  
    stress.s11(i)=sten_h(1,1);
    stress.s12(i)=sten_h(1,2);   
    stress.s22(i)=sten_h(2,2); 
    stress.s33(i)=sten_h(3,3);   
    stress.s1(i)=pstress_h(1);
    stress.s2(i)=pstress_h(2);
    stress.s3(i)=pstress_h(3);  
    stress.vm(i)=VM_h;  
    stress.p(i)=p_h;     
  end    
end

function [stress,strain] = ...
           fem2DT3_stress_and_strain_linelast2d(solution,domainMesh,matProps,config)
  
  fprintf('Postprocessing %s stresses and strains...\n',config.vemlab_method);  
  
  mu=matProps.Ey/(2*(1+matProps.nu));
  kappa=matProps.Ey/(3*(1-2*matProps.nu));
  lambda=kappa-2*mu/3;
  
  % mesh data
  coords=domainMesh.coords;
  connect=domainMesh.connect;    
  num_elem=length(connect); 

  % loop over elements
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
    elseif strcmp(matProps.plane_state,'plane_strain')||strcmp(matProps.plane_state,'plane_strain_stokes')
      eten_h=[str_h(1),str_h(3)/2,0;...
            str_h(3)/2,str_h(2),0;...
            0,0,0];      
    else
      throw_error('Error in plot_numerical_solution_linelast2d.m --> plot_stress_and_strain --> fem_stress_and_strain_linelast2d: plane_state\n');
    end
    pstrain_h=sort(eig(eten_h),'descend'); % principal strain: pstrain1 > pstrain2 > pstrain3      
    strain.e11(i)=eten_h(1,1);
    strain.e12(i)=eten_h(1,2);   
    strain.e22(i)=eten_h(2,2); 
    strain.e33(i)=eten_h(3,3);  
    strain.e1(i)=pstrain_h(1);
    strain.e2(i)=pstrain_h(2);
    strain.e3(i)=pstrain_h(3);    
    % FEM stresses
    D=matProps.D;    % This D is defined as per VEM for strain = [e11 e22 e12].
    D(3,3)=D(3,3)/4; % To avoid troubles, come back to the standard FEM definition of D 
                     % to obtain the stress by multiplying with [e11 e22 2*e12]
    svec_h=D*[eten_h(1,1);eten_h(2,2);2*eten_h(1,2)];
    if strcmp(matProps.plane_state,'plane_stress')
      sten_h=[svec_h(1),svec_h(3),0;svec_h(3),svec_h(2),0;0,0,0];
    elseif strcmp(matProps.plane_state,'plane_strain')||strcmp(matProps.plane_state,'plane_strain_stokes')
      sten_h=[svec_h(1),svec_h(3),0;svec_h(3),svec_h(2),0;0,0,matProps.nu*(svec_h(1)+svec_h(2))];      
    else
      throw_error('Error in plot_numerical_solution_linelast2d.m --> plot_stress_and_strain --> fem_stress_and_strain_linelast2d: plane_state\n');      
    end
    pstress_h = sort(eig(sten_h),'descend');  % principal stresses: pstress1 > pstress2 > pstress3  
    sxx_h=sten_h(1,1); syy_h=sten_h(2,2); szz_h=sten_h(3,3);
    sxy_h=sten_h(1,2); sxz_h=0.0; syz_h=0.0;
    VM_h=(1/sqrt(2))*sqrt((sxx_h-syy_h)^2+(syy_h-szz_h)^2+(szz_h-sxx_h)^2+6*(sxy_h^2+syz_h^2+sxz_h^2));  
    % p_h=-lambda*(eten_h(1,1)+eten_h(2,2));
    p_h=-1/3*(sten_h(1,1)+sten_h(2,2)+sten_h(3,3)); 
    stress.s11(i)=sten_h(1,1);
    stress.s12(i)=sten_h(1,2);   
    stress.s22(i)=sten_h(2,2); 
    stress.s33(i)=sten_h(3,3);   
    stress.s1(i)=pstress_h(1);
    stress.s2(i)=pstress_h(2);
    stress.s3(i)=pstress_h(3);  
    stress.vm(i)=VM_h;
    stress.p(i)=p_h;        
  end
end

function [stress,strain,gp_list,h_min,xmin,xmax,ymin,ymax] = ...
           fem2DQ4_stress_and_strain_linelast2d(solution,domainMesh,matProps,config)
  
  fprintf('Postprocessing %s stresses and strains at Gauss points...\n',...
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
        elseif strcmp(matProps.plane_state,'plane_strain')||strcmp(matProps.plane_state,'plane_strain_stokes')
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
        elseif strcmp(matProps.plane_state,'plane_strain')||strcmp(matProps.plane_state,'plane_strain_stokes')
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
        Ey=matProps.Ey;
        nu=matProps.nu;
        lam=Ey*nu/((1+nu)*(1-2*nu));
        % stress.p(gp)=-lam*(eten_h(1,1)+eten_h(2,2));  
        stress.p(gp)=-1/3*(sten_h(1,1)+sten_h(2,2)+sten_h(3,3));             
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
    %h_size=max_edge_size(elem_coord);
    h_size=polyarea(elem_coord(:,1),elem_coord(:,2));
    if h_size<h_min
      h_min=h_size;
    end        
  end
  h_min=sqrt(h_min); % square root of the minumum element's area among the elements of the mesh    
end

function patch_plot_VEM2D_NIVEM2D_FEM2DT3_linelast2d(domainMesh,solution,stress,strain,config)  

  solutionType=config.vemlab_method;
  nodes=domainMesh.coords;
  nodesNumber=size(nodes,1);  
  polygons=domainMesh.connect;
  range=1:nodesNumber;
  displacementsX=solution(2*range-1);
  displacementsY=solution(2*range);  
  displacementsNorm=sqrt(displacementsX.*displacementsX+displacementsY.*displacementsY);

  fprintf('\n'); 
  fprintf('Plotting %s displacement solution...\n',solutionType);  
  
  if strcmp(config.linelast2d_plot_deformed_domain,'yes') 
    scale=config.linelast2d_scale_for_plotting_deformed_domain;
    nodes_dispx=[nodes(:,1)+displacementsX*scale, nodes(:,2)];
    nodes_dispy=[nodes(:,1), nodes(:,2)+displacementsY*scale];
    nodes=[nodes(:,1)+displacementsX*scale, nodes(:,2)+displacementsY*scale];
    if strcmp(config.linelast2d_plot_displacement.unorm,'yes')
      mytitleclb='$\scriptstyle{\|}\displaystyle{\mathbf{u}_h}\scriptstyle{\|}$';
      mytitlefigfile='uh';   
      colorbar_limits=config.linelast2d_plot_displacement.clim.unorm;       
      patch_plot_figure(config,'on_nodes',polygons,nodes,displacementsNorm,mytitleclb,mytitlefigfile,colorbar_limits,1.09,'numerical');
    end

    if strcmp(config.linelast2d_plot_displacement.ux,'yes')  
      mytitleclb='$u_{1,h}$';
      mytitlefigfile='u1h';  
      colorbar_limits=config.linelast2d_plot_displacement.clim.ux;        
      patch_plot_figure(config,'on_nodes',polygons,nodes_dispx,displacementsX,mytitleclb,mytitlefigfile,colorbar_limits,1.05,'numerical');
    end

    if strcmp(config.linelast2d_plot_displacement.uy,'yes')
      mytitleclb='$u_{2,h}$';
      mytitlefigfile='u2h'; 
      colorbar_limits=config.linelast2d_plot_displacement.clim.uy;        
      patch_plot_figure(config,'on_nodes',polygons,nodes_dispy,displacementsY,mytitleclb,mytitlefigfile,colorbar_limits,1.05,'numerical');
    end    
  else
    if strcmp(config.linelast2d_plot_displacement.unorm,'yes')
      mytitleclb='$\scriptstyle{\|}\displaystyle{\mathbf{u}_h}\scriptstyle{\|}$';
      mytitlefigfile='uh';
      colorbar_limits=config.linelast2d_plot_displacement.clim.unorm;        
      patch_plot_figure(config,'on_nodes',polygons,nodes,displacementsNorm,mytitleclb,mytitlefigfile,colorbar_limits,1.09,'numerical');
    end

    if strcmp(config.linelast2d_plot_displacement.ux,'yes')  
      mytitleclb='$u_{1,h}$';
      mytitlefigfile='u1h';
      colorbar_limits=config.linelast2d_plot_displacement.clim.ux;        
      patch_plot_figure(config,'on_nodes',polygons,nodes,displacementsX,mytitleclb,mytitlefigfile,colorbar_limits,1.05,'numerical');
    end

    if strcmp(config.linelast2d_plot_displacement.uy,'yes')
      mytitleclb='$u_{2,h}$';
      mytitlefigfile='u2h'; 
      colorbar_limits=config.linelast2d_plot_displacement.clim.uy;        
      patch_plot_figure(config,'on_nodes',polygons,nodes,displacementsY,mytitleclb,mytitlefigfile,colorbar_limits,1.05,'numerical');
    end     
  end
  
  if ~isempty(stress)
    fprintf('Plotting %s stress and strain solutions...\n',solutionType);   

    if strcmp(config.linelast2d_plot_stress.s11,'yes')
      mytitleclb='$\sigma_{11,h}$';
      mytitlefigfile='s11h';  
      colorbar_limits=config.linelast2d_plot_stress.clim.s11;       
      patch_plot_figure(config,'on_faces',polygons,nodes,stress.s11',mytitleclb,mytitlefigfile,colorbar_limits,1.23,'numerical');
    end

    if strcmp(config.linelast2d_plot_stress.s22,'yes')  
      mytitleclb='$\sigma_{22,h}$';
      mytitlefigfile='s22h'; 
      colorbar_limits=config.linelast2d_plot_stress.clim.s22;       
      patch_plot_figure(config,'on_faces',polygons,nodes,stress.s22',mytitleclb,mytitlefigfile,colorbar_limits,1.23,'numerical');
    end 

    if strcmp(config.linelast2d_plot_stress.s12,'yes')  
      mytitleclb='$\sigma_{12,h}$';
      mytitlefigfile='s12h'; 
      colorbar_limits=config.linelast2d_plot_stress.clim.s12;       
      patch_plot_figure(config,'on_faces',polygons,nodes,stress.s12',mytitleclb,mytitlefigfile,colorbar_limits,1.23,'numerical');
    end       

    if strcmp(config.linelast2d_plot_stress.s33,'yes')  
      mytitleclb='$\sigma_{33,h}$';
      mytitlefigfile='s33h';
      colorbar_limits=config.linelast2d_plot_stress.clim.s33;       
      patch_plot_figure(config,'on_faces',polygons,nodes,stress.s33',mytitleclb,mytitlefigfile,colorbar_limits,1.23,'numerical');
    end  

    if strcmp(config.linelast2d_plot_stress.s1,'yes')  
      mytitleclb='$\sigma_{1,h}$';
      mytitlefigfile='s1h'; 
      colorbar_limits=config.linelast2d_plot_stress.clim.s1;       
      patch_plot_figure(config,'on_faces',polygons,nodes,stress.s1',mytitleclb,mytitlefigfile,colorbar_limits,1.02,'numerical');
    end   

    if strcmp(config.linelast2d_plot_stress.s2,'yes')  
      mytitleclb='$\sigma_{2,h}$';
      mytitlefigfile='s2h';
      colorbar_limits=config.linelast2d_plot_stress.clim.s2;       
      patch_plot_figure(config,'on_faces',polygons,nodes,stress.s2',mytitleclb,mytitlefigfile,colorbar_limits,1.02,'numerical');
    end 

    if strcmp(config.linelast2d_plot_stress.s3,'yes')  
      mytitleclb='$\sigma_{3,h}$';
      mytitlefigfile='s3h';
      colorbar_limits=config.linelast2d_plot_stress.clim.s3;       
      patch_plot_figure(config,'on_faces',polygons,nodes,stress.s3',mytitleclb,mytitlefigfile,colorbar_limits,1.02,'numerical');
    end    

    if strcmp(config.linelast2d_plot_stress.vm,'yes')  
      mytitleclb='$\sigma_{\mathrm{vm},h}$';
      mytitlefigfile='vmh'; 
      colorbar_limits=config.linelast2d_plot_stress.clim.vm;       
      patch_plot_figure(config,'on_faces',polygons,nodes,stress.vm',mytitleclb,mytitlefigfile,colorbar_limits,1.42,'numerical');
    end 

    if strcmp(config.linelast2d_plot_stress.p,'yes')  
      mytitleclb='$p_h$';
      mytitlefigfile='ph';
      colorbar_limits=config.linelast2d_plot_stress.clim.p;       
      patch_plot_figure(config,'on_faces',polygons,nodes,stress.p',mytitleclb,mytitlefigfile,colorbar_limits,0.64,'numerical');
    end  

    if strcmp(config.linelast2d_plot_strain.e11,'yes')  
      mytitleclb='$\varepsilon_{11,h}$';
      mytitlefigfile='e11h';  
      colorbar_limits=config.linelast2d_plot_strain.clim.e11;      
      patch_plot_figure(config,'on_faces',polygons,nodes,strain.e11',mytitleclb,mytitlefigfile,colorbar_limits,1.14,'numerical');
    end      

    if strcmp(config.linelast2d_plot_strain.e12,'yes')  
      mytitleclb='$\varepsilon_{12,h}$';
      mytitlefigfile='e12h';
      colorbar_limits=config.linelast2d_plot_strain.clim.e12;      
      patch_plot_figure(config,'on_faces',polygons,nodes,strain.e12',mytitleclb,mytitlefigfile,colorbar_limits,1.14,'numerical');
    end    

    if strcmp(config.linelast2d_plot_strain.e22,'yes')  
      mytitleclb='$\varepsilon_{22,h}$';
      mytitlefigfile='e22h'; 
      colorbar_limits=config.linelast2d_plot_strain.clim.e22;      
      patch_plot_figure(config,'on_faces',polygons,nodes,strain.e22',mytitleclb,mytitlefigfile,colorbar_limits,1.14,'numerical');
    end  

    if strcmp(config.linelast2d_plot_strain.e33,'yes')  
      mytitleclb='$\varepsilon_{33,h}$';
      mytitlefigfile='e33h'; 
      colorbar_limits=config.linelast2d_plot_strain.clim.e33;      
      patch_plot_figure(config,'on_faces',polygons,nodes,strain.e33',mytitleclb,mytitlefigfile,colorbar_limits,1.14,'numerical');
    end      

    if strcmp(config.linelast2d_plot_strain.e1,'yes')  
      mytitleclb='$\varepsilon_{1,h}$';
      mytitlefigfile='e1h'; 
      colorbar_limits=config.linelast2d_plot_strain.clim.e1;      
      patch_plot_figure(config,'on_faces',polygons,nodes,strain.e1',mytitleclb,mytitlefigfile,colorbar_limits,0.95,'numerical');
    end  

    if strcmp(config.linelast2d_plot_strain.e2,'yes')  
      mytitleclb='$\varepsilon_{2,h}$';
      mytitlefigfile='e2h'; 
      colorbar_limits=config.linelast2d_plot_strain.clim.e2;      
      patch_plot_figure(config,'on_faces',polygons,nodes,strain.e2',mytitleclb,mytitlefigfile,colorbar_limits,0.95,'numerical');
    end  

    if strcmp(config.linelast2d_plot_strain.e3,'yes')  
      mytitleclb='$\varepsilon_{3,h}$';
      mytitlefigfile='e3h';
      colorbar_limits=config.linelast2d_plot_strain.clim.e3;      
      patch_plot_figure(config,'on_faces',polygons,nodes,strain.e3',mytitleclb,mytitlefigfile,colorbar_limits,0.95,'numerical');
    end  
  end

end

function patch_plot_FEM2DQ4_linelast2d(domainMesh,solution,stress,strain,gp_list,...
                                       h_min,xmin,xmax,ymin,ymax,config) 
                                         
  solutionType=config.vemlab_method;
  nodes=domainMesh.coords;
  nodesNumber=size(nodes,1);  
  polygons=domainMesh.connect;
  range=1:nodesNumber;
  displacementsX=solution(2*range-1);
  displacementsY=solution(2*range);  
  displacementsNorm=sqrt(displacementsX.*displacementsX+displacementsY.*displacementsY);

  fprintf('\n'); 
  fprintf('Plotting %s displacement solution...\n',solutionType);  
  if strcmp(config.linelast2d_plot_deformed_domain,'yes') 
    nodes_disp=nodes;
    scale=config.linelast2d_scale_for_plotting_deformed_domain;
%     nodes_disp=[nodes_disp(:,1)+displacementsX*scale, nodes_disp(:,2)+displacementsY*scale];
    throw_warning('For FEM2DQ4, plot using deformed mesh is available only for displacements');
    
    if strcmp(config.linelast2d_plot_displacement.unorm,'yes')
      nodes_disp=[nodes_disp(:,1)+displacementsX*scale, nodes_disp(:,2)+displacementsY*scale];
      mytitleclb='$\scriptstyle{\|}\displaystyle{\mathbf{u}_h}\scriptstyle{\|}$';
      mytitlefigfile='uh';
      colorbar_limits=config.linelast2d_plot_displacement.clim.unorm;      
      patch_plot_figure(config,'on_nodes',polygons,nodes_disp,displacementsNorm,mytitleclb,mytitlefigfile,colorbar_limits,1.09,'numerical');
    end

    if strcmp(config.linelast2d_plot_displacement.ux,'yes') 
      nodes_disp(:,1)=nodes_disp(:,1)+displacementsX*scale;
      mytitleclb='$u_{1,h}$';  
      mytitlefigfile='u1h';  
      colorbar_limits=config.linelast2d_plot_displacement.clim.ux;       
      patch_plot_figure(config,'on_nodes',polygons,nodes_disp,displacementsX,mytitleclb,mytitlefigfile,colorbar_limits,1.05,'numerical');
    end

    if strcmp(config.linelast2d_plot_displacement.uy,'yes')
      nodes_disp(:,2)=nodes_disp(:,2)+displacementsY*scale;
      mytitleclb='$u_{2,h}$';    
      mytitlefigfile='u2h';     
      colorbar_limits=config.linelast2d_plot_displacement.clim.uy;       
      patch_plot_figure(config,'on_nodes',polygons,nodes_disp,displacementsY,mytitleclb,mytitlefigfile,colorbar_limits,1.05,'numerical');
    end      
  else
    if strcmp(config.linelast2d_plot_displacement.unorm,'yes')
      mytitleclb='$\scriptstyle{\|}\displaystyle{\mathbf{u}_h}\scriptstyle{\|}$';
      mytitlefigfile='uh';   
      colorbar_limits=config.linelast2d_plot_displacement.clim.unorm;       
      patch_plot_figure(config,'on_nodes',polygons,nodes,displacementsNorm,mytitleclb,mytitlefigfile,colorbar_limits,1.09,'numerical');
    end

    if strcmp(config.linelast2d_plot_displacement.ux,'yes')  
      mytitleclb='$u_{1,h}$';     
      mytitlefigfile='u1h';   
      colorbar_limits=config.linelast2d_plot_displacement.clim.ux;       
      patch_plot_figure(config,'on_nodes',polygons,nodes,displacementsX,mytitleclb,mytitlefigfile,colorbar_limits,1.05,'numerical');
    end

    if strcmp(config.linelast2d_plot_displacement.uy,'yes')
      mytitleclb='$u_{2,h}$';   
      mytitlefigfile='u2h'; 
      colorbar_limits=config.linelast2d_plot_displacement.clim.uy;       
      patch_plot_figure(config,'on_nodes',polygons,nodes,displacementsY,mytitleclb,mytitlefigfile,colorbar_limits,1.05,'numerical');
    end     
  end
         
  if ~isempty(stress)  
    mt=0.25; %0.1;
    dx=mt*h_min;
    dy=dx;  
    xs=xmin:dx:xmax;
    ys=ymin:dy:ymax;         
    [xq,yq]=meshgrid(xs,ys);    

    fprintf('Plotting %s stress and strain solutions...\n',solutionType);
    
    if strcmp(config.linelast2d_plot_stress.s11,'yes')
      mytitleclb='$\sigma_{11,h}$';
      mytitlefigfile='s11h';   
      colorbar_limits=config.linelast2d_plot_stress.clim.s11;      
      mesh_plot_figure_fem2dQ4(config,polygons,nodes,gp_list,xq,yq,stress.s11,mytitleclb,mytitlefigfile,colorbar_limits,1.23,'numerical');
    end

    if strcmp(config.linelast2d_plot_stress.s22,'yes')  
      mytitleclb='$\sigma_{22,h}$';
      mytitlefigfile='s22h'; 
      colorbar_limits=config.linelast2d_plot_stress.clim.s22;      
      mesh_plot_figure_fem2dQ4(config,polygons,nodes,gp_list,xq,yq,stress.s22,mytitleclb,mytitlefigfile,colorbar_limits,1.23,'numerical');
    end 

    if strcmp(config.linelast2d_plot_stress.s12,'yes')  
      mytitleclb='$\sigma_{12,h}$';
      mytitlefigfile='s12h';  
      colorbar_limits=config.linelast2d_plot_stress.clim.s12;      
      mesh_plot_figure_fem2dQ4(config,polygons,nodes,gp_list,xq,yq,stress.s12,mytitleclb,mytitlefigfile,colorbar_limits,1.23,'numerical');
    end       

    if strcmp(config.linelast2d_plot_stress.s33,'yes')  
      mytitleclb='$\sigma_{33,h}$';
      mytitlefigfile='s33h';
      colorbar_limits=config.linelast2d_plot_stress.clim.s33;      
      mesh_plot_figure_fem2dQ4(config,polygons,nodes,gp_list,xq,yq,stress.s33,mytitleclb,mytitlefigfile,colorbar_limits,1.23,'numerical');
    end  

    if strcmp(config.linelast2d_plot_stress.s1,'yes')
      mytitleclb='$\sigma_{1,h}$';
      mytitlefigfile='s1h'; 
      colorbar_limits=config.linelast2d_plot_stress.clim.s1;      
      mesh_plot_figure_fem2dQ4(config,polygons,nodes,gp_list,xq,yq,stress.s1,mytitleclb,mytitlefigfile,colorbar_limits,1.02,'numerical');
    end   

    if strcmp(config.linelast2d_plot_stress.s2,'yes')  
      mytitleclb='$\sigma_{2,h}$';
      mytitlefigfile='s2h'; 
      colorbar_limits=config.linelast2d_plot_stress.clim.s2;      
      mesh_plot_figure_fem2dQ4(config,polygons,nodes,gp_list,xq,yq,stress.s2,mytitleclb,mytitlefigfile,colorbar_limits,1.02,'numerical');    
    end 

    if strcmp(config.linelast2d_plot_stress.s3,'yes')  
      mytitleclb='$\sigma_{3,h}$';
      mytitlefigfile='s3h';
      colorbar_limits=config.linelast2d_plot_stress.clim.s3;      
      mesh_plot_figure_fem2dQ4(config,polygons,nodes,gp_list,xq,yq,stress.s3,mytitleclb,mytitlefigfile,colorbar_limits,1.02,'numerical');    
    end    

    if strcmp(config.linelast2d_plot_stress.vm,'yes')  
      mytitleclb='$\sigma_{\mathrm{vm},h}$';
      mytitlefigfile='vmh';
      colorbar_limits=config.linelast2d_plot_stress.clim.vm;      
      mesh_plot_figure_fem2dQ4(config,polygons,nodes,gp_list,xq,yq,stress.vm,mytitleclb,mytitlefigfile,colorbar_limits,1.42,'numerical');    
    end 

    if strcmp(config.linelast2d_plot_stress.p,'yes')  
      mytitleclb='$p_h$';
      mytitlefigfile='ph';
      colorbar_limits=config.linelast2d_plot_stress.clim.p;      
      mesh_plot_figure_fem2dQ4(config,polygons,nodes,gp_list,xq,yq,stress.p,mytitleclb,mytitlefigfile,colorbar_limits,0.64,'numerical');    
    end  

    if strcmp(config.linelast2d_plot_strain.e11,'yes')  
      mytitleclb='$\varepsilon_{11,h}$';
      mytitlefigfile='e11h'; 
      colorbar_limits=config.linelast2d_plot_strain.clim.e11;       
      mesh_plot_figure_fem2dQ4(config,polygons,nodes,gp_list,xq,yq,strain.e11,mytitleclb,mytitlefigfile,colorbar_limits,1.14,'numerical');    
    end      

    if strcmp(config.linelast2d_plot_strain.e12,'yes')  
      mytitleclb='$\varepsilon_{12,h}$';
      mytitlefigfile='e12h';
      colorbar_limits=config.linelast2d_plot_strain.clim.e12;       
      mesh_plot_figure_fem2dQ4(config,polygons,nodes,gp_list,xq,yq,strain.e12,mytitleclb,mytitlefigfile,colorbar_limits,1.14,'numerical');     
    end    

    if strcmp(config.linelast2d_plot_strain.e22,'yes')  
      mytitleclb='$\varepsilon_{22,h}$';
      mytitlefigfile='e22h'; 
      colorbar_limits=config.linelast2d_plot_strain.clim.e22;       
      mesh_plot_figure_fem2dQ4(config,polygons,nodes,gp_list,xq,yq,strain.e22,mytitleclb,mytitlefigfile,colorbar_limits,1.14,'numerical');     
    end  

    if strcmp(config.linelast2d_plot_strain.e33,'yes')  
      mytitleclb='$\varepsilon_{33,h}$';
      mytitlefigfile='e33h'; 
      colorbar_limits=config.linelast2d_plot_strain.clim.e33;       
      mesh_plot_figure_fem2dQ4(config,polygons,nodes,gp_list,xq,yq,strain.e33,mytitleclb,mytitlefigfile,colorbar_limits,1.14,'numerical');     
    end      

    if strcmp(config.linelast2d_plot_strain.e1,'yes')  
      mytitleclb='$\varepsilon_{1,h}$';
      mytitlefigfile='e1h';
      colorbar_limits=config.linelast2d_plot_strain.clim.e1;       
      mesh_plot_figure_fem2dQ4(config,polygons,nodes,gp_list,xq,yq,strain.e1,mytitleclb,mytitlefigfile,colorbar_limits,0.95,'numerical');     
    end  

    if strcmp(config.linelast2d_plot_strain.e2,'yes')  
      mytitleclb='$\varepsilon_{2,h}$';
      mytitlefigfile='e2h';
      colorbar_limits=config.linelast2d_plot_strain.clim.e2;       
      mesh_plot_figure_fem2dQ4(config,polygons,nodes,gp_list,xq,yq,strain.e2,mytitleclb,mytitlefigfile,colorbar_limits,0.95,'numerical');     
    end  

    if strcmp(config.linelast2d_plot_strain.e3,'yes')  
      mytitleclb='$\varepsilon_{3,h}$';
      mytitlefigfile='e3h';
      colorbar_limits=config.linelast2d_plot_strain.clim.e3;       
      mesh_plot_figure_fem2dQ4(config,polygons,nodes,gp_list,xq,yq,strain.e3,mytitleclb,mytitlefigfile,colorbar_limits,0.95,'numerical');     
    end  
  end

end