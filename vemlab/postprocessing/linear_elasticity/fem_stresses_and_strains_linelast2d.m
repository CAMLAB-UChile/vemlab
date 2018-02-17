%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                     VEMLab
%           Source code  : http://camlab.cl/research/software/vemlab/
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:         compute_VEM_stresses_and_strains_linelast2d
%
% Created by : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
% Updated by :
% Updated by :
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% Compute FEM stresses and strains.
%
% Usage
% =====
% [stresses,strains] = ...
%      fem_stresses_and_strains_linelast2d(uh_global,mesh,matProps)
%
% Input
% =====
% uh_global  : nodal solution (entire mesh)
% mesh       : structure containing the polygonal mesh information
% matProps   : material properties
% config     : structure storing VEMLab configuration options and behavior
%
% Output
% ======
% stresses,strains : structures storing VEM stresses and strains
%
%-------------------------------------------------------------------------------
% References 
% ==========
%
%-------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Dec. 26, 2017: first realease (by A. Ortiz-Bernardin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [stresses,strains] = ...
       fem_stresses_and_strains_linelast2d(uh_global,mesh,matProps,config) 
  fprintf('\n');
  fprintf('Computing %s stresses and strains...\n',config.vemlab_method);
  % mesh data
  coords=mesh.coords;
  connect=mesh.connect;    
  num_elem=length(connect); 
  if strcmp(config.vemlab_method,'FEM2DT3')
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
      uh_elem_column=uh_global(dofs); 
      % area of the element
      area=triangle_area(elem_coord);
      % deformation matrix
      x=elem_coord(:,1); y=elem_coord(:,2);  
      % B is defined for [e11,e22,2*e12]
      B=(1/(2*area))*[y(2)-y(3),0,y(3)-y(1),0,y(1)-y(2),0;...
                      0,x(3)-x(2),0,x(1)-x(3),0,x(2)-x(1);...
                      x(3)-x(2),y(2)-y(3),x(1)-x(3),y(3)-y(1),x(2)-x(1),y(1)-y(2)];      
      % FEM strains
      strain_uh=B*uh_elem_column; % In FEM, B is defined for [e11,e22,2*e12]
      if strcmp(matProps.plane_state,'plane_stress')
        eten=[strain_uh(1),strain_uh(3)/2,0;... % for tensor representation divide 
              strain_uh(3)/2,strain_uh(2),0;... % strain coming from B by 2: (2*e12)/2
              0,0,-matProps.nu*(strain_uh(1)+strain_uh(2))];  
      elseif strcmp(matProps.plane_state,'plane_strain')
        eten=[strain_uh(1),strain_uh(3)/2,0;...
              strain_uh(3)/2,strain_uh(2),0;...
              0,0,0];      
      else
        throw_error('Error in fem_stresses_and_strains_linelast.m: plane_state');
      end
      pstrains=sort(eig(eten),'descend'); % principal strain: pstrain1 > pstrain2 > pstrain3      
      strains.e11(i)=eten(1,1);
      strains.e12(i)=eten(1,2);   
      strains.e22(i)=eten(2,2); 
      strains.e33(i)=eten(3,3);  
      strains.e1(i)=pstrains(1);
      strains.e2(i)=pstrains(2);
      strains.e3(i)=pstrains(3);    
      % FEM stresses
      D=matProps.D;    % This D is defined as per VEM for strain = [e11 e22 e12].
      D(3,3)=D(3,3)/4; % To avoid troubles, come back to the standard FEM definition of D 
                       % to obtain the stress by multiplying with [e11 e22 2*e12]
      svec=D*[eten(1,1);eten(2,2);2*eten(1,2)];
      if strcmp(matProps.plane_state,'plane_stress')
        sten=[svec(1),svec(3),0;svec(3),svec(2),0;0,0,0];
      elseif strcmp(matProps.plane_state,'plane_strain')
        sten=[svec(1),svec(3),0;svec(3),svec(2),0;0,0,matProps.nu*(svec(1)+svec(2))];      
      else
        throw_error('Error in fem_stresses_and_strains_linelast2d.m: plane_state');      
      end
      pstresses = sort(eig(sten),'descend');  % principal stresses: pstress1 > pstress2 > pstress3  
      sxx=sten(1,1); syy=sten(2,2); szz=sten(3,3);
      sxy=sten(1,2); sxz=0.0; syz=0.0;
      VM=(1/sqrt(2))*sqrt((sxx-syy)^2+(syy-szz)^2+(szz-sxx)^2+6*(sxy^2+syz^2+sxz^2));    
      stresses.s11(i)=sten(1,1);
      stresses.s12(i)=sten(1,2);   
      stresses.s22(i)=sten(2,2); 
      stresses.s33(i)=sten(3,3);   
      stresses.s1(i)=pstresses(1);
      stresses.s2(i)=pstresses(2);
      stresses.s3(i)=pstresses(3);  
      stresses.vm(i)=VM; 
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
      uh_elem_column=uh_global(dofs); 
      % compute strains and stresses at Gauss points
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
          strain_uh=B*uh_elem_column; % In FEM, B is defined for [e11,e22,2*e12]
          if strcmp(matProps.plane_state,'plane_stress')
            eten=[strain_uh(1),strain_uh(3)/2,0;... % for tensor representation divide 
                  strain_uh(3)/2,strain_uh(2),0;... % strain coming from B by 2: (2*e12)/2
                  0,0,-matProps.nu*(strain_uh(1)+strain_uh(2))];  
          elseif strcmp(matProps.plane_state,'plane_strain')
            eten=[strain_uh(1),strain_uh(3)/2,0;...
                  strain_uh(3)/2,strain_uh(2),0;...
                  0,0,0];      
          else
            throw_error('Error in fem_stresses_and_strains_linelast.m: plane_state');
          end    
          pstrains=sort(eig(eten),'descend'); % principal strain: pstrain1 > pstrain2 > pstrain3      
          strains.e11(gp)=eten(1,1);
          strains.e12(gp)=eten(1,2);   
          strains.e22(gp)=eten(2,2); 
          strains.e33(gp)=eten(3,3);  
          strains.e1(gp)=pstrains(1);
          strains.e2(gp)=pstrains(2);
          strains.e3(gp)=pstrains(3);    
          % FEM stresses
          D=matProps.D;    % This D is defined as per VEM for strain = [e11 e22 e12].
          D(3,3)=D(3,3)/4; % To avoid troubles, come back to the standard FEM definition of D 
                           % to obtain the stress by multiplying with [e11 e22 2*e12]
          svec=D*[eten(1,1);eten(2,2);2*eten(1,2)];
          if strcmp(matProps.plane_state,'plane_stress')
            sten=[svec(1),svec(3),0;svec(3),svec(2),0;0,0,0];
          elseif strcmp(matProps.plane_state,'plane_strain')
            sten=[svec(1),svec(3),0;svec(3),svec(2),0;0,0,matProps.nu*(svec(1)+svec(2))];      
          else
            throw_error('Error in fem_stresses_and_strains_linelast2d.m: plane_state');      
          end
          pstresses = sort(eig(sten),'descend');  % principal stresses: pstress1 > pstress2 > pstress3  
          sxx=sten(1,1); syy=sten(2,2); szz=sten(3,3);
          sxy=sten(1,2); sxz=0.0; syz=0.0;
          VM=(1/sqrt(2))*sqrt((sxx-syy)^2+(syy-szz)^2+(szz-sxx)^2+6*(sxy^2+syz^2+sxz^2));    
          stresses.s11(gp)=sten(1,1);
          stresses.s12(gp)=sten(1,2);   
          stresses.s22(gp)=sten(2,2); 
          stresses.s33(gp)=sten(3,3);   
          stresses.s1(gp)=pstresses(1);
          stresses.s2(gp)=pstresses(2);
          stresses.s3(gp)=pstresses(3);  
          stresses.vm(gp)=VM;   
        end
      end
    end
  else
    throw_error('Error in fem_stiffness_linelast2d.m: vemlab_method');
  end
end

