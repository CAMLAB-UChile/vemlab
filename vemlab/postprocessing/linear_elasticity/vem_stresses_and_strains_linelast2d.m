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
% Compute VEM stresses and strains.
%
% Usage
% =====
% [stresses,strains] = ...
%      vem_stresses_and_strains_linelast2d(uh_global,mesh,matProps)
%
% Input
% =====
% uh_global  : nodal solution (entire mesh)
% mesh       : structure containing the polygonal mesh information
% natProps   : material properties
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
       vem_stresses_and_strains_linelast2d(uh_global,mesh,matProps,config) 
  fprintf('\n');
  fprintf('Computing %s stresses and strains...\n',config.vemlab_method);
  % mesh data
  coords=mesh.coords;
  connect=mesh.connect;    
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
    uh_elem_column=uh_global(dofs);
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
      eten=[pic_grad_uh(1),pic_grad_uh(3),0;...
            pic_grad_uh(3),pic_grad_uh(2),0;...
            0,0,-matProps.nu*(pic_grad_uh(1)+pic_grad_uh(2))];  
    elseif strcmp(matProps.plane_state,'plane_strain')
      eten=[pic_grad_uh(1),pic_grad_uh(3),0;...
            pic_grad_uh(3),pic_grad_uh(2),0;...
            0,0,0];      
    else
      throw_error('Error in vem_stresses_and_strains_linelast.m: plane_state\n');
    end
    pstrains=sort(eig(eten),'descend'); % principal strain: pstrain1 > pstrain2 > pstrain3      
    strains.e11(i)=eten(1,1);
    strains.e12(i)=eten(1,2);   
    strains.e22(i)=eten(2,2); 
    strains.e33(i)=eten(3,3);  
    strains.e1(i)=pstrains(1);
    strains.e2(i)=pstrains(2);
    strains.e3(i)=pstrains(3);    
    % VEM stresses
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
end

