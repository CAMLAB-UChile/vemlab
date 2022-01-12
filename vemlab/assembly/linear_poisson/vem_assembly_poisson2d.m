%-----------------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Jan. 31, 2020: add a check on matProps to figure out if comes on an
%                element-by-element fashion or not.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [K_global,f_global] = vem_assembly_poisson2d(domainMesh,matProps,source_term_fun_values)
  num_nodes=length(domainMesh.coords(:,1)); 
  numel=length(domainMesh.connect);  
  K_global=zeros(num_nodes,num_nodes);
  f_global=zeros(num_nodes,1);
  size_k=length(matProps.k);
  if size_k==numel
    k_is_unique = false;
  elseif size_k==1
    k_is_unique = true;
  end
  for e=1:numel 
    nodes=domainMesh.connect{e,1};
    verts=domainMesh.coords(nodes,:);
    area=polyarea(verts(:,1),verts(:,2));
%     fprintf('Element %d --> Area = %f\n',e,area);
    % conductivity for isotropic material 
    if k_is_unique    
      mp.k=matProps.k; % conductivity is the same for all the elements
    else
      mp.k=matProps.k{e,1}; % conductivity is particular for the current element
    end
    size_k_e=length(mp.k);    
    % assembly element stiffness matrix and element force vector 
    if size_k_e == 1
      if mp.k>0 % only assembly if the element has a conductivity > 0 (otherwise, it is like a void --- cannot conduct)
        K_global(nodes,nodes)=K_global(nodes,nodes)+vem_stiffness_poisson2d(verts,area,mp.k);
        f_global(nodes)=f_global(nodes)+vem_source_vector_poisson2d(verts,area,source_term_fun_values,mp);         
      end
    else
      throw_error('In vem_assembly_poisson2d.m: either conductivity was not defined or multiple conductivities assigned to the element... element stiffness not implemented for this condition');
    end
  end
end

