function [K_global,f_global] = fem_sparse_assembly_poisson2d(domainMesh,config,matProps,source_term_fun_values)
  num_nodes=length(domainMesh.coords(:,1));    
  numel=length(domainMesh.connect);
  
  size_k=length(matProps.k);
  if size_k==numel
    k_is_unique = false;
  elseif size_k==1
    k_is_unique = true;
  end  
  
%   K_global=sparse(num_nodes,num_nodes);
%   f_global=sparse(num_nodes,1); 

  % Lengths of buffer vectors for sparse assembly
  len_vector_K = 0; 
  len_vector_f = 0;
  % Determine exact lengths of buffer vectors
  % (Elements can have arbitrary numbers of nodes)
  for e=1:numel
      num_eldofs = length(domainMesh.connect{e,1}); % 2*length(domainMesh.connect{e,1});
      len_vector_K = len_vector_K + (num_eldofs)^2; % unraveled size of K_local
      len_vector_f = len_vector_f + num_eldofs;
  end
  % Vectors for sparse assembly of K_global
  vector_K = zeros(len_vector_K,1);
  indrow_K = zeros(len_vector_K,1);
  indcol_K = zeros(len_vector_K,1);
  % Vectors for sparse assembly of f_global
  vector_f = zeros(len_vector_f,1);
  indrow_f = zeros(len_vector_f,1);
  pos_vector_K = 0;
  pos_vector_f = 0;
  
  for e=1:numel 
    nodes=domainMesh.connect{e,1};
    verts=domainMesh.coords(nodes,:);
    ind=nodes;    
    % conductivity for isotropic material 
    if k_is_unique    
      mp.k=matProps.k; % conductivity is the same for all the elements
    else
      mp.k=matProps.k{e,1}; % conductivity is particular for the current element
    end
    size_k_e=length(mp.k);     
    % assembly element stiffness matrix and element force vector 
    if size_k_e == 1
%       K_global(nodes,nodes)=K_global(nodes,nodes)+...
%                         fem_stiffness_poisson2d(verts,k,config);
%       f_global(nodes)=f_global(nodes)+...
%                         fem_source_vector_poisson2d(verts,config,source_term_fun_values); 
        num_eldofs = length(nodes);
        K_local = fem_stiffness_poisson2d(verts,mp.k,config);
        f_local = fem_source_vector_poisson2d(verts,config,source_term_fun_values,mp);
        ind_vector_K = (1:num_eldofs^2) + pos_vector_K;
        ind_vector_f = (1:num_eldofs) + pos_vector_f;
        vector_K(ind_vector_K) = K_local(:);
        vector_f(ind_vector_f) = f_local;
        indrow_K_ = ind(:,ones(length(ind),1));
        indcol_K_ = indrow_K_';
        indrow_K(ind_vector_K) = indrow_K_(:);
        indcol_K(ind_vector_K) = indcol_K_(:);
        indrow_f(ind_vector_f) = ind;
        pos_vector_K = pos_vector_K + num_eldofs^2;
        pos_vector_f = pos_vector_f + num_eldofs;
    else
      throw_error('In fem_assembly_poisson2d.m: either conductivity was not defined or multiple conductivities assigned to the element... element stiffness not implemented for this condition');
    end    
  end
  K_global = sparse(indrow_K,indcol_K,vector_K,num_nodes,num_nodes);
  f_global = sparse(indrow_f,ones(len_vector_f,1),vector_f,num_nodes,1);  
end

