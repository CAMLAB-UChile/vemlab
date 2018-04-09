function vem_norms_poisson2d(exact_solution_handle,uh_global,domainMesh,config,matProps)  
  % initialize some variables
  h_max=0; L2_norm_top=0; L2_norm_bottom=0; 
  H1_seminorm_top=0; H1_seminorm_bottom=0; 
  
  % mesh data
  coords=domainMesh.coords; 
  connect=domainMesh.connect; 
  num_elem=length(connect);   
  
  % loop over elements
  for i = 1:num_elem
    % separate dofs
    nodes=connect{i}; elem_coord=coords(nodes,:);    
    dofs=nodes;
    
    % element VEM nodal solution
    uh_elem_column=uh_global(dofs);
    
    % h = nodal spacing = maximum distance between any two neighbor nodes in the mesh
    verts=elem_coord;
    h_size=max_edge_size(verts);
    if h_size>h_max
      h_max=h_size;
    end     
    
    % VEM matrix
    Wc=vem_Wc_poisson2d(verts);
    
    % VEM gradient
    pic_grad_uh=Wc'*uh_elem_column; 
    
    % Gauss points generation by triangulation of the polygonal element
    gauss_order=2;
    plot_triangulation='no';
    [gauss_points,gauss_weights]=...
            gauss_triangulate(elem_coord,nodes,gauss_order,plot_triangulation);
          
    % loop over element Gauss points
    for ngp=1:length(gauss_weights)
      wgp=gauss_weights(ngp);
      xgp=gauss_points(ngp,:); % integration point in the form [x y z]  
      
      % exact solution at Gauss point      
      [u_exact_xgp,~,gradient_exact_xgp]=...
                    exact_solutions_poisson2d(exact_solution_handle,xgp,matProps);
      grad_exact_xgp=[gradient_exact_xgp.x;gradient_exact_xgp.y]; % [dudx; dudy]  
                      
      % displacement projection at Gauss point
      x_bar=(1/length(verts))*[sum(verts(:,1)),sum(verts(:,2))];
      u_bar=(1/length(verts))*sum(uh_elem_column);
      pip_uh_xgp=[pic_grad_uh(1),pic_grad_uh(2)]*(xgp'-x_bar')+u_bar;
      
      % errors at Gauss point
      u_error_xgp=pip_uh_xgp-u_exact_xgp;
      grad_error_xgp=pic_grad_uh-grad_exact_xgp;
      
      % L2-norm of the error
      L2_norm_top=L2_norm_top+dot(u_error_xgp,u_error_xgp)*wgp;
      L2_norm_bottom=L2_norm_bottom+dot(u_exact_xgp,u_exact_xgp)*wgp;     
      
      % H1-seminorm of the error       
      H1_seminorm_top=H1_seminorm_top+...
                     dot(grad_error_xgp,grad_error_xgp)*wgp;
      H1_seminorm_bottom=H1_seminorm_bottom+...
                     dot(grad_exact_xgp,grad_exact_xgp)*wgp;       
    end 
  end
  
  % Finish computation of the norms of the errors
  L2rel = sqrt(L2_norm_top/L2_norm_bottom);
  L2abs = sqrt(L2_norm_top);  
  H1rel = sqrt(H1_seminorm_top/H1_seminorm_bottom);
  H1abs = sqrt(H1_seminorm_top);
  
  % Print out norms
  fprintf('Relative L2-norm of the error = %.20f\n',L2rel); 
  fprintf('Relative H1-seminorm of the error = %.20f\n',H1rel);  
  fprintf('Absolute L2-norm of the error = %.20f\n',L2abs); 
  fprintf('Absolute H1-seminorm of the error = %.20f\n',H1abs);  
  fprintf('Mesh spacing (h) = %.12f\n',h_max);    
end