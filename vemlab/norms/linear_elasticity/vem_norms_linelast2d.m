function [h_size,L2rel,H1rel,L2prel]=vem_norms_linelast2d(exact_sol,uh_global,domainMesh,config,matProps)  
  % initialize some variables
%   h_max=0; 
  h_size=0;  
  L2_norm_top=0; L2_norm_bottom=0;
  L2_pressure_norm_top=0; L2_pressure_norm_bottom=0;  
  H1_seminorm_top=0; H1_seminorm_bottom=0; 
  
  % mesh data
  coords=domainMesh.coords; 
  connect=domainMesh.connect; 
  num_elem=length(connect);
  
  % loop over elements
  for i = 1:num_elem
    % separate dofs
    nodes=connect{i}; elem_coord=coords(nodes,:);    
    dofx=zeros(length(nodes),1); dofy=zeros(length(nodes),1);
    dofs=zeros(2*length(nodes),1);
    for node_i=1:length(nodes)
      dofx(node_i)=2*nodes(node_i)-1;
      dofs(2*node_i-1)=2*nodes(node_i)-1;
      dofy(node_i)=2*nodes(node_i);
      dofs(2*node_i)=2*nodes(node_i);
    end
    
    % element VEM nodal solution
    uh_elem_array=[uh_global(dofx),uh_global(dofy)];
    uh_elem_column=uh_global(dofs);
    
    % h = nodal spacing = maximum distance between any two neighbor nodes in the mesh
    verts=elem_coord;
%     h_size=max_edge_size(verts);
%     if h_size>h_max
%       h_max=h_size;
%     end  
    h_E=polyarea(verts(:,1),verts(:,2));     
    h_size=h_size+h_E;
    
    % VEM matrix
    Wc=vem_Wc_linelast2d(verts);
    Q=zeros(2*length(nodes),2);
    range1=1:2:(size(Wc,1)-1); range2=2:2:size(Wc,1);
    q2=Wc(range1,3); q1=Wc(range2,3);
    Q(range1,:)=[q2,-q2]; 
    Q(range2,:)=[-q1,q1];
    
    % VEM strains
    pic_grad_uh=Wc'*uh_elem_column; 
    pir_grad_uh=Q'*uh_elem_column;
    
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
      [u_exact_xgp,p_exact_xgp,strainvec_exact_xgp]=...
                    exact_solutions_linelast2d(exact_sol,xgp,matProps);
      strain_exact_xgp=strainvec_exact_xgp';            
                  
      % displacement projection at Gauss point
      x_bar=(1/length(verts))*[sum(verts(:,1)),sum(verts(:,2))];
      u_bar=(1/length(verts))*[sum(uh_elem_array(:,1)),sum(uh_elem_array(:,2))];
      pip_uh_xgp=...
        [pic_grad_uh(1),pic_grad_uh(3);pic_grad_uh(3),pic_grad_uh(2)]*(xgp'-x_bar')...
        +[0,pir_grad_uh(1);pir_grad_uh(2),0]*(xgp'-x_bar')+u_bar';
      
      % pressure at Gauss point
      stress_uh_xgp=matProps.D*pic_grad_uh;
      if strcmp(matProps.plane_state,'plane_stress')
        ph_xgp=-1/3*(stress_uh_xgp(1)+stress_uh_xgp(2));
      elseif strcmp(matProps.plane_state,'plane_strain')||strcmp(matProps.plane_state,'plane_strain_stokes')
        ph_xgp=-1/3*(stress_uh_xgp(1) + stress_uh_xgp(2) + matProps.nu*(stress_uh_xgp(1)+stress_uh_xgp(2)));           
      end

      % errors at Gauss point
      u_error_xgp=pip_uh_xgp-u_exact_xgp;
      p_error_xgp=ph_xgp-p_exact_xgp;      
      strain_error_xgp=pic_grad_uh-strain_exact_xgp;
      
      % L2-norm of the error
      L2_norm_top=L2_norm_top+dot(u_error_xgp,u_error_xgp)*wgp;
      L2_norm_bottom=L2_norm_bottom+dot(u_exact_xgp,u_exact_xgp)*wgp;     
      
      % L2-norm of the pressure error
      L2_pressure_norm_top=L2_pressure_norm_top+dot(p_error_xgp,p_error_xgp)*wgp;      
      L2_pressure_norm_bottom=L2_pressure_norm_bottom+dot(p_exact_xgp,p_exact_xgp)*wgp;      
      
      % H1-seminorm of the error       
      H1_seminorm_top=H1_seminorm_top+...
                     dot(strain_error_xgp,(matProps.D)*strain_error_xgp)*wgp;
      H1_seminorm_bottom=H1_seminorm_bottom+...
                     dot(strain_exact_xgp,(matProps.D)*strain_exact_xgp)*wgp;       
    end 
  end
  
  h_size=sqrt(h_size/num_elem); % square root of the average area of the elements
  
  % Finish computation of the norms of the errors
  L2rel = sqrt(L2_norm_top/L2_norm_bottom);
%   L2abs = sqrt(L2_norm_top);  
  L2prel = sqrt(L2_pressure_norm_top/L2_pressure_norm_bottom);
%   L2pabs = sqrt(L2_pressure_norm_top);  
  H1rel = sqrt(H1_seminorm_top/H1_seminorm_bottom);
%   H1abs = sqrt(H1_seminorm_top);
  
  % Print out norms
  fprintf('Relative L2-norm of the displacement error = %.20f\n',L2rel); 
  fprintf('Relative H1-seminorm of the displacement error = %.20f\n',H1rel);  
  fprintf('Relative L2-norm of the pressure error = %.20f\n',L2prel);  
%   fprintf('Absolute L2-norm of the displacement error = %.20f\n',L2abs); 
%   fprintf('Absolute H1-seminorm of the displacement error = %.20f\n',H1abs);  
%   fprintf('Absolute L2-norm of the pressure error = %.20f\n',L2pabs);  
%   fprintf('Mesh spacing (h) = %.12f\n',h_max);    
  fprintf('Element size (sqrt(average area of the elements)): h = %.12f\n\n',h_size);   
end