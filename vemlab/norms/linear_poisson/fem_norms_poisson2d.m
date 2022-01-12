function [h_size,L2rel,H1rel,L2prel]=fem_norms_poisson2d(exact_solution_handle,uh_global,domainMesh,config,matProps)
  % initialize some variables
%   h_max=0; 
  h_size=0;   
  L2_norm_top=0; L2_norm_bottom=0; 
  H1_seminorm_top=0; H1_seminorm_bottom=0; 
  
  % mesh data
  coords=domainMesh.coords; 
  connect=domainMesh.connect; 
  num_elem=length(connect);  
  
  if strcmp(config.vemlab_method,'FEM2DT3')
    
    % loop over elements
    for i = 1:num_elem
      
      % element dofs
      nodes=connect{i}; 
      elem_coord=coords(nodes,:);    
      dofs=nodes;

      % element VEM nodal solution
      uh_elem_column=uh_global(dofs);

      % h = nodal spacing = maximum distance between any two neighbor nodes in the mesh
%       h_size=max_edge_size(elem_coord);
%       if h_size>h_max
%         h_max=h_size;
%       end 
      h_E=polyarea(elem_coord(:,1),elem_coord(:,2));     
      h_size=h_size+h_E;

      % area of the element
      area=triangle_area(elem_coord);
 
      % deformation matrix
      x=elem_coord(:,1); y=elem_coord(:,2);  
      % In FEM, B is defined for [e11,e22,2*e12]
      B=(1/(2*area))*[y(2)-y(3),y(3)-y(1),y(1)-y(2);...
                      x(3)-x(2),x(1)-x(3),x(2)-x(1)];  

      % FEM strains
      grad_uh=B*uh_elem_column; 

      % Gauss points generation by triangulation of the polygonal element
      gauss_order=2;
      [gauss_points,gauss_weights]=gauss_points_T3(gauss_order,elem_coord); 

      % loop over element Gauss points
      for ngp=1:length(gauss_weights)
        wgp=gauss_weights(ngp);
        xgp=gauss_points(ngp,:); % integration point in the form [x y z]            
        [u_exact_xgp,~,gradient_exact_xgp]=...
                      exact_solutions_poisson2d(exact_solution_handle,xgp,matProps);

        % exact solution at Gauss point
        grad_exact_xgp=[gradient_exact_xgp.x;gradient_exact_xgp.y]; % [dudx; dudy]  

        % scalar field solution at Gauss point
        N1=(1/(2*area))*((x(2)*y(3)-x(3)*y(2))+(y(2)-y(3))*xgp(1)+(x(3)-x(2))*xgp(2));
        N2=(1/(2*area))*((x(3)*y(1)-x(1)*y(3))+(y(3)-y(1))*xgp(1)+(x(1)-x(3))*xgp(2));
        N3=(1/(2*area))*((x(1)*y(2)-x(2)*y(1))+(y(1)-y(2))*xgp(1)+(x(2)-x(1))*xgp(2));
        N=[N1,N2,N3];
        uh_xgp=N*uh_elem_column;

        % errors at Gauss point
        u_error_xgp=uh_xgp-u_exact_xgp;
        grad_error_xgp=grad_uh-grad_exact_xgp;

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
    
  elseif strcmp(config.vemlab_method,'FEM2DQ4')
    
    % loop over elements
    for i = 1:num_elem
      
      % element dofs
      nodes=connect{i}; 
      elem_coord=coords(nodes,:);    
      dofs=nodes;

      % element VEM nodal solution
      uh_elem_column=uh_global(dofs);

      % h = nodal spacing = maximum distance between any two neighbor nodes in the mesh
%       h_size=max_edge_size(elem_coord);
%       if h_size>h_max
%         h_max=h_size;
%       end  
      h_E=polyarea(elem_coord(:,1),elem_coord(:,2));     
      h_size=h_size+h_E;
                             
      % compute norms by numerical integration
      num_gp=2; % number of Gauss points per axis of the element
      xi=gauss_points_1d(num_gp);
      eta=gauss_points_1d(num_gp);
      wxi=gauss_weights_1d(num_gp);
      weta=gauss_weights_1d(num_gp);
      for gpxi=1:length(xi)
        dN1deta=-(1-xi(gpxi))/4;
        dN2deta=-(1+xi(gpxi))/4;
        dN3deta=(1+xi(gpxi))/4;
        dN4deta=(1-xi(gpxi))/4;      
        for gpeta=1:length(eta)
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
          B=[dN1dx,dN2dx,dN3dx,dN4dx;...
             dN1dy,dN2dy,dN3dy,dN4dy];   
          % FEM strains
          grad_uh_xgp=B*uh_elem_column; 
          % transform Guass point to Cartesian coordinate             
          N1=(1-xi(gpxi))*(1-eta(gpeta))/4;
          N2=(1+xi(gpxi))*(1-eta(gpeta))/4;
          N3=(1+xi(gpxi))*(1+eta(gpeta))/4;
          N4=(1-xi(gpxi))*(1+eta(gpeta))/4;
          Nv=[N1; N2; N3; N4];
          xcoord=dot(Nv,elem_coord(:,1));
          ycoord=dot(Nv,elem_coord(:,2));
          xgp=[xcoord,ycoord]; % integration point [ x(xi,eta) y(xi,eta) ]   
          % FEM displacements at Gauss point                 
          N=[N1,N2,N3,N4];  
          uh_xgp=N*uh_elem_column;           
          % exact solutions at Gauss point  
          [u_exact_xgp,~,gradient_exact_xgp]=...
                        exact_solutions_poisson2d(exact_solution_handle,xgp,matProps);    
          grad_exact_xgp=[gradient_exact_xgp.x;gradient_exact_xgp.y]; % [dudx; dudy] 
          % errors at Gauss point
          u_error_xgp=uh_xgp-u_exact_xgp;
          grad_error_xgp=grad_uh_xgp-grad_exact_xgp;
          % L2-norm of the error
          L2_norm_top=L2_norm_top+dot(u_error_xgp,u_error_xgp)*wxi(gpxi)*weta(gpeta)*detJ;
          L2_norm_bottom=L2_norm_bottom+dot(u_exact_xgp,u_exact_xgp)*wxi(gpxi)*weta(gpeta)*detJ;     
          % H1-seminorm of the error       
          H1_seminorm_top=H1_seminorm_top+...
                         dot(grad_error_xgp,grad_error_xgp)*wxi(gpxi)*weta(gpeta)*detJ;
          H1_seminorm_bottom=H1_seminorm_bottom+...
                         dot(grad_exact_xgp,grad_exact_xgp)*wxi(gpxi)*weta(gpeta)*detJ;            
        end
      end      
    
    end    
    
  else
    throw_error('Error in fem_norms_poisson2d.m: vemlab_method');
  end
  
  h_size=sqrt(h_size/num_elem); % square root of the average area of the elements   
  
  % Finish computation of the norms of the errors
  L2rel = sqrt(L2_norm_top/L2_norm_bottom);
%   L2abs = sqrt(L2_norm_top);  
  H1rel = sqrt(H1_seminorm_top/H1_seminorm_bottom);
%   H1abs = sqrt(H1_seminorm_top);
  
  L2prel = []; % L2 norm of the pressure error (used in linear elasticity)
  
  % Print out norms
  fprintf('Relative L2-norm of the error = %.20f\n',L2rel); 
  fprintf('Relative H1-seminorm of the error = %.20f\n',H1rel);  
%   fprintf('Absolute L2-norm of the error = %.20f\n',L2abs); 
%   fprintf('Absolute H1-seminorm of the error = %.20f\n',H1abs);  
%   fprintf('Mesh spacing (h) = %.12f\n',h_max);    
  fprintf('Element size (sqrt(average area of the elements)): h = %.12f\n\n',h_size);
end