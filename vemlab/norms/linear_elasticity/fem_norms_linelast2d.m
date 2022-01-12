function [h_size,L2rel,H1rel,L2prel]=fem_norms_linelast2d(exact_sol,uh_global,domainMesh,config,matProps)
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
      dofs=zeros(2*length(nodes),1);
      for node_i=1:length(nodes)
        dofs(2*node_i-1)=2*nodes(node_i)-1;
        dofs(2*node_i)=2*nodes(node_i);
      end

      % element VEM nodal solution
      uh_elem_column=uh_global(dofs);
      
      % area of the element
      area=triangle_area(elem_coord);        

      % h = nodal spacing = maximum distance between any two neighbor nodes in the mesh
%       h_size=max_edge_size(elem_coord);
%       if h_size>h_max
%         h_max=h_size;
%       end  
      h_E=area;    
      h_size=h_size+h_E;

      % deformation matrix
      x=elem_coord(:,1); y=elem_coord(:,2);  
      % In FEM, B is defined for [e11,e22,2*e12]
      B=(1/(2*area))*[y(2)-y(3),0,y(3)-y(1),0,y(1)-y(2),0;...
                      0,x(3)-x(2),0,x(1)-x(3),0,x(2)-x(1);...
                      x(3)-x(2),y(2)-y(3),x(1)-x(3),y(3)-y(1),x(2)-x(1),y(1)-y(2)];    

      % FEM strains
      strain_uh=B*uh_elem_column; % In FEM, B is defined for [e11,e22,2*e12]
      strain_uh(3)=strain_uh(3)/2; % To be consistent with D defined as in VEM,
                                   % divide strain coming from B by 2

      % Gauss points generation by triangulation of the polygonal element
      gauss_order=2;
      [gauss_points,gauss_weights]=gauss_points_T3(gauss_order,elem_coord); 

      % loop over element Gauss points
      for ngp=1:length(gauss_weights)
        wgp=gauss_weights(ngp);
        xgp=gauss_points(ngp,:); % integration point in the form [x y z]            

        % exact solution at Gauss point
        [u_exact_xgp,p_exact_xgp,strainvec_exact_xgp]=...
                      exact_solutions_linelast2d(exact_sol,xgp,matProps);
        strain_exact_xgp=strainvec_exact_xgp';                           
                        

        % displacement projection at Gauss point
        if strcmp(config.vemlab_method,'FEM2DT3')
          N1=(1/(2*area))*((x(2)*y(3)-x(3)*y(2))+(y(2)-y(3))*xgp(1)+(x(3)-x(2))*xgp(2));
          N2=(1/(2*area))*((x(3)*y(1)-x(1)*y(3))+(y(3)-y(1))*xgp(1)+(x(1)-x(3))*xgp(2));
          N3=(1/(2*area))*((x(1)*y(2)-x(2)*y(1))+(y(1)-y(2))*xgp(1)+(x(2)-x(1))*xgp(2));
          N=[N1,0,N2,0,N3,0;...
             0,N1,0,N2,0,N3];
          uh_xgp=N*uh_elem_column;
        else
          throw_error('Error in fem_norms_linelast2d.m: vemlab_method');
        end

        % errors at Gauss point
        u_error_xgp=uh_xgp-u_exact_xgp;
        strain_error_xgp=strain_uh-strain_exact_xgp;

        % L2-norm of the error
        L2_norm_top=L2_norm_top+dot(u_error_xgp,u_error_xgp)*wgp;
        L2_norm_bottom=L2_norm_bottom+dot(u_exact_xgp,u_exact_xgp)*wgp;     

        % H1-seminorm of the error       
        H1_seminorm_top=H1_seminorm_top+...
                       dot(strain_error_xgp,(matProps.D)*strain_error_xgp)*wgp;
        H1_seminorm_bottom=H1_seminorm_bottom+...
                       dot(strain_exact_xgp,(matProps.D)*strain_exact_xgp)*wgp;       
      end 
      
    end
    
    h_size=sqrt(h_size/num_elem); % square root of the average area of the elements 
    
  elseif strcmp(config.vemlab_method,'FEM2DQ4')
    
    % loop over elements
    for i = 1:num_elem
      
      % element dofs
      nodes=connect{i}; 
      elem_coord=coords(nodes,:);    
      dofs=zeros(2*length(nodes),1);
      for node_i=1:length(nodes)
        dofs(2*node_i-1)=2*nodes(node_i)-1;
        dofs(2*node_i)=2*nodes(node_i);
      end

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
      num_gp=2; % number of Gauss points per axis
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
          B=[dN1dx,0,dN2dx,0,dN3dx,0,dN4dx,0;...
             0,dN1dy,0,dN2dy,0,dN3dy,0,dN4dy;...
             dN1dy,dN1dx,dN2dy,dN2dx,dN3dy,dN3dx,dN4dy,dN4dx];
          % FEM strains
          strain_uh_xgp=B*uh_elem_column; % In FEM, B is defined for [e11,e22,2*e12]
          strain_uh_xgp(3)=strain_uh_xgp(3)/2; % To be consistent with D defined as in VEM,
                                               % divide strain coming from B by 2   
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
          N=[N1,0,N2,0,N3,0,N4,0;...
             0,N1,0,N2,0,N3,0,N4];  
          uh_xgp=N*uh_elem_column;           
          % exact solutions at Gauss point                        
          [u_exact_xgp,p_exact_xgp,strainvec_exact_xgp]=...
                        exact_solutions_linelast2d(exact_sol,xgp,matProps);
          strain_exact_xgp=strainvec_exact_xgp';                             
                                          
          % errors at Gauss point
          u_error_xgp=uh_xgp-u_exact_xgp;
          strain_error_xgp=strain_uh_xgp-strain_exact_xgp;
          % L2-norm of the error
          L2_norm_top=L2_norm_top+dot(u_error_xgp,u_error_xgp)*wxi(gpxi)*weta(gpeta)*detJ;
          L2_norm_bottom=L2_norm_bottom+dot(u_exact_xgp,u_exact_xgp)*wxi(gpxi)*weta(gpeta)*detJ;     
          % H1-seminorm of the error       
          H1_seminorm_top=H1_seminorm_top+...
                         dot(strain_error_xgp,(matProps.D)*strain_error_xgp)*wxi(gpxi)*weta(gpeta)*detJ;
          H1_seminorm_bottom=H1_seminorm_bottom+...
                         dot(strain_exact_xgp,(matProps.D)*strain_exact_xgp)*wxi(gpxi)*weta(gpeta)*detJ;            
        end
      end      
    
    end 
    
    h_size=sqrt(h_size/num_elem); % square root of the average area of the elements    
    
  else
    throw_error('Error in fem_norms_linelast2d.m: vemlab_method');
  end
  
  % Finish computation of the norms of the errors
  L2rel = sqrt(L2_norm_top/L2_norm_bottom);
%   L2abs = sqrt(L2_norm_top);  
  H1rel = sqrt(H1_seminorm_top/H1_seminorm_bottom);
%   H1abs = sqrt(H1_seminorm_top);
  
  L2prel = []; % L2 norm of the pressure error
  
  % Print out norms
  fprintf('Relative L2-norm of the error = %.20f\n',L2rel); 
  fprintf('Relative H1-seminorm of the error = %.20f\n',H1rel);  
%   fprintf('Absolute L2-norm of the error = %.20f\n',L2abs); 
%   fprintf('Absolute H1-seminorm of the error = %.20f\n',H1abs);  
%   fprintf('Mesh spacing (h) = %.12f\n',h_max);    
  fprintf('Element size (sqrt(average area of the elements)): h = %.12f\n\n',h_size);  
  
end