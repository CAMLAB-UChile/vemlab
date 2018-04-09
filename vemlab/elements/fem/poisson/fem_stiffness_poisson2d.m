function Kfem = fem_stiffness_poisson2d(verts,matProps,config)
  if strcmp(config.vemlab_method,'FEM2DT3')
    % area of the element
    area=triangle_area(verts);
    % deformation matrix
    x=verts(:,1); y=verts(:,2);  
    B=(1/(2*area))*[y(2)-y(3),y(3)-y(1),y(1)-y(2);...
                    x(3)-x(2),x(1)-x(3),x(2)-x(1)];    
    % FEM element stiffness
    k=matProps.k;    % isotropic material
    Kfem=k*area*(B'*B);
  elseif strcmp(config.vemlab_method,'FEM2DQ4')
    % material matrix
    k=matProps.k;    % isotropic material   
    % compute stiffness matrix by numerical integration
    Kfem=zeros(4,4);
    int_order=2;
    xi=gauss_points_1d(int_order);
    eta=gauss_points_1d(int_order);
    wxi=gauss_weights_1d(int_order);
    weta=gauss_weights_1d(int_order);
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
        xcoord=verts(:,1); ycoord=verts(:,2);
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
        B=[dN1dx,dN2dx,dN3dx,dN4dx;...
           dN1dy,dN2dy,dN3dy,dN4dy];
        Kfem=Kfem+k*(B'*B)*wxi(gpxi)*weta(gpeta)*detJ;
      end
    end
  else
    throw_error('Error in fem_stiffness_linelast2d.m: vemlab_method');
  end
end