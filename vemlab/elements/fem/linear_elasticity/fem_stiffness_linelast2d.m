function Kfem = fem_stiffness_linelast2d(verts,matProps,config)
  if strcmp(config.vemlab_method,'FEM2DT3')
    % area of the element
    area=triangle_area(verts);
    % deformation matrix
    x=verts(:,1); y=verts(:,2);  
    % In FEM, B is defined for [e11,e22,2*e12]
    B=(1/(2*area))*[y(2)-y(3),0,y(3)-y(1),0,y(1)-y(2),0;...
                    0,x(3)-x(2),0,x(1)-x(3),0,x(2)-x(1);...
                    x(3)-x(2),y(2)-y(3),x(1)-x(3),y(3)-y(1),x(2)-x(1),y(1)-y(2)];    
    % FEM element stiffness
    D=matProps.D;    % This D is defined as per VEM for strain = [e11 e22 e12].
    D(3,3)=D(3,3)/4; % To avoid troubles, come back to the standard FEM definition of D 
                     % to obtain the stress by multiplying with [e11 e22 2*e12]
    Kfem=area*B'*D*B;
  elseif strcmp(config.vemlab_method,'FEM2DQ4')
    % material matrix
    D=matProps.D;    % This D is defined as per VEM for strain = [e11 e22 e12].
    D(3,3)=D(3,3)/4; % To avoid troubles, come back to the standard FEM definition of D 
                     % to obtain the stress by multiplying with [e11 e22 2*e12]   
    % compute stiffness matrix by numerical integration
    Kfem=zeros(8,8);
    num_gp=config.number_of_gauss_points_per_axis_FEM2DQ4;
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
        B=[dN1dx,0,dN2dx,0,dN3dx,0,dN4dx,0;...
           0,dN1dy,0,dN2dy,0,dN3dy,0,dN4dy;...
           dN1dy,dN1dx,dN2dy,dN2dx,dN3dy,dN3dx,dN4dy,dN4dx];
        Kfem=Kfem+B'*D*B*wxi(gpxi)*weta(gpeta)*detJ;
      end
    end
  else
    throw_error('Error in fem_stiffness_linelast2d.m: vemlab_method');
  end
end