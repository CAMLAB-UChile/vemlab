function fb = fem_body_force_linelast2d(verts,config,body_force_fun_value) 
  if strcmp(config.vemlab_method,'FEM2DT3')
    % area of the element 
    area=triangle_area(verts);
    % vector of nodal body forces
    bf=body_force_function_linelast2d(verts,body_force_fun_value);
    bx1=bf(1); bx2=bf(3); bx3=bf(5); % when bf is not constant, these are
    by1=bf(2); by2=bf(4); by3=bf(6); % the nodal values of the components of bf.
                                     % Ref.: J. Fish and T. Belytschko, 
                                     % A First Course in Finite Elements, 
                                     % John Wiley & Sons Ltd, 2007
    fb=(area/12)*[2*bx1+bx2+bx3;...
                  2*by1+by2+by3;...
                  bx1+2*bx2+bx3;...
                  by1+2*by2+by3;...
                  bx1+bx2+2*bx3;...
                  by1+by2+2*by3];
  elseif strcmp(config.vemlab_method,'FEM2DQ4')
    %  compute body force vector by numerical integration
    fb=zeros(8,1);    
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
        N1=(1-xi(gpxi))*(1-eta(gpeta))/4;
        N2=(1+xi(gpxi))*(1-eta(gpeta))/4;
        N3=(1+xi(gpxi))*(1+eta(gpeta))/4;
        N4=(1-xi(gpxi))*(1+eta(gpeta))/4;
        N=[N1,0,N2,0,N3,0,N4,0;...
           0,N1,0,N2,0,N3,0,N4];    
        x=N1*xcoord(1)+N2*xcoord(2)+N3*xcoord(3)+N4*xcoord(4);
        y=N1*ycoord(1)+N2*ycoord(2)+N3*ycoord(3)+N4*ycoord(4);
        bf=body_force_function_linelast2d([x,y],body_force_fun_value);
        fb=fb+N'*bf*wxi(gpxi)*weta(gpeta)*detJ;
      end
    end    
  else
    throw_error('In fem_body_force_linelast2d.m: vemlab_method');
  end
end