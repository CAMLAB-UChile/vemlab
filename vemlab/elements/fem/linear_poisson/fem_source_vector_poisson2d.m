function fb = fem_source_vector_poisson2d(verts,config,source_term_fun_values,matProps) 
  if strcmp(config.vemlab_method,'FEM2DT3')
    % area of the element 
    area=triangle_area(verts);
    % vector of nodal source terms
    b=source_term_function_poisson2d(verts,source_term_fun_values,matProps);
    
    b1=b(1); b2=b(2); b3=b(3);    % when b is not a constant, b is linearly 
                                  % interpolated inside the element and we end
                                  % up with nodal values of b. If the nodal
                                  % values are all equal, then these "nodal" 
                                  % source terms recover tha case b = const.
                                  % Ref.: J. Fish and T. Belytschko, 
                                  % A First Course in Finite Elements, 
                                  % John Wiley & Sons Ltd, 2007
    fb=(area/12)*[2*b1+b2+b3;...
                  b1+2*b2+b3;...
                  b1+b2+2*b3]; 
    % note that if nodal values are all equal, we get fb = 1/3*b*ones(3,1), 
    % where b is a constant source term, that is, the fb vector for a constant
    % source b is recovered.
  elseif strcmp(config.vemlab_method,'FEM2DQ4')
    %  compute body force vector by numerical integration
    fb=zeros(4,1);    
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
        N1=(1-xi(gpxi))*(1-eta(gpeta))/4;
        N2=(1+xi(gpxi))*(1-eta(gpeta))/4;
        N3=(1+xi(gpxi))*(1+eta(gpeta))/4;
        N4=(1-xi(gpxi))*(1+eta(gpeta))/4;
        N=[N1,N2,N3,N4];    
        x=N1*xcoord(1)+N2*xcoord(2)+N3*xcoord(3)+N4*xcoord(4);
        y=N1*ycoord(1)+N2*ycoord(2)+N3*ycoord(3)+N4*ycoord(4);
        b=source_term_function_poisson2d([x,y],source_term_fun_values,matProps);
        fb=fb+N'*b*wxi(gpxi)*weta(gpeta)*detJ;
      end
    end    
  else
    throw_error('In fem_body_force_linelast2d.m: vemlab_method');
  end
end