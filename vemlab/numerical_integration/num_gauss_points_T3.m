function ngauss = num_gauss_points_T3(order)
% the number of Gauss points are according to the quadratures in file:
% gauss_points_T3.m
  switch order
    case 1
      ngauss=1;
    case 2
      ngauss=3;
    case 3
      ngauss=6;
    case 4
      ngauss=6;     
    case 5
      ngauss=7;              
    otherwise
      throw_error('In num_gauss_points_T3.m: order');
  end  
end
