%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                     VEMLab
%           Source code  : http://camlab.cl/research/software/vemlab/
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:                   gauss_points_T3
%
% Created by : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
% Updated by :
% Updated by :
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% Compute gauss points data for 3-node triangles according to the specified 
% quadrature rule. 
%
% Usage
% =====
% [x,w] = gauss_points_T3(order,xyz)
%
% Input
% =====
% order  : quadrature order
% xyz    : element nodal coordinates [x1 y1 z1; x2 y2 z2; ...; x4 y4 z4]
%
% Output
% ======
% x      : Gauss points in cartesian coordinates 
%          [x1 y1 z1; x2 y2 z2; ...; x_ngp y_ngp z_ngp]
% w      : Gauss weights (including area or volumen of the element)
%
%-------------------------------------------------------------------------------
% References 
% ==========
% [1] Zienkiewicz O.C., Taylor R.L., The Finite Element Method. Fifth
% Edition. Volume 1, The Basis. Butterworth-Heinemann, MA, USA.
% [2] Jinyun. Y., Symmetric Gaussian quadrature formulae for tetrahedronal
% regions. Computer Methods in Applied Mechanics and Engineering 1984, 43:349-353
% [3] Keast. P., Moderate degree tetrahedral quadrature formulas. Computer
% Methods in Applied Mechanics and Engineering 1986, 55:339-348
%
%-------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Dec. 26, 2017: first realease (by A. Ortiz-Bernardin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,w] = gauss_points_T3(order,xyz)

  if (order==1) % linear
    tcoord=[1.0/3.0 1.0/3.0 1.0/3.0]; % triangular coordinates
    w(1)=1.0; % weight    
    x=tcoord*xyz; % gauss points in cartesian coordinates
  elseif (order==2) % quadratic
    tcoord=[2.0/3.0 1.0/6.0 1.0/6.0;... 
            1.0/6.0 2.0/3.0 1.0/6.0;...
            1.0/6.0 1.0/6.0 2.0/3.0]; 
    w(1)=1.0/3.0; w(2)=1.0/3.0; w(3)=1.0/3.0; 
    x=tcoord*xyz;
  elseif (order==3) % cubic
    tcoord=[0.659027622374092 0.231933368553031 0.109039009072877;...
            0.109039009072877 0.659027622374092 0.231933368553031;...
            0.231933368553031 0.109039009072877 0.659027622374092;...
            0.109039009072877 0.231933368553031 0.659027622374092;...
            0.231933368553031 0.659027622374092 0.109039009072877;...
            0.659027622374092 0.109039009072877 0.231933368553031];    
    w(1)=1/6; w(2)=1/6; w(3)=1/6; w(4)=1/6; w(5)=1/6; w(6)=1/6;      
    x=tcoord*xyz;  
  elseif (order==4) % quartic
    tcoord=[0.816847572980459 0.091576213509771 0.091576213509771;...
            0.091576213509771 0.816847572980459 0.091576213509771;...
            0.091576213509771 0.091576213509771 0.816847572980459;...
            0.108103018168070 0.445948490915965 0.445948490915965;...
            0.445948490915965 0.108103018168070 0.445948490915965;...
            0.445948490915965 0.445948490915965 0.108103018168070]; 
    w(1)=0.109951743655322; w(2)=0.109951743655322; w(3)=0.109951743655322;
    w(4)=0.223381589678011; w(5)=0.223381589678011; w(6)=0.223381589678011;
    x=tcoord*xyz;  
  elseif (order==5) % quintic
    tcoord=[1.0/3.0 1.0/3.0 1.0/3.0;...
            0.059715871789770 0.470142064105115 0.470142064105115;...
            0.470142064105115 0.059715871789770 0.470142064105115;...
            0.470142064105115 0.470142064105115 0.059715871789770;...
            0.797426985353087 0.101286507323456 0.101286507323456;...
            0.101286507323456 0.797426985353087 0.101286507323456;... 
            0.101286507323456 0.101286507323456 0.797426985353087];
    w(1)=0.225000000000000; w(2)=0.132394152788506; w(3)=0.132394152788506; 
    w(4)=0.132394152788506; w(5)=0.125939180544827; w(6)=0.125939180544827; 
    w(7)=0.125939180544827;
    x=tcoord*xyz;   
  else
    throw_error('In gauss_triangulate.m: order.\n');
  end
  
  % compute area of the triangle and update gauss weights
  area=1.0/2.0*abs(xyz(1,1)*(xyz(3,2)-xyz(2,2))+...
       xyz(3,1)*(xyz(2,2)-xyz(1,2))+...
       xyz(2,1)*(xyz(1,2)-xyz(3,2)));
  w=w*area; 
    
end
