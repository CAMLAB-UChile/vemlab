%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                     VEMLab
%           Source code  : http://camlab.cl/research/software/vemlab/
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:                   plot_solution_linelast 
%
% Created by : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
% Updated by :
% Updated by :
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% Plot displacement solutions on a polygonal mesh. Polygon can be a triangle,
% a square, a pentagon, etc.
%
% Usage
% =====
% plot_solution_linelast2d(mesh,solution,titleSolutionX,...
%                          titleSolutionY,titleSolutionNorm,...
%                          plotMeshOverResults,solutionType)
%
% Input
% =====
% mesh     : structure containing mesh data (coords,connect,etc.)
% solution : nodal displacement solution
% titleSolutionX : title for the plot of solution in the X direction (latex format)
% titleSolutionY : title for the plot of solution in the Y direction (latex format)
% titleSolutionNorm : title for the plot of the norm of solutions (latex format)
% plotMeshOverResults : overlay mesh over result plot? 'yes' or 'no'
% solutionType : 'VEM', 'FEM', 'exact' o whatever is appropriate
%
% Output
% ======
%
%-------------------------------------------------------------------------------
% References 
% ==========
% Part of this code has been taken from the plot_solution.m function provided in
% the source code of:
%
% [1] O. J. Sutton. The virtual element method in 50 lines of MATLAB.
%     Numerical Algorithms 2017; 75(4):1141–1159
%
%-------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Dec. 26, 2017: first realease (by A. Ortiz-Bernardin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_solution_linelast2d(mesh,solution,titleSolutionX,...
                                  titleSolutionY,titleSolutionNorm,...
                                  plotMeshOverResults,solutionType)
 
  fprintf('Plotting %s solution...\n',solutionType);                               
  nodes=mesh.coords;
  nodesNumber=size(nodes,1);  
  polygons=mesh.connect;
  
  range=1:nodesNumber;
  displacementsX=solution(2*range-1);
  displacementsY=solution(2*range);  
  displacementsNorm=sqrt(displacementsX.*displacementsX+displacementsY.*displacementsY);

  figure; title(titleSolutionNorm,'Interpreter','latex','FontSize',18);
  maxNumVertices = max(cellfun(@numel,polygons));
  padFunc = @(vertList) [vertList' NaN(1,maxNumVertices-numel(vertList))];
  elements = cellfun(padFunc,polygons,'UniformOutput',false);
  elements = vertcat(elements{:});
  data = [nodes,displacementsNorm];
  if strcmp(plotMeshOverResults,'yes')
    patch('Faces',elements,'Vertices',data,...
          'FaceColor','interp','CData',displacementsNorm);  
  else
    patch('Faces',elements,'Vertices',data,...
          'FaceColor','interp','EdgeColor','interp','CData',displacementsNorm);
  end  
%   axis('square') 
  axis equal;  
  dx = 0; dy = 0; dz = 0;
  xlim([min(nodes(:, 1)) - dx, max(nodes(:, 1)) + dx])
  ylim([min(nodes(:, 2)) - dy, max(nodes(:, 2)) + dy])  
  zlim([min(displacementsNorm) - dz, max(displacementsNorm) + dz])
  xlabel('$x$','Interpreter','latex','FontSize',18); 
  ylabel('$y$','Interpreter','latex','FontSize',18); 
  zlabel(titleSolutionNorm,'Interpreter','latex','FontSize',18);
  colorbar
  colormap jet
  %set(gcf,'Renderer','painters')    
  set(gcf,'InvertHardcopy','off','Color',[1 1 1])
  set(gca, 'FontSize', 12);
  
  figure; title(titleSolutionX,'Interpreter','latex','FontSize',18);
  maxNumVertices = max(cellfun(@numel,polygons));
  padFunc = @(vertList) [vertList' NaN(1,maxNumVertices-numel(vertList))];
  elements = cellfun(padFunc,polygons,'UniformOutput',false);
  elements = vertcat(elements{:});
  data = [nodes,displacementsX];
  if strcmp(plotMeshOverResults,'yes')
    patch('Faces',elements,'Vertices',data,...
          'FaceColor','interp','CData',displacementsX);  
  else
    patch('Faces',elements,'Vertices',data,...
          'FaceColor','interp','EdgeColor','interp','CData',displacementsX);
  end      
%   axis('square') 
  axis equal;  
  dx = 0; dy = 0; dz = 0;
  xlim([min(nodes(:,1)) - dx, max(nodes(:,1)) + dx])
  ylim([min(nodes(:,2)) - dy, max(nodes(:,2)) + dy])
  zlim([min(displacementsX) - dz, max(displacementsX) + dz])
  xlabel('$x$','Interpreter','latex','FontSize',18); 
  ylabel('$y$','Interpreter','latex','FontSize',18); 
  zlabel(titleSolutionX,'Interpreter','latex','FontSize',18);
  colorbar
  colormap jet
  %set(gcf,'Renderer','painters') 
  set(gcf,'InvertHardcopy','off','Color',[1 1 1])  
  set(gca, 'FontSize', 12);
  
  figure; title(titleSolutionY,'Interpreter','latex','FontSize',18);
  maxNumVertices = max(cellfun(@numel,polygons));
  padFunc = @(vertList) [vertList' NaN(1,maxNumVertices-numel(vertList))];
  elements = cellfun(padFunc,polygons,'UniformOutput',false);
  elements = vertcat(elements{:});
  data = [nodes,displacementsY];
  if strcmp(plotMeshOverResults,'yes')
    patch('Faces',elements,'Vertices',data,...
          'FaceColor','interp','CData',displacementsY);  
  else
    patch('Faces',elements,'Vertices',data,...
          'FaceColor','interp','EdgeColor','interp','CData',displacementsY);
  end   
%   axis('square') 
  axis equal;  
  dx = 0; dy = 0; dz = 0;
  xlim([min(nodes(:,1)) - dx, max(nodes(:,1)) + dx])
  ylim([min(nodes(:,2)) - dy, max(nodes(:,2)) + dy])
  zlim([min(displacementsY) - dz, max(displacementsY) + dz])
  xlabel('$x$','Interpreter','latex','FontSize',18); 
  ylabel('$y$','Interpreter','latex','FontSize',18); 
  zlabel(titleSolutionY,'Interpreter','latex','FontSize',18);
  colorbar
  colormap jet
  %set(gcf,'Renderer','painters')  
  set(gcf,'InvertHardcopy','off','Color',[1 1 1])  
  set(gca, 'FontSize', 12);  
  
end