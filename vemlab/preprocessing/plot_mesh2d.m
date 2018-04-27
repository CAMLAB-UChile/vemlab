%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                     VEMLab
%           Source code  : http://camlab.cl/research/software/vemlab/
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:                       plot_mesh 
%
% Created by : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
% Updated by :
% Updated by :
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% Plot a polygonal mesh. Polygon can be a triangle, square, pentagon, etc.
%
% Usage
% =====
% plot_mesh2d(domainMesh,config)
%
% Input
% =====
% domainMesh : structure containing mesh data (coords,connect,etc.)
% config : structure storing VEMLab configuration options and behavior
%
% Output
% ======
%
%-------------------------------------------------------------------------------
% References 
% ==========
%
%-------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Dec. 26, 2017: first realease (by A. Ortiz-Bernardin)
% Mar. 25, 2018: add config (by A. Ortiz-Bernardin)
% Apr. 19, 2018: improve the plotting of axis and fonts
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_mesh2d(domainMesh,config)  
  if strcmp(config.plot_mesh,'yes')
    fprintf('Plotting mesh...\n');
    points=domainMesh.coords; 
    polygons=domainMesh.connect;  
  %   set(gcf,'Renderer','painters')
    figure; 
    maxNumVertices = max(cellfun(@numel,polygons));
    padFunc = @(vertList) [vertList' NaN(1,maxNumVertices-numel(vertList))];
    elements = cellfun(padFunc,polygons,'UniformOutput',false);
    elements = vertcat(elements{:});
    patch('Faces',elements,'Vertices',points,'FaceColor','w'); 
  %   axis('square')  
    set(gca,'FontName','Times New Roman','FontSize',14,'FontWeight','bold');
    axis equal; %axis off; 
    dx = 0; dy = 0;
    xlim([min(points(:, 1)) - dx, max(points(:, 1)) + dx]);
    ylim([min(points(:, 2)) - dy, max(points(:, 2)) + dy]);
    xlabel('$x$','FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    ylabel('$y$','FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
  end
end