%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                    VEMLab
%           Source code  : http://camlab.cl/research/software/vemlab/
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:                       plot_mesh2d 
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
% Jan. 29, 2020: add control commands plot_mesh_linewidth, plot_mesh_nodes, 
%                plot_mesh_nodesize, plot_mesh_axis
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
    patch('Faces',elements,'Vertices',points,'FaceColor','w','LineWidth',config.plot_mesh_linewidth); 
  %   axis('square')  
    set(gca,'FontName','Times New Roman','FontSize',14,'FontWeight','bold');
    axis equal; 
    if strcmp(config.plot_mesh_axis,'no')
      axis('off')
    end    
    %axis off; 
    dx = 0; dy = 0;
    xlim([min(points(:, 1)) - dx, max(points(:, 1)) + dx]);
    ylim([min(points(:, 2)) - dy, max(points(:, 2)) + dy]);
    xlabel('$x$','FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    ylabel('$y$','FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    if strcmp(config.plot_mesh_nodes,'yes')
      hold on
      x=points(:,1);
      y=points(:,2);
      plot(x,y,'bo','MarkerFaceColor','b','MarkerSize',config.plot_mesh_nodesize);
      hold off
    end
  end
end