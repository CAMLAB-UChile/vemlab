%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                      VEMLab
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
% Jan. 29, 2020: add control commands plot_mesh_linewidth, plot_mesh_nodes, 
%                plot_mesh_nodesize, plot_mesh_axis
% Apr. 19, 2018: improve the plotting of axis and fonts
% Mar. 25, 2018: add config (by A. Ortiz-Bernardin)
% Dec. 26, 2017: first realease (by A. Ortiz-Bernardin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_mesh2d(domainMesh,config)   
  if strcmp(config.plot_mesh,'yes')
    fprintf('Plotting mesh...\n');
    points=domainMesh.coords;
    polygons=domainMesh.connect;  
    
%     vemlab_logo_file=[config.icons_folder_location,'vemlab.png']; 
%     Img = imread(vemlab_logo_file); 
    
  %   set(gcf,'Renderer','painters')
    figure; 
    maxNumVertices = max(cellfun(@numel,polygons));
    padFunc = @(vertList) [vertList' NaN(1,maxNumVertices-numel(vertList))];
    elements = cellfun(padFunc,polygons,'UniformOutput',false);
    elements = vertcat(elements{:});
    patch('Faces',elements,'Vertices',points,'FaceColor','w','LineWidth',config.plot_mesh_linewidth); 
  %   axis('square')  
%     set(gca,'FontName','Times New Roman','FontSize',14,'FontWeight','bold');
%     set(gca,'FontName','Times New Roman','FontSize',13,'FontWeight','bold');    
    set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','bold');  
    axis equal; 
    if strcmp(config.plot_mesh_axis,'no')
      axis('off')
%       if strcmp(config.plot_vemlab_logo,'yes') 
%         Lx=max(points(:,1))-min(points(:,1));
%         Ly=max(points(:,2))-min(points(:,2));
%         fLx=0.391304348*Lx;
%         fLy=0.045000000*Ly;
%         image([max(points(:,1))*0.995-fLx max(points(:,1))*0.995],[min(points(:,2))-0.8*fLy min(points(:,2))-1.8*fLy],Img); 
%         hold on    
%       end        
    else
      xlabel('$x_1$','FontWeight','normal','FontSize',18,'FontName',...
             'Times New Roman','Interpreter','latex');
      ylabel('$x_2$','FontWeight','normal','FontSize',18,'FontName',...
             'Times New Roman','Interpreter','latex');      
      axis normal
      ax = gca;
      if ~isempty(config.axis_xtick)
        ax.XTick=config.axis_xtick;
      end
      if ~isempty(config.axis_ytick)    
        ax.YTick=config.axis_ytick;
      end      
%       if strcmp(config.plot_vemlab_logo,'yes') 
%         hold on
%         Lx=max(points(:,1))-min(points(:,1));
%         Ly=max(points(:,2))-min(points(:,2));
%         fLx=0.391304348*Lx;
%         fLy=0.045000000*Ly;
%         image([max(points(:,1))*0.995-fLx max(points(:,1))*0.995],[min(points(:,2))-0.8*fLy min(points(:,2))-1.8*fLy],Img);    
%         axis([min(points(:,1)) max(points(:,1)) min(points(:,2))-2.65*fLy max(points(:,2))])
%         pbaspect([abs(min(points(:,1))-max(points(:,1))) abs(min(points(:,2))-max(points(:,2))) 1])   
%         hold on 
%       else
        axis([min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))])
        pbaspect([abs(min(points(:,1))-max(points(:,1))) abs(min(points(:,2))-max(points(:,2))) 1])
        hold on     
%       end     
    end
      
    if strcmp(config.plot_mesh_nodes,'yes')
      hold on
      x=points(:,1);
      y=points(:,2);
      plot(x,y,'bo','MarkerFaceColor','b','MarkerSize',config.plot_mesh_nodesize);
      hold off
    end
    
    if strcmp(config.print_figures,'yes')
      set(gcf,'InvertHardcopy','off','PaperType','<custom>','PaperSize',[5.83 4.38],'Color',[1 1 1],'Renderer','painters');
      output_figure=[config.matlab_figures_output_folder_location,config.mesh_filename,'_mesh','.pdf']; 
      print('-r300',output_figure,'-dpdf','-bestfit');    
    end
    if strcmp(config.save_matlab_figures,'yes')  
      set(gcf,'InvertHardcopy','off','PaperType','<custom>','PaperSize',[5.83 4.38],'Color',[1 1 1],'Renderer','painters');
      output_figure=[config.matlab_figures_output_folder_location,config.mesh_filename,'_mesh','.fig'];
      saveas(gcf,output_figure);
    end    
   
  end
end