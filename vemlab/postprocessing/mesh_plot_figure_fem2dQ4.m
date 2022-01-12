function mesh_plot_figure_fem2dQ4(config,polygons,nodescoords,gp_list,xq,yq,fieldvar,...
                                  fieldvartitleclb,fieldvartitlefigfile,clb_limits,dxtitle,soltype)

  % fieldvar: field variable to plot (1d matrix)
  % polygons: elements connectivity (cell structure)
  % nodescoords: nodal coordinates of the mesh (2d matrix)
  % fieldvartitleclb: colorbar title (string)
  % fieldvartitlefigfile: title for the pdf and matlab output figure (string)  
  % clb_limits: colorbar limits [min max]   
  % dxtitle: factor to move fieldvartitle horizontally
  % soltype: 'numerical' or 'exact' 
  
  vq=griddata(gp_list.x(:),gp_list.y(:),fieldvar(:),xq,yq);    
  
  % check that the query points (xq,yq) are inside the domain
  num_elem=length(polygons);
  indomain(1:length(xq(:)))=false;
  for e=1:num_elem
    vcoords=nodescoords(polygons{e},:);
    in=inpolygon(xq(:),yq(:),vcoords(:,1),vcoords(:,2));
    indomain(in)=true; % indices of the entries in xq(:) and yq(:) that are inside the polygon      
  end
  vq(~indomain)=NaN; % fill with NaN entries that are outside the domain

  % plot figure
  figure; 
  mesh(xq,yq,vq,'FaceColor','interp','EdgeColor','interp');
  hold on
  if strcmp(config.plot_mesh_over_results,'yes')
%     maxNumVertices = max(cellfun(@numel,polygons));
%     padFunc = @(vertList) [vertList' NaN(1,maxNumVertices-numel(vertList))];
%     elements = cellfun(padFunc,polygons,'UniformOutput',false);
%     elements = vertcat(elements{:});
%     patch('Faces',elements,'Vertices',nodescoords,...
%           'FaceColor','none','LineWidth',config.plot_mesh_linewidth);  
%     hold off
      throw_warning('For FEM2DQ4, plot of mesh over results not available for secondary variables. No mesh will be plotted.');        
  end

  if strcmp(config.matlab_colormap,'parula')
    if isempty(config.colormap_number_of_colors)
      colormap parula
    else
      colormap(parula(config.colormap_number_of_colors));
    end
  elseif strcmp(config.matlab_colormap,'jet')
    if isempty(config.colormap_number_of_colors)
      colormap jet
    else
      colormap(jet(config.colormap_number_of_colors));
    end    
  elseif strcmp(config.matlab_colormap,'hsv')
    if isempty(config.colormap_number_of_colors)
      colormap hsv
    else
      colormap(hsv(config.colormap_number_of_colors));
    end 
  elseif strcmp(config.matlab_colormap,'hot')
    if isempty(config.colormap_number_of_colors)
      colormap hot
    else
      colormap(hot(config.colormap_number_of_colors));
    end 
  elseif strcmp(config.matlab_colormap,'cool')
    if isempty(config.colormap_number_of_colors)
      colormap cool
    else
      colormap(cool(config.colormap_number_of_colors));
    end 
  elseif strcmp(config.matlab_colormap,'spring')
    if isempty(config.colormap_number_of_colors)
      colormap spring
    else
      colormap(spring(config.colormap_number_of_colors));
    end 
  elseif strcmp(config.matlab_colormap,'summer')
    if isempty(config.colormap_number_of_colors)
      colormap summer
    else
      colormap(summer(config.colormap_number_of_colors));
    end 
  elseif strcmp(config.matlab_colormap,'autumn')
    if isempty(config.colormap_number_of_colors)
      colormap autumn
    else
      colormap(autumn(config.colormap_number_of_colors));
    end     
  elseif strcmp(config.matlab_colormap,'winter')
    if isempty(config.colormap_number_of_colors)
      colormap winter
    else
      colormap(winter(config.colormap_number_of_colors));
    end 
  elseif strcmp(config.matlab_colormap,'gray')
    if isempty(config.colormap_number_of_colors)
      colormap gray
    else
      colormap(gray(config.colormap_number_of_colors));
    end 
  elseif strcmp(config.matlab_colormap,'bone')
    if isempty(config.colormap_number_of_colors)
      colormap bone
    else
      colormap(bone(config.colormap_number_of_colors));
    end 
  elseif strcmp(config.matlab_colormap,'copper')
    if isempty(config.colormap_number_of_colors)
      colormap copper
    else
      colormap(copper(config.colormap_number_of_colors));
    end 
  elseif strcmp(config.matlab_colormap,'pink')
    if isempty(config.colormap_number_of_colors)
      colormap pink
    else
      colormap(pink(config.colormap_number_of_colors));
    end 
  elseif strcmp(config.matlab_colormap,'spectral')    
    if isempty(config.colormap_number_of_colors)
      cmap = cbrewer('div','Spectral',256);
      %colormap(cmap);   
      colormap(flipud(cmap));
    else
      cmap = cbrewer('div','Spectral',config.colormap_number_of_colors);
      %colormap(cmap);   
      colormap(flipud(cmap));      
    end
  elseif strcmp(config.matlab_colormap,'RdYlBu')  
    if isempty(config.colormap_number_of_colors)
      cmap = cbrewer('div','RdYlBu',256);
      %colormap(cmap);   
      colormap(flipud(cmap));
    else
      cmap = cbrewer('div','RdYlBu',config.colormap_number_of_colors);
      %colormap(cmap);   
      colormap(flipud(cmap));      
    end
  else
    if isempty(config.colormap_number_of_colors)
      colormap default
    else
      colormap(parula(config.colormap_number_of_colors));
    end     
  end  

  if strcmp(config.front_matlab_plot_view_orientation,'yes')    
    clb=colorbar('FontName','Segoe UI Semibold','FontSize',11,...
                 'FontWeight','normal','LineWidth',config.colorbar_LineWidth,...
                 'TickLength',config.colorbar_TickLength);
%     clb=colorbar('FontName','Segoe UI','FontSize',10,'FontWeight','bold','LineWidth',1.05,'TickLength',0.02);    
    ax = gca;
    if ~isempty(config.axis_xtick)
      ax.XTick=config.axis_xtick;
    end
    if ~isempty(config.axis_ytick)    
      ax.YTick=config.axis_ytick;
    end    
    ax.XMinorTick = config.figure_XMinorTick;
    ax.YMinorTick = config.figure_YMinorTick;   
    ax.XGrid = config.figure_XGrid;
    ax.YGrid = config.figure_YGrid;    
    ax.XMinorGrid = config.figure_XMinorGrid;
    ax.YMinorGrid = config.figure_YMinorGrid;
    if strcmp(config.figure_GridColor,'default')
      ax.GridColor=[0.15 0.15 0.15];
      ax.MinorGridColor=[0.1 0.1 0.1];      
    else
      ax.GridColor = config.figure_GridColor;
      ax.MinorGridColor = config.figure_MinorGridColor;      
    end    
    if strcmp(config.plot_figure_axis,'yes')
      ax.XColor=ax.XColor;
      ax.YColor=ax.YColor;      
    else
      ax.XColor='none';
      ax.YColor='none';      
    end
    view(0,90);    
%         clb.Position(1)=0.93*clb.Position(1); % horizontal position w/r to the figures's bottom left corner
%         clb.Position(2)=1.78*clb.Position(2); % vertical position w/r to the figure's bottom left corner   
%         clb.Position(2)=1.95*clb.Position(2); % vertical position w/r to the figure's bottom left corner  
%         clb.Position(3)=1.1*clb.Position(3);  % width of the colorbar
%         clb.Position(4)=0.82*clb.Position(4);  % height of the colorbar
%         clb.Position(4)=0.80*clb.Position(4);  % height of the colorbar 
%         clb.Position(1:4)=[clb.Position(1) 1.90*clb.Position(2) clb.Position(3) 0.8*clb.Position(4)];
%     clb.Position(1:4)=[clb.Position(1) clb.Position(2) clb.Position(3) clb.Position(4)];
    clb.Title.String = fieldvartitleclb;
    clb.Title.Interpreter = 'latex';
    clb.Title.FontWeight='bold';
    clb.Title.FontName='Times New Roman';        
    clb.Title.FontSize = 18;
    clb.Title.Units = 'normalized';
%     clb.Title.Position(1:2) = [dxtitle max(clb.Limits)+0.020*(max(clb.Limits)-min(clb.Limits))];   
%     clb.Title.Position(1:2) = [2.0*clb.Title.Position(1) 1.0*clb.Title.Position(2)];
%     clb.Title.Position(1:2) = [dxtitle*clb.Title.Position(1) 1.0*clb.Title.Position(2)];
    clb.Title.Position(1:2) = [dxtitle*clb.Position(1) 1.004*clb.Title.Position(2)]; 
    
    % colorbar limits
    if length(clb_limits)==2
      clb.Limits=clb_limits;   
      if ~isempty(config.colormap_number_of_colors) % we set the ytick variable when using a given number of colors so that the ticks match the color boxes
        nc=config.colormap_number_of_colors;
        set(clb,'ylim',clb_limits,'ytick',[clb_limits(1):(clb_limits(2)-clb_limits(1))/nc:clb_limits(2)]');
      end
    elseif (length(clb_limits)==1)||(length(clb_limits)>2)
      throw_warning('In mesh_plot_figure_fem2dQ4.m: clb_limits must be a two-column row. Using default colorbar limits'); 
    else % default limits are used, but we set the ytick variable when using a given number of colors so that the ticks match the color boxes
      if ~isempty(config.colormap_number_of_colors)
        nc=config.colormap_number_of_colors;
        set(clb,'ylim',clb.Limits,'ytick',[clb.Limits(1):(clb.Limits(2)-clb.Limits(1))/nc:clb.Limits(2)]');
      end      
    end    
    
    if strcmp(config.plot_vemlab_logo,'yes') 
%       ylabel(clb,['\color[rgb]{0.360784 0.4 0.435294}VEM','\color[rgb]{0.850980 0.270588 0.325490}LAB'],...
%              'Rotation',0,'Position',[2.019 min(clb.Limits)-0.04*(max(clb.Limits)-min(clb.Limits))],'FontName','Good Times Rg','FontSize',12)
      ylb=ylabel(clb,['\color[rgb]{0.360784 0.4 0.435294}VEM','\color[rgb]{0.850980 0.270588 0.325490}LAB'],...
             'Rotation',0,'FontName','Good Times Rg','FontSize',12.5);     
%       ylb.Position(2)=min(clb.Limits)-0.04*(max(clb.Limits)-min(clb.Limits));
%       ylb.HorizontalAlignment='right';
      ylb.Units = 'normalized';       
      ylb.Position(1:2) = [2.62*clb.Position(1) -0.1*ylb.Position(2)];
    end
    L=cellfun(@(x)sprintf(config.colorbar_tick_label_notation,x),num2cell(get(clb,'ytick')),'Un',0);    
%     L=cellfun(@(x)sprintf('%0.2f',x),num2cell(get(clb,'ytick')),'Un',0);
    set(clb,'yticklabel',L)   
  else
    ax = gca;
    if ~isempty(config.axis_xtick)
      ax.XTick=config.axis_xtick;
    end
    if ~isempty(config.axis_ytick)    
      ax.YTick=config.axis_ytick;
    end    
    if ~isempty(config.axis_ztick)    
      ax.ZTick=config.axis_ztick;
    end     
    ax.XMinorTick = config.figure_XMinorTick;
    ax.YMinorTick = config.figure_YMinorTick;
    ax.ZMinorTick = config.figure_ZMinorTick;   
    ax.XGrid = config.figure_XGrid;
    ax.YGrid = config.figure_YGrid;
    ax.ZGrid = config.figure_ZGrid;
    ax.XMinorGrid = config.figure_XMinorGrid;
    ax.YMinorGrid = config.figure_YMinorGrid;    
    ax.ZMinorGrid = config.figure_ZMinorGrid;
    if strcmp(config.figure_GridColor,'default')
      ax.GridColor=[0.15 0.15 0.15];
      ax.MinorGridColor=[0.1 0.1 0.1];      
    else
      ax.GridColor = config.figure_GridColor;
      ax.MinorGridColor = config.figure_MinorGridColor;      
    end     
    if strcmp(config.plot_figure_axis,'yes')
      ax.XColor=ax.XColor;
      ax.YColor=ax.YColor;      
    else
      ax.XColor='none';
      ax.YColor='none';      
    end    
    view(-45,25);
  end         

  %set(gcf,'Renderer','painters')   
  %set(gcf,'Renderer','opengl')       
  %set(gcf,'InvertHardcopy','off','Color',[1 1 1])
  if strcmp(config.printer_renderer,'opengl')
%     set(gcf,'InvertHardcopy','off','PaperType','<custom>','PaperSize',[5.83 4.38],'Color',[1 1 1],'Renderer','opengl');
    set(gcf,'InvertHardcopy','off','PaperType','<custom>','PaperSize',[6.5 4.88],'Color',[1 1 1],'Renderer','opengl');
  elseif strcmp(config.printer_renderer,'painters')
%     set(gcf,'InvertHardcopy','off','PaperType','<custom>','PaperSize',[5.83 4.38],'Color',[1 1 1],'Renderer','painters');  
    set(gcf,'InvertHardcopy','off','PaperType','<custom>','PaperSize',[6.5 4.88],'Color',[1 1 1],'Renderer','painters');
  else
    throw_error('In mesh_plot_figure_fem2dQ4.m: config.printer_renderer');
  end  
  %set(gca,'FontName','Times New Roman','FontSize',14,'FontWeight','bold');
  %set(gca,'FontName','Times New Roman','FontSize',13,'FontWeight','bold'); 
  %set(gca,'FontName','Segoe UI Semibold','FontSize',12,'FontWeight','normal');   
  if strcmp(config.front_matlab_plot_view_orientation,'no')
    set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','bold');  
    xlabel('$x_1$','FontWeight','normal','FontSize',18,'FontName',...
           'Times New Roman','Interpreter','latex');
    ylabel('$x_2$','FontWeight','normal','FontSize',18,'FontName',...
           'Times New Roman','Interpreter','latex');
    zlabel(fieldvartitleclb,'FontWeight','normal','FontSize',18,'FontName',...
           'Times New Roman','Interpreter','latex');       
  else
    set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','bold'); 
    xlabel('$x_1$','FontWeight','normal','FontSize',18,'FontName',...
           'Times New Roman','Interpreter','latex');
    ylabel('$x_2$','FontWeight','normal','FontSize',18,'FontName',...
           'Times New Roman','Interpreter','latex');
    zlabel(fieldvartitleclb,'FontWeight','normal','FontSize',18,'FontName',...
           'Times New Roman','Interpreter','latex');    
  end
%       if strcmp(config.front_matlab_plot_view_orientation,'yes')
%         title(titleSolutionY,'FontWeight','bold','FontSize',20,'FontName',...
%               'Times New Roman','Interpreter','latex');
%       end
  if strcmp(config.front_matlab_plot_view_orientation,'no') && strcmp(config.plot_vemlab_logo,'yes') 
%     title(['\color[rgb]{0.360784 0.4 0.435294}VEM','\color[rgb]{0.850980 0.270588 0.325490}LAB'],...
%              'FontName','Good Times Rg','FontSize',11,'Position',[0.8*max(xlim) 0.6*min(ylim) 2.4*max(zlim)]);
    tl=title(['\color[rgb]{0.360784 0.4 0.435294}VEM','\color[rgb]{0.850980 0.270588 0.325490}LAB'],...
             'FontName','Good Times Rg','FontSize',12.5);
    tl.Units = 'normalized';
    tl.Position(1:3) = [1.783*tl.Position(1) tl.Position(2) tl.Position(3)];           
  end
  %axis vis3d
  %axis fill
  %axis square
  %axis tight
  axis normal
%       if strcmp(config.plot_vemlab_logo,'yes') 
%         vemlab_logo_file=[config.icons_folder_location,'vemlab.png']; 
%         Img = imread(vemlab_logo_file); 
%         hold on
%         Lx=max(nodes(:,1))-min(nodes(:,1));
%         Ly=max(nodes(:,2))-min(nodes(:,2));
%         fLx=0.391304348*Lx;
%         fLy=0.045000000*Ly;
%         image([max(nodes(:,1))*0.995-fLx max(nodes(:,1))*0.995],[min(nodes(:,2))-0.8*fLy min(nodes(:,2))-1.8*fLy],Img);    
%         axis([min(nodes(:,1)) max(nodes(:,1)) min(nodes(:,2))-2.65*fLy max(nodes(:,2))])
%         pbaspect([abs(min(nodes(:,1))-max(nodes(:,1))) abs(min(nodes(:,2))-max(nodes(:,2))) 1])   
%         hold on 
%       else
  axis([min(nodescoords(:,1)) max(nodescoords(:,1)) min(nodescoords(:,2)) max(nodescoords(:,2))])
  pbaspect([abs(min(nodescoords(:,1))-max(nodescoords(:,1))) abs(min(nodescoords(:,2))-max(nodescoords(:,2))) 1])
  hold on 
  
  % save figures
  if strcmp(soltype,'numerical')
    if strcmp(config.print_figures,'yes')
      print_figure(config,fieldvartitlefigfile);
    end
    if strcmp(config.save_matlab_figures,'yes')  
      output_figure=[config.matlab_figures_output_folder_location,config.mesh_filename,'_',config.vemlab_method,'_',fieldvartitlefigfile,'.fig'];
      saveas(gcf,output_figure);
    end
  elseif strcmp(soltype,'exact')
    if strcmp(config.print_exact_figures,'yes')
      print_figure(config,fieldvartitlefigfile);   
    end
    if strcmp(config.save_exact_matlab_figures,'yes')  
      output_figure=[config.matlab_figures_output_folder_location,config.mesh_filename,'_',config.vemlab_method,'_',fieldvartitlefigfile,'.fig'];
      saveas(gcf,output_figure);
    end  
  else
    throw_error('In mesh_plot_figure_fem2dQ4.m: solutiontype');
  end  
  
end

