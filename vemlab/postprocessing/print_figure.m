function print_figure(config,fieldvartitlefigfile)

  % painters renderer has many issues for plotting pdf and eps figures
  % (e.g., white lines on patch plots, fonts substitution, etc. The code
  % below can deal with the font substitution by saving a svg file and then
  % using inkscape to convert to pdf. This will resolve the font substitution, 
  % but still remains the problem of white lines. In addition, this workaround 
  % is very costly as large images will take too much time in the process
  % of saving the svg file and then converting them to pdf. So, the default
  % will be to use the 'opengl' renderer and the 'painters' renderer
  % will remain here only for experimental purposes.
  
  if strcmp(config.printer_renderer,'painters')
    % We print to a svg file, which will be converted to a PDF file afterwards using inkscape.
    % This is done because Matlab's implementation to save as eps and pdf
    % using painters renderer only accepts few fonts and when the fonts
    % are not within these fonts they are changed. Using inkscape to
    % convert the svg file to a pdf file will keep the special fonts (in
    % this case, we use Good Times Rg and Segoe UI Semibold).    
    output_figure=[config.matlab_figures_output_folder_location,...
                   config.mesh_filename,'_',config.vemlab_method,...
                   '_',fieldvartitlefigfile,'.svg'];
    print(['-r',num2str(config.print_figures_resolution)],output_figure,'-dsvg');
    current_folder = pwd; % store the current directory
    output_figure_noext=[config.matlab_figures_output_folder_location,...
                         config.mesh_filename,'_',config.vemlab_method,...
                         '_',fieldvartitlefigfile];      
    command = strcat(config.path_to_inkscape, {' '}, '"',...
                     output_figure_noext, '.svg" --export-pdf="',...
                     output_figure_noext,...
                     '.pdf" --export-pdf-version=1.5 --export-area-drawing --export-dpi="',...
                     num2str(config.pdf_figures_resolution));                   
    command = string(command);
    [status,msg] = system(command);   
    cd(config.matlab_figures_output_folder_location); 
    delete *.svg
    cd(current_folder);
  elseif strcmp(config.printer_renderer,'opengl')
    % pdf
    output_figure=[config.matlab_figures_output_folder_location,...
                   config.mesh_filename,'_',config.vemlab_method,'_',...
                   fieldvartitlefigfile,'.pdf']; 
    print(['-r',num2str(config.print_figures_resolution)],output_figure,'-dpdf','-bestfit'); 
%     % png
%     output_figure=[config.matlab_figures_output_folder_location,...
%                    config.mesh_filename,'_',config.vemlab_method,'_',...
%                    fieldvartitlefigfile,'.png']; 
%     print(['-r',num2str(config.print_figures_resolution)],output_figure,'-dpng');     
%     % eps
%     output_figure=[config.matlab_figures_output_folder_location,...
%                    config.mesh_filename,'_',config.vemlab_method,'_',...
%                    fieldvartitlefigfile,'.eps']; 
%     print(['-r',num2str(config.print_figures_resolution)],output_figure,'-depsc');
  else
    throw_error('In print_figure: config.printer_renderer');
  end


end

