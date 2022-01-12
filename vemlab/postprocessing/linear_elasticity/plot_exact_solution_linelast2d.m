%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                     VEMLab
%           Source code  : http://camlab.cl/research/software/vemlab/
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:              plot_exact_solution_linelast2d 
%
% Created by : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
% Updated by :
% Updated by :
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% Plot displacement, stress and strain field exact solutions.
%
% Usage
% =====
% plot_exact_solution_linelast2d(domainMesh,exact_solution_handle,...
%                                              matProps,config)
%
% Input
% =====
% domainMesh : structure containing mesh data (coords,connect,etc.)
% exact_solution_handle : handle to exact solution functions
% matProps : structure containing the material properties
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
% Jan. 04, 2021: plotting functionalities improved
% Aug. 21, 2021: function became under simplification so that to plot only at the
%                nodal locations (by A. Ortiz-Bernardin)
% Oct. 20, 2018: add option to switch off all matlab contour plots (by A. Ortiz-Bernardin)
%                add option to deform the domain in matlab coontour plots (by A. Ortiz-Bernardin)
% Apr. 19, 2018: improve the plotting of axis and fonts (by A. Ortiz-Bernardin)
% Mar. 17, 2018: first realease (by A. Ortiz-Bernardin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_exact_solution_linelast2d(domainMesh,exact_solution_handle,...
                                        matProps,config)
  if strcmp(config.create_matlab_exact_contour_plots,'yes')
    [displacement,stress,strain] =...
        exact_nodal_solutions_linelast2d(exact_solution_handle,domainMesh,matProps,config);
    patch_plot_exact_nodal_solutions(domainMesh,displacement,stress,strain,config);
  end
end

function [displacement,stress,strain] =...
       exact_nodal_solutions_linelast2d(exact_solution_handle,domainMesh,matProps,config) 
  fprintf('\n');
  fprintf('Computing exact solutions at the nodes...\n');
  
  coords=domainMesh.coords;  
  num_nodes=size(coords,1);  
  displacement=zeros(2*num_nodes,1);
  
  % loop over elements
  for I=1:num_nodes
    xy=domainMesh.coords(I,:);
    % exact solution
    [disp_exact,p_exact,strainvec_exact]=exact_solutions_linelast2d(exact_solution_handle,xy,matProps);
    % exact displacement
    range=[2*I-1;2*I];
    displacement(range)=disp_exact;
       
    if strcmp(config.linelast2d_plot_stress_and_strain,'yes')
      % exact strain
      str_xgp=[strainvec_exact(1);strainvec_exact(2);...
               2*strainvec_exact(3)];  % [e11,e22,2*e12]   
      if strcmp(matProps.plane_state,'plane_stress')
        eten_xgp=[str_xgp(1),str_xgp(3)/2,0;... % for tensor representation divide 
                  str_xgp(3)/2,str_xgp(2),0;...   % strain by 2: (2*e12)/2
                  0,0,-matProps.nu*(str_xgp(1)+str_xgp(2))];  
      elseif strcmp(matProps.plane_state,'plane_strain')
        eten_xgp=[str_xgp(1),str_xgp(3)/2,0;...
                  str_xgp(3)/2,str_xgp(2),0;...
                  0,0,0];      
      else
        throw_error('Error in plot_exact_solution_linelast2d.m --> exact_nodal_solutions_linelast2d: plane_state\n');
      end
      % check whether there are NaN in the exact derivatives
      [indi,indj]=find(isnan(eten_xgp));  
      len_indi=length(indi);
      len_indj=length(indj);
      if len_indi>0 || len_indj>0
        throw_warning('In plot_exact_solution_linelast2d.m --> exact_nodal_solutions_linelast2d: Exact strains are not defined at one or more nodes. Interpret results carefully.')
      end
      for ig=1:len_indi
        for ik=1:len_indj
          eten_xgp(indi(ig),indj(ik))=0;
        end
      end
      % store exact strain
      pstrain_xgp=sort(eig(eten_xgp),'descend'); % principal strain: pstrain1 > pstrain2 > pstrain3      
      strain.e11(I)=eten_xgp(1,1);
      strain.e12(I)=eten_xgp(1,2);   
      strain.e22(I)=eten_xgp(2,2); 
      strain.e33(I)=eten_xgp(3,3);  
      strain.e1(I)=pstrain_xgp(1);
      strain.e2(I)=pstrain_xgp(2);
      strain.e3(I)=pstrain_xgp(3);        
      % exact stress
      D=matProps.D;    % This D is defined as per VEM for strain = [e11 e22 e12].
      D(3,3)=D(3,3)/4; % To avoid troubles, come back to the standard FEM definition of D 
                       % to obtain the stress by multiplying with [e11 e22 2*e12]
      svec_xgp=D*[eten_xgp(1,1);eten_xgp(2,2);2*eten_xgp(1,2)];
      if strcmp(matProps.plane_state,'plane_stress')
        sten_xgp=[svec_xgp(1),svec_xgp(3),0;svec_xgp(3),svec_xgp(2),0;0,0,0];
      elseif strcmp(matProps.plane_state,'plane_strain')
        sten_xgp=[svec_xgp(1),svec_xgp(3),0;svec_xgp(3),svec_xgp(2),0;0,0,matProps.nu*(svec_xgp(1)+svec_xgp(2))];      
      else
        throw_error('Error in plot_exact_solution_linelast2d.m --> exact_nodal_solutions_linelast2d: plane_state\n');      
      end
      pstress_h = sort(eig(sten_xgp),'descend');  % principal stresses: pstress1 > pstress2 > pstress3  
      sxx_xgp=sten_xgp(1,1); syy_xgp=sten_xgp(2,2); szz_xgp=sten_xgp(3,3);
      sxy_xgp=sten_xgp(1,2); sxz_xgp=0.0; syz_xgp=0.0;
      VM_xgp=(1/sqrt(2))*sqrt((sxx_xgp-syy_xgp)^2+(syy_xgp-szz_xgp)^2+(szz_xgp-sxx_xgp)^2+6*(sxy_xgp^2+syz_xgp^2+sxz_xgp^2));    
      % store exact stress
      stress.s11(I)=sten_xgp(1,1);
      stress.s12(I)=sten_xgp(1,2);   
      stress.s22(I)=sten_xgp(2,2); 
      stress.s33(I)=sten_xgp(3,3);   
      stress.s1(I)=pstress_h(1);
      stress.s2(I)=pstress_h(2);
      stress.s3(I)=pstress_h(3);  
      stress.vm(I)=VM_xgp;
      %Ey=matProps.Ey;
      %nu=matProps.nu;
      %lam=Ey*nu/((1+nu)*(1-2*nu));
      stress.p(I)=p_exact(1);
    else
      stress=[];
      strain=[]; 
    end
  end
end

function patch_plot_exact_nodal_solutions(domainMesh,displacements,stress,strain,config)
  nodes=domainMesh.coords;
  nodesNumber=size(nodes,1);  
  polygons=domainMesh.connect;
  range=1:nodesNumber;
  displacementsX=displacements(2*range-1);
  displacementsY=displacements(2*range);  
  displacementsNorm=sqrt(displacementsX.*displacementsX+displacementsY.*displacementsY);

  fprintf('\n'); 
  fprintf('Plotting exact displacements...\n');
  
  if strcmp(config.linelast2d_plot_deformed_domain,'yes') 
    scale=config.linelast2d_scale_for_plotting_deformed_domain;
    nodes_dispx=[nodes(:,1)+displacementsX*scale, nodes(:,2)];
    nodes_dispy=[nodes(:,1), nodes(:,2)+displacementsY*scale];
    nodes=[nodes(:,1)+displacementsX*scale, nodes(:,2)+displacementsY*scale];
    
    if strcmp(config.linelast2d_plot_displacement.unorm,'yes')
      mytitleclb='$\scriptstyle{\|}\displaystyle{\mathbf{u}}\scriptstyle{\|}$';
      mytitlefigfile='u';    
      colorbar_limits=config.linelast2d_plot_displacement.clim.unorm;       
      patch_plot_figure(config,'on_nodes',polygons,nodes,displacementsNorm,mytitleclb,mytitlefigfile,colorbar_limits,0.79,'exact');
    end

    if strcmp(config.linelast2d_plot_displacement.ux,'yes')  
      mytitleclb='$u_{1}$';  
      mytitlefigfile='u1';  
      colorbar_limits=config.linelast2d_plot_displacement.clim.ux;       
      patch_plot_figure(config,'on_nodes',polygons,nodes_dispx,displacementsX,mytitleclb,mytitlefigfile,colorbar_limits,0.65,'exact');
    end

    if strcmp(config.linelast2d_plot_displacement.uy,'yes')
      mytitleclb='$u_{2}$';  
      mytitlefigfile='u2';  
      colorbar_limits=config.linelast2d_plot_displacement.clim.uy;       
      patch_plot_figure(config,'on_nodes',polygons,nodes_dispy,displacementsY,mytitleclb,mytitlefigfile,colorbar_limits,0.65,'exact');
    end    
  else
    if strcmp(config.linelast2d_plot_displacement.unorm,'yes')
      mytitleclb='$\scriptstyle{\|}\displaystyle{\mathbf{u}}\scriptstyle{\|}$';
      mytitlefigfile='u'; 
      colorbar_limits=config.linelast2d_plot_displacement.clim.unorm;       
      patch_plot_figure(config,'on_nodes',polygons,nodes,displacementsNorm,mytitleclb,mytitlefigfile,colorbar_limits,0.79,'exact');
    end

    if strcmp(config.linelast2d_plot_displacement.ux,'yes')  
      mytitleclb='$u_{1}$'; 
      mytitlefigfile='u1';  
      colorbar_limits=config.linelast2d_plot_displacement.clim.ux;       
      patch_plot_figure(config,'on_nodes',polygons,nodes,displacementsX,mytitleclb,mytitlefigfile,colorbar_limits,0.65,'exact');
    end

    if strcmp(config.linelast2d_plot_displacement.uy,'yes')
      mytitleclb='$u_{2}$'; 
      mytitlefigfile='u2'; 
      colorbar_limits=config.linelast2d_plot_displacement.clim.uy;       
      patch_plot_figure(config,'on_nodes',polygons,nodes,displacementsY,mytitleclb,mytitlefigfile,colorbar_limits,0.65,'exact');
    end     
  end
  
  if ~isempty(stress)
    fprintf('Plotting exact stresses and strains...\n');   

    if strcmp(config.linelast2d_plot_stress.s11,'yes')
      mytitleclb='$\sigma_{11}$';
      mytitlefigfile='s11';  
      colorbar_limits=config.linelast2d_plot_stress.clim.s11;        
      patch_plot_figure(config,'on_nodes',polygons,nodes,stress.s11',mytitleclb,mytitlefigfile,colorbar_limits,0.87,'exact');
    end

    if strcmp(config.linelast2d_plot_stress.s22,'yes')  
      mytitleclb='$\sigma_{22}$';
      mytitlefigfile='s22'; 
      colorbar_limits=config.linelast2d_plot_stress.clim.s22;        
      patch_plot_figure(config,'on_nodes',polygons,nodes,stress.s22',mytitleclb,mytitlefigfile,colorbar_limits,0.87,'exact');
    end 

    if strcmp(config.linelast2d_plot_stress.s12,'yes')  
      mytitleclb='$\sigma_{12}$';
      mytitlefigfile='s12';
      colorbar_limits=config.linelast2d_plot_stress.clim.s12;        
      patch_plot_figure(config,'on_nodes',polygons,nodes,stress.s12',mytitleclb,mytitlefigfile,colorbar_limits,0.87,'exact');
    end       

    if strcmp(config.linelast2d_plot_stress.s33,'yes')  
      mytitleclb='$\sigma_{33}$';
      mytitlefigfile='s33';
      colorbar_limits=config.linelast2d_plot_stress.clim.s33;        
      patch_plot_figure(config,'on_nodes',polygons,nodes,stress.s33',mytitleclb,mytitlefigfile,colorbar_limits,0.87,'exact');
    end  

    if strcmp(config.linelast2d_plot_stress.s1,'yes')  
      mytitleclb='$\sigma_{1}$';
      mytitlefigfile='s1'; 
      colorbar_limits=config.linelast2d_plot_stress.clim.s1;        
      patch_plot_figure(config,'on_nodes',polygons,nodes,stress.s1',mytitleclb,mytitlefigfile,colorbar_limits,0.65,'exact');
    end   

    if strcmp(config.linelast2d_plot_stress.s2,'yes')  
      mytitleclb='$\sigma_{2}$';
      mytitlefigfile='s2'; 
      colorbar_limits=config.linelast2d_plot_stress.clim.s2;        
      patch_plot_figure(config,'on_nodes',polygons,nodes,stress.s2',mytitleclb,mytitlefigfile,colorbar_limits,0.65,'exact');
    end 

    if strcmp(config.linelast2d_plot_stress.s3,'yes')  
      mytitleclb='$\sigma_{3}$';
      mytitlefigfile='s3'; 
      colorbar_limits=config.linelast2d_plot_stress.clim.s3;        
      patch_plot_figure(config,'on_nodes',polygons,nodes,stress.s3',mytitleclb,mytitlefigfile,colorbar_limits,0.65,'exact');
    end    

    if strcmp(config.linelast2d_plot_stress.vm,'yes')  
      mytitleclb='$\sigma_{\mathrm{vm}}$';
      mytitlefigfile='vm';
      colorbar_limits=config.linelast2d_plot_stress.clim.vm;        
      patch_plot_figure(config,'on_nodes',polygons,nodes,stress.vm',mytitleclb,mytitlefigfile,colorbar_limits,1.04,'exact');
    end 

    if strcmp(config.linelast2d_plot_stress.p,'yes')  
      mytitleclb='$p$';
      mytitlefigfile='p'; 
      colorbar_limits=config.linelast2d_plot_stress.clim.p;        
      patch_plot_figure(config,'on_nodes',polygons,nodes,stress.p',mytitleclb,mytitlefigfile,colorbar_limits,0.42,'exact');
    end  

    if strcmp(config.linelast2d_plot_strain.e11,'yes')  
      mytitleclb='$\varepsilon_{11}$';
      mytitlefigfile='e11';
      colorbar_limits=config.linelast2d_plot_strain.clim.e11;       
      patch_plot_figure(config,'on_nodes',polygons,nodes,strain.e11',mytitleclb,mytitlefigfile,colorbar_limits,0.78,'exact');
    end      

    if strcmp(config.linelast2d_plot_strain.e12,'yes')  
      mytitleclb='$\varepsilon_{12}$';
      mytitlefigfile='e12';
      colorbar_limits=config.linelast2d_plot_strain.clim.e12;       
      patch_plot_figure(config,'on_nodes',polygons,nodes,strain.e12',mytitleclb,mytitlefigfile,colorbar_limits,0.78,'exact');
    end    

    if strcmp(config.linelast2d_plot_strain.e22,'yes')  
      mytitleclb='$\varepsilon_{22}$';
      mytitlefigfile='e22'; 
      colorbar_limits=config.linelast2d_plot_strain.clim.e22;       
      patch_plot_figure(config,'on_nodes',polygons,nodes,strain.e22',mytitleclb,mytitlefigfile,colorbar_limits,0.78,'exact');
    end  

    if strcmp(config.linelast2d_plot_strain.e33,'yes')  
      mytitleclb='$\varepsilon_{33}$';
      mytitlefigfile='e33'; 
      colorbar_limits=config.linelast2d_plot_strain.clim.e33;       
      patch_plot_figure(config,'on_nodes',polygons,nodes,strain.e33',mytitleclb,mytitlefigfile,colorbar_limits,0.78,'exact');
    end      

    if strcmp(config.linelast2d_plot_strain.e1,'yes')  
      mytitleclb='$\varepsilon_{1}$';
      mytitlefigfile='e1';
      colorbar_limits=config.linelast2d_plot_strain.clim.e1;       
      patch_plot_figure(config,'on_nodes',polygons,nodes,strain.e1',mytitleclb,mytitlefigfile,colorbar_limits,0.56,'exact');
    end  

    if strcmp(config.linelast2d_plot_strain.e2,'yes')  
      mytitleclb='$\varepsilon_{2}$';
      mytitlefigfile='e2'; 
      colorbar_limits=config.linelast2d_plot_strain.clim.e2;       
      patch_plot_figure(config,'on_nodes',polygons,nodes,strain.e2',mytitleclb,mytitlefigfile,colorbar_limits,0.56,'exact');
    end  

    if strcmp(config.linelast2d_plot_strain.e3,'yes')  
      mytitleclb='$\varepsilon_{3}$';
      mytitlefigfile='e3';
      colorbar_limits=config.linelast2d_plot_strain.clim.e3;       
      patch_plot_figure(config,'on_nodes',polygons,nodes,strain.e3',mytitleclb,mytitlefigfile,colorbar_limits,0.56,'exact');
    end  
  end
end

