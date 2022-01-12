function Neumann_BCs=...
           vem_compute_normal_pressure_Neumann_BCs_linelast2d(domainMesh,Neumann_boundary_nodes,...
                                                              Neumann_boundary_dofs,...
                                                              Neumann_fun_value,matProps)
  coords=domainMesh.coords;  
  nodes=Neumann_boundary_nodes;
  
  % VEM interpolation matrix
  N_bar=[0.5 0.0;...
         0.0 0.5];
  % nodal coordinates of right boundary nodes
  x=coords(nodes,1);
  y=coords(nodes,2);
  % compute Neumman BCs values using a traction vector
  % NOTE: Nodes need to be ordered consecutively on the Neumman edge. A
  % better approach would be to have the edge connectivities available.
  
  num_nodes=length(nodes);
  num_edges=num_nodes-1;
  Neumann_BCs.values=zeros(2*num_nodes,1);    
  for edge=1:num_edges
    first_node=edge;
    second_node=edge+1;
    x1=x(first_node); y1=y(first_node);
    x2=x(second_node); y2=y(second_node);    
    edge_length=norm([x1-x2; y1-y2]);    
    edg=[x2 y2 0]-[x1 y1 0];
    out_of_screen=[0 0 1]; % this assumes the 2D geometry is on the x-y system
    xp=x1; yp=y1; % these can be any function... for the current implementation normal_pressure is constant, so we assing any value, for instance, the coordinates of the firt node
    normal_pressure=Neumann_fun_value.normal_pressure(xp,yp,matProps);
    [tx,ty]=apply_normal_pressure_to_polygon_edge(edg,out_of_screen,normal_pressure);    

%     tx_first_node=tx;
%     ty_first_node=ty;
%     tx_second_node=tx;      
%     ty_second_node=ty;
%     traction_first_node=[tx_first_node; ty_first_node];
%     traction_second_node=[tx_second_node; ty_second_node];       
%     mean_traction=(traction_first_node+traction_second_node)/2;
    mean_traction=[tx;ty];
    
%     nodes(first_node)
%     nodes(second_node)
%     norm(mean_traction)
    
    range1=(2*first_node-1):(2*first_node);
    range2=(2*second_node-1):(2*second_node);
  % The following is used in, for instance, Beirao da Veiga's Basic VEM paper for k=1 (mean_traction can also be computed as: mean_traction = 1/|e| * int_e f ds)
    Neumann_BCs.values(range1)=Neumann_BCs.values(range1)+edge_length*N_bar'*mean_traction;
    Neumann_BCs.values(range2)=Neumann_BCs.values(range2)+edge_length*N_bar'*mean_traction;  
%       This also works: if shape functions are defined as in 1D FEM, the following coincides with a trapezoidal rule
%       Neumann_BCs.values(range1)=Neumann_BCs.values(range1)+edge_length*N_bar'*traction_first_node;
%       Neumann_BCs.values(range2)=Neumann_BCs.values(range2)+edge_length*N_bar'*traction_second_node;    

%     x1
%     y1
%     traction_first_node
%     x2
%     y2
%     traction_second_node

  end
  % indices of the dofs that have an associated Neumman BC
  Neumann_BCs.indexes=Neumann_boundary_dofs;  
   
%   Neumann_BCs.values
  
end

function [tx,ty] = apply_normal_pressure_to_polygon_edge(edge1,edge2,normal_pressure)

  % unit outward normal to polygon edge
  nx=edge1(2)*edge2(3)-edge1(3)*edge2(2);
  ny=edge1(3)*edge2(1)-edge1(1)*edge2(3);
  normal=[nx ny];
  mag=norm(normal);
  normals=normal/mag;  
  % traction components
  tx=normal_pressure*normals(1);
  ty=normal_pressure*normals(2);
  
end