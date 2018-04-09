function Neumann_BCs=...
           vem_compute_Neumann_BCs_poisson2d(domainMesh,Neumann_boundary_nodes,...
                                              Neumann_boundary_dofs,...
                                              Neumann_fun_values)
%
% A Neumann BC for the Poisson problem is a flux-type condition, which is
% currently not implemented
%

  throw_error('In vem_compute_Neumann_BCs_poisson2d.m: A Neumann BC for the Poisson problem is a flux-type condition, which is currently not implemented\n');

%   coords=mesh.coords;  
%   nodes=Neumann_boundary_nodes;  
%   % VEM interpolation matrix
%   N_bar=[0.5 0.0;...
%          0.0 0.5];
%   % nodal coordinates of right boundary nodes
%   x=coords(nodes,1);
%   y=coords(nodes,2);
%   % compute Neumman BCs values using a traction vector
%   % NOTE: Nodes need to be ordered consecutively on the Neumman edge. A
%   % better approach would be to have the edge connectivities available.
%   num_nodes=length(nodes);
%   num_edges=num_nodes-1;
%   Neumann_BCs.values=zeros(2*num_nodes,1);    
%   for edge=1:num_edges
%     first_node=edge;
%     second_node=edge+1;
%     x1=x(first_node); y1=y(first_node);
%     x2=x(second_node); y2=y(second_node);      
%     edge_length=norm([x1-x2; y1-y2]);
%     tx_first_node=Neumann_fun_values.fx(x1,y1);
%     ty_first_node=Neumann_fun_values.fy(x1,y1);
%     tx_second_node=Neumann_fun_values.fx(x2,y2);      
%     ty_second_node=Neumann_fun_values.fy(x2,y2);
%     traction_first_node=[tx_first_node; ty_first_node];
%     traction_second_node=[tx_second_node; ty_second_node];       
%     mean_traction=(traction_first_node+traction_second_node)/2;
%     range1=(2*first_node-1):(2*first_node);
%     range2=(2*second_node-1):(2*second_node);
%   % The following is used in, for instance, Beirao da Veiga's Basic VEM paper for k=1 (mean_traction can also be computed as: mean_traction = 1/|e| * int_e f ds)
%     Neumann_BCs.values(range1)=Neumann_BCs.values(range1)+edge_length*N_bar'*mean_traction;
%     Neumann_BCs.values(range2)=Neumann_BCs.values(range2)+edge_length*N_bar'*mean_traction;  
% %       This also works: if shape functions are defined as in 1D FEM, the following coincides with a trapezoidal rule
% %       Neumann_BCs.values(range1)=Neumann_BCs.values(range1)+edge_length*N_bar'*traction_first_node;
% %       Neumann_BCs.values(range2)=Neumann_BCs.values(range2)+edge_length*N_bar'*traction_second_node;      
%   end
%   % indices of the dofs that have an associated Neumman BC
  Neumann_BCs.indexes=Neumann_boundary_dofs;
end