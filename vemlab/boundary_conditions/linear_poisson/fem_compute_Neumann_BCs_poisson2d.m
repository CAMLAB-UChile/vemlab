function Neumann_BCs=...
              fem_compute_Neumann_BCs_poisson2d(domainMesh,Neumann_boundary_nodes,...
                                              Neumann_boundary_dofs,...
                                              Neumann_fun_value,matProps)
%
% A Neumann BC for the Poisson problem is a flux-type condition, which is
% currently not implemented
%

%   throw_error('In fem_compute_Neumann_BCs_poisson2d.m: A Neumann BC for the Poisson problem is a flux-type condition, which is currently not implemented\n');
                         
  coords=domainMesh.coords;
  nodes=Neumann_boundary_nodes;
  % FEM interpolation matrix for 1-pt rule: N(xi=0)
  N=[0.5 0.5];
  % nodal coordinates of right boundary nodes
  x=coords(nodes,1);
  y=coords(nodes,2);
  % compute Neumman BCs values using a traction vector
  % NOTE: Nodes need to be ordered consecutively on the Neumman edge. A
  % better approach would be to have the edge connectivities available.
  num_nodes=length(nodes);
  num_edges=num_nodes-1;
  Neumann_BCs.values=zeros(num_nodes,1);
  for edge=1:num_edges
    first_node=edge;
    second_node=edge+1;
    x1=x(first_node); y1=y(first_node);
    x2=x(second_node); y2=y(second_node);      
    edge_length=norm([x1-x2; y1-y2]);
    xg=0.5*x1+0.5*x2; % isoparametric mapping at xi=0: xg=N1(0)*x1+N2(0)*x2
    yg=0.5*y1+0.5*y2; % isoparametric mapping at xi=0: yg=N1(0)*y1+N2(0)*y2     
    qn=Neumann_fun_value.qn(xg,yg,matProps); 
    range1=first_node;
    range2=second_node;
    range=[range1,range2];
    % 1-pt rule: N(xi=0)'*traction(xi=0)*(Le/2)*2
    Neumann_BCs.values(range)=Neumann_BCs.values(range)+edge_length*N'*qn;     
  end
  % indices of the dofs that have an associated Neumman BC
  Neumann_BCs.indexes=Neumann_boundary_dofs;
end