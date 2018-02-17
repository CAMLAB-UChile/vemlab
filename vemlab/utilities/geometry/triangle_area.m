function area = triangle_area(nodes_coords)
  % nodes_coords = coordinates of the nodes that define the 3-node triangle.
  % Nodes are assumed to be ordered counterclockwise or clockwise
  area=1.0/2.0*abs(nodes_coords(1,1)*(nodes_coords(3,2)-nodes_coords(2,2))+...
       nodes_coords(3,1)*(nodes_coords(2,2)-nodes_coords(1,2))+...
       nodes_coords(2,1)*(nodes_coords(1,2)-nodes_coords(3,2)));
end

