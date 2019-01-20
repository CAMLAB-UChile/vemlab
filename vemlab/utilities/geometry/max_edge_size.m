function h_max = max_edge_size(verts)

  mysize=size(verts,1);
  mydist=zeros(mysize*mysize,1);
  
  % This procedure also compute the distance between a node with itself (0); however, it
  % doesn't matter since at the end the maximum distance is selected
  r = 1;
  for node_i = 1:mysize
    for node_j = 1:mysize
      % shifts coordinates with respect to "node i" in the connectivity
      shifted_coords= verts-repmat(verts(node_i,:),mysize,1);
      % compute the distances wrt to "node_i" in the element connectivity
      for node = 1:(mysize-1)
        mydist(r) = norm(shifted_coords(node,:));
        r = r + 1;
      end
    end  
  end
  
  % maximum distance from the first node to any node in the element connectivity 
  h_max=max(mydist);
  
end