function h_max = max_edge_size(verts)

  mysize=size(verts,1);
  mydist=zeros(mysize,1);
  
  % shifts coordinates with respect to the first node in the element connectivity
  shifted_coords= verts-repmat(verts(1,:),mysize,1);
  
  % compute the distances wrt to the first node in the element connectivity
  for node = 1:mysize
    mydist(node) = norm(shifted_coords(node,:));
  end
  
  % maximum distance from the first node to any node in the element connectivity 
  h_max=max(mydist);
  
end