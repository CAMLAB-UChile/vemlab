function h_max = max_edge_size(verts)
  mysize=size(verts,1);
  h_max=max(norm(verts-repmat(verts(1,:),mysize,1)));
end