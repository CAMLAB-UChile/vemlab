function [n1] = normals_to_edges(verts)
  v1=[verts; verts(1,:)];
  for i=1:length(verts)
      w(i,:)=v1(i+1,:)-v1(i,:);
      n=[w(:,2) -w(:,1)];
      n1(i,:)=n(i,:)/norm(n(i,:));
  end
end 