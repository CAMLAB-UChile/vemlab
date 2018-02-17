function a = diameter(verts)

for i = 1:length(verts)
    for j = 1:length(verts)
        n(i,:) = norm(verts(i,:) - verts(j,:));
        m(i) = max(n(i,:));
    end
end

a = max(m);
end