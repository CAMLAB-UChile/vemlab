function [Hc] = vem_Hc_poisson2d(verts)
  verts1=verts(:,1);
  verts2=verts(:,2);
  N=length(verts);
  x_bar=(1/N)*[sum(verts1),sum(verts2)];
  Hc=zeros(N,2);
  for i=1:N
    Hc(i,:)=[verts1(i)-x_bar(1),verts2(i)-x_bar(2)];
  end
end 