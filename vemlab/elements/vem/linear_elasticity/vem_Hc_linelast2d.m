function [Hc] = vem_Hc_linelast2d(verts)
  verts1=verts(:,1);
  verts2=verts(:,2);
  x_bar=(1/length(verts))*[sum(verts1),sum(verts2)];
  Hc=zeros(2*length(verts),3);
  for i=1:length(verts)
      Hc_a=[verts1(i)-x_bar(1),0,verts2(i)-x_bar(2);...
            0,verts2(i)-x_bar(2),verts1(i)-x_bar(1)];
      Hc((2*i-1):2*i,:)=Hc_a;
  end
end 