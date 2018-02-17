function [Hr]=vem_Hr_linelast2d(verts)
  verts1=verts(:,1);
  verts2=verts(:,2);
  x_bar=(1/length(verts))*[sum(verts1),sum(verts2)];
  Hr=zeros(2*length(verts),3);
  for i=1:length(verts)
      Hr_a=[1,0,verts2(i)-x_bar(2);...
            0,1,-verts1(i)+x_bar(1)];
      Hr((2*i-1):2*i,:)=Hr_a;
  end
end 