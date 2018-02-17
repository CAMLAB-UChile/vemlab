function [Wc] = vem_Wc_linelast2d(verts)

  area_components=verts(:,1).*verts([2:end,1],2)-verts([2:end,1],1).* verts(:,2);
  area=0.5*abs(sum(area_components));

  v1=[verts; verts(1,:)];

  n=normals_to_edges(verts);
  n1=n(:,1);
  n2=n(:,2);
  n11=n1([length(verts),1:end-1]);
  n22=n2([length(verts),1:end-1]);

  len=zeros(length(verts),1);
  for i=1:length(verts)
      len(i)=norm(v1(i,:)-v1(i+1,:));
  end

  len1=len([length(verts),1:end-1]);
  Wc=zeros(2*length(verts),3);
  for i=1:length(verts)
      q1a=0.25*(len1(i)*n11(i)+len(i)*n1(i))/area;
      q2a=0.25*(len1(i)*n22(i)+len(i)*n2(i))/area; 
      Wc_a=[2*q1a,0,q2a; 0,2*q2a,q1a];
      Wc((2*i-1):2*i,:)=Wc_a;
  end

end 