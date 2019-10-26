function [Wc] = vem_Wc_poisson2d(verts)

%   area_components=verts(:,1).*verts([2:end,1],2)-verts([2:end,1],1).* verts(:,2);
%   area=0.5*abs(sum(area_components));
  area=polyarea(verts(:,1),verts(:,2));

  v1=[verts; verts(1,:)];

  N=length(verts);
  ne=normals_to_edges(verts);
  n1=ne(:,1);
  n2=ne(:,2);
  n11=n1([N,1:end-1]);
  n22=n2([N,1:end-1]);

  len=zeros(N,1);
  for i=1:N
      len(i)=norm(v1(i,:)-v1(i+1,:));
  end

  len1=len([N,1:end-1]);
  Wc=zeros(N,2);
  for i=1:N
      q1a=0.25*(len1(i)*n11(i)+len(i)*n1(i))/area;
      q2a=0.25*(len1(i)*n22(i)+len(i)*n2(i))/area; 
      Wc(i,:)=[2*q1a,2*q2a];
  end

end 