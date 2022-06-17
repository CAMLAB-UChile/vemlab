
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %              Created by  : Rodrigo Silva                                     
  %              E-mail      : rsilvav@outlook.com                                        
  %              Version     : 2.0                         
  %              Date        : June 27, 2018  
  %
  %           Department of Mechanical Engineering                        
  %           University of Chile, Santiago, CHILE    
  %
  %           PURPOSE: mesh from triangles to polygons
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
 
%ENTRADA: 

%Conectividad de los elementos T3:
%Element = [n1 n2 n3; tag2 n3 n4 n5; ...]

%Coordenadas de los nodos:
%Nodos = [x1 y1; x2 y2; ...]

%Conectividad del borde del dominio:
%borde = [nod1 Nod2; Nod1 Nod2; ...]

%esquinas del dominio
%corners = [x1 y1; x2 y2; ...]

%SALIDA:
%Vertices = [x1 y1; x2 y2; ...] coordenadas de los vertices
% V = [v1; v2; v3; v4;]
%     [v5; v6; v7; v8;]... celdas con conectividad de los poligonos

%En esta version los bordes se adaptan a los lados de los poligonos,
%ignorando los nodos del borde, es decir, los nodos del borde ya no son
%tratamos como vertices en el borde (sí los nodos de las esquinas)


function [V,Vertices,Dirichlet_vertices,Neumann_vertices] =...
     T3toPolyMesh(Element,Nodos,borde,corners,num_dirichlet_entities,num_neumann_entities)

  NumberElement = length(Element);
  NumberNodes = length(Nodos(:,1));

  %CENTRO DE CADA ELEMENTO T3
  %UTILIZADOS COMO VERTICES DE LOS POLIGONOS
  centro = zeros(NumberElement,2);
  for i=1:NumberElement   
    centro(i,:) = ( Nodos(Element{i}(1),:) + Nodos(Element{i}(2),:) + Nodos(Element{i}(3),:) )/3;       
  end

  %GENERA CONECTIVIDAD DE LOS POLIGONOS EN BASE A VERTICES (CENTRO DE CADA T3).
  %EL TAG DE CADA NUEVA CELDA CORRESPONDE AL TAG DEL SU NODO CENTRAL 
  %HASTA AQUI FALTAN LOS VERTICES DEL CONTORNO

%  REPLACED BY THE CODE BELOW BY TONGRUI LIU
%   tic
%   V = cell(NumberNodes,1);
%   for i=1:NumberNodes 
%     c = 1;
%     for j=1:NumberElement
%       if(Element{j}(1)==i || Element{j}(2)==i || Element{j}(3)==i)
%         V{i,1}(c) = j;
%         c=c+1;
%       end       
%     end       
%   end
%   toc
  
%   % By Tongrui Liu
  %tic
  connective = reshape(cell2mat(Element),3,[]);
  V = cell(NumberNodes,1);
  tic
  for i=1:NumberNodes
    index = sum(reshape(reshape(connective,[],1)==i*ones(3*length(connective),1),3,[]));
    V{i}= find(index == 1);
  end  
  %toc
  
  %SE GERENA NUEVOS VERTICES EN EL BORDE, Y SE LE ASIGNA UN TAG
  %CADA NUEVO VERTICE ESTA ASOCIADO AL ELEMENTO T3 CORRESPONDINTE Y A SUS
  %NODOS CENTRALES
  Dirichlet_vertices = zeros(num_dirichlet_entities,2);
  Neumann_vertices = zeros(num_neumann_entities,2);
  k = 1; r = 1;
  Cborde=zeros(size(borde,1),5);
  c = length(centro) +1;
  for i=1:size(borde,1)
    % [NuevoTag Nodo1 Nodos2 Cx Cy]  
    Cborde(i,:) = [ c borde(i,1) borde(i,2) (Nodos(borde(i,1),:)+Nodos(borde(i,2),:))/2  ];
    if borde(i,3)==1
      Dirichlet_vertices(k,1) = c; % nodal indice of the Dirichlet vertex
      Dirichlet_vertices(k,2) = borde(i,4); % Dirichlet tag         
      k = k +1;
    elseif borde(i,3)==2
      Neumann_vertices(r,1) = c; % nodal indice of the Neumann vertex
      Neumann_vertices(r,2) = borde(i,4); % Neumann tag           
      r = r + 1;
    end
    c=c+1;
  end

  %SE GENERA MATRIZ DE COORDENADAS DE LOS VERTICES
  %SE DEBE AGREGAR ADEMAS LAS ESQUINAS COMO VERTICES
%   c = Cborde(end,1) +1;
  c = Cborde(end,1);
  CorVer = [];
%   for i=1:length(corners)  %[NewVerticeTag TagNode Coord-Corners]
%     CorVer = [CorVer; c find(Nodos(:,1)==corners(i,1) & Nodos(:,2)==corners(i,2))  corners(i,:) ];
%     c=c+1;    
%   end
  for i=1:size(corners,1)  %[NewVerticeTag TagNode Coord-Corners]
    
    if ~ismember(corners(i,1),corners(1:i-1,1)) % add the corner node only if it has not been previously added
      c=c+1;        
      CorVer = [CorVer; c find(Nodos(:,1)==Nodos(corners(i,1),1) & Nodos(:,2)==Nodos(corners(i,1),2))  Nodos(corners(i,1),:) ];
    end      
    
    if corners(i,2)==1
      Dirichlet_vertices(k,1) = c; % nodal indice of the Dirichlet vertex
      Dirichlet_vertices(k,2) = corners(i,3); % Dirichlet tag     
      k = k +1;
    elseif corners(i,2)==2
      Neumann_vertices(r,1) = c; % nodal indice of the Neumann vertex
      Neumann_vertices(r,2) = corners(i,3); % Neumann tag      
      r = r + 1;     
    end 
    
  end  

  %vertices = centros de celdas + vertices nuevos en el borde + nodos del borde
  Vertices = [centro;  Cborde(:,4:5); CorVer(:,3:4) ];

  %REPARTE CADA NUEVO VERTICE A LA CONECTIVIDAD DE LOS POLIGONOS
  for i=1:size(borde,1)       
   V{Cborde(i,2),1} = [V{Cborde(i,2),1} Cborde(i,1)  ];
   V{Cborde(i,3),1} = [V{Cborde(i,3),1} Cborde(i,1)  ];           
  end

  %REPARTE CADA NODO DE LAS ESQUINAS A LA CONECTIVIDAD DE LOS POLIGONOS
  for i=1:length(CorVer(:,1))      
    V{CorVer(i,2),1} = [V{CorVer(i,2),1} CorVer(i,1)  ];          
  end


  %ORDENA LOS TAG DE LOS VERTICES EN SENTIDO ANTI-HORARIO
  for i=1:length(V)
    Ver = V{i,1}; %vertices del poligono
    Poly = Vertices(Ver,:); %coordenadas de los vertices del poligono
    centro = [mean(Poly(:,1)) mean(Poly(:,2))  ];
    ClockV = Poly - centro;
    angle = atan2(ClockV(:,2),ClockV(:,1));  %calulo de angulo con atan2
    tags = [ Ver' angle];    
    tags = sortrows(tags,2); %ordena segun angulo
    V{i,1} = tags(:,1);
  end
  
  
%   Nodos = CorrectNodes(Nodos, Vertices, corners, Cborde(:,1:3));
  
  
  
%   % El indice de los nodos no coindice con el indice de las celdas, se cambia
%   % entonces el orden de la matriz de conectividad para que coincida
%   figure;
%   V2 = cell(NumberNodes,1);
%   for i=1:length(Nodos)
%     in = inpolygon(Nodos(:,1),Nodos(:,2),Vertices(V{i,1},1),Vertices(V{i,1},2));    
%     
% %     i
% %     in
% %     Vertices(V{i,1},1)
% %     Vertices(V{i,1},2)
%     
%     if i==3
%     pgon = polyshape(Vertices(V{i,1},1),Vertices(V{i,1},2));
%     plot(pgon)
%     hold on;
%     plot(Nodos(i,1),Nodos(i,2),'*')    
%     end
%     
%     
%     V2{in,1} =  V{i,1};
%   end
%   V=V2;

end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %              Created by  : Rodrigo Silva                                     
  %              E-mail      : rsilvav@outlook.com                                        
  %              Version     : 1.0                         
  %              Date        : June 29, 2018  
  %
  %           Department of Mechanical Engineering                        
  %           University of Chile, Santiago, CHILE    
  %
  % PURPOSE: corrects the nodes that are outside the domain,
  %          after transporfar mesh T3 to poligonos
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
%ENTRADAS:  
% Nodos
% x1 y1
% x2 y2
% x3 y3...
% 
% Vertices
% x1 y1
% x2 y2
% x3 y3...
% 
%esquinas del dominio
%corners = [c1 c1; c2 c2; ...]
%
%conectividad del borde de la malla de poligonos borde0
%[TagCelda   v1 v2]
%[TagCelda   v1 v2]
%[TagCelda   v1 v2]...

%SALIDA:
% Nodos
% x1 y1
% x2 y2
% x3 y3...

function Nodos = CorrectNodes(Nodos, Vertices, corners, borde0)

  Cborde=zeros(length(borde0(:,1)),5);
  for i=1:length(borde0(:,1))

  % [TagNode Ver1 Ver2 NewNodx NewNody]  
    Cborde(i,:) = [ borde0(i,1) borde0(i,2) borde0(i,3) (Vertices(borde0(i,2),:)+Vertices(borde0(i,3),:))/2  ];

  end
  %verifica nuevos nodos del borde creados (con nodos falsos cerca de las esquinas)
  % plot(Cborde(:,4), Cborde(:,5), 'ko')
  %  hold on


  %Busca los tags de los nodos de las esquinas
  TagNodeCorners = [];
  for i=1:length(corners(:,1))
   TagNodeCorners = [ TagNodeCorners; find(Nodos(:,1)==corners(i,1) & Nodos(:,2)==corners(i,2))];  

  end

  %elimina nodos generados cerca de las esquinas.
  for i=1:length(TagNodeCorners)
     aa = find(Cborde(:,1)==TagNodeCorners(i) );
     Cborde(aa,:)=[];   
  end
  %verifica nuevos nodos del borde creados
  % plot(Cborde(:,4), Cborde(:,5), 'ko')
  %  hold on


   %reemplaza nodos del borde
  for i=1:length(Cborde(:,1))  
  Nodos(Cborde(i,1),:) = [Cborde(i,4) Cborde(i,5)];  %Cborde(i,4:5);   
  end
  %verifica la totalidad de los nuevos nodos
  % plot(Nodos(:,1), Nodos(:,2), '*b')
  % hold on

end
  
  
  

 