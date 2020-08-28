
%% INTERPOLACION MALLA

% ESTE CODIGO SIRVE PARA HACER UNA MALLA VARIABLE A PARTIR DE UNA MALLA
% SEMILLA PROVENIENTE DE ARCGIS
% LA MEJOR FORMA PARA PROGRAMARLO ES TRABAJAR EN FUNCION DE LA DISTANCIA
% TOMANDO COMO BASE EL TAMAÑO DE LA CELDA DE LA MALLA SEMILLA

% AUTORES: Sebastian Gomez, Mateo Novoa Parra

% DATOS INICIALES:

% % Coordenadas Limites -- No se usa finalmente
% % # EN GRADOS
% lat_inf=11.169754775;      lat_sup=11.302637525;
% lon_izq=-74.253675225;     lon_der=-74.191513125;
% % # EN METROS
% Top=1249564.55625869;      Bottom=1234852.21585203;
% Right=588274.018099295;    Left=581449.226388492;

% dist_x = Right-Left; dist_y = Top-Bottom;

% Para igual tamaño de malla 

% x=linspace(Left,Right,size(Raster,2));
% y=linspace(Top,Bottom,size(Raster,1));

% % OPCIONES MALLA
%Opcion 1: Malla orientada en los ejes x y y. (Porcentaje uniforme)

clear
clc
close all

% CARGAR BATIMETRIA ASCII QUE PROVIENE DE ARCGIS
% load('Bati_60_metros.mat');
M = importdata('A1_60X60_Matlab.txt');
                              %0.0004504505  0.0005405405
tam = 60.1741 % 60.17416;  100.2882        

Raster = M * -1; %Adecuar la malla con los valores negativos
Raster(Raster == 9999) = nan; %Indicar tierra
f = size(M,1); co = size(M,2); %Filas columnas


dis_x = co*tam; dis_y = f*tam;    %Distancias totales de la matriz original 
                                  %Se multiplica por 50 porque ese es el tamaño de las celdas originales

x=linspace(0,dis_x,size(Raster,2)); %Vector x coordenadas originales (Todas con resolucion de 50m)
y=linspace(0,dis_y,size(Raster,1)); %Vector y coordenadas originales (Todas con resolucion de 50m)

% CELDAS DEL EMISARIO
% # Estas se deben dejar sin modificar #
i = [144 145 146 147]; % filas --> Y    (208 209) 25 (104 105 106 107) 50 (88 89 90 91) 60 (75 76 77 78) Ag (47 48 49 50) Lina
j = [79 80 81 82]; % columnas --> X  (163 164) 25 (81 82 83 84) 50 (66 67 68 69) 60 (87 88 89 90) Ag (43 44 45 46) Lina

%Posición vector x y y emisario.
xem = j*tam; %Para ubicar de nuevas casillas de emisario
yem = i*tam; 


%PORCENTAJES DE CRECIMIENTO X / Y

Px = 0.00;   %Porcentaje de crecimiento,  x
ii1 = 1;        %0.03, 0.04 , 0.05 *2 dado que sean proporcionales ambas
vector1(ii1,1) = tam;
suma1(ii1) = 0; %Vector de suma con el que se tiene como ayuda en construir
                        %la malla variable
                        
Py = 0.00;  %Porcentaje de crecimiento,  y
ii2 = 1;    %0.03, 0.04 , 0.05
vector2(ii2,1) = tam;
suma2(ii2) = 0; %Vector de suma con el que se tiene como ayuda en construir
                           %la malla variable
                                
%Tamaños de malla según sea x o sea y

while roundn(suma1(ii1),2) < dis_x
    ii1=ii1+1;
    vector1(ii1,1) = roundn(vector1(ii1-1)*(1+Px),-64);
    suma1(ii1) = suma1(ii1-1) + vector1(ii1,1);
end

while roundn(suma2(ii2),2) < dis_y
    ii2=ii2+1;
    vector2(ii2,1) = roundn( vector2(ii2-1)*(1+Py),-64);
    suma2(ii2) = suma2(ii2-1) + vector2(ii2,1);
end

%Variables auxiliares para encontrar las posiciones de las casillas antes y
%depues del emisario.
y1 = (i(1)-1)*tam; %Distancias por arriba del emisario
y2 = (f-i(length(i)))*tam; %Distancias por abajo del emisario
x1 = (j(1)-1)*tam; %Distancias por izquierda del emisario 
x2 = (co-j(length(j)))*tam; %Distancias por derecha del emisario

a=find(roundn(suma1,-2) >= roundn(x1,-2)); b=find(roundn(suma1,-2) >= roundn(x2,-2)); c=find(roundn(suma2,-2) >= roundn(y1,-2)); d=find(roundn(suma2,-2) >= roundn(y2,-2));
   %Variables auxiliares para buscar dode es mayor el vector para poner el
    %a(1)
    
% Genera los vectores x y y, con el tamaño de las celdas
dx = [sort(vector1(1:a(1)-1,1),'descend');repmat(tam,4,1);vector1(1:b(1)-1,1)]; %Crear vector en x con el orden que se requiere
dy = [sort(vector2(1:c(1)-1,1),'descend');repmat(tam,4,1);vector2(1:d(1)-1,1)];

clear x1 x2 y1 y2 a b c d f co ii i j

% Creamos un nuevo vector conviertiendolo en metros

% Estos dos vectores tienen la función de definir las coordenadas X y Y de
% todos los puntos que se requieran

xq(1,1) = dx(1);
for k=2:length(dx)
    xq(k,1) =  xq(k-1,1)+dx(k);
end
yq(1,1) = dy(1);
for k=2:length(dy)
    yq(k,1) =  yq(k-1,1)+dy(k);
end

% Para tamaño de malla variable 
[X Y] = meshgrid(x,y); %Meshgrid es una funcion en la cual realiza una matriz 
                       %Con las coordenadas x y y de todos los puntos de un
                       %raster, que se definen según las dimensiones de x y y
                       %Matriz con x columnas y y filas
                    
                    
[Xq Yq] = meshgrid(xq,yq); %Igual al anterior, ubica unos nuevas coordenadas

RasterInterp = interp2(X,Y,Raster,Xq,Yq); %Funcion que interpola segun las 
                                          %Coordenadas de todos los puntos y con
%RasterInterp = interp2(Xq,Yq,Raster,Xq,Yq); %Si ya se tiene la batimetria,
                                         %solo ubicacion


RasterInterp(isnan(RasterInterp)) = 0; %Ubicar tierra



%Plotear Escala grises
pcolor(xq,yq,RasterInterp);%shading flat / grafica con los valores del raster
yourColorMap = gray(256); %jet o gray
yourColorMap(256, :) = [0.9290, 0.6940, 0.1250]; %Se asigna en este caso al ulimo valor el blanco (formato RGB)
yourColorMap(256, :) = [139,69,19]/255; %Se asigna en este caso al ulimo valor el blanco (formato RGB)
colormap(yourColorMap); %Graficas

%Plotear Escala color
pcolor(xq,yq,RasterInterp);%shading flat / grafica con los valores del raster
yourColorMap = jet(256); %jet o gray
yourColorMap(256, :) = [255,255,255]/255; %Se asigna en este caso al ulimo valor el blanco (formato RGB)
colormap(yourColorMap); %Graficas



%Barra de colores
c = colorbar('eastoutside');
ylabel(c,{'INTERPOLATED';' DEPTH (m)'});
set(gca,'yDir','reverse') %Ubicarlo con la orientación correcta

%Cambio de nuevo para manejar tierra en ELCOM AQUIIIII
RasterInterp(RasterInterp==0) = 9999; %Indicar tierra



%% Posición Muestras

%Cargar datos 
Muestras = load('MUESTRAS_TOTAL.txt'); %1 col numero de muestreo
                                 %2 col O
                                 %3 col N
                                
Oini = -74.26079;    %Coordenadas esquina superior izquierda %-74.30305 Antes -74.25439 Después (-74.25895 A4) (-74.25895 A3)
Nini =  11.3368;    %11.33494 Antes 11.30521 Despues (11.34784 A4) (11.34784 A3)
                     %lat_sup= 11.30864; Lina
                     %lon_izq=-74.25768; Lina  
                     
%Pasar de coordenadas geograficas a distancias x y y 
Muestras = [Muestras(:,1),(Muestras(:,2)-Oini)*111.3222*1000,(Nini-Muestras(:,3))*111.3222*1000]; %Transformacion grados a m

%Buscar posiciones que se cumplan 
xmst = []; %Vector x donde estara las casillas de los puntos de muestreo
ymst = []; %Vector y donde estará las casillas de los puntos de muestreo

for i = (Muestras(1,1):1:max(Muestras(:,1)))
    aux1 = find (Xq(1,:)>= Muestras(i,2)); %Buscar posicion
    if (isempty(aux1) == 1)
        aux1 = find (Xq(1,:)<= Muestras(i,2)); %Buscar posicion
        aux1 = aux1(length(aux1));
    end
    xmst = [xmst; aux1(1)];%Adjuntar posiciones tipo col
    aux2 = find (Yq(:,1)>= Muestras(i,3));
    if (isempty(aux2) == 1)
       aux2 = find (Yq(:,1)<= Muestras(i,3)); %Buscar posicion
       aux2 = aux2(length(aux2));
    end
    ymst = [ymst; aux2(1)];
end

mst = [ymst,xmst]; %Fil, Col    

%% PARA ORGANIZAR LA MALLA EN FORMATO ELCOM
RasterInterp(isnan(RasterInterp))=9999; %Indicar tierra en los datos vacios
malla = RasterInterp;                   %Cambio de variable

for i=1:length(malla(:,1)) %longitud de las filas de la matriz
    for j=1:length(malla(1,:)) %longitud de las columnas de la matriz
        malla_1(:,1) = 8888; %frontera abierta a la primera columna
        
        malla_1(i+1,j+1)=malla(i,j); %todos los datos de la matriz
        
        if malla(1,j)==9999 %frontera abierta o cerrada a la primera fila
            malla_1(1,j+1)=9999; %Si es tierra sigue siendo tierra
        else 
            malla_1(1,j+1)=8888; %Si no frontera
        end
        %frontera abierta o cerrada a la ultima fila
        if malla(length(malla(:,1)),j)==9999 %Evalua la ultima fila en cada col
            malla_1(length(malla(:,1))+2,j+1)=9999; %Dejar tierra (2 porque crea una fila adicional y tiene que mover 1)
        else
            malla_1(length(malla(:,1))+2,j+1)=8888; %Poner frontera
        end
    end
    malla_1(:,j+2)=9999; %frontera de tierra a la ultima columna
end
%Nuevas Puntos emisario y muestreos
mst = mst + 1;

%Tamaño celdas
Tam = size(malla_1);

%Cambiar dx y dy para acomodar al archivo elcom
dx=[dx(1);dx;dx(length(dx))];
dy=[dy(1);dy;dy(length(dy))];

% MALLA
dlmwrite('Malla_variable.dat',malla_1,'delimiter','\t','precision','%.2f')

% MUESTREO, DESCARGAS Y EMISARIO
dlmwrite('MuDeEmUbicacion.dat',mst,'delimiter','\t','precision','%.0f')

% Tamaño celdas
dlmwrite('Celdas.dat',Tam,'delimiter','\t','precision','%.0f')

% dx
dlmwrite('dy.dat',dx,'delimiter','\t','precision','%.4f') %Se cambian los dx por dy, por la definicion del elcom

% dy
dlmwrite('dx.dat',dy,'delimiter','\t','precision','%.4f') %Se cambian los dy por dx, por la definicion del elcom

%%
% Histograma relaciones
[Dx Dy] = meshgrid(dx,dy); %Meshgrid de diferenciales
Relaciones = [];


for i=1:length(Dx(:,1)) %longitud de las filas
    
    for j=1:length(Dy(1,:)) %longitud de las columnas
        Relaciones = [Relaciones; Dy(i,j)/Dx(i,j)];
        
    end
end
%hist(Relaciones)

%%
% Malla SOD
M = importdata('Malla_variable.dat'); 
SOD = zeros(length(M(:,1)),length(M(1,:)));

for i = 1:length(M(:,1)) %Filas
   for j = 1:length(M(1,:)) %Columnas
       if M(i,j) == 9999
           SOD(i,j) = 9999;
       elseif M(i,j) == 8888
           SOD(i,j) = 8888;
       elseif M(i,j) < -100
           SOD(i,j) = -0.2825;
       elseif M(i,j) < -50
           SOD(i,j) = -0.565;
       elseif  M(i,j) < -30
           SOD(i,j) = -0.8475;
       elseif  M(i,j) < 0
           SOD(i,j) = -1.13;
       end    
   end
end
ymst = ymst +1;
xmst = xmst +1;
xmst(4) = xmst(4)+ 1
ymst(3) = ymst(3)+ 1
ymst(4) = ymst(4)+ 1
mst = [ymst,xmst]; %Fil, Col  
for i = 1:length(ymst(:,1))
    if i == 1
        SOD(ymst(i,1),xmst(1,1)) = -2.95;
    end
    if i == 2
        SOD(ymst(i,1),xmst(i,1)) = -1.32;
    end
    if i == 3
        SOD(ymst(i,1),xmst(i,1)) = -1.26;
    end
    if i == 4
        SOD(ymst(i,1),xmst(i,1)) = -0.89;
    end  
end
%Guardar variable
dlmwrite('MATRIZ_SOD.dat',SOD,'delimiter','\t','precision','%.4f') 
