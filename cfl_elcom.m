function [ CFL stepCFL infoAll infoAllColum] = cfl_elcom()
%LEE los archivo de CFL que da elcom en formato cvs
%   
[NombreDa, DireccionDat] = uigetfile('*.csv','Escoja los archivos (.csv) que desea abrir','MultiSelect','on');
if ischar(NombreDa) %se usa en caso que sea un solo archivo
    NombreDa = cellstr(NombreDa);
end
NombreDat = NombreDa';
% Se especifican los delimitadores para la lectura de los archivos
DELIMITER = ',';    %Delimitador entre los datos
HEADERLINES = 2;   %Número de líneas que no son numéricas, y que se borran en la lectura
%Imprime el porcentaje de lectura de datos
h = waitbar(0,'Porcentaje de lectura: 0%');
%Se crean una matriz multimatrices donde se almacenan los datos leidos para
%cada archivo (el código es más eficiente).
DataCFL=cell(length(NombreDat),1);
FechaCFL = cell(length(NombreDat),1);
Time = cell(length(NombreDat),1);
stepCFL = cell(length(NombreDat),1);

for i=1:length(NombreDat) %Por cada archivo CFL
    fileToRead= [DireccionDat NombreDat{i}]; % Concatena la ruta y el nombre del archivo
    % Se importa el archivo con los parametros mencionados
    newData= importdata(fileToRead, DELIMITER, HEADERLINES);
    % Se almacenan los datos en una nueva multimatriz
    DataCFL{i}= newData.data; 
    % Se almacenan las fechas tipo LakeESP otra multimatriz (fecha juliana)
    stepCFL{i}= newData.textdata(3:end,1); %Juliano
    FechaCFL{i}= newData.textdata(3:end,2); %Gregoriano
%     Se convierte a fecha matlabiana la fecha asociada a COLOMBIA - Approximate_Local_Date_&_Time (UTM - 5 horas)
    Time{i} = datenum(FechaCFL{i,1}(:,1), 'yyyy mmm dd HH:MM:SS');
    % Se obtiene los nombres (encabezados) de las variables y las unidades
    % respectivas
    Header1 = newData.textdata{1,1}; %Titulo
    %Imprime el porcentaje de lectura de datos
    porcen = i/length(NombreDat).*100;
    waitbar(i/length(NombreDat),h,['Porcentaje de lectura: ' sprintf('%2.0f',porcen) '%'])
end
%Se cierra la impresión de lectura
close(h)

%Separa los headers asociados con los archivos y las variables
Header2 = strrep(Header1, ' ', '');
k = strfind(Header2,',');
ini=1;
for i=1:length(k)
    %Se almacenan todos los header en distintas columnas de la multimatriz
    %y adicionalmente el orden, el cual corresponderá luego a la columna
    %ubicación
    
    %Titulos
    Title{1,i} = num2str(i-2); %el -2 porque luego se borran las dos primeras columnas
    Title{2,i} = Header2(ini:k(i));
    ini = k(i)+1;
end
%Se concatena finalmente todas las matrices: tiempo y datos
CFL = cell(length(DataCFL),1); %Resultados finales CFL
for i=1:1:length(DataCFL)
    CFL{i}=[Time{i,1} DataCFL{i,1}];
end

%
infoAll = char(Title(:,3:end)); %Información de las variables del CFL
infoAllColum = Title(:,3:end); %Información de las variables pero en col

end

%¿QUE HACE?
% 1. Lee multiples archivos, y los referencia 
% 2. Extrae los datos (CFL (u,v,w,Umed, Vmed, Wmed) , fecha (jul y greg)
% 3. Extrae titulos
% 4. Organiza el CFL
% 5. Organiza los titutlos
% 6. Organiza el paso del tiempo

