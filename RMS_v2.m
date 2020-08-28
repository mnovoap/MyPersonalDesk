clear
clc
tic
% muestra las variables
% ncdisp('profile.nc')
% lee los datos usando una funcion los crea tipo cell {}, saca los
% perfiles para cada punto
Punto=11; % SE DEBE ESPECIFICAR EL PUNTO

%Se da 1 a los meses a analizar, y 0 a los que no

ene = 0;
feb = 1;
ago = 0;
nov = 0;


% las columnas representan las variables y las filas el archivo(long. malla)
[ Perfil_elcom_1, tiempo, z] = DatosPerfil_elcom(Punto); 
%% 
%Temperatura
tempe={Perfil_elcom_1{:,1}};  %Extrae variable de todas la primera columna, y la separa en puntos
for i=1:size(Perfil_elcom_1,1) 
    temp=[z{i} tempe{i}]; %Cambia los datos NAN por cero
    temp(isnan(temp))=0;  %Cambia los datos NAN por cero
    for k=1:size(temp,2)
        temp(temp(:,k)==0,:)=[]; %elimina las filas de ceros
    end
    % crea una matriz por cada malla
    v=genvarname(['T_Malla',num2str(i)]); %Y guarda la variable en T_malla_1, T_malla_2... 
    eval([v, '=temp']); %Evalua
    clc
end
% Salinidad
sali={Perfil_elcom_1{:,2}}; %Extrae variable de todas la primera columna, y la separa en puntos
for i=1:size(Perfil_elcom_1,1) % Por cada punto
    sal=[z{i} sali{i}]; %Cambia los datos NAN por cero
    sal(isnan(sal))=0;  %Cambia los datos NAN por cero
    for k=1:size(sal,2)
        sal(sal(:,k)==0,:)=[]; %elimina las filas de ceros
    end
    % crea una matriz por cada malla
    v=genvarname(['S_Malla',num2str(i)]); %Y guarda la variable en S_malla_1, S_malla_2... 
    eval([v, '=sal']); %Evalua
    clc
end

% TIEMPO
for i=1:6
    fecha_1(:,i)=[tiempo{1,i}]; % Tiempo para el Archivo 1
    % fecha_2(:,i)=[tiempo{2,i}]; % Tiempo para el Archivo 2
end

tiempo_1=datenum(fecha_1(:,1:6));
% tiempo_2=datenum(fecha_2(:,1:6));

%FECHAS

    tiempo_m = datenum(fecha_1(:,1:6));
    fecha_e  = (matlab2CWRdate(tiempo_m))';
    
    % agosto
    % 
    % noviembre
    % 
    % enero
    %
    % febrero
    %
toc
%% CARGAR CAMPAÑAS 

% Enero
if ene == 1
    load('J:\1.adriana_puello\Desktop\Mateo\Envolventes\Perfiles_CTD\Ene_corr.mat');
    Sal_Punto_Ene_3 = Ene_p3_2(1:enemax(1),[2,6]);
    Temp_Punto_Ene_3 = Ene_p3_2(1:enemax(1),[2,3]);
    Sal_Punto_Ene_5 = Ene_p5_2(1:enemax(2),[2,6]);
    Temp_Punto_Ene_5 = Ene_p5_2(1:enemax(2),[2,3]);
    Sal_Punto_Ene_6 = Ene_p6_2(1:enemax(3),[2,6]);
    Temp_Punto_Ene_6 = Ene_p6_2(1:enemax(3),[2,3]);
    Sal_Punto_Ene_7 = Ene_p7_r_2(1:enemax(4),[2,6]);
    Temp_Punto_Ene_7 = Ene_p7_r_2(1:enemax(4),[2,3]);
    Sal_Punto_Ene_8 = Ene_p8_2(1:enemax(5),[2,6]);
    Temp_Punto_Ene_8 = Ene_p8_2(1:enemax(5),[2,3]);
    Sal_Punto_Ene_9 = Ene_p9_2_2(1:enemax(7),[2,6]);
    Temp_Punto_Ene_9 = Ene_p9_2_2(1:enemax(7),[2,3]);
    Sal_Punto_Ene_10 = Ene_p10_2(1:enemax(8),[2,6]);
    Temp_Punto_Ene_10 = Ene_p10_2(1:enemax(8),[2,3]);
    [a]=find(fecha_e>=2018026.0 & fecha_e<=2018026.005); %[b]=find(fecha_e>=2018027.0 & fecha_e<=2018027.005);
end
% Febrero
if feb == 1
    load('J:\1.adriana_puello\Desktop\Mateo\Envolventes\Perfiles_CTD\Feb_corr.mat');

    Sal_Punto_Feb_2 = Feb_p2_2(1:febmax(1),[2,6]);
    Temp_Punto_Feb_2 = Feb_p2_2(1:febmax(1),[2,3]);
    Sal_Punto_Feb_3 = Feb_p3_2(1:febmax(2),[2,6]);
    Temp_Punto_Feb_3 = Feb_p3_2(1:febmax(2),[2,3]);
    Sal_Punto_Feb_4 = Feb_p4_2(1:febmax(3),[2,6]);
    Temp_Punto_Feb_4 = Feb_p4_2(1:febmax(3),[2,3]);
    Sal_Punto_Feb_5 = Feb_p5_2(1:febmax(4),[2,6]);
    Temp_Punto_Feb_5 = Feb_p5_2(1:febmax(4),[2,3]);
    Sal_Punto_Feb_6 = Feb_p6_2(1:febmax(5),[2,6]);
    Temp_Punto_Feb_6 = Feb_p6_2(1:febmax(5),[2,3]);
    Sal_Punto_Feb_7 = Feb_p7_2(1:febmax(6),[2,6]);
    Temp_Punto_Feb_7 = Feb_p7_2(1:febmax(6),[2,3]);
    Sal_Punto_Feb_8 = Feb_p8_2(1:febmax(7),[2,6]);
    Temp_Punto_Feb_8 = Feb_p8_2(1:febmax(7),[2,3]);
    Sal_Punto_Feb_9 = Feb_p9_2(1:febmax(8),[2,6]);
    Temp_Punto_Feb_9 = Feb_p9_2(1:febmax(8),[2,3]);
    Sal_Punto_Feb_10 = Feb_p10_2(1:febmax(9),[2,6]);
    Temp_Punto_Feb_10 = Feb_p10_2(1:febmax(9),[2,3]);
    Sal_Punto_Feb_E = Feb_pE_2(1:febmax(10),[2,6]);
    Temp_Punto_Feb_E = Feb_pE_2(1:febmax(10),[2,3]);
    
    [a]=find(fecha_e>=2018035.0 & fecha_e<=2018035.01); %[b]=find(fecha_e>=2018036.0 & fecha_e<=2018036.005);

end
% Agosto
if ago == 1
    load('C:\Users\USER\Desktop\Semestre X\TDG\Estadisticos Finales\Perfiles_CTD\Ago_corr.mat');

    Sal_Punto_Ago_2 = Ago_p2_2(1:agomax(1),[2,6]);
    Temp_Punto_Ago_2 = Ago_p2_2(1:agomax(1),[2,3]);
    Sal_Punto_Ago_3 = Ago_p3_2(1:agomax(2),[2,6]);
    Temp_Punto_Ago_3 = Ago_p3_2(1:agomax(2),[2,3]);
    Sal_Punto_Ago_4 = Ago_p4_2(1:agomax(3),[2,6]);
    Temp_Punto_Ago_4 = Ago_p4_2(1:agomax(3),[2,3]);
    Sal_Punto_Ago_5 = Ago_p5_2(1:agomax(4),[2,6]);
    Temp_Punto_Ago_5 = Ago_p5_2(1:agomax(4),[2,3]);
    Sal_Punto_Ago_7 = Ago_p7_2(1:agomax(5),[2,6]);
    Temp_Punto_Ago_7 = Ago_p7_2(1:agomax(5),[2,3]);
    Sal_Punto_Ago_8 = Ago_p8_2(1:agomax(6),[2,6]);
    Temp_Punto_Ago_8 = Ago_p8_2(1:agomax(6),[2,3]);
    Sal_Punto_Ago_9 = Ago_p9_2(1:agomax(7),[2,6]);
    Temp_Punto_Ago_9 = Ago_p9_2(1:agomax(7),[2,3]);
    Sal_Punto_Ago_10 = Ago_p10_2(1:agomax(8),[2,6]);
    Temp_Punto_Ago_10 = Ago_p10_2(1:agomax(8),[2,3]);
    
    [a]=find(fecha_e>=2017240.0 & fecha_e<=2017240.005); % [b]=find(fecha_e>=2017242.0 & fecha_e<=2017242.005);
end
% Noviembre
if nov == 1
    load('C:\Users\USER\Desktop\Semestre X\TDG\Estadisticos Finales\Perfiles_CTD\Nov_corr.mat');

    Sal_Punto_Nov_2 = Nov_p2_2(1:novmax(1),[2,6]);
    Temp_Punto_Nov_2 = Nov_p2_2(1:novmax(1),[2,3]);
    Sal_Punto_Nov_3 = Nov_p3_2(1:novmax(2),[2,6]);
    Temp_Punto_Nov_3 = Nov_p3_2(1:novmax(2),[2,3]);
    Sal_Punto_Nov_4 = Nov_p4_2(1:novmax(3),[2,6]);
    Temp_Punto_Nov_4 = Nov_p4_2(1:novmax(3),[2,3]);
    Sal_Punto_Nov_5 = Nov_p5_2(1:novmax(4),[2,6]);
    Temp_Punto_Nov_5 = Nov_p5_2(1:novmax(4),[2,3]);
    Sal_Punto_Nov_6 = Nov_p6_2(1:novmax(5),[2,6]);
    Temp_Punto_Nov_6 = Nov_p6_2(1:novmax(5),[2,3]);
    Sal_Punto_Nov_7 = Nov_p7_2(1:novmax(6),[2,6]);
    Temp_Punto_Nov_7 = Nov_p7_2(1:novmax(6),[2,3]);
    Sal_Punto_Nov_8 = Nov_p8_2(1:novmax(7),[2,6]);
    Temp_Punto_Nov_8 = Nov_p8_2(1:novmax(7),[2,3]);
    Sal_Punto_Nov_9 = Nov_p9_2(1:novmax(8),[2,6]);
    Temp_Punto_Nov_9 = Nov_p9_2(1:novmax(8),[2,3]);
    Sal_Punto_Nov_10 = Nov_p10_1_2(1:novmax(9),[2,6]);
    Temp_Punto_Nov_10 = Nov_p10_1_2(1:novmax(9),[2,3]);
    [a]=find(fecha_e>=2017333.0 & fecha_e<=2017332.005); %[b]=find(fecha_e>=2017334.0 & fecha_e<=2017334.005);
end
%% ESTADISTICOS DE RENDIMIENTO
% La malla patrón es la de 100 metros

BE_GEN = []; %Error porcentual porcentaje de error
NASH1_GEN = []; %NASH 1 eficiencia
NASH2_GEN = []; %NASH 2 eficiencia
EAM_GEN = []; %ERROR CUADRATICO MEDIO medida del error
MAPE_GEN = []; %MAPE Desviacion promedio de errores en procentaje
RMSE_GEN = []; %ERROR medio, mide el error
R2_GEN =  []; %R2 general
SPIER_GEN = [];  %Spierman, coeficiente de correlación para independencia

for i = 1:11 %Revisión para cada punto de cada estacion
    
    if ene == 1 %Eleccion de mes
        if i == 3  % Punto 3
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla3(:,1);
            k = 0;
            for j = 1:length(T_Malla3(:,1))
                aux_2 = find(abs(Sal_Punto_Ene_3(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Ene_3(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla3(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Ene_3(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Ene_3(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla3(j,a);
                    
                end
            end
            
            
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) );
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
        end
        if i == 5
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla5(:,1);
            k = 0;
            for j = 1:length(T_Malla5(:,1))
                aux_2 = find(abs(Sal_Punto_Ene_5(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Ene_5(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla5(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Ene_5(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Ene_5(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla5(j,a);
                    
                end
            end
            
            
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) );
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
        end
        if  i == 6
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla6(:,1);
            k = 0;
            for j = 1:length(T_Malla6(:,1))
                aux_2 = find(abs(Sal_Punto_Ene_6(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Ene_6(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla6(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Ene_6(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Ene_6(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla6(j,a);
                    
                end
            end
            
            
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) );
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
        end
        if i == 7
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla7(:,1);
            k = 0;
            for j = 1:length(T_Malla7(:,1))
                aux_2 = find(abs(Sal_Punto_Ene_7(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Ene_7(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla7(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Ene_7(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Ene_7(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla7(j,a);
                    
                end
            end
            
            
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) );
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
        end
        if i == 8
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla8(:,1);
            k = 0;
            for j = 1:length(T_Malla8(:,1))
                aux_2 = find(abs(Sal_Punto_Ene_8(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Ene_8(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla8(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Ene_8(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Ene_8(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla8(j,a);
                    
                end
            end
            
          
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) );
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
        end
        if i == 9
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla9(:,1);
            k = 0;
            for j = 1:length(T_Malla9(:,1))
                aux_2 = find(abs(Sal_Punto_Ene_9(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Ene_9(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla9(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Ene_9(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Ene_9(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla9(j,a);
                    
                end
            end
            
           
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) ); 
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
        end
        if i == 10
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla10(:,1);
            k = 0;
            for j = 1:length(T_Malla10(:,1))
                aux_2 = find(abs(Sal_Punto_Ene_10(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Ene_10(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla10(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Ene_10(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Ene_10(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla10(j,a);
                    
                end
            end
            
        
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) );
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
        end
    elseif feb == 1
		if i == 2
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla2(:,1);
            k = 0;
            for j = 1:length(T_Malla2(:,1))
                aux_2 = find(abs(Sal_Punto_Feb_2(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Feb_2(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla2(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Feb_2(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Feb_2(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla2(j,a);
                    
                end
            end
            
            
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) );
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
        end
        if i == 3  % Punto 3
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla3(:,1);
            k = 0;
            for j = 1:length(T_Malla3(:,1))
                aux_2 = find(abs(Sal_Punto_Feb_3(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Feb_3(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla3(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Feb_3(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Feb_3(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla3(j,a);
                    
                end
            end
            
            
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) );
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
            end
        if i == 4
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla4(:,1);
            k = 0;
            for j = 1:length(T_Malla4(:,1))
                aux_2 = find(abs(Sal_Punto_Feb_4(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Feb_4(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla4(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Feb_4(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Feb_4(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla4(j,a);
                    
                end
            end
            
            
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) );
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
           SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
        end
        if i == 5
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla5(:,1);
            k = 0;
            for j = 1:length(T_Malla5(:,1))
                aux_2 = find(abs(Sal_Punto_Feb_5(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Feb_5(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla5(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Feb_5(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Feb_5(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla5(j,a);
                    
                end
            end
            
           
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) );
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
        end
        if  i == 6
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla6(:,1);
            k = 0;
            for j = 1:length(T_Malla6(:,1))
                aux_2 = find(abs(Sal_Punto_Feb_6(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Feb_6(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla6(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Feb_6(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Feb_6(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla6(j,a);
                    
                end
            end
            
           
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) );
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
            
        end
        if i == 7
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla7(:,1);
            k = 0;
            for j = 1:length(T_Malla7(:,1))
                aux_2 = find(abs(Sal_Punto_Feb_7(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Feb_7(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla7(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Feb_7(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Feb_7(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla7(j,a);
                    
                end
            end
            
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) ); 
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
            
        end
        if i == 8
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla8(:,1);
            k = 0;
            for j = 1:length(T_Malla8(:,1))
                aux_2 = find(abs(Sal_Punto_Feb_8(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Feb_8(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla8(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Feb_8(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Feb_8(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla8(j,a);
                    
                end
            end
            
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) );
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
        end
        if i == 9
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla9(:,1);
            k = 0;
            for j = 1:length(T_Malla9(:,1))
                aux_2 = find(abs(Sal_Punto_Feb_9(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Feb_9(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla9(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Feb_9(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Feb_9(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla9(j,a);
                    
                end
            end
            
            
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) ); 
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
            
        end
        if i == 10
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla10(:,1);
            k = 0;
            for j = 1:length(T_Malla10(:,1))
                aux_2 = find(abs(Sal_Punto_Feb_10(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Feb_10(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla10(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Feb_10(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Feb_10(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla10(j,a);
                    
                end
            end
            
           
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) );
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
            
        end
		if i == 11
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla11(:,1);
            k = 0;
            for j = 1:length(T_Malla11(:,1))
                aux_2 = find(abs(Sal_Punto_Feb_E(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Feb_E(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla11(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Feb_E(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Feb_E(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla11(j,a);
                    
                end
            end
            
            
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) );
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
            
        end
    elseif ago == 1
		if i == 2
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla2(:,1);
            k = 0;
            for j = 1:length(T_Malla2(:,1))
                aux_2 = find(abs(Sal_Punto_Ago_2(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Ago_2(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla2(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Ago_2(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Ago_2(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla2(j,a);
                    
                end
            end
            
            
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) );
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
        end
        if i == 3  % Punto 3
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla3(:,1);
            k = 0;
            for j = 1:length(T_Malla3(:,1))
                aux_2 = find(abs(Sal_Punto_Ago_3(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Ago_3(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla3(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Ago_3(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Ago_3(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla3(j,a);
                    
                end
            end
            
            
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) );
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
            
            end
        if i == 4
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla4(:,1);
            k = 0;
            for j = 1:length(T_Malla4(:,1))
                aux_2 = find(abs(Sal_Punto_Ago_4(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Ago_4(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla4(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Ago_4(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Ago_4(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla4(j,a);
                    
                end
            end
            
            
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) );
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
            
        end
        if i == 5
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla5(:,1);
            k = 0;
            for j = 1:length(T_Malla5(:,1))
                aux_2 = find(abs(Sal_Punto_Ago_5(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Ago_5(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla5(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Ago_5(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Ago_5(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla5(j,a);
                    
                end
            end
            
            
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) );
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
            
        end
        if i == 7
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla7(:,1);
            k = 0;
            for j = 1:length(T_Malla7(:,1))
                aux_2 = find(abs(Sal_Punto_Ago_7(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Ago_7(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla7(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Ago_7(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Ago_7(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla7(j,a);
                    
                end
            end
            
            
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) ); 
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
        end
        if i == 8
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla8(:,1);
            k = 0;
            for j = 1:length(T_Malla8(:,1))
                aux_2 = find(abs(Sal_Punto_Ago_8(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Ago_8(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla8(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Ago_8(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Ago_8(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla8(j,a);
                    
                end
            end
            
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) );
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
        end
        if i == 9
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla9(:,1);
            k = 0;
            for j = 1:length(T_Malla9(:,1))
                aux_2 = find(abs(Sal_Punto_Ago_9(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Ago_9(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla9(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Ago_9(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Ago_9(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla9(j,a);
                    
                end
            end
            
            
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) ); 
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
        end
        if i == 10
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla10(:,1);
            k = 0;
            for j = 1:length(T_Malla10(:,1))
                aux_2 = find(abs(Sal_Punto_Ago_10(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Ago_10(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla10(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Ago_10(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Ago_10(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla10(j,a);
                    
                end
            end
            
            
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) );
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SSPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
        end    
    elseif nov == 1 
		if i == 2
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla2(:,1);
            k = 0;
            for j = 1:length(T_Malla2(:,1))
                aux_2 = find(abs(Sal_Punto_Nov_2(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Nov_2(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla2(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Nov_2(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Nov_2(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla2(j,a);
                    
                end
            end
            
            
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) );
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
        end
        if i == 3  % Punto 3
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla3(:,1);
            k = 0;
            for j = 1:length(T_Malla3(:,1))
                aux_2 = find(abs(Sal_Punto_Nov_3(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Nov_3(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla3(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Nov_3(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Nov_3(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla3(j,a);
                    
                end
            end
            
            
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) );
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
            
            end
        if i == 4
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla4(:,1);
            k = 0;
            for j = 1:length(T_Malla4(:,1))
                aux_2 = find(abs(Sal_Punto_Nov_4(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Nov_4(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla4(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Nov_4(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Nov_4(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla4(j,a);
                    
                end
            end
            
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) );
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
            
        end
        if i == 5
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla5(:,1);
            k = 0;
            for j = 1:length(T_Malla5(:,1))
                aux_2 = find(abs(Sal_Punto_Nov_5(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Nov_5(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla5(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Nov_5(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Nov_5(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla5(j,a);
                    
                end
            end
            
            
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) );
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
            
        end
        if i == 6
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla6(:,1);
            k = 0;
            for j = 1:length(T_Malla6(:,1))
                aux_2 = find(abs(Sal_Punto_Nov_6(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Nov_6(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla6(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Nov_6(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Nov_6(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla6(j,a);
                    
                end
            end
            
            
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) );
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
        end
        if i == 7
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla7(:,1);
            k = 0;
            for j = 1:length(T_Malla7(:,1))
                aux_2 = find(abs(Sal_Punto_Nov_7(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Nov_7(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla7(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Nov_7(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Nov_7(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla7(j,a);
                    
                end
            end
            
            
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) );
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
        end
        if i == 8
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla8(:,1);
            k = 0;
            for j = 1:length(T_Malla8(:,1))
                aux_2 = find(abs(Sal_Punto_Nov_8(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Nov_8(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla8(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Nov_8(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Nov_8(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla8(j,a);
                    
                end
            end
            
            
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) );
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
        end
        if i == 9
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla9(:,1);
            k = 0;
            for j = 1:length(T_Malla9(:,1))
                aux_2 = find(abs(Sal_Punto_Nov_9(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Nov_9(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla9(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Nov_9(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Nov_9(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla9(j,a);
                    
                end
            end
            
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) ); 
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
        end
        if i == 10
            OBS = []; %Serie observada por campo
            SIM = []; %Serie simulada en el ELCOM
            aux = -T_Malla10(:,1);
            k = 0;
            for j = 1:length(T_Malla10(:,1))
                aux_2 = find(abs(Sal_Punto_Nov_10(:,1)-aux(j)) <= 0.5);
                if length(aux_2)>0
                    k = 1 + k;
                    aux_2 = aux_2(1);
                    aux_2 = Sal_Punto_Nov_10(aux_2,2);
                    
                    OBS(k,1) = aux_2;
                    SIM(k,1) = S_Malla10(j,a);
                    
                    
                end
                aux_3 = find(abs(Temp_Punto_Nov_10(:,1)-aux(j)) <= 0.5);
                if length(aux_3)>0
                    aux_3 = aux_3(1);
                    aux_3 = Temp_Punto_Nov_10(aux_3,2);
                    
                    OBS(k,2) = aux_3;
                    SIM(k,2) = T_Malla10(j,a);
                    
                end
            end
            
            
            BE(1,1) =  100*abs( mean(OBS(:,1)) - mean(SIM(:,1)) )/mean(OBS(:,1)); %Sal en %
            BE(1,2) =  100*abs( mean(OBS(:,2)) - mean(SIM(:,2)) )/mean(OBS(:,2)); %Temp en %
            
            NASH1(1,1) = 100*( 1 -  sum( ( OBS(:,1) - SIM(:,1) ).^2 )/sum( (OBS(:,1)-mean(OBS(:,1))).^2 ) ); %Sal en %
            NASH1(1,2) = 100*( 1 -  sum( ( OBS(:,2) - SIM(:,2) ).^2 )/sum( (OBS(:,2)-mean(OBS(:,2))).^2 ) ); %Temp en %

            NASH2(1,1) = 100*( 1 -  sum( ( sqrt(OBS(:,1)) - sqrt(SIM(:,1)) ).^2 )/sum( (sqrt(OBS(:,1))-mean(sqrt(OBS(:,1)))).^2)); %Sal en %
            NASH2(1,2) = 100*( 1 -  sum( ( sqrt(OBS(:,2)) - sqrt(SIM(:,2)) ).^2 )/sum( (sqrt(OBS(:,2))-mean(sqrt(OBS(:,2)))).^2)); %Temp en %
            
            EAM(1,1) = sum(abs( OBS(:,1) - SIM(:,1) ))/length(OBS(:,1));
            EAM(1,2) = sum(abs( OBS(:,2) - SIM(:,2) ))/length(OBS(:,2));
            
            RMSE(1,1) = sqrt( (sum( (OBS(:,1) - SIM(:,1) ).^2 ))/length(OBS(:,1)) );
            RMSE(1,2) = sqrt( (sum( (OBS(:,2) - SIM(:,2) ).^2 ))/length(OBS(:,2)) );
            
            MAPE(1,1) = 100/(length(OBS(:,1))) * sum( abs (( ( OBS(:,1) - SIM(:,1) ))./OBS(:,1)) );
            MAPE(1,2) = 100/(length(OBS(:,2))) * sum( abs (( ( OBS(:,2) - SIM(:,2) ))./OBS(:,2)) );
            
            n = length(OBS(:,1));
            R2(1,1) = (  n*sum(OBS(:,1).*SIM(:,1)) - sum(OBS(:,1))*sum(SIM(:,1)) )/sqrt( ( n*sum(OBS(:,1).^2) - n*(sum(OBS(:,1))^2) )*( n*sum(SIM(:,1).^2) - n*(sum(SIM(:,1))^2) ));
            R2(1,2) = (  n*sum(OBS(:,2).*SIM(:,2)) - sum(OBS(:,2))*sum(SIM(:,2)) )/sqrt( ( n*sum(OBS(:,2).^2) - n*(sum(OBS(:,2))^2) )*( n*sum(SIM(:,2).^2) - n*(sum(SIM(:,2))^2) ));
            
            n = length(OBS(:,1));
            SPIER(1,1) = 1 - 6*( sum(( OBS(:,1) - SIM(:,1)).^2 ) )/( n*(n*n-1) );
            SPIER(1,2) = 1 - 6*( sum(( OBS(:,2) - SIM(:,2)).^2 ) )/( n*(n*n-1) );
            
            BE_GEN = [BE_GEN;BE]; %Error porcentual porcentaje de error
            NASH1_GEN = [NASH1_GEN;NASH1]; %NASH 1 eficiencia
            NASH2_GEN = [NASH2_GEN;NASH2]; %NASH 2 eficiencia
            EAM_GEN = [EAM_GEN;EAM]; %ERROR CUADRATICO MEDIO medida del error
            MAPE_GEN = [MAPE_GEN;MAPE]; %MAPE Desviacion promedio de errores en procentaje
            RMSE_GEN = [RMSE_GEN;RMSE]; %ERROR medio, mide el error
            R2_GEN =  [R2_GEN;R2]; %R2 general
            SPIER_GEN = [SPIER_GEN;SPIER]; %Coeficiente de correlacion spierman
            
        end
		
    end
      
end


BE_GEN 
NASH1_GEN 
NASH2_GEN 
EAM_GEN
MAPE_GEN 
RMSE_GEN 
toc