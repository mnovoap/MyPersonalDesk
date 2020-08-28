#Librerias
if(!require("readxl")) install.packages("readxl")
if(!require("dplyr")) install.packages("dplyr")
if(!require("xlsx")) install.packages("xlsx")
if(!require("ie2misc")) install.packages("ie2misc")

library(readxl)
library(dplyr)
library(xlsx)
library(ie2misc)
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------##

#Datos
data = read.csv(choose.files())

Fechas <- select( data, Fecha ) %>%  unlist(use.names = 0) %>% as.Date("%m/%d/%Y") %>% format.Date("%Y") %>% as.numeric()

Pmed <- filter(data, 1997 <= as.numeric(format.Date(as.Date(data$Fecha,"%Y"), "%Y")) & 2011 >= as.numeric(format.Date(as.Date(data$Fecha,"%Y"), "%Y")))  %>% filter(IdParametro == "PRECIPITACION") %>% select(  Valor )  %>%  unlist(use.names = 0) 
Pmed <- read_excel("C:/Users/USER/Downloads/4 TAREA/DATOS/PRECP.xlsx",sheet = "HOJA1") %>% select(  VALOR )  %>%  unlist(use.names = 0) 
Qmed <- filter(data, 1997 <= as.numeric(format.Date(as.Date(data$Fecha,"%Y"), "%Y")) & 2011 >= as.numeric(format.Date(as.Date(data$Fecha,"%Y"), "%Y")))  %>% filter(IdParametro == "CAUDAL") %>% select(  Valor )  %>%  unlist(use.names = 0) 
Qmed <- c(Qmed [1:105],mean(Qmed) ,Qmed[106:179])

View(Pmed)
View(Qmed)

len <- length(Qmed)*0.7

Qmedpron <- Qmed[(len+2):length(Qmed)]
Pmedpron <- Pmed[(len+2):length(Pmed)]

Qmed <- Qmed [1:len]
Pmed <- Pmed[1:len]


##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------##
#Indices

#Seleccionarlos Excel

#MEI <- read_excel(choose.files())
MEI <- read_excel("C:/Users/USER/Downloads/4 TAREA/DATOS/MEI.xlsx",sheet = "Hoja7")
View(MEI)
NAO <- read_excel("C:/Users/USER/Downloads/4 TAREA/DATOS/NAO.xlsx")
View(NAO)
PDO <- read_excel("C:/Users/USER/Downloads/4 TAREA/DATOS/PDO.xlsx")
View(PDO)
SOI <- read_excel("C:/Users/USER/Downloads/4 TAREA/DATOS/SOI.xlsx")
View(SOI)
SST <- read_excel("C:/Users/USER/Downloads/4 TAREA/DATOS/SST.xlsx")
View(SST)

#Filtrar datos
MEI <- filter(MEI, 1997 <= as.numeric(MEI$AÑO) & 2011 >= as.numeric(MEI$AÑO)) %>% select(  VALOR )  %>%  unlist(use.names = 0) 
NAO <- filter(NAO, 1997 <= as.numeric(NAO$AÑO) & 2011 >= as.numeric(NAO$AÑO)) %>% select(  VALOR )  %>%  unlist(use.names = 0) 
PDO <- filter(PDO, 1997 <= as.numeric(PDO$AÑO) & 2011 >= as.numeric(PDO$AÑO)) %>% select(  VALOR )  %>%  unlist(use.names = 0) 
SOI <- filter(SOI, 1997 <= as.numeric(SOI$AÑO) & 2011 >= as.numeric(SOI$AÑO)) %>% select(  VALOR )  %>%  unlist(use.names = 0) 
SST <- filter(SST, 1997 <= as.numeric(SST$AÑO) & 2011 >= as.numeric(SST$AÑO)) %>% select(  VALOR )  %>%  unlist(use.names = 0) 


MEIPred <- MEI[(len+2):length(MEI)]
NAOPred <- NAO[(len+2):length(NAO)]
PDOPred <- PDO[(len+2):length(PDO)]
SOIPred <- SOI[(len+2):length(SOI)]
SSTPred <- SST[(len+2):length(SST)]

MEI <- MEI [1:len]
NAO <- NAO[1:len]
PDO <- PDO [1:len]
SOI <- SOI[1:len]
SST <- SST [1:len]

Ypron3 <- read_excel("C:/Users/USER/Downloads/4 TAREA/DATOS/DatosRedNeuronal.xlsx",sheet = "Hoja1")
View(Ypron3)

##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------##

Pmed <- filter(data, 1997 <= as.numeric(format.Date(as.Date(data$Fecha,"%Y"), "%Y")) & 2011 >= as.numeric(format.Date(as.Date(data$Fecha,"%Y"), "%Y")))  %>% filter(IdParametro == "PRECIPITACION") %>% select(  Valor )  %>%  unlist(use.names = 0) 
Pmed <- read_excel("C:/Users/USER/Downloads/4 TAREA/DATOS/PRECP.xlsx",sheet = "HOJA1") %>% select(  VALOR )  %>%  unlist(use.names = 0) 
Qmed <- filter(data, 1997 <= as.numeric(format.Date(as.Date(data$Fecha,"%Y"), "%Y")) & 2011 >= as.numeric(format.Date(as.Date(data$Fecha,"%Y"), "%Y")))  %>% filter(IdParametro == "CAUDAL") %>% select(  Valor )  %>%  unlist(use.names = 0) 
Qmed <- c(Qmed [1:105],mean(Qmed) ,Qmed[106:179])


#MATRIZ DE COVARIANZAS 
Qmed_inter = (Qmed - mean (Qmed))/(sd(Qmed))
Pmed_inter = (Pmed - mean (Pmed))/(sd(Pmed))

m = 12 #Rezagos
N = length(Pmed) #Longitud de datos

c_wal_Q = c() #Rezagos
c_wal_P = c() #Rezagos

for (j in 0:(m-1)){
  aux1 = c() #Variables aux
  aux2 = c() 
  
  for (i in 1:(N-j)){
    aux1 = c(aux1, Qmed_inter[i]*Qmed_inter[i+j]) 
    aux2 = c(aux2, Pmed_inter[i]*Pmed_inter[i+j])
  }
  
  c_wal_Q = c(c_wal_Q, 1/(N-j*0)*sum(aux1)) 
  c_wal_P = c(c_wal_P, 1/(N-j*0)*sum(aux2))
}


#A <- matrix(aux,nrow=5,ncol=5) #Matriz nxn
A <- toeplitz (c_wal_Q)
B <- toeplitz (c_wal_P)

#DESCOMPOSICION ESPECTRAL
evA <- eigen(A)
valoresA <- evA$values #n valores propios
vectoresA <- evA$vectors #n vectores propios

evB <- eigen(B)
valoresB <- evB$values #n valores propios
vectoresB <- evB$vectors #n vectores propios

#PORCENTAJE VARIANZA
VarA = valoresA/sum(valoresA)
VarB = valoresB/sum(valoresB)
plot (VarA, main = "PORCENTAJE DE VARIANZA", xlab = "EOF #" , ylab= "Fraccion varianza")
lines(VarA)
plot (VarB, main = "PORCENTAJE DE VARIANZA", xlab = "EOF #" , ylab= "Fraccion varianza")
lines(VarB)

#FILTRAR INFORMACION
aux = c()
por = 0.7
for ( i in 1:length(valoresA)){
  aux = c (valoresA[1:i])
  if (sum(aux)/sum(valoresA)>=por){
    break
  }
}
k_Q = length(aux)
m1 = k_Q

aux = c()
por = 0.7
for ( i in 1:length(valoresB)){
  aux = c (valoresB[1:i])
  if (sum(aux)/sum(valoresB)>=por){
    break
  }
}
  
k_P = length(aux)
m2 = k_P

#CP (PROYECCION DE LOS DATOS)
a = matrix ( c( rep( NA, k_Q*(N-m1+1) ) ), nrow=N-m1+1, ncol=k_Q ) #Matriz con cada columna de cada componente proyectado
b = matrix ( c( rep( NA, k_P*(N-m2+1) ) ), nrow=N-m2+1, ncol=k_P ) #Matriz con cada columna de cada componente proyectado

for (v  in 1:k_Q){ #Controla cada componente
  
  for (i in 0:(N-m1)){ #Controla la serie de tiempo
    
    aux  = c()
    
    for (j in 1:m1){ #Control el producto punto, la sumatoria
      aux = c( aux, Qmed_inter[j+i]*vectoresA[j,v])
    }
    aux = sum(aux)
    a[i+1,v]  = aux
  }
}


plot(a[,1], main = "EOF Q 1", ylab  ="", xlab = "REALIZACIÓN EN EL TIEMPO")
lines(a[,1])

plot(a[,2], main = "EOF Q 2", ylab  ="", xlab = "REALIZACIÓN EN EL TIEMPO")
lines(a[,2])

plot(a[,3], main = "EOF Q 3", ylab  ="", xlab = "REALIZACIÓN EN EL TIEMPO")
lines(a[,3])

plot(a[,4], main = "EOF Q 4", ylab  ="", xlab = "REALIZACIÓN EN EL TIEMPO")
lines(a[,4])


for (v  in 1:k_P){ #Controla cada componente
  
  for (i in 0:(N-m2)){ #Controla la serie de tiempo
    
    aux  = c()
    
    for (j in 1:m2){ #Control el producto punto, la sumatoria
      aux = c( aux, Pmed_inter[j+i]*vectoresB[j,v])
    }
    aux = sum(aux)
    b[i+1,v]  = aux
  }
}


plot(b[,1], main = "EOF P 1", ylab  ="", xlab = "REALIZACIÓN EN EL TIEMPO")
lines(b[,1])

plot(b[,2], main = "EOF P 2", ylab  ="", xlab = "REALIZACIÓN EN EL TIEMPO")
lines(b[,2])

plot(b[,3], main = "EOF P 3", ylab  ="", xlab = "MES")
lines(b[,3])

plot(b[,4], main = "EOF P 4", ylab  ="", xlab = "MES")
lines(b[,4])


#CR

RaX_Q  = c() #Reconstruccion 

for (i in 1:N){ #i controla la cantidad de datos
  
  if (1 <= i & i <= (m1-1)){
    
    aux1 = c() #Variable auxialar para sumar
    
    for ( v in 1:k_Q){ #COntrola cada componente
      for (j in 1:i){#Controla rezagos
        aux1 = c(aux1,a[i-j+1,v]*vectoresA[j,v]) #OJOOOO LE PUSE UN +1 al a
      } 
    }
    aux1 = (1/i)*sum(aux1)
    RaX_Q[i] = aux1
    
  } else if(  (m1) <= i & i <= (N - m1 + 1) ){
    aux1 = c() #Variable auxialar para sumar
    
    for ( v in 1:k_Q){ #COntrola cada componente
      for (j in 1:m1){#Controla rezagos
        aux1 = c(aux1,a[ i-j + 1 , v ]*vectoresA[j,v])
      } 
    }
    aux1 = (1/m1)*sum(aux1)
    RaX_Q[i] = aux1
    
  } else if ((N - m1 + 2) <= i & i <= (N) ){
    aux1 = c() #Variable auxialar para sumar
    
    for ( v in 1:k_Q){ #COntrola cada componente
      for (j in (i-N+m1):(m1)){#Controla rezagos
        aux1 = c(aux1,a[i-j + 1,v]*vectoresA[j,v])
      } 
    }
    aux1 = (1/(N-i+1))*sum(aux1)
    RaX_Q[i] = aux1
  }
}

aux2 = RaX_Q*sd(Qmed) + mean (Qmed)
plot(aux2,Qmed, ylim = c(1,5) ,xlim = c(1,5),xlab = "", ylab ="")
lines((1:5),(1:5))

RaX_P  = c() #Reconstruccion 

for (i in 1:N){ #i controla la cantidad de datos
  
  if (1 <= i & i <= (m2-1)){
    
    aux1 = c() #Variable auxialar para sumar
    
    for ( v in 1:k_P){ #COntrola cada componente
      for (j in 1:i){#Controla rezagos
        aux1 = c(aux1,b[i-j+1,v]*vectoresB[j,v])
      } 
    }
    aux1 = (1/i)*sum(aux1)
    RaX_P[i] = aux1
    
  } else if(  (m2) <= i & i <= (N - m2 + 1) ){
    aux1 = c() #Variable auxialar para sumar
    
    for ( v in 1:k_Q){ #COntrola cada componente
      for (j in 1:m2){#Controla rezagos
        aux1 = c(aux1,b[ i-j + 1, v ]*vectoresB[j,v])
      } 
    }
    aux1 = (1/m2)*sum(aux1)
    RaX_P[i] = aux1
    
  } else if ((N - m2 + 2) <= i & i <= (N) ){
    aux1 = c() #Variable auxialar para sumar
    
    for ( v in 1:k_Q){ #COntrola cada componente
      for (j in (i-N+m2):(m2)){#Controla rezagos
        aux1 = c(aux1,b[i-j +1,v]*vectoresB[j,v])
      } 
    }
    aux1 = (1/(N-i+1))*sum(aux1)
    RaX_P[i] = aux1
  }
}


aux3 = RaX_P*sd(Pmed) + mean (Pmed)
plot(aux3,Pmed,xlab = "", ylab ="")
lines((-5:300),(-5:300))


#SERIES RECONSTRUIDAS
plot(aux2, ylim = c(0, max (Qmed)), ylab  = "CAUDAL  (m3/s)", xlab = "MES")
lines(Qmed)
lines(aux2,col = 2,lty = 2)

plot(aux3, ylim = c(0, max (Pmed)), ylab  = "PRECIPITACION  (mm/mes)", xlab = "MES")
lines(Pmed)
lines(aux3,col = 2,lty = 2)

##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------##

#Data frame regresion general

len <- length(Qmed)*0.7

Qmedpron <- Qmed[(len+2):length(Qmed)]
Pmedpron <- Pmed[(len+2):length(Pmed)]

Qmed <- Qmed [1:len]
Pmed <- Pmed[1:len]


data1 <- data.frame(Qmed[2:length(Qmed)],Qmed[1:length(Qmed)-1],MEI[2:length(Qmed)],NAO[2:length(Qmed)],SOI[2:length(Qmed)],PDO[2:length(Qmed)],SST[2:length(Qmed)],Pmed[2:length(Qmed)])
names(data1)<- c("Q","Qrez","MEI","NAO","SOI","PDO","SST","P")
View(data1)


#PRUEBAS

#https://rpubs.com/Joaquin_AR/226291


reg <- function(x, y, ...) {
  points(x,y, ...)
  abline(lm(y~x)) 
}# made to draw regression line instead of lowess line

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  text(0.5, 0.5, txt, cex = 1.1, font = 4)
}


pairs(data1, lower.panel = panel.cor)

cor.test(data1$Q, data1$Qrez, method = "pearson")
cor.test(data1$Q, data1$MEI, method = "pearson")
cor.test(data1$Q, data1$NAO, method = "pearson")
cor.test(data1$Q, data1$SOI, method = "pearson")
cor.test(data1$Q, data1$PDO, method = "pearson")
cor.test(data1$Q, data1$SST, method = "pearson")
cor.test(data1$Q, data1$P, method = "pearson")

modelo <- lm(data1$Q ~ data1$Qrez + data1$MEI +data1$NAO +data1$SOI +data1$PDO +data1$SST + data1$P , data = data1)
summary(modelo)

modelo <- lm(data1$Q ~  data1$Qrez + data1$MEI +data1$SST, data = data1)
summary(modelo)

modelo <- lm(data1$Q ~  data1$Qrez + data1$MEI + data1$P, data = data1)
summary(modelo)

modelo <- lm(data1$Q ~  data1$Qrez + data1$MEI, data = data1)
summary(modelo)
confint(modelo)


##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------##

#Data frame regresion con datos filtrados

data1 <- data.frame(Qmed,aux2[1:length(Qmed)],MEI,NAO,SOI,PDO,SST,aux3[1:length(Qmed)])
names(data1)<- c("Q","Qfil","MEI","NAO","SOI","PDO","SST","Pfil")
View(data1)


#PRUEBAS

#https://rpubs.com/Joaquin_AR/226291

reg <- function(x, y, ...) {
  points(x,y, ...)
  abline(lm(y~x)) 
}# made to draw regression line instead of lowess line

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  text(0.5, 0.5, txt, cex = 1.1, font = 4)
}

pairs(data1,lower.panel = panel.cor)
cor.test(data1$Q, data1$Qfil, method = "pearson")
cor.test(data1$Q, data1$MEI, method = "pearson")
cor.test(data1$Q, data1$NAO, method = "pearson")
cor.test(data1$Q, data1$SOI, method = "pearson")
cor.test(data1$Q, data1$PDO, method = "pearson")
cor.test(data1$Q, data1$SST, method = "pearson")
cor.test(data1$Q, data1$Pfil, method = "pearson")

modelo2 <- lm(data1$Q ~ data1$Qfil + data1$MEI +data1$NAO +data1$SOI +data1$PDO +data1$SST + data1$Pfil , data = data1)
summary(modelo2)

modelo2 <- lm(data1$Q ~  data1$Qfil + data1$MEI +data1$SST, data = data1)
summary(modelo2)

modelo2 <- lm(data1$Q ~  data1$Qfil + data1$MEI + data1$Pfil, data = data1)
summary(modelo2)

modelo2 <- lm(data1$Q ~  data1$Qfil + data1$MEI , data = data1)
summary(modelo2)
confint(modelo2)

##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------##

#ERRORES

#DATOS

Real <- Qmedpron #Reales


QrezPred <- Qmedpron[1:length(Qmedpron)-1] #Qrezagado
QrezPred <-c(Qmed[length(Qmed)],QrezPred)
MEIPred #MEI
QfilPred <- aux2[(len+2):length(aux2)] #Qfil pred
PfilPred <- aux3[(len+2):length(aux3)] #Qfil pred

YPron1 <- modelo$coefficients[1] + modelo$coefficients[2]*QrezPred + modelo$coefficients[3]*MEIPred
YPron2 <- modelo2$coefficients[1] + modelo2$coefficients[2]*QfilPred + modelo2$coefficients[3]*MEIPred

#TEORIA DE RACHAS
#RACHAS

Rachas <- function(x,k){
  #X Vector de datos
  #k veces la desviaciones estandar
  #i es en meses, debe ser positivo
  
  des <- k*sd(x) #Desviaciones estandar por arriba
  arr <- des + mean(x) #Limite superior
  abj <- -des + mean(x) #Limite inferior
  
  ArArr = NULL #areas de rachas arriba
  ArAbj = NULL #areas de rachas abajo
  
  durArr = NULL #Duracion rachas arriba
  durAbj = NULL #DUracion rachas abajo 
  
  n <- length(x) #Total de iteraciones
  
  i = 1 #empieza el ciclo
  while (i <=(n-1)){ 
    #Ciclo for para las n-1 rectas
    #i son los meses
    y1 = x[i]
    y2 = x[i+1]
    
    
    
    aux  = 0 #Variable auxiliar para el aumento de i
    
    if( i == 1 ){ #Criterio distinto para la primera recta
      if ( y1 > arr || y1 < abj ){
        if (abj < y2 & y2 < arr){
          #Crear recta
          m = (y2-y1)/(i+1-i)
          b = y2-m*(i+1)
          
          x2 =  (arr- b)/m #x donde termina racha 
          
          #Area
          A = abs(y1+arr)/2*(x2) - abs(arr)*(x2) #Area 1 por encima 
          
          ArArr = c( A ) #Aqade la racha
          durArr = c( x2 ) #Aqade la duracion de la racha
          
        } 
        else if (arr > y2 & y2 > abj ){
          #Crear recta
          m = (y2-y1)/(i+1-i)
          b = y2-m*(i+1)
          
          x2 =  (abj- b)/m #x donde termina racha 
          
          #Area
          A = abs(y1+abj)/2*(x2) - abs(abj)*(x2) #Area 1 por encima 
          
          ArAbj = c( A ) #Aqade la racha
          durAbj = c( x2 ) #Aqade la duracion de la racha
          
        } 
        else{#Si sigue en la racha
          aux = 1 #aumente por dentro  y no por fuera
          
          if (y2 > arr) { #Por arriba
            # Area acumulativa
            A = c( abs(y2+y1)/2*(i+1-i) - abs(arr)*(i+1-i)) #Area por encima 
            
            while ( y2 > arr){
              i = i+1 #Aumenta el i o el paso
              y1 = x[i]
              y2 = x[i+1]
              
              if (y2 <= arr){
                #Calcular el x2
                #Crear recta
                m = (y2-y1)/(i+1-i)
                b = y2-m*(i+1)
                
                x2 =  (arr- b)/m #x donde inicia racha
                
                #Calcula el area de forma distinta
                A2 = abs(y1+arr)/2*(x2-i) - abs(arr)*(x2-i) #Area 2 por encima 
                A = c(A,A2)
              }else {
                
                # Area acumulativa
                A = c(A, abs(y2+y1)/2*(i+1-i) - abs(arr)*(i+1-i)) #Area por encima 
              }
              
            } 
            ArArr = c( ArArr , sum (A)) #Aqade la racha
            durArr = c( durArr , x2-0 ) #Aqade la duracion de la racha
            
          } else if (y2 < abj){ #Por abajo
            # Area acumulativa
            A = c( abs(y2+y1)/2*(i+1-i) - abs(abj)*(i+1-i)) #Area por encima 
            
            while ( y2 < abj){
              i = i+1 #Aumenta el i o el paso
              y1 = x[i]
              y2 = x[i+1]
              
              if (y2 >= abj){
                #Calcular el x2
                #Crear recta
                m = (y2-y1)/(i+1-i)
                b = y2-m*(i+1)
                
                x2 =  (abj - b)/m #x donde inicia racha
                
                #Calcula el area de forma distinta
                A2 = abs(y1+abj)/2*(x2-i) - abs(abj)*(x2-i) #Area 2 por encima 
                A = c(A,A2)
                
              }else {
                # Area acumulativa
                A = c(A, abs(y2+y1)/2*(i+1-i) - abs(abj)*(i+1-i)) #Area por debajo 
              }
            }
            ArAbj = c( ArAbj , sum (A)) #Aqade la racha
            durAbj = c( durAbj , x2-0) #Aqade la duracion de la racha
          }
          
          
        } 
        
      }else{
        #No haga absolutamente nada, y siga con el ciclo
      }
    }
    else{
      #En caso que la racha no empiece en 1
      
      #Obsservar si empieza una racha
      if ( y2 > arr ){ 
        #Cambia variable auxiliar no permite que aumente
        aux = 1
        
        #Crear recta
        m = (y2-y1)/(i+1-i)
        b = y2-m*(i+1)
        
        x1 =  (arr- b)/m #x donde inicia racha
        
        A = abs(y2+arr)/2*(i+1-x1) - abs(arr)*(i+1-x1) #Area 1 por encima 
        
        while ( y2 > arr){
          i = i+1 #Aumenta el i o el paso
          y1 = x[i]
          y2 = x[i+1]
          
          if (y2 <= arr){
            #Calcular el x2
            #Crear recta
            m = (y2-y1)/(i+1-i)
            b = y2-m*(i+1)
            
            x2 =  (arr- b)/m #x donde inicia racha
            
            #Calcula el area de forma distinta
            A2 = abs(y1+arr)/2*(x2-i) - abs(arr)*(x2-i) #Area 2 por encima 
            A = c(A,A2)
          }else {
            
            # Area acumulativa
            A = c(A, abs(y2+y1)/2*(i+1-i) - abs(arr)*(i+1-i)) #Area por encima 
          }
          
        } 
        ArArr = c( ArArr , sum (A)) #Aqade la racha
        durArr = c( durArr , x2-x1 ) #Aqade la duracion de la racha
      } 
      else if( y2 < abj ){
        #Cambia variable auxiliar no permite que aumente
        aux = 1
        
        #Crear recta
        m = (y2-y1)/(i+1-i)
        b = y2-m*(i+1)
        
        x1 =  (abj - b)/m #x donde inicia racha
        
        A = abs(y2+abj)/2*(i+1-x1) - abs(abj)*(i+1-x1) #Area 1 por debajo 
        
        while ( y2 < abj){
          i = i+1 #Aumenta el i o el paso
          y1 = x[i]
          y2 = x[i+1]
          
          if (y2 >= abj){
            #Calcular el x2
            #Crear recta
            m = (y2-y1)/(i+1-i)
            b = y2-m*(i+1)
            
            x2 =  (abj - b)/m #x donde inicia racha
            
            #Calcula el area de forma distinta
            A2 = abs(y1+abj)/2*(x2-i) - abs(abj)*(x2-i) #Area 2 por encima 
            A = c(A,A2)
            
          }else {
            # Area acumulativa
            A = c(A, abs(y2+y1)/2*(i+1-i) - abs(abj)*(i+1-i)) #Area por debajo 
          }
        }
        ArAbj = c( ArAbj , sum (A)) #Aqade la racha
        durAbj = c( durAbj , x2-x1 ) #Aqade la duracion de la racha
        
      }
    }
    
    
    
    if (aux == 0){
      i = i + 1 #Aumento del ciclo
    }
    
  }
  #Igualar numero de filas
  max.len = max(length(ArArr), length(ArAbj))
  durArr = c(durArr, rep(NA, max.len - length(durArr)))
  ArArr = c(ArArr, rep(NA, max.len - length(ArArr)))
  durAbj = c(durAbj, rep(NA, max.len - length(durAbj)))
  ArAbj = c(ArAbj, rep(NA, max.len - length(ArAbj)))
  
  R = data.frame(durArr,ArArr,durAbj,ArAbj) #Resultados
  return(R)
} 

Rachas(Real,2)
Rachas(YPron1,2)
Rachas(YPron2,2)

#Error cuadratico medio
RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}

RMSE(Real,YPron1)
RMSE(Real,YPron2)

#NASH
Nash <- function(R,P){
  #R reales
  #P prediccion
  aux1 = c() #Sumatarias
  aux2 = c()
  n = min(length(R),length(P))
  
  for (i in 1:n){
    aux1 = c( aux1, (R[i]-P[i])^2 )
    aux2 = c( aux2, (R[i]-mean(R))^2 )
  }
  aux1 <- sum(aux1)
  aux2 <- sum(aux2)
  
  NS = 1 - aux1/aux2
  return(NS)
}

Nash(Real,YPron1)
Nash(Real,YPron2)


#CICLO ANUAL ERROR CUADRATICO MEDIO

Qmed <- filter(data, 1997 <= as.numeric(format.Date(as.Date(data$Fecha,"%Y"), "%Y")) & 2011 >= as.numeric(format.Date(as.Date(data$Fecha,"%Y"), "%Y")))  %>% filter(IdParametro == "CAUDAL") %>% select( Fecha, Valor )  

Qmed$Fecha <- as.Date(Qmed$Fecha)
#  Get months
Qmed$Month <- months(Qmed$Fecha)
#  Get years
Qmed$Year <- format(Qmed$Fecha,format="%y")

#  Aggregate 'X2' on months and year and get mean
Ciclo <- aggregate( Valor ~ Month  , Qmed , mean )

RMSE


