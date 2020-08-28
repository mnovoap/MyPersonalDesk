#Cargar Librerias
library(rafalib)
library(dplyr)
library(fANCOVA)
library(forecast)
library(TSA)
library(dplyr)
library(rafalib)


#Cargar funciones 
#Necesario los Datos ya como vector 
#NOTAAAAAAAAA: HAY QUE CORREGIR EL ALGORITMO PARA HALLAR EL LA FRECUENCIA EN LAS TABLAS

chi_norm= function(datos)   {
  #Parametros
  mean=mean(datos)
  desv=sd(datos)
  #Intervalos
  r=max(datos)-min(datos)
  n=length(datos)
  #Criterios para los intervalos
  k1=1+3.3*log(n,base=2) #Sturges
  k2=2*(n^0.5) #Valleman para muestras pequeñas 50
  k3=10*log10(n) #Dixon y Kronmal, n grande 50
  #Escogencia y ancho final de clases
  if (k1<=k2 & k1<=k3){
    k=k1
  } else if (k2<=k3 & k2<=k1){
    k=k2
  } else {
    k=k3
  }
  c=(r/k)
  #TABLA
  for(i in 1:(ceiling(r/c))){
    if(i==1){
      inf=min(datos)
      sup=inf+c
      INT=c(i)
      P=paste(round(inf,2),"-",round(sup,2))
      RAN=c(P)
      aux=sum(datos<=sup & datos>=inf)
      conteo=c(aux)
      f=conteo/n
      fac=f
      z=c((sup-mean)/desv)
      facre=c(pnorm(z[i]))
      LSR=c(sup)
      fre=c(facre[(i)])
      xi=c(n*((f[(i)]-fre[(i)])^2)/fre[(i)])
    } else if (i==(ceiling(r/c))) {
        inf=sup
        sup=max(datos)
        INT=c(INT,i)
        P=paste(round(inf,2),"-",round(sup,2))
        RAN=c(RAN,P)
        aux=sum(datos<=sup & datos>inf)
        conteo=c(conteo,aux)
        f=conteo/n
        fac=c(fac,fac[(i-1)]+f[(i)])
        z=c(z,(sup-mean)/desv)
        facre=c(facre,pnorm(z[i]))
        LSR=c(LSR,sup)
        fre=c(fre,facre[(i)]-facre[(i-1)])
        xi=c(xi,n*((f[(i)]-fre[(i)])^2)/fre[(i)])
        xit=sum(xi)
        
        #Tabla Final
        final=c(rep(0,ceiling(r/c)-1),xit)
        aux=c(conteo,f,fac,LSR,z,facre,fre,xi,final)
        aux=round(aux,2)
        resul=c(INT,RAN,aux)
        tabla=cbind(matrix(resul,ceiling(r/c),11))
        colnames(tabla)=c("INTERVALO","RANGO","n","fo(x)","Fo(x)",
                          "LSR","Z","F(x)","f(x)","xi2","xt2")
        
        
    } else {
        inf=sup
        sup=inf+c
        INT=c(INT,i)
        P=paste(round(inf,2),"-",round(sup,2))
        RAN=c(RAN,P)
        aux=sum(datos<=sup & datos>inf)
        conteo=c(conteo,aux)
        f=conteo/n
        fac=c(fac,fac[(i-1)]+f[(i)])
        z=c(z,(sup-mean)/desv)
        facre=c(facre,pnorm(z[i]))
        LSR=c(LSR,sup)
        fre=c(fre,facre[(i)]-facre[(i-1)])
        xi=c(xi,n*((f[(i)]-fre[(i)])^2)/fre[(i)])
      }
  }
  return(tabla)
}
chi_lognorm=function(datos){
  #Cambio de dato a logaritmos
  datos=log(datos,base=exp(1))
  #Parametros
  mean=mean(datos)
  desv=sd(datos)
  #Intervalos
  r=max(datos)-min(datos)
  n=length(datos)
  #Criterios para los intervalos
  k1=1+3.3*log(n,base=2) #Sturges
  k2=2*(n^0.5) #Valleman para muestras pequeñas 50
  k3=10*log10(n) #Dixon y Kronmal, n grande 50
  #Escogencia y ancho final de clases
  if (k1<=k2 & k1<=k3){
    k=k1
  } else if (k2<=k3 & k2<=k1){
    k=k2
  } else {
    k=k3
  }
  c=(r/k)
  #TABLA
  for(i in 1:(ceiling(r/c))){
    if(i==1){
      inf=min(datos)
      sup=inf+c
      INT=c(i)
      P=paste(round(inf,2),"-",round(sup,2))
      RAN=c(P)
      aux=sum(datos<=sup & datos>=inf)
      conteo=c(aux)
      f=conteo/n
      fac=f
      z=c((sup-mean)/desv)
      facre=c(pnorm(z[i]))
      LSR=c(sup)
      fre=c(facre[(i)])
      xi=c(n*((f[(i)]-fre[(i)])^2)/fre[(i)])
    } else if (i==(ceiling(r/c))) {
      inf=sup
      sup=max(datos)
      INT=c(INT,i)
      P=paste(round(inf,2),"-",round(sup,2))
      RAN=c(RAN,P)
      aux=sum(datos<=sup & datos>inf)
      conteo=c(conteo,aux)
      f=conteo/n
      fac=c(fac,fac[(i-1)]+f[(i)])
      z=c(z,(sup-mean)/desv)
      facre=c(facre,pnorm(z[i]))
      LSR=c(LSR,sup)
      fre=c(fre,facre[(i)]-facre[(i-1)])
      xi=c(xi,n*((f[(i)]-fre[(i)])^2)/fre[(i)])
      xit=sum(xi)
      
      #Tabla Final
      final=c(rep(0,ceiling(r/c)-1),xit)
      aux=c(conteo,f,fac,LSR,z,facre,fre,xi,final)
      aux=round(aux,2)
      resul=c(INT,RAN,aux)
      tabla=cbind(matrix(resul,ceiling(r/c),11))
      colnames(tabla)=c("INTERVALO","RANGO","n","fo(x)","Fo(x)",
                        "LSR","Z","F(x)","f(x)","xi2","xt2")
      
      
    } else {
      inf=sup
      sup=inf+c
      INT=c(INT,i)
      P=paste(round(inf,2),"-",round(sup,2))
      RAN=c(RAN,P)
      aux=sum(datos<=sup & datos>inf)
      conteo=c(conteo,aux)
      f=conteo/n
      fac=c(fac,fac[(i-1)]+f[(i)])
      z=c(z,(sup-mean)/desv)
      facre=c(facre,pnorm(z[i]))
      LSR=c(LSR,sup)
      fre=c(fre,facre[(i)]-facre[(i-1)])
      xi=c(xi,n*((f[(i)]-fre[(i)])^2)/fre[(i)])
    }
  }
  return(tabla)
}

chi_gumbell=function(datos){
  #Parametros
  mx=mean(datos)
  sx=sd(datos)
  #Intervalos
  r=max(datos)-min(datos)
  n=length(datos)
  #Parametros Gumbell
  sy=0.0262*log(n,base=exp(1))+0.4436
  my=0.133*log(n,base=exp(1))+0.7467
  a=sy/sx
  b=mx-my/sx
  
  #Criterios para los intervalos
  k1=1+3.3*log(n,base=2) #Sturges
  k2=2*(n^0.5) #Valleman para muestras pequeñas 50
  k3=10*log10(n) #Dixon y Kronmal, n grande 50
  #Escogencia y ancho final de clases
  if (k1<=k2 & k1<=k3){
    k=k1
  } else if (k2<=k3 & k2<=k1){
    k=k2
  } else {
    k=k3
  }
  c=(r/k)
  #TABLA
  for(i in 1:(ceiling(r/c))){
    if(i==1){
      inf=min(datos)
      sup=inf+c
      INT=c(i)
      P=paste(round(inf,2),"-",round(sup,2))
      RAN=c(P)
      aux=sum(datos<=sup & datos>=inf)
      conteo=c(aux)
      f=conteo/n
      fac=f
      cg=c((sup-b)*a)
      facre=c(exp(-1*exp(-1*cg[i])))
      LSR=c(sup)
      fre=c(facre[(i)])
      xi=c(n*((f[(i)]-fre[(i)])^2)/fre[(i)])
    } else if (i==(ceiling(r/c))) {
      inf=sup
      sup=max(datos)
      INT=c(INT,i)
      P=paste(round(inf,2),"-",round(sup,2))
      RAN=c(RAN,P)
      aux=sum(datos<=sup & datos>inf)
      conteo=c(conteo,aux)
      f=conteo/n
      fac=c(fac,fac[(i-1)]+f[(i)])
      cg=c(cg,(sup-b)*a)
      facre=c(facre,exp(-1*exp(-1*cg[i])))
      LSR=c(LSR,sup)
      fre=c(fre,facre[(i)]-facre[(i-1)])
      xi=c(xi,n*((f[(i)]-fre[(i)])^2)/fre[(i)])
      xit=sum(xi)
      
      #Tabla Final
      final=c(rep(0,ceiling(r/c)-1),xit)
      aux=c(conteo,f,fac,LSR,cg,facre,fre,xi,final)
      aux=round(aux,2)
      resul=c(INT,RAN,aux)
      tabla=cbind(matrix(resul,length(RAN),11))
      colnames(tabla)=c("INTERVALO","RANGO","n","fo(x)","Fo(x)",
                        "LSR","C","F(x)","f(x)","xi2","xt2")
      
      
    } else {
      inf=sup
      sup=inf+c
      INT=c(INT,i)
      P=paste(round(inf,2),"-",round(sup,2))
      RAN=c(RAN,P)
      aux=sum(datos<=sup & datos>inf)
      conteo=c(conteo,aux)
      f=conteo/n
      fac=c(fac,fac[(i-1)]+f[(i)])
      cg=c(cg,(sup-b)*a)
      facre=c(facre,exp(-1*exp(-1*cg[i])))
      LSR=c(LSR,sup)
      fre=c(fre,facre[(i)]-facre[(i-1)])
      xi=c(xi,n*((f[(i)]-fre[(i)])^2)/fre[(i)])
    }
  }
  return(tabla)
  
  
}

chi_logumbell=function(datos){
  datos=log(datos,base=exp(1))
  #Parametros
  mx=mean(datos)
  sx=sd(datos)
  #Intervalos
  r=max(datos)-min(datos)
  n=length(datos)
  #Parametros Gumbell
  sy=0.0262*log(n,base=exp(1))+0.4436
  my=0.133*log(n,base=exp(1))+0.7467
  a=sy/sx
  b=mx-my/sx
  
  #Criterios para los intervalos
  k1=1+3.3*log(n,base=2) #Sturges
  k2=2*(n^0.5) #Valleman para muestras pequeñas 50
  k3=10*log10(n) #Dixon y Kronmal, n grande 50
  #Escogencia y ancho final de clases
  if (k1<=k2 & k1<=k3){
    k=k1
  } else if (k2<=k3 & k2<=k1){
    k=k2
  } else {
    k=k3
  }
  c=(r/k)
  #TABLA
  for(i in 1:(ceiling(r/c))){
    if(i==1){
      inf=min(datos)
      sup=inf+c
      INT=c(i)
      P=paste(round(inf,2),"-",round(sup,2))
      RAN=c(P)
      aux=sum(datos<=sup & datos>=inf)
      conteo=c(aux)
      f=conteo/n
      fac=f
      cg=c((sup-b)*a)
      facre=c(exp(-1*exp(-1*cg[i])))
      LSR=c(sup)
      fre=c(facre[(i)])
      xi=c(n*((f[(i)]-fre[(i)])^2)/fre[(i)])
    } else if (i==(ceiling(r/c))) {
      inf=sup
      sup=max(datos)
      INT=c(INT,i)
      P=paste(round(inf,2),"-",round(sup,2))
      RAN=c(RAN,P)
      aux=sum(datos<=sup & datos>inf)
      conteo=c(conteo,aux)
      f=conteo/n
      fac=c(fac,fac[(i-1)]+f[(i)])
      cg=c(cg,(sup-b)*a)
      facre=c(facre,exp(-1*exp(-1*cg[i])))
      LSR=c(LSR,sup)
      fre=c(fre,facre[(i)]-facre[(i-1)])
      xi=c(xi,n*((f[(i)]-fre[(i)])^2)/fre[(i)])
      xit=sum(xi)
      
      #Tabla Final
      final=c(rep(0,ceiling(r/c)-1),xit)
      aux=c(conteo,f,fac,LSR,cg,facre,fre,xi,final)
      aux=round(aux,2)
      resul=c(INT,RAN,aux)
      tabla=cbind(matrix(resul,length(RAN),11))
      colnames(tabla)=c("INTERVALO","RANGO","n","fo(x)","Fo(x)",
                        "LSR","C","F(x)","f(x)","xi2","xt2")
      
      
    } else {
      inf=sup
      sup=inf+c
      INT=c(INT,i)
      P=paste(round(inf,2),"-",round(sup,2))
      RAN=c(RAN,P)
      aux=sum(datos<=sup & datos>inf)
      conteo=c(conteo,aux)
      f=conteo/n
      fac=c(fac,fac[(i-1)]+f[(i)])
      cg=c(cg,(sup-b)*a)
      facre=c(facre,exp(-1*exp(-1*cg[i])))
      LSR=c(LSR,sup)
      fre=c(fre,facre[(i)]-facre[(i-1)])
      xi=c(xi,n*((f[(i)]-fre[(i)])^2)/fre[(i)])
    }
  }
  return(tabla)
  
  
}
#Smirnov_Kolgomorov
smirnov_normal=function(datos){
  mx=mean(datos)
  sx=sd(datos)
  n=length(datos)
  m=c(1:n)
  xo=sort(datos,decreasing=TRUE)
  x=(xo-mx)/sx
  fo=1-m/(n+1)
  f=NULL
  for (i in 1:n){
    f=c(f,pnorm(x[i]))
  }
  D=abs(f-fo)
  Df=max(D)
  return(paste("Valor critico",Df ))
}

smirnov_lognormal=function(datos){
  datos=log(datos,base=exp(1))
  mx=mean(datos)
  sx=sd(datos)
  n=length(datos)
  m=c(1:n)
  xo=sort(datos,decreasing=TRUE)
  x=(xo-mx)/sx
  fo=1-m/(n+1)
  f=NULL
  for (i in 1:n){
    f=c(f,pnorm(x[i]))
  }
  D=abs(f-fo)
  Df=max(D)
  return(paste("Valor critico",Df ))
}
smirnov_gumbel=function(datos){
  mx=mean(datos)
  sx=sd(datos)
  n=length(datos)
  #Parametros Gumbell
  sy=0.0262*log(n,base=exp(1))+0.4436
  my=0.133*log(n,base=exp(1))+0.7467
  a=sy/sx
  b=mx-my/sx
  
  #Proceso prueba
  m=c(1:n)
  x=sort(datos,decreasing=TRUE)
  fo=1-m/(n+1)
  f=NULL
  cg=NULL
  for (i in 1:n){
    cg=c(cg,((x[i]-b)*a))
    f=c(f,exp(-1*exp(-1*cg[i])))
  }
  D=abs(f-fo)
  Df=max(D)
  return(paste("Valor critico",Df ))
}
smirnov_loggumbel=function(datos){
  datos=log(datos,base=exp(1))
  mx=mean(datos)
  sx=sd(datos)
  n=length(datos)
  #Parametros Gumbell
  sy=0.0262*log(n,base=exp(1))+0.4436
  my=0.133*log(n,base=exp(1))+0.7467
  a=sy/sx
  b=mx-my/sx
  
  #Proceso prueba
  m=c(1:n)
  x=sort(datos,decreasing=TRUE)
  fo=1-m/(n+1)
  f=NULL
  cg=NULL
  for (i in 1:n){
    cg=c(cg,((x[i]-b)*a))
    f=c(f,exp(-1*exp(-1*cg[i])))
  }
  D=abs(f-fo)
  Df=max(D)
  return(paste("Valor critico",Df ))
}

p_rachas_Hidr=function(datos){
  #Estandarización caudales
  med=mean(datos)
  des=sd(datos)
  datos=(datos-med)/des
  n=length(datos)
  #Clasificacion datos
  for (i in 1:n){
    if(datos[i]>0){
      datos[i]=1
    }else if(datos[i]<0){
      datos[i]=-1
    }else{
      datos[i]=NA
    }
  }
  datos=na.omit(datos)
  
  #Organizar rachas
  a = 1
  b = 0
  c = 0
  d = 1
  aux1 = NA
  aux2 = NA
  r = NA
  i = 1
  
  
  while (i<(n+1)) {
    
    if(d == 1){
      aux1[i]=datos[i]
      e=aux1[i]
      d=0
      i=i+1
      if(datos[i] != e){
        r[a]=length(aux1)
        a = a+1
      }    
      }else if (datos[i] == e | length(aux1) == 0){
        
        aux1[i]=datos[i]
        i=i+1
        
        while(b<1){
          if(datos[i]==aux1[i-1] & i<(n+1)){
            aux1=c(aux1,datos[i])
            e=aux1[i]
            i=i+1
          }else{
            aux1 = na.omit(aux1)
            r[a]=length(aux1)
            a = a+1
            b=1
            c=0
            aux2 = NULL
          }
        }
        
      }else{
        aux2[i]=datos[i]
        i=i+1
        
        while (c<1) {
          if(datos[i]==aux2[i-1] & i<(n+1)){
            aux2=c(aux2,datos[i])
            i=i+1
          }else{
            e=datos[i]
            aux2 = na.omit(aux2)
            r[a]=length(aux2)
            a = a+1
            c = 1
            b = 0
            aux1 = NULL
          }
        }
        
      }
  }
  
  racha=length(r)
  
  #n1 y n2
  t=length(r)
  n1=c()
  n2=c()
  for (i in 1:t){
    if(i%%2 != 0){
      n1<-c(n1,r[i])
    }else{
      n2<-c(n2,r[i])
    }
  }
  n1=sum(n1)
  n2=sum(n2)
  l=c(n1,n2)
  l=sort(l,decreasing=TRUE)
  n1=l[1]
  n2=l[2]
  
  
  #Prueba de aleatoriedad
  if(n1<20 | n2<20){
    sal=c("Buscar en las tablas (n1, n2 y r):",n1,n2,racha)
  } else {
    ur=2*n1*n2/(n1+n2)+1
    sr=(2*n2*n1*(2*n1*n2-n1-n2)/((n1+n2)^2*(n1+n2-1)))^0.5
    z=(racha-ur)/sr
    z=abs(z)
    vp=(1-pnorm(z))*2
    if(vp<=0.05){
      sal=c("Se rechaza la aleatoriedad ", vp)
    }else{
      sal=c("Se acepta la aleatoriedad ", vp)
      }
  }
return(sal)
}
 



#HIDROLOGIA

#Datos
datos=scan(file.choose())
datos=na.omit(datos)
Q=ts(datos,freq=12,start=c(2000,2))



#Plotear serie de tiempo
plot(Q,type="l",main="RIO MARINILLA",sub="SERIE DE TIEMPO",xlab="Realización con medición mensual",ylab="Caudal MAX (m3/s)")
abline(mean(datos),b=0,col="blue")
media=mean(datos)
sd=sd(datos, na.rm = FALSE)
abline((mean(datos)-sd),b=0,col="red")
abline((mean(datos)+sd),b=0,col="red")

#Boxplot
boxplot(Q~cycle(Q),names=c("ENE","FEB","MAR","ABR","MAY","JUN","JUL","AGO","SEP","OCT","NOV","DIC"),main="Boxplot Caudales")

#Histograma
hist(Q,main="HISTOGRAMA CAUDALES",xlab="Caudales (m3/s)",ylab="Conteo")

# Descomposición  aditiva de los primeros n datos
descom=decompose(Q,type="additive")
plot(descom)

#Componente estacional estimada según el filtro
Tt=descom$trend
St=descom$seasonal

# Grafica componente estacional
plot(St, main="Componente estacional")

# Grafica componente de tendencia
plot(Tt,main="Componente de tendencia")


#cdf y fdp
#CDF
prop=function(q){
  mean(datos<=q)} #Parte en el cual se suma y se saca la probabilidad
#esto es posible porque es una suma de los trues
qs = seq(from=min(datos), to=max(datos), length=15) #Organiza los datos en cuantiles

props=sapply(qs,prop)

plot(qs, props,type="b",main="CDF-CAUDALES MARINILLA")

#fdp
prop2=function(q){
  mean(datos<=(q+7.3769) & datos>=(q-7.3769) )}
qs2=seq(from=min(datos),to=max(datos),length=20)
qs2
props2=sapply(qs2,prop2)
plot(qs2, props2,type="b",main="FDP-CAUDALES MARINILLA")

#Pruebas de Bondad y ajuste
chi_norm(datos)
chi_lognorm(datos)
chi_gumbell(datos)
#Grados de libertad
#v=m-p-1 #m es el numero de intervalos, p los parametros (media y desv)