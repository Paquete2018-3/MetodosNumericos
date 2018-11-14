#Calculo del error de interpolaci?n
error<-function(x,f,s,n)
{
  err<-0
  sol<-0
  ans<-0
  for(i in 1:n)
  {
    h<-max(x[i+1]-x[i])
    err=h^(3/2)
    ans = f(x[i])-s[i]
    if(!is.null((ans< err))||!is.null((ans==err)))
    {
      sol[i]=ans
    }
  }
  return(sol)
}

#Imprime los resultados pertinentes
tabla<-function(x,y,err)
{
  datos<-data.frame(cbind(x,y,err))
  colnames(datos)<-c('X','Spline','Error')
  print(datos)
}

#Calculo del tama?o de paso "h"
#De la forma que h(i) = x(i+i) - x(i)
calHi<-function(x)
{
  n = length(x)
  h <- (c(x,0) - c(0,x))[2:n]
  return (h)
}

#Calculo de la matriz A que hace parte del sistema de ecuaciones As=B
#Es de la forma que cuando esta en la iteraci?n i
#Cuando (i,j) i = j = 2(hi+hi+1)
#Cuando (i,j-1) = hi
#Cuando (i,j+1) = hi+1
matrizA<-function(h,a,x)
{
  n = length(x)
  alph <- (3/c(1,h,1,1)*(c(a,1,1) - c(1,a,1)) - 3/c(1,1,h,1)*(c(1,a,1)-c(1,1,a)))[3:n]
  A <- c(1,rep(0,times=n-1))
  for (i in 1:(n-2)) {
    A <- rbind(A,c( rep(0,times=i-1) , c(h[i],2*(h[i]+h[i+1]),h[i+1]) , rep(0,times=n-i-2) ) )
  }
  A <- rbind(A,c(rep(0,times=n-1),1))

  return(A)
}

#Calculo de la matriz B que hace parte del sistema de ecuaciones As=B
#Es de la forma que cuando esta en la iteraci?n i
#Cuando (i,1) = 6((yi+2-yi+1)/hi+1-(yi+1-yi)/hi)
matrizB<-function(h,a)
{
  n = length(x)
  alph <- (3/c(1,h,1,1)*(c(a,1,1) - c(1,a,1)) - 3/c(1,1,h,1)*(c(1,a,1)-c(1,1,a)))[3:n]
  b <- c(0,alph,0)
  return (b)
}

#Despeje de la ecuaci?n de la forma AS=B en el cual obtenemos los valores de Spline
#que nos ayudar? a calcular los polinomios de interpolaci?n
calS<-function(A,b)
{
  c <- solve(A, b)
  return(c)
}

#Calculo de los coeficientes del polinomio de interpolaci?n Spline
#Recibe los spline del despeje de la ecuaci?n As=B
#los valores en X, Y y el tama?o de paso (h)
sNormal<-function(a,c,h)
{
  n = length(x)
  b <- ((c(a,0) - c(0,a))/c(1,h,1) - c(1,h,1)/3*(c(c,0) + 2*c(0,c)))[2:n]
  d <- ((c(c,0) - c(0,c))/(3*c(1,h,1)))[2:n]

  ans = rbind(a[1:n-1],b,c[1:n-1],d)
}

#Splie cubico, llama a las respectivas funciones para que la interpolacion funcione
sCubico<-function(x,y)
{
  h<-calHi(x)
  A<-matrizA(h,y,x)
  B<-matrizB(h,y)
  S<-calS(A,B)
  sNormal(y,S,h)
}

#Funcion que llama a las funciones de interpolar
#Graficar y calcular el error
splineCubico<-function(x,y)
{
  plot(x,y, pch=19, cex=1, col = "red", asp = 1, xlab="X", ylab = "Y", main = "Spline Cubico")
  f = splinefun(x,y,method = "natural", ties = mean)
  t = 1:length(x)
  sx = sCubico(t,x)
  sy = sCubico(t,y)
  for (i in 1:(length(t)-1)) {
    dat<- data.frame(t=seq(t[i],t[i+1], by=0.1) )
    fx <- function(x) (sx[1,i] + sx[2,i]*(x-t[i]) + sx[3,i]*(x-t[i])^2 + sx[4,i]*(x-t[i]))
    fy <- function(x) (sy[1,i] + sy[2,i]*(x-t[i]) + sy[3,i]*(x-t[i])^2 + sy[4,i]*(x-t[i]))
    dat$y=fy(dat$t)
    dat$x=fx(dat$t)
    points(dat$x,dat$y,type='l', col='blue')
    err=error(dat$x,f,dat$y,length(dat$x))
    tabla(dat$x,dat$y,err)
  }
}
