# Funcion que encuentra la derivada aproximada en un punto
derivada <- function(x0,xi,y0,yi)
{
  a <- xi-x0
  b <- yi-y0
  return(b/a)
}
#Funcion que calculo el primer termino del polinomio en la formula de hermite
H1 <- function(o,i,x)
{
  u <- (x-o)/(i-o)
  return((2*(u**3))-(3*(u**2))+1)
}
#Funcion que calculo el segundo termino del polinomio en la formula de hermite
H2 <- function(o,i,x)
{
  u <- (x-o)/(i-o)
  return((-2*(u**3))+(3*(u**2)))
}
#Funcion que calculo el tercer termino del polinomio en la formula de hermite
H3 <- function(o,i,x)
{
  u <- (x-o)/(i-o)
  return((u**3)-(2*(u**2))+u)
}
#Funcion que calculo el cuarto termino del polinomio en la formula de hermite
H4 <- function(o,i,x)
{
  u <- (x-o)/(i-o)
  return((u**3)-(u**2))
}

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

#Funcion que grafica los puntos dados y la curva de las interpolaciones
#Recibe los valores de los splines y los vectores que contienen los puntos
grafica <- function(xH,yH,X,Y)
{
  plot(X,Y, pch=19, cex=1, col = "red", asp = 1, xlab="X", ylab = "Y", main = "Spline de Hermite" )
  lines(xH,yH, col="blue")
}

#Funcion principal de splines de Hermite
#Recibe dos vectores con los puntos enviados
#Llama las funciones derivada, H1,H2,H3,H4 y las funciones para imprimir las tablas y la grafica
#Para sacar las derivadas se tiene encuenta si es un punto inicial, final o intermedio
splineHermite <- function(X,Y)
{
  xHermite <- c()
  yHermite <- c()
  xpto <- X[1]
  c <- 1
  d <- 2
  for(i in 2:length(X)-1)
  {
    xHermite <- c(xHermite,X[i])
    yHermite <- c(yHermite,Y[i])

    lIntervalo <- (X[d]-X[c])/5
    if(i == 1)
    {
      d0 <- derivada(X[i],X[i+1],Y[i],Y[i+1])
      di <- derivada(X[i],X[i+2],Y[i],Y[i+2])
    }
    if(i == length(X)-1)
    {
      d0 <- derivada(X[i-1],X[i+1],Y[i-1],Y[i+1])
      di <- derivada(X[i],X[i+2],Y[i],Y[i+2])
    }
    if(i == length(X)-1)
    {
      d0 <- derivada(X[i-1],X[i+1],Y[i-1],Y[i+1])
      di <- derivada(X[i],X[i+1],Y[i],Y[i+1])
    }
    xpto <- X[c]

    for(j in 1:4)
    {
      xpto <- xpto + lIntervalo
      xHermite <- c(xHermite,xpto)
      if(i == 4)
      {
        h1 <- H1(X[c],X[d],xpto)
        h2 <- H2(X[c],X[d],xpto)
        h3 <- H3(X[c],X[d],xpto)
        h4 <- H4(X[c],X[d],xpto)
        ypto <- (Y[c]*h1)+(Y[d]*h2)+(d0*h3)+(di*h4)
        yHermite <- c(yHermite,ypto)

      }
      if(i != 4)
      {
        h1 <- H1(X[c],X[d],xpto)
        h2 <- H2(X[c],X[d],xpto)
        h3 <- H3(X[c],X[d],xpto)
        h4 <- H4(X[c],X[d],xpto)
        ypto <- (Y[c]*h1)+(Y[d]*h2)+(d0*h3)+(di*h4)
        yHermite <- c(yHermite,ypto)
      }

    }
    c <- c+1
    d <-d+1
  }
  xHermite <- c(xHermite,X[length(X)])
  yHermite <- c(yHermite,Y[length(Y)])
  f <- splinefun(X,Y,method = "natural", ties = mean)
  a <- 1
  b <- 5
  for(i in 2:length(X))
  {
    x <- xHermite[c(a:b)]
    y <- yHermite[c(a:b)]
    err <- error(x,f,y,length(x))
    tabla(x,y,err)
    a <- a+5
    b <- b+5
  }
  grafica(xHermite,yHermite,X,Y)

}

