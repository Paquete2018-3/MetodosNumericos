\name{solucion}
\alias{solucion}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Función Solución
%%  ~~function to do ... ~~
}
\description{
Esta función resuelve el sistema de ecuaciones diferenciales, es llamada dentro de la función principal ODE_Sys con el fin de hallar los valores que conforman la solución.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
solucion (d, vars, init, h, n, method, point)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
  Representa al vector de las ecuaciones diferenciales que componen el sistema a solucionar. Ejemplo -> ["Eqn1","Eqn2","Eqn3", ...]
%%     ~~Describe \code{d} here~~
}
  \item{vars}{
  Representa al vector con los nombre de la variables incluyendo la variable independiente (t) al final del vector. Ejemplo -> ["x1","x2","x3", ...,"t"]
%%     ~~Describe \code{vars} here~~
}
  \item{init}{
  Vector con los valores iniciales incluyendo la variable indepediente (t) al final del vector. Ejemplo -> [x1_0,x2_0,x3_0, ...,t_0]
%%     ~~Describe \code{init} here~~
}
  \item{h}{
  Constante.
%%     ~~Describe \code{h} here~~
}
  \item{n}{
  Número de puntos de la solución que se requieren.
%%     ~~Describe \code{n} here~~
}
  \item{method}{
  Método para encontrar la solución. Se encoge entre las siguientes opciones -> ["euler","midpoint",rk4"]
%%     ~~Describe \code{method} here~~
}
  \item{point}{
  Punto deseado para evaluar la solución.
%%     ~~Describe \code{point} here~~
}
}
\details{
La solución claramente dependerá principalmente del método a implantar para resolver el sistema de ecuaciones diferenciales.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
La función retorna una lista de valores correspondientes a la solución numérica para cada valor de la variable independiente.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Pablo Millan y Oscar Fonseca
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Ejemplos:
#sol=solucion(c("x*y","x*y-1"),c("x","y","t"),c(0.5,0.5,0),0.1,500,"euler",0.62)
#sol=solucion(c("x*y","x*y-1"),c("x","y","t"),c(-2,1,0),0.01,400,"midpoint",0.5)
#sol=solucion(c("x*y","x*y-1"),c("x","y","t"),c(-2,2,0),0.01,700,"rk4",0.674)
#sol=solucion(c("x-y+z","x*y-1","z+3*x"),c("x","y","z","t"),c(0,1,1,0),0.01,450,"euler",0.5)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
Solucion <- function(d, vars, init, h, n, method, point){
  if (method=="euler"){
    values=matrix(init,ncol = 1)

    for (i in 1:length(d)){
      assign(vars[i],init[i])
    }

    npoints=1
    values=cbind(values,numeric(length(d)+1))

    while (npoints<n) {
      for (i in 1:length(d)){
        values[i,npoints+1]=values[i,npoints]+h*eval(parse(text=d[i]))
      }

      for (i in 1:length(d)){
        assign(vars[i],values[i,npoints+1])
      }

      values[length(d)+1,npoints+1]=values[length(d)+1,npoints]+h

      assign(vars[length(d)+1],values[length(d)+1,npoints+1])

      values=cbind(values,numeric(length(d)+1))
      npoints=npoints+1

    }

    values=values[,-npoints-1]
  }

  if (method=="midpoint"){
    values=matrix(init,ncol = 1)

    for (i in 1:length(d)){
      assign(vars[i],init[i]+h/2)
    }

    npoints=1
    values=cbind(values,numeric(length(d)+1))

    while (npoints<n) {
      for (i in 1:length(d)){
        values[i,npoints+1]=values[i,npoints]+h*eval(parse(text=d[i]))
      }

      for (i in 1:length(d)){
        assign(vars[i],values[i,npoints+1]+h/2)
      }

      values[length(d)+1,npoints+1]=values[length(d)+1,npoints]+h

      assign(vars[length(d)+1],values[length(d)+1,npoints+1])

      values=cbind(values,numeric(length(d)+1))
      npoints=npoints+1

    }

    values=values[,-npoints-1]
  }

  if (method=="rk4"){
    values=matrix(init,ncol = 1)

    for (i in 1:length(d)){
      assign(vars[i],init[i]+h/2)
    }

    npoints=1
    values=cbind(values,numeric(length(d)+1))

    while (npoints<n) {
      k1=numeric(length(d))
      k2=numeric(length(d))
      k3=numeric(length(d))
      k4=numeric(length(d))

      for (i in 1:length(d)){
        k1[i]=h*eval(parse(text=d[i]))
        assign(vars[i],values[i,npoints]+k1[i]/2)
      }
      for (i in 1:length(d)){
        k2[i]=h*eval(parse(text=d[i]))
        assign(vars[i],values[i,npoints]+k2[i]/2)
      }
      for (i in 1:length(d)){
        k3[i]=h*eval(parse(text=d[i]))
        assign(vars[i],values[i,npoints]+k3[i])
      }
      for (i in 1:length(d)){
        k4[i]=h*eval(parse(text=d[i]))
      }
      for (i in 1:length(d)){
        values[i,npoints+1]=values[i,npoints]+ k1[i]/6 + k2[i]/3 + k3[i]/3 + k4[i]/6
      }
      for (i in 1:length(d)){
        assign(vars[i],values[i,npoints+1])
      }
      values[length(d)+1,npoints+1]=values[length(d)+1,npoints]+h
      assign(vars[length(d)+1],values[length(d)+1,npoints+1])

      values=cbind(values,numeric(length(d)+1))
      npoints=npoints+1
    }
    values=values[,-npoints-1]
  }
  return (values)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
