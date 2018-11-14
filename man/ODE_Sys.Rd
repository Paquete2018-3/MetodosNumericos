\name{ODE_Sys}
\alias{ODE_Sys}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Función ODE_Sys
%%  ~~function to do ... ~~
}
\description{
Esta es la principal función del paquete y la que debe ser llamada por el usuario, ya que a partir de esta función se llama a las demás funciones con el fin de cumplir el objetivo del paquete que es resolver el sistema de ecuaciones diferenciales.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ODE_Sys(d, vars, init, h, lims, n, method, point)
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
  \item{lims}{
  Vector con los límites desados de la solución. Ejemlo -> [x_liminf,x_limsup,y_liminf,y_limsup]
%%     ~~Describe \code{lims} here~~
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
Si el usuario no ingresa la constante h, esta será definida por default como 0.1. Así mismo, si el usuario no ingresa el número de puntos para hallar la solución, este será definido por default como 250.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
La función retorna una lista de valores correspondientes a la solución numérica para cada valor de la variable independiente, así como la solución aproximada para un valor de la variable independiente. Así mismo también se encuentra como salida la gráfica de la solución, así como el campo de pendientes, únicamente si el sistema ingresado por el usuario es de segundo orden.
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
sol = ODE_Sys(c("-x+y+1","x-y"),c("x","y","t"),c(-3,1,0),0.1,c(-4,4,-3,3),500,"rk4",0.62)

## The function is currently defined as
function (d, vars, init, h, lims, n, method, point)
{
    if (missing(h)) {
        h = 0.1
    }
    if (missing(n)) {
        n = 250
    }
    values = Solucion(d, vars, init, h, n, method, point)
    aux = point - values[length(d) + 1, ]
    index = which.min(abs(aux))
    solution_point = values[, index]
    solution <- list(trajectory = values, answer = solution_point)
    plotfield(field, lims, d, vars)
    return(solution)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
