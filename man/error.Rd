\name{error}
\alias{error}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Error Función
%%  ~~function to do ... ~~
}
\description{
Esta función encuentra el error de la solución comparando la solución analítica del sistema que ya se conoce, con el valor obtenido en la solución.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
error(solucion)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{solucion}{
  La solución tanto en el eje x como en el eje y así como el valor t que se evaluó para encontrar dicha solución.
%%     ~~Describe \code{solucion} here~~
}
}
\details{
La solución analítica que se encuentra en la función hace referencia al sistema de "-x+y+1" y "x-y".
%%  ~~ If necessary, more details than the description above ~~
}
\value{
El valor retornado es un arreglo de tamaño 2, que contiene el error tanto en el eje x como el error en el eje y.
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (solucion)
{
    t = solucion[3]
    xanalitico = (1/4) * (2 * t - 9 * exp(-2 * t) - 3)
    yanalitico = (1/4) * (2 * t + 9 * exp(-2 * t) - 5)
    xerror = abs(xanalitico - solucion[1])
    yerror = abs(yanalitico - solucion[2])
    errores <- c(xerror, yerror)
    return(errores)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
