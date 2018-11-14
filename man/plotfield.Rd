\name{plotfield}
\alias{plotfield}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plotfield función
%%  ~~function to do ... ~~
}
\description{
Función que granfica el campo de pendientes, así como la solución del sistema de ecuaciones diferenciales.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotfield(values, lims, d, vars)
}
\arguments{
%- maybe also 'usage' for other objects documented here.
\item{values}{
  Valores que conforman la solución del sistema.
%%     ~~Describe \code{lims} here~~
}
  \item{lims}{
  Vector con los límites desados de la solución. Ejemlo -> [x_liminf,x_limsup,y_liminf,y_limsup]
%%     ~~Describe \code{lims} here~~
}
  \item{d}{
  Representa al vector de las ecuaciones diferenciales que componen el sistema a solucionar. Ejemplo -> ["Eqn1","Eqn2","Eqn3", ...]
%%     ~~Describe \code{d} here~~
}
  \item{vars}{
  Representa al vector con los nombre de la variables incluyendo la variable independiente (t) al final del vector. Ejemplo -> ["x1","x2","x3", ...,"t"]
%%     ~~Describe \code{vars} here~~
}
}
\details{
Esta función únicamente se llevará a cabo solo si el sistema de ecuaciones diferenciales es de segundo orden.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
La función tiene como salida la grafica de la solución del sistema junto con el campo de pendientes. Similar a la herramienta enviada por la profesora.
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
#plotfield(c(-4,4,-3,3),c("x*y","x*y-1"),c("x","y","t"))
#plotfield(c(-4,4,-3,3),c("x*y","x*y-1"),c("x","y","t"))
#plotfield(c(-4,4,-3,3),c("x*y","x*y-1"),c("x","y","t"))
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as

plotfield <- function(values, lims, d, vars){
  if (length(d)==2){

    Field <- function(t, y_in, parameters){
      dy = numeric(length(y_in))

      for (i in 1:length(dy)){
        for (j in 1:length(dy)){
          assign(vars[j],y_in[j])
        }
        dy[i]=eval(parse(text=parameters[i]))
      }
      return(list(dy))
    }
    Field.flowField <- flowField(Field, xlim=c(lims[1],lims[2]),
                                 ylim = c(lims[3],lims[4]), parameters = d,
                                 points = 20,
                                 add = FALSE, xlab = vars[1], ylab = vars[2] ,
                                 main = "Solution")

    grid()
    lines(values[1,],values[2,], col='red', type = 'l',xlim=c(lims[1],lims[2]),ylim=c(lims[3],lims[4]))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
