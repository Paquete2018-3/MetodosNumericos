\name{SMRoot}
\alias{SMRoot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{SMRoot}
\description{Obtener las raices de la función dada}
\usage{
SMRoot(func = NULL, xo = NULL, E = 10^-6, maxiter = 100, method = "steffensen")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{func}{
 Función matemática dada
}
  \item{xo}{
Valor inicial dado
}
  \item{E}{
Tolerancia
}
  \item{maxiter}{
 Número máximo de iteraciones
}
  \item{method}{
Steffenson o Muller
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
Lista de resultados con: Valor de la raíz, Error and número de iteraciones necesitado
  \item{Root }{Valor de la raíz}
  \item{Error }{Error del método númerico}
  \item{Iterations }{número de iteraciones necesitado}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Maria Fernanda Garces, Juan Sebastian Leon
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
