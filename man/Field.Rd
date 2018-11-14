\name{Field}
\alias{Field}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Field Función
%%  ~~function to do ... ~~
}
\description{
Esta función representa el campo que se va a graficar.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Field(t, y_in, parameters)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
  Variable independiente.
%%     ~~Describe \code{t} here~~
}
  \item{y_in}{
  Variables
%%     ~~Describe \code{y_in} here~~
}
  \item{parameters}{
  Parámetros
%%     ~~Describe \code{parameters} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Retorna una lista que conforma al campo que se va a graficar.
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
function (t, y_in, parameters)
{
    dy = numeric(length(y_in))
    for (i in 1:length(dy)) {
        for (j in 1:length(dy)) {
            assign(vars[j], y_in[j])
        }
        dy[i] = eval(parse(text = parameters[i]))
    }
    return(list(dy))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
