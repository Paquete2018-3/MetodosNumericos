% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Raices.R
\name{secante}
\alias{secante}
\title{Secante}
\usage{
secante(funcion, x0, x1)
}
\arguments{
\item{funcion}{Ecuacion no lineal definida anteriormente como: nombre_funcion <- function(x)(Ecuacion no lineal en terminos de x)}

\item{x0}{Limite inferior}

\item{x1}{Limite superior}
}
\value{
Una tabla con las aproximaciones de la raiz de la ecuacion no lineal, con los diferentes errores desde 10^-1 a 10^-6
}
\description{
Para ecuciones no lineales donde la funcion es doblemente diferenciable, es decir que su segunda derivada existe y es continua, ademas debe tener una raiz unica para generar la raiz.
}
