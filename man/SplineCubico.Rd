\name{SplineCubico}
\alias{SplineCubico}

\title{
Interpolacion Spline Cubico.
}
\description{
Dados unos puntos ( coodenadas x, y ),
se interpola, grafica y presenta una tabla con los datos pertinenetes ( x, spline, error).
}
\usage{
SplineCubico(x,y)
}
\arguments{
  \item{x}{
Vector de puntos en el eje x.
}
  \item{y}{
Vector de puntos en el eje y.
}
}
\details{
La longitud de los vectores debe ser igual.
}
\value{
No tiene retorno, pero imprime la grafica con la interpolacion,
la tabla con los valores de x, el spline y el error.
}
\references{
Barrera Sanchez , P., Hernandez Mederos , V., & Duran Sanchez, C. (1996). El ABC de los splines. Mexico: Sociedad Matematica Mexicana.

Pantoja Carhuavilca, H. (s.f.). Spline Cúbicos. Obtenido de http://hermes22.yolasite.com/resources/Archivos_Antiguos/Spline_cubico.pdf.

}
\author{
Camilo Hoyos.

Catalina Morales.
}
\examples{
x<-seq(1:50)
y<-exp(x)
splineCubico(x,y)

x<-c(0,1,2,2.5,3,4)
y<-c(1.4,0.6,1.0,0.65,0.6,1.0)
splineCubico(x,y)

}

