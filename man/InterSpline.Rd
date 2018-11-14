\name{InterpSpline}
\alias{InterpSpline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
InterpSpline
}
\description{
calcula la interpolación en 3D utilizando el metodo de Spline Cubico
}
\usage{
InterpSpline(x,y,n, xlim = NULL, ylim = NULL, zlim = NULL,
                      xlab = NULL, ylab = NULL, zlab = NULL, add = FALSE, aspect = !add,
                      forceClipregion = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,y}{
  son vectores de los puntos por donde deve pasar la interpolación
}
  \item{n}{
  La longitud de la interpolación deseada.
}
  \item{xlim,ylim,zlim}{
  Límites X, Y y Z. Si está presente, la trama se recorta a esta     región.
}
  \item{xlab,ylab,zlab}{
Títulos para los ejes. nótese bien Estas deben ser cadenas de caracteres; No se aceptan expresiones. Los números serán obligados a cadenas de caracteres.
}
  \item{add}{
  si se desea añadir puntos a una trama existente
  }

  \item{aspect}{
  ya sea una indicación lógica de si se debe ajustar la relación de   aspecto, o una nueva relación.
  }

  \item{forceClippregion}{
  forzar la utilización de una región de recorte, se den o no límites.
  }

  \item{...}{
  Parámetros de material adicionales para pasar a surface3d y decorate3d
  }
}
\value{
retorna una matriz de la función generada.
}
\author{
Javier Marin, Andres Mariño
}
\examples{
x=rnorm(10)
y=rnorm(10)
z=InterpSpline(x,y,10,col="green",xlab="X")
print(z)
z[3,2]
{
  }
}
