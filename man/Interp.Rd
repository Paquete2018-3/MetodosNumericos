\name{interpSplineC}
\alias{interpSplineC}
\title{
interpSplineC(x, y=NULL, z, xo=seq(min(x), max(x), length = nx),
       yo=seq(min(y), max(y), length = ny),
       linear = TRUE, extrap=FALSE, duplicate = "error", dupfun = NULL,
       nx = 40, ny = 40,
       jitter = 10^-12, jitter.iter = 6, jitter.random = FALSE)
}
\description{
Estas funciones implementan la interpolación bivariada en una grilla.
  para datos de entrada espaciados irregularmente. Estriado bilineal o bicúbico
  La interpolación se aplica utilizando diferentes versiones de algoritmos de
  Akima.
  se ha modificado de la función original
}
\usage{
interpSplineC(x, y=NULL, z, xo=seq(min(x), max(x), length = nx),
       yo=seq(min(y), max(y), length = ny),
       linear = TRUE, extrap=FALSE, duplicate = "error", dupfun = NULL,
       nx = 40, ny = 40,
       jitter = 10^-12, jitter.iter = 6, jitter.random = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    Vector de las coordenadas en x de los puntos o
    \code{SpatialPointsDataFrame} objeto.
    Valores faltantes no son aceptados.
  }
  \item{y}{
    Vector de las coordenadas en y de los puntos.
    Valores faltantes no son aceptados.

    Si se deja como NULL indica que \code{x} deberia ser un
  \code{SpatialPointsDataFrame} y \code{z} nombra la variable
  de interés es el banco de datos.
  }
  \item{z}{
    Vector de las coordenadas en z de los puntos o un carácter variable
    que nombra la variable de interés en
    \code{SpatialPointsDataFrame} \code{x}.

    Valores faltantes no son aceptados.

    \code{x}, \code{y}, y \code{z} deben tener la misma longitud
    (excepto si \code{x} es un \code{SpatialPointsDataFrame}) y puede contener
    no menos de 4 puntos. Los puntos de \code{x} y
    \code{y} no deben ser colineales, por ejemplo, no deben caer en la misma
    linea
    (dos vectores \code{x} y \code{y} tales que \code{y = ax + b} para algún
    \code{a}, \code{b} no produciran resultados relevantes). Algunas
    heuristicas se construyen para evitar este caso añadiendo un pequeño de
    jitter de
    \code{x} y \code{y} cuando el numero de valores \code{NA} en el resultado
    excede el 10\%.


    \code{interp} es para casos en los que se tenga n los valores
    \code{x}, \code{y} esparcidos sobre un plano y un valor \code{z}
    para cada uno.  Si, en cambio, estas tratando de evaluar una funcion
    matemática, o tener una interpretación gráfica de relaciones descritas
    por un polinomio, intenta \code{outer()}.
  }
  \item{xo}{
    Vector de coordenadas de x en una grilla de salida. El predeterminado
    es 40 puntos, espaciados regularmente sobre el rango de \code{x}.  si la
    extrapolación no se está usando (\code{extrap=FALSE}, el predeterminado),
    \code{xo} debería tener un rango que es cercano o está dentro del rango de
    \code{x} para resultados significantes.
  }
  \item{yo}{Vector de coordenadas yen la grilla de salida; análogo a
    \code{xo}, ver arriba.}
  \item{linear}{logico -- indicando si es lineal o spline
    interpolacion que debe ser usada.}
  \item{extrap}{
    bandera logica: La extrapolacion debe ser usada afuera de el casco convexo
    determinado por los puntos?}
  \item{duplicate}{caracter string indicando como manejar duplicados en los
  puntos. los valores posibles son:
    \describe{
      \item{\code{"error"}}{Produce mensaje de error,}
      \item{\code{"strip"}}{Elimina los valores z duplicados,}
      \item{ \code{"mean"},\code{"median"},\code{"user"}}{calcula
	mean , median or funcion definida por el usuario (\code{dupfun}) de los
	valores z duplicados.}
  }}
  \item{dupfun}{Una funcion aplicada a los puntos duplicados si,
    \code{duplicate= "user"}.}
  \item{nx}{Dimension de la grilla de salida en direccion de x}
  \item{ny}{Dimension de la grilla de salida en direccion de x}
  \item{jitter}{Jitter de cantidad de \code{diff(range(XX))*jitter} (XX=x
    or y) sera añadido a las coordenadas si se detectan puntos colineares.
    Despues la interpolación se vuelve a intentar.

    Notese que jitter no se genera aleatoriamente a menos que
    \code{jitter.random} se ponga en \code{TRUE}. esto asegura un resultado
    reproducible. \code{\link[tripack]{tri.mesh}} de paquete
    \code{tripack} esa el mismo mecanismo jitter. Eso significa que puedes
    gráficar la triangulación en la cima de la interpolación y ver la misma
    triangulación como la usada para interpolar, ver ejemplos abajo.
  }
  \item{jitter.iter}{número de reintentos de iteración con jitter, cantidad
  que será multiplicada en cada iteración por \code{iter^1.5}}
  \item{jitter.random}{lógico, ver \code{jitter}, predeterminado en
    \code{FALSE}
  }
}

}
\references{
Akima, H. (1978). A Method of Bivariate Interpolation and Smooth Surface Fitting for Irregularly Distributed Data Points. ACM Transactions on Mathematical Software 4, 148-164.

Akima, H. (1996). Algorithm 761: scattered-data surface fitting that has the accuracy of a cubic polynomial. ACM Transactions on Mathematical Software 22, 362–371.

R. J. Renka (1996). Algorithm 751: TRIPACK: a constrained two-dimensional Delaunay triangulation package. ACM Transactions on Mathematical Software. 22, 1-8.
}
\author{
akima modificado por Javier marin, Andres Mariño, traducido, Leonardo Gómez
}

