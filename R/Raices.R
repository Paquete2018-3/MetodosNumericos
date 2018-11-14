#' Punto fijo
#' @description Para ecuaciones no lineales donde el punto fijo enviado por el usuario existe en el dominio de la derivada de la funcion, aproximando este valor hasta llegar a una raiz con un error hasta de 10^-6
#' @param funcion Ecuacion no lineal definida anteriormente como: nombre_funcion <- function(x)(Ecuacion no lineal en terminos de x)
#' @param x0 Punto enviado por el usuario el cual pertenece a la segunda derivada de la funcion donde posiblemente se encuentre la raiz
#'
#' @return Una tabla con las aproximaciones de la raiz de la ecuacion no lineal, con los diferentes errores desde 10^-1 a 10^-6
#' @export
#'
puntofijo =function(funcion, x0){
  print("------------TABLA DE RAICES------------")
  print("X  | VALOR      | ERROR")
  error=0.1
  contador=1
  while(error >(1/1000000)){

    x1 = funcion(x0)
    error = abs(x1 - x0)
    x0 = x1
    cat("   x", contador, " | ", x1, " |",error,"\n")
    contador = contador+1

  }
}
#' Posicion falsa
#' @description Para una ecuacion no lineal que cumpla la condicion de f(x0)*f(x1)<0 recordando que este metodo es un derivado del metodo de secante.
#' @param funcion Ecuacion no lineal definida anteriormente como: nombre_funcion <- function(x)(Ecuacion no lineal en terminos de x)
#' @param x0 Limite inferior
#' @param x1 Limite superior
#'
#' @return Una tabla con las aproximaciones de la raiz de la ecuacion no lineal, con los diferentes errores desde 10^-1 a 10^-6
#' @export
#'
posicionFalsa <- function(funcion,x0,x1) {
  contador=1
  x<-seq(x0,x1,0.1)
  error<-0.1
  print("------------TABLA DE RAICES------------")
  print("X  | VALOR      | ERROR")
  while (error > (1/1000000)) {
    x<-(funcion(x1)*x0-funcion(x0)*x1)/(funcion(x1)-funcion(x0))
    if (funcion(x) == 0) break
    if (funcion(x)*funcion(x0) < 0) {x1 <- x}
    else {x0 <- x}
    error<-abs(funcion(x)-funcion(x1-x0))
    cat("   x", contador, " | ", x, " |",error,"\n")
    contador=contador+1

  }
}
#' Secante
#' @description Para ecuciones no lineales donde la funcion es doblemente diferenciable, es decir que su segunda derivada existe y es continua, ademas debe tener una raiz unica para generar la raiz.
#' @param funcion Ecuacion no lineal definida anteriormente como: nombre_funcion <- function(x)(Ecuacion no lineal en terminos de x)
#' @param x0 Limite inferior
#' @param x1 Limite superior
#'
#' @return Una tabla con las aproximaciones de la raiz de la ecuacion no lineal, con los diferentes errores desde 10^-1 a 10^-6
#' @export
#'
secante <- function(funcion,x0,x1) {
  x<-(funcion(x1)*x0-funcion(x0)*x1)/(funcion(x1)-funcion(x0))
  error <-0.1
  contador=1
  print("------------TABLA DE RAICES------------")
  print("X  | VALOR      | ERROR")
  while (error > 1.e-64) {
    x0<-x1
    x1<-x
    x<-(funcion(x1)*x0-funcion(x0)*x1)/(funcion(x1)-funcion(x0))
    if (funcion(x) == 0) break
    error<-abs(funcion(x0)-funcion(x1))
    cat("   x", contador, " | ", x, " |",error,"\n")
    contador=contador+1
  }
}
#' Raices con metodo de punto fijo, posicion falsa, secante
#' @description Desarrollado para ecuaciones no lineales, por medio de una funcion , y ya sea un punto o dos
#' se generara raices para esa ecuacion, teniendo en cuenta las restricciones de cada funcion descrita anteriormente,
#' debe escribir la ecuacion en forma de funcion y demas parametros y se calculara evaluando las restricciones
#' de los metodos para esa ecuacion.
#' @param funcion Ecuacion no lineal definida anteriormente como: nombre_funcion <- function(x)(Ecuacion no lineal en terminos de x)
#' @param x0 Limite menor del intervalo para Secante y PosicionFalsa, para PuntoFijo es el valor inicial del usuario
#' @param x1 Limite mayor del intervalo para Secante y PosicionFalsa
#'
#' @return Una tabla con las aproximaciones de la raiz de la ecuacion no lineal, con los diferentes errores desde 10^-1 a 10^-6
#' @export
raizPuFijPosFSec= function (funcion, x0 = NULL, x1 = NULL){
  if(is.null(x1)){
    print("Punto fijo")
    return(puntofijo(funcion,x0))
  }
  else{
    if(funcion(x0)*funcion(x1) < 0){
      print("Posicion falsa")
      return(posicionFalsa(funcion,x0,x1))
    }
    else{
      print("secante")
      return(secante(funcion,x0,x1))

    }



  }

}



