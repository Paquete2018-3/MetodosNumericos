#' @title Runge-Kutta orden 4 (RK4)
#' @description Función que permite solucionar una ecuación diferencial ordinaria de primer orden
#' usando el método Runge-Kutta de cuarto orden, además gráfica y entrega el error de la solución.
#' @param dy Function: f(x, y) en una ecuación diferencial de la forma y'=f(x, y).
#' @param ti Número real: t inicial para la solución.
#' @param tf Número real: t final para la solución.
#' @param y0 Número real: valor inicial, es decir y(ti)=y0.
#' @param h Número real: tamaño del paso, no puede ser menor a 10^-4 y debe permitir hacer al menos tres puntos en el intervalo de la solución.
#' @param graficar Valor lógico: si es verdadero grafíca, sino no (por defecto es verdadero).
#' @param numpendientes Número entero: es la raíz cuadrada del número de pendientes de la gráfica (por defecto es 10).
#' @return
#' t: Vector con los valores en x de la solución.
#' w: Vector con los valores en y de la solución.
#' error: Vector con los errores de truncamiento de la solución.
#' @references Método desarrollado por C. Runge y M. W. Kutta en 1900.
#' @references Análisis Numérico. 10a Ed. Richard L. Burden, J. Douglas Faires y Annette M. Burden. Cengage p. 209.
#' @author Jhonny Parra y Laura Donado
#' @export rungekutta4
#' @examples
#' r<-rungekutta4(function(x, y){x-y}, 0, 2, 1, 0.1)
#' data.frame (x=r$t, y=r$w, "Error truncamiento"=r$error)
rungekutta4<-function(dy, ti, tf, y0, h, graficar=TRUE, numpendientes=10){
  library(deSolve)
  library(pracma)
  solucionODE<-function(dy, t, h, y0){
    params <- list(a=1)
    fn <- function(x, y, params) with(params, list(a*dy(x, y)))
    out <- ode(y0, t, fn, params)
    return (out)
  }

  graficarCampoPendiente<-function(x0, xn, y0, yn, func, numpendientes, metodo){
    vectorfield (func, c(x0, xn), c(y0, yn), n=numpendientes, scale = 0.1, col="gray")
    mtext(side=3, metodo,font=2,cex=2)
  }

  graficarSolucionNumerica<-function (x, y){
    points (x, y, pch=20, col="blue")
    for (i in 2:length(x)){
      segments(x[i-1], y[i-1], x[i], y[i], col="red")
    }
  }

  if (ti>=tf){
    cat ("Error: tf debe ser mayor que ti.\n");
    return()
  }
  if (h>(tf-ti)*0.5){
    cat ("Error: El h ingresado es muy grande para el intervalo.\n");
    return()
  }
  if (h<10^-4){
    cat ("Error: El h ingresado es muy pequeño, no puede ser menor que 10^-5.\n");
    return()
  }
  t<-seq(ti, tf, h)
  y<-c(y0)
  error<-c(0)
  solucion=solucionODE (dy, t, h, y0)
  for(i in 2:length(t)){
    k1=h*dy(t[i-1], y[i-1])
    k2=h*dy(t[i-1]+h/2, y[i-1]+k1*(0.5))
    k3=h*dy( t[i-1]+h/2, y[i-1]+k2*(0.5))
    k4=h*dy( t[i-1]+h, y[i-1]+k3)
    y<-c(y, y[i-1]+1/6*(k1+2*k2+2*k3+k4))
    error<-c(error, abs(y[i]-solucion [i, 2]))
  }
  cat (length(y))
  if (graficar){
    graficarCampoPendiente(min(t), max(t), min(y), max(y), dy, numpendientes, "RK4")
    graficarSolucionNumerica(t, y)
  }
  rta<-list(w=y, t=t, error=error)
}

#' @title Runge-Kutta orden 3 (RK3)
#' @description Función que permite solucionar una ecuación diferencial ordinaria de primer orden
#' usando el método Runge-Kutta de tercer orden, además gráfica y entrega el error de la solución.
#' @param dy Function: f(x, y) en una ecuación diferencial de la forma y'=f(x, y).
#' @param ti Número real: t inicial para la solución.
#' @param tf Número real: t final para la solución.
#' @param y0 Número real: valor inicial, es decir y(ti)=y0.
#' @param h Número real: tamaño del paso, no puede ser menor a 10^-4 y debe permitir hacer al menos tres puntos en el intervalo de la solución.
#' @param graficar Valor lógico: si es verdadero grafíca, sino no (por defecto es verdadero).
#' @param numpendientes Número entero: es la raíz cuadrada del número de pendientes de la gráfica (por defecto es 10).
#' @return
#' t: Vector con los valores en x de la solución.
#' w: Vector con los valores en y de la solución.
#' error: Vector con los errores de truncamiento de la solución.
#' @references Método desarrollado por C. Runge y M. W. Kutta en 1900.
#' @references Análisis Numérico. 10a Ed. Richard L. Burden, J. Douglas Faires y Annette M. Burden. Cengage p. 209.
#' @author Jhonny Parra y Laura Donado
#' @export rungekutta3
#' @examples
#' r2<-rungekutta3(function(x, y){x-y}, 0, 2, 1, 0.1)
#' data.frame (x=r2$t, y=r2$w, "Error truncamiento"=r2$error)
rungekutta3<-function(dy, ti, tf, y0, h, graficar=TRUE, numpendientes=10){
  library(deSolve)
  library(pracma)
  solucionODE<-function(dy, t, h, y0){
    params <- list(a=1)
    fn <- function(x, y, params) with(params, list(a*dy(x, y)))
    out <- ode(y0, t, fn, params)
    return (out)
  }

  graficarCampoPendiente<-function(x0, xn, y0, yn, func, numpendientes, metodo){
    vectorfield (func, c(x0, xn), c(y0, yn), n=numpendientes, scale = 0.1, col="gray")
    mtext(side=3, metodo,font=2,cex=2)
  }

  graficarSolucionNumerica<-function (x, y){
    points (x, y, pch=20, col="blue")
    for (i in 2:length(x)){
      segments(x[i-1], y[i-1], x[i], y[i], col="red")
    }
  }

  if (ti>=tf){
    cat ("Error: tf debe ser mayor que ti.\n");
    return()
  }
  if (h>(tf-ti)*0.5){
    cat ("Error: El h ingresado es muy grande para el intervalo.\n");
    return()
  }
  if (h<10^-4){
    cat ("Error: El h ingresado es muy pequeño, no puede ser menor que 10^-5.\n");
    return()
  }
  t<-seq(ti, tf, h)
  y<-c(y0)
  solucion=solucionODE (dy, t, h, y0)
  error<-c(0)
  for(i in 2:length(t)){
    k1=h*dy( t[i-1], y[i-1])
    k2=h*dy(t[i-1]+h/2, y[i-1]+k1*(0.5))
    k3=h*dy(t[i-1]+h, y[i-1]-k1+2*k2)
    y<-c(y, y[i-1]+1/6*(k1+4*k2+k3))
    error<-c(error, abs(y[i]-solucion [i, 2]))
  }
  if (graficar){
    graficarCampoPendiente(min(t), max(t), min(y), max(y), dy, numpendientes, "RK3")
    graficarSolucionNumerica(t, y)
  }
  rta<-list(w=y, t=t, error=error)
}
