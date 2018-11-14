as.function <- function(formula) {
  function(x) eval(expr = parse(text=formula), list(x))
}

lagrange_constructor <- function(t, cor) {
  n = length(t)
  pol <- seq(from = 0, to = 0, len = n)
  for(i in (1:n))
  {
    polytemp <- c(1)
    factor = 1
    for(j in (1:n))
    {
      if(j == i)
        next
      size = length(polytemp)
      new <- seq(from = 0, to = 0, len = size)
      for(k in (1:size))
      {
        temp = polytemp[k]
        otro = -t[j]
        new[k] = new[k] + temp*otro
        new[k+1] = temp
      }
      factor = factor*(t[i] - t[j])
      polytemp <- new
    }
    factor = cor[i]/factor
    for(i in (0:n))
    {
      pol[i] = pol[i] + polytemp[i]*factor
    }
  }
  ##Aqu? ya tenemos contruido el polinomio pol (arreglo de coeficientes) as? que ahora creamos
  ## la equaci?n eq:
  eq = ""
  if(pol[1]!=0)
    eq <- paste(pol[1], "", sep = "")
  if(pol[2]<0)
  {
    ##Negativo
    eq <- paste(eq, paste(-pol[2], "x", sep = "*"), sep = " - ")
  }
  else if(pol[2]>0)
  {
    ##Positivo
    eq <- paste(eq, paste(pol[2], "x", sep = "*"), sep = " + ")
  }
  for(i in (3:n))
  {
    if(pol[i]<0)
    {
      ##Negativo
      eq <- paste(eq, paste(-pol[i], paste("x^",i-1 , sep = ""), sep = "*"), sep = " - ")
    }
    else if(pol[i]>0)
    {
      ##Positivo
      eq <- paste(eq, paste(pol[i], paste("x^",i-1 , sep = ""), sep = "*"), sep = " + ")
    }
  }
  if(substr(eq, 1, 1) == " ")
  {
    if(eq[1]=="-")
      eq = paste("-", substr(eq, 4, nchar(eq)), sep = "")
    else
      eq = substr(eq, 4, nchar(eq))
  }
  return(eq)
}

#' @title lagrange_parametric
#' @description Esta funcion devuelve las ecuaciones \code{x(t)} y \code{y(t)} para una curva parametrica.
#' @param x Vector y Vector t Vector. Los tres vectores deben tener el mismo numero de elementos.
#' @details Esta funcion es parte del paquete de analisis numerico creado por los estudiantes de la
#' Pontificia Universidad Javeriana para este curso en el semestre 1830.
#' @examples
#' t = c(0,0.25,0.5,0.75,1)
#' x = c(-1,0,1,0,1)
#' y = c(0,1,0.5,0,-1)
#' formulas = lagrange_parametric(t,x,y)
#' @export
lagrange_parametric <- function(t, x, y){

  x_cor = lagrange_constructor(t, x)
  y_cor = lagrange_constructor(t, y)
  return(list("fx" = as.function(x_cor), "fy" = as.function(y_cor)))
}
#' @title lagrange_evaluator
#' @description Esta funcion devuelve el valor de las ecuaciones \code{x(t)} y \code{y(t)} para una curva parametrica para un paramero t.
#' @param t Vector formulas Objeto.
#' @details Esta funcion es parte del paquete de analisis numerico creado por los estudiantes de la
#' Pontificia Universidad Javeriana para este curso en el semestre 1830.
#' @examples
#' t = c(0,0.25,0.5,0.75,1)
#' x = c(-1,0,1,0,1)
#' y = c(0,1,0.5,0,-1)
#' formulas = lagrange_parametric(t,x,y)
#' lagrange_evaluator(t, formulas)
#' @export
lagrange_evaluator <- function(t, formulas)
{
  x = formulas$fx(t)
  y = formulas$fy(t)
  return(list("x" = x, "y" = y))
}
