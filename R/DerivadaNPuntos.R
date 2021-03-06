#' @title Halla el error del método de dos puntos
#' @description Halla el error de truncamiento al calcular el valor de la derivada
#' de una expresion en un punto dado teniendo en cuenta la magnitud de la variación en x que determina
#' la precisión del cálculo (h)
#' @param h la magnitud de la variación en x
#' @param expresion expresion a la que se le calculará la derivada
#' @param punto punto en el que se evaluará la derivada
#' @return error de truncamiento del cálculo de la derivada con el método de dos puntos
#' @export hallarErrorDosPuntos
#' @examples
#' expresion=expression((x^2)+x)
#' h=1*10^(-5)
#' punto=0
#' hallarErrorDosPuntos(h,expresion,punto)
hallarErrorDosPuntos<-function(h,expresion,punto){
  segundaDerivada<-D(D(expresion,"x"),"x")
  x=punto
  error<-round((-h/2)*eval(segundaDerivada),10)
  return(error)
}

#' @title Halla el error del método de tres puntos
#' @description Halla el error de truncamiento al calcular el valor de la derivada
#' de una expresion en un punto dado teniendo en cuenta la magnitud de la variación en x que determina
#' la precisión del cálculo (h)
#' @param h la magnitud de la variación en x
#' @param expresion expresion a la que se le calculará la derivada
#' @param punto punto en el que se evaluará la derivada
#' @return error de truncamiento del cálculo de la derivada con el método de tres puntos
#' @export hallarErrorTresPuntos
#' @examples
#' expresion=expression((x^2)+x)
#' h=1*10^(-5)
#' punto=0
#' hallarErrorTresPuntos(h,expresion,punto)
hallarErrorTresPuntos<-function(h,expresion, punto){
  terceraDerivada<-D(D(D(expresion,"x"),"x"),"x")
  x=punto

  error<-round(((h^2)/3)*eval(terceraDerivada),digits = 10)
  return(error)
}

#' @title Halla el error del método de cinco puntos
#' @description Halla el error de truncamiento al calcular el valor de la derivada
#' de una expresion en un punto dado teniendo en cuenta la magnitud de la variación en x
#' que determina la precisión del cálculo (h)
#' @param h la magnitud de la variación en x
#' @param expresion expresion a la que se le calculará la derivada
#' @param punto punto en el que se evaluará la derivada
#' @return error de truncamiento del cálculo de la derivada con el método de cinco puntos
#' @export hallarErrorCincoPuntos
#' @examples
#' expresion=expression((x^2)+x)
#' h=1*10^(-5)
#' punto=0
#' hallarErrorCincoPuntos(h,expresion,punto)
hallarErrorCincoPuntos<-function(h,expresion, punto){
  quintaDerivada<-D(D(D(D(D(expresion,"x"),"x"),"x"),"x"),"x")
  x=punto
  error<-round(((h^4)/30)*eval(quintaDerivada),10)
  return(error)
}

#' @title Halla la variación de x del método de dos puntos
#' @description Halla la variación de x necesaria (h), para que se obtenga el error dado
#' como parametro al calcular el valor de la derivada de una expresion en un punto dado.
#' @param error la magnitud del error que determinará el valor de h
#' @param expresion expresion a la que se le calculará la derivada
#' @param punto punto en el que se evaluará la derivada
#' @return la variación de x (h) para el que la derivada tiene el error de truncamiento dado
#' @export hallarHDosPuntos
#' @examples
#' expresion=expression((x^2)+x)
#' error=1*10^(-5)
#' punto=0
#' hallarHDosPuntos(error,expresion,punto)
hallarHDosPuntos<-function(error,expresion, punto){
  segundaDerivada<-D(D(expresion,"x"),"x")
  x=punto
  termino<-round(eval(segundaDerivada),10)
  if(termino==0){
    termino=1*10^(-10)
  }
  h<-round((-2*error)/termino,10)
  return(h)
}

#' @title Halla la variación de x del método de tres puntos
#' @description Halla la variación de x necesaria (h), para que se obtenga el error dado
#' como parametro al calcular el valor de la derivada de una expresion en un punto dado.
#' @param error la magnitud del error que determinará el valor de h
#' @param expresion expresion a la que se le calculará la derivada
#' @param punto punto en el que se evaluará la derivada
#' @return la variación de x (h) para el que la derivada tiene el error de truncamiento dado
#' @export hallarHTresPuntos
#' @examples
#' expresion=expression((x^2)+x)
#' error=1*10^(-5)
#' punto=0
#' hallarHTresPuntos(error,expresion,punto)
hallarHTresPuntos<-function(error,expresion, punto){
  terceraDerivada<-D(D(D(expresion,"x"),"x"),"x")
  x=punto
  termino<-round(eval(terceraDerivada),10)
  if(termino==0){
    termino=1*10^(-10)
  }
  h<-round(sqrt(3*error/termino),digits = 10)
  return(h)
}

#' @title Halla la variación de x del método de cinco puntos
#' @description Halla la variación de x necesaria (h), para que se obtenga el error dado
#' como parametro al calcular el valor de la derivada de una expresion en un punto dado.
#' @param error la magnitud del error que determinará el valor de h
#' @param expresion expresion a la que se le calculará la derivada
#' @param punto punto en el que se evaluará la derivada
#' @return la variación de x (h) para el que la derivada tiene el error de truncamiento dado
#' @export hallarHCincoPuntos
#' @examples
#' expresion=expression((x^2)+x)
#' error=1*10^(-5)
#' punto=0
#' hallarHCincoPuntos(error,expresion,punto)
hallarHCincoPuntos<-function(error,expresion, punto){
  quintaDerivada<-D(D(D(D(D(expresion,"x"),"x"),"x"),"x"),"x")
  x=punto
  termino<-round(eval(quintaDerivada),10)
  if(termino==0){
    termino=1*10^(-10)
  }
  h<-round(sqrt(sqrt(30*error/termino)),10)
  return(h)
}

#' @title Halla el valor de la derivada en un punto
#' @description Halla el valor de la derivada por el método de dos puntos de una expresion en
#' un punto dado, con la restricción de un error o un h específico(según sea el caso).
#' @param expresion expresion a la que se le calculará la derivada
#' @param punto punto en el que se evaluará la derivada
#' @param dato puede representar el valor de la variación x(h) o el error con el que se quiere
#' calcular el valor de la derivada.
#' @param tipo especifíca si el anterior parametro(dato) corresponde a cualquiera de los dos datos
#' que puede representar: TRUE=el dato representa el error, FALSE=el dato representa h.
#' @return el valor de la derivada de la expresion en el punto dado.
#' @export dosPuntos
#' @examples
#' expresion=expression((x^2)+x)
#' punto=0
#' dato=1*10^(-5)
#' dosPuntos(expresion,punto,dato,FALSE)
dosPuntos<-function(expresion, punto,dato, tipo){
  vector<-c(0,0,0)
  if(tipo == TRUE){
    #QUIERE DECIR QUE EL PARAMETRO RECIBIDO ES EL ERROR POR LO TANTO TOCA CALCULAR H
    h<-round(hallarHDosPuntos(dato,expresion,punto),10)
    if(h==0){
      h<-1*10^(-10)
    }
    vector[2]=h
    vector[1]=dato
  }else{

    #EL DATO RECIBIDO ES H POR LO TANTO SOLO ES NECESARIO ADIGNAR EL VALOR
    h<-round(dato,10)
    if(h==0){
      h<-1*10^(-10)
    }
    error<-hallarErrorDosPuntos(h,expresion,punto)
    if(error==0){
      error<-1*10^(-10)
    }
    vector[1]=h
    vector[2]=error

  }
  #---------------------------------------
  x<-punto+h
  primerTermino<-eval(expresion)
  #---------------------------------------
  x<-punto
  segundoTermino<-eval(expresion)
  #---------------------------------------
  resultado<-(primerTermino-segundoTermino)/h
  vector[3]=resultado


  #cat(Error)
  return(vector)
}

#' @title Halla el valor de la derivada en un punto
#' @description Halla el valor de la derivada por el método de tres puntos de una expresion en
#' un punto dado, con la restricción de un error o un h específico(según sea el caso).
#' @param expresion expresion a la que se le calculará la derivada
#' @param punto punto en el que se evaluará la derivada
#' @param dato puede representar el valor de la variación x(h) o el error con el que se quiere
#' calcular el valor de la derivada.
#' @param tipo especifíca si el anterior parametro(dato) corresponde a cualquiera de los dos datos
#' que puede representar: TRUE=el dato representa el error, FALSE=el dato representa h.
#' @return el valor de la derivada de la expresion en el punto dado.
#' @export tresPuntos
#' @examples
#' expresion=expression((x^2)+x)
#' punto=0
#' dato=1*10^(-5)
#' tresPuntos(expresion,punto,dato,FALSE)
tresPuntos<-function(expresion, punto,dato,tipo){
  vector<-c(0,0,0)
  if(tipo == TRUE){
    #QUIERE DECIR QUE EL PARAMETRO RECIBIDO ES EL ERROR POR LO TANTO TOCA CALCULAR H
    h<-round(hallarHTresPuntos(dato,expresion,punto),10)
    if(h==0){
      h<-1*10^(-10)
    }
    vector[2]=h
    vector[1]=dato
  }else{
    #EL DATO RECIBIDO ES H POR LO TANTO SOLO ES NECESARIO ADIGNAR EL VALOR
    h<-round(dato,10)
    if(h==0){
      h<-1*10^(-10)
    }
    error<-hallarErrorTresPuntos(h,expresion,punto)
    if(error==0){
      error<-1*10^(-10)
    }
    vector[1]=h
    vector[2]=error
  }

  #--------------------------------------
  x<-punto
  primerTermino<-round(3*eval(expresion),digits = 10)
  #--------------------------------------
  x<-punto+h
  segundoTermino<-round(4*eval(expresion),digits = 10)
  #--------------------------------------
  x<-punto+2*h
  tercerTermino<-round(eval(expresion),digits = 10)
  resultado<-round((-primerTermino+segundoTermino-tercerTermino)/(2*h),digits = 10)
  vector[3]=resultado
  return(vector)
}
#' @title Halla el valor de la derivada en un punto
#' @description Halla el valor de la derivada por el método de cinco puntos de una expresion en
#' un punto dado, con la restricción de un error o un h específico(según sea el caso).
#' @param expresion expresion a la que se le calculará la derivada
#' @param punto punto en el que se evaluará la derivada
#' @param dato puede representar el valor de la variación x(h) o el error con el que se quiere
#' calcular el valor de la derivada.
#' @param tipo especifíca si el anterior parametro(dato) corresponde a cualquiera de los dos datos
#' que puede representar: TRUE=el dato representa el error, FALSE=el dato representa h.
#' @return el valor de la derivada de la expresion en el punto dado.
#' @export cincoPuntos
#' @examples
#' expresion=expression((x^2)+x)
#' punto=0
#' dato=1*10^(-5)
#' cincoPuntos(expresion,punto,dato,FALSE)
cincoPuntos<-function(funcion, punto,dato,tipo){
  vector<-c(0,0,0)
  if(tipo == TRUE){
    #QUIERE DECIR QUE EL PARAMETRO RECIBIDO ES EL ERROR POR LO TANTO TOCA CALCULAR H
    h<-round(hallarHCincoPuntos(dato,expresion,punto),10)
    if(h==0){
      h<-1*10^(-10)
    }
    vector[2]=h
    vector[1]=dato
  }else{
    #EL DATO RECIBIDO ES H POR LO TANTO SOLO ES NECESARIO ADIGNAR EL VALOR
    h<-round(dato,10)

    if(h==0){
      h<-1*10^(-10)
    }
    error<-hallarErrorCincoPuntos(h,expresion,punto)

    vector[1]=h
    vector[2]=error
  }
  #--------------------------------------------
  x<-punto-2*h
  primerTermino<-round(eval(expresion),10)
  #--------------------------------------------
  x<-punto-h
  segundoTermino<-round(8*eval(expresion),10)
  #--------------------------------------------
  x<-punto+h
  tercerTermino<-round(8*eval(expresion),10)
  #--------------------------------------------
  x<-punto+2*h
  cuartoTermino<-round(eval(expresion),10)
  resultado<-round((1/(12*h))*(primerTermino-segundoTermino+tercerTermino-cuartoTermino),10)
  vector[3]=resultado
  return(vector)

}

#'@title Tabla de derivadas por dos puntos
#'@description Imprime una tabla de la derivada evaluada f'(a) a partir del error o de h
#'@param h La magnitud de la variación en x
#'@param funcion Expresion a la que se le calculará la derivada
#'@param punto Punto en el que se evaluará la derivada
#'@param dato Aqui se contiene h o el error dependiendo del parametro tipo
#'@param tipo Este es un valor booleano que sirve para saber a partir de que calcular la derivada:(True es a partir del error)(False a partir de h)
#'@param cantidadIteraciones Cuantas veces dividiremos entre dos el valor del dato. Con mas iteraciones, el dato que se encuentra en la ultima tupla es mas exacto, o mas cercano a la respuesta verdadera.
#'@return No retorna nada, imprime la matriz del error, el h y el valor de la funcion derivada en un punto
#'@export tablaDosPuntos
#'@examples
#'expresion=expression(x^5)
#'punto=2
#'dato=90
#'tipo=FALSE
#'iteraciones=20
#'tablaDosPuntos(expresion,punto,dato,tipo,iteraciones)
#''expresion=expression(x^5)
#'punto=2
#'dato=90
#'tipo=TRUE
#'iteraciones=20
#'tablaDosPuntos(expresion,punto,dato,tipo,iteraciones)
#'
tablaDosPuntos<-function(funcion, punto,dato,tipo,cantidadIteraciones){
  cont<-1
  vectorGeneral<-c()
  while(cont<=cantidadIteraciones*3){
    vectorGeneral<-c(vectorGeneral,0)
    cont<-cont+1

  }

  cont<-0

#cat("Error------------h-----------f'(x)------------")
    while(cont<cantidadIteraciones){
      vect<-dosPuntos(funcion, punto,dato,tipo)
      #print(vect)
      vectorGeneral[(cont*3)+1]=vect[1]
      vectorGeneral[(cont*3)+2]=vect[2]
      vectorGeneral[(cont*3)+3]=vect[3]
      #print(vectorGeneral[cont])
      #print(vectorGeneral[cont+1])
      #print(vectorGeneral[cont+2])
      dato=dato/2
      cont<-cont+1
    }
    #print(vectorGeneral)
    #matriz<-matrix(vectorGeneral,cantidadIteraciones)
    matriz = matrix(
         vectorGeneral, # the data elements
        nrow=cantidadIteraciones,              # number of rows
       ncol=3,              # number of columns
         byrow = TRUE)
    #print(matriz)
    if(tipo==TRUE){
      colnames(matriz) <- c("Error--------","H------------","f'(x)-------------")
    }else{
      colnames(matriz) <- c("H------------","Error--------","f'(x)-------------")
    }

    vectorEspacios<-c()
    cont<-1
    while (cont<=cantidadIteraciones) {
      vectorEspacios<-c(vectorEspacios,"")
      cont<-cont+1
    }
    rownames(matriz)<-vectorEspacios
    print(matriz)


}

#'@title Tabla de derivadas por tres puntos
#'@description Imprime una tabla de la derivada evaluada f'(a) a partir del error o de h y la formula de derivada por tres puntos
#'@param h La magnitud de la variacion en x
#'@param funcion Expresion a la que se le calculará la derivada
#'@param punto Punto en el que se evaluará la derivada
#'@param dato Aqui se contiene h o el error dependiendo del parametro tipo
#'@param tipo Este es un valor booleano que sirve para saber a partir de que calcular la derivada:(True es a partir del error)(False a partir de h)
#'@param cantidadIteraciones Cuantas veces dividiremos entre dos el valor del dato. Con mas iteraciones, el dato que se encuentra en la ultima tupla es mas exacto, o mas cercano a la respuesta verdadera.
#'@return No retorna nada, imprime la matriz del error, el h y el valor de la funcion derivada en un punto
#'@export tablaTresPuntos
#'@examples
#'expresion=expression(x^5)
#'punto=2
#'dato=90
#'tipo=FALSE
#'iteraciones=20
#'tablaTresPuntos(expresion,punto,dato,tipo,iteraciones)
#'expresion=expression(x^5)
#'punto=2
#'dato=90
#'tipo=TRUE
#'iteraciones=20
#'tablaTresPuntos(expresion,punto,dato,tipo,iteraciones)
#'
tablaTresPuntos<-function(funcion, punto,dato,tipo,cantidadIteraciones){
  cont<-1
  vectorGeneral<-c()
  while(cont<=cantidadIteraciones*3){
    vectorGeneral<-c(vectorGeneral,0)
    cont<-cont+1

  }

  cont<-0

  #cat("Error------------h-----------f'(x)------------")
  while(cont<cantidadIteraciones){
    vect<-tresPuntos(funcion, punto,dato,tipo)
    #print(vect)
    vectorGeneral[(cont*3)+1]=vect[1]
    vectorGeneral[(cont*3)+2]=vect[2]
    vectorGeneral[(cont*3)+3]=vect[3]
    #print(vectorGeneral[cont])
    #print(vectorGeneral[cont+1])
    #print(vectorGeneral[cont+2])
    dato=dato/2
    cont<-cont+1
  }
  #print(vectorGeneral)
  #matriz<-matrix(vectorGeneral,cantidadIteraciones)
  matriz = matrix(
    vectorGeneral, # the data elements
    nrow=cantidadIteraciones,              # number of rows
    ncol=3,              # number of columns
    byrow = TRUE)
  #print(matriz)
  if(tipo==TRUE){
    colnames(matriz) <- c("Error--------","H------------","f'(x)-------------")
  }else{
    colnames(matriz) <- c("H------------","Error--------","f'(x)-------------")
  }

  vectorEspacios<-c()
  cont<-1
  while (cont<=cantidadIteraciones) {
    vectorEspacios<-c(vectorEspacios,"")
    cont<-cont+1
  }
  rownames(matriz)<-vectorEspacios
  print(matriz)
}

#'@title Tabla de derivadas por cinco puntos
#'@description Imprime una tabla de la derivada evaluada f'(a) a partir del error o de h y la formula de derivada por cinco puntos
#'@param h La magnitud de la variacion en x
#'@param funcion Expresion a la que se le calculará la derivada
#'@param punto Punto en el que se evaluará la derivada
#'@param dato Aqui se contiene h o el error dependiendo del parametro tipo
#'@param tipo Este es un valor booleano que sirve para saber a partir de que calcular la derivada:(True es a partir del error)(False a partir de h)
#'@param cantidadIteraciones Cuantas veces dividiremos entre dos el valor del dato. Con mas iteraciones, el dato que se encuentra en la ultima tupla es mas exacto, o mas cercano a la respuesta verdadera.
#'@return No retorna nada, imprime la matriz del error, el h y el valor de la funcion derivada en un punto
#'@export tablaCincoPuntos
#'@examples
#'expresion=expression(x^5)
#'punto=2
#'dato=90
#'tipo=FALSE
#'iteraciones=20
#'tablaCincoPuntos(expresion,punto,dato,tipo,iteraciones)
#'punto=2
#'dato=90
#'tipo=TRUE
#'iteraciones=20
#'tablaCincoPuntos(expresion,punto,dato,tipo,iteraciones)
tablaCincoPuntos<-function(funcion, punto,dato,tipo,cantidadIteraciones){
  cont<-1
  vectorGeneral<-c()
  while(cont<=cantidadIteraciones*3){
    vectorGeneral<-c(vectorGeneral,0)
    cont<-cont+1
    #print(cont)
  }
  #print(vectorGeneral)
  cont<-0

  #cat("Error------------h-----------f'(x)------------")
  while(cont<cantidadIteraciones){
    vect<-cincoPuntos(funcion, punto,dato,tipo)
    #print(vect)
    vectorGeneral[(cont*3)+1]=vect[1]
    vectorGeneral[(cont*3)+2]=vect[2]
    vectorGeneral[(cont*3)+3]=vect[3]
    #print(vectorGeneral[cont])
    #print(vectorGeneral[cont+1])
    #print(vectorGeneral[cont+2])
    dato=dato/2
    cont<-cont+1
  }
  #print(vectorGeneral)
  #matriz<-matrix(vectorGeneral,cantidadIteraciones)
  matriz = matrix(
    vectorGeneral, # the data elements
    nrow=cantidadIteraciones,              # number of rows
    ncol=3,              # number of columns
    byrow = TRUE)
  #print(matriz)
  if(tipo==TRUE){
    colnames(matriz) <- c("Error--------","H------------","f'(x)-------------")
  }else{
    colnames(matriz) <- c(" H------------","Error--------","f'(x)-------------")
  }

  vectorEspacios<-c()
  cont<-1
  while (cont<=cantidadIteraciones) {
    vectorEspacios<-c(vectorEspacios,"")
    cont<-cont+1
  }
  rownames(matriz)<-vectorEspacios
  print(matriz)
}
