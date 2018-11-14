bisec <- function(f , a , b ){
  epsilon = 0.0001

  c<-(a+b)/2;
  while(f(c)!=0 && b-a> epsilon)
  {

    if(f(c)<0)
    {
      a<-c
    }
    else
    {
      b=c
    }
    {
      c<-(a+b)/2;
    }
  }

  if(abs(f(c) ) <= epsilon ){
    return(c)
  }
  else{
    stop('no hay una solucion en ese intervalo')
  }
}

f<-function(x){
  return ( sin(x) )
}

#plot(f,-6,6)
#x <- (bisec(f , -pi/2 , pi/2));
#print(x)
