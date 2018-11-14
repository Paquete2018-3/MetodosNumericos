sec <- function(f , x1 , x2){
  res <- ((x2-x1) / (f(x2)-f(x1)))

  return ( res )

}


f <- function(x){
  return (as.double(sin(x)));

}

newton <- function(f , ant , xi){
  it = 0;
  epsilon = 0.0001;
  h <- f(xi) * sec(f ,ant , xi);

  while(it < 500 ){

    h <- f(xi) * sec(f , ant , xi);
    ant <- xi;
    xi <- xi - h;
    #print(sec(f , xi , ant));

    if( abs (f(xi) ) < epsilon){
      return (xi);
    }
    it <- it+1;

  }
  print("no converge");

}


#plot(f,-6,6)
#x <- (newton(f , 0.5 , 3));
#print(x)
