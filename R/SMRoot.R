SMRoot<-function(func=NULL,xo=NULL,E=10**-6,maxiter=100,method="steffensen")
{
  if(is.null(func) || !is.function(func))
    stop("Parameter 'func' must be a function.")
  if(is.null(xo) || (!is.numeric(xo) && !is.array(xo)))
    stop("Parameter 'xo' must be a number or array.")
  if(maxiter<1 || !is.numeric(maxiter))
    stop("Parameter 'Maxiter' must be a number greater than 0.")
  if(tolower(method)=="muller")
  {
    if(length(xo)!=3)
      stop("Parameter 'xo' must have 3 values.")

    x3<-xo[3]
    x2<-xo[2]
    x1<-xo[1]
    for(i in c(1:maxiter))
    {
      x12<-((func(x2)-func(x1))/(x2-x1))
      x13<-((func(x3)-func(x1))/(x3-x1))
      x23<-((func(x3)-func(x2))/(x3-x2))
      x123<-(x23-x12)/(x3-x1)
      w<-x12+x13-x23
      q1<-w+sqrt(as.complex((w**2)-(4*(func(x1))*x123)))
      q2<-w-sqrt(as.complex((w**2)-(4*(func(x1))*x123)))
      a<-q2
      if(Mod(q1)>Mod(q2))
        a<-q1
      den<-as.complex(a)
      xo<-as.complex(x1-((2*func(x1))/(den)))
      x3<-x2
      x2<-x1
      x1<-xo
      d<-as.complex((abs(x1-x2)))
      if(Mod(d)< Mod(E))
        break
    }
    if(Mod(d) > Mod(E) && maxiter>=i)
     warning("\n Number of iterations exceeded")
    res<-c(xo,as.integer(Mod(d)),as.integer(Re(i)))
    names(res)<-c("Root","Error","Iterations")
    return (res)
  }
  else if(tolower(method)=="steffensen")
  {
    if(E<=0 || !is.numeric(E))
    {
      stop("Parameter 'E'(error) must be a number greater than 0.")
    }
    for(i in c(1:maxiter))
    {
      x<-xo-((func(xo)**2)/(func(xo+func(xo))-func(xo)))
      d<-abs(x-xo)
      if(is.nan(d))
      {
        d<-Inf
        break
      }

      if(d < E)
        break

      #cat(x,"   ",xo,"\n")
      xo<-x
    }
    if(d > E && maxiter>=i)
      warning("\n  Number of iterations exceeded")

    res<-c(x,d,i)
    names(res)<-c("Root","Error","Iterations")
    return(res)
  }
  else
  {
    stop("Non-existent Method")
  }
}

