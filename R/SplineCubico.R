interpSplineC <- function(x, y=NULL, z,
                   xo = seq(min(x), max(x), length = nx),
                   yo = seq(min(y), max(y), length = ny), linear = TRUE,
                   extrap = FALSE, duplicate = "error", dupfun = NULL,
                   nx=40, ny=40,
                   jitter = 10^-12, jitter.iter = 5, jitter.random = FALSE)
{
  is.sp <- FALSE
  sp.coord <- NULL
  sp.z <- NULL
  if(!(all(is.finite(x)) && all(is.finite(y)) && all(is.finite(z))))
    stop("missing values and Infs not allowed")

  drx <- diff(range(x))
  dry <- diff(range(y))
  if(drx == 0 || dry == 0)
    stop("all data collinear")    # other cases caught in Fortran code
  if(drx/dry > 10000 || drx/dry < 0.0001)
    stop("scales of x and y are too dissimilar")
  n <- length(x)
  nx <- length(xo)
  ny <- length(yo)
  if(length(y) != n || length(z) != n)
    stop("Lengths of x, y, and z do not match")

  xy <- paste(x, y, sep = ",")# trick for 'duplicated' (x,y)-pairs
  miss <- !extrap             # if not extrapolating, set missing values

  ans <- .Fortran("sdsf3p",
                  md = as.integer(1),
                  ndp = as.integer(n),
                  xd = as.double(x),
                  yd = as.double(y),
                  zd = as.double(z),
                  nx = as.integer(nx),
                  x = as.double(xo),
                  ny=as.integer(ny),
                  y = as.double(yo),
                  z = as.double(matrix(0,nx,ny)),
                  ier = integer(1),
                  wk = double(36 * n),
                  iwk = integer(25 * n),
                  extrap = as.logical(matrix(extrap,nx,ny)),
                  near = integer(n),
                  nxt = integer(n),
                  dist = double(n),
                  linear = as.logical(linear),
                  PACKAGE = "akima")
  if(miss)
    ans$z[ans$extrap] <- NA

  if(ans$ier==10){
    warning("collinear points, trying to add some jitter to avoid colinearities!")
    jitter.trials <- 1
    success <- FALSE
    while(jitter.trials<jitter.iter & !success){
      if(jitter.random){
        j <- list()
        j[[1]] <- rep(c(-1,0,1),length.out=length(x))
        j[[2]] <- rep(c(0,1,-1),length.out=length(x))
        j[[3]] <- rep(c(1,-1,0),length.out=length(x))
        jx <- sample(1:3,1)
        jy <- sample(1:3,1)
        xj <- x+j[[jx]]*diff(range(x))*jitter*jitter.trials^1.5
        yj <- y+j[[jy]]*diff(range(y))*jitter*jitter.trials^1.5
      } else {
        xj <- x+rep(c(-1,0,1),length.out=length(x))*diff(range(x))*jitter*jitter.trials^1.5
        yj <- y+rep(c(0,1,-1),length.out=length(y))*diff(range(y))*jitter*jitter.trials^1.5
      }
      ans <- .Fortran("sdsf3p",
                      as.integer(1),
                      as.integer(n),
                      xd=as.double(xj),
                      yd=as.double(yj),
                      as.double(z),
                      as.integer(nx),
                      x = as.double(xo),
                      as.integer(ny),
                      y = as.double(yo),
                      z = as.double(matrix(0,nx,ny)),
                      ier = integer(1),
                      double(36 * n),
                      integer(25 * n),
                      extrap = as.logical(matrix(extrap,nx,ny)),
                      near = integer(n),
                      nxt = integer(n),
                      dist = double(n),
                      linear = as.logical(linear),
                      PACKAGE = "akima")
    }
  }
  ret <- list(x=ans$x,y=ans$y,z=matrix(ans$z,nx,ny))

  ret
}
InterpSpline=function(x,y,n, xlim = NULL, ylim = NULL, zlim = NULL,
                      xlab = NULL, ylab = NULL, zlab = NULL, add = FALSE, aspect = !add,
                      forceClipregion = FALSE, ...){
  list.of.packages <- c("akima")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  list.of.packages <- c("rgl")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  library(akima)
  library(rgl)
  h=c()
  b=c()
  u=c()
  v=c()
  z=c()

  for (i in 1:n-1){
    h[i]=x[i+1]-x[i]
    b[i]=6*(y[i+1]-y[i])/h[i]
  }
  u[1]=2*(h[1]+h[2])

  v[1]=b[2]-b[1]
  for (j in 1:n-1){
    if(j!=1){
      u[j]=2*(h[j]+h[j-1])-(h[j-1])^2/u[j-1]
      v[j]=b[j]-b[j-1]-h[j-1]*v[j-1]/u[j-1]
    }
  }
  z[n]=0
  for (i in (n-1):1){
    z[i]=(v[i]-h[i]*z[i+1])/u[i]
  }
  z[0]=0
  spline_interpolated <- interpSplineC(x, y, z,
                                xo=seq(min(x), max(x), length = n),
                                yo=seq(min(y), max(y), length = n),
                                linear = FALSE, extrap = TRUE)

  x.si <- spline_interpolated$x
  y.si <- spline_interpolated$y
  z.si <- spline_interpolated$z

  persp3d(x.si, y.si, z.si,xlim, ylim, zlim,
          xlab, ylab, zlab, add, aspect, forceClipregion = FALSE, ...)
  return(z.si)
}

