#' moment function generator
#' @param type the type of distribution. It take a value from "Poisson", "power law", "Gamma", "Beta", mixture of "2Poissons", and "2Exponential"
#' @description If we need to estimate n parameters, we can use Moment Generating Function to get from E[x] to E[X^n] in order to obtain n equations
#' The value meanb is our estimate in the currently iteration of the mean of the beta distribution,
#' while x is our data. The gmm() function, in calling our g(), will calculate the average of m1 (m2, m3) over
#' all our data, and then in setting that average to 0, we are equating our parametric estimate of the
#' population mean to our sample mean, exactly what MM is supposed to do.
#'
g<-function(type){
  if(type == "Beta"){
    g1 = function(th, x){
      t1 <- th [1]
      t2 <- th [2]
      t12 <- t1 + t2
      meanb <- t1 / t12
      m1 <- meanb - x
      m2 <- t1*t2 / ( t12^2 * ( t12 +1)) - ( x - meanb)^2
      f <- cbind (m1,m2)
      return ( f )
    }
    return (g1)
  }
  else if(type == "Poisson"){
    g1 = function(th, x){
      lambda <- th [1]
      meanb <- lambda
      m1 <- meanb - x
      m2 <- lambda - (x-meanb)^2
      f <- cbind(m1)
      return (f)
    }
    return (g1)
  }
  else if (type == "Gamma"){
    g1 = function(th, x){
      alpha <- th [1]
      beta <- th [2]
      meanb <- alpha/beta
      m1 <- meanb - x
      m2 <- alpha/(beta^2) - (x-meanb)^2
      f <- cbind(m1, m2)
      return ( f )
    }
    return (g1)
  }
  else if (type == "2Poisson"){
    g1 = function(th, x){
      lambda1 <- th [1]
      lambda2 <- th [2]
      p <- th[3]
      meanb <- p*lambda1 +(1-p)*lambda2
      m1 <- meanb - x
      m2 <- p*((lambda1^2)+lambda1)+(1-p)*((lambda2^2)+lambda2)-(x^2)
      m3 <-  (meanb+3*p*(lambda1^2)+p*(lambda1^3)+3*(1-p)*(lambda2^2)+(1-p)*(lambda2^3)) - (x^3)
      f <- cbind(m1,m2,m3)
      return ( f )
    }
    return (g1)
  }

  else if (type == "2Exponential") {
    g1 = function(th, x){
      lambda1 <- th [1]
      lambda2 <- th [2]
      p <- th[3]
      meanb <- p*(lambda1^(-1)) +(1-p)*(lambda2^(-1))
      m1 <- meanb - x
      m2 <- 2*p*(lambda1^(-2))+2*(1-p)*(lambda2^(-2)) - (x^2)
      m3 <- (6*p*(lambda1^(-3))+6*(1-p)*(lambda2^(-3))) - (x^3)
      f <- cbind(m1,m2,m3)
      return ( f )
    }
    return (g1)
  }
  else { #  power law
    g1 = function(th, x){
      alpha <- th [1]
      meanb <- (alpha-1)/(alpha-2)
      m1 <- meanb-x
      f <-cbind(m1)
      return ( f )
    }
    return (g1)
  }
}
