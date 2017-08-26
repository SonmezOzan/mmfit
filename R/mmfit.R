#' mmfit (method of moment function)
#' @param x a set of sample data
#' @param g moment function
#' @param start initial values to start with
#' @param type distribution type. It take a value from "Poisson", "power law", "Gamma", "Beta", mixture of "2Poissons", and "2Exponential".
#'
#' @examples
#' \dontrun{
#' #############################################
#' #Beta test
#' n = 1000
#' x <- rbeta(n, 0.5, 0.5)
#' summary.mmfit(x, g, c(0.1, 0.1), "Beta")
#' plot.mmfit(x, g, c(0.1, 0.1), "Beta")
#' #############################################
#' #Poisson test
#' n = 1000
#' x <- rpois(n, 4)
#' summary.mmfit(x, g, 1, "Poisson")
#' plot.mmfit(x, g, 1, "Poisson")
#' #############################################
#' # Gamma test
#' n = 1000
#' x <- rgamma(n, 9, rate = 2)
#' summary.mmfit(x, g, c(5, 5), "Gamma")
#' plot.mmfit(x, g, c(5, 5), "Gamma")
#' #############################################
#' # 2Poisson test
#' n = 10000
#' x <- ifelse(runif(n) < 0.3, rpois(n, 4), rpois(n, 5))
#' summary.mmfit(x, g, c(3, 3, 0.1), "2Poisson")
#' plot.mmfit(x, g, c(3, 3, 0.1), "2Poisson")
#' #############################################
#' # 2Exponential test
#' n = 1000
#' x <- ifelse(runif(n) < 0.2, rexp(n, 0.3), rexp(n, 0.9))
#' summary.mmfit(x, g, c(0.1, 0.1, 0.1), "2Exponential")
#' plot.mmfit(x, g, c(0.2, 0.2, 0.2), "2Exponential")
#' #############################################
#' # powerlaw test
#' n = 1000
#' x = poweRlaw::rpldis(n, 1, 2)
#' summary.mmfit(x, g, 1.9, "powerlaw")
#' plot.mmfit(x, g, 1.9, "powerlaw")
#'#############################################
#' }
#' @description mmfit function can estimate the parameters of following distributions, including Poisson, power law, gamma, beta, mixture of 2 Poissons, and mixture of 2 Exponentials.
#' We also do calculation for CDFband in this function
#' For powerlaw distribution, simply keep summing k^(-gamma), k = 1,2,3,... until the sum doesn't change much, and then take c to be the reciprocal of the sum.
#' @return an s3 class object including the estimated distribution parameters, standard error, fitting curves, and empirical CDF with K-S confidence band

mmfit = function(x, g, start, type){

  a = summary(gmm::gmm(g(type), x, start))
  l = list()
  #################################################
  l$thetahat = a$coefficients[,1]
  l$thetahatses = a$coefficients[,2]
  #################################################
  l$denscomp = ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::labs(x="Data", y="Density") +
    ggplot2::labs(title = "Parametric and nonparametric density estimates")

  if(type == "Beta"){
    l$denscomp = l$denscomp + ggplot2::geom_histogram(ggplot2::aes(x = x, y = ..density..),
                                                      binwidth = 0.05,
                                                      color = "black",
                                                      fill = "blue",
                                                      alpha = .2
    )
    l$denscomp = l$denscomp + ggplot2::scale_y_continuous(expand = c(0,0),
                                                          limits=c(0,max(ggplot2::ggplot_build(l$denscomp)$data[[1]]$density)*1.1)) +
      ggplot2::geom_line(ggplot2::aes(x, dbeta(x, l$thetahat[1], l$thetahat[2])), na.rm =TRUE)
  }
  else if(type == "Poisson"){
    l$denscomp = l$denscomp + ggplot2::geom_histogram(ggplot2::aes(x = x, y = ..density..),
                                                      binwidth = 1,
                                                      color = "black",
                                                      fill = "blue",
                                                      alpha = .2
    )
    l$denscomp = l$denscomp + ggplot2::scale_y_continuous(expand = c(0,0),
                                                          limits=c(0,max(ggplot2::ggplot_build(l$denscomp)$data[[1]]$density)*1.1)) +
      ggplot2::geom_line(ggplot2::aes(x, dpois(x, l$thetahat[1])), na.rm =TRUE)
  }
  else if(type == "Gamma"){
    l$denscomp = l$denscomp + ggplot2::geom_histogram(ggplot2::aes(x = x, y = ..density..),
                                                      binwidth = 1,
                                                      color = "black",
                                                      fill = "blue",
                                                      alpha = .2
    )
    l$denscomp = l$denscomp + ggplot2::scale_y_continuous(expand = c(0,0),
                                                          limits=c(0,max(ggplot2::ggplot_build(l$denscomp)$data[[1]]$density)*1.1)) +
      ggplot2::geom_line(ggplot2::aes(x, dgamma(x, l$thetahat[1], rate = l$thetahat[2])), na.rm =TRUE)
  }
  else if(type == "2Poisson"){
    l$denscomp = l$denscomp + ggplot2::geom_histogram(ggplot2::aes(x = x, y = ..density..),
                                                      binwidth = 1,
                                                      color = "black",
                                                      fill = "blue",
                                                      alpha = .2
    )
    l$denscomp = l$denscomp + ggplot2::scale_y_continuous(expand = c(0,0),
                                                          limits=c(0,max(ggplot2::ggplot_build(l$denscomp)$data[[1]]$density)*1.1)) +
      ggplot2::geom_line(ggplot2::aes(x,l$thetahat[3]*dpois(x, l$thetahat[1])+(1-l$thetahat[3])*dpois(x, l$thetahat[2])) , na.rm =TRUE)
  }
  else if(type == "2Exponential"){
    l$denscomp = l$denscomp + ggplot2::geom_histogram(ggplot2::aes(x = x, y = ..density..),
                                                      binwidth = 1,
                                                      color = "black",
                                                      fill = "blue",
                                                      alpha = .2
    )
    l$denscomp = l$denscomp + ggplot2::scale_y_continuous(expand = c(0,0),
                                                          limits=c(0,max(ggplot2::ggplot_build(l$denscomp)$data[[1]]$density)*1.1)) +
      ggplot2::geom_line(ggplot2::aes(x,l$thetahat[3]*dexp(x, l$thetahat[1])+(1-l$thetahat[3])*dexp(x, l$thetahat[2])) , na.rm =TRUE)
  }
  else{
    #calculate c for powerlaw distribution
    alpha <- l$thetahat[1]
    sum <- 0
    k <- 1
    while(TRUE){
      if(k^-alpha <0.001 && k!=1)
        break
      else{
        sum<- sum+(k^-alpha)
        k <- k+1
      }
    }
    c <- 1/sum
    l$denscomp = l$denscomp + ggplot2::geom_histogram(ggplot2::aes(x = x, y = ..density..),
                                                      binwidth = 1,
                                                      color = "black",
                                                      fill = "blue",
                                                      alpha = .2
    )
    l$denscomp = l$denscomp + ggplot2::scale_y_continuous(expand = c(0,0),
                                                          limits=c(0,max(ggplot2::ggplot_build(l$denscomp)$data[[1]]$density)*1.1)) +
      ggplot2::xlim(c(0,50)) +
      ggplot2::geom_line(ggplot2::aes(x, c*poweRlaw::dpldis(x, 1, l$thetahat[1])) , na.rm =TRUE)
  }

  ##################################################
  sv = sort(x)
  y = seq(1/length(x), 1, by = 1/length(x))
  xval = vector()
  yval = vector()
  xval = c(xval, 0)
  yval = c(yval, 0)
  for(i in 1:length(sv)){
    xval = c(xval, sv[i])
    yval = c(yval, yval[length(yval)])
    xval = c(xval, sv[i])
    yval = c(yval, y[i])
  }
  xval = c(xval, 1.01*xval[length(xval)])
  yval = c(yval, 1)
  l$cdfband = ggplot2::ggplot(data.frame(x=xval, y=yval)) +
    ggplot2::theme_bw() +
    ggplot2::labs(x="Data", y="CDF") +
    ggplot2::labs(title = "Empirical CDF and K-S confidence band ") +
    ggplot2::geom_line(mapping=ggplot2::aes(x=xval, y=yval)) +
    ggplot2::geom_ribbon(mapping=ggplot2::aes(x=xval, ymin=yval-1.358*(length(sv)^-0.5),ymax=yval+1.358*(length(sv)^-0.5)), alpha=0.2, fill = "blue")
  ##################################################
  class(l) <- "mmf"
  return(l)
  #################################################
}





























