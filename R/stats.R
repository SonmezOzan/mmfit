#' summary.mmfit function and plot.mmfit function
#'
#' @param x a set of sample data
#' @param g moment function
#' @param start initial values to start with
#' @param type distribution type. It take a value from "Poisson", "power law", "Gamma", "Beta", mixture of "2Poissons", and "2Exponential".
#'
#' @description We can use summary.mmfit(x, g, start, type) to obtain estimated distribution parameters.
#' We use plot.mmfit(x, g, start, type) to plot fitting curve and empirical CDF with K-S confidence band.

summary.mmfit = function(x, g, start, type){
  result = mmfit(x, g, start, type)
  if(type == "Beta"){
    cat("alpha:", result$thetahat[1], "\t", "std. error:", result$thetahatses[1], "\n")
    cat("beta:", result$thetahat[2], "\t", "std. error:", result$thetahatses[2], "\n")
  }
  else if(type == "Poisson"){
    cat("lambda:", result$thetahat[1], "\t", "std. error:", result$thetahatses[1], "\n")
  }
  else if(type == "Gamma"){
    cat("alpha:", result$thetahat[1], "\t", "std. error:", result$thetahatses[1], "\n")
    cat("beta:", result$thetahat[2], "\t", "std. error:", result$thetahatses[2], "\n")
  }
  else if(type == "2Poisson"){
    cat("lambda1:", result$thetahat[1], "\t", "std. error:", result$thetahatses[1], "\n")
    cat("lambda2:", result$thetahat[2], "\t", "std. error:", result$thetahatses[2], "\n")
    cat("p:", result$thetahat[3], "\t", "std. error:", result$thetahatses[3], "\n")
  }
  else if(type == "2Exponential"){
    cat("lambda1:", result$thetahat[1], "\t", "std. error:", result$thetahatses[1], "\n")
    cat("lambda2:", result$thetahat[2], "\t", "std. error:", result$thetahatses[2], "\n")
    cat("p:", result$thetahat[3], "\t", "std. error:", result$thetahatses[3], "\n")
  }
  else{
    cat("alpha:", result$thetahat[1], "\t", "std. error:", result$thetahatses[1], "\n")
  }

}

plot.mmfit = function(x, g, start, type){
  result = mmfit(x, g, start, type)
  gridExtra::grid.arrange(result$denscomp, result$cdfband, ncol=2)

}
