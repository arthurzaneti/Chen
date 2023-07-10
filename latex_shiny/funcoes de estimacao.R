#__________________________________PACKAGES_____________________________________
#install.packages(c("progress","flextable", "knitr", "pacman"))
library(pacman)
pacman::p_load(flextable, knitr, progress)

#___________________________BASIC_FUNCTIONS_____________________________________

rchen_rpr <- function(lambda, mds, tau=0.5){
  nelementos <- length(mds)
  vetor_rchen <- (log(1 - log(1-runif(nelementos))* ((1-exp(mds^lambda))/log(1-tau))))^(1/lambda) #função quantílica da Chen reparametrizada em termos do quantil
  return (vetor_rchen)
}

log_vero_chen <- function(theta, tau, covariables, y){
  lambda <- theta[1]
  betas <- theta[2:length(theta)]
  mds <- exp(covariables%*% betas)
  
  log_vero <- suppressWarnings(log(log(1 - tau) / (1 - exp(mds^lambda))) +
                                 log(lambda) + (lambda - 1) * log(y) +
                                 (log(1 - tau) / (1 - exp(mds^lambda))) * (1 - exp(y^lambda)) + (y^lambda))
  return(sum(log_vero))
}

#_______________________________ESTIMATION______________________________________

estim <- function(vetor_random_chen, covariables, tau=0.5, full = F){
  theta_guesses <- c(rep(1, ncol(covariables) +1))
  tryCatch({suppressWarnings(estimacao <- optim(par = theta_guesses, 
                                                fn = log_vero_chen,
                                                tau = tau,
                                                covariables=covariables,
                                                y = vetor_random_chen,
                                                method ="BFGS",
                                                hessian = full, 
                                                control = list(fnscale=-1)))
    if(full){
      return(estimacao)
    }else{
      return(estimacao$par)
    }
  },error=function(e){
    return(NULL)
  })
}

eval_estim <- function(n_rvalues, monte_carlo_iterations, theta, hist=F){
  
  covariables <- matrix(rep(1, n_rvalues))
  for(i in 2:(length(theta) - 1)){
    covariables <- cbind(covariables, matrix(runif(n_rvalues)))
  }
  lambda <- theta[1]
  betas <- theta[2:length(theta)]
  mds <- exp(covariables%*% betas)
  
  estims <- matrix(0, nrow= monte_carlo_iterations, ncol=length(theta))
  pb <- progress_bar$new(format = "Time left: :eta [:bar] :percent", 
                         total = monte_carlo_iterations)
  
  errors <- 0
  for (i in 1:monte_carlo_iterations){
    y <- rchen_rpr(lambda, mds)
    par_hats <- suppressWarnings(estim(y, covariables, tau = 0.5))
    if(is.null(par_hats)){
      errors <- errors + 1
      if(i !=1){
        estims[i,] = estims[i-1,]
      }
    }else{
      estims[i,] = par_hats
    }
    pb$tick()
  }
  eval <- matrix(0, nrow=4, ncol=length(theta))
  eval[1,] <- colMeans(estims)
  eval[2,] <- eval[1,] - theta
  eval[3,] <- apply(estims, 2, sd)
  eval[4,] <- eval[2,]^2 + apply(estims, 2, var)
  
  cols <- c("\u03BB hat")
  for (i in 1:length(betas)){
    cols <- c(cols, paste0("β", i, " hat"))
  }
  colnames(eval) <- cols
  rownames(eval) <- c("Mean", "Bias", "Standard Error", "MSE")
  
  if(hist){
    opar <- par()
    rows <- floor(sqrt(length(theta) * (16/9)))
    cols <- ceiling(length(theta) / rows)
    par(mfrow = c(rows, cols), family = "serif")
    hist(estims[, 1], main = "\u03BB hat", xlab = "", ylab = "")
    for(i in 2:length(theta)){
      hist(estims[, i], main = paste0("β", i-1, " hat"), xlab = "", ylab = "")
    }
    suppressWarnings(par(opar))
  }
  
  return (eval)
}
#_____________________________CONFIDENCE_INTERVALS______________________________

ci <- function(n_rvalues, theta, alpha, monte_carlo = FALSE, monte_carlo_hist = FALSE) {
  covariables <- matrix(rep(1, n_rvalues))
  for(i in 2:(length(theta) - 1)){
    covariables <- cbind(covariables, matrix(runif(n_rvalues)))
  }
  lambda <- theta[1]
  betas <- theta[2:length(theta)]
  mds <- exp(covariables%*% betas)
  
  y <- rchen_rpr(lambda, mds)
  
  tryCatch({
    estimation <- estim(y, covariables, full = TRUE)
    
    inf <- solve(-estimation$hessian)
    par <- estimation$par
    
    lb <- par - qnorm(1 - alpha/2) * sqrt(diag(inf))
    ub <- par + qnorm(1 - alpha/2) * sqrt(diag(inf))
    
    if (monte_carlo_hist) {
      results <- list(lb[1] <= theta[1] && theta[1] <= ub[1])
      for (i in 2:length(theta)){
        results <- append(results,(lb[i] <= theta[i] && theta[i] <= ub[i]))
      }
      results <- append(results, lb)
      results <- append(results, ub)
      return(results)
    }
    else if (monte_carlo) {
      results <- c(lb[1] <= theta[1] && theta[1] <= ub[1])
      for (i in 2:length(theta)){
        results <- append(results,(lb[i] <= theta[i] && theta[i] <= ub[i]))
      }
      return(results)
    } 
    else {
      return(matrix(c(lb,ub), ncol=2))
    }
  }, error = function(e) {
    return(NULL)
  })
}



eval_ci <- function(n_rvalues, monte_carlo_iterations, theta, alpha, hist = F) {
  covariables <- matrix(rep(1, n_rvalues))
  for(i in 2:(length(theta) - 1)){
    covariables <- cbind(covariables, matrix(runif(n_rvalues)))
  }
  lambda <- theta[1]
  betas <- theta[2:length(theta)]
  mds <- exp(covariables%*% betas)
  
  TF_table <- matrix(0, nrow= monte_carlo_iterations, ncol=length(theta))
  errors <- 0
  
  if (!hist) {
    pb <- progress_bar$new(format = "Time left: :eta [:bar] :percent", total = monte_carlo_iterations)
    
    for (i in 1:monte_carlo_iterations) {
      results <- ci(n_rvalues, theta, alpha, monte_carlo = T)
      
      if (is.null(results)) {
        if (i != 1) {
          errors <- errors + 1
          TF_table[i-1,] <- TF_table[i-1,]
        }
      } else {
        TF_table[i,] <- results
      }
      pb$tick()
    }
    
  } else {
    lbs <- matrix(0, monte_carlo_iterations, length(theta))
    ubs <- matrix(0, monte_carlo_iterations, length(theta))
    errors <- 0
    pb <- progress_bar$new(format = "Time left: :eta [:bar] :percent", total = monte_carlo_iterations)
    
    for (i in 1:monte_carlo_iterations) {
      results <- ci(n_rvalues, theta, alpha, monte_carlo_hist = T)
      
      if (is.null(results)) {
        if (i != 1) {
          TF_table[i,] <- TF_table[i-1,]
          lbs[i,] <- lbs[i-1,]
          ubs[i,] <- ubs[i-1,]
        }
      } else {
        s <- length(theta)
        TF_table[i,] <- unlist(results[1:s])
        lbs[i, ] <- unlist(results[(s+1): (2*s)])
        ubs[i, ] <- unlist(results[(2*s +1): length(results)])
      }
      pb$tick()
    }
    
    opar <- par()
    rows <- floor(sqrt(length(theta) * (16/9)))
    cols <- ceiling(length(theta) / rows)
    par(mfrow = c(rows, cols), family = "serif")
    
    hist(lbs[, 1], main = "Confidence interval distributions for λ hat",
         xlim = c(theta[1] - 0.4, theta[1] + 0.4), col = "green")
    hist(ubs[, 1], col = "red", add = TRUE)
    abline(v=theta[1])
    
    for(i in 2:length(theta)){
      hist(lbs[, i], main = paste0("Confidence interval distributions for β", i-2, " hat"),
           xlim = c(theta[i] - 0.4, theta[i] + 0.4), col = "green")
      hist(ubs[, i], col = "red", add = TRUE)
      abline(v=theta[i])
    }
    suppressWarnings(par(opar))
  }
  
  return_vector <- vector()
  for(i in 1:length(theta)){
    return_vector <- append(return_vector, sum(TF_table[,i] / monte_carlo_iterations))
  }
  return(return_vector)
}