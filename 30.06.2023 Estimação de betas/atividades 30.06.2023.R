# 0) Para tabelas no formato LaTeX usar o pacote knirt
# Exemplo:
#   print(kable(mresults, format="latex", digits=3))
# O arquivo mresults é uma covariables e a função usada é a kable
# 1) Considerar a estrutura de regressão para mu, dada na eq 5 do TCC da Giovana.
# Para essa estrutura gerar números pseudo-aleatórios (y), considerando que o mu é diferente para cada t.
# Você pode criar uma matriz para a estrutura dos x´s.
# X<-cbind(rep(1, n), runif(n)) # essa matriz deve ser fixa, ou seja, ser calculada fora do loop de Monte Carlo. 
# beta<-(beta0, beta1), e assim 
# mu<-exp(X%*%beta) #esse %*% é usado para multiplicação matricial
# Os valores de beta você escolhe, podendo ser positivo, negativo...
# Mas é interessante que forneçam valores médios de mu semelhantes aos exemplos que já tem feito.  Não dá para escolher aleatoriamente os valores do parâmetro
# algum sentido precisam fazer. 
# 2) Ajustar a função de log-verossimilhança e a parte da estimação para que consiga estimar lambda, beta0 e beta1.
# 3) Avaliar o desempenho dos estimadores lambda, beta0 e beta1, usando simulação de Monte Carlo. 
# 4) Considerar outras variações do modelo, com mais covariáveis (colunas em X e consequentemente mais betas) 
# 5) Estudar a parte 3.3 Intervalos de confiança e testes de hipóteses do TCC da Giovana, mais especificamente a de intervalo de confiança, implementar os intervalos 
# para os parâmetros do modelo estudado. 
# 6) Estudar a avaliação numérica de intervalos de confiança (ver 6.4 Avaliação numérica de intervalos de confiança da apostila) e implementar essa avaliação
# - Coloquei um exemplo de programação em que os intervalos de confiança para a distribuição Kumaraswamy foram obtidos, a ideia é a mesma. 
# Usar o optim e a hessiana, que é a matriz de segundas derivadas, fazer o negativo dela e inverter, usando a função solve. Esse será o J^{-1}(theta_hat) descrito na seção 3.3 do TCC da Giovana.

#__________________________________PACKAGES_____________________________________
install.packages(c("progress","flextable", "knitr", "pacman"))
library(pacman)
pacman::p_load(flextable, knitr, pacman)

#___________________________BASIC_FUNCTIONS_____________________________________

rchen_rpr <- function(lambda, mds, tau=0.5){
  nelementos <- length(mds)
  romegas <- runif(nelementos) # omega é o quantil da Chen reparametrizada
  vetor_rchen <- (log(1 - log(1-romegas)* ((1-exp(mds^lambda))/log(1-tau))))^(1/lambda) #função quantílica da Chen reparametrizada em termos do quantil
  return (vetor_rchen)
}

log_vero_chen <- function(par, tau, covariables, y){
    lambda <- par[1]
    beta0 <- par[2]
    beta1 <- par[3]
    
    mds <- exp(covariables%*% c(beta0,beta1))
    
    log_vero <- suppressWarnings(log(log(1 - tau) / (1 - exp(mds^lambda))) +
                             log(lambda) + (lambda - 1) * log(y) +
                             (log(1 - tau) / (1 - exp(mds^lambda))) * (1 - exp(y^lambda)) + (y^lambda))
    return(sum(log_vero))
}

#_______________________________ESTIMATION______________________________________

estim <- function(vetor_random_chen, covariables, tau=0.5, full = F){
    parametros_inciais <- c(1,1,1)
    tau <- 0.5
    tryCatch({suppressWarnings(estimacao <- optim(par = parametros_inciais, 
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
  
  covariables <- cbind(rep(1, n_rvalues), runif(n_rvalues))
  
  lambda <- theta[1]
  beta0 <- theta[2]
  beta1 <- theta[3]
  mds <- exp(covariables%*% c(beta0,beta1))
  
  lambda_hats <- rep(0, monte_carlo_iterations)
  beta0_hats <- rep(0, monte_carlo_iterations)
  beta1_hats <- rep(0, monte_carlo_iterations)
  
  pb <- progress_bar$new(format = "Time left: :eta [:bar] :percent", 
                         total = monte_carlo_iterations)
  
  errors <- 0
  for (i in 1:monte_carlo_iterations){
    y <- rchen_rpr(lambda, mds)
    par_hats <- suppressWarnings(estim(y, covariables, tau = 0.5))
    if(is.null(par_hats)){
      errors <- errors + 1
      if(i !=1){
        lambda_hats[i] = lambda_hats[i-1]
        beta0_hats[i] = beta0_hats[i-1]
        beta1_hats[i] = beta1_hats[i-1]
      }
    }else{
      lambda_hats[i] = par_hats[1]
      beta0_hats[i] = par_hats[2]
      beta1_hats[i] = par_hats[3]
    }
    pb$tick()
  }
 eval <- matrix(0, nrow=4, ncol=3)
 eval[1,] <- colMeans(cbind(lambda_hats, beta0_hats, beta1_hats))
 eval[2,] <- eval[1,] - c(lambda, beta0, beta1)
 eval[3,] <- apply(cbind(lambda_hats, beta0_hats, beta1_hats), 2, sd)
 eval[4,] <- eval[2,]^2 + apply(cbind(lambda_hats, beta0_hats, beta1_hats), 2, var)

 colnames(eval) = c("\u03BB hat", "\u03B2\u2080 hat", "\u03B2\u2081 hat")
 rownames(eval) = c("Mean", "Bias", "Standard Error", "MSE")
 
 if(hist){
   opar <- par()
   par(mfrow = c(2,2), cex.main = 2, family = "serif")
   hist(lambda_hats, main = "\u03BB hat", xlab = "", ylab = "")
   ?hist
   hist(beta0_hats, main = "\u03B2\u2080 hat", xlab = "", ylab = "")
   hist(beta1_hats, main = "\u03B2\u2081 hat", xlab = "", ylab = "")
   par(opar)
 }
 
 return (eval)
}
#_____________________________CONFIDENCE_INTERVALS______________________________

ci <- function(n_rvalues, theta, alpha, monte_carlo = F, monte_carlo_hist = F){
  
  covariables <- cbind(rep(1, n_rvalues), runif(n_rvalues))
  
  lambda <- theta[1]
  beta0 <- theta[2]
  beta1 <- theta[3]
  
  mds <- exp(covariables %*% c(beta0, beta1))
  
  y <- rchen_rpr(lambda, mds)
  tryCatch({estimation <- estim(y, covariables, full=T)
    inf <- solve(-estimation$hessian)
    
    lb <- theta - qnorm(1 - alpha/2)*sqrt(diag(inf))
    ub <- theta + qnorm(1 - alpha/2)*sqrt(diag(inf))
    par <- estimation$par
    if(monte_carlo_hist){
      return(list(lambda_in_ic = (lb[1] <= par[1] && par[1] <= ub[1]), 
                  beta0_in_ic = (lb[2] <= par[2] && par[2]<= ub[2]), 
                  beta1_in_ic = (lb[3] <= par[3] && par[3]<= ub[3]),
                  lb = lb, ub = ub))
    }
    if(monte_carlo){
      return(c(lb[1] <= par[1] && par[1] <= ub[1], 
               lb[2] <= par[2] && par[2]<= ub[2], 
               lb[3] <= par[3] && par[3]<= ub[3]))
    }else{
      return(par)
    }
    }, error = function(e){ return(NULL)})
}

eval_ci <- function(n_rvalues, monte_carlo_iterations, theta, alpha, hist = F) {
  covariables <- cbind(rep(1, n_rvalues), runif(n_rvalues))
  mds <- exp(covariables %*% c(theta[2], theta[3]))
  lambda_in_ic <- rep(F, monte_carlo_iterations)
  beta0_in_ic <- rep(F, monte_carlo_iterations)
  beta1_in_ic <- rep(F, monte_carlo_iterations)
  errors <- 0
  
  if (!hist) {
    pb <- progress_bar$new(format = "Time left: :eta [:bar] :percent", total = monte_carlo_iterations)
    
    for (i in 1:monte_carlo_iterations) {
      results <- ci(n_rvalues, theta, alpha, monte_carlo = T)
      
      if (is.null(results)) {
        if (i != 1) {
          errors <- errors + 1
          lambda_in_ic[i] <- lambda_in_ic[i - 1]
          beta0_in_ic[i] <- beta0_in_ic[i - 1]
          beta1_in_ic[i] <- beta1_in_ic[i - 1]
        }
      } else {
        lambda_in_ic[i] <- results[1]
        beta0_in_ic[i] <- results[2]
        beta1_in_ic[i] <- results[3]
      }
      pb$tick()
    }
    
  } else {
    lbs <- matrix(0, monte_carlo_iterations, 3)
    ubs <- matrix(0, monte_carlo_iterations, 3)
    errors <- 0
    pb <- progress_bar$new(format = "Time left: :eta [:bar] :percent", total = monte_carlo_iterations)
    
    for (i in 1:monte_carlo_iterations) {
      results <- ci(n_rvalues, theta, alpha, monte_carlo_hist = T)
      
      if (is.null(results)) {
        if (i != 1) {
          lambda_in_ic[i] <- lambda_in_ic[i - 1]
          beta0_in_ic[i] <- beta0_in_ic[i - 1]
          beta1_in_ic[i] <- beta1_in_ic[i - 1]
          lbs[i, ] <- lbs[i - 1, ]
          ubs[i, ] <- ubs[i - 1, ]
        }
      } else {
        lambda_in_ic[i] <- results$lambda_in_ic
        beta0_in_ic[i] <- results$beta0_in_ic
        beta1_in_ic[i] <- results$beta1_in_ic
        lbs[i, ] <- results$lb
        ubs[i, ] <- results$ub
      }
      pb$tick()
    }
    
    opar <- par()
    par(mfrow = c(3, 1), cex.main = 2, family = "serif")
    hist(lbs[, 1], main = "Confidence interval distributions for λ hat",
         xlim = c(theta[1] - theta[1] / 5, theta[1] + theta[1] / 5), col = "green")
    hist(ubs[, 1], col = "red", add = TRUE)
    
    hist(lbs[, 2], main = "Confidence interval distributions for β₀ hat",
         xlim = c(theta[2] - theta[2] / 5, theta[2] + theta[2] / 5), col = "green")
    hist(ubs[, 2], col = "red", add = TRUE)
    
    hist(lbs[, 3], main = "Confidence interval distributions for β₁ hat",
         xlim = c(theta[3] - theta[3] / 5, theta[3] + theta[3] / 5), col = "green")
    hist(ubs[, 3], col = "red", add = TRUE)
  }
  
  return(c(sum(lambda_in_ic, na.rm = TRUE) / monte_carlo_iterations, 
           sum(beta0_in_ic, na.rm = TRUE) / monte_carlo_iterations, 
           sum(beta1_in_ic, na.rm = TRUE) / monte_carlo_iterations))
}

#_________________________________PRINTING______________________________________
print_as_kable <- function(eval, latex = F){
  if(latex){
    kable(eval, format="latex")
  }else{
    kable(eval)
  }
}

print_as_flextable <- function(eval){
  flextable(as.data.frame(eval))
}
#____________________________________MAIN_______________________________________

#eval_estim(n_rvalues, monte_carlo_iterations, theta(lambda, beta0 and beta1), hist = False) is the formatting for eval_estim
#print as flextable or kable

print_as_kable(eval_estim(50, 500, c(1.1, 2, 1)), latex=T)
print_as_kable(eval_estim(50, 5000, c(1.1, 2, 1)))
print_as_flextable(eval_estim(50, 50000, c(1.1, 2, 1)))

print_as_kable(eval_estim(50, 50, c(0.7, 1.5, 2)))
print_as_kable(eval_estim(50, 5000, c(0.7, 1.5, 2)))
print_as_kable(eval_estim(50, 50000,  c(0.7, 1.5, 2)))

print(ci1 <- eval_ci(50, 500, c(1.1, 2, 1), 0.05, hist=T))
ci2 <- eval_ci(50, 5000, c(1.1, 2, 1), 0.05)
ci3 <- eval_ci(50, 50000, c(1.1, 2, 1), 0.05)
