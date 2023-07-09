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

#print as flextable or kable
#eval_estim(n_rvalues, monte_carlo_iterations, theta(lambda, beta0,...), hist = False) is the formatting for eval_estim

print_as_kable(eval_estim(50, 500, c(1.1, 2, -2)), latex=T)
print_as_kable(eval_estim(50, 5000, c(0.7, 1, 2)))
print_as_flextable(eval_estim(50, 50000, c(1.1, 2, 1)))


#eval_ci(n_rvalues, monte_carlo_iterations, theta, hist=F) is the formatting for eval_ci

ci1 <- eval_ci(200, 100, c(1.1, 2, 1), 0.05)
ci2 <- eval_ci(200, 500, c(1.1, 2, 1), 0.05, hist=T)
ci3 <- eval_ci(200, 5000, c(1.1, 2, 1), 0.05)

eval_ci_matrix <- matrix(c(ci1, ci2, ci3), nrow=3, ncol=3)
colnames(eval_ci_matrix) <- c("λ hat", "β₀ hat", "β₁ hat")
rownames(eval_ci_matrix) <- c("500 iterations", "5000 iterations", "50000 iterations")
print_as_kable(eval_ci_matrix)


#multiple betas

print_as_kable(eval_estim(500, 500, c(1.1, 2, -2, 1, 1.2, 1, 1), hist=T)) 

print_as_kable(eval_estim(30, 5000, c(0.7, 1.1, 0.5)), latex=T)
print_as_kable(eval_estim(100, 5000, c(1.0, 1, 1.2, -1.5, -2, 0.3)), latex=T)

print_as_kable(eval_estim(50, 500, c(1.1, 2, -2, 1, -1, 0.5, 2, 0.7, 0.8)))

print(eval_ci(100, 500, c(0.7,1,2,0.3), 0.05, hist=T))
print(eval_ci(100, 500, c(0.7,1,1,2,0.3,1), 0.05, hist=T))
print(eval_ci(100, 5000, c(0.7,1,1,1), 0.05))
