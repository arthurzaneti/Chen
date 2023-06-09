rshen <- function(lambda, delta, nelementos){
  valores_aleatorios_tau <- runif(nelementos)
  valores_shen <- log((1-(log(1-valores_aleatorios_tau))/delta))^(1/lambda)
  return (valores_shen)
}

quantilica <- function(lambda, delta, tau){
  return (log((1-(log(1-tau))/delta))^(1/lambda))
}

lambda <- 1;
delta <- 4;
valores_shen <- rshen(lambda, delta, 1000);
hist(valores_shen)
