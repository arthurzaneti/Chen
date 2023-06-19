# 1) Estudar a parte de simulação de Monte Carlo
# 2) Avaliar o desempenho dos Estimadores de máxima verossimilhança mu e lambda por meio de simulação de Monte Carlo. Para isso usar:
#   - Média das estimativas (mean(vetor_estimativas))
# - Viés (vies<-mean(vetor_estimativas)-valor_parâmetro_fixado)
# - Erro quadrático médio (EQM<-var(vetor_estimativas)+vies^2 )
# - Assimetria e Curtose - Pode usar o pacote moments do R, tem as funções skewness e kurtosis
# - Fazer gráficos para ver se a distribuição do vetor das estimativas se parece com a distribuição normal.
# 3) Variar o tamanho amostral, sugiro usar os tamanhos n=15, 20, 30, 50, 100, 500. Já deve dar uma ideia do desempenho do estimador nos diferentes n´s.

library(glue)
library(progress)

#_______________________________________________
rchen <- function(lambda, delta, nelementos){
  valores_aleatorios_tau <- runif(nelementos)
  valores_chen <- log((1-(log(1-valores_aleatorios_tau))/delta))^(1/lambda)
  return (valores_chen)
}

calcula_delta <- function(lambda, mu, tau){
  return (log(1-tau)/(1- exp(mu^lambda)))
}

chen_reparametrizada <- function(y, lambda, mu, tau){
  exp1 <- (log(1-tau)/(1-exp(mu^lambda)))*lambda*(y^(lambda-1))
  exp2 <- exp((log(1-tau)/(1-exp(mu^lambda)))* (1-exp(y^lambda)) + y^lambda)
  return (exp1 * exp2)
}

estima_parametros <- function(vetor_random_chen){
  tryCatch({
    suppressWarnings(estimacao <- optim(par=c(2,2), fn=ll_chen, y=vetor_random_chen, method ="BFGS", hessian =TRUE, control = list(fnscale=-1)))
    return(estimacao$par)
    },error = function(e){
      return(c(NULL,NULL))
    })
}

ll_chen <- function(y, parametros){
  n<- length(y)
  lambda <- parametros[1]
  mu <- parametros[2]
  delta<-calcula_delta(lambda, mu, 0.5)
  log_verossimilhanca<- n*log(delta)+n*log(lambda) + sum((lambda-1)*log(y)) + delta*(sum(1-exp(y^lambda)))+ sum(y^lambda)
  return (log_verossimilhanca)
}

print_resultados <- function(vetor_estimacoes_lambda, vetor_estimacoes_mu){
  media_lambda_chapeu <- mean(vetor_estimacoes_lambda)
  media_mu_chapeu <- mean(vetor_estimacoes_mu)
  
  glue("MEDIA:
     lambda chapeu: {media_lambda_chapeu}
     mu chapeu: {media_mu_chapeu}
     
     VÍES:
     lambda chapeu: {media_lambda_chapeu - lambda}
     mu chapeu: {media_mu_chapeu - mu}
     
     ERRO-PADRÃO:
     lambda chapeu: {sd(vetor_estimacoes_lambda)}
     mu chapeu: {sd(vetor_estimacoes_mu)}
     
     ERRO QUADRADO MÉDIO:
     lambda chapeu: {(media_lambda_chapeu - lambda)**2 + var(vetor_estimacoes_lambda)}
     mu chapeu: {(media_mu_chapeu - mu)**2 + var(vetor_estimacoes_mu)}")
}

lambda <- 0.5
mu <- 2
tau <- 0.5

iteracoes_monte_carlo <- 500000
valores_aleatorios_gerados <- 5000

vetor_estimacoes_lambda <- vector()
vetor_estimacoes_mu <- vector()

#__________________MONTE CARLO PESADA (CUIDADO AO RODAR)________________________
pb <- progress_bar$new(total = iteracoes_monte_carlo)
for (i in 1:iteracoes_monte_carlo){
  random_chen <- rchen(lambda, calcula_delta(lambda, mu, tau), valores_aleatorios_gerados)
  estimacoes <- estima_parametros(random_chen)
  vetor_estimacoes_lambda <- append(vetor_estimacoes_lambda, estimacoes[1])
  vetor_estimacoes_mu <- append(vetor_estimacoes_mu, estimacoes[2])
  pb$tick()
}
#_________________________FIM MONTE CARLO PESADA_______________________________________

opar <-par()
par(family="serif", cex = 2)
print_resultados(vetor_estimacoes_lambda, vetor_estimacoes_mu)

# Sem histograma porque com 5 x 10⁵ valores com variancias baixa ele fica estranho
# Caso queira conferir deixei as linhas prontas abaixo, apenas descomentá-las
#hist(vetor_estimacoes_lambda, breaks=1000)
#hist(vetor_estimacoes_mu, breaks=1000)

#A grande simulação pesada foi feita para chegar o mais perto possível do verdadeiro viés, erro-padrao e EQM do estimador

iteracoes_monte_carlo_2 <- 500
valores_aleatorios_gerados_2 <- 5000

vetor_estimacoes_lambda_2 <- vector()
vetor_estimacoes_mu_2 <- vector()

pb <- progress_bar$new(total = iteracoes_monte_carlo_2)
for (i in 1:iteracoes_monte_carlo_2){
  random_chen <- rchen(lambda , calcula_delta(lambda, mu, tau), valores_aleatorios_gerados_2)
  estimacoes <- estima_parametros(random_chen)
  vetor_estimacoes_lambda_2 <- append(vetor_estimacoes_lambda_2, estimacoes[1])
  vetor_estimacoes_mu_2 <- append(vetor_estimacoes_mu_2, estimacoes[2])
  pb$tick()
}

hist(vetor_estimacoes_lambda_2, ann = FALSE)
title(main= expression(paste("Distribuição do estimador de ", lambda)), sub="500 iterações e vetor tamanho 5000")
hist(vetor_estimacoes_mu_2, ann=FALSE)
title(main= expression(paste("Distribuição do estimador de ", mu)), sub= "500 iterações e vetor tamanho 5000")

print_resultados(vetor_estimacoes_lambda_2, vetor_estimacoes_mu_2)




iteracoes_monte_carlo_3 <- 500
valores_aleatorios_gerados_3 <- 1000

vetor_estimacoes_lambda_3 <- vector()
vetor_estimacoes_mu_3 <- vector()

pb <- progress_bar$new(total = iteracoes_monte_carlo_3)
for (i in 1:iteracoes_monte_carlo_3){
  random_chen <- rchen(lambda, calcula_delta(lambda, mu, tau), valores_aleatorios_gerados_3)
  estimacoes <- estima_parametros(random_chen)
  vetor_estimacoes_lambda_3 <- append(vetor_estimacoes_lambda_3, estimacoes[1])
  vetor_estimacoes_mu_3 <- append(vetor_estimacoes_mu_3, estimacoes[2])
  pb$tick()
}

hist(vetor_estimacoes_lambda_3, ann = FALSE)
title(main= expression(paste("Distribuição do estimador de ", lambda)), sub="500 iterações e vetor tamanho 1000")
hist(vetor_estimacoes_mu_3, ann=FALSE)
title(main= expression(paste("Distribuição do estimador de ", mu)), sub= "500 iterações e vetor tamanho 1000")

print_resultados(vetor_estimacoes_lambda_3, vetor_estimacoes_mu_3)

lambda <- 0.7
mu <- 7

iteracoes_monte_carlo_4 <- 500
valores_aleatorios_gerados_4 <- 1000

vetor_estimacoes_lambda_4 <- vector()
vetor_estimacoes_mu_4 <- vector()

pb <- progress_bar$new(total = iteracoes_monte_carlo_4)
for (i in 1:iteracoes_monte_carlo_4){
  random_chen <- rchen(lambda, calcula_delta(lambda, mu, tau), valores_aleatorios_gerados_4)
  estimacoes <- estima_parametros(random_chen)
  vetor_estimacoes_lambda_4 <- append(vetor_estimacoes_lambda_4, estimacoes[1])
  vetor_estimacoes_mu_4 <- append(vetor_estimacoes_mu_4, estimacoes[2])
  pb$tick()
}

hist(vetor_estimacoes_lambda_4, ann = FALSE)
title(main= expression(paste("Distribuição do estimador de ", lambda)), sub="500 iterações e vetor tamanho 1000
      com lambda = 0.7 e mu = 7")
hist(vetor_estimacoes_mu_4, ann=FALSE)
title(main= expression(paste("Distribuição do estimador de ", mu)), sub= "500 iterações e vetor tamanho 1000
      com lambda = 0.7 e mu = 7")

print_resultados(vetor_estimacoes_lambda_4, vetor_estimacoes_mu_4)


