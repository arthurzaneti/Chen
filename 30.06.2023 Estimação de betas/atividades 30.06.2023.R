# 0) Para tabelas no formato LaTeX usar o pacote knirt
# Exemplo:
#   print(kable(mresults, format="latex", digits=3))
# O arquivo mresults é uma matriz e a função usada é a kable
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

library(progress)

chen_reparametrizada <- function(y, lambda, mu, tau){
    p1 <- log(1-tau)/(1-exp(mu^lambda))
    p2 <- lambda*(y^(lambda-1))
    p3 <- exp((log(1-tau)/(1-exp(mu^lambda)))*(1-exp(y^lambda)) +y^lambda)
    return(p1*p2*p3)
}

# exatamente como está no artigo da Giovanna, pensei em simplificar mais mas assim está mais organizado
log_vero_betas <- function(mus, lambda, tau){
    p1 <- log(log(1-tau)) - log(1-exp(mus^lambda))
    p2 <- (lambda-1)*log(y) + log(lambda)
    p3 <- (log(1-tau)*(1-exp(y^lambda))) / (1 - exp(mus^lambda)) + y
    return (p1 - p2 + p3)
}

estima_parametros <- function(vetor_random_chen){
  tryCatch({
    suppressWarnings(estimacao <- optim(par=c(2,2), fn=ll_chen, y=vetor_random_chen, method ="BFGS", hessian =TRUE, control = list(fnscale=-1)))
    return(estimacao$par)
  },error = function(e){
    return(c(NULL,NULL))
  })
}

chen_reparametrizada(7, 0.7, 7, 0.5)

tamanho_matriz <- 100
matriz <- cbind(rep(1, tamanho_matriz), runif(tamanho_matriz))

matriz
beta <- c(1.5, 1.5)
mus <- exp(matriz%*%beta)
View(mus)
