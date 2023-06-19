library(RColorBrewer)
rchen <- function(lambda, delta, nelementos){
  valores_aleatorios_tau <- runif(nelementos)
  valores_chen <- log((1-(log(1-valores_aleatorios_tau))/delta))^(1/lambda)
  return (valores_chen)
}

quantilica <- function(lambda, delta, tau){
  return (log((1-(log(1-tau))/delta))^(1/lambda))
}

calcula_delta <- function(lambda, mu, tau){
  return (log(1-tau)/(1- exp(mu^lambda)))
}

chen_reparametrizada <- function(y, lambda, mu, tau){
  exp1 <- (log(1-tau)/(1-exp(mu^lambda)))*lambda*(y^(lambda-1))
  exp2 <- exp((log(1-tau)/(1-exp(mu^lambda)))* (1-exp(y^lambda)) + y^lambda)
  return (exp1 * exp2)
}


col1 <- brewer.pal(4, "Dark2")
col2 <- brewer.pal(4, "Set1")
opar <- par(no.readonly=TRUE)
par(lwd = 2, font.main=4, cex.main=2)

tau <- 0.5
lambdas <- c(0.5,0.7,0.9,1.1)
mus <- c(2,7,3,8)


curve(chen_reparametrizada(x, lambdas[1], mus[1], tau) ,from=0, to=20, ylim=c(0,0.5), col = col1[1], ann=FALSE)
title(main="Função densidade \n probabilidade reparametrizada", xlab = "y", ylab="Densidade de probabilidade" )
curve(chen_reparametrizada(x, lambdas[2], mus[2], tau), add = TRUE, col = col1[2])
curve(chen_reparametrizada(x, lambdas[3], mus[3], tau), add = TRUE, col = col1[3])
curve(chen_reparametrizada(x, lambdas[4], mus[4], tau), add = TRUE, col = col1[4])

legends <- c(expression(paste(lambda,"=0.5, ",mu,"=2")),
             expression(paste(lambda,"=0.7, ",mu,"=7")),
             expression(paste(lambda,"=0.9, ",mu,"=3")),
             expression(paste(lambda,"=1.1, ",mu,"=8")))
legend(x="topright",fill = col1, legend = legends)

curve(chen_reparametrizada(x, 0.7, 7, 0.7), from=0, to=20, xlab="Densidade de propabilidade", col = col2[1])
curve(chen_reparametrizada(x, 0.7, 7, 0.5),  add =TRUE, col = col2[2])
curve(chen_reparametrizada(x, 0.7, 7, 0.3),  add =TRUE, col = col2[3])
curve(chen_reparametrizada(x, 0.7, 7, 0.1),  add =TRUE, col = col2[4])


#note o eixo Y, ele varia de exemplo pra exemplo
random_chen1 <- rchen(0.5, calcula_delta(0.5, 2, 0.5), 1000)
hist(random_chen1)
random_chen2 <- rchen(1.1, calcula_delta(0.7, 7, 0.5), 1000000)
hist(random_chen2)
random_chen3 <- rchen(1.1, calcula_delta(0.9, 3, 0.5), 1000000)
hist(random_chen3)
random_chen4 <- rchen(1.1, calcula_delta(1.1, 8, 0.5), 1000000)
hist(random_chen4)


negativo_ll_chen <- function(y, parametros){
  n<- length(y)
  lambda <- parametros[1]
  mu <- parametros[2]
  delta<-calcula_delta(lambda, mu, 0.5)
  log_verossimilhanca<- n*log(delta)+n*log(lambda) + sum((lambda-1)*log(y)) + delta*(sum(1-exp(y^lambda)))+ sum(y^lambda)
  return (-log_verossimilhanca)
}
random_chen1 <- rchen(0.5, calcula_delta(0.5, 2, 0.5), 1000)
suppressWarnings(estimacao <- optim(par=c(2,2), fn=negativo_ll_chen, y=random_chen1, method ="BFGS", hessian =TRUE))
estimacao$par
