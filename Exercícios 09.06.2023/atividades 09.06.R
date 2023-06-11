library("RColorBrewer")

rchen <- function(lambda, delta, nelementos){
  valores_aleatorios_tau <- runif(nelementos)
  valores_chen <- log((1-(log(1-valores_aleatorios_tau))/delta))^(1/lambda)
  return (valores_chen)
}

quantilica <- function(lambda, delta, tau){
  return (log((1-(log(1-tau))/delta))^(1/lambda))
}

calcula_omega <- function(lambda, mu, tau){
  return (log(1-tau)/(1- exp(mu^lambda)))
}

chen_reparametrizada <- function(y, lambda, mu, tau){
  exp1 <- (log(1-tau)/(1-exp(mu^lambda)))*lambda*(y^(lambda-1))
  exp2 <- exp((log(1-tau)/(1-exp(mu^lambda)))* (1-exp(y^lambda)) + y^lambda)
  return (exp1 * exp2)
}

tau <- 0.5

dark2 <- brewer.pal(4, "Dark2")
curve(chen_reparametrizada(x, 0.5, 2, tau) ,from=0, to=20, xlab = "y", ylab = "Densidade de probabilidade", col = dark2[1])
curve(chen_reparametrizada(x, 0.7, 7, tau) ,from=0, to=20, xlab = "y", add = TRUE, col = dark2[2])
curve(chen_reparametrizada(x, 0.9, 3, tau) ,from=0, to=20, xlab = "y", add = TRUE, col = dark2[3])
curve(chen_reparametrizada(x, 1.1, 8, tau) ,from=0, to=20, xlab = "y", add = TRUE, col = dark2[4])

set1 <- brewer.pal(4, "Set1")
curve(chen_reparametrizada(x, 0.7, 7, 0.7), from=0, to=20, xlab="Densidade de propabilidade", col = set1[1])
curve(chen_reparametrizada(x, 0.7, 7, 0.5), from=0, to=20, add =TRUE, col = set1[2])
curve(chen_reparametrizada(x, 0.7, 7, 0.3), from=0, to=20, add =TRUE, col = set1[3])
curve(chen_reparametrizada(x, 0.7, 7, 0.1), from=0, to=20, add =TRUE, col = set1[4])

#note o eixo Y, ele varia de exemplo pra exemplo
random_chen1 <- rchen(0.5, calcula_omega(0.5, 2, 0.5), 1000000)
hist(random_chen1)
random_chen2 <- rchen(1.1, calcula_omega(0.7, 7, 0.5), 1000000)
hist(random_chen2)
random_chen3 <- rchen(1.1, calcula_omega(0.9, 3, 0.5), 1000000)
hist(random_chen3)
random_chen4 <- rchen(1.1, calcula_omega(1.1, 8, 0.5), 1000000)
hist(random_chen4)
