chen_reparametrizada <- function(y, lambda, mu, tau){
  exp1 <- (log(1-tau)/(1-exp(mu^lambda)))*lambda*(y^(lambda-1))
  exp2 <- exp((log(1-tau)/(1-exp(mu^lambda)))* (1-exp(y^lambda)) + y^lambda)
  return (exp1 * exp2)
}
curve(chen_reparametrizada(x, lambdas[1], mus[1], tau) ,from=0, to=20, ylim=c(0,0.5), col = col1[1], ann=FALSE)
