library("glue")

fatorial <- function(x){
  if(x <=1){
    return(1)
  }else{
    return (x * fatorial(x-1))
  }
}

normal <- function(x) {
  (1/(sqrt(2*pi)))*exp(-x^2/2)
}

poisson <- function(x){
  3**x * exp(-3)/fatorial(x)
}
numero_estimacoes <- 1000
a_normal <- -20 # escolhido arbitrariamente
b_normal <-  20 # escolhido arbitrariamente

somatorio_normal <- 0 
for (i in 1:numero_estimacoes){
  somatorio_normal <- somatorio_normal + normal(runif(n=1, a_normal, b_normal))
}

area_estimada_normal <- (b_normal-a_normal)*(somatorio_normal/numero_estimacoes)
somatorio_normal <- 0
glue("A área estimada da distribuição normal é: {area_estimada_normal}")

a_poisson <- 0
b_poisson <-  60 # esacolhido arbitrariamente 

somatorio_poisson <- 0 
for (i in 1:numero_estimacoes){
  somatorio_poisson <- somatorio_poisson + poisson(runif(n=1, a_poisson, b_poisson))
}
area_estimada_poisson <- (b_poisson-a_poisson)*(somatorio_poisson/numero_estimacoes)
somatorio_poisson <- 0
glue("A área estimada da distribuição poisson é: {area_estimada_poisson}")

#



