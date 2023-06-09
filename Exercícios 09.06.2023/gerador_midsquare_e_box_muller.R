library(tidyverse)

gerador_midsquare <- function(n_valores){
  random <- sample(0:10000, 1)
  print(glue("Semente: {random}"))  
  vec <- c(random/10000)
  for (i in 1:(n_valores-1)){
    random_squared <- random^2
    random_squared_formatado <- sprintf("%08d", random_squared)
    random_squared_4digitos <- as.integer(substr(random_squared_formatado, 4,7))
    vec <- append(vec, random_squared_4digitos/10000)
    random <- random_squared_4digitos
  }
  return (vec)
}

valores_aleatorios <- gerador_midsquare(1000)
hist(valores_aleatorios)

gerador_normal_box_muller <- function(n_valores){
  n_pares <- ceiling(n_valores)
  vec<- c()
  for (i in 1:(n_pares-1)){
    random_uniforme <- runif(2)
    vec <- append(vec, sqrt(-2*log(random_uniforme[1]))*cos(2*pi*random_uniforme[2]))
    if (length(vec)==n_valores){
      break
    }
    vec <- append(vec, sqrt(-2*log(random_uniforme[1]))*sin(2*pi*random_uniforme[2]))
  }
  return (vec)
}
rnormal <- gerador_normal_box_muller(1000)
hist(rnormal)

