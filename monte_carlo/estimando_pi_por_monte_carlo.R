library(glue)
raio <- 1
iteracoes <- 10000
gera_aleatoria <- function(raio){
  return(runif(1, -1, 1))
}

pontos_dentro_circulo <- 0
for (i in 1:iteracoes){
  abcissa <- gera_aleatoria(raio)
  coordenada <- gera_aleatoria(raio)
  if (abcissa**2 + coordenada**2 < 1){
    pontos_dentro_circulo <- pontos_dentro_circulo +1
  }
}
glue("Pi é: {(pontos_dentro_circulo/iteracoes)*(2*raio)**2}")
