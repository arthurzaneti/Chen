# 10 primeiras potências de 3

l <- list()
for (i in 0:9){
  l <- append(l, 3**i)
}
l

# obtendo maior e menor valores de uma amostra gerada (distribuição normal)

vetor_aleatorio <- rnorm(100, mean= pi)
min(vetor_aleatorio)
max(vetor_aleatorio)

# calculando fatorial com loop

number_sorted <- sample.int(50,1)
fatorial <- 1
for (i in 2:number_sorted){
  fatorial <- fatorial*i
}
fatorial
number_sorted

# calculando fatorial com função recursiva, mais elegante ;)

factorial_func = function(n){
  if (n==1 || n==0){
    return (1)
  }else{
    return (factorial_func(n-1) * n)
  }
}
factorial_func(6)

# programas para calcular numericamente a probabilidade de ver 3 em um dado
# --------1--------
cont<-1
for(i in 1: 10000){
  face<-sample(6,1)
  if (3==face) cont<-cont+1
}
print(cont/10000)
# --------2--------

cont<-1
new_vector<-vector()
for (j in 1:200){
  cont<-0
  for(i in 1: 100000){
    face<-round(6*runif(1))
    if (6 ==face) cont<-cont+1
  }
  average <- cont/100000
  new_vector<- append(new_vector, average - 0.16666)
  print(average)
  remove(average)
}
hist(new_vector)
# A primeira é melhor, isso acontece porque, quando usamos o método de gerar um número aleatório entre 0 e 1 e multiplica-lo por 6 temos menos probabilidade de tirar 6 que os outros númeors


# aproximação de e usando for

x <- as.integer(readline(prompt=">>> "))
enax <-0
for (i in 0:9){
  enax <- enax + (x**i/factorial_func(i))
}
enax

# aproximação de e usando função

x <- as.integer(readline(prompt="Enter x: "))
n <- as.integer(readline(prompt="Enter the number of terms you want to approximate e**x: "))
aprox_vec <- vector()
aprox_e <- function(x, n){
  enax <-0
  for (i in 0:(n-1)){
    enax <- enax + (x**i/factorial_func(i))
    aprox_vec <- append(aprox_vec,exp(x)-enax)
  }
  plot(aprox_vec)
  return (enax)
}
valor <- aprox_e(x,n)

# calculando a função de máxima verossimilhança da distribuição poisson

calcula_verossimilhança_Poission = function(lambda, vetor){
  return (mean(vetor))
}
