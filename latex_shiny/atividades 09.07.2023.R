source("./09.07.2023 Latex/funcoes de estimacao.R")
#________________________________________________________________________________

theta1 <- c(0.7, 2, -1, 0.8)
n1 <- 30
estimacoes1 <- eval_estim(n1, 5000, theta1)
ci1 <- ci(n1, theta1, 0.05)

limites_inferiores1 = c(ci1[,1])
limites_superiores1 = c(ci1[,2])

taxa_de_cobertura1 = c(eval_ci(100, 5000, theta1, 0.05))
cenario1 = rbind(estimacoes1, limites_inferiores1, limites_superiores1, taxa_de_cobertura1)
rownames(cenario1) <- c("Média", "Viés", "Erro padrão", "EQM", "Limites inferiores", "Limites superiores", "Taxa de cobertura")
print_as_kable(cenario1, latex=T)

#________________________________________________________________________________

theta2 <- c(1.1, 0.4, 2, 1.1, 1.6, -3)
n2 <- 100
estimacoes2 <- eval_estim(n2, 5000, theta2)
ci2 <- ci(n2, theta2, 0.01)

limites_inferiores2 = c(ci2[,1])
limites_superiores2 = c(ci2[,2])

taxa_de_cobertura2 = c(eval_ci(n2, 5000, theta2, 0.05))
cenario2 = rbind(estimacoes2, limites_inferiores2, limites_superiores2, taxa_de_cobertura2)
cenario2
rownames(cenario2) <- c("Média", "Viés", "Erro padrão", "EQM", "Limites inferiores", "Limites superiores", "Taxa de cobertura")
print_as_kable(cenario2, latex=T)
