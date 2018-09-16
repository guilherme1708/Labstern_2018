#Laboratorio de Computacao e Simulacao
# Guilherme Navarro NºUSP: 8943160

#Importance Sampling 
# Método 2

rg = 50057666 #numero do RG
cpf = 43396344847 #numero do CPF
nusp = 8943160 #numero USP

a = rg / (10^(nchar(rg))) #a = 0.RG
b = cpf / (10^(nchar(cpf))) #b = 0.CPF
c = 1 + (nusp / (10^(nchar(nusp)))) #c = 1.NUSP

f <- function(x) ((x+abs(sin(a*x+b)))^(c-1))*exp(-((x+abs(sin(a*x+b))))) #funcao f(x)

importance.sampling <- function(i, f, B){ # função importance sampling
  x <- rexp(B, 1)
  f(x) / dexp(x,1)
}

# inicialização de varíavaeis
med <- vector()
sd <- vector()
resultado <- vector()
erro <- vector()
valor_est <- 0
prec <- 0.005 #precisao
N <- 10 # Número de Simulações
B <- 10000 # Tamanho da Amostra

for(i in 1:N){
  I <- importance.sampling(i, f, B)
  med[i] <- mean(I)
  sd[i] <- sd(I)
  erro[i] <- sd(I)/sqrt(B)
  if (erro[i] < prec){
    resultado[i] <- mean(med[i])
  }
}
valor_est <- mean(resultado)
message("Usando a Distribuição Exponencial com paramentro lambda = 1, temos que o valor da integral é: ", valor_est)

df <- data.frame( resultado, sd, erro) # para construção da tabela