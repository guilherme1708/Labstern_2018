#Laboratorio de Computacao e Simulacao
# Guilherme Navarro NºUSP: 8943160

#Metropolis
#nucleo <- 'Normal'

rg <- 50057666 #numero do RG
cpf <- 43396344847 #numero do CPF
nusp <- 8943160 #numero USP

a <- rg / (10^(nchar(rg))) #a = 0.RG
b <- cpf / (10^(nchar(cpf))) #b = 0.CPF
c <- 1 + (nusp / (10^(nchar(nusp)))) #c = 1.NUSP

f <-function(x) ((x+abs(sin(a*x+b)))^(c-1))*exp(-((x+abs(sin(a*x+b))))) #funcao f(x)
g <- function(x) (x^(c-1)*exp(-x)/gamma(c)) # função de desidade

MCMC <- function(Nsim){
  X <- rep(NA,Nsim) # inicia a cadeia
  X[1] <- 0.1 # chute inicial
  
  for (i in 2:Nsim){
    xprop <- abs(rnorm(1,2,10)) # cadidato proposto com distr Normal
    numer <- g(xprop) # função densidade aplicada no candidato
    denom <- g(X[i-1]) # função densidade aplicada no termo anterior da cadeia
    R <- numer/denom # Razão de Hastings
    alpha <- min(1,R) # probabilidade de aceitação/ rejeição
    if(runif(1)< alpha){
      X[i] <- xprop
    }else{
      X[i] <- X[i-1]
    }
  }
  return(X[20001:Nsim]) # vetor com numeros gerados de acordo com a distribuição g(x) descartando os 20000 primeiros numeros
}

importance.sampling <- function(i, f, B){  #função importance sampling
  X <- MCMC(B) # valores gerados com MCMC
  f(X)/(X^(c-1)*exp(-X)/gamma(c)) # função auxiliar que cobre a f(x)
}

# inicialização de varíavaeis
med <- vector() # vetor de médias
sd <- vector() # vetor de desvios padrões
resultado <- vector()
erro <- vector()
prec <- 0.005 #precisao
N <- 50 # Número de Simulações
B <- 50000 # Tamanho da Amostra

for(i in 1:N){
  I <- importance.sampling(i,f,B)
  med[i] <- mean(I)
  sd[i] <- sd(I)
  erro[i] <- sd(I)/sqrt(B)
  if (erro[i] < prec){
    resultado[i] <- mean(med[i])
  }
} 
valor_est <- mean(resultado, na.rm = T)
message("Implentando um gerador MCMC com núcleo Normal, temos que o valor da integral é: ", round(valor_est,5))
#--------------------------------------------------------------------#

#Gráficos e tabelas

teste <- MCMC(B)

hist(teste, breaks = 50, main = "Metropolis Kernel Normal", las = 1, fre = F, xlim = c(0,8), xlab = "Valores gerados pelo MCMC", ylab = "Densidade", col="grey")
curve(g, col = "blue", add = T)
legend("topright", c('Valores aleatórios',"Densidade Gamma "), lty = 1, col = c(1,4))

df <- data.frame(mean(resultado, na.rm = T), mean(sd), mean(erro))

