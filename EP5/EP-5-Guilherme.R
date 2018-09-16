#Laboratorio de Computação e Simulação
# Guilherme Navarro NºUSP: 8943160

#Dados do tempo de falha:
tempos <- c(0.01, 0.19, 0.51, 0.57, 0.70, 0.73, 0.75, 0.75, 1.11, 1.16,
            1.21, 1.22, 1.24, 1.48, 1.54, 1.59, 1.61, 1.61, 1.62, 1.62,
            1.71, 1.75, 1.77, 1.79, 1.88, 1.9, 1.93, 2.01, 2.16, 2.18,
            2.3, 2.3, 2.41, 2.44, 2.57, 2.61, 2.62, 2.72, 2.76, 2.84,
            2.96, 2.98, 3.19, 3.25, 3.31, 1.19, 3.5, 3.5, 3.5, 3.5)

ro <- 0.3 # ro dado pelos chineses

#x[1] = alfa*, x[2] = beta* e x[3] = gamma*

weibull <- function(x,t) ((x[2]*(t+x[1])^(x[2]-1)/(x[3]^x[2]))*exp(-((t+x[1])/x[3])^x[2]))/exp(-(x[1]/x[3])^x[2]) # densidade da distribuição weibull truncada

reability <- function(x,t) (exp(-((t+x[1])/x[3])^x[2]))/exp(-(x[1]/x[3])^x[2]) # função de sobrevivencia

vero_ot <- function(x) { # Função verosimilhança p/ otimização
  prod_w <- prod(weibull(x,tempos[1:45]))
  prod_r <- prod(reability(x,tempos[46:50]))
  vero <- (-1)*(prod_w*prod_r)
  return (vero)
}

vero <- function(x) { # Função verosimilhança p/ mcmc
  prod_w <- prod(weibull(x,tempos[1:45]))
  prod_r <- prod(reability(x,tempos[46:50]))
  vero <- (prod_w*prod_r)
  return (vero)
}

#Processo de otimização 
#install.packages("alabama") # intala o pacote alabama

library(numDeriv)
library(alabama) # Pacote de otimização não linear

restr_beta <- function(x){ # restrição do intervalo de beta, beta tem que estar entre 3 e 4
  h <- 3 <= x[2] &&  x[2] <= 4
  return(h)
}

heq <- function(x){ # Restrição para otimização = 0
  h1 <- ro*x[3]*gamma(1+(1/x[2]))-x[1]
  return(h1)
}

# ro dado pelo fabricante dos leds
otimiza <- constrOptim.nl(par = c(1.1,3,3.1), fn = vero_ot, hin = restr_beta, heq=heq)

#Metropolis-Hastings

theta <- otimiza$par # vetor theta* com os parâmetros otimizados anteiormente
mi <- theta[3]*gamma(1+1/theta[2]) # média da distribuição weibull truncada
var <- theta[3]^2*(gamma(1+2/theta[2]))+(gamma(1+1/theta[2]))^2 # variância da distribuição weibull truncada
ro <- theta[1]/mi # Valor de ro da amostra

#install.packages("MASS")
library(MASS) # Para a normal multivariada 

MCMC <- function(Nsim){
  sample <- array(dim = c(Nsim+1,3))
  m_covar <- diag(3) #matriz das covariancas
  sample[1,] <- c(1.3,3.3,3.5) # Pontos iniciais
  
  for (i in 1:Nsim){
    prop <- sample[i,] + mvrnorm(1, c(0,0,0), m_covar) #Proposta de incremento o chute inicial com valores da normal multivariada.
    
    if((prop[1] > 0) & (prop[2] >= 3) & (prop[2] <= 4) & (prop[3] >= 0)){ # Se a proposta pertence ao espaço paramétrico continuo o processo.
      R <-  vero(prop)/vero(sample[i,]) # Razão de Hastings
      if (runif(1) <  R & !is.na(R)){ # probabilidade de aceitação/ rejeição
        sample[i+1,] <- prop
      }else{
        sample[i+1,] <- sample[i,]
      }
    }else{
      sample[i+1,] <- sample[i,]
    }
  }
  return(sample[4000:Nsim,])
}

# Calculo da probabilidade a posteriori 

Nsim <- 10000 # Tamanho da amostra total
amostras<- MCMC(Nsim) # Amostra estabilizada
t <- vector() # vetor T*

for (i in 1:length(amostras[,1])){ # conjunto T*
  if (vero(amostras[i,]) >= vero(theta)){
    t[i] <- 1
  }else{
    t[i] <- 0 
  }
} 

contador <- 0

for (i in 1:length(t)){
  if (t[i] == 1){
    contador <- contador +1
  }
}

k <- contador/length(t) # a credibilidade é proporção de do conjunto T*

ev_h <- 1-k #cálculo da evidência
message("Valor aproximado do e-value: ", round(ev_h,3))

