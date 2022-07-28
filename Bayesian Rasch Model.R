###===--- Precision ---===###

### Limpar tudo (Tome cuidado para n?o apagar nada importante!)====
rm(list=ls()) # Excluir todos os objetos do ambiente de trabalho
gc()          # Liberar mem?ria RAM inutilizada
dev.off()     # Apagar todos os gr?ficos
cat("\014")   # Limpar o console

### Carregar pacote====
require(jagsUI)
require(energy)
KL <- function(x,y) {
  .5 * {{{sd(x)/sd(y)}^2} + {{{mean(y)-mean(x)}^2}/var(y)} - 1 + {2 * log(sd(y)/sd(x))}}
}
FI <- function(theta) {
  1/{plogis(theta) * {1-plogis(theta)}}
}

### Dados Simulados====
set.seed(1234) #  Estabelecer um seed para sempre gerar os mesmos resultados
V <- 20; N <- 200 # V: tamanho do teste; N: tamanho amostral
f <- rnorm(N)     # f: escore verdadeiro
diffs <- seq(-2,2,len=V) # diffs: dificuldades do itens
lambdas <- runif(V, .5, .8) # lambdas: cargas fatoriais
X <- sapply(1:V, function(g) diffs[g] + {lambdas[g]*f} + rnorm(N)) # Modelo fatorial linear
D <- {X > 0} * 1 # Dicotomização dos escores
colnames(D) <- paste("D", 1:V, sep="_") # Gerar nomes para as variáveis
dataList <- list( D=D, N=N, V=V ) # Agregar dados em uma lista

##############============ Modelo de Rasch====
RM <- " ### Rasch Model
model {
  # Priors
  for(j in 1:V) {
    delta[j] ~ dnorm(0,1)
  }
  for(i in 1:N) {
    theta[i] ~ dnorm(0,1)
  }

  # VerossimilhanÃÂ§a
  for (i in 1:N) {
    for(j in 1:V) {
      D[i,j] ~ dbern(ilogit(theta[i] - delta[j]))
      D_hat[i,j] <- ifelse(ilogit(theta[i] - delta[j]) > .5, 1, 0)
    }
  }
}"

### Run the chains====
# Name the parameters to be monitored
params <- c("delta","theta","D_hat")
# Random initial values
inits <- NULL
# Define some MCMC parameters for JAGS
nthin    = 1    # How Much Thinning?
nchains  = 5    # How Many Chains?
nburnin  = 1500 # How Many Burn-in Samples?
nsamples = 6500 # How Many Recorded Samples?
nadapt   = 2500 # How Many adaptation Samples?
# Calling jagsUI
set.seed(666)
model = textConnection(RM)
fit1 <- jagsUI::jags(dataList, inits=inits, params, model.file=model,
                     n.chains=nchains, n.adapt=nadapt, n.iter=nsamples,
                     n.burnin=nburnin, n.thin=nthin, DIC=T)

### Precision====
# Ajustar modelo de Rasch pelo método de máxima verossimilhança
rasch <- mirt::mirt(D, 1, "Rasch")
# Escores estimados pelo método Bayesiano
theta <- fit1$sims.list$theta
# Escores estimados pelo método de ML
score <- c(mirt::fscores(rasch, method="ML"))
# Informação de Fisher Bayesiano
FIB     <- apply(theta, 1, FI)
# Informação de Fisher Bayesiano
FisherB <- rowMeans(FIB)
# Informação de Fisher de máxima verossimilhança
FisherF <- FI(score)
# KLD entre escores estimados e distribuição a priori
KLD   <- sapply(1:ncol(theta), function(g) KL(theta[,g], c(scale(theta[,g]))))
# Precisão das estimativas Bayesianas dos escores
Prec  <- apply(theta, 2, sd)
# Acurácia das estimativas Bayesianas
Acc <- rowMeans(abs(colMeans(fit1$sims.list$D_hat) - D))
# Agregar resultados
Results <- cbind(Prec, KLD, FisherB, FisherF,
                 Acc, mirt::personfit(rasch)[,-c(2,4)],
                 EAP=colMeans(theta), ML=score, true=f)
# Gerar a Figura 1 do Editorial
#jpeg("Figure01.jpeg", height=20, width=40, units="cm", res=1200, pointsize=18)
par(mfrow=c(ncol(Results), ncol(Results)), mar=c(1,1,1,1))
for(i in 1:ncol(Results)) {
  for(j in 1:ncol(Results)) {
    if(i == j) {
      hist(Results[,i], main=colnames(Results)[j])
    } else if (j > i) {
      plot(Results[,i] ~ Results[,j], xaxt = 'n', yaxt = 'n',
           bty = 'n', pch = '', ylab = '', xlab = '')
      text(sum(range(Results[,j]))/2, sum(range(Results[,i]))/2, cex=1.5,
           format(round(energy::dcor(Results[,i], Results[,j]),3), nsmall = 3),
           col=grey(1-energy::dcor(Results[,i], Results[,j])))
    } else {
      if({{j == 1} & {i == 09}} | {{j == 2} & {i == 09}} |
         {{j == 1} & {i == 10}} | {{j == 2} & {i == 10}} |
         {{j == 1} & {i == 11}} | {{j == 2} & {i == 11}} |
         {{j == 3} & {i == 09}} | {{j == 4} & {i == 09}} |
         {{j == 3} & {i == 10}} | {{j == 4} & {i == 10}} |
         {{j == 3} & {i == 11}} | {{j == 4} & {i == 11}} ) {
        plot(Results[,i] ~ Results[,j], col="gray")
        temp <- smooth.spline(Results[,j] ~ Results[,i])
        lines(temp$y, temp$x, lwd=2)
      } else {
        plot(Results[,i] ~ Results[,j], col="gray")
        temp <- smooth.spline(Results[,i] ~ Results[,j])
        lines(temp$x, temp$y, lwd=2)
      }
    }
  }
}
#dev.off()
psych::pairs.panels(Results)

### Bias====
colMeans(cbind(abs(scale(colMeans(theta)) - scale(f)),
               abs(scale(score) - scale(f))))

####====---- THE END ----====####
