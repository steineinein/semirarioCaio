#------------------------------------------------------------#
# Para rodar este programa  deixe no objeto fit.model a saída 
# do  ajuste  da  regressão com erros n.inversa. Deixe  
# os dados disponíveis  através do comando attach(...). 
# Depois  use  o comando source(...) no R ou S-Plus para
# executar o programa. A sequência de comandos é a seguinte:
#
#        > fit.model <- ajuste
#        > attach(dados)
#        > source("diag_ninv")
#
# A saída terá quatro gráficos: de pontos de alacanca, de
# pontos influentes  e  dois de resíduos. Para identificar 
# os pontos que  mais  se destacam usar o comando 
# identify(...). Se por exemplo se destacam três pontos 
# no plot(fitted(fit.model),h,...),  após esse comando coloque
#     
#        > identify(fitted(fit.model),h,n=3)
#
# O mesmo pode ser feito nos demais gráficos. Nos gráficos de 
# resíduos foram colocados os limites ylim=c(a-1,b+1), em que a
# é o menor valor e b o maior valor para o resíduo. 
#------------------------------------------------------------#
X <- model.matrix(fit.model)
n <- nrow(X)
p <- ncol(X)
w <- fit.model$weights
W <- diag(w)
H <- solve(t(X)%*%W%*%X)
H <- sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
h <- diag(H)
soma <- resid(fit.model, type="pearson")
soma <- sum(soma^2)
fi <- (n-p)/soma
ts <- resid(fit.model,type="pearson")*sqrt(fi/(1-h))
td <- resid(fit.model,type="deviance")*sqrt(fi/(1-h))
par(mfrow=c(2,2))
di <- (h/(1-h))*(ts^2)
a <- max(td)
b <- min(td)
#
plot(fitted(fit.model),h,xlab="Valor Ajustado", ylab="Medida h", pch=16)
#identify(fitted(fit.model), h, n=3)
#
plot(fitted(fit.model), di,xlab="Valor Ajustado", ylab="Distância de Cook", pch=16)
#identify(fitted(fit.model),di, n=6)
#
plot(fitted(fit.model),td,xlab="Valor Ajustado", ylab="Resíduo Componente do Desvio",
     ylim=c(b-1,a+1),pch=16)
abline(2,0,lty=2)
abline(-2,0,lty=2)
#identify(fitted(fit.model),td, n=1)
eta <- predict(fit.model)
z <- eta + resid(fit.model, type="pearson")/sqrt(w)
plot(predict(fit.model),z,xlab="Preditor Linear", 
     ylab="Variavel z", pch=16)
lines(smooth.spline(predict(fit.model), z, df=2))
par(mfrow=c(1,1))
#------------------------------------------------------------#