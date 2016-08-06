
require(tseries)
require(forecast)

dataT <- read.csv("CS (2).csv", header=FALSE, sep=",")
dataT <- t(dataT)
dataTlog <- log(dataT + 1)
modelo <- Arima(dataTlog,c(7,0,7))

#-------------------------------------------
plot(dataT,type="l",xlab="Time",ylab="",xlim = c(0,150),ylim=c(0,16000))
lines(exp(fitted(modelo)),col="red")

prediccion <- predict(modelo,10)
lines(exp(prediccion$pred),type="o")
lines(exp(prediccion$pred-prediccion$se),col="blue",type="o")
lines(exp(prediccion$pred+prediccion$se),col="blue",type="o")
