library(shiny)
library(xts)
library(knitr)
source("Reliability.R")

shinyServer(function(input, output) {
   
    impX <- function(){as.numeric(input$x)}
    impP <- function(){as.numeric(input$p)}
    impR <- function(){as.numeric(input$r)}
    impI <- function(){as.numeric(input$pil)}
    
    impmusaP <- function(){as.numeric(input$musaP)}
    impduaneP <- function(){as.numeric(input$duaneP)}
    implvP <- function(){as.numeric(input$lvP)}
    impokumotoP <- function(){as.numeric(input$okumotoP)}
    impgeoP <- function(){as.numeric(input$geoP)}
    impjmP <- function(){as.numeric(input$jmP)}
    
    impT <- function(){   
      inFile <- input$file
      if (is.null(inFile))
        return(NULL)
      data <- read.csv(inFile$datapath, header=FALSE, sep=",")
      data <- as.numeric(data)
      return(data)
    }
  
    
    output$total_plot <- renderPlot({
      
      inFile <- input$file
      if (is.null(inFile))
        return(NULL)
      t <- impT()
      musa.basic.par3 <- impX()
      jelinski.moranda.par3 <- impP()
      jelinski.moranda.par4 <- impR()
      
      duane.par1 <- duane(t)$rho
      duane.par2 <- duane(t)$theta
      
      lit.par1 <- littlewood.verall(t, linear = T)$theta0
      lit.par2 <- littlewood.verall(t, linear = T)$theta1
      lit.par3 <- littlewood.verall(t, linear = T)$rho
      
      mor.par1 <- moranda.geometric(t)$D
      mor.par2 <- moranda.geometric(t)$theta
      
      musa.par1 <- musa.okumoto(t)$theta0
      musa.par2 <- musa.okumoto(t)$theta1
      
      musa.basic.par1 <- musa.basic(t,musa.basic.par3)$beta0
      musa.basic.par2 <- musa.basic(t,musa.basic.par3)$beta1
      
      jelinski.moranda.par1 <- jelisnki.moranda(t, jelinski.moranda.par3, jelinski.moranda.par4)$N0
      jelinski.moranda.par2 <- jelisnki.moranda(t, jelinski.moranda.par3, jelinski.moranda.par4)$phi
      
      total.plot(jelinski.moranda.par1, jelinski.moranda.par2, jelinski.moranda.par3, jelinski.moranda.par4,
                 musa.basic.par1, musa.basic.par2,duane.par1, duane.par2, lit.par1, lit.par2, lit.par3, mor.par1,
                 mor.par2, musa.par1, musa.par2, t, linear = T, xlab = "Time ", main = "All models")
    })
    
    output$rel_plot <- renderPlot({
       
       inFile <- input$file
       if (is.null(inFile))
          return(NULL)
       t <- impT()
       musa.basic.par3 <- impX()
       jelinski.moranda.par3 <- impP()
       jelinski.moranda.par4 <- impR()
       
       duane.par1 <- duane(t)$rho
       duane.par2 <- duane(t)$theta
       
       lit.par1 <- littlewood.verall(t, linear = T)$theta0
       lit.par2 <- littlewood.verall(t, linear = T)$theta1
       lit.par3 <- littlewood.verall(t, linear = T)$rho
       
       mor.par1 <- moranda.geometric(t)$D
       mor.par2 <- moranda.geometric(t)$theta
       
       musa.par1 <- musa.okumoto(t)$theta0
       musa.par2 <- musa.okumoto(t)$theta1
       
       musa.basic.par1 <- musa.basic(t,musa.basic.par3)$beta0
       musa.basic.par2 <- musa.basic(t,musa.basic.par3)$beta1
       
       jelinski.moranda.par1 <- jelisnki.moranda(t, jelinski.moranda.par3, jelinski.moranda.par4)$N0
       jelinski.moranda.par2 <- jelisnki.moranda(t, jelinski.moranda.par3, jelinski.moranda.par4)$phi
       
       rel.plot(jelinski.moranda.par1, jelinski.moranda.par2, jelinski.moranda.par3, 
                jelinski.moranda.par4,
                musa.basic.par1, musa.basic.par2,
                duane.par1, duane.par2,
                lit.par1, lit.par2, lit.par3,
                mor.par1, mor.par2,
                musa.par1, musa.par2,
                t, linear = T, ymin = -1,
                ymax = 2.5, xlab = "Time", main = "Relative error")
    })
    
    output$rank <- renderDataTable({
       inFile <- input$file
       if (is.null(inFile))
          return(NULL)
       t <- impT()
       musa.basic.par3 <- impX()
       jelinski.moranda.par3 <- impP()
       jelinski.moranda.par4 <- impR()
       
       duane.par1 <- duane(t)$rho
       duane.par2 <- duane(t)$theta
       
       lit.par1 <- littlewood.verall(t, linear = T)$theta0
       lit.par2 <- littlewood.verall(t, linear = T)$theta1
       lit.par3 <- littlewood.verall(t, linear = T)$rho
       
       mor.par1 <- moranda.geometric(t)$D
       mor.par2 <- moranda.geometric(t)$theta
       
       musa.par1 <- musa.okumoto(t)$theta0
       musa.par2 <- musa.okumoto(t)$theta1
       
       musa.basic.par1 <- musa.basic(t,musa.basic.par3)$beta0
       musa.basic.par2 <- musa.basic(t,musa.basic.par3)$beta1
       
       jelinski.moranda.par1 <- jelisnki.moranda(t, jelinski.moranda.par3, jelinski.moranda.par4)$N0
       jelinski.moranda.par2 <- jelisnki.moranda(t, jelinski.moranda.par3, jelinski.moranda.par4)$phi
       
       return(rank.rel(jelinski.moranda.par1, jelinski.moranda.par2, jelinski.moranda.par3, 
                          jelinski.moranda.par4,
                          musa.basic.par1, musa.basic.par2, duane.par1, duane.par2, lit.par1, 
                          lit.par2, lit.par3, mor.par1,
                          mor.par2, musa.par1, musa.par2, t, linear = T))
    }, options = list(pageLength = 10))
    
    output$pred_plot <- renderPlot({
       
       inFile <- input$file
       if (is.null(inFile))
          return(NULL)
       t <- impT()
       musa.basic.par3 <- impX()
       jelinski.moranda.par3 <- impP()
       jelinski.moranda.par4 <- impR()
       
       duane.par1 <- duane(t)$rho
       duane.par2 <- duane(t)$theta
       
       lit.par1 <- littlewood.verall(t, linear = T)$theta0
       lit.par2 <- littlewood.verall(t, linear = T)$theta1
       lit.par3 <- littlewood.verall(t, linear = T)$rho
       
       mor.par1 <- moranda.geometric(t)$D
       mor.par2 <- moranda.geometric(t)$theta
       
       musa.par1 <- musa.okumoto(t)$theta0
       musa.par2 <- musa.okumoto(t)$theta1
       
       musa.basic.par1 <- musa.basic(t,musa.basic.par3)$beta0
       musa.basic.par2 <- musa.basic(t,musa.basic.par3)$beta1
       
       jelinski.moranda.par1 <- jelisnki.moranda(t, jelinski.moranda.par3, jelinski.moranda.par4)$N0
       jelinski.moranda.par2 <- jelisnki.moranda(t, jelinski.moranda.par3, jelinski.moranda.par4)$phi

       t_aux <- cumsum(t)
       tmax <- t_aux[length(t)]
       interval <- impI()
       new_t <- tmax + 0:interval
       
       musaP <- impmusaP()
       duaneP <- impduaneP()
       lvP <- implvP()
       okumotoP <- impokumotoP()
       geoP <- impgeoP()
       jmP <- impjmP()
       
       prediction.plot(musaP, duaneP, lvP, okumotoP, geoP, jmP,
                       jelinski.moranda.par1, jelinski.moranda.par2, jelinski.moranda.par3, jelinski.moranda.par4,
                       musa.basic.par1, musa.basic.par2,duane.par1, duane.par2, lit.par1, lit.par2, lit.par3, mor.par1,
                       mor.par2, musa.par1, musa.par2, rep(1,interval), linear = T, xlab = "Time ", main = "All models")
    })
    
    
    output$Time_series <- renderPlot({
    
       require(tseries)
       require(forecast)
       
       inFile <- input$file
       if (is.null(inFile))
          return(NULL)
       dataT <- read.csv(inFile$datapath, header=FALSE, sep=",")
       dataT <- t(dataT)
       
       # plot(dataT)
       
       dataTlog <- log(dataT + 1)
       modelo <- Arima(dataTlog,c(7,0,7))
       
       #-------------------------------------------
       plot(dataT,type="l",xlab="Time",ylab="",xlim = c(0,length(dataT)+input$prediccion+1),ylim=c(0,16000))
       lines(exp(fitted(modelo)),col="red")
       
       prediccion <- predict(modelo,input$prediccion)
       lines(exp(prediccion$pred),type="o")
       lines(exp(prediccion$pred-prediccion$se),col="blue",type="o")
       lines(exp(prediccion$pred+prediccion$se),col="blue",type="o")
       
       
   })
    
    output$Time_series2 <- renderDygraph({
       
       require(tseries)
       require(forecast)
       
       inFile <- input$file
       if (is.null(inFile))
          return(NULL)
       dataT <- read.csv(inFile$datapath, header=FALSE, sep=",")
       dataT <- t(dataT)
       
       # plot(dataT)
       
       dataTlog <- log(dataT + 1)
       modelo <- Arima(dataTlog,c(7,0,7))
       #-------------------------------------------
       prediccion <- predict(modelo,input$prediccion)
       
       
       dPrediccion <- data.frame(LowLimit = exp(prediccion$pred-prediccion$se),
                                 Prediction = exp(prediccion$pred),
                                 UpperLimit = exp(prediccion$pred+prediccion$se))
       
       nas <- matrix(rep(NA,length(dataT)*3),ncol=3)
       colnames(nas) <- names(dPrediccion)
       dataT2 <- rbind(nas,dPrediccion)
       
       nas2 <- matrix(rep(NA,length(prediccion$pred)),ncol=1)
       colnames(nas2) <- names(dataT)
       dataT3 <- rbind(dataT,nas2)
       
       dataFinal <- cbind(Historic = dataT3,dataT2)
       Tiempo= as.Date(1:length(dataT3))
       
       dataFinal<-xts(dataFinal, order.by = Tiempo)
       
       dygraph(dataFinal, "Tecnology Failtures") %>%
          dySeries("Historic", label = "Actual") %>%
          dySeries(c("LowLimit", "Prediction", "UpperLimit"), label = "Predicted")
       
       
    })

    
    # the values
    sliderValues <- reactive({
       
       
       require(tseries)
       require(forecast)
       
       inFile <- input$file
       if (is.null(inFile))
          return(NULL)
       dataT <- read.csv(inFile$datapath, header=FALSE, sep=",")
       dataT <- t(dataT)
       
       # plot(dataT)
       
       dataTlog <- log(dataT + 1)
       modelo <- Arima(dataTlog,c(7,0,7))
       #-------------------------------------------
       prediccion <- predict(modelo,input$prediccion)
       
       
       res <- data.frame(LowLimit = exp(prediccion$pred-prediccion$se),
                  Prediction = exp(prediccion$pred),
                  UpperLimit = exp(prediccion$pred+prediccion$se))
       #kable(res,format = "html" )
    }) 
    
    # Show the values using an HTML table
    output$Prediction_values <- renderTable({
       sliderValues()
    })
    
  }
)
