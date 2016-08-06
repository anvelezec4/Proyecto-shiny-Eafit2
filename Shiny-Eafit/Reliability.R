# This file includes following functions
# duane, littlewood.verall, moranda.geometric, musa.okumoto, musa.basic, jelinski.moranda


# Function for estimating the parameters of the Duane model 
"duane" <- function(t, init = c(1, 1), method = "Nelder-Mead", maxit = 10000, ...) 
{
  eq1 <- function(rho, theta, t) 
  {
    n <- length(t)
    t <- cumsum(t)
    i <- seq(along = t)
    tn <- t[length(t)]
    rho - n / (tn^theta)
  }
  
  eq2 <- function(rho, theta, t) 
  {
    i <- seq(along = t[1:length(t) - 1])
    n <- length(t)
    t <- cumsum(t)
    tn <- t[length(t)]
    theta - n / (sum(log(tn / t[i])))  
  }
  
  # Merging the two parts
  global <- function(vec, t) 
  {
    rho <- vec[1]
    theta <- vec[2]
    v1 <- eq1(rho, theta, t)^2
    v2 <- eq2(rho, theta, t)^2
    v1 + v2
  }
  
  # Minimum search for rho and theta
  res <- optim(init, global, t = t, method = method, 
               control = list(maxit = maxit, ...))
  rho <- res$par[1]
  theta <- res$par[2]
  
  return(list(rho = rho, theta = theta))
}


# Function for plotting the mean value function of the Duane model

"duane.plot" <- function(rho, theta, t, xlab = "time", 
                         ylab = "Cumulated failures and estimated mean value function", 
                         main = NULL) 
{
  plot(cumsum(t), 1:length(t), type = "o", pch = 3, cex = 0.8, xlab = xlab,      
       ylab = ylab)                                         
  lines(cumsum(t), mvf.duane(rho, theta, cumsum(t)), lty = 2, type = "o", 
        pch = 4, cex = 0.8, col = "mediumblue")   
  legend("bottomright", c("Data", "Duane"), col = c("black", "mediumblue"), pch = c(3, 4), 
         lty = c(1, 2))    
  
  if(!is.null(main))
  {
    title(main)
  }                                                                                 
}


# Function for estimating the parameters of the Littlewood-Verall model 

"littlewood.verall" <- function(t, linear = T, init = c(1, 1, 1), method = "Nelder-Mead", 
                                maxit = 10000, ...) 
{
  if(linear == T) 
  {
    eq1 <- function(theta0, theta1, rho, t) 
    {
      n <- length(t)
      i <- seq(along = t)
      s1 <- sum(log(theta0 + theta1 * i))
      s2 <- sum(log(theta0 + theta1 * i + t[i]))
      (n / rho) + s1 - s2
    }
    
    eq2 <- function(theta0, theta1, rho, t) 
    {
      i <- seq(along = t)
      s1 <- rho * sum(1 / (theta0 + theta1 * i))
      s2 <- (rho + 1) * sum(1 / (theta0 + theta1 * i + t[i]))
      s1 - s2
    }
    
    eq3 <- function(theta0, theta1, rho, t) 
    {
      i <- seq(along = t)
      s1 <- rho * sum(i / (theta0 + theta1 * i))
      s2 <- (rho + 1) * sum(i / (theta0 + theta1 * i + t[i]))
      s1 - s2
    }
    
    # Merging the three parts
    global.linear <- function(vec, t) 
    {
      theta0 <- vec[1]
      theta1 <- vec[2]
      rho <- vec[3]
      v1 <- eq1(theta0, theta1, rho, t)^2
      v2 <- eq2(theta0, theta1, rho, t)^2
      v3 <- eq3(theta0, theta1, rho, t)^2
      v1 + v2 + v3
    }
    
    # Minimum search for theta0, theta1 and rho
    res <- optim(init, global.linear, t = t, method = method, 
                 control = list(maxit = maxit, ...))
    theta0 <- res$par[1]
    theta1 <- res$par[2]
    rho <- res$par[3]
    
    # Repeated minimum search for theta0, theta1 and rho
    res <- optim(c(theta0, theta1, rho), global.linear, t = t, method = method,
                 control = list(maxit = maxit, ...))
  }
  else 
  {
    eq1 <- function(theta0, theta1, rho, t) 
    {
      n <- length(t)
      i <- seq(along = t)
      print(theta0 + theta1 * i^2)
      s1 <- sum(log(theta0 + theta1 * i^2))
      s2 <- sum(log(theta0 + theta1 * i^2 + t[i]))
      (n / rho) + s1 - s2
    }
    
    eq2 <- function(theta0, theta1, rho, t) 
    {
      i <- seq(along = t)
      s1 <- rho * sum(1 / (theta0 + theta1 * i^2))
      s2 <- (rho + 1) * sum(1 / (theta0 + theta1 * i^2 + t[i]))
      s1 - s2
    }
    
    eq3 <- function(theta0, theta1, rho, t) 
    {
      i <- seq(along = t)
      s1 <- rho * sum((i^2) / (theta0 + theta1 * i^2))
      s2 <- (rho + 1) * sum((i^2) / (theta0 + theta1 * i^2 + t[i]))
      s1 - s2
    }
    
    # Merging the three parts
    global.quadratic <- function(vec, t) 
    {
      theta0 <- vec[1]
      theta1 <- vec[2]
      rho <- vec[3]
      v1 <- eq1(theta0, theta1, rho, t)^2
      v2 <- eq2(theta0, theta1, rho, t)^2
      v3 <- eq3(theta0, theta1, rho, t)^2
      v1 + v2 + v3
    }
    
    # Minimum search for theta0, theta1 and rho
    res <- optim(init, global.quadratic, t = t, method = method, 
                 control = list(maxit = maxit, ...))
    theta0 <- res$par[1]
    theta1 <- res$par[2]
    rho <- res$par[3]
    
    # Repeated minimum search for theta0, theta1 and rho
    res <- optim(c(theta0, theta1, rho), global.quadratic, t = t, method = method,
                 control = list(maxit = maxit, ...))
  }
  
  theta0 <- res$par[1]
  theta1 <- res$par[2]
  rho <- res$par[3]
  
  return(list(theta0 = theta0, theta1 = theta1, rho = rho))
}


# Function for plotting the mean value function of the Littlewood-Verall model

"littlewood.verall.plot" <- function(theta0, theta1, rho, t, linear = T, xlab = "time", 
                                     ylab = "Cumulated failures and estimated mean value function", main = NULL) 
{
  if(linear == T) 
  {
    mvf.ver <- mvf.ver.lin(theta0, theta1, rho, cumsum(t))
    
    plot(cumsum(t), 1:length(t), type = "o", pch = 3, cex = 0.8, xlab=xlab,
         ylab = ylab)
  }
  else 
  {
    mvf.ver <- mvf.ver.quad(theta0, theta1, rho, cumsum(t))
    
    plot(cumsum(t), 1:length(t), type = "o", pch = 3, cex = 0.8, xlab = xlab,
         ylab = ylab)
  }
  
  lines(cumsum(t), mvf.ver, lty=2, type = "o", pch = 4, cex = 0.8, col="mediumblue")
  legend("bottomright", c("Data", "Littlewood-Verall"), col=c("black", "mediumblue"), 
         pch=c(3, 4), lty=c(1, 2))  
  
  if(!is.null(main))
  {
    title(main)
  }  
}


# Function for estimating the parameters of the Moranda-Geometric model 

"moranda.geometric" <- function(t, init = c(0, 1), tol = .Machine$double.eps^0.25) 
{
  # Function for computation of Dhat
  Dhat <- function(phi, t) 
  {
    n <- length(t)
    i <- seq(along = t)
    return(phi * n / (sum(phi^i * t[i])))
  }
  
  # Squared Maximum Likelihood estimate for phihat
  phihat <- function(phi, t) 
  {
    i <- seq(along = t)
    n <- length(t)
    return((sum(i * (phi^i) * t[i]) / sum((phi^i) * t[i]) - (n + 1) / 2)^2)
  }
  
  # Minimum search for phihat
  min <- optimize(phihat, init, tol = tol, t = t)  
  phi <- min$minimum
  D <- Dhat(phi, t)
  theta<- -log(phi)
  
  return(list(D = D, theta = theta))
}


# Function for plotting the mean value function of the Moranda-Geometric model

"moranda.geometric.plot" <- function(D, theta, t, xlab = "time", 
                                     ylab = "Cumulated failures and estimated mean value function", main = NULL) 
{
  plot(cumsum(t), 1:length(t), type = "o", pch = 3, cex = 0.8, xlab = xlab,      
       ylab = ylab)                                         
  lines(cumsum(t), mvf.mor(D, theta, cumsum(t)), lty = 2, type = "o", pch = 4, 
        cex = 0.8, col = "mediumblue")  
  legend("bottomright", c("Data", "Moranda-Geometric"), col=c("black", "mediumblue"), 
         pch=c(3, 4), lty=c(1, 2)) 
  
  if(!is.null(main))
  {
    title(main)
  }  
}


# Function for estimating the parameters of the Musa-Okumoto model 

"musa.okumoto" <- function(t, init = c(0, 1), tol = .Machine$double.eps^0.25) 
{
  # Function for computation of theta0hat
  theta0hat <- function(theta1, t) 
  {
    n <- length(t)
    t <- cumsum(t)
    tn <- t[length(t)]
    return(n / log(1 + theta1 * tn))
  }
  
  # Squared Maximum Likelihood estimate for theta1hat
  theta1hat <- function(theta1, t) 
  {
    i <- seq(along = t)
    n <- length(t)
    t <- cumsum(t)
    tn <- t[length(t)]
    return((1 / theta1 * sum(1 / (1 + theta1 * t[i])) - n * tn / 
              ((1 + theta1 * tn) * log(1 + theta1 * tn)))^2)
  }
  
  # Minimum search for theta1hat
  min <- optimize(theta1hat, init, tol = tol, t = t)
  theta1 <- min$minimum
  theta0 <- theta0hat(theta1, t)
  
  return(list(theta0 = theta0, theta1 = theta1))
}


# Function for plotting the mean value function of the Musa-Okumoto model

"musa.okumoto.plot" <- function(theta0, theta1, t, xlab = "time", 
                                ylab = "Cumulated failures and estimated mean value function", main = NULL) 
{
  plot(cumsum(t), 1:length(t), type = "o", pch = 3, cex = 0.8, xlab = xlab,
       ylab = ylab)
  lines(cumsum(t), mvf.musa(theta0, theta1, cumsum(t)), lty = 2, type = "o", 
        pch = 4, cex = 0.8, col = "mediumblue")
  legend("bottomright", c("Data", "Musa-Okumoto"), col = c("black", "mediumblue"), 
         pch = c(3, 4), lty = c(1, 2)) 
  
  if(!is.null(main))
  {
    title(main)
  }  
}


# Function for plotting all mean value functions of all models in one plot.

"total.plot" <- function(jelinski.moranda.par1, jelinski.moranda.par2, jelinski.moranda.par3, jelinski.moranda.par4, 
                         musa.basic.par1, musa.basic.par2, duane.par1, duane.par2, lit.par1, lit.par2, lit.par3, 
                         mor.par1, mor.par2, musa.par1, musa.par2, t, linear = T, xlab = "Time",
                         ylab = "Cumulated failures and estimated mean value functions", main = NULL)
{
  if(linear == T) 
  {
    mvf.ver <- mvf.ver.lin(lit.par1, lit.par2, lit.par3, cumsum(t))
  }
  else 
  {
    mvf.ver <- mvf.ver.quad(lit.par1, lit.par2, lit.par3, cumsum(t))
  }
  
  plot(cumsum(t), 1:length(t), type = "o", pch = 3, cex = 0.8, xlab = xlab,      
       ylab = ylab)
  lines(cumsum(t), mvf.duane(duane.par1, duane.par2, cumsum(t)), lty = 2,          
        type = "o", pch = 4, cex = 0.8, col = "mediumblue")     
  lines(cumsum(t), mvf.ver, lty = 3, type = "o", pch = 8, cex = 0.8, 
        col = "steelblue")    
  lines(cumsum(t), mvf.mor(mor.par1, mor.par2, cumsum(t)), lty = 4,          
        type = "o", pch = 0, cex = 0.8, col = "skyblue")   
  lines(cumsum(t), mvf.musa(musa.par1, musa.par2, cumsum(t)), lty = 5,          
        type = "o", pch = 15, cex = 0.8, col = "royalblue")   
  lines(cumsum(t), mvf.musa.basic(musa.basic.par1,musa.basic.par2, cumsum(t)), lty = 6,
        type = "o", pch = 16, cex = 0.8, col = "lightseagreen")   
  lines(cumsum(t), mvf.jelinski.moranda(jelinski.moranda.par1, jelinski.moranda.par2, jelinski.moranda.par3,
                                        jelinski.moranda.par4, cumsum(t)), lty = 7,
        type = "o", pch = 17, cex = 0.8, col = "cornflowerblue")   
  legend("bottomright", c("Data", "D", "LV", "MG", "MO", "BM", "JM"), 
         col = c("black", "mediumblue", "steelblue", "skyblue", "royalblue", "lightseagreen", "cornflowerblue" ),
         pch = c(3, 4, 8, 0, 15, 16, 17),lty = c(1, 2, 3, 4, 5, 6, 7), horiz=T)
  
  if(!is.null(main))
  {
    title(main)
  }  
}



# Function for calculating ranking based on mean absolute relative errors.
"rank.rel" <- function(jelinski.moranda.par1, jelinski.moranda.par2, jelinski.moranda.par3, jelinski.moranda.par4,
                       musa.basic.par1, musa.basic.par2, duane.par1, duane.par2, lit.par1, lit.par2, lit.par3, mor.par1,
                       mor.par2, musa.par1, musa.par2, t, linear = T)
{
  l <- 1:length(t)
  Rrel.duane <- mean(abs((mvf.duane(duane.par1, duane.par2, cumsum(t)) - l) / l))
  
  if(linear == T) 
  {
    Rrel.ver <- mean(abs((mvf.ver.lin(lit.par1, lit.par2, lit.par3, cumsum(t)) - l) / l))
  }
  else 
  {
    Rrel.ver <- mean(abs((mvf.ver.quad(lit.par1, lit.par2, lit.par3, cumsum(t)) - l) / l))
  } 
  
  Rrel.mor <- mean(abs((mvf.mor(mor.par1, mor.par2, cumsum(t)) - l) / l))
  Rrel.musa <- mean(abs((mvf.musa(musa.par1, musa.par2, cumsum(t)) - l) / l))
  Rrel.musa.basic <- mean(abs((mvf.musa.basic(musa.basic.par1, musa.basic.par2, cumsum(t)) - l) / l))
  Rrel.jelinski.moranda <- mean(abs((mvf.jelinski.moranda(jelinski.moranda.par1, jelinski.moranda.par2,
                                                          jelinski.moranda.par3, jelinski.moranda.par4, cumsum(t)) - l) / l))
  
  mae <-c(Rrel.duane, Rrel.ver, Rrel.mor, Rrel.musa, Rrel.musa.basic, Rrel.jelinski.moranda) 
  meths <- c("Duane", "Littlewood-Verall", "Moranda Geometric", "Musa-Okumoto", "Musa Basic","Jelinski-Moranda")
  m2 <- cbind(meths, mae)
  colnames(m2, do.NULL = FALSE)
  colnames(m2) <- c("Method","Mean Absolute Relative Error")
  return(m2)	
}

# Function for plotting all relative errors in one plot.
"rel.plot" <- function(jelinski.moranda.par1, jelinski.moranda.par2, jelinski.moranda.par3, jelinski.moranda.par4,
                       musa.basic.par1, musa.basic.par2, duane.par1, duane.par2, lit.par1, lit.par2, lit.par3, mor.par1,
                       mor.par2, musa.par1, musa.par2, t, linear = T, ymin, ymax, xlab = "Time", ylab = "Relative error", main = NULL)
{
  l <- 1:length(t)
  rel.duane <- (mvf.duane(duane.par1, duane.par2, cumsum(t)) - l) / l
  
  if(linear == T) 
  {
    rel.ver <- (mvf.ver.lin(lit.par1, lit.par2, lit.par3, cumsum(t)) - l) / l
  }
  else 
  {
    rel.ver <- (mvf.ver.quad(lit.par1, lit.par2, lit.par3, cumsum(t)) - l) / l
  } 
  
  rel.mor <- (mvf.mor(mor.par1, mor.par2, cumsum(t)) - l) / l
  rel.musa <- (mvf.musa(musa.par1, musa.par2, cumsum(t)) - l) / l
  rel.musa.basic <- (mvf.musa.basic(musa.basic.par1, musa.basic.par2, cumsum(t)) - l) / l
  rel.jelinski.moranda <- (mvf.jelinski.moranda(jelinski.moranda.par1, jelinski.moranda.par2,
                                                jelinski.moranda.par3, jelinski.moranda.par4, cumsum(t)) - l) / l
  
  if(missing(ymin))
  {
    ymin <- min(c(rel.duane, rel.ver, rel.mor, rel.musa, rel.musa.basic))
  }
  
  if(missing(ymax))
  {
    ymax <- max(c(rel.duane, rel.ver, rel.mor, rel.musa, rel.musa.basic))
  }
  
  plot(cumsum(t), rel.duane, type = "o", pch = 4, cex = 0.8, col = "mediumblue",xlab = xlab, ylim = c(ymin, ymax), lty = 2, ylab = ylab)
  lines(cumsum(t), rel.ver, lty = 3, type = "o", pch = 8, cex = 0.8, col = "steelblue")     
  lines(cumsum(t), rel.mor, lty = 4, type = "o", pch = 0, cex = 0.8, col = "skyblue")    
  lines(cumsum(t), rel.musa, lty = 5, type = "o", pch = 15, cex = 0.8, col = "royalblue")   
  lines(cumsum(t), rel.musa.basic, lty = 6, type = "o", pch = 16, cex = 0.8, col = "lightseagreen")
  lines(cumsum(t), rel.jelinski.moranda, lty = 7, type = "o", pch = 17, cex = 0.8, col = "cornflowerblue")
  segments(-max(cumsum(t)) * 1.05, 0, max(cumsum(t)) * 1.05, 0, lty = 3)
  legend("topright", c("D", "LV", "MG", "MO", "BM", "JM"), 
         col = c("mediumblue", "steelblue", "skyblue", "royalblue", "lightseagreen", "cornflowerblue"), pch = c(4, 8, 0, 15, 16, 17), 
         lty = c(2, 3, 4, 5, 6, 7), horiz=T)
  
  if(!is.null(main))
  {
    title(main)
  }  
}


# Function for plotting the predictions of all models in one plot.

"prediction.plot" <- function(musaP, duaneP, lvP, okumotoP, geoP, jmP,
                              jelinski.moranda.par1, jelinski.moranda.par2, jelinski.moranda.par3, jelinski.moranda.par4, 
                         musa.basic.par1, musa.basic.par2, duane.par1, duane.par2, lit.par1, lit.par2, lit.par3, 
                         mor.par1, mor.par2, musa.par1, musa.par2, t, linear = T, xlab = "Time",
                         ylab = "Cumulated failures and estimated mean value functions", main = NULL)
{
   if(linear == T) 
   {
      mvf.ver <- mvf.ver.lin(lit.par1, lit.par2, lit.par3, cumsum(t))
   }
   else 
   {
      mvf.ver <- mvf.ver.quad(lit.par1, lit.par2, lit.par3, cumsum(t))
   }
   
   if (lvP==F) mvf.ver <- 0
   
   if (duaneP==T) mvf.duane <- mvf.duane(duane.par1, duane.par2, cumsum(t))
   else mvf.duane <- 0
   
   if (geoP == T )    mvf.mor <- mvf.mor(mor.par1, mor.par2, cumsum(t))
   else mvf.mor <- 0
   
   if (okumotoP == T) mvf.musa <- mvf.musa(musa.par1, musa.par2, cumsum(t))
   else mvf.musa <- 0
   
   if (musaP == T) mvf.musa.basic <- mvf.musa.basic(musa.basic.par1, musa.basic.par2, cumsum(t))
   else  mvf.musa.basic <- 0
   
   if (jmP == T) mvf.jelinski.moranda <- mvf.jelinski.moranda(jelinski.moranda.par1, jelinski.moranda.par2, jelinski.moranda.par3, jelinski.moranda.par4, cumsum(t))
   else  mvf.jelinski.moranda <- 0
   
   
   ymin <- min(c(mvf.duane, mvf.ver, mvf.mor, mvf.musa, mvf.musa.basic))
   ymax <- max(c(mvf.duane, mvf.ver, mvf.mor, mvf.musa, mvf.musa.basic))
   
   plot(cumsum(t), 1:length(t), type = "o", pch = 3, cex = 0.8, col="white", xlab = xlab, ylab = ylab, ylim = c(ymin, ymax))
   
   if (duaneP == T) lines(cumsum(t), mvf.duane(duane.par1, duane.par2, cumsum(t)), lty = 2, type = "o", pch = 4, cex = 0.8, col = "mediumblue")     
   
   if (lvP == T) lines(cumsum(t), mvf.ver, lty = 3, type = "o", pch = 8, cex = 0.8, col = "steelblue")
   
   if (geoP == T) lines(cumsum(t), mvf.mor(mor.par1, mor.par2, cumsum(t)), lty = 4, type = "o", pch = 0, cex = 0.8, col = "skyblue")   
   
   if (okumotoP == T) lines(cumsum(t), mvf.musa(musa.par1, musa.par2, cumsum(t)), lty = 5, type = "o", pch = 15, cex = 0.8, col = "royalblue")
   
   if (musaP == T) lines(cumsum(t), mvf.musa.basic(musa.basic.par1,musa.basic.par2, cumsum(t)), lty = 6,type = "o", pch = 16, cex = 0.8, col = "lightseagreen")   
   
   if (jmP == T) lines(cumsum(t), mvf.jelinski.moranda(jelinski.moranda.par1, jelinski.moranda.par2, jelinski.moranda.par3, jelinski.moranda.par4, cumsum(t)), lty = 7, type = "o", pch = 17, cex = 0.8, col = "cornflowerblue")
   
   legend("top", c("R", "D", "LV", "MG", "MO", "BM", "JM"), 
          col = c("black", "mediumblue", "steelblue", "skyblue", "royalblue", "lightseagreen", "cornflowerblue" ),
          pch = c(3, 4, 8, 0, 15, 16, 17),lty = c(1, 2, 3, 4, 5, 6, 7),horiz = T)
   
   if(!is.null(main))
   {
      title(main)
   }  
}

# Mean value function for Duane model

"mvf.duane" <- function(rho, theta, t) 
{
  if(rho <= 0 || theta <= 0)
  {
    stop("rho and theta should be larger than 0")
  }
  
  if(length(rho) != 1 || length(theta) != 1)
  {
    stop("rho and theta should have length 1")
  }
  
  return(rho * t^theta)
}


# Mean value function Moranda-Geometric model

"mvf.mor" <- function(D, theta, t) 
{
  if(theta <= 0)
  {
    stop("theta should be larger than 0")
  }
  
  if(length(D) != 1 || length(theta) != 1)
  {
    stop("D and theta should have length 1")
  }
  
  return(1 / theta * log((D * theta * exp(theta)) * t + 1))
}


# Mean Value Funktion Musa-Okumoto model

"mvf.musa" <- function(theta0, theta1, t) 
{
  if(length(theta0) != 1 || length(theta1) != 1)
  {
    stop("theta0 and theta1 should have length 1")
  }
  
  return(theta0 * log(theta1 * t + 1))
}


# Mean value function Littlewood-Verall model linear

"mvf.ver.lin" <- function(theta0, theta1, rho, t) 
{
  if(theta1 == 0)
  {
    stop("theta1 should not be equal 0")
  }
  
  if(length(theta0) != 1 || length(theta1) != 1 || length(rho) != 1)
  {
    stop("theta0, theta1 and rho should have length 1")
  }
  
  return(1 / theta1 * sqrt(theta0^2 + 2 * theta1 * t * rho))
}


# Mean value function Littlewood-Verall model quadratic 

"mvf.ver.quad" <- function(theta0, theta1, rho, t) 
{
  if(theta1 == 0)
  {
    stop("theta1 should not be equal 0")
  }
  
  if(length(theta0) != 1 || length(theta1) != 1 || length(rho) != 1)
  {
    stop("theta0, theta1 and rho should have length 1")
  }
  
  v1 <- (rho - 1)^(1 / 3) / ((18 * theta1)^(1 / 3))
  v2 <- 4 * (theta0^3) / (9 * (rho - 1)^2 * theta1)
  Q1 <- (cumsum(t) + (cumsum(t)^2 + v2)^(1 / 2))^(1 / 3)
  Q2 <- (cumsum(t) - (cumsum(t)^2 + v2)^(1 / 2))^(1 / 3)
  
  return(3 * v1 * (Q1 + Q2))
}

# Function for estimating the parameters of the Basic Musa model 

"musa.basic" <- function(t, x, init = c(0, 1), tol = .Machine$double.eps^0.25) 
{
  # Function for computation of beta0hat
  beta0hat <- function(beta1, t, x) 
  {
    n <- length(t)
    t <- cumsum(t)
    tn <- t[length(t)]
    return(n /( 1- exp(-beta1 * (tn + x))))
  }
  
  # Squared Maximum Likelihood estimate for beta1hat
  beta1hat <- function(beta1, t) 
  {
    i <- seq(along = t)
    n <- length(t)
    t <- cumsum(t)
    tn <- t[length(t)]
    return((n / beta1 - n*(tn + x)/(exp(beta1 * (tn + x))-1) - sum(t[i]))^2)
  }
  
  # Minimum search for theta1hat
  min <- optimize(beta1hat, init, tol = tol, t = t)
  beta1 <- min$minimum
  beta0 <- beta0hat(beta1, t, x)
  
  return(list(beta0 = beta0, beta1 = beta1))
}


# Function for plotting the mean value function of the Musa-Okumoto model

"musa.basic.plot" <- function(beta0, beta1, t, xlab = "time", 
                              ylab = "Cumulated failures and estimated mean value function", main = NULL) 
{
  plot(cumsum(t), 1:length(t), type = "o", pch = 3, cex = 0.8, xlab = xlab,
       ylab = ylab)
  lines(cumsum(t), mvf.musa.basic(beta0, beta1, cumsum(t)), lty = 2, type = "o", 
        pch = 4, cex = 0.8, col = "mediumblue")
  legend("bottomright", c("Data", "Basic Musa"), col = c("black", "mediumblue"), 
         pch = c(3, 4), lty = c(1, 2)) 
  
  if(!is.null(main))
  {
    title(main)
  }  
}

# Mean Value Funktion Basic Musa model

"mvf.musa.basic" <- function(beta0, beta1, t) 
{
  if(length(beta0) != 1 || length(beta1) != 1)
  {
    stop("beta0 and beta1 should have length 1")
  }
  
  return(beta0 *(1 - exp(-beta1 * t)))
}

# Function for estimating the parameters of the Jelinski Moranda model 

"jelisnki.moranda" <- function(t, p, r, init = c(length(t), 150), tol = .Machine$double.eps^0.25) 
{
  # Function for computation of phihat
  phihat <- function(N0, t, p, r) 
  {
    n <- length(t)
    i <- seq(along = t)
    s1 <- sum(t[i])
    s2 <- sum((i-1)*(p-r)*t[i])
    return( n/( N0*s1 - s2 ))
  }
  
  # Maximum Likelihood estimate for N0hat
  N0hat <- function(t, p, r) 
  {
    n <- length(t)
    s <- n
    u <- n + 500 * n 
    i <- seq(along = t)
    s1 <- sum(t[i])
    s2 <- sum((i-1)*(p-r)*t[i])
    rhsN <- n * s1
    rho <- 1/(p-r)
    
    while ((u-s)>0.01){
      N0 <- (s+u)/2;  
      lhs <- sum( 1 / (N0 - (i - 1) / rho ))
      rhsD <- sum( ( N0 - (i - 1) / rho ) * t[i] )
      dif <- lhs - rhsN / rhsD;
      if (dif < 0){
        u <- N0
      } else {
        s <- N0
      }
    }
    
    return(ceiling(N0))
  }
  
  # Search for N0hat
  N0 <- N0hat(t, p, r)
  phi <- phihat(N0, t, p, r)
  
  return(list(N0 = N0, phi = phi))
}


# Function for plotting the mean value function of the Jelinski-Moranda model

"jelinski.moranda.plot" <- function(N0, phi, p, r, t, xlab = "Time", 
                                    ylab = "Cumulated failures and estimated mean value function", main = NULL) 
{
  plot(cumsum(t), 1:length(t), type = "o", pch = 3, cex = 0.8, xlab = xlab,
       ylab = ylab)
  lines(cumsum(t), mvf.jelinski.moranda(N0, phi, p, r, cumsum(t)), lty = 2, type = "o", 
        pch = 4, cex = 0.8, col = "mediumblue")
  legend("bottomright", c("Data", "Jelinski-Moranda"), col = c("black", "mediumblue"), 
         pch = c(3, 4), lty = c(1, 2)) 
  
  if(!is.null(main))
  {
    title(main)
  }  
}

# Mean Value Funktion Jelinski-Moranda model

"mvf.jelinski.moranda" <- function(N0, phi, p, r, t) 
{
  if(length(N0) != 1 || length(phi) != 1 || length(p) != 1 || length(r) != 1)
  {
    stop("N0, phi, p and r should have length 1")
  }
  
  return(N0 *(1 - exp(-phi * ( p - r ) * t)))
}