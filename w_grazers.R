library(deSolve)
library(fields)



# Parameters
k <- 0.25
w <- 0.4
a <- 0.2
s <- 0.95
n <- 0.5
g <- 0.4
b <- 0.2
m <- 0.05
z <- 0.64
o <- 4
h <- 0.1
e <- 0     # ranges between 0 and 1


# Initial conditions
R <- 0.2
C <- 0.3
M <- 0.3
G <- 

times <- seq(0,1000,1)

derivative = function(t,y,param){
  R = y[1]
  C = y[2]
  M = y[3]
  G = y[4]
  
  dR = k*(w*C + e)*(1 - R - C - M) - a*R - s*R*M - n*R
  dC = a*R + g*C*(1 - R - C - M) - b*C*M - m*C
  dM = s*M*(1 - C - M) + b*C*M - z*M*G - h*M
  dG = z*M*G*(1- (G/(beta*K))) - f*G 
  
  list(c(dR,dC,dM,dG))
}

init <- c(R,C,M,G)
all.out <- ode(init, times, derivative, parms=NULL)


#RESULTS
plot(x = times, y = all.out[,2], type = "l", main = "Seabed coverage",
     ylab = "Proportion of occupied seabed", xlab = "Time (years)", ylim=c(0,1),
     pch = 3, lwd = 4, col = '#F6C64A', cex.main = 1.5, cex.axis=1.3, cex.lab=1.3,
     panel.first = grid())

lines(x = times, y = all.out[,3], type = "l",
      pch = 3, lwd = 4, col = '#E95951')
lines(x = times, y = all.out[,4], type = "l",
      pch = 3, lwd = 4, col = '#39656E')

legend(x = "topright",
       legend = c("Coral recruits", "Coral adults", "Macroalgae"),
       lty = 1,
       col = c('#F6C64A', '#E95951', "#39656E"),
       lwd = 4)


##################################################################################################################################
