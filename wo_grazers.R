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
e <- 0.2       # ranges between 0 and 1


# Initial conditions
C <- 0.2
R <- 0.2
M <- 0.1

init <- c(C,R,M)
times <- seq(0,1000,1)


derivative = function(t,init,params){
  
  
  dR = k*(w*C + e)*(1 - R - C - M) - a*R - s*R*M - n*R
  dC = a*R + g*C*(1 - R - C - M) - b*C*M - m*C
  dM = s*M*(1 - C - M) + b*C*M - z*M*((o*C)/(1+o*C)) - h*M
  
  list(c(dR,dC,dM))
}


algae.out <- ode(init, times, derivative, parms=NULL)
tail(algae.out)






