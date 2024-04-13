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
C <- 0.1
M <- 0.3

times <- seq(0,1000,1)
num_scenarios <- 1000

derivative = function(t,y,param){
  R = y[1]
  C = y[2]
  M = y[3]
  
  dR = k*(w*C + e)*(1 - R - C - M) - a*R - s*R*M - n*R
  dC = a*R + g*C*(1 - R - C - M) - b*C*M - m*C
  dM = s*M*(1 - C - M) + b*C*M - z*M*((o*C)/(1+o*C)) - h*M
  
  list(c(dR,dC,dM))
}

algae.out <- list()

# Random scenario
for (i in 1:(num_scenarios)) {
  init <- runif(4)
  init <- init / sum(init)
  init <- init[1:3]
  algae.out <- append(algae.out, list(ode(init, times, derivative, parms=NULL)))
}


# RESULTS
plot(x = times, y = algae.out[[7]][, 2], type = "l", main = "Seabed coverage",
     ylab = "Proportion of occupied seabed", xlab = "Time (years)", ylim=c(0,1),
     pch = 3, lwd = 4, col = '#F6C64A', cex.main = 1.5, cex.axis=1.3, cex.lab=1.3)

lines(x = times, y = algae.out[[7]][,3], type = "l",
     pch = 3, lwd = 4, col = '#E95951')
lines(x = times, y = algae.out[[7]][,4], type = "l",
      pch = 3, lwd = 4, col = '#39656E')

legend(x = "topright",         
       legend = c("Coral recruits", "Coral adults", "Macroalgae"), 
       lty = 1,           
       col = c('#F6C64A', '#E95951', "#39656E"),
       lwd = 4) 
grid()

# Plot phase plane: Corals vs. Macroalgae
plot(algae.out[[1]][,4], algae.out[[1]][,3], type="l", xlab="Macroalgae", ylab="Corals", main="Phase plane: Corals vs. Macroalgae", 
     ylim=c(0,1), xlim=c(0,1), pch = 3, lwd = 1, cex.main = 1.5, cex.axis=1.3, cex.lab=1.3, col="#39656E")
for (i in 2:num_scenarios) {
  lines(algae.out[[i]][,4], algae.out[[i]][,3], type="l", pch = 3, lwd = 1, col="#39656E")
}



