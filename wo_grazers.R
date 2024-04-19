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

derivative = function(t,y,param){
  R = y[1]
  C = y[2]
  M = y[3]
  
  dR = k*(w*C + e)*(1 - R - C - M) - a*R - s*R*M - n*R
  dC = a*R + g*C*(1 - R - C - M) - b*C*M - m*C
  dM = s*M*(1 - C - M) + b*C*M - z*M*((o*C)/(1+o*C)) - h*M
  
  list(c(dR,dC,dM))
}

init <- c(R,C,M)
all.out <- ode(init, times, derivative, parms=NULL)

all.out[,2]


# RESULTS
plot(x = times, y = all.out[,2], type = "l", main = "Seabed coverage",
     ylab = "Proportion of occupied seabed", xlab = "Time (years)", ylim=c(0,1),
     pch = 3, lwd = 4, col = '#F6C64A', cex.main = 1.5, cex.axis=1.3, cex.lab=1.3)

lines(x = times, y = all.out[,3], type = "l",
     pch = 3, lwd = 4, col = '#E95951')
lines(x = times, y = all.out[,4], type = "l",
      pch = 3, lwd = 4, col = '#39656E')

legend(x = "topright",         
       legend = c("Coral recruits", "Coral adults", "Macroalgae"), 
       lty = 1,           
       col = c('#F6C64A', '#E95951', "#39656E"),
       lwd = 4) 
grid()


##################################################################################################################################

initcond.out <- list()
num_scenarios <- 1000

# Random scenario
for (i in 1:(num_scenarios)) {
  init <- runif(4)
  init <- init / sum(init)
  init <- init[1:3]
  initcond.out <- append(initcond.out, list(ode(init, times, derivative, parms=NULL)))
}

# Color coding
corals.out <- list()
algae.out <- list()
for (i in 1:num_scenarios) {
  init <- runif(4)
  init <- init / sum(init)
  init <- init[1:3]
  init_values <- append(init_values, list(c(R,C,M)))
  temp.out <- ode(init, times, derivative, parms=NULL)
  
  if (temp.out[999,3] > 0.8){
    corals.out <- append(corals.out, list(temp.out))
  }
  else {
    algae.out <- append(algae.out, list(temp.out))
  }
}

#Not so random scenario
values <- seq(0,1,0.01)
init_values <- list()
for (c in 1:100){
  C <- values[c]
  
  for (m in 1:(100*(1-C))) {
    M <- values[m] 
    R <- runif(1, 0, 1-C-M)
    
    init_values <- append(init_values, list(c(R,C,M)))
  }
}


corals.out <- list()
algae.out <- list()
init <- c(0,0,0)

for (i in 0:length(init_values)){
  init <- init_values[i]
  print(init)
  
}

len = length(init_values)
for (i in 0:len) {
  init <- init_values[i]
  print(init_values[i])
  temp.out <- ode(init, times, derivative, parms=NULL)
  if (temp.out[999,3] >= 0.8){
    corals.out <- append(corals.out, list(temp.out))
  }
  else {
    algae.out <- append(algae.out, list(temp.out))
  }
}




# Plot phase plane: Corals vs. Macroalgae [one matrix]
plot(initcond.out[[1]][,4], initcond.out[[1]][,3], type="l", xlab="Macroalgae", ylab="Corals", main="Phase plane: Corals vs. Macroalgae", 
     ylim=c(0,0.9), xlim=c(0,0.9), pch = 3, lwd = 1, cex.main = 1.5, cex.axis=1.3, cex.lab=1.3, col="#39656E")
for (i in 2:num_scenarios) {
  lines(initcond.out[[i]][,4], initcond.out[[i]][,3], type="l", pch = 3, lwd = 1, col="#39656E")
}


# Plot phase plane: Corals vs. Macroalgae [two matrix]
plot(corals.out[[1]][,4], corals.out[[1]][,3], type="l", xlab="Macroalgae", ylab="Corals", main="Phase plane: Corals vs. Macroalgae", 
     ylim=c(0,0.9), xlim=c(0,0.9), pch = 3, lwd = 1, cex.main = 1.5, cex.axis=1.3, cex.lab=1.3, col="#E95951")
lines(algae.out[[1]][,4], algae.out[[1]][,3], type="l", pch = 3, lwd = 1, col="#39656E")

len <- max(length(corals.out),length(algae.out))
for (i in 2:len) {
  if (i < length(corals.out)){
    lines(corals.out[[i]][,4], corals.out[[i]][,3], type="l", pch = 3, lwd = 1, col="#E95951")
  }
  if (i < length(algae.out)) {
    lines(algae.out[[i]][,4], algae.out[[i]][,3], type="l", pch = 3, lwd = 1, col="#39656E")
  }
}


##################################################################################################################################

e_values <- seq(0,1,10)
e_values <- 

# Random scenario
for (i in 1:(num_scenarios)) {
  init <- runif(4)
  init <- init / sum(init)
  init <- init[1:3]
  initcond.out <- append(initcond.out, list(ode(init, times, derivative, parms=NULL)))
}


# Coral cover and external supply
plot(algae.out[[1]][,4], algae.out[[1]][,3], type="l", xlab="Macroalgae", ylab="Corals", main="Phase plane: Corals vs. Macroalgae", 
     ylim=c(0,1), xlim=c(0,1), pch = 3, lwd = 1, cex.main = 1.5, cex.axis=1.3, cex.lab=1.3, col="#39656E")
for (i in 2:num_scenarios) {
  lines(algae.out[[i]][,4], algae.out[[i]][,3], type="l", pch = 3, lwd = 1, col="#39656E")
}



