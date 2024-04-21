library(deSolve)
library(fields)
library(colorspace)


  
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
e <- 0.2     # ranges between 0 and 1


# Initial conditions
R <- 0.2
C <- 0.3
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
times <- seq(0,50,0.1)
R <- 0.1
M <- 0.5
init_values <- list(c(R, 0.10, M ), 
                    c(R, 0.15, M ), 
                    c(R, 0.165,M ), 
                    c(R, 0.20, M ), 
                    c(R, 0.25, M ))


corals.out <- list()
for (i in 1:length(init_values)) {
  temp.out <- ode(init_values[[i]], times, derivative, parms=NULL)
  corals.out <- append(corals.out, list(temp.out))
  print(tail(temp.out))
}



#RESULTS
coral_palette <- colorRampPalette(c("#FFD1DC", "#FFB6C1", "#FF98A8", "#FF7A8F", "#FF5C75", "#FF3E5C", "#FF2042", "#FF001E"))
plot(x=times, y=corals.out[[1]][,3]*100, type="l", main="Coral coverage",
     ylab="Coral coverage (%)", xlab="Time (yrs)", ylim=c(0, 100),
     pch=1, lwd=5, col=coral_palette(length(init_values)), cex.main=1.5, cex.axis=1.3, cex.lab=1.3,
     panel.first=grid(lty=1, col='lightgrey'), axes=0)

for (i in 2:length(init_values)) {
  lines(x = times, y = corals.out[[i]][,3]*100, type = "l",
        pch = 1, lwd = 5, col = coral_palette(length(init_values))[i])
}


for (x in seq(0, 50, 1)) {symbols(x, 15.65, circles = 0.5, inches = FALSE, add = TRUE, fg = "black", lty = NULL, lwd = 1)}
for (x in seq(0, 50, 1)) {symbols(x, 86.57404, circles = 0.5, inches = FALSE, add = TRUE, fg = "black", bg = "black", lty = NULL, lwd = 2)}
for (x in seq(0, 50, 1)) {symbols(x, 0, circles = 0.5, inches = FALSE, add = TRUE, fg = "black", bg = "black", lty = NULL, lwd = 2)}
box(lwd=3)
axis(side = 2, at = seq(0, 100, 10), lwd = 2, xaxs = "i", las = 1, cex.axis = 1.3)
axis(side = 1, at = seq(0, 50, 10), lwd = 2, yaxs = "i", cex.axis = 1.3)

##################################################################################################################################

num_scenarios <- 1000
init_values <- list()

# Random scenario
for (i in 1:num_scenarios) {
  init <- runif(4)
  init <- init / sum(init)
  init <- init[1:3]
  init_values <- append(init_values, list(init))
}

# #Not so random scenario
# values <- seq(0,1,0.1)
# init_values <- list()
# for (c in 1:10){
#   C <- values[c]
#   
#   for (m in 1:(10*(1-C))) {
#     M <- values[m] 
#     R <- runif(1, 0, 1-C-M)
#     
#     init <- c(R,C,M)
#     
#     init_values <- append(init_values, list(init))
#   }
# }
# 
# num_scenarios_per_value <- 20  # Number of scenarios per value of C or M
# C_values <- seq(0.1, 0.9, 0.1)  # Values of C from 0.1 to 0.9 with spacing 0.1
# M_values <- seq(0.1, 0.9, 0.1)  # Values of M from 0.1 to 0.9 with spacing 0.1
# 
# init_values <- list()
# 
# for (C in C_values) {
#   for (M in M_values) {
#     if (C+M >= 1){
#       break
#     }
#     for (i in 1:num_scenarios_per_value) {
#       # Calculate the remaining free space
#       R_free_space <- 1 - C - M
# 
#       # Randomly assign the remaining free space to R
#       R <- runif(1, 0, R_free_space)
# 
#       # Store the initial conditions
#       init_values <- append(init_values, list(c(R, C, M)))
#     }
#   }
# }
# 


corals.out <- list()
algae.out <- list()
initcond.out <- list()

for (i in 1:length(init_values)) {
  
  temp.out <- ode(init_values[[i]], times, derivative, parms=NULL)
  initcond.out <- append(initcond.out, temp.out)
  
  if (tail(temp.out[,3], 1) >= 0.8) {
    corals.out <- append(corals.out, list(temp.out))
  } else {
    algae.out <- append(algae.out, list(temp.out))
  }
}




# Plot phase plane: Corals vs. Macroalgae [one matrix]
# plot(initcond.out[[1]][,4], initcond.out[[1]][,3], type="l", xlab="Macroalgae", ylab="Corals", main="Phase plane: Corals vs. Macroalgae", 
#      ylim=c(0,0.9), xlim=c(0,0.9), pch = 3, lwd = 1, cex.main = 1.5, cex.axis=1.3, cex.lab=1.3, col="#39656E")
# for (i in 2:num_scenarios) {
#   lines(initcond.out[[i]][,4], initcond.out[[i]][,3], type="l", pch = 3, lwd = 1, col="#39656E")
# }


# Plot phase plane: Corals vs. Macroalgae [color coded]
plot(corals.out[[1]][,4]*100, corals.out[[1]][,3]*100, type="l", xlab="Macroalgae", ylab="Corals", main="Phase plane: Corals vs. Macroalgae", 
     ylim=c(0,90), xlim=c(0,90), pch = 3, lwd = 1, cex.main = 1.5, cex.axis=1.3, cex.lab=1.3, col="#E95951", axes=0)
lines(algae.out[[1]][,4]*100, algae.out[[1]][,3]*100, type="l", pch = 3, lwd = 1, col="#39656E")

len <- max(length(corals.out),length(algae.out))
for (i in 2:len) {
  if (i < length(corals.out)){
    lines(corals.out[[i]][,4]*100, corals.out[[i]][,3]*100, type="l", pch = 3, lwd = 1, col="#E95951")
  }
  if (i < length(algae.out)) {
    lines(algae.out[[i]][,4]*100, algae.out[[i]][,3]*100, type="l", pch = 3, lwd = 1, col="#39656E")
  }
}
box(lwd=3)
axis(side = 2, at = seq(0, 100, 10), lwd = 2, xaxs = "i", las = 1, cex.axis = 1.3)
axis(side = 1, at = seq(0, 100, 10), lwd = 2, yaxs = "i", cex.axis = 1.3)

##################################################################################################################################
times <- seq(0, 200, 1)
e_values <- seq(0, 1, 0.05)
num_scenarios <- 10000
e_all <- numeric(length(e_values))

# Random scenario
init_values <- list()
for (i in 1:num_scenarios) {
  init <- runif(4)
  init <- init / sum(init)  # Ensure the sum is 1
  init <- init[1:3]
  init_values[[i]] <- init
}

# Solve for different e values
for (i in seq_along(e_values)) {
  tail_temp <- 0
  e <- e_values[i]
  for (j in 1:num_scenarios) {
    init <- init_values[[j]]
    temp.out <- ode(init, times, derivative, parms = NULL)
    tail_temp <- tail_temp + temp.out[length(times), 3]
  }
  e_all[i] <- tail_temp / num_scenarios
}

e_all

# Coral cover and external supply
plot(e_values, e_all, type="p", xlab="External larval supply", ylab="Coral cover", main="Effect of external larval supply", 
     ylim=c(0,1), xlim=c(0,1), pch = 16, lwd = 5, cex.main = 1.5, cex.axis=1.3, cex.lab=1.3, col="black", cex=2, panel.first = grid())



# # Initialize results vector
# e_values <- seq(0, 1, 0.1)
# init_values <- expand.grid(C = e_values, M = e_values)
# init_values <- init_values[init_values$C + init_values$M <= 1, ]
# 
# init_values
# e_all <- numeric(nrow(init_values))
# 
# # Solve for different initial conditions
# for (i in seq_len(nrow(init_values))) {
#   C <- init_values[i, "C"]
#   M <- init_values[i, "M"]
#   R <- runif(1, 0, 1 - C - M)  
#   init <- c(R, C, M)
#   result.out <- ode(init, times, derivative, parms = NULL)
#   e_all[i] <- result.out[length(times), 3]
# }
# 
# # Coral cover and external supply
# plot(e_values, e_all, type = "p", xlab = "External larval supply", ylab = "Coral cover", 
#      main = "Effect of external larval supply", ylim = c(0, 1), xlim = c(0, 1), 
#      pch = 16, lwd = 5, cex.main = 1.5, cex.axis = 1.3, cex.lab = 1.3, col = "black", cex = 2, panel.first = grid())
# 
