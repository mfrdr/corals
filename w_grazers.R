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
e <- 0.2     # ranges between 0 and 1
beta <- 10
f <- 0.2 # from Blackwood et al 2010


# Initial conditions
R <- 0.2
C <- 0.4
M <- 0.3
G <- 1
B <- 1 - R - C - M

times <- seq(0,1000,1)

derivative = function(t,y,param){
  R = y[1]
  C = y[2]
  M = y[3]
  G = y[4]
  
  dR = k*(w*C + e)*(1 - R - C - M) - a*R - s*R*M - n*R
  dC = a*R + g*C*(1 - R - C - M) - b*C*M - m*C
  dM = s*M*(1 - C - M) + b*C*M - z*M*G - h*M
  dG = (z*M*G*(1 - (G/(beta*(C)))) - f*G)
  
  list(c(dR,dC,dM,dG))
}

init <- c(R,C,M,G)
all.out <- ode(init, times, derivative, parms=NULL, method = vode )

# blank space
B <- 1 - all.out[,2] - all.out[,3] - all.out[,4]

#RESULTS upto the timestep at which the integration failed
maxtime <- times[1:length(all.out[,1])]

plot(x = maxtime, y = all.out[,2], type = "l", main = "Seabed coverage",
     ylab = "Proportion of occupied seabed", xlab = "Time (years)", ylim=c(0,1),
     pch = 3, lwd = 4, col = '#F6C64A', cex.main = 1.5, cex.axis=1.3, cex.lab=1.3,
     panel.first = grid())

lines(x = maxtime, y = all.out[,3], type = "l",
      pch = 3, lwd = 4, col = '#E95951')
lines(x = maxtime, y = all.out[,4], type = "l",
      pch = 3, lwd = 4, col = '#39656E')
lines(x = maxtime, y = all.out[,5], type = "l",
      pch = 3, lwd = 4, col = "black", lty = 2)
lines(x = maxtime, y = B, type = "l",
      pch = 3, lwd = 4, col = "green")


legend(x = "topright",
       legend = c("Coral recruits", "Coral adults", "Macroalgae", "Parrotfish"),
       lty = 1,
       col = c('#F6C64A', '#E95951', "#39656E", "black"),
       lwd = 4)

# PHASE PLOT
plot(all.out[,4]*100, all.out[,3]*100, type = "l")
points(all.out[1,4]*100, all.out[1,3]*100, col = "green")

# #RESULTS
# plot(x = times, y = all.out[,2], type = "l", main = "Seabed coverage",
#      ylab = "Proportion of occupied seabed", xlab = "Time (years)", ylim=c(0,1),
#      pch = 3, lwd = 4, col = '#F6C64A', cex.main = 1.5, cex.axis=1.3, cex.lab=1.3,
#      panel.first = grid())
# 
# lines(x = times, y = all.out[,3], type = "l",
#       pch = 3, lwd = 4, col = '#E95951')
# lines(x = times, y = all.out[,4], type = "l",
#       pch = 3, lwd = 4, col = '#39656E')
# lines(x = times, y = all.out[,5], type = "l",
#       pch = 3, lwd = 4, col = "black", lty = 2)
# 
# 
# legend(x = "topright",
#        legend = c("Coral recruits", "Coral adults", "Macroalgae", "Parrotfish"),
#        lty = 1,
#        col = c('#F6C64A', '#E95951', "#39656E", "black"),
#        lwd = 4)


##################################################################################################################################
# CORAL COVERAGE EQUILIBRIUM AND INSTABILITY ANALYSIS

times <- seq(0,150,0.1)
e <- 0.2
num_scenarios <- 1000
c_values <- c(0.10, 0.15, 0.166, 0.20, 0.25)

R <- 0.2
M <- 0.3
G <- 1
init_values <- list(c(R, 0.10, M, G ), 
                    c(R, 0.15, M, G ), 
                    c(R, 0.1777,M , G ), 
                    c(R, 0.20, M, G ), 
                    c(R, 0.9, M, G ))

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
# PHASE PLOT (CORAL vs MACROALGAE) e=0.2
e <- 0.2
num_scenarios <- 1000
init_values <- list()

# Random scenarios
for (i in 1:num_scenarios) {
  init <- runif(4)
  init <- init / sum(init)
  init <- init[1:3]
  init[4] <- 1
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


# ISOCLINE ZERO NET GROWTH M
c_values <- seq(0,0.9,0.01) # initial values for coral
m_zg <- numeric(length(c_values))
for (i in seq_along(c_values)){
  c <- c_values[i]
  m_zg[i] = 1 - c + b*c/s - (z/s)*G - h/s
  print(m_zg[i])
}

# ISOCLINE ZERO NET GROWTH C
# m_values <- seq(0,0.9,0.01)
# c_zg <- numeric(length(c_values))
# for (i in seq_along(c_values)){
#   m <- m_values[i]
#   m_zg[i] = 
# }

# Plot phase plane: Corals vs. Macroalgae [color coded]
plot(corals.out[[1]][,4]*100, corals.out[[1]][,3]*100, type="l", xlab="Macroalgae", ylab="Corals", main="Phase plane: Corals vs. Macroalgae (e=0.2)", 
     ylim=c(0,90), xlim=c(0,90), pch = 3, lwd = 1, cex.main = 1.5, cex.axis=1.3, cex.lab=1.3, col="#E95951", axes=0)
lines(algae.out[[1]][,4]*100, algae.out[[1]][,3]*100, type="l", pch = 3, lwd = 1, col="#39656E")

len <- max(length(corals.out),length(algae.out))
for (i in 2:len) {
  if (i < length(algae.out)) {
    lines(algae.out[[i]][,4]*100, algae.out[[i]][,3]*100, type="l", pch = 3, lwd = 1, col="#39656E")
  }
  if (i < length(corals.out)){
    lines(corals.out[[i]][,4]*100, corals.out[[i]][,3]*100, type="l", pch = 3, lwd = 1, col="#E95951")
  }
}
lines(y=c_values*100, x=m_zg*100, type="l", pch = 3, lwd = 3, col="black")
# lines(y=c_zg*100, x=m_values*100, type="l", pch = 3, lwd = 3, col="grey")

box(lwd=3)
axis(side = 2, at = seq(0, 100, 10), lwd = 2, xaxs = "i", las = 1, cex.axis = 1.3)
axis(side = 1, at = seq(0, 100, 10), lwd = 2, yaxs = "i", cex.axis = 1.3)

