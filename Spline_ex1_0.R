load("wmap.RData")
library(RColorBrewer)
library(splines)

# 0.3 SPLINES ------------------------------------------------------------------

# Define a function to get the GCV respect to different values of *spar*

smoothing.splines = function(spar){
  "This function returns the GCV of a specific SS. The argument of the function is:
  - spar : is a number, in particular lambda is a monotonic function of spar."
  # Perform the SS model
  mod = smooth.spline(wmap$x, wmap$y, spar = spar)
  # Get the GCV
  gcv = mod$cv.crit  
  
  return (gcv)
}

# Define a grid of spars values
grid.spars = seq(0, 1.5, by = 0.001)

# Initialize the vector of SS's GCV
ss.gcv = rep(0, length(grid.spars))

# For each value of *spar*
for (i in 1:length(grid.spars)){
  # Get the GCV and fill in into the initialize vector
  ss.gcv[i] = smoothing.splines(grid.spars[i])
}

# Visualize the GCV values respect to the value of spar
plot(spars, ss, 'l', xlab = 'spar', ylab = 'Generalized Cross-Validation' , main = 'GCV vs Spar')

# Compute the spar.optim
spar.optim = spars[which.min(ss.gcv)]
# Plot the point on the previous graph
points(spar.optim , min(ss), col = 'red', pch = 19)
text(spar.optim, min(ss), expression(spar[opt]), pos = 3, cex = .7)

# Fit the model with the spar.optim
ss.fit = smooth.spline(wmap$x, wmap$y, spar = spar.optim)
# Plot the rusults
plot(wmap$y, pch = 19, col = rgb(0.5,0.5,.6,.3), cex = .7, main = "CMB data (WMAP)", xlab = "Multipole", ylab = "Power")
points(wmap$x, ss.fit$y, col = 'black', cex = 0.2, pch = 19)



# AGGIUNGERE TRUNCATED BASISE SPLINES?????