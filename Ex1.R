
# TARGET FUNCTION ---------------------------------------------------------

# Define the target function

# Doppler function scaled in [0,1]
doppler.fun <-  function(x){
  "This function returns the Doppler function in [0,1]:
  x : function argument in [0,1]"
  out = sqrt(x*(1 - x))*sin( (2.1*pi)/(x + 0.05) )
  return(out)
}

# Doppler function scaled in [-1,1]
h_x = function(x){
  "This function returns the Doppler function scaled in [-1,1]:
  x : function argument in [-1,1]"
  return(doppler.fun(0.5*(x+1)))
}
# Plot the scaled Doppler
par(mfrow=c(1,2))
curve(doppler.fun, ylab = expression(g(x)), main = 'Doppler function in [0,1]', col = 'blue')
curve(h_x(x), xlim = c(-1,1), ylab = expression(h(x)), main = 'Doppler function in [-1,1]', col = 'blue')


# LEGENDRE BASIS DEFINITION -----------------------------------------------

# Set up the function to construct recursively the Legendre basis
leg.basis = function(x, j) {
  "The function returns the list of the orthogonal and orthonormal basis vectors."
  if (j==0)
    p=1
  else if (j==1)
    p=x
  else {
    # Define the Legendre polynomial p_(j)
    p = x
    # Define the Legendre polynomial p_(j-1)
    p_1 = 1
    for (k in 1:(j-1)){
      # Compute the Legendre polynomial p_(j+1)
      p_next = ((2 * k + 1) * x * p - k * p_1) / (k+1)
      p_1 = p
      p = p_next
    }
  }
  # Orthonormalize the basis vectors
  q = sqrt((2*j+1)/2) * p
  return (list(p,q))
}
# Vectorize the function
leg.basis.vec = Vectorize(leg.basis)

# Plot first five Legendre polynomials
library(RColorBrewer)
mycol = brewer.pal(5, "Blues")
# Declare the lwd of the plot
plot_lty=c(1,2,3,4,1)
# Plot the first file Legendre the polynomials 
for (j in 0 : 4){
  curve(leg.basis.vec(x,j)[2,], n = 501,  ylim = c(-3,3), xlim=c(-1,1), add=(j!=0),
        col = mycol[j+1], lty=plot_lty[j+1], lwd = 2, xlab = 'x', ylab = expression(phi[j](x)),
        main = 'Legendre polynomials in [-1,1]')
}
abline(v=0, h=0)
legend('bottomright', c("Q_0(x)","Q_1(x)","Q_2(x)","Q_3(x)","Q_4(x)"), lty=c(1,2,3,4,1), col=mycol,
       bty="n", cex=0.75, inset=0.05, horiz=F)


# BASIS ORTHONORMALITY ----------------------------------------------------

grade=5
# Declare a NA matrix
orth_leg.basis = matrix(nrow=grade, ncol=grade)
# Fill in the matrix : m(j,k) = p_(j)*p_(k)
for (j in 0 : (grade-1)){
  for (k in 0 : (grade-1)){
    # Compute the integral in [-1,1] to verify that the basis elements satisfy the requirements
    orth_leg.basis[(j+1),(k+1)] = round( integrate(function(x, j1, j2) 
      unlist(leg.basis.vec(x, j1)[2,]) * 
        unlist(leg.basis.vec(x,j2)[2,]),
      lower = -1, upper = 1, j1 = k, j2 = j)$value, 3)
  }
}
# Check the orthogonality, you expect a diagonal matrix where d(i,i) = 1
orth_leg.basis


# FOURIER COEFFICIENTS ----------------------------------------------------

# Define the coefficients function
beta_integral = function(j) {
  "This function return the vector of Fourier coefficients. It's lenght depend on:
  - j : length of the vector."
  beta=c()
  # Compute Beta as the integral in [-1,1] of the inner product btw g and phi
  for (i in 0:(j-1)) {
    beta = c(beta, integrate(function(x) h_x(x) * unlist(leg.basis.vec(x,i)[2,]),
                             lower=-1, upper=1)$value) 
  }
  return(beta)
}

# Evaluate the first 200 Fourier coefficients of h in the Legendre basis
beta = beta_integral(200)
```

# Plot beta 
par(mfrow=c(1,1))
plot(seq(1,200, by=1),beta, type="h", main = 'Fourier Coefficients : j in {0, 1, .., 200}', xlab = 'j', ylab = expression(beta[j]))


# LINEAR APPROXIMATION -----------------------------------------------------------

# Built two linear approximations for the target function.
# Set up the function to compute the approximation of the target function
h_j = function (x, j, b){
  "This function returns the value of the approximation function given:
  - x : the value in which I compute the function;
  - j : the j-th term for the non-linear approximation;
  - b : the vector of Beta's related to j."
  out = 0
  for (k in 0 : (j-1)){
    out = out + b[k+1] * unlist(leg.basis.vec(x, k)[2,])
  }
  return (out)
}


#Plot the linear approximation on the target function
par(mfrow=c(1,1))
curve(h_x(x), xlim = c(-1,1), ylab = expression(h(x)), main = 'Doppler function in [-1,1] :: Approximate g, J = 50', col = 'green', lwd = 5)
curve(h_j(x, 50, beta[1:50]), xlim = c(-1,1), add = T, lty = 2, col = 'black')
legend('bottomright', legend = c('Doppler', 'g_J = 50'), col = c('green', 'black'), lty = c(1, 2), text.width = 0.5)

curve(h_x(x), xlim = c(-1,1), ylab = expression(h(x)), main = 'Doppler function in [-1,1] :: Approximate g, J = 50', col = 'pink', lwd = 5)
curve(h_j(x, 100, beta[1:100]), xlim = c(-1,1), add = T, lty = 2, col = 'black')
legend('bottomright', legend = c('Doppler', 'g_J = 100'), col = c('pink', 'black'), lty = c(1, 2), text.width = 0.5)
```

# EVALUATE THE LINEAR APPROXIMATION ----------------------------------------------

# The evaluation consists in compute integral in [-1,1] of the squared differences
# btw the target and the approximated function.
# One where J = 50
integrate(function(x) (h_x(x) - h_j(x, 50, beta))^2, lower=-1, upper=1)$value
# The other for J = 100
integrate(function(x) (h_x(x) - h_j(x, 100, beta))^2, lower=-1, upper=1)$value


# NON-LINEAR APPROXIMATION ------------------------------------------------

# In order to compute the greedy approximation define the following function
h_j_star = function (x, j, b){
  "This function returns the value of the approximation function given:
  - x : the value in which I compute the function;
  - j : the j-th term for the non-linear approximation;
  - b : the vector of Beta's related to j."
  out = 0
  # Put equal to zero all the beta's smaller than the j-th term.
  b[abs(b) < sort(abs(b), decreasing = TRUE)[j]] = 0
  for (k in 0:(length(b)-1)){
    out = out + b[k+1] * unlist(leg.basis.vec(x,k)[2,])
  }
  return (out)
}


# Plot of the non-linear approximation on the target functino
par(mfrow=c(1,1))
curve(h_x(x), xlim = c(-1,1), ylab = expression(h(x)), main = 'Doppler function in [-1,1] :: Approximate g*, J = 50', col = 'green', lwd = 5)
curve(h_j_star(x, 50, beta[1:50]), xlim = c(-1,1), add = T, lty = 2, col = 'black')
legend('bottomright', legend = c('Doppler', 'g*_J = 50'), col = c('green', 'black'), lty = c(1, 2), text.width = 0.5)

curve(h_x(x), xlim = c(-1,1), ylab = expression(h(x)), main = 'Doppler function in [-1,1] :: Approximate g*, J = 50', col = 'pink', lwd = 5)
curve(h_j_star(x, 100, beta[1:100]), xlim = c(-1,1), add = T, lty = 2, col = 'black')
legend('bottomright', legend = c('Doppler', 'g*_J = 100'), col = c('pink', 'black'), lty = c(1, 2), text.width = 0.5)


# EVALUATE THE NON-LINEAR APPROXIMATION -----------------------------------

# Evaluate the approximation for j = 50
integrate(function(x) (h_x(x) - h_j_star(x, 50, beta))^2, lower=-1, upper=1)$value
# Do the same for j = 100
integrate(function(x) (h_x(x) - h_j_star(x, 100, beta))^2, lower=-1, upper=1)$value
