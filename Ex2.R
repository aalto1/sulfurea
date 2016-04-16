
# Define a target function -----------------------------------------
library(plot3D)
library(rgl)

# Define the target function that we want to approximate
target.fun = function(x_1, x_2){
  z = outer(x_1, x_2, FUN = function (x,y) {out = x + cos(y); return(out)})
  return(z)
}
# Declare the range of values (x1,x2) used to plot the target function
x1 = seq(0, 1, by = 0.01)
x2 = seq(0, 1, by = 0.01)
f = target.fun(x1,x2)
# Plot the target function
persp3d(f,col = 'red', border = "yellow")

# Define the basis functions ----------------------------------------------

# Cosine-basis
cos.basis = function(x, j){
  return(1*(j == 0) + sqrt(2)*cos(pi*j*x)*(j > 0))
}

# Since we know that the cosine basis is an orthonormal basis for L2[0,1], 
# we affirm that the tensor basis, obtained as : phi(x1,x2) = phi(x1)*phi(x_2), 
# is orthonormal too.

# Tensor product basis

tensor.basis = function(x_1, x_2, j1, j2){
  return(cos.basis(x_1, j1) * cos.basis(x_2, j2))
}

# Fourier coefficients ----------------------------------------------------

library(cubature)

# Define the argument that we should integrate to compute the Fourier coefficients
integrande.fun = function(x, j_1, j_2){
  x_1 = x[1]
  x_2 = x[2]
  return( tensor.basis(x_1, x_2, j_1, j_2)*target.fun(x_1, x_2))
}
# Then we vectorize the function( for j_1 and j_2) in order to give vectors as input
integrande.fun.vec= Vectorize(integrande.fun, vectorize.args=c("j_1","j_2"))

# We create the Fourier coefficients matrix

# At first we define a function to compute a single Fourier coefficient

fourier.coef = function(j1, j2){
  out = adaptIntegrate(function(x, j_1, j_2) integrande.fun(x, j_1, j_2),
                       lowerLimit = c(0,0), upperLimit = c(1,1), 
                       j_1 = j1, j_2 = j2, maxEval = 50000, 
                       absError = 1e-4)$integral
  return(out)
}

# Then we vectorize the function
fourier.vec = Vectorize(fourier.coef)
# And finally create the Fourier coefficient matrix
fourier.matrix = function(j_1, j_2) {return(outer(0:j_1, 0:j_2, FUN = fourier.vec ))}

# Target function approximation ------------------------------------------

# Define the function that computes the approximation of g given x1, x2 and j1 and j2
approximate.fun = function(x.1, x.2, j1, j2){
  vec_1 = c()
  for(k in 0 : j1){
    vec_2 = c()
    for (i in 0 : j2){
      vec_2[i+1] = coeff.matrix[k+1, i+1]*tensor.basis(x.1, x.2, k, i)
    }
    vec_1[k+1] = sum(vec_2)
  }
  out = sum(vec_1)
  return(out)
}  
# Vectorize the approximation function
approximate.vec = Vectorize(approximate.fun)

# Define the matrix of values necessary to plot the approximated function
approximate.matrix = function(x_1, x_2, j_1, j_2){
  out = outer(x1, x2, FUN = approximate.vec, j1 = j_1 , j2 = j_2 )
  return(out)
}

# Evaluate three different approximation ----------------------------------
# Define the functin to compute the loss related to one approximation
loss = function(j.1,j.2){
  adaptIntegrate(function(x, j1, j2) (target.fun(x[1],x[2]) - approximate.vec(x[1],x[2],j1,j2))^2,
                 lowerLimit = c(0,0), upperLimit = c(1,1), j1 = j.1, j2 = j.2, 
                 maxEval = 50000, absError = 1e-4)$integral
}
# Define the matrix of the losses for different combination of beta's
loss.matrix = function(j1, j2){
  out = matrix(0, j1, j2)
  for (k in 0:j1-1) {
    for (i in 0:j2-1){
      out[k+1,i+1] = loss(k,i)
    }
  }
  return(out)
}

# Define error matrix
j_1=j_2=10
coeff.matrix = fourier.matrix(j_1, j_2)
err = loss.matrix(j_1,j_2)

# First
persp3d(x = x1, y = x2, f, col = 'red')
z = approximate.matrix(x1, x2, 2, 2)
persp3d(x = x1, y = x2, z,col = 'blue', add = T)
loss(3,3)

# Second
persp3d(f,col = 'red')
z = approximate.matrix(x1, x2, 9, 9)
persp3d(x = x1, y = x2, z,col = 'blue', add = T)
loss(10,10)

# Third
persp3d(x = x1, y = x2, f,col = 'red')
z = approximate.matrix(x1, x2, 1, 2)
persp3d(x = x1, y = x2, z,col = 'blue', add = T)
loss(2,3)