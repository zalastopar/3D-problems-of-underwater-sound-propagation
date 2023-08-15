library(sdetorus)
library(matlib)
library(cmna)
library(tidyr)
library(plotly)

setwd("~/AMAT/IZMENJAVA/nov_poskus")

old <- Sys.time()

library(readr)
manjsiH <- read_csv("h.csv")

PML <- TRUE

# Nx <- nrow(manjsiH) - 1
# dx <- 1;
# x_points <- seq(0, Nx, by = dx);
# x_steps <- length(x_points);
# 
# Ny <- ncol(manjsiH) - 1;
# dy <- 1;
# y_lower <- - 500;
# L = 100; #------------------------------------------------------- L
# y_points <- seq(y_lower, y_lower + Ny, by = dy);
# y_steps <- length(y_points);
# 
# h_0 <- 50
# alpha <- 2
# c <- 1500
# f <- 100
# omega <- 2*pi*f
# sigma <- 300 #-------------------------------------------------- sigma
# 
# j <- 5
# H <- manjsiH
# #H00 <- h_0 + tan(alpha*pi/180);
# 
# H00 <- H[,1]
# umesni1 <- omega**2/c**2 - (pi*j/H00)**2
# K00j <- sqrt(apply(umesni1, 1, as.complex));
# #K00j <- sqrt(omega**2/c**2 - (pi*j/H00)**2);
# 
# umesni <- omega**2/c**2 - (pi*j)**2/(H*H)
# Kj <- t(sqrt(apply(umesni, 1, as.complex)))
# K00j <- Kj[,1]


#-------------------------------------------------------------------------------------------------------
Nx <- 1000;
dx <- 1;
x_points <- seq(0, Nx, by = dx);
x_steps <- length(x_points);

Ny <- 1000;
dy <- 1;
y_lower <- - 500;
L = 100; #------------------------------------------------------- L
y_points <- seq(y_lower, y_lower + Ny, by = dy);
y_steps <- length(y_points);

h_0 <- 50
alpha <- 2
c <- 1500
f <- 100
omega <- 2*pi*f
sigma <- 300 #-------------------------------------------------- sigma


#j <- 5

H00 <- h_0 + tan(alpha*pi/180)*0;

A <- matrix(0, y_steps, x_steps);

for (j in 1:5){
  
  
  
K00j <- sqrt(omega**2/c**2 - (pi*j/H00)**2);


function_H <- function(x, y, h_0, alpha) {
  solution <- h_0 + tan(alpha *pi/180)*y
  return(solution)}



# xy_combinatioons <- expand.grid(x_points, y_points)
# 
# function_H <- function(combo, h_0, alpha) {
# 
#   combo[,3] <- h_0 + tan(alpha *pi/180)*combo[,2]
# 
#   return(combo)}
# 
# 
# new_H <- function_H(xy_combinatioons, h_0, alpha)
# H <- spread(new_H, Var2, V3)[,-1]


H1 <- h_0 + tan(alpha *pi/180)*y_points #function_H(x_points, y_points, h_0, alpha);

H <-t(matrix(H1, nrow=x_steps, ncol=length(H1), byrow=TRUE));


Kj1 <- sqrt(as.complex(omega**2/c**2 - (pi*j)**2/(H1*H1)))
#Kj <- matrix(Kj1, nrow=x_steps, ncol=length(Kj1), byrow=TRUE);
umesni <- omega**2/c**2 - (pi*j)**2/(H*H)
Kj <- sqrt(apply(umesni, 1, as.complex))



K <- Kj**2 - K00j**2;
funk <- K/(2)


# Q_pml ------------------------------------------------------
make_Q_pml <- function(y_points, y_steps, is_PML, dy, Ny, y_lower, L, sigma) {
  Q_pml <- matrix(0, y_steps, 1);
  
  if (is_PML) {
    Q_pml <- y_points
    Q_pml[Q_pml >= y_lower + L & Q_pml <= y_lower + Ny - L] <- 0
    Q_pml[1: (floor(L/dy))] <-  rev( - sigma*1i*c(0:(floor(L/dy)-1))**3/((floor(L/dy) - 1)**3))
    Q_pml[-(length(Q_pml) - (floor(L/dy)): length(Q_pml))] <- - sigma*1i*c(0:(floor(L/dy)-1))**3/((floor(L/dy) - 1)**3)
  }
  
  return(Q_pml)
  
}

# Q_pml <- matrix(0, y_steps, 1);
# 
# if (PML) {
#   Q_pml <- y_points
#   # y_l <- y_lower + L
#   # y_r <- y_lower + Ny - L
#   Q_pml[Q_pml >= y_lower + L & Q_pml <= y_lower + Ny - L] <- 0
#   Q_pml[1: (floor(L/dy))] <-  rev( - sigma*1i*c(0:(floor(L/dy)-1))**3/((floor(L/dy) - 1)**3))
#   Q_pml[-(length(Q_pml) - (floor(L/dy)): length(Q_pml))] <- - sigma*1i*c(0:(floor(L/dy)-1))**3/((floor(L/dy) - 1)**3)
# }

Q_pml <- make_Q_pml(y_points, y_steps, PML, dy, Ny, y_lower, L, sigma)


# (1i*K00j/(1i*K00j + sigma*(y_l - Q_pml[1: floor(L/dy)])**2))**2
# (1i*K00j/(1i*K00j + sigma*(Q_pml[-(length(Q_pml)- floor(L/dy) : length(Q_pml))] - y_r)**2))**2

#-------------------------------------------------------------


# matrix Q -----------------------------------------------------
# B <- matrix(0, y_steps, 1);
# 
# 
# 
# diag(Q) <- 2*1i*K00j/dx - B/dy**2
# new_diag_Q <- matrix(0, y_steps, y_steps);
# diag(new_diag_Q) <- B/(2*dy**2)
# 
# lower_Q <- rbind(rep(0, y_steps-1), new_diag_Q[2:length(new_diag_Q)])

# Q <- matrix(0, y_steps, y_steps);
# diag(Q) <- 2*1i*K00j/dx - (1 + Q_pml)/dy**2
# 
# new_diag_Q <- matrix(0, y_steps, y_steps);
# diag(new_diag_Q) <- (1 + Q_pml)/(2*dy**2)
# 
# lower_diag_Q <- rbind(replicate(y_steps, 0), cbind(new_diag_Q[1:y_steps-1, 1:y_steps-1], replicate(y_steps-1, 0))) # od zacetka
# upper_diag_Q <- rbind(cbind(replicate(y_steps-1, 0), new_diag_Q[2:y_steps, 2:y_steps]), replicate(y_steps, 0))
# Q <- Q + lower_diag_Q + upper_diag_Q


n <- (1 + Q_pml)/(2*dy**2)
# lower_diagonal <- n
# lower_diagonal[1] <- 0
# upper_diagonal <- n
# upper_diagonal[length(upper_diagonal)] <- 0

lower_diagonal1 <- n[2:length(n)]
upper_diagonal1 <- n[1:length(n)-1]



# matrix P -----------------------------------------------------


P <- matrix(0, y_steps, y_steps);
diag(P) <- 2*1i*K00j/dx + (1 + Q_pml)/dy**2;

new_diag_P <- matrix(0, y_steps, y_steps);
diag(new_diag_P) <- - (1 + Q_pml)/(2*dy**2)

lower_diag_P <- rbind(replicate(y_steps, 0), cbind(new_diag_P[1:y_steps-1, 1:y_steps-1], replicate(y_steps-1, 0))) # od zacetka
upper_diag_P <- rbind(cbind(replicate(y_steps-1, 0), new_diag_P[2:y_steps, 2:y_steps]), replicate(y_steps, 0))
P <- P + lower_diag_P + upper_diag_P


# initial condition -------------------------------------------
a0 <- 1/(2 * sqrt(pi))* exp(as.numeric(-y_points**2 * K00j**2));
a0[!is.finite(a0)] <- 0


B <- matrix(0, y_steps, x_steps);
B[, 1] <- a0



make_phi_j <- function(j, z, H) {
  new <- sqrt(2/H)*sin(pi*j*z/H)
  return(new)
}

# create system 


for (x in 2:x_steps){ 
  umesna <- matrix(0, y_steps, y_steps);
  diag(umesna) <- funk[x,];
  
  P_update <- P - umesna
  d <- P_update%*%B[,x-1]
  
  nasa_diag <- 2*1i*K00j/dx - (1 + Q_pml)/dy**2 + funk[x,]
  

  
  resitev <- tridiagmatrix(lower_diagonal1, nasa_diag, upper_diagonal1, d)   
  B[, x] <- resitev
  
}
phi_j <- make_phi_j(j, z_source, H)
A <- A + B*matrix(exp(K00j*1i*x_points), nrow=x_steps, ncol=length(x_points), byrow=TRUE)

}
new <- Sys.time()
print(new-old)


manjsa_risba <- function(m){
  m[m >= min(m)+ 0.04] <- min(m)+ 0.04
  return(m)
}

to_decibels <- function(m) {
  s <- 20*log10(abs(m))[,-c(1)]

  return(s)
}

#image(y_points, x_points, to_decibels(A), col = rainbow(100))
plot_ly(z = to_decibels(A), y = y_points, x = x_points, colorscale = "Jet", 
        zauto = F,zmin = max(max(to_decibels(A)))-50, zmax = max(max(to_decibels(A)))-10,type = "heatmap") %>% layout(yaxis = list(autorange = "reversed"))






