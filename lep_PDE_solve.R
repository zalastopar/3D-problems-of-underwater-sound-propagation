library(sdetorus)
library(matlib)
library(cmna)
library(tidyr)
library(plotly)

setwd("~/AMAT/IZMENJAVA/nov_poskus")


library(readr)


# ----------------------------------------------------------------------------------------------------
# PODATKI
# ----------------------------------------------------------------------------------------------------

# 
# H <- t(read_csv("h.csv"))
# 
# PML <- TRUE
# 
# Nx <- nrow(H) - 1
# dx <- 1;
# x_points <- seq(0, Nx, by = dx);
# x_steps <- length(x_points);
# 
# Ny <- ncol(H) - 1;
# dy <- 1;
# y_lower <- - 500;
# L = 100; #------------------------------------------------------------------------------- L
# y_points <- seq(y_lower, y_lower + Ny, by = dy);
# y_steps <- length(y_points);
# 
# 
# c <- 1500
# f <- 100 # -------------------------------------------------------------------------------f
# omega <- 2*pi*f
# sigma <- 100 #-------------------------------------------------- sigma
# 
# 
# #H00 <- H[,1];
# H00 <- as.numeric(H['0', 1])#h_0 + tan(alpha*pi/180);
# 


z_source <- 20 #------------------------------------------------------------------------- z



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
alpha <- 1
c <- 1500
f <- 100 #----------------------------------------------------------------- f
omega <- 2*pi*f
sigma <- 100 #-------------------------------------------------- sigma


#j <- 5

H00 <- h_0 + tan(alpha*pi/180)*0;
H1 <- h_0 + tan(alpha *pi/180)*y_points #function_H(x_points, y_points, h_0, alpha)#y_points*0 + h_0 ;

H <-t(matrix(H1, nrow=x_steps, ncol=length(H1), byrow=TRUE));
#---------------------------------------------------------------------------------------------------------------------
#%5555555555555555555555555555555555555

# PML <- TRUE
# #-------------------------------------------------------------------------------------------------------
# aNx <- 1000;
# adx <- 1;
# ax_points <- seq(0, Nx, by = dx);
# ax_steps <- length(x_points);
# 
# aNy <- 1000;
# ady <- 1;
# ay_lower <- - 500;
# aL = 100; #------------------------------------------------------- L
# ay_points <- seq(y_lower, y_lower + Ny, by = dy);
# ay_steps <- length(y_points);
# 
# ah_0 <- 50
# aalpha <- 2
# ac <- 1500
# af <- 100
# aomega <- 2*pi*f
# asigma <- 300 #-------------------------------------------------- sigma
# 
# 
# aH00 <- h_0 + tan(alpha*pi/180)*0;
# 
# aH1 <- h_0 + tan(alpha *pi/180)*y_points
# aH <-t(matrix(H1, nrow=x_steps, ncol=length(H1), byrow=TRUE));


# ----------------------------------------------------------------------------------------------------
# FUNKCIJE
# ----------------------------------------------------------------------------------------------------

make_f <- function(omega, j, H, K00j) {

  umesni <- omega**2/c**2 - (pi*j)**2/(H*H)
  Kj <- sqrt(apply(umesni, 1, as.complex))
  K <- Kj**2 - K00j**2;
  f <- K/2
  return(f)
  
}

make_phi_j <- function(j, z, H) {
  new <- sqrt(2/H)*sin(pi*j*z/H)
  return(new)
}

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



make_P <- function(y_steps, K00j, Q_pml, dy, dx) {
  
  P <- matrix(0, y_steps, y_steps);
  diag(P) <- 2*1i*K00j/dx + (1 + Q_pml)/dy**2;
  
  new_diag_P <- matrix(0, y_steps, y_steps);
  diag(new_diag_P) <- - (1 + Q_pml)/(2*dy**2)
  
  lower_diag_P <- rbind(replicate(y_steps, 0), cbind(new_diag_P[1:y_steps-1, 1:y_steps-1], replicate(y_steps-1, 0))) # od zacetka
  upper_diag_P <- rbind(cbind(replicate(y_steps-1, 0), new_diag_P[2:y_steps, 2:y_steps]), replicate(y_steps, 0))
  P <- P + lower_diag_P + upper_diag_P
  
  return(P)
}

make_and_solve_system <- function(y_points, K00j, x_steps, f, P, Q_pml, dx, dy, lower_diagonal1, upper_diagonal1) {
  
  # initial condition -------------------------------------------
  a0 <- 1/(2 * sqrt(pi))* exp(-y_points**2 * K00j**2);
  
  B <- matrix(0, y_steps, x_steps);
  B[, 1] <- a0
  
  for (x in 2:x_steps){ 
    umesna <- matrix(0, y_steps, y_steps);
    diag(umesna) <- f[x,];
    
    P_update <- P - umesna
    d <- P_update%*%B[,x-1]
    
    nasa_diag <- 2*1i*K00j/dx - (1 + Q_pml)/dy**2 + f[x,]
    
    resitev <- tridiagmatrix(lower_diagonal1, nasa_diag, upper_diagonal1, d)   
    B[, x] <- resitev
    
  }
  return(B)
}

to_decibels <- function(m) {
  s <- 20*log10(abs(m))[,-c(1)]
  return(s)
}
# ----------------------------------------------------------------------------------------------------
# SKUPI
# ----------------------------------------------------------------------------------------------------
A <- matrix(0, 1001, 1001)

for (j in 1:5){
  
  K00j <- sqrt(omega**2/c**2 - (pi*j/H00[1])**2)


  funk <- make_f(omega, j, H, K00j)
  Q_pml <- make_Q_pml(y_points, y_steps, PML, dy, Ny, y_lower, L, sigma)
  P <- make_P(y_steps, K00j, Q_pml, dy, dx)

  n <- (1 + Q_pml)/(2*dy**2)
  lower_diagonal1 <- n[2:length(n)]
  upper_diagonal1 <- n[1:length(n)-1]

  B <- make_and_solve_system(y_points, K00j, x_steps, funk, P, Q_pml, dx, dy, lower_diagonal1, upper_diagonal1)

  phi_j <- make_phi_j(j, z_source, H)

  A <- A + B*matrix(exp(K00j*1i*x_points), nrow=y_steps, ncol=length(x_points), byrow=TRUE)*phi_j
  j = j+1
}






plot_ly(z = to_decibels(A), y = y_points, x = x_points, colorscale = "Jet", 
        zauto = F,
        zmin = max(max(to_decibels(A))) - 30,
        zmax = max(max(to_decibels(A))) - 8,
        type = "heatmap") %>% layout(yaxis = list(autorange = "reversed", title = 'y'), 
                                     title = 'z20_L100_f100_alfa1', 
                                     xaxis = list(title = 'x'))
