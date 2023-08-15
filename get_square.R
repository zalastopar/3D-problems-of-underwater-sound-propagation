
library(plotly)
library(dplyr)

library(pracma) # za dot



file <- t(volcano)

# za pike --------------------------

x_r <- 10
y_r <- 50

x_s <- 50
y_s <- 30

z_r <- 180

x <- c(x_r, x_s)
y <- c(y_r, y_s)
z <- c(z_r, z_r)

df <- data.frame(x, y, z)





find_square <- function(y, x, x_s, y_s, x_r, y_r, E){
  k = (y_s - y_r)/(x_s - x_r)


  if (y_s <= y_r){ # receiver in the right
  if (y >= -1/k*x + y_s + 1/k*x_s & # desno od B
      (y <= -1/k*x + y_r + 1/k*x_r) & # levo od D
      (x <= 1/k * (y - (y_s + 20) + k*(x_s + 20))) & # pod C
      (x >= 1/k * (y - (y_s - 20) + k*(x_s - 20))) # nad A

  ){return(TRUE)} else {
    return(FALSE)
  }
  } else{
    if (y <= -1/k * x + y_s + 1/k*x_s & # desno od B
        y >= -1/k * x + y_r + 1/k*x_r & # levo od D
        (x <= 1/k * (y - (y_s + 20) + k*(x_s + 20))) & # nad C
        x >= 1/k * (y - (y_s - 20) + k*(x_s - 20)) # pod A
    ){
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
}

dodaj_vrednosti <- function(op, vr){
  n <- vr
  n[n == n] <- NA
  for (i in 1:nrow(op)){
    if (op[i,3] == TRUE){
      n[op[i, 1], op[i, 2]] <- vr[op[i, 1], op[i, 2]]}
  }
  return(n)
}


# https://math.stackexchange.com/questions/2157931/how-to-check-if-a-point-is-inside-a-square-2d-plane
nova_f <- function(x, y, x_s, y_s, x_r, y_r, E){
  AB <- c(-2*E, 2*E)
  AD <- c(x_r - x_s, y_r - y_s)
  AM <- c(x - x_s - E, y - y_s + E)
  
  if ( (0 <= dot(AM, AB)) & (dot(AM, AB) <= dot(AB, AB)) &
       (0 <= dot(AM, AD)) & (dot(AM, AD) <= dot(AD, AD))) {
    return(TRUE)
  } else {return(FALSE)}
  
}

nova_2 <- function(x, y, x_s, y_s, x_r, y_r, E){
  h <- sqrt((x_r-x_s)**2 + (y_r-y_s)**2)
  w <- 2*E
  
  a <- x - x_r + E
  b <- y - y_r - E
  
  if (0 <= a*2*E/w + b*(-2)*E/h & a*2*E/w + b*(-2)*E/h <= w &
      0 <= a*(x_s-x_r)/w + b*(y_s - y_r)/h & a*(x_s-x_r)/w + b*(y_s - y_r)/h <= h){
    return(TRUE)
  } else {return(FALSE)}
}



options[,3] <- apply(options, 1, function(x) nova_f(x[1],x[2], x_s, y_s, x_r, y_r, E)) # dodamo true/false
options <- options %>% filter(V3 == TRUE) # filtriramo v une s true
# to string options number
#options[,1] <- as.character(options[,1])
#options[,2] <- as.character(options[,2])
m <- dodaj_vrednosti(options, file)

m[!is.na(m)] <- z_r



plot_ly(z = file, type = 'surface') %>%
  add_trace(data = df, x = x, y = y, z = z, mode = "markers", type = "scatter3d",
            marker = list(size = 4, color = c('dodgerblue', 'deeppink'), symbol = 104),
            text=c('receiver', 'source'), mode="text", inherit=FALSE, showlegend = F) %>%
  add_trace(z = m, showlegend = F)




