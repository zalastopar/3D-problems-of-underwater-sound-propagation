setwd("~/AMAT/IZMENJAVA/nov_poskus")
library(shiny)
library(plotly)
library(dplyr)
library(geometry)
library(shinythemes)
library("grDevices")
library(tidyr)
library(reshape2)

options(shiny.maxRequestSize=30*1024^2)
# -------------------------------------------------------------------------------
# Function for relief plot
# -------------------------------------------------------------------------------
add_true_false <- function(op, vr){
  n <- vr
  n[n == n] <- NA
  for (i in 1:nrow(op)){
    if (op[i,3] == TRUE){
      n[op[i, 1], op[i, 2]] <- vr[op[i, 1], op[i, 2]]}
  }
  return(n)
}

find_square <- function(y, x, xA, yA, xB, yB, xD, yD){
  # function checks if point (x, y) is inside the square, where 
  # (xA, yA), (xB, yB) and (xD, yD) are 3 corner points
  AB <- c(xB - xA, yB - yA)
  AD <- c(xD - xA, yD - yA)
  AM <- c(x - xA, y - yA)
  
  if (
    0 <= dot(AM, AB) & dot(AM, AB) <= dot(AB, AB) &
    0 <= dot(AM, AD) & dot(AM, AD) <= dot(AD, AD)
  ) {return(TRUE)} else {return(FALSE)}
}

get_square_points <- function(file, xA, yA, xB, yB, xD, yD, z_rec){#sq_size, x_rec, y_rec, x_so, y_so, z_rec){
  # Find a rectangle between receiver and source on height z_rec
  # INPUT: file = relief of the bottom of the water body, sq_size = 1/2 length of the rectangle
  # x_, y_, z_ = coordinates of the source and receiver
  # OUTPUT: data frame, points that are in the square keep their value (become z_rec), other NaN
  
  
  # # ----------------------------------------
  # # Get square
  # # ----------------------------------------
  # # https://www.quora.com/How-do-I-find-the-coordinates-of-a-point-at-a-given-distance-from-the-point-on-a-line
  # k_2 <- - (x_rec - x_so)/(y_rec - y_so) 
  # 
  # # first point
  # x_E  <- x_rec - sq_size/sqrt(1 + k_2**2)
  # y_E <- y_rec - k_2 * sq_size/sqrt(1 + k_2**2)
  # # second point
  # x_E2  <- x_rec + sq_size/sqrt(1 + k_2**2)
  # y_E2 <- y_rec + k_2 * sq_size/sqrt(1 + k_2**2)
  
  # ----------------------------------------
  # All points for the plain
  # ----------------------------------------
  
  x_ex <- c(1:nrow(file))
  y_ex <- c(1:ncol(file))
  options <- expand.grid(x_ex, y_ex)
  
  # ----------------------------------------
  # Extract points we need
  # ----------------------------------------
  
  options[,3] <- apply(options, 1, function(x) find_square(x[1],x[2], xA, yA, xB, yB, xD, yD)) # add true/false
  options <- options %>% filter(V3 == TRUE) # take just true
  
  return(options)}

umesna <- function(our_points, file, file2, ravnina, z=0){
  m <- add_true_false(our_points, file)
  #a <- interp2(unique(sort(as.numeric(file2[,2]))), unique(sort(as.numeric(file2[,1]))), as.numeric(file), as.numeric(our_points[,1]), as.numeric(our_points[,2]))
  # a <- interp2(unique(sort(file_melted[,2])), unique(sort(file_melted[,1])), file, square[,2], square[,1])
  # b <- cbind(square, a)
  # c <- spread(b, Var2, a)
  # c <- c[,-c(1,2)]
  
  if (ravnina){
    m[!is.na(m)] <- z # m[!is.na(m)]#z_rec#m[!is.na(m)] + 20#m[!is.na(m)] + 50#z_rec
  }
  
  return(m)
}

plot_3d_basic <- function(file){
  p <- plot_ly(z = as.matrix(file), type = 'surface', colors = 'PuRd')
  return(p)
}

plot_3d <- function(p_b, square_points, x_coo, y_coo, z_coo, points_col, points_text){
  
  df <- data.frame(x_coo, y_coo, z_coo)
  
  # plot
  p <- p_b %>%
    add_trace(data = df, x = x_coo, y = y_coo, z = z_coo, mode = "markers", type = "scatter3d",
              marker = list(size = 4, color = points_col, symbol = 104),
              text = points_text, mode="text", inherit=FALSE, showlegend = F) %>%
    add_trace(z = as.matrix(square_points), showscale = FALSE) 
  
  return(p)}

# -------------------------------------------------------------------------------
# Function for wave plot
# -------------------------------------------------------------------------------
make_f <- function(omega, j, H, K00j) {
  
  umesni <- omega**2/c**2 - (pi*j)**2/(H*H)
  Kj <- t(sqrt(apply(umesni, 1, as.complex)))
  K <- Kj**2 - K00j**2;
  f <- K/2
  return(f)
  
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

manjsa_risba <- function(m){
  m[m >= min(m) + 0.04] <- min(m) + 0.04
  return(m)
}

to_decibels <- function(m) {
  s <- 20*log10(abs(m))[,-c(1)]
  s[s <= max(s) - 150] <- max(s) - 150
  
  return(s)
}

# -------------------------------------------------------------------------------
# SERVER
# -------------------------------------------------------------------------------

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

    # get relief of the bottom
    mydata <- reactive(
      
      {
        bottom <- input$table_H
        ext <- tools::file_ext(bottom$datapath)
        req(bottom)
        validate(need(ext == "csv", "Please upload a csv file")) # we want csv
        df <- read.csv(bottom$datapath)
      }
    )
    
    output$primer <- renderTable(mydata())
    

    #---------------------------------------------------------------------
    
    # source coordinates
    x_coord_source <- reactive({input$x_position_source})
    y_coord_source <- reactive({input$y_position_source})
    z_coord_source <- reactive({input$z_position})

    # receiver coordinates
    x_coord_receiver <- reactive({input$x_position_receiver})
    y_coord_receiver <- reactive({input$y_position_receiver})
    z_coord_receiver <- reactive({input$z_position})
    
    # size of the square
    sq_size <- reactive({input$sq_size})
    
    # square coordinates
    k_2 <- reactive({if (y_coord_receiver() == y_coord_source()){
      0
    } else {
      -(x_coord_receiver() - x_coord_source())/(y_coord_receiver() - y_coord_source())}})
    
    # First - A
    x_A  <- reactive({
      if (y_coord_receiver() == y_coord_source()){x_coord_receiver()} else {
      x_coord_receiver() + sq_size()/sqrt(1 + k_2()**2)}
      
      })
    y_A <- reactive({
      if (y_coord_receiver() == y_coord_source()){y_coord_receiver() + sq_size()} else {
      y_coord_receiver() + k_2() * sq_size()/sqrt(1 + k_2()**2)}
      })
    
    # Second - B
    x_B  <- reactive({
      if (y_coord_receiver() == y_coord_source()){ x_coord_receiver()} else {
      x_coord_receiver() - sq_size()/sqrt(1 + k_2()**2)}
        })
    y_B <- reactive({
      if (y_coord_receiver() == y_coord_source()){y_coord_receiver() - sq_size()} else {
      y_coord_receiver() - k_2() * sq_size()/sqrt(1 + k_2()**2)}
      })
    
    # Third - D 
    x_D <- reactive({
      if (y_coord_receiver() == y_coord_source()){x_coord_source()} else {
      x_coord_source() + sq_size()/sqrt(1 + k_2()**2)}
      })
    y_D <- reactive({
      if (y_coord_receiver() == y_coord_source()){y_coord_source() + sq_size()} else {
      y_coord_source() + k_2() * sq_size()/sqrt(1 + k_2()**2)}
      })


    square_points <- reactive({
        get_square_points(mydata(), x_A(), y_A(), x_B(), y_B(), x_D(), y_D(), z_coord_receiver())
    })

    
    # plot bottom
    #output$plot <- renderPlotly({
    plot <- reactive({
      p <- plot_3d_basic(t(mydata()))
    })
    
    ravnina <- reactive({ 
      file <- as.matrix(mydata())
      file_melted <- melt(file) # transfor matrix to data frame with 3 columns
      square <- get_square_points(file, x_A(), y_A(), x_B(), y_B(), x_D(), y_D(), z_coord_receiver())
      
      umesna(square, file, file_melted, TRUE, z_coord_receiver())})
    
    izrez <- reactive({
      file <- as.matrix(mydata())
      
      file_melted <- melt(file) # transfor matrix to data frame with 3 columns
      square <- get_square_points(file, x_A(), y_A(), x_B(), y_B(), x_D(), y_D(), z_coord_receiver())
      umesna(square, file, file_melted, FALSE)
    })

    output$myplot <- renderPlotly({
      
      x <- c(x_coord_source(), x_coord_receiver(), x_A(), x_B(), x_D())
      y <- c(y_coord_source(), y_coord_receiver(), y_A(), y_B(), y_D())
      z <- c(z_coord_source(), z_coord_receiver(), z_coord_receiver(), z_coord_receiver(), z_coord_receiver())
      
      colors <- c('dodgerblue', 'deeppink', 'red', 'blue', 'green')
      names <- c('receiver', 'source', 'A', 'B', 'D')
      
      df <- data.frame(x, y, z)


      p <- plot_3d(plot(), ravnina(), x, y, z, colors, names)
    })
    
    
    #---------------------------------------------------------------------
    
    plot_waves <- renderImage({
    manjsiH <- mydata()
    
    PML <- TRUE
    
    Nx <- nrow(manjsiH) - 1
    dx <- 1;
    x_points <- seq(0, Nx, by = dx);
    x_steps <- length(x_points);
    
    Ny <- ncol(manjsiH) - 1;
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
    
    j <- 5
    H <- manjsiH
    H00 <-  h_0 + tan(alpha*pi/180);
    K00j <- sqrt(omega**2/c**2 - (pi*j/H00)**2 + 0i);

    f <- make_f(omega, j, H, K00j)
    Q_pml <- make_Q_pml(y_points, y_steps, PML, dy, Ny, y_lower, L, sigma)
    P <- make_P(y_steps, K00j, Q_pml, dy, dx)
    
    n <- (1 + Q_pml)/(2*dy**2)
    lower_diagonal1 <- n[2:length(n)]
    upper_diagonal1 <- n[1:length(n)-1]
    
    B <- make_and_solve_system(y_points, K00j, x_steps, f, P, Q_pml, dx, dy, lower_diagonal1, upper_diagonal1)
    
    
    image(y_points, x_points, to_decibels(abs(B)), col = rainbow(120))
    
    })
})
    
    


