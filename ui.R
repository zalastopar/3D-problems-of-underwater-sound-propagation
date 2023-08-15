setwd("~/AMAT/IZMENJAVA/nov_poskus")
library(shiny)
library(shinythemes)
library(plotly)





# Define UI for application that draws a histogram
# shinyUI(fluidPage(
# 
#   theme = shinytheme("sandstone"),
# 
#   # Application title
#   titlePanel("Underwater sound propagation"),
# 
# 
#   # Sidebar with a slider input for number of bins
#   sidebarLayout(
# 
#     sidebarPanel(
#       splitLayout(
#         fileInput("table_H", "Upload bottom relief", accept = ".csv"),
#         actionButton("button", "Upload", style = "margin-top: +25px;"),
#         cellWidths = c("85%", "15%")),
#       #fileInput("table_H", "Upload bottom relief", accept = ".csv"),
# 
#       hr(style = "border-top: 1px #DFC7B5;margin-top: -10px;"),
# 
#       splitLayout(
#         h3('Source coordinates'),
#         h3('Receiver coordinates')),
# 
#       splitLayout(
#         numericInput('x_position_source', 'x coordinate', value = 0, min = 0, max = 100, width = '80%'),
#         numericInput('x_position_receiver', 'x coordinate', value = 30, min = 0, max = 100, width = '80%')),
# 
#       splitLayout(
#         numericInput('y_position_source', 'y coordinate', value = 0, min = 0, max = 100, width = '80%'),
#         numericInput('y_position_receiver', 'y coordinate', value = 20, min = 0, max = 100, width = '80%')),
# 
#       splitLayout(
#         numericInput('z_position', 'z coordinate (same for both, source and receiver)',
#                      value = 0, min = 0, max = 100, width = '100%')),
# 
#       splitLayout(
#         numericInput('sq_size', 'Size of the square',
#                      min = 0, max = 100, value = 20, width = '100%'))
#     ),
# 
#     # Show a plot of the generated distribution
#     mainPanel(
#       #tableOutput("primer")
#       plotlyOutput('myplot'),
#       imageOutput('plot_waves')
#     )
#   )
# ))


shinyUI(fluidPage(

theme = shinytheme("sandstone"),
titlePanel("Underwater sound propagation"),
  
tabsetPanel(
  tabPanel("Relief", fluid = TRUE,
           sidebarLayout(
             sidebarPanel(
               splitLayout(
                 fileInput("table_H", "Upload bottom relief", accept = ".csv"),
                 actionButton("button", "Upload", style = "margin-top: +25px;"),
                 cellWidths = c("85%", "15%")),
               #fileInput("table_H", "Upload bottom relief", accept = ".csv"),
               
               hr(style = "border-top: 1px #DFC7B5;margin-top: -10px;"),
               
               splitLayout(
                 # h3('Source coordinates'),
                 # h3('Receiver coordinates')),
               
               # splitLayout(
               #   numericInput('x_position_source', 'x coordinate', value = 0, min = 0, max = 100, width = '80%'),
               #   numericInput('x_position_receiver', 'x coordinate', value = 30, min = 0, max = 100, width = '80%')),
               # 
               # splitLayout(
               #   numericInput('y_position_source', 'y coordinate', value = 0, min = 0, max = 100, width = '80%'),
               #   numericInput('y_position_receiver', 'y coordinate', value = 20, min = 0, max = 100, width = '80%')),
               
               splitLayout(
                 numericInput('z_position', 'z coordinate (same for both, source and receiver)',
                              value = 0, min = 0, max = 100, width = '100%'))
               
               # splitLayout(
               #   numericInput('sq_size', 'Size of the square',
               #                min = 0, max = 100, value = 20, width = '100%'))
             )),
             mainPanel(
               plotlyOutput('myplot')
             )
           )
  ),
  tabPanel("Waves", fluid = TRUE,
           sidebarLayout(
             sidebarPanel(

               
               splitLayout(
                 numericInput('L_size', 'Choose L', value = 0, min = 0, max = 100, width = '80%'),
                 numericInput('omega', 'Choose omega', value = 0, min = 0, max = 100, width = '80%')),
               
               splitLayout(
                 numericInput('j_step', 'Choose j', value = 0, min = 0, max = 100, width = '80%'),
                 numericInput('c_size', 'Choose c', value = 0, min = 0, max = 100, width = '80%')),
               
               splitLayout(
                 numericInput('sigma', 'Choose sigma', value = 0, min = 0, max = 100, width = '80%'),
                 numericInput('c_size', 'Choose c', value = 0, min = 0, max = 100, width = '80%'))
               
             ),
             mainPanel(
               imageOutput('plot_waves')
             )

           )
  )
)
))

