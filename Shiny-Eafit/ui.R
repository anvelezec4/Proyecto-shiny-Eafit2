library(shiny)
library(dygraphs)
library(knitr)

shinyUI(pageWithSidebar(
   headerPanel("Software Reliability Statistical Models"),
   sidebarPanel(
#       img(src="http://www.maquilempaques.com.co/ESW/Images/LOGO-BANCOLOMBIA.png", width = 350, align="top"),
#       h4("\t"),
      br(" "),
      img(src="http://www.eafit.edu.co/SiteCollectionImages/logo_eafit_55.png", width = 350, align="top"),
      h5(" "),
      
      helpText("Tool for fitting software reliability statistical models given time-between-failures data.", align="center"),
      fileInput("file", label = h4("File input")),
      
      h4("Model Properties"),
      checkboxInput("musaB","Musa Basic"),
      conditionalPanel(condition = "input.musaB == true",numericInput("x", h6("Musa Basic final failure-free time:"), 0 ,min = 0)),
      checkboxInput("jmB","Jelinski-Moranda Imperfect Debugging"),
      conditionalPanel(condition = "input.jmB == true",
                       numericInput("p", h6("Jelinski-Moranda perfect debugging probability:"), 1),
                       numericInput("r", h6("Jelinski-Moranda imperfect debugging probability:"), 0)),
      h1(" "),
      h4("Prediction Parameters"),
      numericInput("pil", h6("Prediction interval length (in hours):"), 24 ,min = 0),
      checkboxInput("musaP","Musa Basic", value = T),
      checkboxInput("duaneP","Duane", value = T),
      checkboxInput("lvP","Littlewood-Verall", value = T),
      checkboxInput("okumotoP","Musa-Okumoto", value = T),
      checkboxInput("geoP","Moranda Geometric", value = T),
      checkboxInput("jmP","Jelinski-Moranda", value = T),
      
      
      sliderInput("prediccion",
                  "Time Serie Prediction",
                  min = 1,  max = 15,  value = 5)
      ),
   
   mainPanel(
      tabsetPanel(
         tabPanel("Ranking",textOutput("a"),dataTableOutput("rank")),
         tabPanel("Prediction",plotOutput("pred_plot", width="100%")),
         tabPanel("Relative Error Plot",plotOutput(outputId = "rel_plot", width = "100%")),
         tabPanel("MVF Plot",plotOutput("total_plot", width="100%")),
         tabPanel("Time Series",plotOutput("Time_series", width="100%")),
         tabPanel("Time Series2",dygraphOutput("Time_series2", width="100%")),
         tabPanel("Time Series Prediction",tableOutput("Prediction_values"),dataTableOutput("rank 100%"))
      )
   )
))