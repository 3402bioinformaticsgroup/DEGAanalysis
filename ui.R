# ui.R

library(shiny)

shinyUI(fluidPage(
 titlePanel("Differential Expression Pt2"),
 sidebarLayout(
   sidebarPanel(

      fileInput("NanoStringData", label = h3("Normalized Expression Profile(.csv)")),
      fileInput("PanelReference", label = h3("Panel Reference(.csv)")),
      textInput("Condition1", "Index of condition1:", "1,2"),
      textInput("Condition2", "Index of condition2:", "3,4"),
      textInput("LogFC.cutoff", "LogFC cutoff value:", "1"), #not used
      textInput("outputfile", "Output file name","Test Result"),
      #actionButton("outputreport", "Generate PDF report"),
      downloadButton('downloadData', 'Download Test Result')
      #uiOutput("pdfview")

   ),

 mainPanel(
   tableOutput("input.view"),
   fluidRow(
       column(6, plotOutput("plotPower")),
       column(6, plotOutput("plotDEG"))
)

)

              )
           )
        )




