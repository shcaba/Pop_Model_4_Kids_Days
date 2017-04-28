library(shiny)

shinyUI(
  fluidPage(
    titlePanel("Fish Population Dynamics"),
    h5(p(em("This tool shows you how a population changes over time under fishing."))),
    br(),
    
    headerPanel("Let's go fishing"),
    sidebarLayout(
        sidebarPanel
       (
        wellPanel(
          h4("Fishing target: Rockfish"),
          numericInput("years.rockfish", "How many years shall we go fishing for this rockfish?", value=10,min=1, max=1000, step=1),    
          br(),
          sliderInput("F.in.rockfish","How hard should we fish the rockfish?", value=0,min=0, max=0.1, step=0.01)
         ),
        br(),
        br(),
        br(),
        br(),
        br(),
        br(),
         wellPanel(
           h4("Fishing target: Flatfish"),
           numericInput("years.flatfish", "How many years shall we go fishing for this flatfish?", value=10,min=1, max=1000, step=1),    
           br(),
           sliderInput("F.in.flatfish","How hard should we fish the flatfish?", value=0,min=0, max=0.3, step=0.01)
         ),
         br(),
         br(),
         br(),
         br(),
         br(),
         br(),
         wellPanel(
           h4("Fishing target: shark"),
           numericInput("years.shark", "How many years shall we go fishing for this shark?", value=10,min=1, max=1000, step=1),    
           br(),
           sliderInput("F.in.shark","How hard should we fish the shark?", value=0,min=0, max=0.15, step=0.01)
         )
         ),
         
       mainPanel(
          h4("Rockfish Population"),
          plotOutput("rockfish_pop"),
          h4("Flatfish Population"),
          plotOutput("flatfish_pop"),
          h4("Shark Population"),
          plotOutput("shark_pop")
        )
    ) 
)
)