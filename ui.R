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
         ),
      
       br(),
       br(),
       br(),
       br(),
       br(),
       br(),
       br(),
       br(),
       wellPanel(
         h4("Fishing target: Custom stock"),
         numericInput("years.custom", "How many years shall we go fishing for this shark?", value=10,min=1, max=1000, step=1),    
         br(),
         sliderInput("F.in.custom","How hard should we fish the shark?", value=0,min=0, max=0.15, step=0.01),
         br(),
         fluidRow(column(6,numericInput("Max_age_cs", "Maximum age", value=65,min=1, max=1000, step=1)),column(6,numericInput("M_cs", "Natural mortality", value=0.05,min=0.00001, max=5, step=0.01))),
         fluidRow(column(4,numericInput("Linf_cs", "VBGF Linf", value=61,min=1, max=5000, step=0.5)),column(4, numericInput("k_cs", "VBGF k", value=0.13,min=0.00001, max=5, step=0.01)),column(4, numericInput("t0_cs", "VBGF t0", value=0,min=-15, max=15, step=0.01))),
         fluidRow(column(4,numericInput("R0_cs", "Init. Rec.", value=1000,min=1, max=10000000000, step=10)),column(4, numericInput("h_cs", "Steepness", value=0.7,min=0.2, max=1, step=0.01)),column(4, numericInput("a_i_cs", "Age step", value=1,min=0.1, max=100, step=0.5))),
         fluidRow(column(6,numericInput("LW_a_cs", "L-W a", value=1.18e-05,min=0.00000000001, max=1, step=0.0000001)),column(6,numericInput("LW_b_cs", "L-W b", value=3.094,min=0.00001, max=5, step=0.0001))),
         fluidRow(column(6,numericInput("Mat_a_cs", "Maturity slope", value=-0.5,min=-10, max=10, step=0.001)),column(6,numericInput("Mat_L50_cs", "L50%", value=40,min=0.00001, max=5000, step=0.1))),
         fluidRow(column(6,numericInput("Sel_a_cs", "Sel50%", value=10,min=10, max=5000, step=0.1)),column(6,numericInput("Sel_b_cs", "Selectivity slope", value=0.5,min=-10, max=10, step=0.01)))
    )
    ),
    
       mainPanel(
          h4("Rockfish Population"),
          plotOutput("rockfish_pop"),
          h4("Flatfish Population"),
          plotOutput("flatfish_pop"),
          h4("Shark Population"),
          plotOutput("shark_pop"),
          h4("Custom Stock Population"),
          plotOutput("custom_pop")
       )
    ) 
)
)