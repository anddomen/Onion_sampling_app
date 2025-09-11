library(shiny)
library(bslib)
library(gamlss)

source("00_setup.R")

# Define UI ----
ui <- fluidPage(
    
  ## App title ----
    page_sidebar(
      title = "Onion Sampling App",
      
      sidebar = sidebar(
        width = 450,
        
        
        ## Scenario 1 options ----
        card(
          card_header("Scenario 1"),
          sliderInput("scenario1.prev",
                      "Incoming contamination prevalence proportion",
                      min = 0.001,
                      max = 1,
                      value = 0.017
          ),
          tags$ul(
            tags$li("A low contamination prevalence is 0.001"),
            tags$li("A normal contamination prevalence is 0.017"),
            tags$li("A high contamination prevalence is 0.363")
          ),
          
          numericInput("scenario1.lot",
                       "Enter lot size",
                       value = 22000000,
                       min = 1000,
                       max = 30000000,
                       step = 1000000),
          
          numericInput("scenario1.sample", # Fixed duplicate ID
                       "How many onions do you want to sample?",
                       value = 60,
                       min   = 1,
                       max   = 500),
        ),
        
        ## Scenario 2 options ----
        card(
          card_header("Scenario 2"),
          sliderInput("scenario2.prev",
                      "Incoming contamination prevalence",
                      min = 0.001,
                      max = 1,
                      value = 0.017
          ),
          
          numericInput("scenario2.lot",
                       "Enter lot size",
                       value = 4250000,
                       min = 1000,
                       max = 30000000,
                       step = 1000000),
          
          numericInput("scenario2.sample", # Fixed duplicate ID
                       "How many onions do you want to sample?",
                       value = 60,
                       min   = 1,
                       max   = 500),
          
          actionButton("run_sim", "Go!", class = "btn-success", width = "100%")
        ),
      ),
        

    ## Scenario results ----
    layout_columns(
      col_widths = c(6, 6),
      
      card(
        card_header("Scenario 1 Results"),
        card_body(
          uiOutput("scenario1_output")
        )
      ),
      card(
        card_header("Scenario 2 Results"),
        card_body(
          uiOutput("scenario2_output")
        )
      )
    )
    )
)


# Define server ----
server <- function(input, output) {
  # define number of iterations
  n_sim <- 100
  
  # Create reactive values to store results
  results <- reactiveValues(
    scen1_results = NULL,
    scen2_results = NULL
  )
  
  observeEvent(input$run_sim, {
    # Prep vectors
    scen1.stor <- numeric(length = n_sim)
    scen2.stor <- numeric(length = n_sim)
    
    for (i in 1:n_sim){
      ## Generate distributions based off user input ----
      scenario1.sim <- generate_filtered_ZAGA(input$scenario1.sample,    # Units: logCFU/onion
                                              mu    = norm.incom.contam.mu,
                                              sigma = norm.incom.contam.sigma,
                                              nu    = input$scenario1.prev)
      
      scenario2.sim <- generate_filtered_ZAGA(input$scenario2.sample,    # Units: logCFU/onion
                                              mu    = norm.incom.contam.mu,
                                              sigma = norm.incom.contam.sigma,
                                              nu    = input$scenario2.prev)
      
      ## Update storage vectors ----
      # across iterations, sum up the ones that catch a positive
      scen1.stor[i] <- sum(scenario1.sim > 0)
      scen2.stor[i] <- sum(scenario2.sim > 0)
    }
    
    # Store results in reactive values
    results$scen1_results <- scen1.stor
    results$scen2_results <- scen2.stor
  })
  
  ## Render results for scenario 1 ----
  output$scenario1_output <- renderUI({
    if (is.null(results$scen1_results)) {
      p("Click Go to see results")
    } else {
      results$scen1_results
    }
  })

  ## Render results for scenario 2 ----
  output$scenario2_output <- renderUI({
    if (is.null(results$scen2_results)) {
      p("Click Go to see results")
    } else {
      results$scen2_results
    }
  })

}

# Run the application 
shinyApp(ui = ui, server = server)
