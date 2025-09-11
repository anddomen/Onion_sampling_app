library(shiny)
library(bslib)
library(gamlss)
library(tidyverse)

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
  n_sim <- 10000
  
  # Create reactive values to store results
  results <- reactiveValues(
    pos.lots.scen1_results  = NULL,
    pos.lots.scen2_results  = NULL,
    pos.onions.scen1_results = NULL,
    pos.onions.scen2_results = NULL,
    scen1.sim_results        = NULL,
    scen2.sim_results        = NULL
  )
  
  observeEvent(input$run_sim, {
    # Prep vectors
    # positive lots
    pos.lots.scen1.stor <- numeric(length = n_sim)
    pos.lots.scen2.stor <- numeric(length = n_sim)
    
    # positive onions
    pos.onions.scen1.stor <- numeric(length = n_sim)
    pos.onions.scen2.stor <- numeric(length = n_sim)
    
    # actual simulation values
    scen1.sim.stor <- tibble()
    scen2.sim.stor <- tibble()
    
    for (i in 1:n_sim){
      ## Generate distributions based off user input ----
      scenario1.sim <- generate_filtered_ZAGA(input$scenario1.sample,    # Units: logCFU/onion
                                              mu    = norm.incom.contam.mu,
                                              sigma = norm.incom.contam.sigma,
                                              nu    = 1-input$scenario1.prev)
      
      scenario2.sim <- generate_filtered_ZAGA(input$scenario2.sample,    # Units: logCFU/onion
                                              mu    = norm.incom.contam.mu,
                                              sigma = norm.incom.contam.sigma,
                                              nu    = 1-input$scenario2.prev)
      
      ## Update storage vectors ----
      # across iterations (lots), sum up any time the any of the samples are pos
      pos.lots.scen1.stor[i] <- sum(any(scenario1.sim > 0))
      pos.lots.scen2.stor[i] <- sum(any(scenario2.sim > 0))
      
      # across iterations, sum up the number of positive onions
      pos.onions.scen1.stor[i] <- sum(scenario1.sim > 0)
      pos.onions.scen2.stor[i] <- sum(scenario2.sim > 0)
      
      # save the actual simulations
      # scenario 1
      iteration.data.scen1 <- tibble(
        iteration = i,
        sample    = 1:length(scenario1.sim),
        value     = scenario1.sim
      )
      scen1.sim.stor <- bind_rows(scen1.sim.stor, iteration.data.scen1)
      
      # # scenario 2
      # iteration.data.scen2 <- tibble(
      #   iteration = i,
      #   sample    = 1:length(scenario2.sim),
      #   value     = scenario2.sim
      # )
      # scen2.sim.stor <- bind_rows(scen2.sim.stor, iteration.data.scen2)
    }
    
    # Store results in reactive values
    results$pos.lots.scen1_results   <- pos.lots.scen1.stor
    results$pos.lots.scen2_results   <- pos.lots.scen2.stor
    results$pos.onions.scen1_results <- pos.onions.scen1.stor 
    results$pos.onions.scen2_results <- pos.onions.scen2.stor  
    # results$scen1.sim_results        <- scen1.sim.stor
    # results$scen2.sim_results        <- scen2.sim.stor
  })
  
  ## Render results for scenario 1 ----
  output$scenario1_output <- renderUI({
    if (is.null(results$pos.lots.scen1_results)) {
      p("Click Go to see results")
    } else {
      div(
        h5("Summary Statistics"),
        h6(paste("Across", n_sim, "simulated lots:")),
        p(paste("Number of positive lots caught:", sum(results$pos.lots.scen1_results))),
        p(paste("Number of positive onions caught:", sum(results$pos.onions.scen1_results), "out of", input$scenario1.sample*n_sim, "total onions sampled.")),
        plotOutput("scenario1_plot", height = "300px")
      )
    }
  })
  
  ## Render results for scenario 2 ----
  output$scenario2_output <- renderUI({
    # Check for the correct variable name
    if (is.null(results$pos.lots.scen2_results)) {
      p("Click Go to see results")
    } else {
      div(
        h5("Summary Statistics"),
        h6(paste("Across", n_sim, "simulated lots:")),
        p(paste("Number of positive lots caught:", sum(results$pos.lots.scen2_results))),
        p(paste("Number of positive onions caught:", sum(results$pos.onions.scen2_results))),

        # plotOutput("scenario2_plot", height = "300px")
      )
    }
  })
  
  # plot outputs
  # output$scenario1_plot <- renderPlot({
  #   if (!is.null(results$scen1.sim_results)) {
  #     results$scen1.sim_results |> 
  #       ggplot(aes(x = value)) +
  #       geom_density() +
  #       labs(title = "Distribution of Simulation Values",
  #            x = "Contamination Level (log CFU/onion)",
  #            y = "Number of onions")
  #     # hist(results$pos.onions.scen1_results,
  #     #      main = "Distribution of Positive Onions (Scenario 1)",
  #     #      xlab = "Number of Positive Onions",
  #     #      col = "lightblue")
  #   }
  # })
  # 
  # output$scenario2_plot <- renderPlot({
  #   if (!is.null(results$pos.onions.scen2_results)) {
  #     hist(results$pos.onions.scen2_results, 
  #          main = "Distribution of Positive Onions (Scenario 2)",
  #          xlab = "Number of Positive Onions")
  #   }
  # })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
