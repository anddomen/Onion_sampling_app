library(shiny)
library(bslib)
library(gamlss)
library(tidyverse)
library(wesanderson)

source("00_setup.R")


# Define UI ----
ui <- fluidPage(
  theme = bs_theme(version = 5, bootswatch = "lux"),
  
  div(class = "navbar navbar-expand-lg bg-dark",
      div(class = "container-fluid",
          h1("Product Sampling App", 
             class = "navbar-brand mb-0 h1",
             style = "color: white !important;")
      )
  ),
  
  
  ## App title ----
    page_sidebar(
      title = NULL,
      
      sidebar = sidebar(
        width = 450,
        
        
        ## Scenario 1 options ----
        card(
          card_header(class = "bg-dark text-white",
            "Scenario 1"),

          numericInput("scenario1.lot",
                       "Enter lot size",
                       value = 22000000,
                       min   = 1000,
                       max   = 30000000,
                       step  = 1000000),
          
          numericInput("scenario1.sample", 
                       "How many onions do you want to sample?",
                       value = 60,
                       min   = 1,
                       max   = 500),
          
          sliderInput("scenario1.prev",
                      "Percent of lot that's contaminated",
                      min   = 0.01,
                      max   = 10,
                      value = 0.1,
                      step  = 0.01,
                      post  = "%",
                      ticks = FALSE
          ),
          
        ),
        
        ## Scenario 2 options ----
        card(
          card_header(class = "bg-dark text-white",
                      "Scenario 2"),
          
          numericInput("scenario2.lot",
                       "Enter lot size",
                       value = 4250000,
                       min   = 1000,
                       max   = 30000000,
                       step  = 1000000),
          
          numericInput("scenario2.sample", # Fixed duplicate ID
                       "How many onions do you want to sample?",
                       value = 60,
                       min   = 1,
                       max   = 500),
          
          sliderInput("scenario2.prev",
                      "Percent of lot that's contaminated",
                      min   = 0.01,
                      max   = 10,
                      value = 0.1,
                      step  = 0.01,
                      post  = "%",
                      ticks = FALSE
          ),
          
          actionButton("run_sim", "Go!", class = "btn-success", width = "100%")
        ),
      ),
        

    ## Scenario results ----
    layout_columns(
      col_widths = c(6, 6, 12),
      row_heights = c(1,2), 
      
      card(
        card_header(class="card text-white bg-info mb-3",
          "Scenario 1 Results"),
        card_body(
          uiOutput("scenario1_output")
        )
      ),
      card(
        card_header(class = "card text-white bg-warning mb-3",
          "Scenario 2 Results"),
        card_body(
          uiOutput("scenario2_output")
        )
      ),
      card(
        card_body(
          plotOutput("summary_plot")
        )
      )
    )
    ),
  div(
    style = "position: fixed; bottom: 20px; right: 20px; z-index: 1000;",
    div(
      style = "display: flex; align-items: center; gap: 15px;",
      # Logo 1 - Replace with your actual logo path/URL
      tags$img(
        src = "labLogo_orange+black.png",
        # alt = "Logo 1",
        style = "height: 100px; width: auto;"
      ),
      # Logo 2 - Replace with your actual logo path/URL  
      tags$img(
        src = "OSU_horizontal_2C_O_over_B.png", 
        # alt = "Logo 2", 
        style = "height: 50px; width: auto;"
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

    for (i in 1:n_sim){
      ## Generate distributions based off user input ----
      scenario1.sim <- generate_filtered_ZAGA(input$scenario1.sample,    # Units: logCFU/onion
                                              mu    = norm.incom.contam.mu,
                                              sigma = norm.incom.contam.sigma,
                                              nu    = 1-(input$scenario1.prev/100))
      
      scenario2.sim <- generate_filtered_ZAGA(input$scenario2.sample,    # Units: logCFU/onion
                                              mu    = norm.incom.contam.mu,
                                              sigma = norm.incom.contam.sigma,
                                              nu    = 1-(input$scenario2.prev/100))
      
      ## Update storage vectors ----
      # across iterations (lots), sum up any time the any of the samples are pos
      pos.lots.scen1.stor[i] <- sum(any(scenario1.sim > 0))
      pos.lots.scen2.stor[i] <- sum(any(scenario2.sim > 0))
      
      # across iterations, sum up the number of positive onions
      pos.onions.scen1.stor[i] <- sum(scenario1.sim > 0)
      pos.onions.scen2.stor[i] <- sum(scenario2.sim > 0)
      
    }
    
    # Store results in reactive values
    results$pos.lots.scen1_results   <- pos.lots.scen1.stor
    results$pos.lots.scen2_results   <- pos.lots.scen2.stor
    results$pos.onions.scen1_results <- pos.onions.scen1.stor 
    results$pos.onions.scen2_results <- pos.onions.scen2.stor  
  })
  
  ## Create the plot data ----
  create_summary_data <- function(results, n_sim) {
    
    # Calculate proportions for lots
    lots_data <- data.frame(
      scenario = c("Scenario 1", "Scenario 2"),
      positive = c(
        sum(results$pos.lots.scen1_results) / n_sim * 100,
        sum(results$pos.lots.scen2_results) / n_sim * 100
      ),
      metric = "Lots"
    )
    
    # Calculate proportions for onions
    total_onions_scen1 <- length(results$pos.onions.scen1_results) * input$scenario1.sample
    total_onions_scen2 <- length(results$pos.onions.scen2_results) * input$scenario2.sample
    
    onions_data <- data.frame(
      scenario = c("Scenario 1", "Scenario 2"),
      positive = c(
        sum(results$pos.onions.scen1_results) / total_onions_scen1 * 100,
        sum(results$pos.onions.scen2_results) / total_onions_scen2 * 100
      ),
      metric = "Onions"
    )
    
    # Combine and add negative percentages
    combined_data <- rbind(lots_data, onions_data)
    combined_data$negative <- 100 - combined_data$positive
    
    # Reshape to long format for stacking
    plot_data <- combined_data %>%
      pivot_longer(cols = c(positive, negative), 
                   names_to = "result", 
                   values_to = "percentage") %>%
      mutate(
        result = factor(result, levels = c("negative", "positive")),
        metric = factor(metric, levels = c("Lots", "Onions"))
      )
    
    return(plot_data)
  }
  
  ## Render results for scenario 1 ----
  output$scenario1_output <- renderUI({
    if (is.null(results$pos.lots.scen1_results)) {
      p("Click Go to see results")
    } else {
      div(
        h5("Summary Statistics"),
        h6(paste("Across", 
                 format(n_sim, big.mark = ","), 
                 "simulated lots:")),
        p(paste("Number of positive lots caught:",
                format(sum(results$pos.lots.scen1_results), big.mark = ","))),
        p(paste("Number of positive onions caught:",
                format(sum(results$pos.onions.scen1_results), big.mark = ","),
                "out of",
                format(input$scenario1.sample * n_sim, big.mark = ",", scientific = FALSE),
                "total onions sampled."))
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
        h6(paste("Across", 
                 format(n_sim, big.mark = ","), 
                 "simulated lots:")),
        p(paste("Number of positive lots caught:",
                format(sum(results$pos.lots.scen2_results), big.mark = ","))),
        p(paste("Number of positive onions caught:",
                format(sum(results$pos.onions.scen2_results), big.mark = ","),
                "out of",
                format(input$scenario2.sample * n_sim, big.mark = ",", scientific = FALSE),
                "total onions sampled."))
      )
      
    }
  })
  
  ## Generate summary plot ----
  output$summary_plot <- renderPlot({
    if (!is.null(results$pos.lots.scen1_results)) {
      
      plot_data <- create_summary_data(results, n_sim)
      
      ggplot(plot_data, aes(x = scenario, 
                            y = percentage, 
                            fill = result)) +
        geom_col(position = "stack", width = 0.7) +
        facet_grid(rows = vars(metric), 
                   labeller = labeller(metric = c("Lots" = "Positive Lots Detected", 
                                                  "Onions" = "Positive Onions Detected"))) +
        scale_fill_manual(values = c("negative" = "#9fc2b2", "positive" = "#a90636"),
                          labels = c("negative" = "Negative", "positive" = "Positive")) +
        scale_y_continuous(labels = function(x) paste0(x, "%"), 
                           limits = c(0, 100)) +
        labs(
          title = "Detection Results by Scenario",
          subtitle = paste("Based on", format(n_sim, big.mark = ","), "simulated lots"),
          x = "Scenario",
          y = "Percentage",
          fill = "Result"
        ) +
        theme_linedraw() +
        theme(
          text = element_text(size = 18),
          strip.text = element_text(face = "bold"),
          strip.background = element_rect(fill = "#343a40"),
          axis.title.x = element_blank(),
          legend.position = "bottom"
        )
      
    }
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
