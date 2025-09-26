library(shiny)
library(shinycssloaders)
library(bslib)
library(gamlss)
library(tidyverse)

source("00_setup.R")


# Define UI ----
ui <- fluidPage(
  theme = bs_theme(version = 5, bootswatch = "lux"),
  
  div(class = "navbar navbar-expand-lg bg-dark",
      div(class = "container-fluid",
          h1("Product Sampling", 
             class = "navbar-brand mb-0 h1",
             style = "color: white !important;"),
          
          # Add logos to the right side
          div(class = "navbar-nav ms-auto",
              div(class = "d-flex align-items-center gap-3",
                  tags$img(
                    src = "labLogo_orange+white.png",
                    style = "height: 100px; width: auto;"
                  ),
                  tags$img(
                    src = "OSU_horizontal_2C_O_over_W.png", 
                    style = "height: 50px; width: auto;"
                  )
              )
          )
      )
  ),
  
  
  ## App title ----
  page_sidebar(
    title = NULL,
    
    sidebar = sidebar(
      width = 450,
      
      
      ## Scenario 1 options ----
      card(
        card_header(class="card text-white bg-info mb-3",
                    "Scenario 1"),
        
        numericInput("scenario1.lot",
                     "Enter lot size",
                     value = 22000000,
                     min   = 1000,
                     max   = 30000000,
                     step  = 1000000),
        
        numericInput("scenario1.sample", 
                     "How many samples do you want to take?",
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
        card_header(class = "card text-white bg-warning mb-3",
                    "Scenario 2"),
        
        numericInput("scenario2.lot",
                     "Enter lot size",
                     value = 4250000,
                     min   = 1000,
                     max   = 30000000,
                     step  = 1000000),
        
        numericInput("scenario2.sample", # Fixed duplicate ID
                     "How many samples do you want to take?",
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
      col_widths = c(6, 6),
      
      ### Scenario 1 ----
      card(
        card_header(class="card text-white bg-info mb-3",
                    "Scenario 1 Results"),
        card_body(
          withSpinner(uiOutput("scenario1_output"), type = 4, color = "steelblue", proxy.height = "200px"),
          withSpinner(plotOutput("scen1_posLot_plot"), type = 0, proxy.height = "0px"),
          withSpinner(uiOutput("scenario1_graph_title"), type = 0, proxy.height = "0px"),
          withSpinner(plotOutput("scen1_posOnions.posLot_plot"), type = 0, proxy.height = "0px"),
        )
      ),
      
      ### Scenario 2 ----
      card(
        card_header(class = "card text-white bg-warning mb-3",
                    "Scenario 2 Results"),
        card_body(
          withSpinner(uiOutput("scenario2_output"), type = 4, color = "steelblue", proxy.height = "200px"),
          withSpinner(plotOutput("scen2_posLot_plot"), type = 0, proxy.height = "0px"),
          withSpinner(uiOutput("scenario2_graph_title"), type = 0, proxy.height = "0px"),
          withSpinner(plotOutput("scen2_posOnions.posLot_plot"), type = 0, proxy.height = "0px"),
        )
      )
    )
  )
)


# Define server ----
server <- function(input, output) {
  # define number of iterations
  n_sim <- 100000
  

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
    showSpinner("scenario1_output")
    showSpinner("scen1_posLot_plot")
    showSpinner("scenario1_graph_title")
    showSpinner("scen1_posOnions.posLot_plot")
    showSpinner("scenario2_output")
    showSpinner("scen2_posLot_plot")
    showSpinner("scenario2_graph_title")
    showSpinner("scen2_posOnions.posLot_plot")
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
    
    hideSpinner("scenario1_output")
    hideSpinner("scen1_posLot_plot")
    hideSpinner("scenario1_graph_title")
    hideSpinner("scen1_posOnions.posLot_plot")
    hideSpinner("scenario2_output")
    hideSpinner("scen2_posLot_plot")
    hideSpinner("scenario2_graph_title")
    hideSpinner("scen2_posOnions.posLot_plot")

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
  
  ## Scenario 1 result render ----
  ### Text output ----
  output$scenario1_output <- renderUI({
    if (is.null(results$pos.lots.scen1_results)) {
      p("Click Go to see results")
    } else {
      div(
        h5("Summary Statistics"),
        h6(paste("Across", 
                 format(n_sim, big.mark = ","), 
                 "simulated lots:")),
        p("Number of positive lots caught:",
          tags$span(class = "h2",
                    format(sum(results$pos.lots.scen1_results), big.mark = ","))),
        p("That's ", 
          tags$span(class = "h2",
                    paste0(sum(results$pos.lots.scen1_results)/n_sim * 100, "%")))
        
      )
    }
  })
  
  output$scenario1_graph_title <- renderUI({
    if (is.null(results$pos.lots.scen1_results)) {
      NULL
    } else {
      h5("Even in positive lots, few individual samples test positive.")
    }
  })
  
  ### Positive lot graph ----
  output$scen1_posLot_plot <- renderPlot({
    showSpinner()
    if (!is.null(results$pos.lots.scen1_results)) {
      plot_data <- create_summary_data(results, n_sim) 
      
      plot_data |> 
        filter(metric == "Lots",
               scenario == "Scenario 1") |> 
        mutate(pct_formatted = paste0(round(percentage, digits = 1), "%")) |> 
        ggplot(aes(x = scenario, 
                   y = percentage, 
                   fill = result,
                   label = pct_formatted)) +
        geom_col(position = "stack", width = 0.7) +
        geom_text(
          position = position_stack(vjust = 0.5),
          size = 10,
          fontface = "bold",
          color = "gray0"
        ) +
        scale_fill_manual(values = c("negative" = "#9fc2b2", "positive" = "#a90636"),
                          labels = c("negative" = "Negative", "positive" = "Positive")) +
        scale_y_continuous(labels = function(x) paste0(x, "%"), 
                           limits = c(0, 100)) +
        labs(
          title = "Scenario 1 — Detection Results by Lot",
          subtitle = paste("Based on", format(n_sim, big.mark = ","), "simulated lots"),
          y = "Percentage",
          fill = "Result"
        ) +
        theme_classic()+
        coord_flip()+
        theme(
          text = element_text(size = 18),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "bottom"
        ) 
      
    }
  })
  
  ### Positive samples graph ----
  output$scen1_posOnions.posLot_plot <- renderPlot({
    if (!is.null(results$pos.lots.scen1_results)) {
      
      onions_data <- data.frame(
        scenario = c("Scenario 1"),
        positive = c(
          results$pos.onions.scen1_results)) |>
        filter(positive != 0)
      
      onions_data |>
        ggplot(aes(x = positive)) +
        geom_bar(fill = "#1f9bcf")+
        geom_text(
          stat = "count",
          position = position_stack(vjust = 0.5),
          size = 10,
          fontface = "bold",
          color = "gray0",
          aes(label = ..count..))+
        scale_x_continuous(breaks = seq(min(onions_data$positive), max(onions_data$positive), by = 1)) +
        labs(title = paste("Distribution across the", format(nrow(onions_data), big.mark = ",", scientific = FALSE), "simulated lots that tested positive"),
             y = "Number of lots that test positive",
             x = "Number of positive samples per positive lot") +
        theme_classic()+
        theme(text = element_text(size = 18))
    }
  })
  
  
  
  
  
  
  
  
  ## Scenario 2 result render ----
  ### Text output ----
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
        p("Number of positive lots caught:",
          tags$span(class = "h2",
                    format(sum(results$pos.lots.scen2_results), big.mark = ","))),
        p("That's ", 
          tags$span(class = "h2",
                    paste0(sum(results$pos.lots.scen2_results)/n_sim * 100, "%")))
      )
      
    }
  })
  
  output$scenario2_graph_title <- renderUI({
    if (is.null(results$pos.lots.scen2_results)) {
      NULL
    } else {
      h5("Even in positive lots, few individual samples test positive.")
    }
  })
  
  ### Positive lot graph ----
  output$scen2_posLot_plot <- renderPlot({
    if (!is.null(results$pos.lots.scen2_results)) {
      plot_data <- create_summary_data(results, n_sim) 
      
      plot_data |> 
        filter(metric == "Lots",
               scenario == "Scenario 2") |> 
        mutate(pct_formatted = paste0(round(percentage, digits = 1), "%")) |> 
        ggplot(aes(x = scenario, 
                   y = percentage, 
                   fill = result,
                   label = pct_formatted)) +
        geom_col(position = "stack", width = 0.7) +
        geom_text(
          position = position_stack(vjust = 0.5),
          size = 10,
          fontface = "bold",
          color = "gray0"
        ) +
        scale_fill_manual(values = c("negative" = "#9fc2b2", "positive" = "#a90636"),
                          labels = c("negative" = "Negative", "positive" = "Positive")) +
        scale_y_continuous(labels = function(x) paste0(x, "%"), 
                           limits = c(0, 100)) +
        labs(
          title = "Scenario 2 — Detection Results by Lot",
          subtitle = paste("Based on", format(n_sim, big.mark = ","), "simulated lots"),
          y = "Percentage",
          fill = "Result"
        ) +
        theme_classic()+
        coord_flip()+
        theme(
          text = element_text(size = 18),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "bottom"
        ) 
      
    }
  })
  
  ### Positive samples graph ----
  output$scen2_posOnions.posLot_plot <- renderPlot({
    if (!is.null(results$pos.lots.scen2_results)) {
      
      onions_data <- data.frame(
        scenario = c("Scenario 2"),
        positive = c(
          results$pos.onions.scen2_results)) |>
        filter(positive != 0)
      
      onions_data |>
        ggplot(aes(x = positive)) +
        geom_bar(fill = "#f0ad4e")+
        geom_text(
          stat = "count",
          position = position_stack(vjust = 0.5),
          size = 10,
          fontface = "bold",
          color = "gray0",
          aes(label = ..count..))+
        scale_x_continuous(breaks = seq(min(onions_data$positive), max(onions_data$positive), by = 1)) +
        labs(title = paste("Distribution across the", format(nrow(onions_data), big.mark = ",", scientific = FALSE), "simulated lots that tested positive"),
             y = "Number of lots that test positive",
             x = "Number of positive samples per positive lot") +
        theme_classic()+
        theme(text = element_text(size = 18))
    }
  })
  
  
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
