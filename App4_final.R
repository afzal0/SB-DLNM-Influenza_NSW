##############################
# app.R
# Aesthetic 2-tab Shiny Dashboard for SB-DLNM
# (In the name of Al)
##############################

# -- 1) Required Libraries
library(shiny)
library(bslib)
library(ggplot2)
library(ggrepel)
library(plotly)
library(leaflet)
library(sf)
library(splines)
library(dplyr)
library(coda)
library(dlnm)
library(DT)
library(spdep)
library(viridis)

##############################
# 2) Load Model, Data, and Shapefile
##############################

# Load Model 3 posterior samples
if (file.exists("data/model3_results.RData")) {
  load("data/model3_results.RData")  # loads model3_samples
  model3_results <- model3_samples
  samples_mat <- as.matrix(model3_results)
} else {
  stop("Cannot find 'data/model3_results.RData'")
}

# Load time-series data (data_ts)
if (file.exists("data/data_timeseries.RData")) {
  load("data/data_timeseries.RData")  # loads data_ts
} else {
  stop("Cannot find 'data/data_timeseries.RData'")
}

# If the expected climate variables are missing, create dummy data for them
climate_vars <- c("temperature_2m", "relative_humidity_2m", "rain")
for (var in climate_vars) {
  if (!var %in% names(data_ts)) {
    message("Column '", var, "' missing in data_ts. Creating dummy data for demonstration.")
    if (var == "temperature_2m") {
      data_ts[[var]] <- runif(nrow(data_ts), 15, 35)
    } else if (var == "relative_humidity_2m") {
      data_ts[[var]] <- runif(nrow(data_ts), 30, 100)
    } else if (var == "rain") {
      data_ts[[var]] <- runif(nrow(data_ts), 0, 50)
    }
  }
}

# Ensure data_ts has region_index
if (!"region_index" %in% names(data_ts)) {
  data_ts$region_index <- as.numeric(factor(data_ts$LHD))
}

# Create region names and mapping
region_names <- sort(unique(data_ts$LHD))
region_mapping_ts <- data_ts %>% select(LHD, region_index) %>% distinct()

# Load the shapefile
if (!file.exists("data/NSW_LHD_Boundaries.shp")) {
  stop("Cannot find 'data/NSW_LHD_Boundaries.shp'")
}
shapefile_nsw <- st_read("data/NSW_LHD_Boundaries.shp")
shapefile_nsw <- st_make_valid(shapefile_nsw)

# Suppose the shapefile has a column 'lhd_name' matching data_ts$LHD. Rename to 'LHD'.
if ("lhd_name" %in% names(shapefile_nsw)) {
  shapefile_nsw <- shapefile_nsw %>% rename(LHD = lhd_name)
} else {
  stop("We expected a column 'lhd_name' in the shapefile to match data_ts$LHD. Adjust the rename() step if needed.")
}

##############################
# 3) Define Helper Functions
##############################

safe_ns <- function(x, df = NULL, knots = NULL, Boundary.knots) {
  unique_x <- sort(unique(x))
  if (length(unique_x) < 2) {
    return(matrix(1, nrow = length(x), ncol = 1))
  }
  if (!is.null(knots)) {
    knots <- knots[knots > Boundary.knots[1] & knots < Boundary.knots[2]]
  }
  tryCatch(
    ns(x, df = df, knots = knots, Boundary.knots = Boundary.knots),
    error = function(e) {
      message("ns() error: ", e$message, " => returning constant basis.")
      matrix(1, nrow = length(x), ncol = 1)
    }
  )
}

get_median_coeff <- function(climate_var, region, samples_matrix) {
  prefix <- switch(
    climate_var,
    "temperature_2m" = "b_temp",
    "relative_humidity_2m" = "b_hum",
    "rain" = "b_rain"
  )
  pattern <- paste0("^", prefix, "\\[", region, ",")
  param_names <- grep(pattern, colnames(samples_matrix), value = TRUE)
  median_coeff <- sapply(param_names, function(name) median(samples_matrix[, name]))
  unname(as.vector(median_coeff))
}

predict_prob_occurrence <- function(x_val, climate_var, region, data, samples_matrix) {
  median_coeff <- get_median_coeff(climate_var, region, samples_matrix)
  L <- length(median_coeff)
  b_knots <- range(data[[climate_var]], na.rm = TRUE)
  if (length(x_val) == 1) {
    x_eval <- c(x_val, mean(b_knots))
    basis_mat <- safe_ns(x_eval, df = L, Boundary.knots = b_knots)
    basis_vec <- as.vector(basis_mat[1, ])
  } else {
    basis_vec <- as.vector(safe_ns(x_val, df = L, Boundary.knots = b_knots))
  }
  log_effect <- sum(basis_vec * median_coeff)
  lambda <- exp(log_effect)
  1 - exp(-lambda)
}

predict_rr_curve <- function(climate_var, region, data, samples_matrix) {
  median_coeff <- get_median_coeff(climate_var, region, samples_matrix)
  L <- length(median_coeff)
  x_seq <- seq(
    min(data[[climate_var]], na.rm = TRUE),
    max(data[[climate_var]], na.rm = TRUE),
    length.out = 100
  )
  b_knots <- range(data[[climate_var]], na.rm = TRUE)
  basis_mat <- safe_ns(x_seq, df = L, Boundary.knots = b_knots)
  rr <- exp(as.vector(basis_mat %*% median_coeff))
  data.frame(Exposure = x_seq, RR = rr)
}

compute_region_rr <- function(climate_var, exposure, data, samples_matrix) {
  regs <- sort(unique(data$region_index))
  rr_vals <- sapply(regs, function(reg) {
    median_coeff <- get_median_coeff(climate_var, reg, samples_matrix)
    L <- length(median_coeff)
    b_knots <- range(data[[climate_var]], na.rm = TRUE)
    if (length(exposure) == 1) {
      x_eval <- c(exposure, mean(b_knots))
      basis_mat <- safe_ns(x_eval, df = L, Boundary.knots = b_knots)
      basis_vec <- as.vector(basis_mat[1, ])
    } else {
      basis_vec <- as.vector(safe_ns(exposure, df = L, Boundary.knots = b_knots))
    }
    exp(sum(basis_vec * median_coeff))
  })
  data.frame(region_index = regs, RR = rr_vals)
}

create_percentile_plot <- function(data, model, climate_var, lhd_name) {
  lhd_data <- data[data$LHD == lhd_name, ]
  if (nrow(lhd_data) == 0) return(NULL)
  region_idx <- unique(lhd_data$region_index)[1]
  
  var_knots <- quantile(lhd_data[[climate_var]], c(0.10, 0.50, 0.90), na.rm = TRUE)
  cb <- crossbasis(
    lhd_data[[climate_var]],
    lag = 1,
    argvar = list(fun = "ns", knots = var_knots),
    arglag = list(fun = "ns", df = 2)
  )
  prefix <- switch(climate_var,
                   "temperature_2m" = "b_temp",
                   "relative_humidity_2m" = "b_hum",
                   "rain" = "b_rain"
  )
  pat <- paste0("^", prefix, "\\[", region_idx, ",")
  param_names <- grep(pat, colnames(model), value = TRUE)
  if (length(param_names) == 0) return(NULL)
  median_chain <- sapply(param_names, function(x) median(model[, x]))
  x_range <- range(lhd_data[[climate_var]], na.rm = TRUE)
  x_seq <- seq(x_range[1], x_range[2], length.out = 100)
  ref_val <- switch(
    climate_var,
    "temperature_2m" = 20,
    "relative_humidity_2m" = 60,
    "rain" = 300
  )
  pred <- crosspred(
    cb,
    coef = median_chain,
    vcov = diag(0, ncol(cb), ncol(cb)),
    model.link = "log",
    at = x_seq,
    cen = ref_val
  )
  data.frame(
    Exposure = x_seq,
    fit = pred$allRRfit,
    low = pred$allRRlow,
    high = pred$allRRhigh
  )
}

##############################
# 4. Define the UI
##############################

ui <- navbarPage(
  title = "NSW Health: SB-DLNM Dashboard",
  theme = bslib::bs_theme(
    bootswatch = "flatly",
    primary = "#2c3e50",
    secondary = "#95a5a6"
  ),
  
  # ============== TAB 1: Main Dashboard ==============
  tabPanel(
    "Main Dashboard",
    tags$head(
      tags$style(HTML("
        .well { background-color: #f8f9fa; border-radius: 10px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
        .control-label { color: #2c3e50; font-weight: 600; }
        .plot-container { background: white; padding: 15px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.05); }
        .nav-tabs { margin-bottom: 20px; }
        .title-box { background-color: #ecf0f1; padding: 10px; border-radius: 5px; margin-bottom: 15px; }
        .dashboard-title { color: #2c3e50; font-size: 24px; font-weight: 600; }
      "))
    ),
    
    fluidPage(
      fluidRow(
        column(3,
               wellPanel(
                 style = "position: fixed; width: 23%;",
                 h4("Dashboard Controls", style = "color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px;"),
                 selectInput("dash_region", "Select Region (LHD):", choices = region_names),
                 selectInput("dash_climate", "Select Climate Variable:", choices = climate_vars),
                 sliderInput("dash_exposure", "Exposure Value:",
                             min = 10, max = 50, value = 20,
                             ticks = FALSE
                 ),
                 radioButtons("dash_rr_type", "RR Plot Type:",
                              choices = c("Simple", "Percentile"),
                              selected = "Simple"
                 ),
                 hr(),
                 h5("3D Plot Settings", style = "color: #2c3e50;"),
                 numericInput("dash_3d_npoints", "Surface Resolution:",
                              value = 20, min = 5, max = 50
                 )
               )
        ),
        
        column(9, offset = 3,
               fluidRow(
                 column(6,
                        div(class = "plot-container",
                            h4("Relative Risk Curve", style = "color: #2c3e50;"),
                            plotlyOutput("dash_rrPlot", height = "300px")
                        )
                 ),
                 column(6,
                        div(class = "plot-container",
                            h4("Probability of Occurrence", style = "color: #2c3e50;"),
                            plotlyOutput("dash_probPlot", height = "300px")
                        )
                 )
               ),
               
               fluidRow(
                 column(12,
                        div(class = "plot-container",
                            style = "margin-top: 20px;",
                            h4("Regional Risk Map", style = "color: #2c3e50;"),
                            leafletOutput("dash_map", height = "400px")
                        )
                 )
               ),
               
               fluidRow(
                 column(12,
                        div(class = "plot-container",
                            style = "margin-top: 20px;",
                            h4("3D Risk Surface", style = "color: #2c3e50;"),
                            plotlyOutput("dash_3dPlot", height = "400px")
                        )
                 )
               )
        )
      )
    )
  ),
  
  # ============== TAB 2: Additional Analysis ==============
  tabPanel(
    "Additional Analysis",
    fluidPage(
      fluidRow(
        column(3,
               wellPanel(
                 h4("Analysis Settings", style = "color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px;"),
                 selectInput("analysis_climate", "Climate Variable:", 
                             choices = climate_vars, 
                             selected = "rain"
                 )
               )
        ),
        column(9,
               fluidRow(
                 column(12,
                        div(class = "plot-container",
                            h4("Regional Coefficients", style = "color: #2c3e50;"),
                            DTOutput("coefTable")
                        )
                 )
               ),
               fluidRow(
                 column(12,
                        div(class = "plot-container",
                            style = "margin-top: 20px;",
                            h4("Exposure-Coefficient Relationship", style = "color: #2c3e50;"),
                            plotlyOutput("coefPlot", height = "400px")
                        )
                 )
               )
        )
      )
    )
  )
)

##############################
# 5. Server
##############################

server <- function(input, output, session) {
  
  ### 1) Enhanced RR Plot
  output$dash_rrPlot <- renderPlotly({
    req(input$dash_region, input$dash_climate)
    region_num <- unique(data_ts$region_index[data_ts$LHD == input$dash_region])[1]
    
    if (input$dash_rr_type == "Simple") {
      rr_df <- predict_rr_curve(input$dash_climate, region_num, data_ts, samples_mat)
      p <- ggplot(rr_df, aes(x = Exposure, y = RR)) +
        geom_line(color = "#3498db", size = 1.2) +
        geom_ribbon(aes(ymin = 0.95 * RR, ymax = 1.05 * RR),
                    fill = "#3498db", alpha = 0.2) +
        labs(
          title = paste("Relative Risk:", input$dash_region),
          subtitle = paste("Climate Variable:", input$dash_climate),
          x = input$dash_climate,
          y = "Relative Risk"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 12),
          axis.title = element_text(size = 11),
          panel.grid.minor = element_blank()
        )
      
      ggplotly(p) %>%
        layout(
          hoverlabel = list(bgcolor = "white"),
          plot_bgcolor = "rgba(0,0,0,0)"
        )
    } else {
      pc_df <- create_percentile_plot(data_ts, samples_mat, input$dash_climate, input$dash_region)
      if (is.null(pc_df)) {
        return(ggplotly(ggplot() + ggtitle("No data available")))
      }
      
      p <- ggplot(pc_df, aes(x = Exposure, y = fit)) +
        geom_ribbon(aes(ymin = low, ymax = high),
                    fill = "#3498db", alpha = 0.2) +
        geom_line(color = "#3498db", size = 1.2) +
        scale_y_continuous(trans = "log") +
        labs(
          title = paste("Relative Risk (Percentile):", input$dash_region),
          subtitle = paste("Climate Variable:", input$dash_climate),
          x = input$dash_climate,
          y = "Relative Risk (log scale)"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 12),
          axis.title = element_text(size = 11),
          panel.grid.minor = element_blank()
        )
      
      ggplotly(p) %>%
        layout(
          hoverlabel = list(bgcolor = "white"),
          plot_bgcolor = "rgba(0,0,0,0)"
        )
    }
  })
  
  ### 2) Enhanced Probability Plot
  output$dash_probPlot <- renderPlotly({
    req(input$dash_region, input$dash_climate, input$dash_exposure)
    region_num <- unique(data_ts$region_index[data_ts$LHD == input$dash_region])[1]
    
    x_seq <- seq(10, 50, length.out = 51)
    prob_vals <- sapply(x_seq, function(v) {
      predict_prob_occurrence(v, input$dash_climate, region_num, data_ts, samples_mat)
    })
    df_plot <- data.frame(Exposure = x_seq, Probability = prob_vals)
    
    user_prob <- predict_prob_occurrence(input$dash_exposure, 
                                         input$dash_climate, 
                                         region_num, 
                                         data_ts, 
                                         samples_mat)
    
    p <- ggplot(df_plot, aes(x = Exposure, y = Probability)) +
      geom_ribbon(aes(ymin = 0, ymax = Probability),
                  fill = "#e74c3c", alpha = 0.1) +
      geom_line(color = "#e74c3c", size = 1.2) +
      geom_vline(xintercept = input$dash_exposure,
                 color = "#2980b9",
                 linetype = "dashed",
                 size = 1) +
      geom_point(
        data = data.frame(x = input$dash_exposure, y = user_prob),
        aes(x = x, y = y),
        color = "#2980b9",
        size = 4
      ) +
      annotate("text",
               x = input$dash_exposure,
               y = user_prob,
               label = sprintf("p = %.3f", user_prob),
               color = "#2980b9",
               vjust = -1
      ) +
      labs(
        title = paste("Occurrence Probability:", input$dash_region),
        subtitle = paste("Climate Variable:", input$dash_climate),
        x = input$dash_climate,
        y = "Probability"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 11),
        panel.grid.minor = element_blank()
      )
    
    ggplotly(p) %>%
      layout(
        hoverlabel = list(bgcolor = "white"),
        plot_bgcolor = "rgba(0,0,0,0)"
      )
  })
  
  ### 3) Enhanced Map
  output$dash_map <- renderLeaflet({
    req(input$dash_climate, input$dash_exposure)
    
    rr_df <- compute_region_rr(input$dash_climate, input$dash_exposure, data_ts, samples_mat)
    rr_df <- rr_df %>%
      left_join(region_mapping_ts %>% distinct(), by = "region_index")
    
    map_data <- shapefile_nsw %>%
      left_join(rr_df, by = "LHD")
    
    if (!"RR" %in% names(map_data)) {
      map_data$RR <- 1
    }
    map_data$RR[is.na(map_data$RR)] <- 1
    
    bb <- st_bbox(map_data)
    center_lng <- mean(c(bb["xmin"], bb["xmax"]))
    center_lat <- mean(c(bb["ymin"], bb["ymax"]))
    
    pal <- colorNumeric(
      palette = "RdYlBu",
      domain = map_data$RR,
      reverse = TRUE,
      na.color = "transparent"
    )
    
    leaflet(map_data) %>%
      addProviderTiles(
        providers$CartoDB.Positron,
        options = providerTileOptions(opacity = 0.8)
      ) %>%
      addPolygons(
        fillColor = ~pal(RR),
        weight = 1,
        color = "white",
        dashArray = "3",
        fillOpacity = 0.7,
        highlightOptions = highlightOptions(
          weight = 2,
          color = "#666",
          dashArray = "",
          fillOpacity = 0.9,
          bringToFront = TRUE
        ),
        label = ~sprintf("<strong>%s</strong><br/>RR: %.2f", LHD, RR) %>% lapply(htmltools::HTML)
      ) %>%
      addLegend(
        pal = pal,
        values = ~RR,
        opacity = 0.7,
        title = "Relative Risk",
        position = "bottomright",
        labFormat = labelFormat(digits = 2)
      ) %>%
      setView(lng = center_lng, lat = center_lat, zoom = 6)
  })
  
  ### 4) Enhanced 3D Plot
  output$dash_3dPlot <- renderPlotly({
    req(input$dash_climate, input$dash_3d_npoints)
    
    exp_seq <- seq(10, 50, length.out = input$dash_3d_npoints)
    reg_seq <- sort(unique(data_ts$region_index))
    reg_seq <- head(reg_seq, input$dash_3d_npoints)
    
    z_matrix <- outer(exp_seq, reg_seq, Vectorize(function(e, r) {
      predict_prob_occurrence(e, input$dash_climate, r, data_ts, samples_mat)
    }))
    
    plot_ly(
      x = exp_seq,
      y = reg_seq,
      z = ~z_matrix,
      type = "surface",
      colorscale = "Viridis"
    ) %>%
      layout(
        scene = list(
          xaxis = list(title = "Exposure"),
          yaxis = list(title = "Region Index"),
          zaxis = list(title = "Probability"),
          camera = list(
            eye = list(x = 1.87, y = 0.88, z = 0.64)
          )
        ),
        title = list(
          text = paste("3D Probability Surface -", input$dash_climate),
          font = list(size = 14)
        )
      )
  })
  
  ### Additional Analysis Tab
  
  # Coefficients Table
  output$coefTable <- renderDT({
    req(input$analysis_climate)
    
    coefs <- lapply(region_names, function(rn) {
      rn_idx <- unique(data_ts$region_index[data_ts$LHD == rn])[1]
      cvals <- get_median_coeff(input$analysis_climate, rn_idx, samples_mat)
      data.frame(
        Region = rn,
        'Median Coefficient' = round(mean(cvals), 3),
        'Standard Error' = round(sd(cvals), 3)
      )
    })
    final_df <- do.call(rbind, coefs)
    
    datatable(
      final_df,
      options = list(
        pageLength = 10,
        searching = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel'),
        initComplete = JS(
          "function(settings, json) {",
          "$(this.api().table().header()).css({'background-color': '#f8f9fa', 'color': '#2c3e50'});",
          "}"
        ),
        rowCallback = JS(
          "function(row, data) {",
          "  $(row).css({'border-left': '3px solid #3498db'});",
          "}"
        )
      ),
      class = 'cell-border stripe hover',
      rownames = FALSE,
      style = 'bootstrap4',
      selection = 'single'
    ) %>%
      formatStyle(
        columns = colnames(final_df),
        backgroundColor = '#ffffff',
        color = '#2c3e50'
      )
  })
  
  # Coefficient Plot
  output$coefPlot <- renderPlotly({
    req(input$analysis_climate)
    
    df_scatter <- data.frame(
      Region = region_names,
      MeanExposure = sapply(region_names, function(rn) {
        mean(data_ts[[input$analysis_climate]][data_ts$LHD == rn], na.rm = TRUE)
      }),
      MedianCoef = sapply(region_names, function(rn) {
        rn_idx <- unique(data_ts$region_index[data_ts$LHD == rn])[1]
        mean(get_median_coeff(input$analysis_climate, rn_idx, samples_mat))
      })
    )
    
    # Add pseudo "confidence intervals" for demonstration
    df_scatter$CI_lower <- df_scatter$MedianCoef - sapply(region_names, function(rn) {
      rn_idx <- unique(data_ts$region_index[data_ts$LHD == rn])[1]
      sd(get_median_coeff(input$analysis_climate, rn_idx, samples_mat))
    })
    df_scatter$CI_upper <- df_scatter$MedianCoef + sapply(region_names, function(rn) {
      rn_idx <- unique(data_ts$region_index[data_ts$LHD == rn])[1]
      sd(get_median_coeff(input$analysis_climate, rn_idx, samples_mat))
    })
    
    p <- ggplot(df_scatter, aes(
      x = MeanExposure, y = MedianCoef,
      text = paste(
        "Region:", Region,
        "\nMean Exposure:", round(MeanExposure, 2),
        "\nMedian Coefficient:", round(MedianCoef, 3)
      )
    )) +
      geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                    width = 0.2,
                    color = "#95a5a6",
                    alpha = 0.5) +
      geom_point(color = "#3498db", size = 4, alpha = 0.8) +
      geom_smooth(method = "loess", se = TRUE, color = "#e74c3c", fill = "#e74c3c", alpha = 0.2) +
      geom_text_repel(aes(label = Region), box.padding = 0.5, point.padding = 0.2, force = 2, size = 3) +
      labs(
        title = paste("Exposure-Coefficient Relationship -", input$analysis_climate),
        x = "Mean Exposure",
        y = "Median Coefficient"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 11),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "#ecf0f1"),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white")
      )
    
    ggplotly(p, tooltip = "text") %>%
      layout(
        hoverlabel = list(bgcolor = "white"),
        plot_bgcolor = "rgba(0,0,0,0)"
      )
  })
  
  # A dynamic observer to adjust the exposure slider range based on the chosen climate
  observe({
    climate_var <- input$dash_climate
    range_vals <- range(data_ts[[climate_var]], na.rm = TRUE)
    
    updateSliderInput(session, "dash_exposure",
                      min = floor(range_vals[1]),
                      max = ceiling(range_vals[2]),
                      value = round(mean(range_vals))
    )
  })
}

##############################
# 5. Run the App
##############################
shinyApp(ui = ui, server = server)
