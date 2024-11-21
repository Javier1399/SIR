library(shiny)
library(shinyWidgets)
library(ggplot2)
library(deSolve)
library(DT)  # Para tablas interactivas

# Función para el modelo EMOD SIR
simulate_emod <- function(S0, I0, R0, beta, gamma, timesteps) {
  vaccination_rate = 0.01
  mask_factor = 0.2
  sd_factor = 0.2
  # Initialize vectors to store results
  S <- numeric(timesteps + 1)
  I <- numeric(timesteps + 1)
  R <- numeric(timesteps + 1)
  p_infection <- numeric(timesteps + 1)
  p_recovery <- numeric(timesteps + 1)
  
  # Set initial values
  S[1] <- S0
  I[1] <- I0
  R[1] <- R0
  p_infection[1] <- NA
  p_recovery[1] <- NA
  
  # Run the simulation over discrete time steps
  for (t in 1:timesteps) {
    # Apply random variability in social distancing and mask usage
    distance_effect <- rnorm(1, mean = 1, sd = sd_factor)  # Random effect for social distancing
    mask_effect <- rnorm(1, mean = 1, sd = mask_factor)  # Random effect for mask usage
    
    # Adjust the transmission rate with these factors
    adjusted_beta <- beta * distance_effect * mask_effect
    
    # Calculate probabilities of infection and recovery
    p_infection[t] <- 1 - exp(-adjusted_beta * I[t] / (S[t] + I[t] + R[t]))
    p_recovery[t] <- 1 - exp(-gamma)
    
    # Calculate the number of new infections and recoveries
    new_infections <- rbinom(1, S[t], p_infection[t])
    new_recoveries <- rbinom(1, I[t], p_recovery[t])
    
    vaccinated <- rbinom(1, S[t], vaccination_rate)
    
    # Update the populations
    S[t + 1] <- S[t] - new_infections-vaccinated
    I[t + 1] <- I[t] + new_infections - new_recoveries
    R[t + 1] <- R[t] + new_recoveries
  }
  
  # Return a data frame with results
  data.frame(Time = 0:timesteps, 
             S = S, 
             I = I, 
             R = R, 
             P_Infection = round(p_infection, 4), 
             P_Recovery = round(p_recovery, 4))
}


# Función para el modelo basado en ecuaciones diferenciales
sir <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    beta <- R0 * gamma  # Tasa de infección dependiente de gamma
    dS <- mu * N - beta * I * S / N - mu * S
    dI <- beta * I * S / N - gamma * I - mu * I
    dR <- gamma * I - mu * R
    return(list(c(dS, dI, dR)))
  })
}

# UI de la aplicación Shiny
ui <- fluidPage(
  titlePanel("Modelos SIR Interactivos"),
  tabsetPanel(
    # Pestaña del modelo EMOD SIR
    tabPanel("Modelo EMOD SIR",
             sidebarLayout(
               sidebarPanel(
                 sliderInput("S0_emod_1", "Susceptibles Iniciales (S0 - Caso 1):", 100, 10000, 990, step = 10),
                 sliderInput("I0_emod_1", "Infectados Iniciales (I0 - Caso 1):", 1, 100, 10, step = 1),
                 sliderInput("R0_emod_1", "Recuperados Iniciales (R0 - Caso 1):", 0, 1000, 0, step = 10),
                 sliderInput("beta_emod_1", "Tasa de Transmisión (β - Caso 1):", 0.01, 1, 0.3, step = 0.01),
                 sliderInput("gamma_emod_1", "Tasa de Recuperación (γ - Caso 1):", 0.01, 1, 0.1, step = 0.01),
                 sliderInput("timesteps_emod", "Duración de la Simulación (T):", 10, 365, 100, step = 5),
                 checkboxInput("enable_case2_emod", "Activar Caso 2", FALSE),
                 
                 # Para el segundo caso en el modelo EMOD
                 conditionalPanel(
                   condition = "input.enable_case2_emod == true",
                   sliderInput("S0_emod_2", "Susceptibles Iniciales (S0 - Caso 2):", 100, 10000, 990, step = 10),
                   sliderInput("I0_emod_2", "Infectados Iniciales (I0 - Caso 2):", 1, 100, 10, step = 1),
                   sliderInput("R0_emod_2", "Recuperados Iniciales (R0 - Caso 2):", 0, 1000, 0, step = 10),
                   sliderInput("beta_emod_2", "Tasa de Transmisión (β - Caso 2):", 0.01, 1, 0.3, step = 0.01),
                   sliderInput("gamma_emod_2", "Tasa de Recuperación (γ - Caso 2):", 0.01, 1, 0.1, step = 0.01)
                 ),
                 
                 # Selección para mostrar tablas
                 selectInput("table_case_emod", "Seleccionar Tabla a Mostrar:",
                             choices = c("Caso 1", "Caso 2"),
                             selected = "Caso 1")
               ),
               mainPanel(
                 plotOutput("plot_emod_comparison"),
                 textOutput("summary_emod"),
                 DT::dataTableOutput("table_emod_comparison")
               )
             )
    ),
    # Pestaña del modelo SIR basado en ecuaciones diferenciales
    tabPanel("Modelo Diferencial SIR",
             sidebarLayout(
               sidebarPanel(
                 sliderInput("S0_sir_1", "Susceptibles Iniciales (S0 - Caso 1):", 0, 1, 0.999, step = 0.001),
                 sliderInput("I0_sir_1", "Infectados Iniciales (I0 - Caso 1):", 0, 1, 0.001, step = 0.001),
                 sliderInput("R0_sir_1", "Número Básico de Reproducción (R0 - Caso 1):", 1, 10, 4, step = 0.1),
                 sliderInput("gamma_sir_1", "Tasa de Recuperación (γ - Caso 1):", 0.01, 1, 0.1, step = 0.01),
                 sliderInput("mu_sir_1", "Tasa de Cambio de Población (μ - Caso 1):", 0, 0.1, 0, step = 0.01),
                 sliderInput("timesteps_sir", "Duración de la Simulación (T):", 10, 365, 100, step = 5),
                 checkboxInput("enable_case2_sir", "Activar Caso 2", FALSE),
                 
                 # Para el segundo caso en el modelo SIR
                 conditionalPanel(
                   condition = "input.enable_case2_sir == true",
                   sliderInput("S0_sir_2", "Susceptibles Iniciales (S0 - Caso 2):", 0, 1, 0.999, step = 0.001),
                   sliderInput("I0_sir_2", "Infectados Iniciales (I0 - Caso 2):", 0, 1, 0.001, step = 0.001),
                   sliderInput("R0_sir_2", "Número Básico de Reproducción (R0 - Caso 2):", 1, 10, 4, step = 0.1),
                   sliderInput("gamma_sir_2", "Tasa de Recuperación (γ - Caso 2):", 0.01, 1, 0.1, step = 0.01),
                   sliderInput("mu_sir_2", "Tasa de Cambio de Población (μ - Caso 2):", 0, 0.1, 0, step = 0.01)
                 ),
                 
                 # Selección para mostrar tablas
                 selectInput("table_case_sir", "Seleccionar Tabla a Mostrar:",
                             choices = c("Caso 1", "Caso 2"),
                             selected = "Caso 1")
               ),
               mainPanel(
                 plotOutput("plot_sir_comparison"),
                 textOutput("summary_sir"),
                 DT::dataTableOutput("table_sir_comparison")
               )
             )
    )
  )
)



server <- function(input, output, session) {
  
  # Modelo EMOD SIR - Comparación de los dos casos
  output$plot_emod_comparison <- renderPlot({
    # Caso 1
    results_1 <- simulate_emod(
      S0 = input$S0_emod_1,
      I0 = input$I0_emod_1,
      R0 = input$R0_emod_1,
      beta = input$beta_emod_1,
      gamma = input$gamma_emod_1,
      timesteps = input$timesteps_emod
    )
    results_1$Case <- "Caso 1"
    
    combined_results <- results_1
    
    if (input$enable_case2_emod) {
      # Caso 2
      results_2 <- simulate_emod(
        S0 = input$S0_emod_2,
        I0 = input$I0_emod_2,
        R0 = input$R0_emod_2,
        beta = input$beta_emod_2,
        gamma = input$gamma_emod_2,
        timesteps = input$timesteps_emod
      )
      results_2$Case <- "Caso 2"
      combined_results <- rbind(results_1, results_2)
    }
    
    ggplot(combined_results, aes(x = Time, color = Case)) +
      geom_line(aes(y = S, linetype = "Susceptible")) +
      geom_line(aes(y = I, linetype = "Infected")) +
      geom_line(aes(y = R, linetype = "Recovered")) +
      labs(
        title = "Comparación del Modelo EMOD SIR",
        y = "Población",
        color = "Caso",
        linetype = "Compartimiento"
      ) +
      theme_minimal()
  })
  
  output$table_emod_comparison <- renderDataTable({
    # Caso 1
    results_1 <- simulate_emod(
      S0 = input$S0_emod_1,
      I0 = input$I0_emod_1,
      R0 = input$R0_emod_1,
      beta = input$beta_emod_1,
      gamma = input$gamma_emod_1,
      timesteps = input$timesteps_emod
    )
    results_1$Case <- "Caso 1"
    
    # Caso 2
    results_2 <- simulate_emod(
      S0 = input$S0_emod_2,
      I0 = input$I0_emod_2,
      R0 = input$R0_emod_2,
      beta = input$beta_emod_2,
      gamma = input$gamma_emod_2,
      timesteps = input$timesteps_emod
    )
    results_2$Case <- "Caso 2"
    
    combined_results <- rbind(results_1, results_2)
    
    # Filtrar la tabla según el caso seleccionado
    if (input$table_case_emod == "Caso 1") {
      combined_results <- subset(combined_results, Case == "Caso 1")
    } else {
      combined_results <- subset(combined_results, Case == "Caso 2")
    }
    
    DT::datatable(combined_results, options = list(pageLength = 10))
  })
  
  # Modelo SIR Diferencial - Comparación de los dos casos
  output$plot_sir_comparison <- renderPlot({
    # Caso 1
    init_1 <- c(S = input$S0_sir_1, I = input$I0_sir_1, R = 1 - input$S0_sir_1 - input$I0_sir_1)
    parms_1 <- c(mu = input$mu_sir_1, N = 1, R0 = input$R0_sir_1, gamma = input$gamma_sir_1)
    t <- seq(0, input$timesteps_sir, by = 1)
    sir_output_1 <- ode(y = init_1, times = t, func = sir, parms = parms_1)
    sir_df_1 <- as.data.frame(sir_output_1)
    sir_df_1$Case <- "Caso 1"
    
    combined_sir_results <- sir_df_1
    
    if (input$enable_case2_sir) {
      # Caso 2
      init_2 <- c(S = input$S0_sir_2, I = input$I0_sir_2, R = 1 - input$S0_sir_2 - input$I0_sir_2)
      parms_2 <- c(mu = input$mu_sir_2, N = 1, R0 = input$R0_sir_2, gamma = input$gamma_sir_2)
      sir_output_2 <- ode(y = init_2, times = t, func = sir, parms = parms_2)
      sir_df_2 <- as.data.frame(sir_output_2)
      sir_df_2$Case <- "Caso 2"
      combined_sir_results <- rbind(sir_df_1, sir_df_2)
    }
    
    ggplot(combined_sir_results, aes(x = time, color = Case)) +
      geom_line(aes(y = S, linetype = "Susceptible")) +
      geom_line(aes(y = I, linetype = "Infected")) +
      geom_line(aes(y = R, linetype = "Recovered")) +
      labs(
        title = "Comparación del Modelo Diferencial SIR",
        y = "Proporción de Población",
        color = "Caso",
        linetype = "Compartimiento"
      ) +
      theme_minimal()
  })
  
  output$table_sir_comparison <- renderDataTable({
    # Caso 1
    init_1 <- c(S = input$S0_sir_1, I = input$I0_sir_1, R = 1 - input$S0_sir_1 - input$I0_sir_1)
    parms_1 <- c(mu = input$mu_sir_1, N = 1, R0 = input$R0_sir_1, gamma = input$gamma_sir_1)
    t <- seq(0, input$timesteps_sir, by = 1)
    sir_output_1 <- ode(y = init_1, times = t, func = sir, parms = parms_1)
    sir_df_1 <- as.data.frame(sir_output_1)
    sir_df_1$Case <- "Caso 1"
    
    # Caso 2
    init_2 <- c(S = input$S0_sir_2, I = input$I0_sir_2, R = 1 - input$S0_sir_2 - input$I0_sir_2)
    parms_2 <- c(mu = input$mu_sir_2, N = 1, R0 = input$R0_sir_2, gamma = input$gamma_sir_2)
    sir_output_2 <- ode(y = init_2, times = t, func = sir, parms = parms_2)
    sir_df_2 <- as.data.frame(sir_output_2)
    sir_df_2$Case <- "Caso 2"
    
    combined_sir_results <- rbind(sir_df_1, sir_df_2)
    
    # Filtrar la tabla según el caso seleccionado
    if (input$table_case_sir == "Caso 1") {
      combined_sir_results <- subset(combined_sir_results, Case == "Caso 1")
    } else {
      combined_sir_results <- subset(combined_sir_results, Case == "Caso 2")
    }
    
    DT::datatable(combined_sir_results, options = list(pageLength = 10))
  })
}




shinyApp(ui = ui, server = server)
                                 
