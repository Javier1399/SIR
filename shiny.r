library(shiny)
library(ggplot2)
library(deSolve)

# Función para el modelo EMOD SIR
simulate_emod <- function(S0, I0, R0, beta, gamma, timesteps) {
  S <- numeric(timesteps + 1)
  I <- numeric(timesteps + 1)
  R <- numeric(timesteps + 1)
  
  S[1] <- S0
  I[1] <- I0
  R[1] <- R0
  
  for (t in 1:timesteps) {
    new_infections <- rbinom(1, S[t], 1 - exp(-beta * I[t] / (S[t] + I[t] + R[t])))
    new_recoveries <- rbinom(1, I[t], 1 - exp(-gamma))
    
    S[t + 1] <- S[t] - new_infections
    I[t + 1] <- I[t] + new_infections - new_recoveries
    R[t + 1] <- R[t] + new_recoveries
  }
  
  data.frame(Time = 0:timesteps, S = S, I = I, R = R)
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
                 sliderInput("S0_emod", "Susceptibles Iniciales (S0):", 100, 10000, 990, step = 10),
                 sliderInput("I0_emod", "Infectados Iniciales (I0):", 1, 100, 10, step = 1),
                 sliderInput("R0_emod", "Recuperados Iniciales (R0):", 0, 1000, 0, step = 10),
                 sliderInput("beta_emod", "Tasa de Transmisión (β):", 0.01, 1, 0.3, step = 0.01),
                 sliderInput("gamma_emod", "Tasa de Recuperación (γ):", 0.01, 1, 0.1, step = 0.01),
                 sliderInput("timesteps_emod", "Duración de la Simulación (T):", 10, 365, 100, step = 5)
               ),
               mainPanel(
                 plotOutput("plot_emod"),
                 textOutput("summary_emod")
               )
             )
    ),
    # Pestaña del modelo SIR basado en ecuaciones diferenciales
    tabPanel("Modelo Diferencial SIR",
             sidebarLayout(
               sidebarPanel(
                 sliderInput("S0_sir", "Susceptibles Iniciales (S0):", 0, 1, 0.999, step = 0.001),
                 sliderInput("I0_sir", "Infectados Iniciales (I0):", 0, 1, 0.001, step = 0.001),
                 sliderInput("R0_sir", "Número Básico de Reproducción (R0):", 1, 10, 4, step = 0.1),
                 sliderInput("gamma_sir", "Tasa de Recuperación (γ):", 0.01, 1, 0.1, step = 0.01),
                 sliderInput("mu_sir", "Tasa de Cambio de Población (μ):", 0, 0.1, 0, step = 0.01),
                 sliderInput("timesteps_sir", "Duración de la Simulación (T):", 10, 365, 100, step = 5)
               ),
               mainPanel(
                 plotOutput("plot_sir"),
                 textOutput("summary_sir")
               )
             )
    )
  )
)

# Server de la aplicación Shiny
server <- function(input, output) {
  # Modelo EMOD SIR
  output$plot_emod <- renderPlot({
    results <- simulate_emod(
      S0 = input$S0_emod,
      I0 = input$I0_emod,
      R0 = input$R0_emod,
      beta = input$beta_emod,
      gamma = input$gamma_emod,
      timesteps = input$timesteps_emod
    )
    
    ggplot(results, aes(x = Time)) +
      geom_line(aes(y = S, color = "Susceptible")) +
      geom_line(aes(y = I, color = "Infected")) +
      geom_line(aes(y = R, color = "Recovered")) +
      labs(
        title = "Modelo EMOD SIR",
        y = "Población",
        color = "Compartment"
      ) +
      theme_minimal()
  })
  
  output$summary_emod <- renderText({
    paste("Parámetros EMOD SIR:",
          "\nS0 =", input$S0_emod,
          "| I0 =", input$I0_emod,
          "| R0 =", input$R0_emod,
          "\nβ =", input$beta_emod,
          "| γ =", input$gamma_emod,
          "| T =", input$timesteps_emod)
  })
  
  # Modelo Diferencial SIR
  output$plot_sir <- renderPlot({
    init <- c(S = input$S0_sir, I = input$I0_sir, R = 1 - input$S0_sir - input$I0_sir)
    parms <- c(mu = input$mu_sir, N = 1, R0 = input$R0_sir, gamma = input$gamma_sir)
    t <- seq(0, input$timesteps_sir, by = 1)
    sir_output <- ode(y = init, times = t, func = sir, parms = parms)
    sir_df <- as.data.frame(sir_output)
    
    ggplot(sir_df, aes(x = time)) +
      geom_line(aes(y = S, color = "Susceptible")) +
      geom_line(aes(y = I, color = "Infectados")) +
      geom_line(aes(y = R, color = "Recuperados")) +
      labs(
        title = "Modelo Diferencial SIR",
        y = "Proporción de Población",
        x = "Días",
        color = "Compartment"
      ) +
      theme_minimal()
  })
  
  output$summary_sir <- renderText({
    paste("Parámetros SIR Diferencial:",
          "\nS0 =", input$S0_sir,
          "| I0 =", input$I0_sir,
          "\nR0 =", input$R0_sir,
          "| γ =", input$gamma_sir,
          "| μ =", input$mu_sir,
          "| T =", input$timesteps_sir)
  })
}

# Ejecuta la aplicación Shiny
shinyApp(ui = ui, server = server)
