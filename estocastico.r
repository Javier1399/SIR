simulate_emod <- function(S0, I0, R0, beta, gamma, timesteps, vaccination_rate, sd_factor, mask_factor, 
                          age_factor_mean, heart_disease_factor_mean, diabetes_factor_mean, 
                          lung_disease_factor_mean, obesity_factor_mean, kidney_disease_factor_mean) {
  
  # Initialize vectors to store results
  S <- numeric(timesteps + 1)
  I <- numeric(timesteps + 1)
  R <- numeric(timesteps + 1)
  Deaths <- numeric(timesteps + 1)  # Track the number of deaths
  p_infection <- numeric(timesteps + 1)
  p_recovery <- numeric(timesteps + 1)
  
  # Initialize health conditions for each individual
  age_factor <- rbinom(S0, 1, age_factor_mean)  # Randomly assign age vulnerability
  heart_disease_factor <- rbinom(S0, 1, heart_disease_factor_mean)  # Randomly assign heart disease
  diabetes_factor <- rbinom(S0, 1, diabetes_factor_mean)  # Randomly assign diabetes
  lung_disease_factor <- rbinom(S0, 1, lung_disease_factor_mean)  # Randomly assign lung disease
  obesity_factor <- rbinom(S0, 1, obesity_factor_mean)  # Randomly assign obesity
  kidney_disease_factor <- rbinom(S0, 1, kidney_disease_factor_mean)  # Randomly assign kidney disease
  
  S[1] <- S0
  I[1] <- I0
  R[1] <- R0
  p_infection[1] <- NA
  p_recovery[1] <- NA
  Deaths[1] <- 0  # No deaths initially
  
  # Run the simulation over discrete time steps
  for (t in 1:timesteps) {
    # Apply random variability in social distancing and mask usage
    distance_effect <- rnorm(1, mean = 1, sd = sd_factor)  # Social distancing effect
    mask_effect <- rnorm(1, mean = 1, sd = mask_factor)  # Mask usage effect
    
    # Adjust the transmission rate with these factors
    adjusted_beta <- beta * distance_effect * mask_effect
    
    # Calculate the probability of infection and recovery
    p_infection[t] <- 1 - exp(-adjusted_beta * I[t] / (S[t] + I[t] + R[t]))
    p_recovery[t] <- 1 - exp(-gamma)
    
    # Calculate the number of new infections and recoveries
    new_infections <- rbinom(1, S[t], p_infection[t])
    new_recoveries <- rbinom(1, I[t], p_recovery[t])
    
    vaccinated <- rbinom(1, S[t], vaccination_rate)
    
    # Calculate the number of deaths based on random comorbidity presence
    # Calculate probability of death for infected individuals
    # The sum of all disease risk factors is used to calculate the death probability
    death_prob <- rnorm(new_infections, mean = 0, sd = 0.05)  # Base probability of death (without comorbidities)
    
    # Apply each disease risk factor by adding their effects if the person has the condition
    death_prob <- death_prob + (age_factor[1:new_infections] * 0.3)  # Add risk if in vulnerable age
    death_prob <- death_prob + (heart_disease_factor[1:new_infections] * 0.5)  # Heart disease risk
    death_prob <- death_prob + (diabetes_factor[1:new_infections] * 0.4)  # Diabetes risk
    death_prob <- death_prob + (lung_disease_factor[1:new_infections] * 0.6)  # Lung disease risk
    death_prob <- death_prob + (obesity_factor[1:new_infections] * 0.7)  # Obesity risk
    death_prob <- death_prob + (kidney_disease_factor[1:new_infections] * 0.5)  # Kidney disease risk
    
    # Apply the death probability (with bounds between 0 and 1)
    death_prob <- pmin(pmax(death_prob, 0), 1)  # Ensure it's between 0 and 1
    
    # Calculate the deaths today based on the death probability for each new infection
    deaths_today <- sum(rbinom(new_infections, 1, death_prob))  # Number of deaths today
    
    # Update the populations
    S[t + 1] <- S[t] - new_infections - vaccinated
    I[t + 1] <- I[t] + new_infections - new_recoveries - deaths_today
    R[t + 1] <- R[t] + new_recoveries
    Deaths[t + 1] <- Deaths[t] + deaths_today
  }
  
  # Return a data frame with results
  data.frame(Time = 0:timesteps, 
             S = S, 
             I = I, 
             R = R, 
             Deaths = Deaths, 
             P_Infection = round(p_infection, 4), 
             P_Recovery = round(p_recovery, 4))
}

# Establecer parámetros
S0 <- 990       # Número inicial de individuos susceptibles
I0 <- 10        # Número inicial de infectados
R0 <- 0         # Número inicial de recuperados
beta <- 0.3     # Tasa de transmisión
gamma <- 0.1    # Tasa de recuperación
vaccination_rate <- 0.05  # Tasa de vacunación
timesteps <- 200 # Número de pasos de tiempo

# Factores de riesgo (probabilidades entre 0 y 1)
age_factor_mean <- 0.02  # 20% de la población tiene vulnerabilidad por edad
heart_disease_factor_mean <- 0.03  # 30% de la población tiene enfermedades cardíacas
diabetes_factor_mean <- 0.05  # 25% de la población tiene diabetes
lung_disease_factor_mean <- 0.02  # 20% de la población tiene enfermedades pulmonares
obesity_factor_mean <- 0.05  # 35% de la población tiene obesidad
kidney_disease_factor_mean <- 0.01  # 10% de la población tiene enfermedad renal

# Ejecutar la simulación
results <- simulate_emod(S0, I0, R0, beta, gamma, timesteps, vaccination_rate, 
                         sd_factor = 0.2, mask_factor = 0.2, 
                         age_factor_mean, heart_disease_factor_mean, diabetes_factor_mean, 
                         lung_disease_factor_mean, obesity_factor_mean, kidney_disease_factor_mean)

# Graficar los resultados
library(ggplot2)
ggplot(results, aes(x = Time)) +
  geom_line(aes(y = S, color = "Susceptible")) +
  geom_line(aes(y = I, color = "Infected")) +
  geom_line(aes(y = R, color = "Recovered")) +
  geom_line(aes(y = Deaths, color = "Deaths")) +
  labs(y = "Population", color = "Compartment") +
  theme_minimal() +
  ggtitle("EMOD SIR Model with Vaccination, Random Comorbidities, and Mortality")

