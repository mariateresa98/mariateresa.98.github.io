# Required libraries
library(dplyr)
library(ggplot2)
library(gridExtra)

#Start measuring time
start_time = Sys.time()

#Definition of Heston model's parameters
r = 0.0433 # risk-free interest rate
k = 2 # rate of reversion
theta = 0.04 # long-run volatility
beta = 0.3 #volatility of volatility
rho = -0.5 #correlation parameter between the underlying asset and its volatility
v0 = 0.04 # initial condition for volatility

#Definition of Monte Carlo simulation parameters
Time = 1 # Maturity
N = 100 # Number of steps
dt = Time/N # Time increment
nsim = 5000 #Number of Monte Carlo simulations

# Initialization of the vectors for the trajectories for two sets of paths
S1 = matrix(NA, nsim/2, N+1)
v1 = matrix(NA, nsim/2, N+1)
S2 = matrix(NA, nsim/2, N+1)
v2 = matrix(NA, nsim/2, N+1)

# Setting initial conditions for both sets
S1[,1] = 131.23
v1[,1] = v0
S2[,1] = 131.23
v2[,1] = v0
set.seed(1)

# Trajectories simulations for antithetic variates
for (i in 1:N) {
  dW1 = rnorm(nsim/2) * sqrt(dt)
  dW2 = rho * dW1 + sqrt(1 - rho^2) * rnorm(nsim/2) * sqrt(dt)
  
  dv1 = k * (theta - v1[,i]) * dt + beta * sqrt(v1[,i]) * dW2
  dv2 = k * (theta - v2[,i]) * dt + beta * sqrt(v2[,i]) * (-dW2)  # Antithetic paths
  
  dS1 = r * S1[,i] * dt + sqrt(v1[,i]) * S1[,i] * dW1
  dS2 = r * S2[,i] * dt + sqrt(v2[,i]) * S2[,i] * (-dW1)  # Antithetic paths
  
  v1[,i+1] = pmax(v1[,i] + dv1, 0)
  v2[,i+1] = pmax(v2[,i] + dv2, 0)
  S1[,i+1] = S1[,i] + dS1
  S2[,i+1] = S2[,i] + dS2
}

# Pricing an arithmetic Asian option with K = 131.23 using antithetic variates
K = 131.23
A1 = apply(S1, 1, FUN = mean)
A2 = apply(S2, 1, FUN = mean)
discounted_payoff1 = exp(-r * Time) * pmax(A1 - K, 0) #Regular path
discounted_payoff2 = exp(-r * Time) * pmax(A2 - K, 0) #Antithetic path

# Average the discounted payoffs from both sets
MCP_antithetic = 0.5 * (mean(discounted_payoff1) + mean(discounted_payoff2))
corr= cor(discounted_payoff1, discounted_payoff2)
sd_MCP = sqrt(1/4*(var(discounted_payoff1 + discounted_payoff2)))
error_MCP = sd_MCP/sqrt(nsim/2)

#MC error under the hypothesis of independence
error_ind = sqrt(1/(4*nsim/2)*(var(discounted_payoff1) + var(discounted_payoff2)))

# Calculate the 95% confidence interval
margin = qnorm(0.975)*error_MCP
confidence_interval_lower = MCP_antithetic - margin
confidence_interval_upper = MCP_antithetic + margin

# End measuring time
end_time = Sys.time()

# Calculate the computation time
computation_time = end_time - start_time
cat("Computation Time:", computation_time, "seconds\n")

#Start measuring time
start_time_graph = Sys.time()

#Graphing underlying asset antithetic trajectories
# Create a data frame for the regular and antithetic trajectories of the underlying asset
regular_trajectories = data.frame(
  Time = rep(seq(0, Time, by = dt), each = nsim/2),
  Trajectory = rep(1:nsim/2, times = N + 1),
  Price = as.vector(S1)
)

antithetic_trajectories = data.frame(
  Time = rep(seq(0, Time, by = dt), each = nsim/2),
  Trajectory = rep(1:nsim/2, times = N + 1),
  Price = as.vector(S2)
)

# Create line plots for the regular and antithetic trajectories
regular_plot = ggplot(data = regular_trajectories, aes(x = Time, y = Price, group = Trajectory, color = as.factor(Trajectory))) +
  geom_line() +
  labs(title = "Regular Electricity Price Trajectories (MWh)", x = "Time", y = "Price")+
  theme_minimal() + theme(legend.position = "none")

antithetic_plot = ggplot(data = antithetic_trajectories, aes(x = Time, y = Price, group = Trajectory, color = as.factor(Trajectory))) +
  geom_line() +
  labs(title = "Antithetic Electricity Price Trajectories (MWh)", x = "Time", y = "Price")+
  theme_minimal() + theme(legend.position = "none")

# Print the plots side by side
grid.arrange(regular_plot, antithetic_plot, ncol = 2)

# End measuring time
end_time_graph = Sys.time()

# Calculate the computation time
computation_time_graph = end_time_graph - start_time_graph
cat("Computation Time:", computation_time_graph, "seconds\n")
