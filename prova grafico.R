# Required libraries
library(dplyr)
library(ggplot2)

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
nsim = 100 #Number of Monte Carlo simulations

#Initialization of the vectors for the trajectories
S = matrix(NA, nsim, N+1)
v = matrix(NA, nsim, N+1)

#Setting initial conditions
S[,1] = 131.23 #Initial price of the underlying asset
v[,1] = v0
set.seed(1)

#Trajectories simulations
for (i in 1:N){
  dW1 = rnorm(nsim)*sqrt(dt)
  dW2 = rho*dW1+sqrt(1-rho^2)*rnorm(nsim)*sqrt(dt)
  
  dv = k*(theta-v[,i])*dt+beta*sqrt(v[,i])*dW2
  dS = r*S[,i]*dt+sqrt(v[,i])*S[,i]*dW1
  
  v[,i+1] = pmax(v[,i]+dv, 0)
  S[,i+1] = S[,i]+dS
}

#Pricing an arithmetic Asian option with K = 131.23
K = 131.23
A = apply(S, 1, FUN = mean)
discounted_payoff= exp(-r*Time)*pmax(A-K,0)

#MC price
MCP = mean(discounted_payoff)
sd_MCP= sd(discounted_payoff)
error_MCP = sd_MCP/sqrt(nsim)

# Calculate the 95% confidence interval
margin = qnorm(0.975)*error_MCP
confidence_interval_lower = MCP - margin
confidence_interval_upper = MCP + margin

# End measuring time
end_time = Sys.time()

# Calculate the computation time
computation_time = end_time - start_time
cat("Computation Time:", computation_time, "seconds\n")

#Start measuring time
start_time_graph = Sys.time()

#Graphing underlying asset trajectories
price_simulation_data <- data.frame(
  Time = rep(seq(0, Time, by = dt), each = nsim),
  Trajectory = rep(1:nsim, times = N + 1),
  Price = as.vector(S))

ggplot(data = price_simulation_data, aes(x = Time, y = Price, group = Trajectory, color = as.factor(Trajectory))) +
  geom_line() +
  labs(x = "Time", y = "Price", title = "Simulated Electricity Price Trajectories (MWh)") +
  theme_minimal() + theme(legend.position = "none")

# End measuring time
end_time_graph = Sys.time()

# Calculate the computation time
computation_time_graph = end_time_graph - start_time_graph
cat("Computation Time:", computation_time_graph, "seconds\n")



margin = qt(0.975,df= nsim-1)
margin2 = qnorm(0.975)
