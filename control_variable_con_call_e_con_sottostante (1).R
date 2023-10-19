# Required libraries
library(NMOF)

#Start measuring time
start_time = Sys.time()

#Definition of Heston model's parameters
r = 0.0433 # risk-free interest rate
k = 2 # rate of reversion
theta = 0.04 # long-run volatility
beta = 0.3 #volatility of volatility
rho = -0.5 #correlation parameter between the underlying asset and its volatility
v0 = 0.04 # initial condition for volatility

#Definition of simulation parameters
Time = 1 # Maturity
N = 100 # Number of steps
dt = Time/N # Time increment
nsim = 10000 #Number of Monte Carlo simulations

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

#Pricing an arithmetic Asian call option with K = 131.23
K = 131.23
A = apply(S, 1, FUN = mean)
discounted_payoff_asian= exp(-r*Time)*pmax(A-K,0)

#Pricing a European call option with K = 131.23
discounted_payoff_european= exp(-r*Time)*pmax(S[,N+1]-K,0)
Mcp=mean(discounted_payoff_european)

#Pricing underlying asset
discounted_payoff_under= exp(-r*Time)*S[,N+1]

#Pricing the control variable using the explicit formula for the Heston model
Heston_price_european = callHestoncf(S = 131.23, X = K, tau = Time, r = r, q = 0, v0 = v0, vT = theta, rho = rho, k = k, sigma = beta, implVol = FALSE)

#Error European call
Error = (discounted_payoff_european-Heston_price_european)

# Error  underlying asset
Error_under = discounted_payoff_under-S[,1]

#Computing the adjustment factor lambda for the European call 
cov_asian_european = cov(discounted_payoff_asian,Error)
var_european = var(Error)
lambda = -cov_asian_european / var_european 

#Computing the adjustment factor lambda for the underlying asset 
cov_asian_under = cov(discounted_payoff_asian,Error_under)
var_under = var(Error_under)
lambda_under = -cov_asian_under / var_under 

# Pricing the Asian option using the European call as control variable
MCP_Asian = mean(discounted_payoff_asian+lambda*Error)
sd_MCP_withreduction = sqrt(var(discounted_payoff_asian)-cov_asian_european^2/var_european)
error_MCP_withreduction = sd_MCP_withreduction/sqrt(nsim)

sd_MCP_withoutreduction = sqrt(var(discounted_payoff_asian))
error_MCP_withoutreduction = sd_MCP_withoutreduction/sqrt(nsim)

# Calculate the 95% confidence interval
margin = qnorm(0.975)*error_MCP_withreduction
confidence_interval_lower = MCP_Asian - margin
confidence_interval_upper = MCP_Asian + margin

# Pricing the Asian option using the underlying asset as control variable
MCP_Asian_under = mean(discounted_payoff_asian+lambda_under*Error_under)
sd_MCP_withreduction_under = sqrt(var(discounted_payoff_asian)-cov_asian_under^2/var_under)
error_MCP_withreduction_under = sd_MCP_withreduction_under/sqrt(nsim)

# End measuring time
end_time = Sys.time()

# Calculate the computation time
computation_time = end_time - start_time
cat("Computation Time:", computation_time, "seconds\n")
