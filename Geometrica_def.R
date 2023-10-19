library(NMOF)

#Definition of Heston model's parameters
r = 0.0433 # risk-free interest rate
k = 2 # rate of reversion
theta = 0.04 # long-run volatility
beta = 0.3 # volatility of volatility
rho = -0.5 # correlation parameter between the underlying asset and its volatility
v0 = 0.04 # initial condition for volatility
s0 = 131.23 # initial condition for asset price
K = 131.23 # strike price
Time = 1 # maturity
t = 0 # initial time
tau = 1 # Time-t

# Define function z0
z0 = function(s, w, t, Time, St, r, beta, k, theta, rho, v_t) {
  result = s * (((Time - t) / Time) * log(St) +
    (((r * beta - k * theta * rho) * (Time - t)^2) / (2 * beta * Time)) -
    (rho * (Time - t) / (beta * Time)* v_t)) + w * (log(St)-(rho / beta)*v_t + 
    (r - ((rho * k * theta)/ beta))*(Time - t))
  return(result)
}

# Define functions z1, z2, z3, and z4 for complex numbers s w
z1 = function(s, w, rho, beta, Time) {
  return((s^2 * (1 - rho^2)) / (2 * Time^2))
}

z2 = function(s, w, rho, beta, Time) {
  return( s*(2 * rho * k - beta) / (2 * beta * Time) + s * w * (1 - rho^2) / Time)
}

z3 = function(s, w, rho, beta, Time) {
  return((s * rho) / (beta * Time) + w * (2 * rho * k - beta) / (2 * beta ) + (w^2 * (1 - rho^2)) / 2)
}

z4 = function(s, w, rho, beta) {
  return(w * rho / (beta))
}

# Define the recursive functions f_n 
fn = function(n, s, w, tau, beta, k, rho) {
  if (Re(n) == -2 || Re(n) == -1) {
    return(0)
  } else if (Re(n) == 0) {
    return(1)
  } else if (Re(n) == 1) {
    return((k - z4(s, w, rho, beta)*beta^2) * tau / 2) 
  } else if (Re(n) >= 2) {
    sum_term <- z1(s, w, rho, beta, tau) * tau^2 * fn(n - 4, s, w, tau, beta, k, rho) +
      z2(s, w, rho, beta, tau) * tau * fn(n - 3, s, w, tau, beta, k, rho) +
      (z3(s, w, rho, beta, tau)-k^2/(2*beta^2)) * fn(n - 2, s, w, tau, beta, k, rho)
    return(-((beta^2 * tau^2) / (2 * Re(n) * (Re(n) - 1))) * sum_term)
  } else {
    return(0)
  }
}

# Define the function F_tau with a specified cut-off term "N_cutoff"
F_tau = function(s, w, tau, beta, k, rho, N_cutoff) {
  result = 0
  for (n in 0:N_cutoff) {
    fn_value = fn(n, s, w, tau, beta, k, rho)
    result = result + fn_value
  }
  return(result)
}

# Define the function F_tau_approx with a specified cut-off term "N_cutoff"
F_tau_approx = function(s, w, tau, beta, k, rho, N_cutoff) {
  result = 0
  for (n in 1:N_cutoff) {
    fn_value = fn(n, s, w, tau, beta, k, rho)
    result = result + (n / tau) * fn_value
  }
  return(result)
}

# Vectorization
z0_v = Vectorize(z0, vectorize.args = "s")
F_tau_v = Vectorize(F_tau, vectorize.args = "s")
F_tau_approx_v = Vectorize(F_tau_approx, vectorize.args = "s")

# Calculate the characteristic function
characteristic_function= function(om,S, X, tau, r, v0, theta, rho, k, beta, N_cutoff){
  om <- 1i*om
  z0_value = z0_v(om, 0, 0, tau, S, r, beta, k, theta, rho, v0)
  F_tau_value = F_tau_v(s=om, w=0, tau, beta, k, rho, N_cutoff)
  F_tau_approx_value = F_tau_approx_v(s=om, w=0, tau, beta, k, rho, N_cutoff)
  
  res = exp(z0_value) * exp(
    ((k * v0 + k^2 * theta * tau) / beta^2) -
      (2 * v0 / beta^2) * (F_tau_approx_value / F_tau_value) -
      (2 * k * theta / beta^2) * log(F_tau_value))
  return(res)
}

characteristic_function(om=0,S=131.23, X=0, tau=1, r=0, v0, theta, rho, k, beta, N_cutoff=15)

#Price the geometric Asian option
Geometric_Heston_price = callCF(characteristic_function, S=s0, X=K, tau=Time, r=r, q = 0, v0, theta, rho, k, beta, N_cutoff=15, implVol = FALSE)
