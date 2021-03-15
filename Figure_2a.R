library(Rfast)
library(igraph)
library(ggplot2)
library(reshape2)


###### PARAMETER SETUP -----------------

### Social network parameters

# Population of institution
N = 10000

# Mean number of workplace contacts
mu = 2.3

# Variance in number of workplace contacts
sigma2 = 2.4

# Mean number of peripheral contacts
mu_star = 0.2

# Variance of outside contacts
sigma2_star = 1.9

# Proportion of outside contacts known to testing authorities
omega <- 1


### Viral dynamics parameters

# Secondary attack rate - mean
mu_rho = 0.16

# Secondary attack rate - variance
sigma2_rho = 0.01

# Recovery rate
gamma = 1/10

# Exposed-to-infectious transition rate
delta = 1/5

# Initial infection rate within the institution
v0 = rep(0.01, N)

# Initial infection rate in the periphery
V0_star = 0.06
E0_star <- V0_star * gamma / delta


### Diagnostic testing parameters

# Number of tests per person per day
C = 0.12

# Test sensitivity (E stage)
psi_e = 0.2

# Test sensitivity (V stage)
psi_v = 1

# Proportion of tests used inside the institution (5 possible scenarios)
p = c(0, 0.25, 0.5, 0.75, 1)


###### CONTACT NETWORK GENERATION -----------------



# Solve NBin(r, p) parameters via MoM
NBin_solve <- function(mean, var){
  return(
    c(
      mean^2/(var - mean),
      mean/var
    )
  )
}

target_contacts <- rgamma(N, mu^2/(sigma2 - mu), mu/(sigma2 - mu))
mtx <- outer(target_contacts, target_contacts)/sum(target_contacts)
diag(mtx) = 0

to.assign <- runif(N * (N-1) / 2) < upper_tri(mtx)
mtx <- upper_tri.assign(mtx, to.assign)
mtx <- t(mtx)
mtx <- upper_tri.assign(mtx, to.assign)

mean(colsums(mtx))
var(colsums(mtx))


# Contacts outside the institution
y <- rnbinom(
  N, 
  NBin_solve(mu_star, sigma2_star)[1], 
  NBin_solve(mu_star, sigma2_star)[2]
)

# Known contacts outside the institution
z <- rbinom(
  N,
  y,
  omega
)

known_periphery_size <- sum(z)
total_known_contacts <- colsums(mtx) + z

M <- N + sum(y)

mtx2 <- matrix(0, nrow = N, ncol = M)
mtx2[1:N, 1:N] <- mtx
count <- N + 1
known_periphery <- c()
for (i in 1:N) {
  if(y[i] > 0){
    mtx2[i, count:(count+y[i] - 1)] <- 1
    if(z[i] > 0){
      known_periphery <- c(known_periphery, count:(count+z[i] - 1))
    }
  }
  count <- count + y[i]
}

# SAR 
rho <- rbeta(
  M,
  ((mu_rho * (1 - mu_rho) / sigma2_rho) - 1)*mu_rho,
  ((mu_rho * (1 - mu_rho) / sigma2_rho) - 1)*(1 - mu_rho)
)

# Attack rate (in terms of other parameters)
beta = 1 - exp(-gamma*rho/(1-rho))

###### DIFFERENTIAL EQUATION SOLVING -----------------

# Number of days
n_days = 41

# Lists of states matrices for each replication
S_list <- list()
E_list <- list()
V_list <- list()
U_list <- list()
W_list <- list()
Q_list <- list()
R_list <- list()

# Optimization function used to compute eta
adjust <- function(k, vec){
  return((sum(1 - exp(-k*vec)) - sum(vec))^2)
}


### Solve differential equations under each value of p

for (t in 1:5) {
  
  # Number of tests assigned to the institution
  n_institution_test <- C * N * p[t]
  
  # Number of tests assigned to the periphery
  n_periphery_test <- C * N * (1 - p[t])
  
  # Rate at which periphery members get tested
  eta_star <- min(n_periphery_test/known_periphery_size, 1)
  
  # Steady state probability of being exposed in the periphery
  E_infty_star <- (gamma*V0_star)/(eta_star*psi_e + delta)
  
  # Steady state probability of being infectious in the periphery
  V_infty_star <- (delta*gamma*V0_star)/((eta_star*psi_v + gamma) * (eta_star*psi_e + delta))
  
  
  # States matrices. Entry j,k is probability that on day j, agent k is in the given state
  S <- matrix(nrow = n_days, ncol = N)
  S[1, ] <- 1 - v0
  
  E <- matrix(nrow = n_days, ncol = N)
  E[1, ] <- v0
  
  V <- matrix(nrow = n_days, ncol = N)
  V[1, ] <- rep(0, N)
  
  U <- matrix(nrow = n_days, ncol = N)
  U[1, ] <- rep(0, N)
  
  W <- matrix(nrow = n_days, ncol = N)
  W[1, ] <- rep(0, N)
  
  Q <- matrix(nrow = n_days, ncol = N)
  Q[1, ] <- rep(0,N)
  
  R <- matrix(nrow = n_days, ncol = N)
  R[1, ] <- rep(0,N)
  
  
  ### Populate matrices by looping through all relevant days
  
  for (j in 2:n_days) {
    
    # Compute eta, the probability of being tested inside the institution 
    eta <- total_known_contacts * (S[j-1, ] + E[j-1, ] + V[j-1, ])
    eta <- n_institution_test * eta / sum(eta)
    alpha <- optimize(adjust, c(1,10), tol = 0.01, vec = eta)$minimum
    eta <- 1 - exp(-alpha*eta)
    
    p_infectious <- c(V[j-1, ], rep(V0_star, sum(y)))
    p_infectious[known_periphery] <- V_infty_star
    
    BETA <- 1 - exp(-((mtx2 %*% (p_infectious*beta))))
    GAMMA <- 1 - exp(-((mtx %*% ((E[j-1, ]*psi_e + V[j-1, ]*psi_v) * eta)) + ((E_infty_star*psi_e + V_infty_star*psi_v) * z * eta_star)))
    
    # Update matrices as per differential equations
    S[j, ] <- S[j-1, ] - 
      S[j-1, ] * BETA -
      S[j-1, ] * (
        GAMMA
      ) +
      eta * U[j-1, ]
    
    E[j, ] <- E[j-1, ] - 
      delta * E[j-1, ] -
      eta * psi_e * E[j-1, ] +
      S[j-1, ] * BETA -
      E[j-1, ] * (
        GAMMA
      )
    
    V[j, ] <- V[j-1, ] - 
      gamma * V[j-1, ] -
      eta * psi_v* V[j-1, ] + 
      delta*E[j-1, ] -
      V[j-1, ] * (
        GAMMA
      )
    
    U[j, ] <- U[j-1, ] -
      eta * U[j-1, ] +
      S[j-1, ] * (
        GAMMA
      )
    
    W[j, ] <- W[j-1, ] -
      delta * W[j-1, ] + 
      eta * psi_e * E[j-1, ] + 
      E[j-1, ] * (
        GAMMA
      )
    
    Q[j, ] <- Q[j-1, ] -
      gamma * Q[j-1, ] + 
      delta * W[j-1, ] +
      eta * psi_v * V[j-1, ] +
      V[j-1, ] * (
        GAMMA
      )
    
    R[j, ] <- R[j-1,] +
      gamma * V[j-1, ] + 
      gamma * Q[j-1, ]
  }
  
  # Populate lists
  S_list[[t]] <- S
  E_list[[t]] <- E
  V_list[[t]] <- V
  U_list[[t]] <- U
  W_list[[t]] <- W
  Q_list[[t]] <- Q
  R_list[[t]] <- R
  
}

plot_data <- data.frame(day = 1:n_days)
for (t in 1:5) {
  totals <- rowsums(E_list[[t]] + V_list[[t]] + W_list[[t]] + Q_list[[t]] + R_list[[t]])
  plot_data <- cbind(plot_data,  totals)
}
colnames(plot_data) <- c("day", "V1", "V2", "V3", "V4", "V5")
dd = melt(plot_data, id=c("day"))

colors <- colorRampPalette(c("#3A1C71", "#D76D77", "#FFAF7B"))(5)
labels <- c()
for (i in 1:5) {
  labels[i] <- paste("$p = ", p[i], sep = "") 
}

plot <- ggplot(dd) + geom_line(aes(x=day, y=value, color=variable)) + 
  scale_color_manual(values=colors,
                     name="",
                     labels = labels
  ) + 
  labs(x = "Day", y = "Cumulative Cases, Institution") +
  theme_minimal()

print(plot)

