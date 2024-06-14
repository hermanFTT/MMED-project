############################ mini PROJECT MMED ###############################################

########## Load the Data ############

library(deSolve)
library(tidyverse)

load("WAevddat.Rdata")
# Function to partition the data
partition_data <- function(data) {
  liberia <- data[1:810,]
  sierra <- data[811:1734,]
  guinea <- data[1735:4869,]
  return(list(liberia = liberia, sierra = sierra, guinea = guinea))
}

# Function to summarize the data
summarize_data <- function(data) {
  summarized_data <- data %>%
    group_by(date) %>%
    summarize(total_cases = sum(cases))
  return(summarized_data)
}

partitioned_data = partition_data(evddat)

Liberia <- summarize_data(partitioned_data$liberia) 

Sierra <- summarize_data(partitioned_data$sierra)

Guinea <- summarize_data(partitioned_data$guinea)


####### SEIR for Ebola ( without vaccination )

####################################### Model With Vaccination Component #############################

seir_v <- function(t,y,params){
  S <- y[1]
  E <- y[2]
  I <- y[3]
  R <- y[4]
  
  beta <- params["beta"]
  #N <- params["N"]
  eps <- params["eps"]
  w <-params["w"] 
  sigma <- params["sigma"]
  delta <- params["delta"]
  rho <-params["rho"]
  phi=params["phi"]
  w_i <- params["w_i"]
  alpha <- params["alpha"]
  
  dSdt <- -beta*S*I -eps*S*R-w*S+alpha*R-phi*S+sigma
  dEdt <- beta*S*I + eps*S*R - delta*E - w*E-phi*E
  dIdt <- delta*E - (rho+w_i)*I
  dRdt <- rho*I-alpha*R+phi*S+phi*E
  
  return(list(c(dSdt,dEdt,dIdt,dRdt)))}

## R_o , beta , E* , S* (equilibruim)

  beta2_comp <- function(R_o=1.83,rho=0.1,delta,w=0.073,w_i=0.10,phi=0.5){R_o*(rho+w_i)*(delta+w+phi)/delta}
  S_star2 <- param.vals['sigma']/param.vals['w']
  R_star2 <- param.vals['phi']*param.vals['sigma']/param.vals['w']*param.vals['alpha']
  
  
# parameters setting
  N0=700
  
  # parameters
  # β = 0.0003589, ε = 0.000799, δ = 0.0179
  param.vals <- c(
    delta=0.083, 
    beta=beta_compute(delta=0.083), # 0.71,
    N=N0,
    eps=0.089,    #0.089,
    sigma=1.7,  # 1.7
    rho=0.3,     # 0.3
    w=0.073, # 0.073
    phi=0.6, # vaccination rate ( set it according to the chosen vaccination program)
    w_i=0.1,
    alpha=3 #2.57
  )
  S.star1=param.vals["sigma"]/param.vals["w"]
  # over two years 
  times <- seq(160,2*365,7)
  S0 <- 500; E0 <-50 ; I0 <- 5
  init <- c(sus=S0,exp=E0,inf=I0,rec=0)
  
######### Solve the EDO
  tc2 <- data.frame(lsoda(init,times,seir_v,param.vals))
  
############# models  plot
  
  plot(seq(1,nrow(Liberia),1),Liberia$total_cases,main= "Cases in Liberia",xlab="time(weeks)",ylab="Number of cases")
  lines(lines(tc2$time/7,(tc2$exp+tc2$inf),col="green",lty=1) )
  legend("topright", legend=c( "vaccine (model)"," without vaccine (model)"), col=c("green","red"), lty=c(1,3), bty="n")
  
#############
inf_no_vac <- tc2$inf # save the infectious (without vaccine)
inf_vac <- tc2$inf    # save the infectious ( with vaccine)  
sum_diff <- sum(inf_no_vac-inf_vac) # averted sum ( sum of diff )

################## Adverted infected result for each vaccination program (0.2, 0.6, 0.9) 

Sl <- c(15,40,56)          # results in Sierra Leone
G <- c(13,34,47)           # results in Guinea
L <- c(16,44,61)           # results in Liberia
vac <- c(0.2,0.6,0.9)      # vaccination program 

R_vec <- c(2.02,2.02,2.02,1.71,1.71,1.71,1.83,1.83,1.83)  # concatenation of "Ro's" 
cout_vec <- c(15,40,56,13,34,47,16,44,61)                 # conctenation of results averted infected numbers

########### graph
plot( R_vec,cout_vec, xaxt="n",
     ylab = " infected cases averted", 
     xlab = "Ro (basic reproduction number)", 
     pch = 19,   # Solid circle for points
     cex = 1.5,
     #xlim=c(1.71,1.83,2.02),
     ylim=c(12,67),# Increase point size (1.5 times the default size)
)
axis(1,at=c(1.71,1.83,2.02),labels=c("1.71 (Guinea)","1.83 (Liberia)","2.02 (Sierra Leone)"))

# Add lines to join the points
lines(c(1.71,1.83,2.02), c(47,61,56), type = "b", pch = 19, lty = 1,col="green")
lines(c(1.71,1.83,2.02), c(34,44,40), type = "b", pch = 19, lty = 1,col="orange")
lines(c(1.71,1.83,2.02), c(13,16,15), type = "b", pch = 19, lty = 1,col="red")
legend("topright", legend=c("vacc rate=90%", "vacc rate = 60%", "vacc rate = 20%"), col=c("green", "orange", "red"), lty=1, bty="n")


########################### change of R_o with respect to the vaccination coverage rate #######################

#delta <- c(0.099,0.074,0.066) # for Sierra Leone, Liberia, Guinea respectively
R_o <- function(phi,rho=0.1,beta=0.71,delta=0.066,w=0.07,w_i=0.1){ delta*beta/((rho+w_i)*(delta+w+phi))}

phi_val <- seq(0,1,0.1)
vec_R <- sapply(phi_val,R_o)

  
plot(phi_val, vec_R, 
     xlab = "vaccination rate", 
     ylab = "Ro (basic reproduction number)", 
     pch = 16,   # Solid circle for points
     #cex = 1.5,
     col='blue',
     ylim=c(0,2.5),# Increase point size (1.5 times the default size)
)

# Add lines to join the points
lines(phi_val, vec_R, type = "b", pch = 16, lty = 1,col="blue")
points(phi_val, vec_R, type="b",pch=16,lty=1,col='orange')
points(phi_val, vec_R, type = "b", pch = 16, lty = 1, col='green')

abline(h=1, lty=2, col="red")

# Add thresholds  (vertical dashed line)
abline(v=0.1, lty=3, col="green")
abline(v=0.13, lty=3, col="orange")
abline(v=0.2, lty=3, col="blue")
axis(1,at=c(0.1,0.13),labels=c(0.1,0.13))

legend("topright", legend=c("Sierra", "Guinea", "Liberia"), col=c( "blue", "green","orange"), lty=1, bty="n")

