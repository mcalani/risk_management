#Problem set 3
# Question 4 : 
# generate a set of simulated price paths for a stock with lognormally 
# distributed returns by using standardized discrete random variables z = +/- 1

Nt = 252
Np = 10000
dt = 1/252
S0 = 100
mu = 0.06
sigma = 0.4

#Simulation
r <- matrix(rnorm(Nt*Np,mean=mu*dt,sd=sigma*sqrt(dt)),nrow=Nt)

#Calibration of the binomial tree
x <- sign(r)
#p=sum(x>0)/(Nt*Np)
p=0.5

a = (mu*dt) - (sigma*sqrt(dt)) *sqrt(p/(1-p))
b = (sigma*sqrt(dt)) / sqrt(p*(1-p))

r_c <-  b*x/2 #calibrated shock  
r_c <- mean(r)-mean(r_c)+  b*x/2
mean(r_c)
mean(r)
sd(r_c)
sd(r)

S <- matrix(0,Nt+1,Np)
Sd <- matrix(0,Nt+1,Np)

for (k in 1:Nt) {
  #z ~ N(0,1)
  S[k+1,] <- S[k,]  + r[k,] # x[k,] #
  
  #z~ (+1,-1) 
  Sd[k+1,] <- Sd[k,]  + r_c[k,] # x[k,] #
}

matplot(S[,1:100], type = "l", lty=1,ylab="Price", xlab="Time",
        main="Sample lognormal price path");grid()
matplot(Sd[,1:100], type = "l", lty=1,ylab="Price", xlab="Time",
        main="Sample lognormal price path");grid()

hist(S[Nt+1,]) #normal
hist(Sd[Nt+1,]) #normal


St=S0*exp(S)
Sdt=S0*exp(Sd)

matplot(St[,1:100], type = "l", lty=1,ylab="Price", xlab="Time",
        main="Sample lognormal price path");grid()
matplot(Sdt[,1:100], type = "l", lty=1,ylab="Price", xlab="Time",
        main="Sample lognormal price path");grid()


hist(St[Nt+1,],breaks=50)
hist(Sdt[Nt+1,],breaks=50)

St.terminal =St[Nt+1,]
Sdt.terminal =Sdt[Nt+1,]

#Part A----------------------------------------------
mean(St[Nt+1,]) #right : 116.02
sd(St[Nt+1,])   #right : 48.37

mean(Sdt[Nt+1,]) #
sd(Sdt[Nt+1,])   #


#Part B----------------------------------------------
#how many times on average S_t cross or hit S_0 ?
crossings <- function(P) {
  # Count number of times each path returns to initial value
  # P is a matrix whose columns are time series sample paths
  # x takes value +1 for a crossing occurring between t-1 and t
  # Exact returns to the origin take value 1/2 (and will be counted twice)
  Nt <- nrow(P)
  x <- 0 * P # Initialize matrix of crossing counts
  for (t in 3:Nt) { # Begin after initial step from the origin
    x[t, ] = (1 - sign(P[t, ]-P[1, ])*sign(P[t-1, ]-P[1, ]))/2
  }
  return(x)
}



matplot(St[,1:100], type = "l", lty=1,ylab="Price", xlab="Time",
        main="Sample lognormal price path");grid()

x <- crossings(S)
x.count <- apply(x,2,sum)

mean(x.count)

#Part C----------------------------------------------
K=100 #strike price
C <- (St - K)*(St > K) #Call option
Cd <- (Sdt - K)*(Sdt > K) #Call option
matplot(C[,1:100], type = "l", lty=1,ylab="Price", xlab="Time",
        main="European Call option price path");grid()
matplot(Cd[,1:100], type = "l", lty=1,ylab="Price", xlab="Time",
        main="European Call option price path");grid()

C.terminal=C[Nt+1,]
Cd.terminal=Cd[Nt+1,]

hist(C.terminal)
hist(Cd.terminal)




#Attempt1: 26.79 (wrong)
#Attempt2: 27.501 (wrong)
mean(C[dim(C)[1],])
mean(Cd[dim(Cd)[1],])

#Attempt3: 16.116 (wrong)
mean(C)
mean(Cd)

#Attempt4: 53.793
mean(C[C>0])
mean(Cd[Cd>0])

#Attempt5: 25.73 (right)
mean(C[Nt,])
mean(Cd[Nt,])

#Part D----------------------------------------------
K=100 #strike price
P <- (K-St)*(St < K) #Put option
Pd <- (K-Sdt)*(Sdt < K) #Put option

matplot(P[,1:100], type = "l", lty=1,ylab="Price", xlab="Time",
        main="Sample lognormal price path");grid()
matplot(Pd[,1:100], type = "l", lty=1,ylab="Price", xlab="Time",
        main="Sample lognormal price path");grid()

#Attempt1: 11.10 (wrong)
mean(P[Nt+1,])
mean(Pd[Nt+1,])

hist(P[Nt+1,],breaks = 80)

#Attempt2: 17.294 (wrong)
mean(P[P>0])
mean(Pd[Pd>0])


#Attempt3: 
mean(P[Nt,])
mean(Pd[Nt,])
