#Analysis for Grant, Krishnan, and Bradbury

#Analysis of Teresa Blader 260th street data from 2015 to get stage survival probabilities

library(runjags)
runjags.options(force.summary=TRUE)

#maximum time to pupation,
pupamax = 16
#number of days in the study (46) + pupamax (16) = 64
n = 46 + pupamax
#mean temperature each day of the study period of length n days - from Ames WSW weather station
MT = c(21.9,21.9,20.8,17.5,18.1,19.7,22.2,24.7,17.2,16.4,18.9,20.3,23.3,26.9,28.3,26.4,
       24.4,23.65,26.7,26.7,24.15,23.6,20.8,20.55,21.95,23.6,26.1,24.4,24.75,22.5,21.95,21.95,22.75,23.6,25.25,23.6,20.85,21.7,21.4,23.85,24.45,24.15,
       22.5,21.7,22.2,21.7,24.4,23.6,23.6,24.15,19.75,15.8,17.2,20.25,20.85,17.5,15.25,15.85,16.65,18.05,18.9,18.35)
length(MT)


#calculate accumulated day-degrees for each cohort
ADD = matrix(nrow=n+pupamax, ncol=n); colnames(ADD) = 1:n
for (j in 1:n){
  for (i in (j+1):(n+1)) {
    ADD[i,j] = sum(MT[j:(i-1)])
  }
}
#set upper right off diagonal cells to 0 - these have to be set to 0 for JAGS (this code was originally in the model)
for (j in 1:n) {
  for (i in 1:j){
    ADD[i,j] = 0    
  }
}
#set lower cells to 0
for (j in 1:n) {
  for (i in (n+2):(n+pupamax)){
    ADD[i,j] = 0    
  }
}

#calculate stage durations for each cohort (this code was originally in the JAGS model, which is why it is a string of ifelse() statements instead of something more elegant)
SDC = matrix(nrow = 7, ncol = n); colnames(SDC) = 1:n; rownames(SDC) = c("egg","first","second","third","fourth","fifth","pupa")
for(j in 1:n) {
  SDC[1,j] = 1
  SDC[2,j] = ifelse(ADD[j+2,j] > 45,3,
                    ifelse(ADD[j+3,j] > 45,4,
                           5))
  SDC[3,j] = ifelse(ADD[j+3,j] > 77.3,4,
                    ifelse(ADD[j+4,j] > 77.3,5,
                           ifelse(ADD[j+5,j] > 77.3,6,
                                  7)))
  SDC[4,j] = ifelse(ADD[j+4,j] > 105.1,5,
                    ifelse(ADD[j+5,j] > 105.1,6,
                           ifelse(ADD[j+6,j] > 105.1,7,
                                  ifelse(ADD[j+7,j] > 105.1,8,
                                         9))))
  SDC[5,j] = ifelse(ADD[j+5,j] > 129.6,6,
                    ifelse(ADD[j+6,j] > 129.6,7,
                           ifelse(ADD[j+7,j] > 129.6,8,
                                  ifelse(ADD[j+8,j] > 129.6,9,
                                         10))))
  SDC[6,j] = ifelse(ADD[j+6,j] > 165.3,7,
                    ifelse(ADD[j+7,j] > 165.3,8,
                           ifelse(ADD[j+8,j] > 165.3,9,
                                  ifelse(ADD[j+9,j] > 165.3,10,
                                         ifelse(ADD[j+10,j] > 165.3,11,
                                                ifelse(ADD[j+11,j] > 165.3,12,
                                                       13))))))
  SDC[7,j] = ifelse(ADD[j+9,j] > 231.9,10,
                    ifelse(ADD[j+10,j] > 231.9,11,
                           ifelse(ADD[j+11,j] > 231.9,12,
                                  ifelse(ADD[j+12,j] > 231.9,13,
                                         ifelse(ADD[j+13,j] > 231.9,14,
                                                ifelse(ADD[j+14,j] > 231.9,15,
                                                       ifelse(ADD[j+15,j] > 231.9,16,
                                                              17)))))))
}


#########  calculate stage durations for each cohort and mean stage durations  ###############
#vector of mean stage durations
MSD = vector(length = 6)
#cohorts that finished 5th instar
I5 = ADD[(n+1),] > 231.9
#remove first 16 days before study started
I5[1:16] = FALSE
#stage durations of cohorts that reached pupation
SDC[7,I5]
#fifth instar stage durations for each cohort
SD5 = SDC[7,I5]-SDC[6,I5]
MSD[6] = mean(SD5); range(SD5)
#cohorts that finished 4th instar
I4 = ADD[(n+1),] > 165.3
I4[1:16] = FALSE
MSD[5] = mean(SDC[6,I4]-SDC[5,I4])
#cohorts that finished 3rd instar
I3 = ADD[(n+1),] > 129.6; I3[1:16] = FALSE
MSD[4] = mean(SDC[5,I3]-SDC[4,I3])
#cohorts that finished 2nd instar
I2 = ADD[(n+1),] > 105.1; I2[1:16] = FALSE
MSD[3] = mean(SDC[4,I2]-SDC[3,I2])
#cohorts that finished 1st instar
I1 = ADD[(n+1),] > 77.3; I1[1:16] = FALSE
MSD[2] = mean(SDC[3,I1]-SDC[2,I1])
#cohorts that hatched
I0 = ADD[(n+1),] > 45; I0[1:16] = FALSE
MSD[1] = mean(SDC[2,I0]-SDC[1,I0])
MSD

#Field Counts - counts were made approximately once per week for several weeks
#counts were taken on day 1 and the last count was taken on day 46
C = matrix(nrow = 46, ncol = 6)
C[1,] = c(37,2,4,0,0,0)
C[9,] = c(90,15,9,7,1,1)
C[17,] = c(135,19,5,5,3,7)
C[24,] = c(156,22,7,3,1,3)
C[31,] = c(121,5,2,3,1,3)
C[38,] = c(77,1,7,2,1,1)
C[46,] = c(11,1,0,2,3,2)


#To allow the model to estimate eggs laid before the first day of the study, 16 days are added before the first day of the study
Buffer = matrix(nrow = pupamax, ncol = 6)
CB = rbind(Buffer,C)

#run the model with runjags package
#data inputs for JAGS
data = list(n=n, pupamax=pupamax, MSD = MSD, SDC = SDC, ADD = ADD, C=CB)
#initial vlaues for lambda2 and S
inits = list(list(lambda2 = 200*runif(n), S = runif(6)),
             list(lambda2 = 200*runif(n), S = runif(6)),
             list(lambda2 = 200*runif(n), S = runif(6)),
             list(lambda2 = 200*runif(n), S = runif(6)),
             list(lambda2 = 200*runif(n), S = runif(6)),
             list(lambda2 = 200*runif(n), S = runif(6)),
             list(lambda2 = 200*runif(n), S = runif(6)))

#runjags, monitoring daily survival, stage survival, cumulative survival (SC), and eggs laid B.  
out.260.np = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 7, data = data, inits = inits)
out.260.np.e = extend.jags(out.260.np.e, sample = 100000)
plot(out.260.np.e, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))
write.csv(out.260.np.e$summaries, "out.260.np.e.csv")
out.260.np.e.sum1 = out.260.np.e$summaries
out.260.np.e.sum2 = out.260.np.e$summary
save(out.260.np.e, file = "out.260.np.e.RData")
rm(out.260.np.e)
write.csv(out.260.np.e.sum1, file = "out.260.np.e.sum1.csv")


#run with de Anda priors - using MonarchModelAnda2.bug
out.260 = run.jags(model="MonarchModelAnda2.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 7, data = data, inits = inits)
plot(out.260, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))
#add some more samples
out.260.e = extend.jags(out.260, sample = 90000)
#added 100,000 more
out.260.e = extend.jags(out.260.e, sample = 100000)
plot(out.260.e, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))
write.csv(out.260.e$summaries, "out.260.e.csv")
out.260.e.sum1 = out.260.e$summaries
out.260.e.sum2 = out.260.e$summary













######################### FIGURES AND GRAPHS ############################

#graph priors

#sequence to calculate density at
q = seq(from = 0, to = 1, length.out = 51)
sd = 0.5
var = sd^2
precision = 1/var
precision #precision for JAGS


plot(q, dnorm(q, mean = 0.63, sd = sd), type = "l", xlab = "Survival Probability Density", ylab = "", ylim = c(0,max(dnorm(q, mean = 0.63, sd = sd))), 
     main = "Egg Survival Prior Distribution")

plot(q, dnorm(q, mean = 0.60, sd = sd), type = "l", xlab = "Survival Probability Density", ylab = "", ylim = c(0,max(dnorm(q, mean = 0.60, sd = sd))), 
     main = "First Instar Survival Prior Distribution")

plot(q, dnorm(q, mean = 0.56, sd = sd), type = "l", xlab = "Survival Probability Density", ylab = "", ylim = c(0,max(dnorm(q, mean = 0.56, sd = sd))), 
     main = "Second Instar Survival Prior Distribution")




#graph posteriors and priors

library(coda)
library(mcmcplots)
library(runjags)

out.260.e.MCMC = as.mcmc(out.260.e) #says it combined the 7 chains together

#huge file so load to make graphs then delete from workspace
load("out.260.e.MCMC.RData")

denplot(out.260.e.MCMC, parms = "S[1]", xlim = c(0,1), xlab = "Daily Egg Survival", ylab = "Density", main = "", auto.layout = FALSE, col = "Black") 
lines(q, 3.1*(dnorm(q, mean = 0.63, sd = sd)), lty = 'dotted')

denplot(out.260.e.MCMC, parms = "S[2]", xlim = c(0,1), xlab = "Daily First Instar Survival", ylab = "Density", main = "", auto.layout = FALSE, col = "Black") 
lines(q, 3.5*(dnorm(q, mean = 0.60, sd = sd)), lty = 'dotted')

denplot(out.260.e.MCMC, parms = "S[3]", xlim = c(0,1), xlab = "Daily Second Instar Survival", ylab = "Density", main = "", auto.layout = FALSE, col = "Black") 
lines(q, 2.25*(dnorm(q, mean = 0.56, sd = sd)), lty = 'dotted')

rm(out.260.e.MCMC)


library(ggmcmc)

out.260.gg = ggs(out.260.e.MCMC)
str(out.260.gg)
#ggmcmc(out.260.gg) this command crashed my computer by using up ALL the RAM

save(out.260.e.MCMC, file = "out.260.e.MCMC.RData")
rm(out.260.gg)
rm(out.260.e.MCMC)

#load("out.260.e.MCMC.RData")














