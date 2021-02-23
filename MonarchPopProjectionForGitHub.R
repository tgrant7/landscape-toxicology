#Simulate/Project monarch population from eggs to adults

####### Table of Contents  ##############################
#   1.  Simulate pop to insecticide drift point - L17
#   2.  Calculate size of population in drift zone - L298
#   3.  Population produced with only natural survival - L744
#   4.  Population produced with natural and insecticide mortality - L1028

#########################################################

library(ggplot2)
library(scales)
library(plot3D)

######################## 1.  Simulate population that gets hit by insecticide using survival est code  ##################################

#load survival estimates from StageSurvival.R
out.260.e.sum1 = read.csv("out.260.e.csv")

#simulation parameters
#number of days in simulation
n=20
#nuisance variable for simulation
pupamax = 16
#daily temperature
MT = rep(21,n) #same temp every day for simulations
#daily survival probabilities
S = out.260.e.sum1$Median
#total eggs laid from scenario 3
TE3 = 2379294
#when do eggs from first day reach mid 5th instar stage? try 14 days first.  So 14 days of eggs
d = 12
#number of eggs laid per day.  Should be of length n
B = rep(TE3/d, d);B[(d+1):n] = 0  
B



#code from SimulateandAnalyzeMonarchCounts.R
ADD = matrix(nrow=n+pupamax, ncol=n); colnames(ADD) = 1:n
for (j in 1:n){
  for (i in (j+1):(n+1)) {
    ADD[i,j] = sum(MT[j:(i-1)])
  }
}
#set upper right off diagonal cells to 0 - set to 0 so following sum operations will work smoothly
for (j in 1:n) {
  for (i in 1:j){
    ADD[i,j] = 0    
  }
}
#set lower cells to 0 as well
for (j in 1:n) {
  for (i in (n+2):(n+pupamax)){
    ADD[i,j] = 0    
  }
}

#calculate stage durations for each cohort. 
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


#########  calculate mean stage durations  ###############
#vector of mean stage durations
MSD = vector(length = 6)
#cohorts that finished 5th instar
I5 = ADD[(n+1),] > 231.9
#stage durations of cohorts that reached pupation
SDC[7,I5]
#fifth instar stage durations for each cohort
SD5 = SDC[7,I5]-SDC[6,I5]
MSD[6] = mean(SD5)
#cohorts that finished 4th instar
I4 = ADD[(n+1),] > 165.3
MSD[5] = mean(SDC[6,I4]-SDC[5,I4])
#cohorts that finished 3rd instar
I3 = ADD[(n+1),] > 129.6
MSD[4] = mean(SDC[5,I3]-SDC[4,I3])
#cohorts that finished 2nd instar
I2 = ADD[(n+1),] > 105.1
MSD[3] = mean(SDC[4,I2]-SDC[3,I2])
#cohorts that finished 1st instar
I1 = ADD[(n+1),] > 77.3
MSD[2] = mean(SDC[3,I1]-SDC[2,I1])
#cohorts that hatched
I0 = ADD[(n+1),] > 45
MSD[1] = mean(SDC[2,I0]-SDC[1,I0])
MSD

#calculate stage survival probabilities and cumulative survival probability
S1 = S[1]^MSD[1] 
S2 = S[2]^MSD[2]
S3 = S[3]^MSD[3]
S4 = S[4]^MSD[4]
S5 = S[5]^MSD[5]
S6 = S[6]^MSD[6]
SC = S1*S2*S3*S4*S5*S6 #cumulative survival probability
SN = S1*S2*S3*S4*S5 #survival probability most similar to Nail et al. 2015 estimator
SC


#PEDD = Proportion surviving each day for each cohort.  Rows are days, cols are cohorts.  This is p_ik*pi_j-k-bik from Grant et al.
PEDD = matrix(nrow=n+pupamax, ncol=n); colnames(PEDD) = 1:n
for(j in 1:n){
  PEDD[j,j] = 1
  for (i in 2:SDC[2,j]){
    PEDD[i+j-1,j] = (S[1])^(i-1)
  }
  for (i in (SDC[2,j]+1):SDC[3,j]){
    PEDD[i+j-1,j] = (S[1])^(SDC[2,j]-1)*(S[2])^(i-SDC[2,j])
  }
  for (i in (SDC[3,j]+1):SDC[4,j]){
    PEDD[i+j-1,j] = (S[1])^(SDC[2,j]-1)*(S[2])^(SDC[3,j]-SDC[2,j])*S[3]^(i-SDC[3,j])
  }
  for (i in (SDC[4,j]+1):SDC[5,j]){
    PEDD[i+j-1,j] = (S[1])^(SDC[2,j]-1)*(S[2])^(SDC[3,j]-SDC[2,j])*S[3]^(SDC[4,j]-SDC[3,j])*S[4]^(i-SDC[4,j])
  }
  for (i in (SDC[5,j]+1):SDC[6,j]){
    PEDD[i+j-1,j] = (S[1])^(SDC[2,j]-1)*(S[2])^(SDC[3,j]-SDC[2,j])*S[3]^(SDC[4,j]-SDC[3,j])*S[4]^(SDC[5,j]-SDC[4,j])*S[5]^(i-SDC[5,j])
  }
  for (i in (SDC[6,j]+1):SDC[7,j]){
    PEDD[i+j-1,j] = (S[1])^(SDC[2,j]-1)*(S[2])^(SDC[3,j]-SDC[2,j])*S[3]^(SDC[4,j]-SDC[3,j])*S[4]^(SDC[5,j]-SDC[4,j])*S[5]^(SDC[6,j]-SDC[5,j])*S[6]^(i-SDC[6,j])
  }
}
#set NA's to 0
for (j in 2:n) {
  for (i in 1:(j-1)) {
    PEDD[i,j] = 0
  }
}
for (j in 1:(n-1)) {
  for (i in (SDC[7,j]+j):(n+pupamax)) {
    PEDD[i,j] = 0
  }
}
#write.csv(PEDD, "PEDD.csv")


#Number surviving each day of each cohort.  rows are days, cols are cohorts.  These are the cohort contributions that are summed to get l_ij of Grant et al.  
NDEDD = matrix(nrow = n+pupamax, ncol = n); colnames(NDEDD) = 1:n
for (j in 1:n) {
  for (i in j:(j+SDC[7,j]-1)){
    NDEDD[i,j] = B[j]*PEDD[i,j]
  }
}
for (j in 2:n) {
  for (i in 1:(j-1)) {
    NDEDD[i,j] = 0
  }
}
for (j in 1:(n-1)) {
  for (i in (j+SDC[7,j]):(n+pupamax)){
    NDEDD[i,j] = 0
  }
}

#ST - index indicating what stage each cohort is in each day.  rows are days, cols are cohorts.  99 means it has pupated
#this is essentially the indicator function I(s_kj = i) of Grant et al. 
ST = matrix(nrow=n, ncol=n); colnames(ST) = 1:n
for (j in 1:(n-1)) {
  for (i in (j+1):n) {
    ST[i,j] = ifelse(ADD[i,j] >= 1 && ADD[i,j] < 45, 1, 
                     ifelse(ADD[i,j] >= 45 && ADD[i,j] < 77.3, 2,
                            ifelse(ADD[i,j] >= 77.3 && ADD[i,j] < 105.1, 3,
                                   ifelse(ADD[i,j] >= 105.1 && ADD[i,j] < 129.6, 4,
                                          ifelse(ADD[i,j] >= 129.6 && ADD[i,j] < 165.3, 5,
                                                 ifelse(ADD[i,j] >= 165.3 && ADD[i,j] < 231.9, 6, 99))))))
  }
}
for (j in 1:n){
  for (i in j:j){
    ST[i,j] = 1
  }
}
for (j in 2:n){
  for (i in 1:(j-1)){
    ST[i,j] = 0
  }
}

#A matrix for each stage - copy numbers from NDEDD using ST to tell which numbers go into which matrix
#i.e., if a ST[4,5] = 1, then NDEDD[4,5] goes into matrix M1
#this makes it easier/possible to sum the population across cohorts in each stage on any particular day
M1 = matrix(nrow=n, ncol=n); colnames(M1) = 1:n
M2 = matrix(nrow=n, ncol=n); colnames(M2) = 1:n
M3 = matrix(nrow=n, ncol=n); colnames(M3) = 1:n
M4 = matrix(nrow=n, ncol=n); colnames(M4) = 1:n
M5 = matrix(nrow=n, ncol=n); colnames(M5) = 1:n
M6 = matrix(nrow=n, ncol=n); colnames(M6) = 1:n

for (j in 1:n){
  for (i in 1:n) {
    M1[i,j] = ifelse(ST[i,j]==1,NDEDD[i,j],0)
    M2[i,j] = ifelse(ST[i,j]==2,NDEDD[i,j],0)
    M3[i,j] = ifelse(ST[i,j]==3,NDEDD[i,j],0)
    M4[i,j] = ifelse(ST[i,j]==4,NDEDD[i,j],0)
    M5[i,j] = ifelse(ST[i,j]==5,NDEDD[i,j],0)
    M6[i,j] = ifelse(ST[i,j]==6,NDEDD[i,j],0)
  }
}

#NDD = sum populations from each cohort into a matrix of population of each stage each day.  rows are days, cols are stages.  
#This is l_ij of Grant et al. 
NDD = matrix(nrow=n,ncol=6); colnames(NDD) = c("E","I","II","III","IV","V")

for (i in 1:n) {
  NDD[i,1] =  sum(M1[i,])
}
for (i in 1:n) {
  NDD[i,2] =  sum(M2[i,])
}
for (i in 1:n) {
  NDD[i,3] =  sum(M3[i,])
}
for (i in 1:n) {
  NDD[i,4] =  sum(M4[i,])
}
for (i in 1:n) {
  NDD[i,5] =  sum(M5[i,])
}
for (i in 1:n) {
  NDD[i,6] =  sum(M6[i,])
}
NDD


#proportion of total population for each stage
PropNDD = matrix(nrow=n,ncol=6); colnames(PropNDD) = c("E","I","II","III","IV","V")
for (i in 1:n){
  PropNDD[i,] = NDD[i,]/(sum(NDD[i,]))
}
PropNDD


#extract day 12 population
D12 = NDD[12,]
P12 = sum(D12) #479,957
TE3 #2,379,294
P12/TE3 #20.17%


#read in populations by landcover type
Eggs = read.csv("EggsLaid.csv")
Eggs = Eggs[-c(18:25),]

#Scenario 3 population by landcover type
MPop3 = data.frame(matrix(nrow = 18, ncol = 7))
colnames(MPop3) = c("Landcover","E","I","II","III","IV","V")
MPop3$Landcover = Eggs$Buffer
#day 12 proportions of starting total for each stage
Prop12 = NDD[12,]/(sum(B))
for (i in 1:18){
    for (j in 2:7){
      MPop3[i,j] = Eggs[i,4]*Prop12[j-1]
    }
}
write.csv(MPop3, file = "MPop3.csv")





############################   2.  Calculate size of population in drift zone  #############################



###### Drift on only N side of soybean fields - first test ############

MedAugDBF = read.csv("MedAugUnionNBuff4_Export.txt")

#empty vector for area
Area3 = c()
#new column for proportion
MedAugDBF$Proportion = NA
#new column for sum of areas for each ORIGINAL poly, as a check
MedAugDBF$OrigArea = NA

#FID_StoryC is FID of original poly

for(i in 1:52952) {
  
  #FID for original poly
  FID = MedAugDBF[i,2]
  #reset Area vector
  Area3 = c()
  #reset index for area vector
  w=1
  
  #get areas into a vector
  for (j in 1:52952){
    if (MedAugDBF[j,2] == FID){
      Area3[w] = MedAugDBF[j,18]
      w = w+1
    }
  }
  
  #put sum of area3 into df
  MedAugDBF[i,20] = sum(Area3)
  
  #find correct Area for the numerator and calc the proportion
  for (m in 1:length(Area3)){
    if (MedAugDBF[i,18]==Area3[m]){
      MedAugDBF[i,19] = Area3[m]/sum(Area3)
    }
    
  }
  
}

#calculate eggs in each new poly - eggs in original poly times proportional area
MedAugDBF$PropEggs = MedAugDBF$Proportion*MedAugDBF$Scen3Egg_2
sum(MedAugDBF$PropEggs) #check, should be 2,379,294 and it is

#eggs inside drift buffer
IN = sum(MedAugDBF[which(MedAugDBF[,17]=="Drift"),21]) #128,399.8
#eggs outside drift buffer
OUT = sum(MedAugDBF[which(MedAugDBF[,17]==" "),21]) #2,250,894
IN+OUT #2,379,294
IN/(IN+OUT) #5.4%
OUT/(IN+OUT) #94.6%

#eggs in the drift zone by habitat type
DZEggs3 = data.frame(matrix(nrow = 18, ncol = 2))
colnames(DZEggs3) = c("Landcover","Eggs")
DZEggs3$Landcover = levels(MedAugDBF[,3])

ind = levels(MedAugDBF[,3])
for (i in 1:18){
  DZEggs3[i,2] = sum(MedAugDBF[which(MedAugDBF[,3]==ind[i] & MedAugDBF[,17]=="Drift"),21])
}
sum(DZEggs3$Eggs) #128,399.8 is correct
DZEggs3$Prop = DZEggs3$Eggs/(sum(DZEggs3$Eggs))
sum(DZEggs3$Prop)
write.csv(DZEggs3, file = "DZEggs3.csv")


####### Drift on NW side of soybean fields Scenario 3 #########

MedAugNWDBF = read.csv("MedAugUnionNWBuff_Export.txt")

#empty vector for area
Area4 = c()
#new column for proportion
MedAugNWDBF$Proportion = NA
#new column for sum of areas for each ORIGINAL poly, as a check
MedAugNWDBF$OrigArea = NA

#FID_StoryC is FID of original poly

time = system.time(
for(i in 1:58751) {
  
  #FID for original poly
  FID = MedAugNWDBF[i,2]
  #reset Area vector
  Area4 = c()
  #reset index for area vector
  w=1
  
  #get areas into a vector
  for (j in 1:58751){
    if (MedAugNWDBF[j,2] == FID){
      Area4[w] = MedAugNWDBF[j,17]
      w = w+1
    }
  }
  
  #put sum of area3 into df
  MedAugNWDBF[i,20] = sum(Area4)
  
  #find correct Area for the numerator and calc the proportion
  for (m in 1:length(Area4)){
    if (MedAugNWDBF[i,17]==Area4[m]){
      MedAugNWDBF[i,19] = Area4[m]/sum(Area4)
    }
    
  }
  
}

)

#calculate eggs in each new poly - eggs in original poly times proportional area
MedAugNWDBF$PropEggs = MedAugNWDBF$Proportion*MedAugNWDBF$Scen3Egg_2
sum(MedAugNWDBF$PropEggs) #check, should be 2,379,294 and it is

#eggs inside drift buffer
IN = sum(MedAugNWDBF[which(MedAugNWDBF[,18]=="Drift"),21]) #232,002.4
#eggs outside drift buffer
OUT = sum(MedAugNWDBF[which(MedAugNWDBF[,18]==" "),21]) #2,147,292
IN+OUT #2,379,294
IN/(IN+OUT) #9.8%
OUT/(IN+OUT) #90.2%

#eggs in the drift zone by habitat type
NWBeanDZEggs3 = data.frame(matrix(nrow = 18, ncol = 2))
colnames(NWBeanDZEggs3) = c("Landcover","Eggs")
NWBeanDZEggs3$Landcover = levels(MedAugNWDBF[,3])

ind = levels(MedAugNWDBF[,3])
for (i in 1:18){
  NWBeanDZEggs3[i,2] = sum(MedAugNWDBF[which(MedAugNWDBF[,3]==ind[i] & MedAugNWDBF[,18]=="Drift"),21])
}
sum(NWBeanDZEggs3$Eggs) #232,002.4 is correct
NWBeanDZEggs3$Prop = NWBeanDZEggs3$Eggs/(sum(NWBeanDZEggs3$Eggs))
sum(NWBeanDZEggs3$Prop)
write.csv(NWBeanDZEggs3, file = "NWBeanDZEggs3.csv")



####### Drift on NW side of soybean fields Scenario 1 #########

BaselineNWDBF = read.csv("BaselineUnionNWBuff_Export.txt")

#empty vector for area
Area4 = c()
#new column for proportion
BaselineNWDBF$Proportion = NA
#new column for sum of areas for each ORIGINAL poly, as a check
BaselineNWDBF$OrigArea = NA

#FID_StoryC is FID of original poly

time = system.time(
  for(i in 1:58751) {
    
    #FID for original poly
    FID = BaselineNWDBF[i,2]
    #reset Area vector
    Area4 = c()
    #reset index for area vector
    w=1
    
    #get areas into a vector
    for (j in 1:58751){
      if (BaselineNWDBF[j,2] == FID){
        Area4[w] = BaselineNWDBF[j,17]
        w = w+1
      }
    }
    
    #put sum of area3 into df
    BaselineNWDBF[i,20] = sum(Area4)
    
    #find correct Area for the numerator and calc the proportion
    for (m in 1:length(Area4)){
      if (BaselineNWDBF[i,17]==Area4[m]){
        BaselineNWDBF[i,19] = Area4[m]/sum(Area4)
      }
      
    }
    
  }
  
)


#calculate eggs in each new poly - eggs in original poly times proportional area
BaselineNWDBF$PropEggs = BaselineNWDBF$Proportion*BaselineNWDBF$Scen1Egg_2
sum(BaselineNWDBF$PropEggs) #check, should be 2,176,354 and it is

#eggs inside drift buffer
IN = sum(BaselineNWDBF[which(BaselineNWDBF[,18]=="Drift"),21]) #186,639.1
#eggs outside drift buffer
OUT = sum(BaselineNWDBF[which(BaselineNWDBF[,18]==" "),21]) #1,989,715
IN+OUT #2,176,354
IN/(IN+OUT) #8.6%
OUT/(IN+OUT) #91.4%

#eggs in the drift zone by habitat type
NWBeanDZEggs1 = data.frame(matrix(nrow = 18, ncol = 2))
colnames(NWBeanDZEggs1) = c("Landcover","Eggs")
NWBeanDZEggs1$Landcover = levels(BaselineNWDBF[,3])

ind = levels(BaselineNWDBF[,3])
for (i in 1:18){
  NWBeanDZEggs1[i,2] = sum(BaselineNWDBF[which(BaselineNWDBF[,3]==ind[i] & BaselineNWDBF[,18]=="Drift"),21])
}
sum(NWBeanDZEggs1$Eggs) #186,639.1 is correct
NWBeanDZEggs1$Prop = NWBeanDZEggs1$Eggs/(sum(NWBeanDZEggs1$Eggs))
sum(NWBeanDZEggs1$Prop)
write.csv(NWBeanDZEggs1, file = "NWBeanDZEggs1.csv")



####### Drift on NW side of soybean fields Scenario 2 #########

MaxAugNWDBF = read.csv("MaxAugUnionNWBuff_Export.txt")

#empty vector for area
Area4 = c()
#new column for proportion
MaxAugNWDBF$Proportion = NA
#new column for sum of areas for each ORIGINAL poly, as a check
MaxAugNWDBF$OrigArea = NA

#FID_StoryC is FID of original poly

time = system.time(
  for(i in 1:58751) {
    
    #FID for original poly
    FID = MaxAugNWDBF[i,2]
    #reset Area vector
    Area4 = c()
    #reset index for area vector
    w=1
    
    #get areas into a vector
    for (j in 1:58751){
      if (MaxAugNWDBF[j,2] == FID){
        Area4[w] = MaxAugNWDBF[j,18]
        w = w+1
      }
    }
    
    #put sum of area3 into df
    MaxAugNWDBF[i,21] = sum(Area4)
    
    #find correct Area for the numerator and calc the proportion
    for (m in 1:length(Area4)){
      if (MaxAugNWDBF[i,18]==Area4[m]){
        MaxAugNWDBF[i,20] = Area4[m]/sum(Area4)
      }
      
    }
    
  }
  
)

#calculate eggs in each new poly - eggs in original poly times proportional area
MaxAugNWDBF$PropEggs = MaxAugNWDBF$Proportion*MaxAugNWDBF$Scen2Egg_1
sum(MaxAugNWDBF$PropEggs) #check, should be 2,713,414 and it is

#eggs inside drift buffer
IN = sum(MaxAugNWDBF[which(MaxAugNWDBF[,19]=="Drift"),22]) #292,971.4
#eggs outside drift buffer
OUT = sum(MaxAugNWDBF[which(MaxAugNWDBF[,19]==" "),22]) #2,420,443
IN+OUT #2,713,414
IN/(IN+OUT) #10.8%
OUT/(IN+OUT) #89.2%

#eggs in the drift zone by habitat type
NWBeanDZEggs2 = data.frame(matrix(nrow = 18, ncol = 2))
colnames(NWBeanDZEggs2) = c("Landcover","Eggs")
NWBeanDZEggs2$Landcover = levels(MaxAugNWDBF[,3])

ind = levels(MaxAugNWDBF[,3])
for (i in 1:18){
  NWBeanDZEggs2[i,2] = sum(MaxAugNWDBF[which(MaxAugNWDBF[,3]==ind[i] & MaxAugNWDBF[,19]=="Drift"),22])
}
sum(NWBeanDZEggs2$Eggs) #292,971.4 is correct
NWBeanDZEggs2$Prop = NWBeanDZEggs2$Eggs/(sum(NWBeanDZEggs2$Eggs))
sum(NWBeanDZEggs2$Prop)
write.csv(NWBeanDZEggs2, file = "NWBeanDZEggs2.csv")



####### Drift on NW side of soybean fields Scenario 4 #########

AugOutNWDBF = read.csv("AugOutsideUnionNWBuff_Export.txt")

#empty vector for area
Area4 = c()
#new column for proportion
AugOutNWDBF$Proportion = NA
#new column for sum of areas for each ORIGINAL poly, as a check
AugOutNWDBF$OrigArea = NA

#FID_StoryC is FID of original poly

time = system.time(
  for(i in 1:58751) {
    
    #FID for original poly
    FID = AugOutNWDBF[i,6]
    #reset Area vector
    Area4 = c()
    #reset index for area vector
    w=1
    
    #get areas into a vector
    for (j in 1:58751){
      if (AugOutNWDBF[j,6] == FID){
        Area4[w] = AugOutNWDBF[j,4]
        w = w+1
      }
    }
    
    #put sum of area3 into df
    AugOutNWDBF[i,20] = sum(Area4)
    
    #find correct Area for the numerator and calc the proportion
    for (m in 1:length(Area4)){
      if (AugOutNWDBF[i,4]==Area4[m]){
        AugOutNWDBF[i,19] = Area4[m]/sum(Area4)
      }
      
    }
    
  }
  
)

#calculate eggs in each new poly - eggs in original poly times proportional area
AugOutNWDBF$PropEggs = AugOutNWDBF$Proportion*AugOutNWDBF$Scen4Egg_2
sum(AugOutNWDBF$PropEggs) #check, should be 2,253,138 and it is

#eggs inside drift buffer
IN = sum(AugOutNWDBF[which(AugOutNWDBF[,5]=="Drift"),21]) #172,365.1
#eggs outside drift buffer
OUT = sum(AugOutNWDBF[which(AugOutNWDBF[,5]==" "),21]) #2,080,773
IN+OUT #2,253,138
IN/(IN+OUT) #7.65%
OUT/(IN+OUT) #92.35%

#eggs in the drift zone by habitat type
NWBeanDZEggs4 = data.frame(matrix(nrow = 18, ncol = 2))
colnames(NWBeanDZEggs4) = c("Landcover","Eggs")
NWBeanDZEggs4$Landcover = levels(AugOutNWDBF[,7])

ind = levels(AugOutNWDBF[,7])
for (i in 1:18){
  NWBeanDZEggs4[i,2] = sum(AugOutNWDBF[which(AugOutNWDBF[,7]==ind[i] & AugOutNWDBF[,5]=="Drift"),21])
}
sum(NWBeanDZEggs4$Eggs) #172,365.1 is correct
NWBeanDZEggs4$Prop = NWBeanDZEggs4$Eggs/(sum(NWBeanDZEggs4$Eggs))
sum(NWBeanDZEggs4$Prop)
write.csv(NWBeanDZEggs4, file = "NWBeanDZEggs4.csv")


#### Drift on NW side of army cornworm fields (random 4% of corn and bean fields) #############

ArmyBaselineNWDBF = read.csv("ArmyBaselineUnionNWBuff_Export2.txt")

#empty vector for area
Area4 = c()
#new column for proportion
ArmyBaselineNWDBF$Proportion = NA
#new column for sum of areas for each ORIGINAL poly, as a check
ArmyBaselineNWDBF$OrigArea = NA

#FID_StoryC is FID of original poly

time = system.time(
  for(i in 1:43124) {
    
    #FID for original poly
    FID = ArmyBaselineNWDBF[i,2]
    #reset Area vector
    Area4 = c()
    #reset index for area vector
    w=1
    
    #get areas into a vector
    for (j in 1:43124){
      if (ArmyBaselineNWDBF[j,2] == FID){
        Area4[w] = ArmyBaselineNWDBF[j,18]
        w = w+1
      }
    }
    
    #put sum of area3 into df
    ArmyBaselineNWDBF[i,20] = sum(Area4)
    
    #find correct Area for the numerator and calc the proportion
    for (m in 1:length(Area4)){
      if (ArmyBaselineNWDBF[i,18]==Area4[m]){
        ArmyBaselineNWDBF[i,19] = Area4[m]/sum(Area4)
      }
      
    }
    
  }
  
)

#calculate eggs in each new poly - eggs in original poly times proportional area
ArmyBaselineNWDBF$PropEggs = ArmyBaselineNWDBF$Proportion*ArmyBaselineNWDBF$Scen1Egg_2
sum(ArmyBaselineNWDBF$PropEggs) #check, should be 2,176,354 and it is

#eggs inside drift buffer
IN = sum(ArmyBaselineNWDBF[which(ArmyBaselineNWDBF[,17]=="Drift"),21]) #15,799
#eggs outside drift buffer
OUT = sum(ArmyBaselineNWDBF[which(ArmyBaselineNWDBF[,17]==" "),21]) #2,160,555
IN+OUT #2,176,354
IN/(IN+OUT) #0.73%
OUT/(IN+OUT) #99.27%

#eggs in the drift zone by habitat type
ArmyNWDZEggs = data.frame(matrix(nrow = 18, ncol = 2))
colnames(ArmyNWDZEggs) = c("Landcover","Eggs")
ArmyNWDZEggs$Landcover = levels(ArmyBaselineNWDBF[,3])

ind = levels(ArmyBaselineNWDBF[,3])
for (i in 1:18){
  ArmyNWDZEggs[i,2] = sum(ArmyBaselineNWDBF[which(ArmyBaselineNWDBF[,3]==ind[i] & ArmyBaselineNWDBF[,17]=="Drift"),21])
}
sum(ArmyNWDZEggs$Eggs) #15,798.77 is correct
ArmyNWDZEggs$Prop = ArmyNWDZEggs$Eggs/(sum(ArmyNWDZEggs$Eggs))
sum(ArmyNWDZEggs$Prop)
write.csv(ArmyNWDZEggs, file = "ArmyNWDZEggs.csv")





##############  4. Add to Surv Est Code to simulate adults produced  ########################

#code below copied from section 1 above then modified

#simulation parameters
#number of days in simulation - need at least 12+13=25 to get everything to adult stage
n=30
#nuisance variable for simulation
pupamax = 16
#daily temperature
MT = rep(21,n) #same temp every day for simulations
#daily survival probabilities
S = out.260.e.sum1$Median
#total eggs laid from scenario 3
TE3 = 2379294
#when do eggs from first day reach mid 5th instar stage? try 14 days first.  So 14 days of eggs
d = 12
#number of eggs laid per day.  Should be of length n
B = rep(TE3/d, d);B[(d+1):n] = 0  
B



#code from SimulateandAnalyzeMonarchCounts.R
ADD = matrix(nrow=n+pupamax, ncol=n); colnames(ADD) = 1:n
for (j in 1:n){
  for (i in (j+1):(n+1)) {
    ADD[i,j] = sum(MT[j:(i-1)])
  }
}
#set upper right off diagonal cells to 0 - set to 0 so following sum operations will work smoothly
for (j in 1:n) {
  for (i in 1:j){
    ADD[i,j] = 0    
  }
}
#set lower cells to 0 as well
for (j in 1:n) {
  for (i in (n+2):(n+pupamax)){
    ADD[i,j] = 0    
  }
}

#calculate stage durations for each cohort. 
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
                                                       ifelse(ADD[j+15,j] > 231.9,16,  ####  Why are there two > 231.9 's here
                                                              17)))))))
}


#########  calculate mean stage durations  ###############
#vector of mean stage durations
MSD = vector(length = 6)
#cohorts that finished 5th instar
I5 = ADD[(n+1),] > 231.9
#stage durations of cohorts that reached pupation
SDC[7,I5]
#fifth instar stage durations for each cohort
SD5 = SDC[7,I5]-SDC[6,I5]
MSD[6] = mean(SD5)
#cohorts that finished 4th instar
I4 = ADD[(n+1),] > 165.3
MSD[5] = mean(SDC[6,I4]-SDC[5,I4])
#cohorts that finished 3rd instar
I3 = ADD[(n+1),] > 129.6
MSD[4] = mean(SDC[5,I3]-SDC[4,I3])
#cohorts that finished 2nd instar
I2 = ADD[(n+1),] > 105.1
MSD[3] = mean(SDC[4,I2]-SDC[3,I2])
#cohorts that finished 1st instar
I1 = ADD[(n+1),] > 77.3
MSD[2] = mean(SDC[3,I1]-SDC[2,I1])
#cohorts that hatched
I0 = ADD[(n+1),] > 45
MSD[1] = mean(SDC[2,I0]-SDC[1,I0])
MSD

#calculate stage survival probabilities and cumulative survival probability
S1 = S[1]^MSD[1] 
S2 = S[2]^MSD[2]
S3 = S[3]^MSD[3]
S4 = S[4]^MSD[4]
S5 = S[5]^MSD[5]
S6 = S[6]^MSD[6]
SC = S1*S2*S3*S4*S5*S6 #cumulative survival probability
SC


#PEDD = Proportion surviving each day for each cohort.  Rows are days, cols are cohorts.  This is p_ik*pi_j-k-bik from Grant et al.
PEDD = matrix(nrow=n+pupamax+1, ncol=n); colnames(PEDD) = 1:n #added 1 to pupamax, seems all it needed was 1 extra row
for(j in 1:n){
  PEDD[j,j] = 1
  for (i in 2:SDC[2,j]){
    PEDD[i+j-1,j] = (S[1])^(i-1)
  }
  for (i in (SDC[2,j]+1):SDC[3,j]){
    PEDD[i+j-1,j] = (S[1])^(SDC[2,j]-1)*(S[2])^(i-SDC[2,j])
  }
  for (i in (SDC[3,j]+1):SDC[4,j]){
    PEDD[i+j-1,j] = (S[1])^(SDC[2,j]-1)*(S[2])^(SDC[3,j]-SDC[2,j])*S[3]^(i-SDC[3,j])
  }
  for (i in (SDC[4,j]+1):SDC[5,j]){
    PEDD[i+j-1,j] = (S[1])^(SDC[2,j]-1)*(S[2])^(SDC[3,j]-SDC[2,j])*S[3]^(SDC[4,j]-SDC[3,j])*S[4]^(i-SDC[4,j])
  }
  for (i in (SDC[5,j]+1):SDC[6,j]){
    PEDD[i+j-1,j] = (S[1])^(SDC[2,j]-1)*(S[2])^(SDC[3,j]-SDC[2,j])*S[3]^(SDC[4,j]-SDC[3,j])*S[4]^(SDC[5,j]-SDC[4,j])*S[5]^(i-SDC[5,j])
  }
  for (i in (SDC[6,j]+1):SDC[7,j]){
    PEDD[i+j-1,j] = (S[1])^(SDC[2,j]-1)*(S[2])^(SDC[3,j]-SDC[2,j])*S[3]^(SDC[4,j]-SDC[3,j])*S[4]^(SDC[5,j]-SDC[4,j])*S[5]^(SDC[6,j]-SDC[5,j])*S[6]^(i-SDC[6,j])
  }
  #pupa who make it to adult
  PEDD[(SDC[7,j]+1)+j-1,j] = (S[1])^(SDC[2,j]-1)*(S[2])^(SDC[3,j]-SDC[2,j])*S[3]^(SDC[4,j]-SDC[3,j])*S[4]^(SDC[5,j]-SDC[4,j])*S[5]^(SDC[6,j]-SDC[5,j])*S[6]^(i-SDC[6,j])*0.76
}
#set NA's to 0
#for (j in 2:n) {
#  for (i in 1:(j-1)) {
#    PEDD[i,j] = 0
#  }
#}
#for (j in 1:(n-1)) {
#  for (i in (SDC[7,j]+1+j):(n+pupamax)) {
#    PEDD[i,j] = 0
#  }
#}


#Number surviving each day of each cohort.  rows are days, cols are cohorts.  These are the cohort contributions that are summed to get l_ij of Grant et al.  
NDEDD = matrix(nrow = n+pupamax+1, ncol = n); colnames(NDEDD) = 1:n
for (j in 1:n) {
  for (i in j:(j+SDC[7,j]-1)){
    NDEDD[i,j] = B[j]*PEDD[i,j]
  }
}
#for (j in 2:n) {
#  for (i in 1:(j-1)) {
#    NDEDD[i,j] = 0
#  }
#}
#for (j in 1:(n-1)) {
#  for (i in (j+SDC[7,j]):(n+pupamax)){
#    NDEDD[i,j] = 0
#  }
#}

#ST - index indicating what stage each cohort is in each day.  rows are days, cols are cohorts.  99 means it has pupated
#this is essentially the indicator function I(s_kj = i) of Grant et al. 
ST = matrix(nrow=n, ncol=n); colnames(ST) = 1:n
for (j in 1:(n-1)) {
  for (i in (j+1):n) {
    ST[i,j] = ifelse(ADD[i,j] >= 1 && ADD[i,j] < 45, 1, 
                     ifelse(ADD[i,j] >= 45 && ADD[i,j] < 77.3, 2,
                            ifelse(ADD[i,j] >= 77.3 && ADD[i,j] < 105.1, 3,
                                   ifelse(ADD[i,j] >= 105.1 && ADD[i,j] < 129.6, 4,
                                          ifelse(ADD[i,j] >= 129.6 && ADD[i,j] < 165.3, 5,
                                                 ifelse(ADD[i,j] >= 165.3 && ADD[i,j] < 231.9, 6, 99))))))
  }
}
for (j in 1:n){
  for (i in j:j){
    ST[i,j] = 1
  }
}
for (j in 2:n){
  for (i in 1:(j-1)){
    ST[i,j] = 0
  }
}

#A matrix for each stage - copy numbers from NDEDD using ST to tell which numbers go into which matrix
#i.e., if a ST[4,5] = 1, then NDEDD[4,5] goes into matrix M1
#this makes it easier/possible to sum the population across cohorts in each stage on any particular day
M1 = matrix(nrow=n, ncol=n); colnames(M1) = 1:n
M2 = matrix(nrow=n, ncol=n); colnames(M2) = 1:n
M3 = matrix(nrow=n, ncol=n); colnames(M3) = 1:n
M4 = matrix(nrow=n, ncol=n); colnames(M4) = 1:n
M5 = matrix(nrow=n, ncol=n); colnames(M5) = 1:n
M6 = matrix(nrow=n, ncol=n); colnames(M6) = 1:n

for (j in 1:n){
  for (i in 1:n) {
    M1[i,j] = ifelse(ST[i,j]==1,NDEDD[i,j],0)
    M2[i,j] = ifelse(ST[i,j]==2,NDEDD[i,j],0)
    M3[i,j] = ifelse(ST[i,j]==3,NDEDD[i,j],0)
    M4[i,j] = ifelse(ST[i,j]==4,NDEDD[i,j],0)
    M5[i,j] = ifelse(ST[i,j]==5,NDEDD[i,j],0)
    M6[i,j] = ifelse(ST[i,j]==6,NDEDD[i,j],0)
  }
}

#NDD = sum populations from each cohort into a matrix of population of each stage each day.  rows are days, cols are stages.  
#This is l_ij of Grant et al. 
NDD = matrix(nrow=n,ncol=6); colnames(NDD) = c("E","I","II","III","IV","V")

for (i in 1:n) {
  NDD[i,1] =  sum(M1[i,])
}
for (i in 1:n) {
  NDD[i,2] =  sum(M2[i,])
}
for (i in 1:n) {
  NDD[i,3] =  sum(M3[i,])
}
for (i in 1:n) {
  NDD[i,4] =  sum(M4[i,])
}
for (i in 1:n) {
  NDD[i,5] =  sum(M5[i,])
}
for (i in 1:n) {
  NDD[i,6] =  sum(M6[i,])
}
NDD


#proportion of total population for each stage
PropNDD = matrix(nrow=n,ncol=6); colnames(PropNDD) = c("E","I","II","III","IV","V")
for (i in 1:n){
  PropNDD[i,] = NDD[i,]/(sum(NDD[i,]))
}
PropNDD


#extract day 12 population
D12 = NDD[12,]
P12 = sum(D12) #479,957
TE3 #2,379,294
P12/TE3 #20.17%


#read in populations by landcover type
Eggs = read.csv("EggsLaid.csv")
Eggs = Eggs[-c(18:25),]

#Scenario 3 population by landcover type
MPop3 = data.frame(matrix(nrow = 18, ncol = 7))
colnames(MPop3) = c("Landcover","E","I","II","III","IV","V")
MPop3$Landcover = Eggs$Buffer
#day 12 proportions of starting total for each stage
Prop12 = NDD[12,]/(sum(B))
for (i in 1:18){
  for (j in 2:7){
    MPop3[i,j] = Eggs[i,4]*Prop12[j-1]
  }
}
write.csv(MPop3, file = "MPop3.csv")







############  Finalized code for combined Natural and Insectide Survival  ####################

#code except for new TISM and PEDD and after copied from Section 1 and 4 above and maybe slightly modified. 


#load survival estimates from StageSurvival.R
out.260.e.sum1 = read.csv("out.260.e.csv")

#simulation parameters
#number of days in simulation
n=30
#nuisance variable for simulation
pupamax = 27
#daily temperature
MT = rep(21,n) #same temp every day for simulations
#daily survival probabilities
S = out.260.e.sum1$Median
#proportion exposed to spray drift (equals 1-proportion under leaf)
#decided to just hard code the numbers because makes it easier to understand
#E = 0.40
#delay between exposure and actual mortality.  In my sims, exposure happens on day 13, change in pop/death happens day 17
#decided not to use the delay, caused too many problems
#del = 4

#total eggs laid for different scenarios, in order that they were analyzed

#Scenario 3 total eggs laid
TE3 = 2379294
#Scenario 3 eggs laid outside soybean drift zone
TE3O = 2147292
#Scenario 3 eggs laid inside soybean drift zone
TE3I = 232002

#Scenario 1 total eggs laid = from section 2 above
TE1 = 2176354
#Scenario 1 eggs laid outside soybean drift zone
TE1O = 1989715
#Scenario 1 eggs laid in soybean drift zone
TE1I = 186639

#Scenario 2 total eggs laid - from section 2 above
TE2 = 2713414
#Scenario 1 eggs laid outside soybean drift zone
TE2O = 2420443
#Scenario 1 eggs laid in soybean drift zone
TE2I = 292971

#Scenario 4 total eggs laid - from section 2 above
TE4 = 2253138
#Scenario 4 eggs laid outside soybean drift zone
TE4O = 2080773
#Scenario 4 eggs laid in soybean drift zone
TE4I = 172365


#number of days for eggs to be laid
d = 12

#number of eggs laid per day.  Should be of length n
#change variable TE* to variable for correct scenario above
B = rep(TE4I/d, d);B[(d+1):n] = 0  
B

####### be sure to define TISM correctly, it is at the bottom though  ###########


#code from SimulateandAnalyzeMonarchCounts.R
ADD = matrix(nrow=n+pupamax, ncol=n); colnames(ADD) = 1:n
for (j in 1:n){
  for (i in (j+1):(n+1)) {
    ADD[i,j] = sum(MT[j:(i-1)])
  }
}
#set upper right off diagonal cells to 0 - set to 0 so following sum operations will work smoothly
for (j in 1:n) {
  for (i in 1:j){
    ADD[i,j] = 0    
  }
}
#set lower cells to 0 as well
for (j in 1:n) {
  for (i in (n+2):(n+pupamax)){
    ADD[i,j] = 0    
  }
}

#calculate stage durations for each cohort. 
SDC = matrix(nrow = 8, ncol = n); colnames(SDC) = 1:n; rownames(SDC) = c("egg","first","second","third","fourth","fifth","pupa","adult")
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
  SDC[8,j] = ifelse(ADD[j+12,j] > 351.8,13,
                    ifelse(ADD[j+13,j] > 351.8,14,
                           ifelse(ADD[j+14,j] > 351.8,15,
                                  ifelse(ADD[j+15,j] > 351.8,16,
                                         ifelse(ADD[j+16,j] > 351.8,17,
                                                ifelse(ADD[j+17,j] > 351.8,18,
                                                       ifelse(ADD[j+18,j] > 351.8,19,
                                                              ifelse(ADD[j+19,j] > 351.8,20,
                                                                     ifelse(ADD[j+21,j] > 351.8,22,
                                                                            ifelse(ADD[j+22,j] > 351.8,23,
                                                                                   ifelse(ADD[j+23,j] > 351.8,24,
                                                                                          ifelse(ADD[j+24,j] > 351.8,25,
                                                                                                 ifelse(ADD[j+25,j] > 351.8,26,
                                                                                                        27)))))))))))))
}



#########  calculate mean stage durations  ###############
#vector of mean stage durations
MSD = vector(length = 7)
#cohorts that finished pupation
I6 = ADD[(n+1),] > 351.8
#stage durations of cohorts that reached adult
SDC[8,I6]
#cohorts that finished 5th instar
I5 = ADD[(n+1),] > 231.9
MSD[7] = mean(SDC[8,I6]-SDC[7,I6])
#stage durations of cohorts that reached pupation
SDC[7,I5]
#fifth instar stage durations for each cohort
SD5 = SDC[7,I5]-SDC[6,I5]
MSD[6] = mean(SD5)
#cohorts that finished 4th instar
I4 = ADD[(n+1),] > 165.3
MSD[5] = mean(SDC[6,I4]-SDC[5,I4])
#cohorts that finished 3rd instar
I3 = ADD[(n+1),] > 129.6
MSD[4] = mean(SDC[5,I3]-SDC[4,I3])
#cohorts that finished 2nd instar
I2 = ADD[(n+1),] > 105.1
MSD[3] = mean(SDC[4,I2]-SDC[3,I2])
#cohorts that finished 1st instar
I1 = ADD[(n+1),] > 77.3
MSD[2] = mean(SDC[3,I1]-SDC[2,I1])
#cohorts that hatched
I0 = ADD[(n+1),] > 45
MSD[1] = mean(SDC[2,I0]-SDC[1,I0])
MSD

#calculate stage survival probabilities and cumulative survival probability
S1 = S[1]^MSD[1] 
S2 = S[2]^MSD[2]
S3 = S[3]^MSD[3]
S4 = S[4]^MSD[4]
S5 = S[5]^MSD[5]
S6 = S[6]^MSD[6]
SC = S1*S2*S3*S4*S5*S6 #cumulative survival probability to pupation
SC
#pupal survival is 0.76 total, daily survival over msd below
S[7] = 0.76^(1/MSD[7])


#moved ST here from between NDEDD and M's, change nrow from n to n+pupamax+1
#ST - index indicating what stage each cohort is in each day.  rows are days, cols are cohorts.  
#modified to include adults - 7 is pupa, 8 is adult
ST = matrix(nrow=n+pupamax+1, ncol=n); colnames(ST) = 1:n
for (j in 1:(n-1)) {
  for (i in (j+1):n) {
    ST[i,j] = ifelse(ADD[i,j] >= 1 && ADD[i,j] < 45, 1, 
                     ifelse(ADD[i,j] >= 45 && ADD[i,j] < 77.3, 2,
                            ifelse(ADD[i,j] >= 77.3 && ADD[i,j] < 105.1, 3,
                                   ifelse(ADD[i,j] >= 105.1 && ADD[i,j] < 129.6, 4,
                                          ifelse(ADD[i,j] >= 129.6 && ADD[i,j] < 165.3, 5,
                                                 ifelse(ADD[i,j] >= 165.3 && ADD[i,j] < 231.9, 6, 
                                                        ifelse(ADD[i,j] >= 231.9 && ADD[i,j] < 351.8, 7, 8)))))))
  }
}
for (j in 1:n){
  for (i in j:j){
    ST[i,j] = 1
  }
}
#set upper right off diagonal cells to 1 - will use egg stage for days before eggs laid - don't need now
#for (j in 1:n) {
#  for (i in 1:j){
#    ST[i,j] = 1    
#  }
#}
#add buffer of size del rows to ST
#STB = matrix(1, nrow = del, ncol = n)
#ST = rbind(STB,ST)



##############  Insecticide Mortality/Survival  #############

#matrix of insecticide survival.  Put the survival in the day it happens, not the next day.  The survival is used to calculate the pop
#for the next day but the code takes care of it.  To repeat, put the survival in the day the drift/exposure occurs
#number should be a survival probality ranging from 0 to 1.  
#TISM = matrix(1,nrow = n+pupamax+1, ncol = 7); rownames(TISM) = 1:(n+pupamax+1); colnames(TISM) = c("egg","first","second","third","fourth","fifth","pupa")

#scenarios for TISM are below, after NDD


#TISM buffer for delay in mortality - don't need now, 2nd v. of PEDD doesn't work
#TISMB = matrix(1, nrow = del, ncol = 7)
#TISM = rbind(TISMB, TISM); rownames(TISM) = 1:(n+pupamax+1+del)


#NEW version of PEDD
#PEDD = Proportion surviving each day for each cohort.  Rows are days, cols are cohorts.  This is p_ik*pi_j-k-bik from Grant et al.
PEDD = matrix(nrow=n+pupamax+1, ncol=n); colnames(PEDD) = 1:n #added 1 to pupamax, seems all it needed was 1 extra row
for(j in 1:n){
  PEDD[j,j] = 1
  for (i in 2:SDC[2,j]){
    PEDD[i+j-1,j] = S[1]^(i-1)*(prod(TISM[j:(i+j-2),1]))
  }
  for (i in (SDC[2,j]+1):SDC[3,j]){   
    PEDD[i+j-1,j] = S[1]^(SDC[2,j]-1)*(prod(TISM[j:(j+SDC[2,j]-2),1]))*S[2]^(i-SDC[2,j])*(prod(TISM[(j+SDC[2,j]-1):(i+j-2),2]))
  }
  for (i in (SDC[3,j]+1):SDC[4,j]){
    PEDD[i+j-1,j] = S[1]^(SDC[2,j]-1)*(prod(TISM[j:(j+SDC[2,j]-2),1]))*S[2]^(SDC[3,j]-SDC[2,j])*(prod(TISM[(j+SDC[2,j]-1):(j+SDC[3,j]-2),2]))*S[3]^(i-SDC[3,j])*(prod(TISM[(j+SDC[3,j]-1):(i+j-2),3]))
  }
  
  #so for j=10, day 17, i = 8, TISM should be TISM[17,X], where X is stage at day 13, use ST and move it up
  #test - PEDD[17,10] = ...*S[4]^(i-SDC[4,j])*(prod(TISM[(j+SDC[4,j]-1):(i+j-2),4])) = ...*S[4]^1*TISM[(j+SDC[4,j]):(i+j-1),ST[i+j-1,j]].
  
  for (i in (SDC[4,j]+1):SDC[5,j]){
    PEDD[i+j-1,j] = S[1]^(SDC[2,j]-1)*(prod(TISM[j:(j+SDC[2,j]-2),1]))*S[2]^(SDC[3,j]-SDC[2,j])*(prod(TISM[(j+SDC[2,j]-1):(j+SDC[3,j]-2),2]))*S[3]^(SDC[4,j]-SDC[3,j])*(prod(TISM[(j+SDC[3,j]-1):(j+SDC[4,j]-2),3]))*S[4]^(i-SDC[4,j])*(prod(TISM[(j+SDC[4,j]-1):(i+j-2),4]))
  }
  for (i in (SDC[5,j]+1):SDC[6,j]){
    PEDD[i+j-1,j] = S[1]^(SDC[2,j]-1)*(prod(TISM[j:(j+SDC[2,j]-2),1]))*S[2]^(SDC[3,j]-SDC[2,j])*(prod(TISM[(j+SDC[2,j]-1):(j+SDC[3,j]-2),2]))*S[3]^(SDC[4,j]-SDC[3,j])*(prod(TISM[(j+SDC[3,j]-1):(j+SDC[4,j]-2),3]))*S[4]^(SDC[5,j]-SDC[4,j])*(prod(TISM[(j+SDC[4,j]-1):(j+SDC[5,j]-2),4]))*S[5]^(i-SDC[5,j])*(prod(TISM[(j+SDC[5,j]-1):(i+j-2),5]))
  }
  for (i in (SDC[6,j]+1):SDC[7,j]){
    PEDD[i+j-1,j] = S[1]^(SDC[2,j]-1)*(prod(TISM[j:(j+SDC[2,j]-2),1]))*S[2]^(SDC[3,j]-SDC[2,j])*(prod(TISM[(j+SDC[2,j]-1):(j+SDC[3,j]-2),2]))*S[3]^(SDC[4,j]-SDC[3,j])*(prod(TISM[(j+SDC[3,j]-1):(j+SDC[4,j]-2),3]))*S[4]^(SDC[5,j]-SDC[4,j])*(prod(TISM[(j+SDC[4,j]-1):(j+SDC[5,j]-2),4]))*S[5]^(SDC[6,j]-SDC[5,j])*(prod(TISM[(j+SDC[5,j]-1):(j+SDC[6,j]-2),5]))*S[6]^(i-SDC[6,j])*(prod(TISM[(j+SDC[6,j]-1):(i+j-2),6]))
  }
  #for pupal survival, the overall rate of 0.76 is transformed into a daily rate for however long the pupal stage is.  daily rate calculated above under MSD
  for (i in (SDC[7,j]+1):SDC[8,j]){
  PEDD[i+j-1,j] = S[1]^(SDC[2,j]-1)*(prod(TISM[j:(j+SDC[2,j]-2),1]))*S[2]^(SDC[3,j]-SDC[2,j])*(prod(TISM[(j+SDC[2,j]-1):(j+SDC[3,j]-2),2]))*S[3]^(SDC[4,j]-SDC[3,j])*(prod(TISM[(j+SDC[3,j]-1):(j+SDC[4,j]-2),3]))*S[4]^(SDC[5,j]-SDC[4,j])*(prod(TISM[(j+SDC[4,j]-1):(j+SDC[5,j]-2),4]))*S[5]^(SDC[6,j]-SDC[5,j])*(prod(TISM[(j+SDC[5,j]-1):(j+SDC[6,j]-2),5]))*S[6]^(SDC[7,j]-SDC[6,j])*(prod(TISM[(j+SDC[6,j]-1):(j+SDC[7,j]-2),6]))*S[7]^(i-SDC[7,j])*(prod(TISM[(j+SDC[7,j]-1):(i+j-2),7]))
  }
}


#Number surviving each day of each cohort.  rows are days, cols are cohorts.  These are the cohort contributions that are summed to get l_ij of Grant et al.  
#modified here to add adults to the end - just had to remove -1 from j+SDC[7,j]-1
NDEDD = matrix(nrow = n+pupamax+1, ncol = n); colnames(NDEDD) = 1:n
for (j in 1:n) {
  for (i in j:(j+SDC[8,j])){
    NDEDD[i,j] = B[j]*PEDD[i,j]
  }
}


#A matrix for each stage - copy numbers from NDEDD using ST to tell which numbers go into which matrix
#i.e., if a ST[4,5] = 1, then NDEDD[4,5] goes into matrix M1
#this makes it easier/possible to sum the population across cohorts in each stage on any particular day
M1 = matrix(nrow=n, ncol=n); colnames(M1) = 1:n
M2 = matrix(nrow=n, ncol=n); colnames(M2) = 1:n
M3 = matrix(nrow=n, ncol=n); colnames(M3) = 1:n
M4 = matrix(nrow=n, ncol=n); colnames(M4) = 1:n
M5 = matrix(nrow=n, ncol=n); colnames(M5) = 1:n
M6 = matrix(nrow=n, ncol=n); colnames(M6) = 1:n
M7 = matrix(nrow=n, ncol=n); colnames(M7) = 1:n
M8 = matrix(nrow=n, ncol=n); colnames(M8) = 1:n

for (j in 1:n){
  for (i in 1:n) {
    M1[i,j] = ifelse(ST[i,j]==1,NDEDD[i,j],0)
    M2[i,j] = ifelse(ST[i,j]==2,NDEDD[i,j],0)
    M3[i,j] = ifelse(ST[i,j]==3,NDEDD[i,j],0)
    M4[i,j] = ifelse(ST[i,j]==4,NDEDD[i,j],0)
    M5[i,j] = ifelse(ST[i,j]==5,NDEDD[i,j],0)
    M6[i,j] = ifelse(ST[i,j]==6,NDEDD[i,j],0)
    M7[i,j] = ifelse(ST[i,j]==7,NDEDD[i,j],0)
    M8[i,j] = ifelse(ST[i,j]==8,NDEDD[i,j],0)
  }
}

#NDD = sum populations from each cohort into a matrix of population of each stage each day.  rows are days, cols are stages.
#note that for landscape tox sims that carry out to adults, adults are not cumulative.  Other stages are cumulative each day, but 
#for adults, its the number of adults produced that day.  So total adults produced is the sum of the adults column.  
#This is l_ij of Grant et al. 
NDD = matrix(nrow=n,ncol=8); colnames(NDD) = c("E","I","II","III","IV","V","P","A")

for (i in 1:n) {
  NDD[i,1] =  sum(M1[i,], na.rm = TRUE)
}
for (i in 1:n) {
  NDD[i,2] =  sum(M2[i,], na.rm = TRUE)
}
for (i in 1:n) {
  NDD[i,3] =  sum(M3[i,], na.rm = TRUE)
}
for (i in 1:n) {
  NDD[i,4] =  sum(M4[i,], na.rm = TRUE)
}
for (i in 1:n) {
  NDD[i,5] =  sum(M5[i,], na.rm = TRUE)
}
for (i in 1:n) {
  NDD[i,6] =  sum(M6[i,], na.rm = TRUE)
}
for (i in 1:n) {
  NDD[i,7] =  sum(M7[i,], na.rm = TRUE)
}
for (i in 1:n) {
  NDD[i,8] =  sum(M8[i,], na.rm = TRUE)
}
NDD

#total adults produced
sum(NDD[,8])


plot(NDD[,3])
lines(NDD[,1])
lines(NDD[,2])
lines(NDD[,3])
lines(NDD[,4])
lines(NDD[,5])
lines(NDD[,6])
lines(NDD[,7])
lines(NDD[,8])





######## RESULTS - Soybean Aphid with Aerial Drift  ##################

#pesticide exposure scenarios
#create TISM and then run code above from PEDD down 
#cuticular exposure happens on day 13, survival effect happens on day 14
#oral exposure survival happens on day 13 also, mortality effect calculated on day 17

#for cuticular survival, used first for second and 3rd for fourth
#for dietary survival, used second for first and 3rd for fourth
#for now, 100% insecticide survival for eggs and pupa


####### SCENARIO 3 ########################

#Scenario 3 Soybean Aphid Monarch Production - 2379294 total eggs laid, 232,002.4 outside, 2,147,292 inside
SC3SoyProd = data.frame(matrix(nrow = 3, ncol = 6))
colnames(SC3SoyProd) = c("No Drift","BCF","CTR","CFS","IMI","TMX")
rownames(SC3SoyProd) = c("Outside Drift Zone Prod","Inside Drift Zone Prod","Total Production")

#production
SC3SoyProd[3,1] = AdsProd[1] #I orginally used AdsProd but then switched to SC3SoyProd, so put AdsProd into SC3SoyProd
#production inside drift zone - just subtraction - outside calc'd below
SC3SoyProd[2,1] = SC3SoyProd[3,1]-SC3SoyProd[1,1]

#Baseline - no insecticide, TISM all = 1
TISM = matrix(1,nrow = n+pupamax+1, ncol = 7); rownames(TISM) = 1:(n+pupamax+1); colnames(TISM) = c("egg","first","second","third","fourth","fifth","pupa")
TISM[13,]=c(1,1,1,1,1,1,1)
sum(NDD[,8]) #38,038.46 when TE3 = 2379294
#vector for final adults produced in each scenario
AdsProd = c()
AdsProd[1] = sum(NDD[,8])


#production out side soybean drift zone
SC3SoyProd[1,1:6] = sum(NDD[,8])


######### TISM ##############

#cuticular and dietary exposure
#numbers in parentheses are cuticular exposure times dietary exposure
#eggs have no cuticular exposure, but we decided to use dietary exposure for 1st instar for eggs
#pupal survival is 100% unless hits spiracle, which maybe never happens

#the pattern is (0.40*S.CU + 0.60)*S.DE
#0.40 is the proportion of individuals exposed to spray drift
#where S.CU is cuticular survival rate and S.DE is dietary exposure survival rate

#BCF
TISM[13,] = c(((0.40*1+0.60)*0.13),((0.40*0+0.60)*0.13),((0.40*0+0.60)*0.13),((0.40*0.01+0.60)*0.44),((0.40*0.01+0.60)*0.44),((0.40*0.06+0.60)*0.39),1)
SC3SoyProd[2,2] = sum(NDD[,8])
SC3SoyProd[3,2] = SC3SoyProd[1,2]+SC3SoyProd[2,2]

#CTR
TISM[13,]=c(((0.40*1+0.60)*0.01),((0.40*0+0.60)*0.01),((0.40*0+0.60)*0.01),((0.40*0.06+0.60)*0.12),((0.40*0.06+0.60)*0.12),((0.40*0.09+0.60)*0.35),1)
SC3SoyProd[2,3] = sum(NDD[,8])
SC3SoyProd[3,3] = SC3SoyProd[1,3]+SC3SoyProd[2,3]

#CFS
TISM[13,]=c(((0.40*1+0.60)*0.13),((0.40*0.37+0.60)*0.13),((0.40*0.37+0.60)*0.13),((0.40*0.56+0.60)*0.15),((0.40*0.56+0.60)*0.15),((0.40*0.82+0.60)*0.19),1)
SC3SoyProd[2,4] = sum(NDD[,8])
SC3SoyProd[3,4] = SC3SoyProd[1,4]+SC3SoyProd[2,4]

#IMI
TISM[13,]=c(((0.40*1+0.60)*0.75),((0.40*0.77+0.60)*0.75),((0.40*0.77+0.60)*0.75),((0.40*0.89+0.60)*0.93),((0.40*0.89+0.60)*0.93),((0.40*0.98+0.60)*0.70),1)
SC3SoyProd[2,5] = sum(NDD[,8])
SC3SoyProd[3,5] = SC3SoyProd[1,5]+SC3SoyProd[2,5]

#TMX
TISM[13,]=c(((0.40*1+0.60)*0.73),((0.40*0.78+0.60)*0.73),((0.40*0.78+0.60)*0.73),((0.40*0.96+0.60)*0.81),((0.40*0.96+0.60)*0.81),((0.40*1+0.60)*0.96),1)
SC3SoyProd[2,6] = sum(NDD[,8])
SC3SoyProd[3,6] = SC3SoyProd[1,6]+SC3SoyProd[2,6]

#CDN
TISM[13,]=c(((0.40*1+0.60)*0.47),((0.40*0.01+0.60)*0.47),((0.40*0.01+0.60)*0.47),((0.40*0.19+0.60)*0.93),((0.40*0.19+0.60)*0.93),((0.40*0.78+0.60)*0.27),1)
SC3SoyProdCDN = c()
SC3SoyProdCDN[1] = SC3SoyProd[1,6]
SC3SoyProdCDN[2] = sum(NDD[,8])
SC3SoyProdCDN[3] = SC3SoyProdCDN[1]+SC3SoyProdCDN[2]
SC3SoyProd$CDN = SC3SoyProdCDN

write.csv(SC3SoyProd, file = "SC3SoyProd.csv")

#bar graphs
library(ggplot2)
library(scales)
SC3SoyProdBar = t(SC3SoyProd)
DriftCol = c("NoDrift","BCF","CTR","CFS","IMI","TMX")
SC3SoyProdBar = cbind(SC3SoyProdBar,DriftCol)
typeof(SC3SoyProdBar)
SC3SoyProdBar = data.frame(SC3SoyProdBar)
#doesn't work, type is all wrong

SC3SoyProdBar = data.frame(matrix(nrow = 7, ncol = 4))
colnames(SC3SoyProdBar) = c("Insecticide","Out","In","Total")
SC3SoyProdBar$Insecticide = c("No Drift","BCF","CTR","CFS","IMI","TMX","CDN")
SC3SoyProdBar$In = c(SC3SoyProd[2,])
SC3SoyProdBar$Out = c(SC3SoyProd[1,])
SC3SoyProdBar$Total = c(SC3SoyProd[3,])
SC3SoyProdBar$Insecticide = factor(SC3SoyProdBar$Insecticide, levels = c("No Drift","BCF","CTR","CFS","IMI","TMX","CDN"))
#for some reason SC3SoyProdBar is all lists, so data labels don't work
SC3SoyProdBar$Total = round(unlist(SC3SoyProdBar$Total))
SC3SoyProdBar$In = round(unlist(SC3SoyProdBar$In))


p = ggplot(data = SC3SoyProdBar, aes(x = Insecticide, y = Total)) + geom_bar(stat="identity") + ylab("Adult Monarch Production") +
  geom_text(aes(label = comma(Total)), vjust = 1.5, color = "white")
p

p2 = ggplot(data = SC3SoyProdBar, aes(x = Insecticide, y = In)) + geom_bar(stat="identity") + ylab("Adult Monarch Production") +
  geom_text(aes(label = comma(In)), vjust = 1.5, color = "white")
p2



####### SCENARIO 1 ########################

#Scenario 1 Soybean Aphid Monarch Production
SC1SoyProd = data.frame(matrix(nrow = 3, ncol = 6))
colnames(SC1SoyProd) = c("No Drift","BCF","CTR","CFS","IMI","TMX")
rownames(SC1SoyProd) = c("Outside Drift Zone Prod","Inside Drift Zone Prod","Total Production")

#Baseline - no insecticide, TISM all = 1
TISM[13,]=c(1,1,1,1,1,1,1)
#total eggs
SC1SoyProd[3,1] = sum(NDD[,8]) 
#eggs laid out side drift zone
SC1SoyProd[1,1:6] = sum(NDD[,8])
#eggs laid inside drift zone
SC1SoyProd[2,1] = SC1SoyProd[3,1]-SC1SoyProd[1,1]

#BCF - used TISM from above
SC1SoyProd[2,2] = sum(NDD[,8]) 
SC1SoyProd[3,2] = SC1SoyProd[1,2]+SC1SoyProd[2,2]

#CTR
SC1SoyProd[2,3] = sum(NDD[,8]) 
SC1SoyProd[3,3] = SC1SoyProd[1,3]+SC1SoyProd[2,3]

#CFS
SC1SoyProd[2,4] = sum(NDD[,8]) 
SC1SoyProd[3,4] = SC1SoyProd[1,4]+SC1SoyProd[2,4]

#IMI
SC1SoyProd[2,5] = sum(NDD[,8]) 
SC1SoyProd[3,5] = SC1SoyProd[1,5]+SC1SoyProd[2,5]

#TMX
SC1SoyProd[2,6] = sum(NDD[,8]) 
SC1SoyProd[3,6] = SC1SoyProd[1,6]+SC1SoyProd[2,6]

#CDN
SC1SoyProdCDN = c()
SC1SoyProdCDN[1] = SC1SoyProd[1,6]
SC1SoyProdCDN[2] = sum(NDD[,8])
SC1SoyProdCDN[3] = SC1SoyProdCDN[1]+SC1SoyProdCDN[2]
SC1SoyProd$CDN = SC1SoyProdCDN


#figures
SC1SoyProdBar = data.frame(matrix(nrow = 7, ncol = 4))
colnames(SC1SoyProdBar) = c("Insecticide","Out","In","Total")
SC1SoyProdBar$Insecticide = c("No Drift","BCF","CTR","CFS","IMI","TMX","CDN")
SC1SoyProdBar$In = c(SC1SoyProd[2,])
SC1SoyProdBar$Out = c(SC1SoyProd[1,])
SC1SoyProdBar$Total = c(SC1SoyProd[3,])
SC1SoyProdBar$Insecticide = factor(SC1SoyProdBar$Insecticide, levels = c("No Drift","BCF","CTR","CFS","IMI","TMX","CDN"))
#for some reason SC1SoyProdBar is all lists, so data labels don't work
SC1SoyProdBar$Total = round(unlist(SC1SoyProdBar$Total))
SC1SoyProdBar$In = round(unlist(SC1SoyProdBar$In))

p3 = ggplot(data = SC1SoyProdBar, aes(x = Insecticide, y = Total)) + geom_bar(stat="identity") + ylab("Adult Monarch Production") +
  geom_text(aes(label = comma(Total)), vjust = 1.5, color = "white")
p3

p4 = ggplot(data = SC1SoyProdBar, aes(x = Insecticide, y = In)) + geom_bar(stat="identity") + ylab("Adult Monarch Production") +
  geom_text(aes(label = comma(In)), vjust = 1.5, color = "white")
p4



####### SCENARIO 2 ########################

#Scenario 2 Soybean Aphid Monarch Production
SC2SoyProd = data.frame(matrix(nrow = 3, ncol = 6))
colnames(SC2SoyProd) = c("No Drift","BCF","CTR","CFS","IMI","TMX")
rownames(SC2SoyProd) = c("Outside Drift Zone Prod","Inside Drift Zone Prod","Total Production")

#Baseline - no insecticide, TISM all = 1
TISM[13,]=c(1,1,1,1,1,1,1)
#total eggs
SC2SoyProd[3,1] = sum(NDD[,8]) 
#eggs laid out side drift zone
SC2SoyProd[1,1:6] = sum(NDD[,8])
#eggs laid inside drift zone
SC2SoyProd[2,1] = SC2SoyProd[3,1]-SC2SoyProd[1,1]

#BCF - used TISM from above
SC2SoyProd[2,2] = sum(NDD[,8]) 
SC2SoyProd[3,2] = SC2SoyProd[1,2]+SC2SoyProd[2,2]

#CTR
SC2SoyProd[2,3] = sum(NDD[,8]) 
SC2SoyProd[3,3] = SC2SoyProd[1,3]+SC2SoyProd[2,3]

#CFS
SC2SoyProd[2,4] = sum(NDD[,8]) 
SC2SoyProd[3,4] = SC2SoyProd[1,4]+SC2SoyProd[2,4]

#IMI
SC2SoyProd[2,5] = sum(NDD[,8]) 
SC2SoyProd[3,5] = SC2SoyProd[1,5]+SC2SoyProd[2,5]

#TMX
SC2SoyProd[2,6] = sum(NDD[,8]) 
SC2SoyProd[3,6] = SC2SoyProd[1,6]+SC2SoyProd[2,6]

#CDN
SC2SoyProdCDN = c()
SC2SoyProdCDN[1] = SC2SoyProd[1,6]
SC2SoyProdCDN[2] = sum(NDD[,8])
SC2SoyProdCDN[3] = SC2SoyProdCDN[1]+SC2SoyProdCDN[2]
SC2SoyProd$CDN = SC2SoyProdCDN


#figures
SC2SoyProdBar = data.frame(matrix(nrow = 7, ncol = 4))
colnames(SC2SoyProdBar) = c("Insecticide","Out","In","Total")
SC2SoyProdBar$Insecticide = c("No Drift","BCF","CTR","CFS","IMI","TMX","CDN")
SC2SoyProdBar$In = c(SC2SoyProd[2,])
SC2SoyProdBar$Out = c(SC2SoyProd[1,])
SC2SoyProdBar$Total = c(SC2SoyProd[3,])
SC2SoyProdBar$Insecticide = factor(SC2SoyProdBar$Insecticide, levels = c("No Drift","BCF","CTR","CFS","IMI","TMX","CDN"))
#for some reason SC2SoyProdBar is all lists, so data labels don't work
SC2SoyProdBar$Total = round(unlist(SC2SoyProdBar$Total))
SC2SoyProdBar$In = round(unlist(SC2SoyProdBar$In))

p5 = ggplot(data = SC2SoyProdBar, aes(x = Insecticide, y = Total)) + geom_bar(stat="identity") + ylab("Adult Monarch Production") +
  geom_text(aes(label = comma(Total)), vjust = 1.5, color = "white")
p5

p6 = ggplot(data = SC2SoyProdBar, aes(x = Insecticide, y = In)) + geom_bar(stat="identity") + ylab("Adult Monarch Production") +
  geom_text(aes(label = comma(In)), vjust = 1.5, color = "white")
p6



####### SCENARIO 4 ########################

#Scenario 4 Soybean Aphid Monarch Production
SC4SoyProd = data.frame(matrix(nrow = 3, ncol = 6))
colnames(SC4SoyProd) = c("No Drift","BCF","CTR","CFS","IMI","TMX")
rownames(SC4SoyProd) = c("Outside Drift Zone Prod","Inside Drift Zone Prod","Total Production")

#Baseline - no insecticide, TISM all = 1
TISM[13,]=c(1,1,1,1,1,1,1)
#total eggs
SC4SoyProd[3,1] = sum(NDD[,8]) 
#eggs laid out side drift zone
SC4SoyProd[1,1:6] = sum(NDD[,8])
#eggs laid inside drift zone
SC4SoyProd[2,1] = SC4SoyProd[3,1]-SC4SoyProd[1,1]

#BCF - used TISM from above
SC4SoyProd[2,2] = sum(NDD[,8]) 
SC4SoyProd[3,2] = SC4SoyProd[1,2]+SC4SoyProd[2,2]

#CTR
SC4SoyProd[2,3] = sum(NDD[,8]) 
SC4SoyProd[3,3] = SC4SoyProd[1,3]+SC4SoyProd[2,3]

#CFS
SC4SoyProd[2,4] = sum(NDD[,8]) 
SC4SoyProd[3,4] = SC4SoyProd[1,4]+SC4SoyProd[2,4]

#IMI
SC4SoyProd[2,5] = sum(NDD[,8]) 
SC4SoyProd[3,5] = SC4SoyProd[1,5]+SC4SoyProd[2,5]

#TMX
SC4SoyProd[2,6] = sum(NDD[,8]) 
SC4SoyProd[3,6] = SC4SoyProd[1,6]+SC4SoyProd[2,6]

#CDN
SC4SoyProdCDN = c()
SC4SoyProdCDN[1] = SC4SoyProd[1,6]
SC4SoyProdCDN[2] = sum(NDD[,8])
SC4SoyProdCDN[3] = SC4SoyProdCDN[1]+SC4SoyProdCDN[2]
SC4SoyProd$CDN = SC4SoyProdCDN


#figures
SC4SoyProdBar = data.frame(matrix(nrow = 7, ncol = 4))
colnames(SC4SoyProdBar) = c("Insecticide","Out","In","Total")
SC4SoyProdBar$Insecticide = c("No Drift","BCF","CTR","CFS","IMI","TMX","CDN")
SC4SoyProdBar$In = c(SC4SoyProd[2,])
SC4SoyProdBar$Out = c(SC4SoyProd[1,])
SC4SoyProdBar$Total = c(SC4SoyProd[3,])
SC4SoyProdBar$Insecticide = factor(SC4SoyProdBar$Insecticide, levels = c("No Drift","BCF","CTR","CFS","IMI","TMX","CDN"))
#for some reason SC4SoyProdBar is all lists, so data labels don't work
SC4SoyProdBar$Total = round(unlist(SC4SoyProdBar$Total))
SC4SoyProdBar$In = round(unlist(SC4SoyProdBar$In))

p7 = ggplot(data = SC4SoyProdBar, aes(x = Insecticide, y = Total)) + geom_bar(stat="identity") + ylab("Adult Monarch Production") +
  geom_text(aes(label = comma(Total)), vjust = 1.5, color = "white")
p7

p8 = ggplot(data = SC4SoyProdBar, aes(x = Insecticide, y = In)) + geom_bar(stat="identity") + ylab("Adult Monarch Production") +
  geom_text(aes(label = comma(In)), vjust = 1.5, color = "white")
p8


#figure for no drift monarch production for scenarios 1-4
SC1.4SoyProdBar = data.frame(matrix(nrow = 4, ncol = 2))
colnames(SC1.4SoyProdBar) = c("Scenario","Total")
SC1.4SoyProdBar$Scenario = c("Scenario 1","Scenario 2","Scenario 3","Scenario 4")
SC1.4SoyProdBar$Total[1] = SC1SoyProdBar[1,4]
SC1.4SoyProdBar$Total[2] = SC2SoyProdBar[1,4]
SC1.4SoyProdBar$Total[3] = SC3SoyProdBar[1,4]
SC1.4SoyProdBar$Total[4] = SC4SoyProdBar[1,4]
SC1.4SoyProdBar$Total = round(unlist(SC1.4SoyProdBar$Total))

p9 = ggplot(data = SC1.4SoyProdBar, aes(x = Scenario, y = Total)) + geom_bar(stat="identity") + ylab("Adult Monarch Production") +
  geom_text(aes(label = comma(Total)), vjust = 1.5, color = "white")
p9




#3D graph for scenarios, pesticides and production

library(plot3D)
data("VADeaths")

SC1SoyProdBar2 = SC1SoyProdBar
SC1SoyProdBar2$Scenario = rep("Scenario 1",7)
SC1SoyProdBar2 = SC1SoyProdBar2[, c(1,3,5)]
SC2SoyProdBar2 = SC2SoyProdBar
SC2SoyProdBar2$Scenario = rep("Scenario 2",7)
SC2SoyProdBar2 = SC2SoyProdBar2[, c(1,3,5)]
SC3SoyProdBar2 = SC3SoyProdBar
SC3SoyProdBar2$Scenario = rep("Scenario 3",7)
SC3SoyProdBar2 = SC3SoyProdBar2[, c(1,3,5)]
SC4SoyProdBar2 = SC4SoyProdBar
SC4SoyProdBar2$Scenario = rep("Scenario 4",7)
SC4SoyProdBar2 = SC4SoyProdBar2[, c(1,3,5)]
SoyProd3D = rbind(SC1SoyProdBar2,SC2SoyProdBar2,SC3SoyProdBar2,SC4SoyProdBar2)
SoyProd3D = SoyProd3D[, c(3,1,2)]
#oops, needs to be matrix, not dataframe
#this code from here to 10 year doesn't work after adding CDN, would be easy to fix but not necessary
SoyProd3DM = matrix(c(SoyProd3D[1:6,3],SoyProd3D[7:12,3],SoyProd3D[13:18,3],SoyProd3D[19:24,3]), ncol = 4)
colnames(SoyProd3DM) = c("Scenario 1","Scenario 2","Scenario 3","Scenario 4")
rownames(SoyProd3DM) = SC1SoyProdBar2$Insecticide

hist3D(z = SoyProd3DM, expand = 0.6, border = "black")

#reorder rows
SoyProd3DMR = SoyProd3DM[c(1,6,5,2,4,3),]
hist3D(x = 1:6, y = 1:4, z = SoyProd3DMR, expand = 0.6, border = "black", space = 0.2)

text3D(x = 1:6, y = 1:4, z = SoyProd3DMR,
       labels = colnames(SoyProd3DMR),
       add = TRUE, adj = 0)





#############  10 YEAR Aerial Spraying Soybean Aphid Figures ###############################

#Only within drift zone
#Scenario 3

#for 10 years, x = pesticide, y = production, 
#fill = North, Central, South
#North = 5 spray applications over 10 yrs, Central = 3 over 10 yrs, South = 1 over 10 yrs, Worst-case = 10 of 10

SoyProd3D[15:21,] #scen 3

SoyProd.10yr.S3 = data.frame(matrix(nrow = 28, ncol = 4))
colnames(SoyProd.10yr.S3) = c("Region","Insecticide","1Year","TenYear")
SoyProd.10yr.S3$Region = c(rep("North",7),rep("Central",7),rep("South",7),rep("Worst-Case",7))
SoyProd.10yr.S3$Region = factor(SoyProd.10yr.S3$Region, levels = c("Worst-Case","North","Central","South"))
SoyProd.10yr.S3$Insecticide = SoyProd3D[1:28,2]
SoyProd.10yr.S3$`1Year` = rep(SoyProd3D[15:21,3],4)
#note that SoyProd.10yr.S3 has the same numbers repeated 4x in column `1Year`, which is why code below seems inconsistent
SoyProd.10yr.S3[1:7,4] = (5*SoyProd.10yr.S3[8:14,3])+5*SoyProd.10yr.S3[8,3]
SoyProd.10yr.S3[8:14,4] = (3*SoyProd.10yr.S3[8:14,3])+7*SoyProd.10yr.S3[8,3]
SoyProd.10yr.S3[15:21,4] = (1*SoyProd.10yr.S3[8:14,3])+9*SoyProd.10yr.S3[8,3]
SoyProd.10yr.S3[22:28,4] = 10*SoyProd.10yr.S3[1:7,3]

p18 = ggplot(data = SoyProd.10yr.S3, aes(x = Insecticide, y = TenYear, fill = Region)) + 
  geom_bar(stat="identity", position = "dodge", color = "black") +
  scale_y_continuous(name = "Ten-Year Adult Monarch Production", labels = comma) +
  scale_fill_brewer(palette = "RdYlBu") +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title = element_text(size = 20), axis.text = element_text(size = 16), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black")
  )
p18


#Scenario 1

SoyProd3D[1:7,] #scen 1

SoyProd.10yr.S1 = data.frame(matrix(nrow = 28, ncol = 4))
colnames(SoyProd.10yr.S1) = c("Region","Insecticide","1Year","TenYear")
SoyProd.10yr.S1$Region = c(rep("North",7),rep("Central",7),rep("South",7),rep("Worst-Case",7))
SoyProd.10yr.S1$Region = factor(SoyProd.10yr.S1$Region, levels = c("Worst-Case","North","Central","South"))
SoyProd.10yr.S1$Insecticide = SoyProd3D[1:28,2]
SoyProd.10yr.S1$`1Year` = rep(SoyProd3D[1:7,3],4)
SoyProd.10yr.S1[1:7,4] = (5*SoyProd.10yr.S1[8:14,3])+5*SoyProd.10yr.S1[8,3]
SoyProd.10yr.S1[8:14,4] = (3*SoyProd.10yr.S1[8:14,3])+7*SoyProd.10yr.S1[8,3]
SoyProd.10yr.S1[15:21,4] = (1*SoyProd.10yr.S1[8:14,3])+9*SoyProd.10yr.S1[8,3]
SoyProd.10yr.S1[22:28,4] = 10*SoyProd.10yr.S1[1:7,3]

p19 = ggplot(data = SoyProd.10yr.S1, aes(x = Insecticide, y = TenYear, fill = Region)) + 
  geom_bar(stat="identity", position = "dodge", color = "black") +
  scale_y_continuous(name = "Ten-Year Adult Monarch Production", labels = comma) +
  scale_fill_brewer(palette = "RdYlBu") +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title = element_text(size = 20), axis.text = element_text(size = 16), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black")
  )
p19


#Scenario 2

SoyProd3D[8:14,] #scen 2

SoyProd.10yr.S2 = data.frame(matrix(nrow = 28, ncol = 4))
colnames(SoyProd.10yr.S2) = c("Region","Insecticide","1Year","TenYear")
SoyProd.10yr.S2$Region = c(rep("North",7),rep("Central",7),rep("South",7),rep("Worst-Case",7))
SoyProd.10yr.S2$Region = factor(SoyProd.10yr.S2$Region, levels = c("Worst-Case","North","Central","South"))
SoyProd.10yr.S2$Insecticide = SoyProd3D[1:28,2]
SoyProd.10yr.S2$`1Year` = rep(SoyProd3D[8:14,3],4)
SoyProd.10yr.S2[1:7,4] = (5*SoyProd.10yr.S2[8:14,3])+5*SoyProd.10yr.S2[8,3]
SoyProd.10yr.S2[8:14,4] = (3*SoyProd.10yr.S2[8:14,3])+7*SoyProd.10yr.S2[8,3]
SoyProd.10yr.S2[15:21,4] = (1*SoyProd.10yr.S2[8:14,3])+9*SoyProd.10yr.S2[8,3]
SoyProd.10yr.S2[22:28,4] = 10*SoyProd.10yr.S2[1:7,3]

p20 = ggplot(data = SoyProd.10yr.S2, aes(x = Insecticide, y = TenYear, fill = Region)) + 
  geom_bar(stat="identity", position = "dodge", color = "black") +
  scale_y_continuous(name = "Ten-Year Adult Monarch Production", labels = comma) +
  scale_fill_brewer(palette = "RdYlBu") +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title = element_text(size = 20), axis.text = element_text(size = 16), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black")
  )
p20


#Scenario 4

SoyProd3D[22:28,] #scen 4

SoyProd.10yr.S4 = data.frame(matrix(nrow = 28, ncol = 4))
colnames(SoyProd.10yr.S4) = c("Region","Insecticide","1Year","TenYear")
SoyProd.10yr.S4$Region = c(rep("North",7),rep("Central",7),rep("South",7),rep("Worst-Case",7))
SoyProd.10yr.S4$Region = factor(SoyProd.10yr.S4$Region, levels = c("Worst-Case","North","Central","South"))
SoyProd.10yr.S4$Insecticide = SoyProd3D[1:28,2]
SoyProd.10yr.S4$`1Year` = rep(SoyProd3D[22:28,3],4)
SoyProd.10yr.S4[1:7,4] = (5*SoyProd.10yr.S4[1:7,3])+5*SoyProd.10yr.S4[8,3]
SoyProd.10yr.S4[8:14,4] = (3*SoyProd.10yr.S4[1:7,3])+7*SoyProd.10yr.S4[8,3]
SoyProd.10yr.S4[15:21,4] = (1*SoyProd.10yr.S4[1:7,3])+9*SoyProd.10yr.S4[8,3]
SoyProd.10yr.S4[22:28,4] = 10*SoyProd.10yr.S4[1:7,3]

p21 = ggplot(data = SoyProd.10yr.S4, aes(x = Insecticide, y = TenYear, fill = Region)) + 
  geom_bar(stat="identity", position = "dodge", color = "black") +
  scale_y_continuous(name = "Ten-Year Adult Monarch Production", labels = comma) +
  scale_fill_brewer(palette = "RdYlBu") +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title = element_text(size = 20), axis.text = element_text(size = 16), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black")
  )
p21


############# 3D graphs per region ####################

#North
SoyProd.10yr.S1[1:7,4]
SoyProd.10yr.S2[1:7,4]
SoyProd.10yr.S3[1:7,4]
SoyProd.10yr.S4[1:7,4]

SoyProd.10yr.3D.North = matrix(c(SoyProd.10yr.S1[1:7,4],SoyProd.10yr.S2[1:7,4],SoyProd.10yr.S3[1:7,4],SoyProd.10yr.S4[1:7,4]), ncol = 4)
colnames(SoyProd.10yr.3D.North) = c("Scenario 1","Scenario 2","Scenario 3","Scenario 4")
rownames(SoyProd.10yr.3D.North) = SC1SoyProdBar2$Insecticide
SoyProd.10yr.3D.North

hist3D(z = SoyProd.10yr.3D.North, expand = 0.6, border = "black")

SoyProd.10yr.3D.North.R = SoyProd.10yr.3D.North[c(1,6,5,7,2,4,3),] #reorder pesticides
SoyProd.10yr.3D.North.R = SoyProd.10yr.3D.North.R[,c(4,1,3,2)]#reorder scenarios

#all commands I know
hist3D(x = seq(1, nrow(SoyProd.10yr.3D.North.R), by = 1), y = seq(1, ncol(SoyProd.10yr.3D.North.R), by = 1), 
       z = SoyProd.10yr.3D.North.R, expand = 0.6, border = "black", space = 0.2, ticktype = "detailed", colkey = FALSE, bty = "g",
       xlab = "Pesticide", ylab = "Scenario", zlab = "Adult Monarchs", main = "Adult Monarch Production", zlim = c(0,max(SoyProd.10yr.3D.North.R)),
       shade = 0.2, col = rev(jet.col())
       )

#current version
hist3D(x = seq(1, nrow(SoyProd.10yr.3D.North.R), by = 1), y = seq(1, ncol(SoyProd.10yr.3D.North.R), by = 1), 
       z = SoyProd.10yr.3D.North.R, expand = 0.6, border = "black", space = 0.2, bty = "g", colkey = FALSE, ticktype = "detailed",
       xlab = "", ylab = "", zlab = "", main = "Adult Monarch Production", zlim = c(0,max(SoyProd.10yr.3D.North.R)),
       col = rev(jet.col()))

#for some reason appears x,y, and z need to be same length as labels vector
#somehow, the vectors for x and y specify the location of the labels, apparently they are the x,y,z coords for each label
#6 rownames
text3D(x = 1:7, y = rep(0.2,7), z = rep(0,7), labels = rownames(SoyProd.10yr.3D.North.R), adj = 1.9, add = TRUE)
#4 colnames
text3D(x = rep(8.3,4), y = c(0.92,1.92,2.92,3.92), z = rep(0,4), labels = colnames(SoyProd.10yr.3D.North.R), adj = 0, add = TRUE)

write.csv(SoyProd.10yr.3D.North.R, file = "SoyProd.10yr.3D.North.R.csv")


#Central
SoyProd.10yr.S1[8:14,4]
SoyProd.10yr.S2[8:14,4]
SoyProd.10yr.S3[8:14,4]
SoyProd.10yr.S4[8:14,4]

SoyProd.10yr.3D.Cent = matrix(c(SoyProd.10yr.S1[8:14,4],SoyProd.10yr.S2[8:14,4],SoyProd.10yr.S3[8:14,4],SoyProd.10yr.S4[8:14,4]), ncol = 4)
colnames(SoyProd.10yr.3D.Cent) = c("Scenario 1","Scenario 2","Scenario 3","Scenario 4")
rownames(SoyProd.10yr.3D.Cent) = SC1SoyProdBar2$Insecticide
SoyProd.10yr.3D.Cent
#reord
SoyProd.10yr.3D.Cent.R = SoyProd.10yr.3D.Cent[c(1,6,5,7,2,4,3),] #reorder pesticides
SoyProd.10yr.3D.Cent.R = SoyProd.10yr.3D.Cent.R[,c(4,1,3,2)]#reorder scenarios
SoyProd.10yr.3D.Cent.R

hist3D(x = seq(1, nrow(SoyProd.10yr.3D.Cent.R), by = 1), y = seq(1, ncol(SoyProd.10yr.3D.Cent.R), by = 1), 
       z = SoyProd.10yr.3D.Cent.R, expand = 0.6, border = "black", space = 0.2, bty = "g", colkey = FALSE, ticktype = "detailed",
       xlab = "", ylab = "", zlab = "", main = "Adult Monarch Production", zlim = c(0,max(SoyProd.10yr.3D.Cent.R)),
       col = rev(jet.col()))
text3D(x = 1:7, y = rep(0.2,7), z = rep(0,7), labels = rownames(SoyProd.10yr.3D.Cent.R), adj = 1.9, add = TRUE)
text3D(x = rep(8.3,4), y = c(0.92,1.92,2.92,3.92), z = rep(0,4), labels = colnames(SoyProd.10yr.3D.Cent.R), adj = 0, add = TRUE)



#South
SoyProd.10yr.S1[15:21,4]
SoyProd.10yr.S2[15:21,4]
SoyProd.10yr.S3[15:21,4]
SoyProd.10yr.S4[15:21,4]

SoyProd.10yr.3D.South = matrix(c(SoyProd.10yr.S1[15:21,4],SoyProd.10yr.S2[15:21,4],SoyProd.10yr.S3[15:21,4],SoyProd.10yr.S4[15:21,4]), ncol = 4)
colnames(SoyProd.10yr.3D.South) = c("Scenario 1","Scenario 2","Scenario 3","Scenario 4")
rownames(SoyProd.10yr.3D.South) = SC1SoyProdBar2$Insecticide
SoyProd.10yr.3D.South
#reord
SoyProd.10yr.3D.South.R = SoyProd.10yr.3D.South[c(1,6,5,7,2,4,3),] #reorder pesticides
SoyProd.10yr.3D.South.R = SoyProd.10yr.3D.South.R[,c(4,1,3,2)]#reorder scenarios
SoyProd.10yr.3D.South.R

hist3D(x = seq(1, nrow(SoyProd.10yr.3D.South.R), by = 1), y = seq(1, ncol(SoyProd.10yr.3D.South.R), by = 1), 
       z = SoyProd.10yr.3D.South.R, expand = 0.6, border = "black", space = 0.2, bty = "g", colkey = FALSE, ticktype = "detailed",
       xlab = "", ylab = "", zlab = "", main = "Adult Monarch Production", zlim = c(0,max(SoyProd.10yr.3D.South.R)),
       col = rev(jet.col()))
text3D(x = 1:7, y = rep(0.2,7), z = rep(0,7), labels = rownames(SoyProd.10yr.3D.South.R), adj = 1.9, add = TRUE)
text3D(x = rep(8.3,4), y = c(0.92,1.92,2.92,3.92), z = rep(0,4), labels = colnames(SoyProd.10yr.3D.South.R), adj = 0, add = TRUE)



#Worst
SoyProd.10yr.S1[22:28,4]
SoyProd.10yr.S2[22:28,4]
SoyProd.10yr.S3[22:28,4]
SoyProd.10yr.S4[22:28,4]

SoyProd.10yr.3D.Worst = matrix(c(SoyProd.10yr.S1[22:28,4],SoyProd.10yr.S2[22:28,4],SoyProd.10yr.S3[22:28,4],SoyProd.10yr.S4[22:28,4]), ncol = 4)
colnames(SoyProd.10yr.3D.Worst) = c("Scenario 1","Scenario 2","Scenario 3","Scenario 4")
rownames(SoyProd.10yr.3D.Worst) = SC1SoyProdBar2$Insecticide
SoyProd.10yr.3D.Worst
#reord
SoyProd.10yr.3D.Worst.R = SoyProd.10yr.3D.Worst[c(1,6,5,7,2,4,3),] #reorder pesticides
SoyProd.10yr.3D.Worst.R = SoyProd.10yr.3D.Worst.R[,c(4,1,3,2)]#reorder scenarios
SoyProd.10yr.3D.Worst.R

hist3D(x = seq(1, nrow(SoyProd.10yr.3D.Worst.R), by = 1), y = seq(1, ncol(SoyProd.10yr.3D.Worst.R), by = 1), 
       z = SoyProd.10yr.3D.Worst.R, expand = 0.6, border = "black", space = 0.2, bty = "g", colkey = FALSE, ticktype = "detailed",
       xlab = "", ylab = "", zlab = "", main = "Adult Monarch Production", zlim = c(0,max(SoyProd.10yr.3D.Worst.R)),
       col = rev(jet.col()))
text3D(x = 1:7, y = rep(0.2,7), z = rep(0,7), labels = rownames(SoyProd.10yr.3D.Worst.R), adj = 1.9, add = TRUE)
text3D(x = rep(8.3,4), y = c(0.92,1.92,2.92,3.92), z = rep(0,4), labels = colnames(SoyProd.10yr.3D.Worst.R), adj = 0, add = TRUE)











######## RESULTS - Soybean Aphid with High Boom Spraying  ##################

#TISM generated from Table S9

#the pattern is (0.40*S.CU + 0.60)*S.DE
#0.40 is the proportion of individuals exposed to spray drift
#where S.CU is cuticular survival rate and S.DE is dietary exposure survival rate

#BCF
TISM[13,] = c(((0.40*1+0.60)*0.47),((0.40*0+0.60)*0.47),((0.40*0+0.60)*0.47),((0.40*0.10+0.60)*0.90),((0.40*0.10+0.60)*0.90),((0.40*0.56+0.60)*0.69),1)

#CTR
TISM[13,] = c(((0.40*1+0.60)*0.04),((0.40*0.02+0.60)*0.04),((0.40*0.02+0.60)*0.04),((0.40*0.31+0.60)*0.23),((0.40*0.31+0.60)*0.23),((0.40*0.71+0.60)*0.57),1)

#CFS
TISM[13,] = c(((0.40*1+0.60)*0.38),((0.40*0.99+0.60)*0.38),((0.40*0.99+0.60)*0.38),((0.40*0.98+0.60)*0.36),((0.40*0.98+0.60)*0.36),((0.40*1.00+0.60)*0.47),1)

#IMI
TISM[13,] = c(((0.40*1+0.60)*0.98),((0.40*0.99+0.60)*0.98),((0.40*0.99+0.60)*0.98),((0.40*0.98+0.60)*0.99),((0.40*0.98+0.60)*0.99),((0.40*1.00+0.60)*0.87),1)

#TMX
TISM[13,] = c(((0.40*1+0.60)*0.99),((0.40*0.98+0.60)*0.99),((0.40*0.98+0.60)*0.99),((0.40*1.00+0.60)*0.97),((0.40*1.00+0.60)*0.97),((0.40*1.00+0.60)*0.99),1)

#CDN
TISM[13,] = c(((0.40*1+0.60)*1.00),((0.40*0.07+0.60)*1.00),((0.40*0.07+0.60)*1.00),((0.40*0.51+0.60)*1.00),((0.40*0.51+0.60)*1.00),((0.40*0.98+0.60)*0.47),1)



########################################## SCENARIO 1 HIGH BOOM ################################################

#Scenario 1 Soybean Aphid Monarch Production with HIGH BOOM spraying
SC1SoyProdHB = data.frame(matrix(nrow = 3, ncol = 6))
colnames(SC1SoyProdHB) = c("No Drift","BCF","CTR","CFS","IMI","TMX")
rownames(SC1SoyProdHB) = c("Outside Drift Zone Prod","Inside Drift Zone Prod","Total Production")

#Baseline - no insecticide, TISM all = 1
TISM[13,]=c(1,1,1,1,1,1,1)
#total eggs
SC1SoyProdHB[3,1] = sum(NDD[,8]) 
#eggs laid out side drift zone
SC1SoyProdHB[1,1:6] = sum(NDD[,8])
#eggs laid inside drift zone
SC1SoyProdHB[2,1] = SC1SoyProdHB[3,1]-SC1SoyProdHB[1,1]

#BCF - used TISM from above
SC1SoyProdHB[2,2] = sum(NDD[,8]) 
SC1SoyProdHB[3,2] = SC1SoyProdHB[1,2]+SC1SoyProdHB[2,2]

#CTR
SC1SoyProdHB[2,3] = sum(NDD[,8]) 
SC1SoyProdHB[3,3] = SC1SoyProdHB[1,3]+SC1SoyProdHB[2,3]

#CFS
SC1SoyProdHB[2,4] = sum(NDD[,8]) 
SC1SoyProdHB[3,4] = SC1SoyProdHB[1,4]+SC1SoyProdHB[2,4]

#IMI
SC1SoyProdHB[2,5] = sum(NDD[,8]) 
SC1SoyProdHB[3,5] = SC1SoyProdHB[1,5]+SC1SoyProdHB[2,5]

#TMX
SC1SoyProdHB[2,6] = sum(NDD[,8]) 
SC1SoyProdHB[3,6] = SC1SoyProdHB[1,6]+SC1SoyProdHB[2,6]

#CDN
SC1SoyProdHBCDN = c()
SC1SoyProdHBCDN[1] = SC1SoyProdHB[1,6]
SC1SoyProdHBCDN[2] = sum(NDD[,8])
SC1SoyProdHBCDN[3] = SC1SoyProdHBCDN[1]+SC1SoyProdHBCDN[2]
SC1SoyProdHB$CDN = SC1SoyProdHBCDN


#figures
library(ggplot2)
library(scales)
SC1SoyProdHBBar = data.frame(matrix(nrow = 7, ncol = 4))
colnames(SC1SoyProdHBBar) = c("Insecticide","Out","In","Total")
SC1SoyProdHBBar$Insecticide = c("No Drift","BCF","CTR","CFS","IMI","TMX","CDN")
SC1SoyProdHBBar$In = c(SC1SoyProdHB[2,])
SC1SoyProdHBBar$Out = c(SC1SoyProdHB[1,])
SC1SoyProdHBBar$Total = c(SC1SoyProdHB[3,])
SC1SoyProdHBBar$Insecticide = factor(SC1SoyProdHBBar$Insecticide, levels = c("No Drift","BCF","CTR","CFS","IMI","TMX","CDN"))
#for some reason SC1SoyProdHBBar is all lists, so data labels don't work
SC1SoyProdHBBar$Total = round(unlist(SC1SoyProdHBBar$Total))
SC1SoyProdHBBar$In = round(unlist(SC1SoyProdHBBar$In))

p10 = ggplot(data = SC1SoyProdHBBar, aes(x = Insecticide, y = Total)) + geom_bar(stat="identity") +
  geom_text(aes(label = comma(Total)), vjust = 1.5, color = "white") +
  scale_y_continuous(name = "Adult Monarch Production", labels = comma) +
  theme(axis.title.x = element_text(size = 20, margin = margin(t=10)),
        axis.title.y = element_text(size = 20, margin = margin(r=10)),
        axis.title = element_text(size = 20), axis.text = element_text(size = 16), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black")
  )
p10

p11 = ggplot(data = SC1SoyProdHBBar, aes(x = Insecticide, y = In)) + geom_bar(stat="identity") + ylab("Adult Monarch Production") +
  geom_text(aes(label = comma(In)), vjust = 1.5, color = "white") +
  scale_y_continuous(name = "Adult Monarch Production", labels = comma) +
  theme(axis.title.x = element_text(size = 20, margin = margin(t=10)),
        axis.title.y = element_text(size = 20, margin = margin(r=10)),
        axis.title = element_text(size = 20), axis.text = element_text(size = 16), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black")
  )
p11


####### SCENARIO 2 HIGH BOOM ########################

#Scenario 2 Soybean Aphid Monarch Production with HIGH BOOM spraying
SC2SoyProdHB = data.frame(matrix(nrow = 3, ncol = 6))
colnames(SC2SoyProdHB) = c("No Drift","BCF","CTR","CFS","IMI","TMX")
rownames(SC2SoyProdHB) = c("Outside Drift Zone Prod","Inside Drift Zone Prod","Total Production")

#Baseline - no insecticide, TISM all = 1
TISM[13,]=c(1,1,1,1,1,1,1)
#total eggs
SC2SoyProdHB[3,1] = sum(NDD[,8]) 
#eggs laid out side drift zone
SC2SoyProdHB[1,1:6] = sum(NDD[,8])
#eggs laid inside drift zone
SC2SoyProdHB[2,1] = SC2SoyProdHB[3,1]-SC2SoyProdHB[1,1]

#BCF - used TISM from above
SC2SoyProdHB[2,2] = sum(NDD[,8]) 
SC2SoyProdHB[3,2] = SC2SoyProdHB[1,2]+SC2SoyProdHB[2,2]

#CTR
SC2SoyProdHB[2,3] = sum(NDD[,8]) 
SC2SoyProdHB[3,3] = SC2SoyProdHB[1,3]+SC2SoyProdHB[2,3]

#CFS
SC2SoyProdHB[2,4] = sum(NDD[,8]) 
SC2SoyProdHB[3,4] = SC2SoyProdHB[1,4]+SC2SoyProdHB[2,4]

#IMI
SC2SoyProdHB[2,5] = sum(NDD[,8]) 
SC2SoyProdHB[3,5] = SC2SoyProdHB[1,5]+SC2SoyProdHB[2,5]

#TMX
SC2SoyProdHB[2,6] = sum(NDD[,8]) 
SC2SoyProdHB[3,6] = SC2SoyProdHB[1,6]+SC2SoyProdHB[2,6]

#CDN
SC2SoyProdHBCDN = c()
SC2SoyProdHBCDN[1] = SC2SoyProdHB[1,6]
SC2SoyProdHBCDN[2] = sum(NDD[,8])
SC2SoyProdHBCDN[3] = SC2SoyProdHBCDN[1]+SC2SoyProdHBCDN[2]
SC2SoyProdHB$CDN = SC2SoyProdHBCDN

#figures
SC2SoyProdHBBar = data.frame(matrix(nrow = 7, ncol = 4))
colnames(SC2SoyProdHBBar) = c("Insecticide","Out","In","Total")
SC2SoyProdHBBar$Insecticide = c("No Drift","BCF","CTR","CFS","IMI","TMX","CDN")
SC2SoyProdHBBar$In = c(SC2SoyProdHB[2,])
SC2SoyProdHBBar$Out = c(SC2SoyProdHB[1,])
SC2SoyProdHBBar$Total = c(SC2SoyProdHB[3,])
SC2SoyProdHBBar$Insecticide = factor(SC2SoyProdHBBar$Insecticide, levels = c("No Drift","BCF","CTR","CFS","IMI","TMX","CDN"))
#for some reason SC2SoyProdHBBar is all lists, so data labels don't work
SC2SoyProdHBBar$Total = round(unlist(SC2SoyProdHBBar$Total))
SC2SoyProdHBBar$In = round(unlist(SC2SoyProdHBBar$In))

p12 = ggplot(data = SC2SoyProdHBBar, aes(x = Insecticide, y = Total)) + geom_bar(stat="identity") +
  geom_text(aes(label = comma(Total)), vjust = 1.5, color = "white") +
  scale_y_continuous(name = "Adult Monarch Production", labels = comma) +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title = element_text(size = 20), axis.text = element_text(size = 16), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black")
  )
p12

p13 = ggplot(data = SC2SoyProdHBBar, aes(x = Insecticide, y = In)) + geom_bar(stat="identity") + ylab("Adult Monarch Production") +
  geom_text(aes(label = comma(In)), vjust = 1.5, color = "white") +
  scale_y_continuous(name = "Adult Monarch Production", labels = comma) +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title = element_text(size = 20), axis.text = element_text(size = 16), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black")
  )
p13


####### SCENARIO 3 HIGH BOOM ########################

#Scenario 3 Soybean Aphid Monarch Production with HIGH BOOM spraying
SC3SoyProdHB = data.frame(matrix(nrow = 3, ncol = 6))
colnames(SC3SoyProdHB) = c("No Drift","BCF","CTR","CFS","IMI","TMX")
rownames(SC3SoyProdHB) = c("Outside Drift Zone Prod","Inside Drift Zone Prod","Total Production")

#Baseline - no insecticide, TISM all = 1
TISM[13,]=c(1,1,1,1,1,1,1)
#total eggs
SC3SoyProdHB[3,1] = sum(NDD[,8]) 
#eggs laid out side drift zone
SC3SoyProdHB[1,1:6] = sum(NDD[,8])
#eggs laid inside drift zone
SC3SoyProdHB[2,1] = SC3SoyProdHB[3,1]-SC3SoyProdHB[1,1]

#BCF - used TISM from above
SC3SoyProdHB[2,2] = sum(NDD[,8]) 
SC3SoyProdHB[3,2] = SC3SoyProdHB[1,2]+SC3SoyProdHB[2,2]

#CTR
SC3SoyProdHB[2,3] = sum(NDD[,8]) 
SC3SoyProdHB[3,3] = SC3SoyProdHB[1,3]+SC3SoyProdHB[2,3]

#CFS
SC3SoyProdHB[2,4] = sum(NDD[,8]) 
SC3SoyProdHB[3,4] = SC3SoyProdHB[1,4]+SC3SoyProdHB[2,4]

#IMI
SC3SoyProdHB[2,5] = sum(NDD[,8]) 
SC3SoyProdHB[3,5] = SC3SoyProdHB[1,5]+SC3SoyProdHB[2,5]

#TMX
SC3SoyProdHB[2,6] = sum(NDD[,8]) 
SC3SoyProdHB[3,6] = SC3SoyProdHB[1,6]+SC3SoyProdHB[2,6]

#CDN
SC3SoyProdHBCDN = c()
SC3SoyProdHBCDN[1] = SC3SoyProdHB[1,6]
SC3SoyProdHBCDN[2] = sum(NDD[,8])
SC3SoyProdHBCDN[3] = SC3SoyProdHBCDN[1]+SC3SoyProdHBCDN[2]
SC3SoyProdHB$CDN = SC3SoyProdHBCDN

#figures
SC3SoyProdHBBar = data.frame(matrix(nrow = 7, ncol = 4))
colnames(SC3SoyProdHBBar) = c("Insecticide","Out","In","Total")
SC3SoyProdHBBar$Insecticide = c("No Drift","BCF","CTR","CFS","IMI","TMX","CDN")
SC3SoyProdHBBar$In = c(SC3SoyProdHB[2,])
SC3SoyProdHBBar$Out = c(SC3SoyProdHB[1,])
SC3SoyProdHBBar$Total = c(SC3SoyProdHB[3,])
SC3SoyProdHBBar$Insecticide = factor(SC3SoyProdHBBar$Insecticide, levels = c("No Drift","BCF","CTR","CFS","IMI","TMX","CDN"))
#for some reason SC3SoyProdHBBar is all lists, so data labels don't work
SC3SoyProdHBBar$Total = round(unlist(SC3SoyProdHBBar$Total))
SC3SoyProdHBBar$In = round(unlist(SC3SoyProdHBBar$In))

p14 = ggplot(data = SC3SoyProdHBBar, aes(x = Insecticide, y = Total)) + geom_bar(stat="identity") +
  geom_text(aes(label = comma(Total)), vjust = 1.5, color = "white") +
  scale_y_continuous(name = "Adult Monarch Production", labels = comma) +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title = element_text(size = 20), axis.text = element_text(size = 16), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black")
  )
p14

p15 = ggplot(data = SC3SoyProdHBBar, aes(x = Insecticide, y = In)) + geom_bar(stat="identity") + ylab("Adult Monarch Production") +
  geom_text(aes(label = comma(In)), vjust = 1.5, color = "white") +
  scale_y_continuous(name = "Adult Monarch Production", labels = comma) +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title = element_text(size = 20), axis.text = element_text(size = 16), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black")
  )
p15



####### SCENARIO 4 HIGH BOOM ########################

#Scenario 4 Soybean Aphid Monarch Production with HIGH BOOM spraying
SC4SoyProdHB = data.frame(matrix(nrow = 3, ncol = 6))
colnames(SC4SoyProdHB) = c("No Drift","BCF","CTR","CFS","IMI","TMX")
rownames(SC4SoyProdHB) = c("Outside Drift Zone Prod","Inside Drift Zone Prod","Total Production")

#Baseline - no insecticide, TISM all = 1
TISM[13,]=c(1,1,1,1,1,1,1)
#total eggs
SC4SoyProdHB[3,1] = sum(NDD[,8]) 
#eggs laid out side drift zone
SC4SoyProdHB[1,1:6] = sum(NDD[,8])
#eggs laid inside drift zone
SC4SoyProdHB[2,1] = sum(NDD[,8])
SC4SoyProdHB[3,1]-SC4SoyProdHB[1,1]

#BCF - used TISM from above
SC4SoyProdHB[2,2] = sum(NDD[,8]) 
SC4SoyProdHB[3,2] = SC4SoyProdHB[1,2]+SC4SoyProdHB[2,2]

#CTR
SC4SoyProdHB[2,3] = sum(NDD[,8]) 
SC4SoyProdHB[3,3] = SC4SoyProdHB[1,3]+SC4SoyProdHB[2,3]

#CFS
SC4SoyProdHB[2,4] = sum(NDD[,8]) 
SC4SoyProdHB[3,4] = SC4SoyProdHB[1,4]+SC4SoyProdHB[2,4]

#IMI
SC4SoyProdHB[2,5] = sum(NDD[,8]) 
SC4SoyProdHB[3,5] = SC4SoyProdHB[1,5]+SC4SoyProdHB[2,5]

#TMX
SC4SoyProdHB[2,6] = sum(NDD[,8]) 
SC4SoyProdHB[3,6] = SC4SoyProdHB[1,6]+SC4SoyProdHB[2,6]

#CDN
SC4SoyProdHBCDN = c()
SC4SoyProdHBCDN[1] = SC4SoyProdHB[1,6]
SC4SoyProdHBCDN[2] = sum(NDD[,8])
SC4SoyProdHBCDN[3] = SC4SoyProdHBCDN[1]+SC4SoyProdHBCDN[2]
SC4SoyProdHB$CDN = SC4SoyProdHBCDN

#figures
SC4SoyProdHBBar = data.frame(matrix(nrow = 7, ncol = 4))
colnames(SC4SoyProdHBBar) = c("Insecticide","Out","In","Total")
SC4SoyProdHBBar$Insecticide = c("No Drift","BCF","CTR","CFS","IMI","TMX","CDN")
SC4SoyProdHBBar$In = c(SC4SoyProdHB[2,])
SC4SoyProdHBBar$Out = c(SC4SoyProdHB[1,])
SC4SoyProdHBBar$Total = c(SC4SoyProdHB[3,])
SC4SoyProdHBBar$Insecticide = factor(SC4SoyProdHBBar$Insecticide, levels = c("No Drift","BCF","CTR","CFS","IMI","TMX","CDN"))
#for some reason SC4SoyProdHBBar is all lists, so data labels don't work
SC4SoyProdHBBar$Total = round(unlist(SC4SoyProdHBBar$Total))
SC4SoyProdHBBar$In = round(unlist(SC4SoyProdHBBar$In))

p16 = ggplot(data = SC4SoyProdHBBar, aes(x = Insecticide, y = Total)) + geom_bar(stat="identity") +
  geom_text(aes(label = comma(Total)), vjust = 1.5, color = "white") +
  scale_y_continuous(name = "Adult Monarch Production", labels = comma) +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title = element_text(size = 20), axis.text = element_text(size = 16), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black")
  )
p16

p17 = ggplot(data = SC4SoyProdHBBar, aes(x = Insecticide, y = In)) + geom_bar(stat="identity") + ylab("Adult Monarch Production") +
  geom_text(aes(label = comma(In)), vjust = 1.5, color = "white") +
  scale_y_continuous(name = "Adult Monarch Production", labels = comma) +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title = element_text(size = 20), axis.text = element_text(size = 16), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black")
  )
p17





######################  10-year 3D graphs of High Boom Soybean Aphid Spray Drift #############################

#insecticide by scenario, per region

#scenario 1 10-yr calcs
SoyProdHB.10yr.S1 = data.frame(matrix(nrow = 28, ncol = 4))
colnames(SoyProdHB.10yr.S1) = c("Region","Insecticide","1Year","TenYear")
SoyProdHB.10yr.S1$Region = c(rep("North",7),rep("Central",7),rep("South",7),rep("Worst-Case",7))
SoyProdHB.10yr.S1$Region = factor(SoyProdHB.10yr.S1$Region, levels = c("Worst-Case","North","Central","South"))
SoyProdHB.10yr.S1$Insecticide = SoyProd3D[1:28,2]
SoyProdHB.10yr.S1$`1Year` = rep(SC1SoyProdHBBar[1:7,3],4)
SoyProdHB.10yr.S1[1:7,4] = (5*SoyProdHB.10yr.S1[1:7,3])+5*SoyProdHB.10yr.S1[1,3]
SoyProdHB.10yr.S1[8:14,4] = (3*SoyProdHB.10yr.S1[1:7,3])+7*SoyProdHB.10yr.S1[1,3]
SoyProdHB.10yr.S1[15:21,4] = (1*SoyProdHB.10yr.S1[1:7,3])+9*SoyProdHB.10yr.S1[1,3]
SoyProdHB.10yr.S1[22:28,4] = 10*SoyProdHB.10yr.S1[1:7,3]

#scenario 2 10-yr calcs
SoyProdHB.10yr.S2 = data.frame(matrix(nrow = 28, ncol = 4))
colnames(SoyProdHB.10yr.S2) = c("Region","Insecticide","1Year","TenYear")
SoyProdHB.10yr.S2$Region = c(rep("North",7),rep("Central",7),rep("South",7),rep("Worst-Case",7))
SoyProdHB.10yr.S2$Region = factor(SoyProdHB.10yr.S2$Region, levels = c("Worst-Case","North","Central","South"))
SoyProdHB.10yr.S2$Insecticide = SoyProd3D[1:28,2]
SoyProdHB.10yr.S2$`1Year` = rep(SC2SoyProdHBBar[1:7,3],4)
SoyProdHB.10yr.S2[1:7,4] = (5*SoyProdHB.10yr.S2[1:7,3])+5*SoyProdHB.10yr.S2[1,3]
SoyProdHB.10yr.S2[8:14,4] = (3*SoyProdHB.10yr.S2[1:7,3])+7*SoyProdHB.10yr.S2[1,3]
SoyProdHB.10yr.S2[15:21,4] = (1*SoyProdHB.10yr.S2[1:7,3])+9*SoyProdHB.10yr.S2[1,3]
SoyProdHB.10yr.S2[22:28,4] = 10*SoyProdHB.10yr.S2[1:7,3]

#scenario 3 10-yr calcs
SoyProdHB.10yr.S3 = data.frame(matrix(nrow = 28, ncol = 4))
colnames(SoyProdHB.10yr.S3) = c("Region","Insecticide","1Year","TenYear")
SoyProdHB.10yr.S3$Region = c(rep("North",7),rep("Central",7),rep("South",7),rep("Worst-Case",7))
SoyProdHB.10yr.S3$Region = factor(SoyProdHB.10yr.S3$Region, levels = c("Worst-Case","North","Central","South"))
SoyProdHB.10yr.S3$Insecticide = SoyProd3D[1:28,2]
SoyProdHB.10yr.S3$`1Year` = rep(SC3SoyProdHBBar[1:7,3],4)
SoyProdHB.10yr.S3[1:7,4] = (5*SoyProdHB.10yr.S3[1:7,3])+5*SoyProdHB.10yr.S3[1,3]
SoyProdHB.10yr.S3[8:14,4] = (3*SoyProdHB.10yr.S3[1:7,3])+7*SoyProdHB.10yr.S3[1,3]
SoyProdHB.10yr.S3[15:21,4] = (1*SoyProdHB.10yr.S3[1:7,3])+9*SoyProdHB.10yr.S3[1,3]
SoyProdHB.10yr.S3[22:28,4] = 10*SoyProdHB.10yr.S3[1:7,3]

#scenario 4 10-yr calcs
SoyProdHB.10yr.S4 = data.frame(matrix(nrow = 28, ncol = 4))
colnames(SoyProdHB.10yr.S4) = c("Region","Insecticide","1Year","TenYear")
SoyProdHB.10yr.S4$Region = c(rep("North",7),rep("Central",7),rep("South",7),rep("Worst-Case",7))
SoyProdHB.10yr.S4$Region = factor(SoyProdHB.10yr.S4$Region, levels = c("Worst-Case","North","Central","South"))
SoyProdHB.10yr.S4$Insecticide = SoyProd3D[1:28,2]
SoyProdHB.10yr.S4$`1Year` = rep(SC4SoyProdHBBar[1:7,3],4)
SoyProdHB.10yr.S4[1:7,4] = (5*SoyProdHB.10yr.S4[1:7,3])+5*SoyProdHB.10yr.S4[1,3]
SoyProdHB.10yr.S4[8:14,4] = (3*SoyProdHB.10yr.S4[1:7,3])+7*SoyProdHB.10yr.S4[1,3]
SoyProdHB.10yr.S4[15:21,4] = (1*SoyProdHB.10yr.S4[1:7,3])+9*SoyProdHB.10yr.S4[1,3]
SoyProdHB.10yr.S4[22:28,4] = 10*SoyProdHB.10yr.S4[1:7,3]


#North
SoyProdHB.10yr.S1[1:7,4]
SoyProdHB.10yr.S2[1:7,4]
SoyProdHB.10yr.S3[1:7,4]
SoyProdHB.10yr.S4[1:7,4]

SoyProdHB.10yr.3D.North = matrix(c(SoyProdHB.10yr.S1[1:7,4],SoyProdHB.10yr.S2[1:7,4],SoyProdHB.10yr.S3[1:7,4],SoyProdHB.10yr.S4[1:7,4]), ncol = 4)
colnames(SoyProdHB.10yr.3D.North) = c("Scenario 1","Scenario 2","Scenario 3","Scenario 4")
rownames(SoyProdHB.10yr.3D.North) = SC1SoyProdBar2$Insecticide
SoyProdHB.10yr.3D.North

SoyProdHB.10yr.3D.North.R = SoyProdHB.10yr.3D.North[c(1,6,5,7,2,4,3),]
SoyProdHB.10yr.3D.North.R = SoyProdHB.10yr.3D.North.R[,c(4,1,3,2)]
SoyProdHB.10yr.3D.North.R

hist3D(x = seq(1, nrow(SoyProdHB.10yr.3D.North.R), by = 1), y = seq(1, ncol(SoyProdHB.10yr.3D.North.R), by = 1), 
       z = SoyProdHB.10yr.3D.North.R, expand = 0.6, border = "black", space = 0.2, bty = "g", colkey = FALSE, ticktype = "detailed",
       xlab = "", ylab = "", zlab = "", main = "Adult Monarch Production", zlim = c(0,max(SoyProdHB.10yr.3D.North.R)),
       col = rev(jet.col()))
text3D(x = c(0.92,1.92,2.92,3.92,4.92,5.92,6.92), y = rep(0,7), z = rep(0,7), labels = rownames(SoyProdHB.10yr.3D.North.R), adj = 1.9, add = TRUE)
text3D(x = rep(8.3,4), y = c(0.92,1.92,2.92,3.92), z = rep(0,4), labels = colnames(SoyProdHB.10yr.3D.North.R), adj = 0, add = TRUE)


#Central
SoyProdHB.10yr.S1[8:14,4]
SoyProdHB.10yr.S2[8:14,4]
SoyProdHB.10yr.S3[8:14,4]
SoyProdHB.10yr.S4[8:14,4]

SoyProdHB.10yr.3D.Cent = matrix(c(SoyProdHB.10yr.S1[8:14,4],SoyProdHB.10yr.S2[8:14,4],SoyProdHB.10yr.S3[8:14,4],SoyProdHB.10yr.S4[8:14,4]), ncol = 4)
colnames(SoyProdHB.10yr.3D.Cent) = c("Scenario 1","Scenario 2","Scenario 3","Scenario 4")
rownames(SoyProdHB.10yr.3D.Cent) = SC1SoyProdBar2$Insecticide
SoyProdHB.10yr.3D.Cent

SoyProdHB.10yr.3D.Cent.R = SoyProdHB.10yr.3D.Cent[c(1,6,5,7,2,4,3),]
SoyProdHB.10yr.3D.Cent.R = SoyProdHB.10yr.3D.Cent.R[,c(4,1,3,2)]
SoyProdHB.10yr.3D.Cent.R

hist3D(x = seq(1, nrow(SoyProdHB.10yr.3D.Cent.R), by = 1), y = seq(1, ncol(SoyProdHB.10yr.3D.Cent.R), by = 1), 
       z = SoyProdHB.10yr.3D.Cent.R, expand = 0.6, border = "black", space = 0.2, bty = "g", colkey = FALSE, ticktype = "detailed",
       xlab = "", ylab = "", zlab = "", main = "Adult Monarch Production", zlim = c(0,max(SoyProdHB.10yr.3D.Cent.R)),
       col = rev(jet.col()))
text3D(x = 1:7, y = rep(0.2,7), z = rep(0,7), labels = rownames(SoyProdHB.10yr.3D.Cent.R), adj = 1.5, add = TRUE)
text3D(x = rep(7.2,4), y = c(0.92,1.92,2.92,3.92), z = rep(0,4), labels = colnames(SoyProdHB.10yr.3D.Cent.R), adj = 0, add = TRUE)


#South
SoyProdHB.10yr.S1[15:21,4]
SoyProdHB.10yr.S2[15:21,4]
SoyProdHB.10yr.S3[15:21,4]
SoyProdHB.10yr.S4[15:21,4]

SoyProdHB.10yr.3D.South = matrix(c(SoyProdHB.10yr.S1[15:21,4],SoyProdHB.10yr.S2[15:21,4],SoyProdHB.10yr.S3[15:21,4],SoyProdHB.10yr.S4[15:21,4]), ncol = 4)
colnames(SoyProdHB.10yr.3D.South) = c("Scenario 1","Scenario 2","Scenario 3","Scenario 4")
rownames(SoyProdHB.10yr.3D.South) = SC1SoyProdBar2$Insecticide
SoyProdHB.10yr.3D.South

SoyProdHB.10yr.3D.South.R = SoyProdHB.10yr.3D.South[c(1,6,5,7,2,4,3),]
SoyProdHB.10yr.3D.South.R = SoyProdHB.10yr.3D.South.R[,c(4,1,3,2)]
SoyProdHB.10yr.3D.South.R

hist3D(x = seq(1, nrow(SoyProdHB.10yr.3D.South.R), by = 1), y = seq(1, ncol(SoyProdHB.10yr.3D.South.R), by = 1), 
       z = SoyProdHB.10yr.3D.South.R, expand = 0.6, border = "black", space = 0.2, bty = "g", colkey = FALSE, ticktype = "detailed",
       xlab = "", ylab = "", zlab = "", main = "Adult Monarch Production", zlim = c(0,max(SoyProdHB.10yr.3D.South.R)),
       col = rev(jet.col()))
text3D(x = 1:7, y = rep(0.2,7), z = rep(0,7), labels = rownames(SoyProdHB.10yr.3D.South.R), adj = 1.5, add = TRUE)
text3D(x = rep(8.3,4), y = c(0.92,1.92,2.92,3.92), z = rep(0,4), labels = colnames(SoyProdHB.10yr.3D.South.R), adj = 0, add = TRUE)


#Worst-Case Scenario
SoyProdHB.10yr.S1[22:28,4]
SoyProdHB.10yr.S2[22:28,4]
SoyProdHB.10yr.S3[22:28,4]
SoyProdHB.10yr.S4[22:28,4]

SoyProdHB.10yr.3D.Worst = matrix(c(SoyProdHB.10yr.S1[22:28,4],SoyProdHB.10yr.S2[22:28,4],SoyProdHB.10yr.S3[22:28,4],SoyProdHB.10yr.S4[22:28,4]), ncol = 4)
colnames(SoyProdHB.10yr.3D.Worst) = c("Scenario 1","Scenario 2","Scenario 3","Scenario 4")
rownames(SoyProdHB.10yr.3D.Worst) = SC1SoyProdBar2$Insecticide
SoyProdHB.10yr.3D.Worst

SoyProdHB.10yr.3D.Worst.R = SoyProdHB.10yr.3D.Worst[c(1,6,5,7,2,4,3),]
SoyProdHB.10yr.3D.Worst.R = SoyProdHB.10yr.3D.Worst.R[,c(4,1,3,2)]
SoyProdHB.10yr.3D.Worst.R

hist3D(x = seq(1, nrow(SoyProdHB.10yr.3D.Worst.R), by = 1), y = seq(1, ncol(SoyProdHB.10yr.3D.Worst.R), by = 1), 
       z = SoyProdHB.10yr.3D.Worst.R, expand = 0.6, border = "black", space = 0.2, bty = "g", colkey = FALSE, ticktype = "detailed",
       xlab = "", ylab = "", zlab = "", main = "Adult Monarch Production", zlim = c(0,max(SoyProdHB.10yr.3D.Worst.R)),
       col = rev(jet.col()))
text3D(x = 1:7, y = rep(0.2,7), z = rep(0,7), labels = rownames(SoyProdHB.10yr.3D.Worst.R), adj = 1.5, add = TRUE)
text3D(x = rep(8.3,4), y = c(0.92,1.92,2.92,3.92), z = rep(0,4), labels = colnames(SoyProdHB.10yr.3D.Worst.R), adj = 0, add = TRUE)






###################### 10-year Graphs when 100% Mortality from Spray Drift in Some Years - Aerial then High Boom ############################

#first going to add it to regular graph as another row in pesticide scenarios
#not redoing above code, taking above code and remaking it with the 100% mortality pesticide scenario

#attempt to make the color ramp constant across graphs (successful)

#maximum range between max and min production
#max is always Scen 4, no drift = 46,840
#min is 0 in worst case, 100% mort
#min for the 8 graphs
min(SoyProd100.10yr.3D.North.R) #13780
min(SoyProd100.10yr.3D.Cent.R) #19292
min(SoyProd100.10yr.3D.South.R) #24804
min(SoyProd100.10yr.3D.Worst.R) #0
min(SoyProd100HB.10yr.3D.North.R) #13780
min(SoyProd100HB.10yr.3D.Cent.R) #19292
min(SoyProd100HB.10yr.3D.South.R) #24804
min(SoyProd100HB.10yr.3D.Worst.R) #0
#range divided by max = proportion of colors should be used
(max(SoyProd100.10yr.3D.North.R)-min(SoyProd100.10yr.3D.North.R))/max(SoyProd100.10yr.3D.North.R) #range = 33060, prop = 0.71
(max(SoyProd100.10yr.3D.Cent.R)-min(SoyProd100.10yr.3D.Cent.R))/max(SoyProd100.10yr.3D.Cent.R) #range = 27548, prop = 0.59
(max(SoyProd100.10yr.3D.South.R)-min(SoyProd100.10yr.3D.South.R))/max(SoyProd100.10yr.3D.South.R) #range = 22036, prop = 0.47
(max(SoyProd100.10yr.3D.Worst.R)-min(SoyProd100.10yr.3D.Worst.R))/max(SoyProd100.10yr.3D.Worst.R) #range = 46840, prop = 1

#base color ramp for maximum range
col.base = rev(jet.col(100))

#"col = " for first figure below - aerial, north
col1 = col.base[(100-71):100]
col2 = col.base[(100-59):100]
col3 = col.base[(100-47):100]
col4 = col.base[1:100] 
col5 = col.base[(100-71):100]
col6 = col.base[(100-59):100]
col7 = col.base[(100-47):100]
col8 = col.base[1:100]


#AERIAL SPRAYING

#modify this to SoyProd100.10yr.S1
SoyProdHB.10yr.S1

#scenario 1 10-yr calcs
SoyProd100.10yr.S1 = data.frame(matrix(nrow = 32, ncol = 4))
colnames(SoyProd100.10yr.S1) = c("Region","Insecticide","1Year","TenYear")
SoyProd100.10yr.S1$Region = c(rep("North",8),rep("Central",8),rep("South",8),rep("Worst-Case",8))
SoyProd100.10yr.S1$Region = factor(SoyProd100.10yr.S1$Region, levels = c("Worst-Case","North","Central","South"))
SoyProd100.10yr.S1$Insecticide = rep(c(as.character(SC1SoyProdBar$Insecticide),"100%"),4)
SoyProd100.10yr.S1$`1Year` = rep(c(SC1SoyProdBar[1:7,3],0),4)
SoyProd100.10yr.S1[1:8,4] = (5*SoyProd100.10yr.S1[1:8,3])+5*SoyProd100.10yr.S1[1,3]
SoyProd100.10yr.S1[9:16,4] = (3*SoyProd100.10yr.S1[1:8,3])+7*SoyProd100.10yr.S1[1,3]
SoyProd100.10yr.S1[17:24,4] = (1*SoyProd100.10yr.S1[1:8,3])+9*SoyProd100.10yr.S1[1,3]
SoyProd100.10yr.S1[25:32,4] = 10*SoyProd100.10yr.S1[1:8,3]

#scenario 2 10-yr calcs
SoyProd100.10yr.S2 = data.frame(matrix(nrow = 32, ncol = 4))
colnames(SoyProd100.10yr.S2) = c("Region","Insecticide","1Year","TenYear")
SoyProd100.10yr.S2$Region = c(rep("North",8),rep("Central",8),rep("South",8),rep("Worst-Case",8))
SoyProd100.10yr.S2$Region = factor(SoyProd100.10yr.S2$Region, levels = c("Worst-Case","North","Central","South"))
SoyProd100.10yr.S2$Insecticide = rep(c(as.character(SC1SoyProdBar$Insecticide),"100%"),4)
SoyProd100.10yr.S2$`1Year` = rep(c(SC2SoyProdBar[1:7,3],0),4)
SoyProd100.10yr.S2[1:8,4] = (5*SoyProd100.10yr.S2[1:8,3])+5*SoyProd100.10yr.S2[1,3]
SoyProd100.10yr.S2[9:16,4] = (3*SoyProd100.10yr.S2[1:8,3])+7*SoyProd100.10yr.S2[1,3]
SoyProd100.10yr.S2[17:24,4] = (1*SoyProd100.10yr.S2[1:8,3])+9*SoyProd100.10yr.S2[1,3]
SoyProd100.10yr.S2[25:32,4] = 10*SoyProd100.10yr.S2[1:8,3]

#scenario 3 10-yr calcs
SoyProd100.10yr.S3 = data.frame(matrix(nrow = 32, ncol = 4))
colnames(SoyProd100.10yr.S3) = c("Region","Insecticide","1Year","TenYear")
SoyProd100.10yr.S3$Region = c(rep("North",8),rep("Central",8),rep("South",8),rep("Worst-Case",8))
SoyProd100.10yr.S3$Region = factor(SoyProd100.10yr.S3$Region, levels = c("Worst-Case","North","Central","South"))
SoyProd100.10yr.S3$Insecticide = rep(c(as.character(SC1SoyProdBar$Insecticide),"100%"),4)
SoyProd100.10yr.S3$`1Year` = rep(c(SC3SoyProdBar[1:7,3],0),4)
SoyProd100.10yr.S3[1:8,4] = (5*SoyProd100.10yr.S3[1:8,3])+5*SoyProd100.10yr.S3[1,3]
SoyProd100.10yr.S3[9:16,4] = (3*SoyProd100.10yr.S3[1:8,3])+7*SoyProd100.10yr.S3[1,3]
SoyProd100.10yr.S3[17:24,4] = (1*SoyProd100.10yr.S3[1:8,3])+9*SoyProd100.10yr.S3[1,3]
SoyProd100.10yr.S3[25:32,4] = 10*SoyProd100.10yr.S3[1:8,3]

#scenario 4 10-yr calcs
SoyProd100.10yr.S4 = data.frame(matrix(nrow = 32, ncol = 4))
colnames(SoyProd100.10yr.S4) = c("Region","Insecticide","1Year","TenYear")
SoyProd100.10yr.S4$Region = c(rep("North",8),rep("Central",8),rep("South",8),rep("Worst-Case",8))
SoyProd100.10yr.S4$Region = factor(SoyProd100.10yr.S4$Region, levels = c("Worst-Case","North","Central","South"))
SoyProd100.10yr.S4$Insecticide = rep(c(as.character(SC1SoyProdBar$Insecticide),"100%"),4)
SoyProd100.10yr.S4$`1Year` = rep(c(SC4SoyProdBar[1:7,3],0),4)
SoyProd100.10yr.S4[1:8,4] = (5*SoyProd100.10yr.S4[1:8,3])+5*SoyProd100.10yr.S4[1,3]
SoyProd100.10yr.S4[9:16,4] = (3*SoyProd100.10yr.S4[1:8,3])+7*SoyProd100.10yr.S4[1,3]
SoyProd100.10yr.S4[17:24,4] = (1*SoyProd100.10yr.S4[1:8,3])+9*SoyProd100.10yr.S4[1,3]
SoyProd100.10yr.S4[25:32,4] = 10*SoyProd100.10yr.S4[1:8,3]


#North
SoyProd100.10yr.S1[1:8,4]
SoyProd100.10yr.S2[1:8,4]
SoyProd100.10yr.S3[1:8,4]
SoyProd100.10yr.S4[1:8,4]

SoyProd100.10yr.3D.North = matrix(c(SoyProd100.10yr.S1[1:8,4],SoyProd100.10yr.S2[1:8,4],SoyProd100.10yr.S3[1:8,4],SoyProd100.10yr.S4[1:8,4]), ncol = 4)
colnames(SoyProd100.10yr.3D.North) = c("Scenario 1","Scenario 2","Scenario 3","Scenario 4")
rownames(SoyProd100.10yr.3D.North) = c(as.character(SC1SoyProdBar$Insecticide),"100%")
SoyProd100.10yr.3D.North

SoyProd100.10yr.3D.North.R = SoyProd100.10yr.3D.North[c(1,6,5,7,2,4,3,8),]
SoyProd100.10yr.3D.North.R = SoyProd100.10yr.3D.North.R[,c(1,4,3,2)]
SoyProd100.10yr.3D.North.R

hist3D(x = seq(1, nrow(SoyProd100.10yr.3D.North.R), by = 1), y = seq(1, ncol(SoyProd100.10yr.3D.North.R), by = 1), 
       z = SoyProd100.10yr.3D.North.R, expand = 0.6, border = "black", space = 0.2, bty = "g", colkey = FALSE, ticktype = "detailed",
       xlab = "", ylab = "", zlab = "", main = "Adult Monarch Production\nSpray Drift Occurs 5 in 10 Years", zlim = c(0,max(SoyProd100.10yr.3D.North.R)),
       col = col1
       #col = rev(jet.col())
       )
text3D(x = c(0.6,1.6,2.6,3.6,4.6,5.6,6.6,7.6), y = rep(0,8), z = rep(0,8), labels = rownames(SoyProd100.10yr.3D.North.R), adj = 0.9, add = TRUE)
text3D(x = rep(8.3,4), y = c(0.5,1.5,2.5,3.5), z = rep(0,4), labels = colnames(SoyProd100.10yr.3D.North.R), adj = -0.7, add = TRUE)

#Spray.Aer.N

#Central
SoyProd100.10yr.S1[9:16,4]
SoyProd100.10yr.S2[9:16,4]
SoyProd100.10yr.S3[9:16,4]
SoyProd100.10yr.S4[9:16,4]

SoyProd100.10yr.3D.Cent = matrix(c(SoyProd100.10yr.S1[9:16,4],SoyProd100.10yr.S2[9:16,4],SoyProd100.10yr.S3[9:16,4],SoyProd100.10yr.S4[9:16,4]), ncol = 4)
colnames(SoyProd100.10yr.3D.Cent) = c("Scenario 1","Scenario 2","Scenario 3","Scenario 4")
rownames(SoyProd100.10yr.3D.Cent) = c(as.character(SC1SoyProdBar$Insecticide),"100%")
SoyProd100.10yr.3D.Cent

SoyProd100.10yr.3D.Cent.R = SoyProd100.10yr.3D.Cent[c(1,6,5,7,2,4,3,8),]
SoyProd100.10yr.3D.Cent.R = SoyProd100.10yr.3D.Cent.R[,c(1,4,3,2)]
SoyProd100.10yr.3D.Cent.R

hist3D(x = seq(1, nrow(SoyProd100.10yr.3D.Cent.R), by = 1), y = seq(1, ncol(SoyProd100.10yr.3D.Cent.R), by = 1), 
       z = SoyProd100.10yr.3D.Cent.R, expand = 0.6, border = "black", space = 0.2, bty = "g", colkey = FALSE, ticktype = "detailed",
       xlab = "", ylab = "", zlab = "", main = "Adult Monarch Production\nSpray Drift Occurs 3 in 10 Years", zlim = c(0,max(SoyProd100.10yr.3D.Cent.R)),
       #col = rev(jet.col())
       col = col2
       )
text3D(x = c(0.6,1.6,2.6,3.6,4.6,5.6,6.6,7.6), y = rep(0,8), z = rep(0,8), labels = rownames(SoyProd100.10yr.3D.Cent.R), adj = 0.9, add = TRUE)
text3D(x = rep(8.3,4), y = c(0.5,1.5,2.5,3.5), z = rep(0,4), labels = colnames(SoyProd100.10yr.3D.Cent.R), adj = -0.7, add = TRUE)


#South
SoyProd100.10yr.S1[17:24,4]
SoyProd100.10yr.S2[17:24,4]
SoyProd100.10yr.S3[17:24,4]
SoyProd100.10yr.S4[17:24,4]

SoyProd100.10yr.3D.South = matrix(c(SoyProd100.10yr.S1[17:24,4],SoyProd100.10yr.S2[17:24,4],SoyProd100.10yr.S3[17:24,4],SoyProd100.10yr.S4[17:24,4]), ncol = 4)
colnames(SoyProd100.10yr.3D.South) = c("Scenario 1","Scenario 2","Scenario 3","Scenario 4")
rownames(SoyProd100.10yr.3D.South) = c(as.character(SC1SoyProdBar$Insecticide),"100%")
SoyProd100.10yr.3D.South

SoyProd100.10yr.3D.South.R = SoyProd100.10yr.3D.South[c(1,6,5,7,2,4,3,8),]
SoyProd100.10yr.3D.South.R = SoyProd100.10yr.3D.South.R[,c(1,4,3,2)]
SoyProd100.10yr.3D.South.R

hist3D(x = seq(1, nrow(SoyProd100.10yr.3D.South.R), by = 1), y = seq(1, ncol(SoyProd100.10yr.3D.South.R), by = 1), 
       z = SoyProd100.10yr.3D.South.R, expand = 0.6, border = "black", space = 0.2, bty = "g", colkey = FALSE, ticktype = "detailed",
       xlab = "", ylab = "", zlab = "", main = "Adult Monarch Production\nSpray Drift Occurs 1 in 10 Years", zlim = c(0,max(SoyProd100.10yr.3D.South.R)),
       #col = rev(jet.col())
       col = col3
       )
text3D(x = c(0.6,1.6,2.6,3.6,4.6,5.6,6.6,7.6), y = rep(0,8), z = rep(0,8), labels = rownames(SoyProd100.10yr.3D.South.R), adj = 0.9, add = TRUE)
text3D(x = rep(8.3,4), y = c(0.5,1.5,2.5,3.5), z = rep(0,4), labels = colnames(SoyProd100.10yr.3D.South.R), adj = -0.7, add = TRUE)


#Worst-Case Scenario
SoyProd100.10yr.S1[25:32,4]
SoyProd100.10yr.S2[25:32,4]
SoyProd100.10yr.S3[25:32,4]
SoyProd100.10yr.S4[25:32,4]

SoyProd100.10yr.3D.Worst = matrix(c(SoyProd100.10yr.S1[25:32,4],SoyProd100.10yr.S2[25:32,4],SoyProd100.10yr.S3[25:32,4],SoyProd100.10yr.S4[25:32,4]), ncol = 4)
colnames(SoyProd100.10yr.3D.Worst) = c("Scenario 1","Scenario 2","Scenario 3","Scenario 4")
rownames(SoyProd100.10yr.3D.Worst) = c(as.character(SC1SoyProdBar$Insecticide),"100%")
SoyProd100.10yr.3D.Worst

SoyProd100.10yr.3D.Worst.R = SoyProd100.10yr.3D.Worst[c(1,6,5,7,2,4,3,8),]
SoyProd100.10yr.3D.Worst.R = SoyProd100.10yr.3D.Worst.R[,c(1,4,3,2)]
SoyProd100.10yr.3D.Worst.R

hist3D(x = seq(1, nrow(SoyProd100.10yr.3D.Worst.R), by = 1), y = seq(1, ncol(SoyProd100.10yr.3D.Worst.R), by = 1), 
       z = SoyProd100.10yr.3D.Worst.R, expand = 0.6, border = "black", space = 0.2, bty = "g", colkey = FALSE, ticktype = "detailed",
       xlab = "", ylab = "", zlab = "", main = "Adult Monarch Production\nSpray Drift Occurs 10 in 10 Years", zlim = c(0,max(SoyProd100.10yr.3D.Worst.R)),
       col = col4
       #col = rev(jet.col())
       )
text3D(x = c(0.6,1.6,2.6,3.6,4.6,5.6,6.6,7.6), y = rep(0,8), z = rep(0,8), labels = rownames(SoyProd100.10yr.3D.Worst.R), adj = 0.9, add = TRUE)
text3D(x = rep(8.3,4), y = c(0.5,1.5,2.5,3.5), z = rep(0,4), labels = colnames(SoyProd100.10yr.3D.Worst.R), adj = -0.7, add = TRUE)




###########HIGH BOOM###########

#scenario 1 10-yr calcs
SoyProd100HB.10yr.S1 = data.frame(matrix(nrow = 32, ncol = 4))
colnames(SoyProd100HB.10yr.S1) = c("Region","Insecticide","1Year","TenYear")
SoyProd100HB.10yr.S1$Region = c(rep("North",8),rep("Central",8),rep("South",8),rep("Worst-Case",8))
SoyProd100HB.10yr.S1$Region = factor(SoyProd100HB.10yr.S1$Region, levels = c("Worst-Case","North","Central","South"))
SoyProd100HB.10yr.S1$Insecticide = rep(c(as.character(SC1SoyProdBar$Insecticide),"100%"),4)
SoyProd100HB.10yr.S1$`1Year` = rep(c(SC1SoyProdHBBar[1:7,3],0),4)
SoyProd100HB.10yr.S1[1:8,4] = (5*SoyProd100HB.10yr.S1[1:8,3])+5*SoyProd100HB.10yr.S1[1,3]
SoyProd100HB.10yr.S1[9:16,4] = (3*SoyProd100HB.10yr.S1[1:8,3])+7*SoyProd100HB.10yr.S1[1,3]
SoyProd100HB.10yr.S1[17:24,4] = (1*SoyProd100HB.10yr.S1[1:8,3])+9*SoyProd100HB.10yr.S1[1,3]
SoyProd100HB.10yr.S1[25:32,4] = 10*SoyProd100HB.10yr.S1[1:8,3]

#scenario 2 10-yr calcs
SoyProd100HB.10yr.S2 = data.frame(matrix(nrow = 32, ncol = 4))
colnames(SoyProd100HB.10yr.S2) = c("Region","Insecticide","1Year","TenYear")
SoyProd100HB.10yr.S2$Region = c(rep("North",8),rep("Central",8),rep("South",8),rep("Worst-Case",8))
SoyProd100HB.10yr.S2$Region = factor(SoyProd100HB.10yr.S2$Region, levels = c("Worst-Case","North","Central","South"))
SoyProd100HB.10yr.S2$Insecticide = rep(c(as.character(SC1SoyProdBar$Insecticide),"100%"),4)
SoyProd100HB.10yr.S2$`1Year` = rep(c(SC2SoyProdHBBar[1:7,3],0),4)
SoyProd100HB.10yr.S2[1:8,4] = (5*SoyProd100HB.10yr.S2[1:8,3])+5*SoyProd100HB.10yr.S2[1,3]
SoyProd100HB.10yr.S2[9:16,4] = (3*SoyProd100HB.10yr.S2[1:8,3])+7*SoyProd100HB.10yr.S2[1,3]
SoyProd100HB.10yr.S2[17:24,4] = (1*SoyProd100HB.10yr.S2[1:8,3])+9*SoyProd100HB.10yr.S2[1,3]
SoyProd100HB.10yr.S2[25:32,4] = 10*SoyProd100HB.10yr.S2[1:8,3]

#scenario 3 10-yr calcs
SoyProd100HB.10yr.S3 = data.frame(matrix(nrow = 32, ncol = 4))
colnames(SoyProd100HB.10yr.S3) = c("Region","Insecticide","1Year","TenYear")
SoyProd100HB.10yr.S3$Region = c(rep("North",8),rep("Central",8),rep("South",8),rep("Worst-Case",8))
SoyProd100HB.10yr.S3$Region = factor(SoyProd100HB.10yr.S3$Region, levels = c("Worst-Case","North","Central","South"))
SoyProd100HB.10yr.S3$Insecticide = rep(c(as.character(SC1SoyProdBar$Insecticide),"100%"),4)
SoyProd100HB.10yr.S3$`1Year` = rep(c(SC3SoyProdHBBar[1:7,3],0),4)
SoyProd100HB.10yr.S3[1:8,4] = (5*SoyProd100HB.10yr.S3[1:8,3])+5*SoyProd100HB.10yr.S3[1,3]
SoyProd100HB.10yr.S3[9:16,4] = (3*SoyProd100HB.10yr.S3[1:8,3])+7*SoyProd100HB.10yr.S3[1,3]
SoyProd100HB.10yr.S3[17:24,4] = (1*SoyProd100HB.10yr.S3[1:8,3])+9*SoyProd100HB.10yr.S3[1,3]
SoyProd100HB.10yr.S3[25:32,4] = 10*SoyProd100HB.10yr.S3[1:8,3]

#scenario 4 10-yr calcs
SoyProd100HB.10yr.S4 = data.frame(matrix(nrow = 32, ncol = 4))
colnames(SoyProd100HB.10yr.S4) = c("Region","Insecticide","1Year","TenYear")
SoyProd100HB.10yr.S4$Region = c(rep("North",8),rep("Central",8),rep("South",8),rep("Worst-Case",8))
SoyProd100HB.10yr.S4$Region = factor(SoyProd100HB.10yr.S4$Region, levels = c("Worst-Case","North","Central","South"))
SoyProd100HB.10yr.S4$Insecticide = rep(c(as.character(SC1SoyProdBar$Insecticide),"100%"),4)
SoyProd100HB.10yr.S4$`1Year` = rep(c(SC4SoyProdHBBar[1:7,3],0),4)
SoyProd100HB.10yr.S4[1:8,4] = (5*SoyProd100HB.10yr.S4[1:8,3])+5*SoyProd100HB.10yr.S4[1,3]
SoyProd100HB.10yr.S4[9:16,4] = (3*SoyProd100HB.10yr.S4[1:8,3])+7*SoyProd100HB.10yr.S4[1,3]
SoyProd100HB.10yr.S4[17:24,4] = (1*SoyProd100HB.10yr.S4[1:8,3])+9*SoyProd100HB.10yr.S4[1,3]
SoyProd100HB.10yr.S4[25:32,4] = 10*SoyProd100HB.10yr.S4[1:8,3]


#North
SoyProd100HB.10yr.S1[1:8,4]
SoyProd100HB.10yr.S2[1:8,4]
SoyProd100HB.10yr.S3[1:8,4]
SoyProd100HB.10yr.S4[1:8,4]

SoyProd100HB.10yr.3D.North = matrix(c(SoyProd100HB.10yr.S1[1:8,4],SoyProd100HB.10yr.S2[1:8,4],SoyProd100HB.10yr.S3[1:8,4],SoyProd100HB.10yr.S4[1:8,4]), ncol = 4)
colnames(SoyProd100HB.10yr.3D.North) = c("Scenario 1","Scenario 2","Scenario 3","Scenario 4")
rownames(SoyProd100HB.10yr.3D.North) = c(as.character(SC1SoyProdBar$Insecticide),"100%")
SoyProd100HB.10yr.3D.North

SoyProd100HB.10yr.3D.North.R = SoyProd100HB.10yr.3D.North[c(1,6,5,7,2,4,3,8),]
SoyProd100HB.10yr.3D.North.R = SoyProd100HB.10yr.3D.North.R[,c(1,4,3,2)]
SoyProd100HB.10yr.3D.North.R

hist3D(x = seq(1, nrow(SoyProd100HB.10yr.3D.North.R), by = 1), y = seq(1, ncol(SoyProd100HB.10yr.3D.North.R), by = 1), 
       z = SoyProd100HB.10yr.3D.North.R, expand = 0.6, border = "black", space = 0.2, bty = "g", colkey = FALSE, ticktype = "detailed",
       xlab = "", ylab = "", zlab = "", main = "Adult Monarch Production\nSpray Drift Occurs 5 in 10 Years", zlim = c(0,max(SoyProd100HB.10yr.3D.North.R)),
       #col = rev(jet.col())
       col = col5
       )
text3D(x = c(0.6,1.6,2.6,3.6,4.6,5.6,6.6,7.6), y = rep(0,8), z = rep(0,8), labels = rownames(SoyProd100HB.10yr.3D.North.R), adj = 0.9, add = TRUE)
text3D(x = rep(8.3,4), y = c(0.5,1.5,2.5,3.5), z = rep(0,4), labels = colnames(SoyProd100HB.10yr.3D.North.R), adj = -0.7, add = TRUE)

#Spray.HB.N

#Central
SoyProd100HB.10yr.S1[9:16,4]
SoyProd100HB.10yr.S2[9:16,4]
SoyProd100HB.10yr.S3[9:16,4]
SoyProd100HB.10yr.S4[9:16,4]

SoyProd100HB.10yr.3D.Cent = matrix(c(SoyProd100HB.10yr.S1[9:16,4],SoyProd100HB.10yr.S2[9:16,4],SoyProd100HB.10yr.S3[9:16,4],SoyProd100HB.10yr.S4[9:16,4]), ncol = 4)
colnames(SoyProd100HB.10yr.3D.Cent) = c("Scenario 1","Scenario 2","Scenario 3","Scenario 4")
rownames(SoyProd100HB.10yr.3D.Cent) = c(as.character(SC1SoyProdBar$Insecticide),"100%")
SoyProd100HB.10yr.3D.Cent

SoyProd100HB.10yr.3D.Cent.R = SoyProd100HB.10yr.3D.Cent[c(1,6,5,7,2,4,3,8),]
SoyProd100HB.10yr.3D.Cent.R = SoyProd100HB.10yr.3D.Cent.R[,c(1,4,3,2)]
SoyProd100HB.10yr.3D.Cent.R

hist3D(x = seq(1, nrow(SoyProd100HB.10yr.3D.Cent.R), by = 1), y = seq(1, ncol(SoyProd100HB.10yr.3D.Cent.R), by = 1), 
       z = SoyProd100HB.10yr.3D.Cent.R, expand = 0.6, border = "black", space = 0.2, bty = "g", colkey = FALSE, ticktype = "detailed",
       xlab = "", ylab = "", zlab = "", main = "Adult Monarch Production\nSpray Drift Occurs 3 in 10 Years", zlim = c(0,max(SoyProd100HB.10yr.3D.Cent.R)),
       #col = rev(jet.col())
       col = col6
       )
text3D(x = c(0.6,1.6,2.6,3.6,4.6,5.6,6.6,7.6), y = rep(0,8), z = rep(0,8), labels = rownames(SoyProd100HB.10yr.3D.Cent.R), adj = 0.9, add = TRUE)
text3D(x = rep(8.3,4), y = c(0.5,1.5,2.5,3.5), z = rep(0,4), labels = colnames(SoyProd100HB.10yr.3D.Cent.R), adj = -0.7, add = TRUE)


#South
SoyProd100HB.10yr.S1[17:24,4]
SoyProd100HB.10yr.S2[17:24,4]
SoyProd100HB.10yr.S3[17:24,4]
SoyProd100HB.10yr.S4[17:24,4]

SoyProd100HB.10yr.3D.South = matrix(c(SoyProd100HB.10yr.S1[17:24,4],SoyProd100HB.10yr.S2[17:24,4],SoyProd100HB.10yr.S3[17:24,4],SoyProd100HB.10yr.S4[17:24,4]), ncol = 4)
colnames(SoyProd100HB.10yr.3D.South) = c("Scenario 1","Scenario 2","Scenario 3","Scenario 4")
rownames(SoyProd100HB.10yr.3D.South) = c(as.character(SC1SoyProdBar$Insecticide),"100%")
SoyProd100HB.10yr.3D.South

SoyProd100HB.10yr.3D.South.R = SoyProd100HB.10yr.3D.South[c(1,6,5,7,2,4,3,8),]
SoyProd100HB.10yr.3D.South.R = SoyProd100HB.10yr.3D.South.R[,c(1,4,3,2)]
SoyProd100HB.10yr.3D.South.R

hist3D(x = seq(1, nrow(SoyProd100HB.10yr.3D.South.R), by = 1), y = seq(1, ncol(SoyProd100HB.10yr.3D.South.R), by = 1), 
       z = SoyProd100HB.10yr.3D.South.R, expand = 0.6, border = "black", space = 0.2, bty = "g", colkey = FALSE, ticktype = "detailed",
       xlab = "", ylab = "", zlab = "", main = "Adult Monarch Production\nSpray Drift Occurs 1 in 10 Years", zlim = c(0,max(SoyProd100HB.10yr.3D.South.R)),
       #col = rev(jet.col())
       col = col7
       )
text3D(x = c(0.6,1.6,2.6,3.6,4.6,5.6,6.6,7.6), y = rep(0,8), z = rep(0,8), labels = rownames(SoyProd100HB.10yr.3D.South.R), adj = 0.9, add = TRUE)
text3D(x = rep(8.3,4), y = c(0.5,1.5,2.5,3.5), z = rep(0,4), labels = colnames(SoyProd100HB.10yr.3D.South.R), adj = -0.7, add = TRUE)


#Worst-Case Scenario
SoyProd100HB.10yr.S1[25:32,4]
SoyProd100HB.10yr.S2[25:32,4]
SoyProd100HB.10yr.S3[25:32,4]
SoyProd100HB.10yr.S4[25:32,4]

SoyProd100HB.10yr.3D.Worst = matrix(c(SoyProd100HB.10yr.S1[25:32,4],SoyProd100HB.10yr.S2[25:32,4],SoyProd100HB.10yr.S3[25:32,4],SoyProd100HB.10yr.S4[25:32,4]), ncol = 4)
colnames(SoyProd100HB.10yr.3D.Worst) = c("Scenario 1","Scenario 2","Scenario 3","Scenario 4")
rownames(SoyProd100HB.10yr.3D.Worst) = c(as.character(SC1SoyProdBar$Insecticide),"100%")
SoyProd100HB.10yr.3D.Worst

SoyProd100HB.10yr.3D.Worst.R = SoyProd100HB.10yr.3D.Worst[c(1,6,5,7,2,4,3,8),]
SoyProd100HB.10yr.3D.Worst.R = SoyProd100HB.10yr.3D.Worst.R[,c(1,4,3,2)]
SoyProd100HB.10yr.3D.Worst.R

hist3D(x = seq(1, nrow(SoyProd100HB.10yr.3D.Worst.R), by = 1), y = seq(1, ncol(SoyProd100HB.10yr.3D.Worst.R), by = 1), 
       z = SoyProd100HB.10yr.3D.Worst.R, expand = 0.6, border = "black", space = 0.2, bty = "g", colkey = FALSE, ticktype = "detailed",
       xlab = "", ylab = "", zlab = "", main = "Adult Monarch Production\nSpray Drift Occurs 10 in 10 Years", zlim = c(0,max(SoyProd100HB.10yr.3D.Worst.R)),
       #col = rev(jet.col())
       col = col8
       )
text3D(x = c(0.6,1.6,2.6,3.6,4.6,5.6,6.6,7.6), y = rep(0,8), z = rep(0,8), labels = rownames(SoyProd100HB.10yr.3D.Worst.R), adj = 0.9, add = TRUE)
text3D(x = rep(8.3,4), y = c(0.5,1.5,2.5,3.5), z = rep(0,4), labels = colnames(SoyProd100HB.10yr.3D.Worst.R), adj = -0.7, add = TRUE)




############## 3D Graphs of full Pop, rather than just spray drift zone pop ##########################


#colors
#range divided by max = proportion of colors should be used
(max(SoyProd100.10yr.3D.North.R.TOT)-min(SoyProd100.10yr.3D.North.R.TOT))/(max(SoyProd100.10yr.3D.Worst.R.TOT)-min(SoyProd100.10yr.3D.Worst.R.TOT)) #range = 101,922, prop = 0.88
(max(SoyProd100.10yr.3D.Cent.R.TOT)-min(SoyProd100.10yr.3D.Cent.R.TOT))/(max(SoyProd100.10yr.3D.Worst.R.TOT)-min(SoyProd100.10yr.3D.Worst.R.TOT)) #range = 96,410, prop = 0.83
(max(SoyProd100.10yr.3D.South.R.TOT)-min(SoyProd100.10yr.3D.South.R.TOT))/(max(SoyProd100.10yr.3D.Worst.R.TOT)-min(SoyProd100.10yr.3D.Worst.R.TOT)) #range = 90,898, prop = 0.79
(max(SoyProd100.10yr.3D.Worst.R.TOT)-min(SoyProd100.10yr.3D.Worst.R.TOT))/(max(SoyProd100.10yr.3D.Worst.R.TOT)-min(SoyProd100.10yr.3D.Worst.R.TOT)) #range = 115,702, prop = 1

(max(SoyProd100HB.10yr.3D.North.R.TOT)-min(SoyProd100HB.10yr.3D.North.R.TOT))/(max(SoyProd100HB.10yr.3D.Worst.R.TOT)-min(SoyProd100HB.10yr.3D.Worst.R.TOT)) #range = 101,922, prop = 0.88
(max(SoyProd100HB.10yr.3D.Cent.R.TOT)-min(SoyProd100HB.10yr.3D.Cent.R.TOT))/(max(SoyProd100HB.10yr.3D.Worst.R.TOT)-min(SoyProd100HB.10yr.3D.Worst.R.TOT)) #range = 96,410, prop = 0.83
(max(SoyProd100HB.10yr.3D.South.R.TOT)-min(SoyProd100HB.10yr.3D.South.R.TOT))/(max(SoyProd100HB.10yr.3D.Worst.R.TOT)-min(SoyProd100HB.10yr.3D.Worst.R.TOT)) #range = 90,898, prop = 0.79
(max(SoyProd100HB.10yr.3D.Worst.R.TOT)-min(SoyProd100HB.10yr.3D.Worst.R.TOT))/(max(SoyProd100HB.10yr.3D.Worst.R.TOT)-min(SoyProd100HB.10yr.3D.Worst.R.TOT)) #range = 115,702, prop = 1

#base color ramp for maximum range
col.base2 = rev(jet.col(100))

#"col = " for first figure below - aerial, north
col9 = col.base2[(100-88):100]
col10 = col.base2[(100-83):100]
col11 = col.base2[(100-79):100]
col12 = col.base2[1:100] 
col13 = col.base2[(100-88):100]
col14 = col.base2[(100-83):100]
col15 = col.base2[(100-79):100]
col16 = col.base2[1:100]


#10 yr aerial


#North
SoyProd100.10yr.3D.North.R.TOT = SoyProd100.10yr.3D.North.R
SoyProd100.10yr.3D.North.R.TOT[,1] = SoyProd100.10yr.3D.North.R[,1]+(SC1SoyProd[1,1]*10)
SoyProd100.10yr.3D.North.R.TOT[,2] = SoyProd100.10yr.3D.North.R[,2]+(SC4SoyProd[1,1]*10)
SoyProd100.10yr.3D.North.R.TOT[,3] = SoyProd100.10yr.3D.North.R[,3]+(SC3SoyProd[1,1]*10)
SoyProd100.10yr.3D.North.R.TOT[,4] = SoyProd100.10yr.3D.North.R[,4]+(SC2SoyProd[1,1]*10)

hist3D(x = seq(1, nrow(SoyProd100.10yr.3D.North.R.TOT), by = 1), y = seq(1, ncol(SoyProd100.10yr.3D.North.R.TOT), by = 1), 
       z = SoyProd100.10yr.3D.North.R.TOT, expand = 0.6, border = "black", space = 0.2, bty = "g", colkey = FALSE, ticktype = "detailed",
       xlab = "", ylab = "", zlab = "", main = "Adult Monarch Production\nSpray Drift Occurs 5 in 10 Years", zlim = c(0,max(SoyProd100.10yr.3D.North.R.TOT)),
       col = col9)
text3D(x = c(0.6,1.6,2.6,3.6,4.6,5.6,6.6,7.6), y = rep(0,8), z = rep(0,8), labels = rownames(SoyProd100.10yr.3D.North.R.TOT), adj = 0.9, add = TRUE)
text3D(x = rep(8.3,4), y = c(0.5,1.5,2.5,3.5), z = rep(0,4), labels = colnames(SoyProd100.10yr.3D.North.R.TOT), adj = -0.7, add = TRUE)

#Full.Aer.N

#Central
SoyProd100.10yr.3D.Cent.R.TOT = SoyProd100.10yr.3D.Cent.R
SoyProd100.10yr.3D.Cent.R.TOT[,1] = SoyProd100.10yr.3D.Cent.R[,1]+(SC1SoyProd[1,1]*10)
SoyProd100.10yr.3D.Cent.R.TOT[,2] = SoyProd100.10yr.3D.Cent.R[,2]+(SC4SoyProd[1,1]*10)
SoyProd100.10yr.3D.Cent.R.TOT[,3] = SoyProd100.10yr.3D.Cent.R[,3]+(SC3SoyProd[1,1]*10)
SoyProd100.10yr.3D.Cent.R.TOT[,4] = SoyProd100.10yr.3D.Cent.R[,4]+(SC2SoyProd[1,1]*10)

hist3D(x = seq(1, nrow(SoyProd100.10yr.3D.Cent.R.TOT), by = 1), y = seq(1, ncol(SoyProd100.10yr.3D.Cent.R.TOT), by = 1), 
       z = SoyProd100.10yr.3D.Cent.R.TOT, expand = 0.6, border = "black", space = 0.2, bty = "g", colkey = FALSE, ticktype = "detailed",
       xlab = "", ylab = "", zlab = "", main = "Adult Monarch Production\nSpray Drift Occurs 3 in 10 Years", zlim = c(0,max(SoyProd100.10yr.3D.Cent.R.TOT)),
       col = col10)
text3D(x = c(0.6,1.6,2.6,3.6,4.6,5.6,6.6,7.6), y = rep(0,8), z = rep(0,8), labels = rownames(SoyProd100.10yr.3D.Cent.R.TOT), adj = 0.9, add = TRUE)
text3D(x = rep(8.3,4), y = c(0.5,1.5,2.5,3.5), z = rep(0,4), labels = colnames(SoyProd100.10yr.3D.Cent.R.TOT), adj = -0.7, add = TRUE)


#South
SoyProd100.10yr.3D.South.R.TOT = SoyProd100.10yr.3D.South.R
SoyProd100.10yr.3D.South.R.TOT[,1] = SoyProd100.10yr.3D.South.R[,1]+(SC1SoyProd[1,1]*10)
SoyProd100.10yr.3D.South.R.TOT[,2] = SoyProd100.10yr.3D.South.R[,2]+(SC4SoyProd[1,1]*10)
SoyProd100.10yr.3D.South.R.TOT[,3] = SoyProd100.10yr.3D.South.R[,3]+(SC3SoyProd[1,1]*10)
SoyProd100.10yr.3D.South.R.TOT[,4] = SoyProd100.10yr.3D.South.R[,4]+(SC2SoyProd[1,1]*10)

hist3D(x = seq(1, nrow(SoyProd100.10yr.3D.South.R.TOT), by = 1), y = seq(1, ncol(SoyProd100.10yr.3D.South.R.TOT), by = 1), 
       z = SoyProd100.10yr.3D.South.R.TOT, expand = 0.6, border = "black", space = 0.2, bty = "g", colkey = FALSE, ticktype = "detailed",
       xlab = "", ylab = "", zlab = "", main = "Adult Monarch Production\nSpray Drift Occurs 1 in 10 Years", zlim = c(0,max(SoyProd100.10yr.3D.South.R.TOT)),
       col = col11)
text3D(x = c(0.6,1.6,2.6,3.6,4.6,5.6,6.6,7.6), y = rep(0,8), z = rep(0,8), labels = rownames(SoyProd100.10yr.3D.South.R.TOT), adj = 0.9, add = TRUE)
text3D(x = rep(8.3,4), y = c(0.5,1.5,2.5,3.5), z = rep(0,4), labels = colnames(SoyProd100.10yr.3D.South.R.TOT), adj = -0.7, add = TRUE)


#worst case scenario
SoyProd100.10yr.3D.Worst.R.TOT = SoyProd100.10yr.3D.Worst.R
SoyProd100.10yr.3D.Worst.R.TOT[,1] = SoyProd100.10yr.3D.Worst.R[,1]+(SC1SoyProd[1,1]*10)
SoyProd100.10yr.3D.Worst.R.TOT[,2] = SoyProd100.10yr.3D.Worst.R[,2]+(SC4SoyProd[1,1]*10)
SoyProd100.10yr.3D.Worst.R.TOT[,3] = SoyProd100.10yr.3D.Worst.R[,3]+(SC3SoyProd[1,1]*10)
SoyProd100.10yr.3D.Worst.R.TOT[,4] = SoyProd100.10yr.3D.Worst.R[,4]+(SC2SoyProd[1,1]*10)

hist3D(x = seq(1, nrow(SoyProd100.10yr.3D.Worst.R.TOT), by = 1), y = seq(1, ncol(SoyProd100.10yr.3D.Worst.R.TOT), by = 1), 
       z = SoyProd100.10yr.3D.Worst.R.TOT, expand = 0.6, border = "black", space = 0.2, bty = "g", colkey = FALSE, ticktype = "detailed",
       xlab = "", ylab = "", zlab = "", main = "Adult Monarch Production\nSpray Drift Occurs 10 in 10 Years", zlim = c(0,max(SoyProd100.10yr.3D.Worst.R.TOT)),
       col = col12)
text3D(x = c(0.6,1.6,2.6,3.6,4.6,5.6,6.6,7.6), y = rep(0,8), z = rep(0,8), labels = rownames(SoyProd100.10yr.3D.Worst.R.TOT), adj = 0.9, add = TRUE)
text3D(x = rep(8.3,4), y = c(0.5,1.5,2.5,3.5), z = rep(0,4), labels = colnames(SoyProd100.10yr.3D.Worst.R.TOT), adj = -0.7, add = TRUE)


#10 year high boom


#North
SoyProd100HB.10yr.3D.North.R.TOT = SoyProd100HB.10yr.3D.North.R
SoyProd100HB.10yr.3D.North.R.TOT[,1] = SoyProd100HB.10yr.3D.North.R[,1]+(SC1SoyProdHB[1,1]*10)
SoyProd100HB.10yr.3D.North.R.TOT[,2] = SoyProd100HB.10yr.3D.North.R[,2]+(SC4SoyProdHB[1,1]*10)
SoyProd100HB.10yr.3D.North.R.TOT[,3] = SoyProd100HB.10yr.3D.North.R[,3]+(SC3SoyProdHB[1,1]*10)
SoyProd100HB.10yr.3D.North.R.TOT[,4] = SoyProd100HB.10yr.3D.North.R[,4]+(SC2SoyProdHB[1,1]*10)

hist3D(x = seq(1, nrow(SoyProd100HB.10yr.3D.North.R.TOT), by = 1), y = seq(1, ncol(SoyProd100HB.10yr.3D.North.R.TOT), by = 1), 
       z = SoyProd100HB.10yr.3D.North.R.TOT, expand = 0.6, border = "black", space = 0.2, bty = "g", colkey = FALSE, ticktype = "detailed",
       xlab = "", ylab = "", zlab = "", main = "Adult Monarch Production\nSpray Drift Occurs 5 in 10 Years", zlim = c(0,max(SoyProd100HB.10yr.3D.North.R.TOT)),
       col = col13)
text3D(x = c(0.6,1.6,2.6,3.6,4.6,5.6,6.6,7.6), y = rep(0,8), z = rep(0,8), labels = rownames(SoyProd100HB.10yr.3D.North.R.TOT), adj = 0.9, add = TRUE)
text3D(x = rep(8.3,4), y = c(0.5,1.5,2.5,3.5), z = rep(0,4), labels = colnames(SoyProd100HB.10yr.3D.North.R.TOT), adj = -0.7, add = TRUE)

#Full.HB.N

#Central
SoyProd100HB.10yr.3D.Cent.R.TOT = SoyProd100HB.10yr.3D.Cent.R
SoyProd100HB.10yr.3D.Cent.R.TOT[,1] = SoyProd100HB.10yr.3D.Cent.R[,1]+(SC1SoyProdHB[1,1]*10)
SoyProd100HB.10yr.3D.Cent.R.TOT[,2] = SoyProd100HB.10yr.3D.Cent.R[,2]+(SC4SoyProdHB[1,1]*10)
SoyProd100HB.10yr.3D.Cent.R.TOT[,3] = SoyProd100HB.10yr.3D.Cent.R[,3]+(SC3SoyProdHB[1,1]*10)
SoyProd100HB.10yr.3D.Cent.R.TOT[,4] = SoyProd100HB.10yr.3D.Cent.R[,4]+(SC2SoyProdHB[1,1]*10)

hist3D(x = seq(1, nrow(SoyProd100HB.10yr.3D.Cent.R.TOT), by = 1), y = seq(1, ncol(SoyProd100HB.10yr.3D.Cent.R.TOT), by = 1), 
       z = SoyProd100HB.10yr.3D.Cent.R.TOT, expand = 0.6, border = "black", space = 0.2, bty = "g", colkey = FALSE, ticktype = "detailed",
       xlab = "", ylab = "", zlab = "", main = "Adult Monarch Production\nSpray Drift Occurs 3 in 10 Years", zlim = c(0,max(SoyProd100HB.10yr.3D.Cent.R.TOT)),
       col = col14)
text3D(x = c(0.6,1.6,2.6,3.6,4.6,5.6,6.6,7.6), y = rep(0,8), z = rep(0,8), labels = rownames(SoyProd100HB.10yr.3D.Cent.R.TOT), adj = 0.9, add = TRUE)
text3D(x = rep(8.3,4), y = c(0.5,1.5,2.5,3.5), z = rep(0,4), labels = colnames(SoyProd100HB.10yr.3D.Cent.R.TOT), adj = -0.7, add = TRUE)


#South
SoyProd100HB.10yr.3D.South.R.TOT = SoyProd100HB.10yr.3D.South.R
SoyProd100HB.10yr.3D.South.R.TOT[,1] = SoyProd100HB.10yr.3D.South.R[,1]+(SC1SoyProdHB[1,1]*10)
SoyProd100HB.10yr.3D.South.R.TOT[,2] = SoyProd100HB.10yr.3D.South.R[,2]+(SC4SoyProdHB[1,1]*10)
SoyProd100HB.10yr.3D.South.R.TOT[,3] = SoyProd100HB.10yr.3D.South.R[,3]+(SC3SoyProdHB[1,1]*10)
SoyProd100HB.10yr.3D.South.R.TOT[,4] = SoyProd100HB.10yr.3D.South.R[,4]+(SC2SoyProdHB[1,1]*10)

hist3D(x = seq(1, nrow(SoyProd100HB.10yr.3D.South.R.TOT), by = 1), y = seq(1, ncol(SoyProd100HB.10yr.3D.South.R.TOT), by = 1), 
       z = SoyProd100HB.10yr.3D.South.R.TOT, expand = 0.6, border = "black", space = 0.2, bty = "g", colkey = FALSE, ticktype = "detailed",
       xlab = "", ylab = "", zlab = "", main = "Adult Monarch Production\nSpray Drift Occurs 1 in 10 Years", zlim = c(0,max(SoyProd100HB.10yr.3D.South.R.TOT)),
       col = col15)
text3D(x = c(0.6,1.6,2.6,3.6,4.6,5.6,6.6,7.6), y = rep(0,8), z = rep(0,8), labels = rownames(SoyProd100HB.10yr.3D.South.R.TOT), adj = 0.9, add = TRUE)
text3D(x = rep(8.3,4), y = c(0.5,1.5,2.5,3.5), z = rep(0,4), labels = colnames(SoyProd100HB.10yr.3D.South.R.TOT), adj = -0.7, add = TRUE)


#Worst-case Scenario
SoyProd100HB.10yr.3D.Worst.R.TOT = SoyProd100HB.10yr.3D.Worst.R
SoyProd100HB.10yr.3D.Worst.R.TOT[,1] = SoyProd100HB.10yr.3D.Worst.R[,1]+(SC1SoyProdHB[1,1]*10)
SoyProd100HB.10yr.3D.Worst.R.TOT[,2] = SoyProd100HB.10yr.3D.Worst.R[,2]+(SC4SoyProdHB[1,1]*10)
SoyProd100HB.10yr.3D.Worst.R.TOT[,3] = SoyProd100HB.10yr.3D.Worst.R[,3]+(SC3SoyProdHB[1,1]*10)
SoyProd100HB.10yr.3D.Worst.R.TOT[,4] = SoyProd100HB.10yr.3D.Worst.R[,4]+(SC2SoyProdHB[1,1]*10)

hist3D(x = seq(1, nrow(SoyProd100HB.10yr.3D.Worst.R.TOT), by = 1), y = seq(1, ncol(SoyProd100HB.10yr.3D.Worst.R.TOT), by = 1), 
       z = SoyProd100HB.10yr.3D.Worst.R.TOT, expand = 0.6, border = "black", space = 0.2, bty = "g", colkey = FALSE, ticktype = "detailed",
       xlab = "", ylab = "", zlab = "", main = "Adult Monarch Production\nSpray Drift Occurs 10 in 10 Years", zlim = c(0,max(SoyProd100HB.10yr.3D.Worst.R.TOT)),
       col = col16)
text3D(x = c(0.6,1.6,2.6,3.6,4.6,5.6,6.6,7.6), y = rep(0,8), z = rep(0,8), labels = rownames(SoyProd100HB.10yr.3D.Worst.R.TOT), adj = 0.9, add = TRUE)
text3D(x = rep(8.3,4), y = c(0.5,1.5,2.5,3.5), z = rep(0,4), labels = colnames(SoyProd100HB.10yr.3D.Worst.R.TOT), adj = -0.7, add = TRUE)

write.csv(SoyProd100.10yr.3D.Worst.R.TOT, file = "SoyProd100.10yr.3D.Worst.R.TOT.csv")

save(SoyProd100.10yr.3D.North.R.TOT, SoyProd100.10yr.3D.Cent.R.TOT, SoyProd100.10yr.3D.South.R.TOT, SoyProd100.10yr.3D.Worst.R.TOT,
     SoyProd100HB.10yr.3D.North.R.TOT, SoyProd100HB.10yr.3D.Cent.R.TOT, SoyProd100HB.10yr.3D.South.R.TOT, SoyProd100HB.10yr.3D.Worst.R.TOT,
     file = "FullMonProd.RData")





#########  Break-Even Milkweed Vs. Pesticide Mortality ##################
# if pesticide mortality is 100%, how does remove all gains from milkweed planting?
#eggs and adult production without drift, or, outside eggs/adults = 100% mortality in drift zone
#all single year prod numbers

BreakEven = data.frame(matrix(nrow = 37, ncol = 4))
colnames(BreakEven) = c("Scenario 1","Scenario 2","Scenario 3","Scenario 4")
rownames(BreakEven) = c("Outside Drift Zone Eggs","Inside Drift Zone Eggs",
                        "Total Eggs","Outside Drift Zone Prod","Inside Drift Zone Prod","Total Production","Natural Survival Rate",
                        "BCF Prod","CTR Prod","CFS Prod","IMI Prod","TMX Prod","CDN Prod",
                        "BCF Survival","CTR Survival","CFS Survival","IMI Survival","TMX Survival","CDN Survival",
                        "BCF Survival Delta","CTR Survival Delta","CFS Survival Delta","IMI Survival Delta","TMX Survival Delta","CDN Survival Delta",
                        "BCF 4-3 Surv","CTR 4-3 Surv","CFS 4-3 Surv","IMI 4-3 Surv","TMX 4-3 Prod","CDN 4-3 Surv",
                        "BCF 4-3 Delta","CTR 4-3 Delta","CFS 4-3 Delta","IMI 4-3 Delta","TMX 4-3 Delta","CDN 4-3 Delta"
                        )

BreakEven[1,1] = sum(BaselineNWDBF[which(BaselineNWDBF[,18]==" "),21]) #1,989,715
BreakEven[2,1] = sum(BaselineNWDBF[which(BaselineNWDBF[,18]=="Drift"),21]) #186,639.1
BreakEven[3,1] = sum(BaselineNWDBF$PropEggs)

BreakEven[1,2] = sum(MaxAugNWDBF[which(MaxAugNWDBF[,19]==" "),22]) #2,420,443
BreakEven[2,2] = sum(MaxAugNWDBF[which(MaxAugNWDBF[,19]=="Drift"),22]) #292,971.4
BreakEven[3,2] = sum(MaxAugNWDBF$PropEggs)

BreakEven[1,3] = sum(MedAugNWDBF[which(MedAugNWDBF[,18]==" "),21]) #2,147,292
BreakEven[2,3] = sum(MedAugNWDBF[which(MedAugNWDBF[,18]=="Drift"),21]) #232,002.4
BreakEven[3,3] = sum(MedAugNWDBF$PropEggs)

BreakEven[1,4] = sum(AugOutNWDBF[which(AugOutNWDBF[,5]==" "),21]) #2,080,773
BreakEven[2,4] = sum(AugOutNWDBF[which(AugOutNWDBF[,5]=="Drift"),21]) #172,365.1
BreakEven[3,4] = sum(AugOutNWDBF$PropEggs)

BreakEven[4:6,1] = SC1SoyProd[1:3,1]
BreakEven[4:6,2] = SC2SoyProd[1:3,1]
BreakEven[4:6,3] = SC3SoyProd[1:3,1]
BreakEven[4:6,4] = SC4SoyProd[1:3,1]

#natural morality rate - no pesticide drift
BreakEven[7,] = BreakEven[5,]/BreakEven[2,]

#bring in soybean aphid aerial spray application adult production numbers
BreakEven[8:13,1] = SC1SoyProdBar[2:7,3]
BreakEven[8:13,2] = SC2SoyProdBar[2:7,3]
BreakEven[8:13,3] = SC3SoyProdBar[2:7,3]
BreakEven[8:13,4] = SC4SoyProdBar[2:7,3]

#survival rate when spray drift occurs
BreakEven[14:19,1] = BreakEven[8:13,1]/BreakEven[2,1]
BreakEven[14:19,2] = BreakEven[8:13,2]/BreakEven[2,2]
BreakEven[14:19,3] = BreakEven[8:13,3]/BreakEven[2,3]
BreakEven[14:19,4] = BreakEven[8:13,4]/BreakEven[2,4]

#mortality/survival due to the insecticide - natural subtracted out
BreakEven[20:25,1] = BreakEven[7,1] - BreakEven[8:13,1]/BreakEven[2,1]
BreakEven[20:25,2] = BreakEven[7,2] - BreakEven[8:13,2]/BreakEven[2,2]
BreakEven[20:25,3] = BreakEven[7,3] - BreakEven[8:13,3]/BreakEven[2,3]
BreakEven[20:25,4] = BreakEven[7,4] - BreakEven[8:13,4]/BreakEven[2,4]

#overall survival in scenario 3 given scenario 4 prod numbers
BreakEven[26,3] = BreakEven[8,4]/BreakEven[2,3]
BreakEven[26:31,3] = BreakEven[8:13,4]/BreakEven[2,3]

#survival delta for scenario 3 that would create production numbers like scenario 4
#described in another way, the survival rate for scenario 3 that would create production numbers like scenario 4
BreakEven[32,3] = (BreakEven[8,3]-BreakEven[8,4])/BreakEven[2,3]
BreakEven[32:37,3] = (BreakEven[8:13,3]-BreakEven[8:13,4])/BreakEven[2,3]

#check - should be difference between scen 3 and scen 4 BCF prod = 238
BreakEven[32,3]*BreakEven[2,3]
BreakEven[32:37,3]*BreakEven[2,3]

write.csv(BreakEven, "Breakeven.csv")




