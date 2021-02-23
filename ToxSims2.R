#Analysis for Landscape-scale insecticide risk to monarch butterflies

library(dplyr)


### read in RS output  ###########

scen1EZ = read.csv("CumEggsPerZone.2019.Oct.24.09_52_48.txt")
sum(scen1EZ$Eggs) #2,176,354
scen2EZ = read.csv("CumEggsPerZone.2019.Oct.24.10_33_58.txt")
sum(scen2EZ$Eggs) #2,713,414
scen3EZ = read.csv("CumEggsPerZone.2019.Oct.24.11_16_49.txt")
sum(scen3EZ$Eggs) #2,379,294
scen4EZ = read.csv("CumEggsPerZone.2019.Oct.24.14_39_14.txt")
sum(scen4EZ$Eggs) #2,253,138

length(scen1EZ$run)/50 #41,326, the number of polgyons in the shapefile


#combine instances - run code below for each map
Denresults = data.frame(matrix(nrow=41326, ncol=5))
colnames(Denresults) = c("PolygonID", "HabType", "CumEggs","PolygonArea","ProbEggs")

dens = scen1EZ
dens = scen2EZ
dens = scen3EZ
dens = scen4EZ
nrow(dens) #should be 41,326*50=2,066,300

#loop through dens to create object with densities for each polygon - 33-45 min
system.time(
  for(i in 1:41326)
  {
    densi = filter(dens, ID == i) #combine all 20 instances for each polygon
    Denresults[i,1] = densi[1,3] #polygon ID number
    Denresults[i,2] = as.character(densi[1,4]) #habitat type
    Denresults[i,3] = sum(densi$Eggs) #sum across the 20 instances to get total eggs for per polygon
    Denresults[i,4] = densi[1,6] #polygon area - sometimes lat/long, depends on shapefile
    Denresults[i,5] = densi[1,5] #probEggs
  }
)

Denresults1 = Denresults
sum(Denresults1$CumEggs) #2,176,354, same as above
filter(dens, ID == 12) #just seeing how consistent results are for 1 poly - not too consistent because only a few butterflies per instance
Denresults2 = Denresults
sum(Denresults2$CumEggs) #2,713,414, same as above
Denresults3 = Denresults
sum(Denresults3$CumEggs) #2,379,294, same as above
Denresults4 = Denresults
sum(Denresults4$CumEggs) #2,253,138, same as above

#calculate eggs density as eggs/ha
Denresults1$Density = (Denresults1$CumEggs/Denresults1$PolygonArea)*10000
Denresults2$Density = round((Denresults2$CumEggs/Denresults2$PolygonArea)*10000,5)
Denresults3$Density = (Denresults3$CumEggs/Denresults3$PolygonArea)*10000
Denresults4$Density = (Denresults4$CumEggs/Denresults4$PolygonArea)*10000

#set PolygonID to 0-41325
Denresults1$PolygonID = 0:41325
Denresults2$PolygonID = 0:41325
Denresults3$PolygonID = 0:41325
Denresults4$PolygonID = 0:41325

#save a .csv to join in ArcMap for making maps
write.csv(Denresults1, "Scen1EggDens.csv")
write.csv(Denresults2, "Scen2EggDens.csv")
write.csv(Denresults3, "Scen3EggDens.csv")
write.csv(Denresults4, "Scen4EggDens.csv")

#calculate eggs laid in different habitat types
Eggs = data.frame(matrix(nrow=17, ncol=5))
colnames(Eggs) = c("HabType","Scen1", "Scen2", "Scen3","Scen4")
Eggs$HabType = unique(Denresults1$HabType)

for(i in Eggs$HabType){
  eggsi = filter(Denresults1, HabType == i) #combine all polygons for each habitat type
  Eggs[which(Eggs$HabType==i),2] = sum(eggsi$CumEggs) #polygon ID number
}

for(i in Eggs$HabType){
  eggsi = filter(Denresults2, HabType == i) #combine all polygons for each habitat type
  Eggs[which(Eggs$HabType==i),3] = sum(eggsi$CumEggs) #polygon ID number
}

for(i in Eggs$HabType){
  eggsi = filter(Denresults3, HabType == i) #combine all polygons for each habitat type
  Eggs[which(Eggs$HabType==i),4] = sum(eggsi$CumEggs) #polygon ID number
}

for(i in Eggs$HabType){
  eggsi = filter(Denresults4, HabType == i) #combine all polygons for each habitat type
  Eggs[which(Eggs$HabType==i),5] = sum(eggsi$CumEggs) #polygon ID number
}

Eggs[18,1] = "Total"
Eggs[18,2] = sum(Eggs[1:17,2])
Eggs[18,3] = sum(Eggs[1:17,3])
Eggs[18,4] = sum(Eggs[1:17,4])
Eggs[18,5] = sum(Eggs[1:17,5])

write.csv(Eggs, "EggsLaid.csv")


#get dummy variable for buffer zones
MedAugTable = read.csv("Export_MedAug.txt")
table(MedAugTable$CLASS2)
for(i in 1:41326){
  if(MedAugTable[i,6] != "Buffer"){
    MedAugTable[i,6] = "Outside"
  }
}
#well, that turned all the blank cells into NA, which isn't so bad I guess
#need to define "Outside" as one of the factors, or it won't use it
MedAugTable$CLASS2 = factor(MedAugTable$CLASS2, levels = c("Buffer","Outside",NA))
for(i in 1:41326){
  if(is.na(MedAugTable[i,6])){
    MedAugTable[i,6] = "Outside"
  }
}
table(MedAugTable$CLASS2)

Denresults1$Buffer = MedAugTable$CLASS2
Denresults2$Buffer = MedAugTable$CLASS2
Denresults3$Buffer = MedAugTable$CLASS2
Denresults4$Buffer = MedAugTable$CLASS2


#eggs within and outside buffer zone
BufferEggs = data.frame(matrix(nrow=2, ncol=5))
colnames(BufferEggs) = c("Buffer","Scen1", "Scen2", "Scen3","Scen4")
BufferEggs[1,1] = "Total-Inside"; BufferEggs[2,1] = "Total-Outside"

#total eggs inside and outside buffer regardless of landcover type - for Denresults1,2,3, or 4
densB = filter(Denresults4, Buffer == "Buffer") 
BufferEggs[1,5] = sum(densB$CumEggs)

densC = filter(Denresults4, Buffer == "Outside") 
BufferEggs[2,5] = sum(densC$CumEggs)

sum(BufferEggs$Scen1) #2,176,354
sum(BufferEggs$Scen2) #2,713,414
sum(BufferEggs$Scen3) #2,379,294
sum(BufferEggs$Scen4) #2,253,138

#calculate roadsides, low intensity dev, and grass/pasture inside/outside
BufferEggs[3,1] = "G/P-In"; BufferEggs[4,1] = "G/P-Out"

densB = filter(Denresults4, Buffer == "Buffer" & HabType == "Grass/Pasture") 
BufferEggs[3,5] = sum(densB$CumEggs)

densC = filter(Denresults4, Buffer == "Outside" & HabType == "Grass/Pasture") 
BufferEggs[4,5] = sum(densC$CumEggs)

BufferEggs[5,1] = "Roadsides-In"; BufferEggs[6,1] = "Roadsides-Out"

densB = filter(Denresults4, Buffer == "Buffer" & (HabType == "MWROW0" | HabType == "MWROW1-5" | HabType == "MWROW5-20" | HabType == "MWROW20-60" | HabType == "MWROW60-100"))
BufferEggs[5,5] = sum(densB$CumEggs)

densC = filter(Denresults4, Buffer == "Outside" & (HabType == "MWROW0" | HabType == "MWROW1-5" | HabType == "MWROW5-20" | HabType == "MWROW20-60" | HabType == "MWROW60-100"))
BufferEggs[6,5] = sum(densC$CumEggs)

BufferEggs[7,1] = "LowInt-In"; BufferEggs[8,1] = "LowInt-Out"

densB = filter(Denresults4, Buffer == "Buffer" & HabType == "Developed/Open Space/Developed/Low Intensity") 
BufferEggs[7,5] = sum(densB$CumEggs)

densC = filter(Denresults4, Buffer == "Outside" & HabType == "Developed/Open Space/Developed/Low Intensity") 
BufferEggs[8,5] = sum(densC$CumEggs)

write.csv(BufferEggs, "BufferEggsLaid.csv")



















