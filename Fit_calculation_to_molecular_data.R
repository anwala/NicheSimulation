#!/usr/bin/env Rscript
##Selective Mutation Accumulation: A Compuational Model of the Paternal Age Effect
##Eoin C Whelan*, Alexander C Nwala, Christopher Osgood, Stephan Olariu
##Old Dominion University Norfolk VA USA 
##*Contact: ewhel001@odu.edu

##Supplementary codes
##R script for generating mutation curves and matching them to sequencing data.

#Calculates number of mutations per million cells

library(expm) #required for exponential multiplication of matricies

options(scipen=999)#removes scientific notation as an output

n<-324                       #number of cells per niche
states<-n+1                  #number of rows/columns in the matrix (i.e. mutation states including 0)
maxAge <- 80
divisionsEveryNumberOfDays <- 16
divisionsPerCellPerYear<-365.25/divisionsEveryNumberOfDays

mutationsPerNumber<- 1000000  #this is the metric of mutants per X number of cells

maxN <- c(52000, 52000, 52000, 52000) #default maximum lifetime number of cell divisions per niche

p.list <- c(2*10^-10, 1*10^-10, 6*10^-10, 2*10^-10) #mutation rate including adjustment for CpG status

numberOfDataPoints = c(7, 7, 8, 9) #number of age categories, based on the data used.

basic.output.data <- data.frame(Age=1:numberOfDataPoints[4], MutantsProportion=1:numberOfDataPoints[4], MutantsPerMillion=1:numberOfDataPoints[4])
basic.output.data[] <- 0

FGFR2_C755G <- c(3, 12, 15, 23, 37, 19, 53, 0, 0)  
FGFR2_C758G <- c(4, 5, 11, 15, 14, 14.5, 21, 0, 0)
HRAS_G34A <- c(10, 7, 13.5, 15, 22, 22.5, 27.5, 39, 0) 
FGFR3_A1948G <- c(10.87, 18.97, 18.32, 14.26, 26.16, 31.11, 28.38, 58.45, 38.20) 
#output needs to have the numberOfDataPoints equal to the highest value, 0s will be ignored

#formula for age = MX + C, this is the average age for each age category. Each data set has a different set of values. 
Age.M<-c(5, 5, 5, 5)
Age.C<-c(17.5, 17.5, 15, 17.5) 

diseaseData <- cbind(FGFR2_C755G, FGFR2_C758G,  HRAS_G34A, FGFR3_A1948G) 
#combines sequence data all together in one data-frame. 
#Note that in this case "disease data" is disease alleles in father's sperm

alleleName <- c("C755G", "C758G", "G34A", "A1948G")
diseaseName <- c("Apert", "Apert", "Costello", "Thanatophoric_dysplasia")
geneName <- c("FGFR2", "FGFR2", "HRAS", "FGFR3")

firstDisease=4 #can specify which disease to start/end with
lastDisease=4

matrix.function<-function(r){
  
  maxN[diseaseIdentifier] <- round(divisionsPerCellPerYear*maxAge*n)               #max number of cell divisions
  divisionsPerYear <- maxN[diseaseIdentifier]/maxAge                             #how many cell divisions per year, total, in one niche
  
  p_i<-function(i){
    (((n-i)/n)*(p+((i*r)/(n+1))))} #Equation 2.6
  
  BasicMatrix <- matrix(0, nrow=states, ncol=states) #empty matrix
  
  for (a in 1:states)                     #this makes the matrix in equation 2.7
  {BasicMatrix[a,a] <- (1-(p_i(a-1)))
   if(a<states){
     BasicMatrix[a, (a+1)]<- p_i(a-1)} 
  }
  
  for (counterN in 1:numberOfDataPoints[diseaseIdentifier]){ #this loop works out every value for age category used 
        
    currentN <- floor((counterN*Age.M[diseaseIdentifier]+Age.C[diseaseIdentifier]) * divisionsPerYear) #number of cell divisions to this age category
    
    output.data[counterN,1]<-counterN*Age.M[diseaseIdentifier]+Age.C[diseaseIdentifier] #age for this age category
    
    ExponentialMatrix <- BasicMatrix %^% currentN #This is T^K
    
    toprow<-ExponentialMatrix[1,] #equation 2.9
    
    for (c in 1:states)
    {toprow[c] <-(toprow[c] * (c-1))}
    
    sum.toprow<-(sum(toprow)) #equation 2.10
    
    output.data[counterN,2]<-(sum.toprow) #number of mutants per niche
    output.data[counterN,3]<- (output.data[counterN,2])*(mutationsPerNumber/n) #this gives # mutants/million cells (or whatever value).
  }
  Distances <- (output.data[,4]-output.data[,3]) 
  
  DistancesSquared <- (Distances^2) 
  
  sumOfSquaredDistances <- sum(DistancesSquared)
  cat(alleleName[diseaseIdentifier], "r", r, "SOS", sumOfSquaredDistances, "\n")
  
  if(run.through==1){
    return (sumOfSquaredDistances)
  }else{
    return (output.data)
  }

}  

for (diseaseIdentifier in firstDisease:lastDisease){
  run.through=1
  p=p.list[diseaseIdentifier] #sets the probability of mutation
  
  output.data <- cbind(basic.output.data, diseaseData[,diseaseIdentifier])
  
  Optimum_r <- optim(0.005, matrix.function, method="Brent", lower=0, upper=1) #this is the main optimize function
    
  run.through=2 #sets the return from the function to be the output data frame
  r=Optimum_r$par #this sets r to the optimum value found by function Optimum_r
  
  output.data <- matrix.function(r) #runs the funtion again with the optimum function
  
  information.names <- c("", "Disease", "Gene", "Mutation", "p", "r (extrapolated)", "Niche size", "Total cell divisions", "")
  information.data <- c(diseaseIdentifier, diseaseName[diseaseIdentifier], geneName[diseaseIdentifier], alleleName[diseaseIdentifier], p, r, n, maxN[diseaseIdentifier], "")
  
  output.data2<-matrix.function(0) #runs the function again with a r value of 0
  
  output.data <- cbind(information.names, information.data, output.data) 
    
    if(diseaseIdentifier == firstDisease){
      final.output.data <- output.data
      final.output.data2 <- output.data2
    }else{
      final.output.data <- rbind(final.output.data, output.data) # combines the output into a single file
      final.output.data2 <- rbind(final.output.data2, output.data2)
    }
}

final.output.data <- cbind(final.output.data, final.output.data2[,3]) #adds the r=0 data with no positive selection.
colnames(final.output.data)[7]<-"O/E r0"

write.csv(final.output.data, file=paste("1Molecular_data_n_", n, "_OE_comparison_one_variable.csv", sep="")  )

