##Selective Mutation Accumulation: A Compuational Model of the Paternal Age Effect
##Eoin C Whelan*, Alexander C Nwala, Christopher Osgood, Stephan Olariu
##Old Dominion University Norfolk VA USA 
##*Contact: ewhel001@odu.edu

##Supplementary codes
##R script for generating O/E curves and matching them to existing data.

library(expm) #required for exponentially multiplying matricies

options(scipen=999)#removes scientific notation as an output

n<-324                  #number of cells per niche
states<-n+1             #number of rows/columns in the matrix    
maxAge <- 80            #maximum age (in years)
divisionRate <- 16     #how many days between cell divisions. I.e. each cell divides once per 16 days
divisionsPerCellPerYear<-365.25/divisionRate

maxN <- round(divisionsPerCellPerYear*maxAge*n)    #max number of extractions/cell divisions
divisionsPerYear <- maxN/maxAge                  #how many divisions PER NICHE per year

baseline_p<-(4)*(10^-11)  #muation rate per nucleotide per cell division

numberOfDataPoints = 7    #This is however many data points are in the final output

basic.output.data <- data.frame(Years=1:numberOfDataPoints, Mutants_Prop=1:numberOfDataPoints, Observed=1:numberOfDataPoints, Expected=1:numberOfDataPoints, O_E=1:numberOfDataPoints)

Apert <- c(0.26, 1.55, 0.94, 2.33, 1.42, 4.3, 5.2)
FOP <- c(0.43, 0.53, 1.23, 1.18, 2.08, 3.01, 3.21)
Crouzon <-c(0.59,	0.71,	1.77,	0.9,	1.04,	2.86,	5.56)
Pfeiffer <- c(0.33,	0.31,	1.46,	1.49,	4.6,	6.25,	5.88)
Costello <- c(0, 0.206, 1.878, 1.655, 2.443, 9.38, 7.354)
Marfan <- c(0, 0.43, 1.90, 1.31, 2.07, 0, 11.11)
ACH <-c(.4, .4, 1.1, .9, 2.8, 7.4, 4.8)
TD <- c(.4, .7, 1, 2.9, 2, 3.4, 5.3)
#These are literature O/E values corresponding to ages 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5

diseaseData <- cbind(Apert, FOP, Crouzon, Pfeiffer, Costello, Marfan, ACH, TD)

birthNumbers <-c(1098528, 1091918, 719933, 408643, 184490, 72113, 30648) #1966 birth numbers.

diseaseName <- c("Apert", "FOP", "Crouzon", "Pfeiffer", "Costello", "Marfan", "Achondroplasia", "TD")

incidence <- c(100000, 2000000, 125000, 100000, 300000, 70000, 26000, 40000)

alleles <- c(2, 1, 11, 12, 14, 50, 1, 1)
p.list <- 1-((1-baseline_p)^alleles) #this adjusts the p-value of a mutation occurring by the number of alleles (see equation 3.1)

p.list[1] <-1-((1-baseline_p*5)*(1-baseline_p)) #this is to adjust for two mutable sites for Apert's syndrome, where one is a transversion at a CpG dinucleotide
p.list[2] <-baseline_p*15                       #transition at CpG dinucleotide
p.list[8] <-baseline_p*15                       #transition at CpG dinucleotide

firstDisease=1  #This specifies the disease to start or end with, in order above. 
lastDisease=8   #Setting firstDisease =1 and lastDisease=8 runs through all of the 8 diseases

matrix.function<-function(r){
  BasicMatrix <- matrix(0, nrow=states, ncol=states) #empty matrix

  p_i<-function(i){
    (((n-i)/n)*(p+((i*r)/(n+1))))   #Equation 2.6
    }
    
  for (a in 1:states)                       #fills the empty matrix
  {BasicMatrix[a,a] <- (1-(p_i(a-1)))
   if(a<states){
     BasicMatrix[a, (a+1)]<- p_i(a-1)}      #this makes the matrix in equation 2.7
  }
  
  for (counterN in 1:numberOfDataPoints){ #this loop works out every value for each age category used 
        
    currentN <- floor((counterN*5+17.5) * divisionsPerYear) #This calculates a single N value (number of cell divisions total in the niche at this specific age)
    
    output.data[counterN,1]<-counterN*5+17.5 #this calculates equation 2.8
    
    ExponentialMatrix <- BasicMatrix %^% currentN   #This is T^K
    
    toprow<-ExponentialMatrix[1,] #equation 2.9
    
    for (c in 1:states)
    {toprow[c] <-(toprow[c] * (c-1))}
    
    sum.toprow<-(sum(toprow)) #equation 2.10
    
    output.data[counterN,2]<-(sum.toprow/n) #SO IT GIVES A PROPORTION OF MUTANT:TOTAL
    output.data[counterN,3]<- (output.data[counterN,2])*(output.data[counterN,6]) #This is proportion x number of births
  }
  sumBirths <- sum(output.data[,6])
  sumMutants <- sum(output.data[,3])
    
  for (i in 1:numberOfDataPoints){
    output.data[i, 4] <- (output.data[i, 6]/sumBirths) * sumMutants
    output.data[i, 5] <- output.data[i, 3]/output.data[i, 4]  #this calculates the S/P value equivalent to O/E
  }
  Distances <- (output.data[,5]-output.data[,7]) #this is the raw distance between calculated and actual data
  
  DistancesSquared <- (Distances^2) 
  
  sumOfSquaredDistances <- sum(DistancesSquared)
  if(run.through==0){
    return (sumOfSquaredDistances)  #initially returns sumOfSquaredDistances to the optim function
  }else{
    return (output.data)            #returns output.data when calculating final output values
  }
}  

for (diseaseIdentifier in firstDisease:lastDisease){
  p=p.list[diseaseIdentifier]
  run.through=0
  
  output.data <- cbind(basic.output.data, birthNumbers, diseaseData[,diseaseIdentifier])
    
  Optimum_r <- optim(0.005, matrix.function, method="Brent", lower=0, upper=1)  #MAIN FUNCTION
  
  r=Optimum_r$par
  
  run.through=1
  output.data<-matrix.function(r)   #Generates data with optimum r value
  output.data2<-matrix.function(0)  #This runs through with r=0
  
  sumBirths <- sum(output.data[,6])
  sumMutants <- sum(output.data[,3])       
  PredictedIncidenceRate <- (sumBirths/sumMutants)
  
  cat("Final values: \nDisease: ", diseaseName[diseaseIdentifier], "\np-value: ", p, "\nr value: ", r, 
      "\nPredicted incidence rate: 1 in", PredictedIncidenceRate, 
      "\nActual incidence rate: 1 in", (incidence[diseaseIdentifier]), "\n")
  
  information.names <- c("Disease", "Actual incidence", "Predicted incidence", "p", "r (extrapolated)", "Niche size", "Total cell divisions")
  information.data <- c(diseaseName[diseaseIdentifier], incidence[diseaseIdentifier], PredictedIncidenceRate, p, r, n, maxN)
  
  output.data <- cbind(information.names, information.data, output.data)
    
  if(diseaseIdentifier == firstDisease){
    final.output.data <- output.data
    final.output.data2 <- output.data2
  }else{
    final.output.data <- rbind(final.output.data, output.data)
    final.output.data2 <- rbind(final.output.data2, output.data2)
  }
}

final.output.data <- cbind(final.output.data, final.output.data2[,5])
colnames(final.output.data)[7]<-"O/E calc"
colnames(final.output.data)[9]<-"O/E disease"
colnames(final.output.data)[10]<-"O/E r=0"

write.csv(final.output.data, file=paste("Census_data_n_", n, "_OE_comparison.csv", sep=""))
