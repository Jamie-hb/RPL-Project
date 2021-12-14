setwd("C:/Users/jamie/OneDrive/Documents/MMath Project/Summer2021")

library("coda")

################
#              #
# Simple Model #
#              #   
################

# data <- read.csv("simpleoutput.csv",header=FALSE)
# 
# output <- matrix(0,10001,2)
# output[,1] <- t(data[1,])
# output[,2] <- t(data[2,])
# 
# Output <- mcmc(data=output)
# 
# autocorr(Output)
# effectiveSize(Output)
# gelman.diag(Output)

#############################
#                           #
# Time Covariates Simulated #
#                           #     
#############################

# data <- read.csv("TimeCovariates.csv",header=FALSE)
# 
# Output = mcmc.list()
# 
# for(k in 1:10){
#  OutputTemp <- matrix(0,10001,5)
#  
#  for(i in  1:5){
#    OutputTemp[,i] <- t(data[5*(k-1)+i,])
#  }
#  
#  Output[[k]] <- mcmc(data=OutputTemp)
# }
# 
# autocorr(Output)
# lapply(Output,effectiveSize)
# gelman.diag(Output)

########################
#                      #
# Time Covariates Data #
#                      #     
########################

data <- read.csv("TimeCovariatesData.csv",header=FALSE)

OutputTemp <- matrix(0,10001,5)

for(i in 1:5){
  OutputTemp[,i] <- t(data[i,])
}

Output <- mcmc(data=OutputTemp)

effectiveSize(Output)

###############################
#                             #
# Time Covariates Diagnostics #
#                             #     
###############################

# data <- read.csv("TimeCovariatesDataDiagnostic.csv",header=FALSE)
# 
# Output = mcmc.list()
# 
# for(k in 1:10){
#  OutputTemp <- matrix(0,10001,5)
# 
#  for(i in  1:5){
#    OutputTemp[,i] <- t(data[5*(k-1)+i,])
#  }
# 
#  Output[[k]] <- mcmc(data=OutputTemp[5000:10001,])
# }
# 
# autocorr(Output)
# lapply(Output,effectiveSize)
# gelman.diag(Output)

################################
#                              #
# Fertile Covariates Simulated #
#                              #     
################################

# data <- read.csv("FertileCovariatesTest.csv",header=FALSE)
# 
# Output <- mcmc.list()
# 
# for(k in 1:10){
#  OutputTemp <- matrix(0,10001,5)
# 
#  for(i in  1:5){
#    OutputTemp[,i] <- t(data[5*(k-1)+i,])
#  }
# 
#  Output[[k]] <- mcmc(data=OutputTemp)
# }
# 
# autocorr(Output)
# lapply(Output,effectiveSize)
# gelman.diag(Output)

##################################
#                                #
# Fertile Covariates Diagnostics #
#                                #     
##################################

# data <- read.csv("FertileCovariatesDataDiagnostic.csv",header=FALSE)
# 
# Output = mcmc.list()
# 
# for(k in 1:10){
#  OutputTemp <- matrix(0,10001,5)
# 
#  for(i in  1:5){
#    OutputTemp[,i] <- t(data[5*(k-1)+i,])
#  }
# 
#  Output[[k]] <- mcmc(data=OutputTemp[5000:10001,])
# }
# 
# autocorr(Output)
# lapply(Output,effectiveSize)
# gelman.diag(Output)

###########################
#                         #
# Fertile Covariates Data #
#                         #     
###########################

# data <- read.csv("EventClassMCMC.csv",header=FALSE)
# 
# OutputTemp <- matrix(0,10001,6)
# 
# for(i in 1:6){
#   OutputTemp[,i] <- t(data[i,])
# }
# 
# Output <- mcmc(data=OutputTemp)
# 
# autocorr(Output)
# effectiveSize(Output)

############################
#                          #
# All Covariates Simulated #
#                          #     
############################

# data <- read.csv("FullModelTest.csv",header=FALSE)
# 
# Output <- mcmc.list()
# 
# K <- c(1,2,4,7,8,9,10)
# 
# for(k in 1:10){
#  OutputTemp <- matrix(0,10001,8)
# 
#  for(i in 1:8){
#    OutputTemp[,i] <- t(data[8*(k-1)+i,])
#  }
# 
#  Output[[k]] <- mcmc(data=OutputTemp)
# }
# 
# autocorr(Output)
# lapply(Output,effectiveSize)
# gelman.diag(Output)

#######################
#                     #
# All Covariates Data #
#                     #     
#######################

# data <- read.csv("FullModelData.csv",header=FALSE)
# 
# OutputTemp <- matrix(0,10001,8)
# 
# for(i in 1:8){
#   OutputTemp[,i] <- t(data[i,])
# }
# 
# Output <- mcmc(data=OutputTemp)
# 
# effectiveSize(Output)

##############################
#                            #
# All Covariates Diagnostics #
#                            #     
##############################

# data <- read.csv("FullModelDataDiagnostic.csv",header=FALSE)
# 
# Output = mcmc.list()
# 
# for(k in 1:10){
#  OutputTemp <- matrix(0,10001,5)
# 
#  for(i in  1:8){
#    OutputTemp[,i] <- t(data[5*(k-1)+i,])
#  }
# 
#  Output[[k]] <- mcmc(data=OutputTemp)
# }
# 
# autocorr(Output)
# lapply(Output,effectiveSize)
# gelman.diag(Output)
