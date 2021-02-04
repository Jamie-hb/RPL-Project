numberinprogramme <- function(t){
  Ids <- demographics$IDMAT[demographics$dayssincefirstconsultation >= t & !is.na(demographics$dayssincefirstconsultation)]
  return(length(Ids))
}

# event OR in programme >= t days
f0 <- function(t){
  # find those without an event but >= t days
  Ids <- demographics$IDMAT[demographics$dayssincefirstconsultation >= t & 
                              !(demographics$IDMAT %in% firstpregnancy$IDMAT) & 
                              !is.na(demographics$dayssincefirstconsultation)]
  
  # add on those with an event and return
  return(length(Ids) + length(firstpregnancy$IDMAT))
}

# viable pregnancy with no miscarriage
f1 <- function(t){
  res <- 0
  for(i in 1:length(firstpregnancy$IDMAT)){
    # extract the information for woman i
    subframe <- postpregnancies[postpregnancies$IDMAT==firstpregnancy$IDMAT[i],]
    outcomes <- subframe$outcome
    conceptiondate <- subframe$ConceptionDate
    dayssincefirstconsultation <- firstpregnancy$dayssincefirstconsultation[i]
    DateOfFirstConsultation <- firstpregnancy$DateOfFirstConsultation[i]
    
    # see if woman i achieves within t days a viable pregnancy before a miscarriage, res +=1 if so
    L <- length(outcomes)
    
    for(j in 1:L){
      # check for viable pregnancy
      if(!is.na(outcomes[j]) & outcomes[j]==1){
        
       # check to see if the viable pregnancy was achieved within t days
       daystillpregnant <- as.numeric(conceptiondate[j] - DateOfFirstConsultation)
       if(!is.na(daystillpregnant) & daystillpregnant <= t){
          res <- res + 1
          break
        }
      }
      
      # check for miscarriage, break if so  
      if(!is.na(outcomes[j]) & outcomes[j] %in% c(2,3,4,5,10)){
        break
      }
    }
    
  }
  
  return(res)
}


# <24 weeks with no miscarriage       
f2 <- function(t){
  res <- 0
  for(i in 1:length(firstpregnancy$IDMAT)){
    # extract the information for woman i 
    subframe <- postpregnancies[postpregnancies$IDMAT==firstpregnancy$IDMAT[i],]
    outcomes <- subframe$outcome
    conceptiondate <- subframe$ConceptionDate
    dayssincefirstconsultation <- firstpregnancy$dayssincefirstconsultation[i]
    DateOfFirstConsultation <- firstpregnancy$DateOfFirstConsultation[i]
    
    # see if woman i achieves within t days a <24 weeks pregnancy before a miscarriage, res +=1 if so
    L <- length(outcomes)
    for(j in 1:L){
      # check for < 24 weeks pregnancy
      if(!is.na(outcomes[j]) & outcomes[j]==10){
        # check if this pregnancy was achieved within t days
        daystillpregnant <- as.numeric(conceptiondate[j] - DateOfFirstConsultation)
        if(!is.na(daystillpregnant) & daystillpregnant <= t){
         res <- res + 1
         break
        }
      }
      
      # check for miscarriage or viable pregnancy, break if so  
      if(!is.na(outcomes[j]) & outcomes[j] %in% c(1,2,3,4,5)){
        break
      }
    }
    
  }
  
  return(res)
}


# viable pregnancy after miscarriage
f3 <- function(t){
  res <- 0
  for(i in 1:length(firstpregnancy$IDMAT)){
    # extract the information for woman i
    subframe <- postpregnancies[postpregnancies$IDMAT==firstpregnancy$IDMAT[i],]
    outcomes <- subframe$outcome
    conceptiondate <- subframe$ConceptionDate
    dayssincefirstconsultation <- firstpregnancy$dayssincefirstconsultation[i]
    DateOfFirstConsultation <- firstpregnancy$DateOfFirstConsultation[i]
    
    # create a miscarriage indicator 
    mis <- 0
    
    L <- length(outcomes)
    # see if woman i achieves within t days a viable pregnancy following a miscarriage, if so res+=1
    for(j in 1:L){
      # check for miscarriage, mis += 1 if so
      if(!is.na(outcomes[j]) & outcomes[j] %in% c(2,3,4,5)){
        mis <- 1
      }
        
      if(!is.na(outcomes[j]) & outcomes[j]==1){
        # check for previous miscarriage
        if(mis==1){
          # check to see if viable pregnancy was within t days
          daystillpregnant <- as.numeric(conceptiondate[j] - DateOfFirstConsultation)
          if(!is.na(daystillpregnant) & daystillpregnant <= t){
            res <- res + 1
            break
          }
        }
      }
        
    }
    
  }
  
  return(res)
}


# < 24 weeks after miscarriage
f4 <- function(t){
  res <- 0
  for(i in 1:length(firstpregnancy$IDMAT)){
    # extract the information for woman i
    subframe <- postpregnancies[postpregnancies$IDMAT==firstpregnancy$IDMAT[i],]
    outcomes <- subframe$outcome
    conceptiondate <- subframe$ConceptionDate
    dayssincefirstconsultation <- firstpregnancy$dayssincefirstconsultation[i]
    DateOfFirstConsultation <- firstpregnancy$DateOfFirstConsultation[i]
    
    # create a miscarriage indicator for woman i 
    mis <- 0

    L <- length(outcomes)
    # see if woman i achieves within t days a viable pregnancy following a miscarriage, if so res+=1
    # check for miscarriage, mis += 1 if so
    for(j in 1:L){
      if(!is.na(outcomes[j]) & outcomes[j] %in% c(2,3,4,5)){
        mis <- 1
      }
      
      # check for < 24 weeks pregnancy  
      if(!is.na(outcomes[j]) & outcomes[j]==10){
        # check for previous miscarriage
        if(mis==1){
          # check to see if pregnancy was within t days
          daystillpregnant <- as.numeric(conceptiondate[j] - DateOfFirstConsultation)
          if(!is.na(daystillpregnant) & daystillpregnant <= t){
            res <- res + 1
            break
          }
        }
      }
        
    }
    
  }
  
  return(res)
}


# miscarriage with no successful pregnancy
f5 <- function(t){
  res <- 0
  for(i in 1:length(firstpregnancy$IDMAT)){
    # extract information for woman i 
    subframe <- postpregnancies[postpregnancies$IDMAT==firstpregnancy$IDMAT[i],]
    outcomes <- subframe$outcome
    conceptiondate <- subframe$ConceptionDate
    dayssincefirstconsultation <- firstpregnancy$dayssincefirstconsultation[i]
    DateOfFirstConsultation <- firstpregnancy$DateOfFirstConsultation[i]
    
    # create miscarriage indicator for woman i
    mis <- 0
    
    # check for miscarriage
    if(1 %in% outcomes | 2 %in% outcomes | 4 %in% outcomes | 5 %in% outcomes){
      # check for lack of pregnancy
      if(!(1 %in% outcomes) & !(10 %in% outcomes)){
        # check for conception within t days
        if(!is.na(conceptiondate[1]) & !is.na(DateOfFirstConsultation) & as.numeric(conceptiondate[1] - DateOfFirstConsultation) <= t){
          res <- res + 1          
        }
      }
    }
        
  }
  
  return(res)
}


# ======= Make the plot =======

t <- seq(0,1200,100)
y0 <- sapply(t, f0)
y1 <- sapply(t, f1)
y2 <- sapply(t, f2)
y3 <- sapply(t, f3)
y4 <- sapply(t, f4)
y5 <- sapply(t, f5)

# plot
plot(t, y0, type="l", lwd=2, xlab="days since first consultation", ylab="number of women",
     main="Outcomes After t Days")
lines(t, y1, col="darkgreen", lwd=2)
lines(t, y1 + y2, col="blue", lwd=2)
lines(t, y1 + y2 + y3, col="green", lwd=2)
lines(t, y1 + y2 + y3 + y4, col="red", lwd=2)
lines(t, y1 + y2 + y3 + y4 + y5, col="magenta", lwd=2)

# shade in 
polygon(c(t,rev(t)), c(rep(0,length(y1)),rev(y1)), col="darkgreen", density=20, border=NA)
polygon(c(t,rev(t)), c(y1,rev(y1+y2)), col="blue", density=20, border=NA)
polygon(c(t,rev(t)), c(y1+y2,rev(y1+y2+y3)), col="green", density=20, border=NA)
polygon(c(t,rev(t)), c(y1+y2+y3,rev(y1+y2+y3+y4)), col="red", density=20, border=NA)
polygon(c(t,rev(t)), c(y1+y2+y3+y4,rev(y1+y2+y3+y4+y5)), col="magenta", density=20, border=NA)

# add a legend
legend("topright",legend=c("number in programme >= t days", "viable pregnancy with no miscarriage", "< 24 weeks pregnancy with no miscarriage",
                "viable pregnancy after a miscarriage", "< 24 weeks pregnancy after a miscarriage",
                "miscarriage but no viable or <24 weeks pregnancy"), 
       col=c("black","darkgreen","blue","green","red","purple"), lty=1, lwd=2)

