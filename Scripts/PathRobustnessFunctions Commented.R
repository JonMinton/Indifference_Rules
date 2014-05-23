

#################################################################################################################
####################################FUNCTIONS####################################################################
#################################################################################################################

# The functions expect a data frame as an input with three columns:
# trt: treatment name (e.g. A, B. C)
# cost: cost
# qaly : qaly (or other health benefit)

# the dataframe has one row per option



CalcAngleBlock <- function(Data=NULL){
  
  # Product an empty square matrix whose dimensions are equal to the number of options being compared
  
  AngleBlock <- matrix(nrow=dim(Data)[1], ncol=dim(Data)[1])
  
  # name the row names and column names according to the names of the treatments e.g. 
  
  #      A      B       C
  # A  [1, 1]   [1, 2]  [1, 3]
  # B  [2, 1]   [2, 2]  [2, 3]
  # C  [3, 1]   [3, 2]  [3, 3]
  
  
  # For each row, and for each column, find dcost and dqaly 
  # Use this information (cartesian coordinates)
  # to calculate the angle from the horizontal going in an anticlockwise direction, 
  # arranged so all values 
  colnames(AngleBlock) <- rownames(AngleBlock) <- Data$trt
  for (i in 1:dim(Data)[1]){
    for (j in (1:dim(Data)[1])[-i]){ # note the [-i], so the algorithm won't calculate the angle between an option and itself
      c2 <- as.numeric(Data$cost[j]) # 'to' cost
      c1 <- as.numeric(Data$cost[i]) # 'from' cost
      
      q2 <- as.numeric(Data$qaly[j]) # 'to' qaly
      q1 <- as.numeric(Data$qaly[i]) # 'from' qaly
      
      this.angle <- atan2(c2 - c1, q2 - q1) # angle in radians; angles below the horizontal are negative
      if (this.angle < 0){ this.angle <- 2*pi + this.angle } # turning negative angles into angles between 180 and 360 degrees
      AngleBlock[i,j] <- this.angle
    }  
  }
  return(AngleBlock)
}

# Distance matrix: not used but calculated for completeness
CalcDistBlock <- function(Data=NULL){
  DistBlock <- matrix(nrow=dim(Data)[1], ncol=dim(Data)[1])
  colnames(DistBlock) <- rownames(DistBlock) <- Data$trt
  for (i in 1:dim(Data)[1]){
    for (j in (1:dim(Data)[1])[-i]){
      
      c2 <- as.numeric(Data$cost[j])
      c1 <- as.numeric(Data$cost[i])
      
      q2 <- as.numeric(Data$qaly[j])
      q1 <- as.numeric(Data$qaly[i])
      
      this.dist <- (  (c2 - c1)^2   + (q2 - q1)^2 )^0.5 # hypotenuse 
      DistBlock[i,j] <- this.dist
    }  
  }
  return(DistBlock)
}


# Calculating efficiency frontier  - algorithm to be copied in vba in excel
CalcFrontier <- function(Data=NULL){
  
  AngleBlock <- CalcAngleBlock(Data)
  DistBlock <- CalcDistBlock(Data)
  # to see as degrees do :
  #AngleBlock * 180 / pi
  
  # The algorithm:
  # 1) Find option with min Y
  costs <- as.numeric(Data$cost)
  start.position <- which.min(costs)
  PathName <- as.character(Data$trt[start.position])
  PathPosition <- start.position

  # set flag to continue
  continue <- T
  # start at position 1
  this.position <- start.position
  # create empty vector for storing frontier angles as they are identified
  FrontierAngles <- c()
  
  # while loop: continues until all options in frontier have been identified
  while(continue==T){
    # find the row with the angles from the option of interest to all other options
    # check if any of these angles are less than 90 degrees (pi/2 radians)
    # If not, then the full frontier has been identified, so set continue to false to exit while loop
    # If yes, the other options in the frontier remain to be identified, 
    # and the else part of the logic statement is run
    if(all(AngleBlock[this.position,]> pi / 2, na.rm=T) ){
      # end of the frontier
      continue <- F}
    else {
      # find the next option by looking for the smallest angle on the 
      # row corresponding to the current option of interest
      next.position <- which.min(AngleBlock[this.position,])
      # add new position to current vector of positions of options on the frontier
      PathPosition <- c(PathPosition, next.position)
      # add to a character vector giving the name of the options on the frontier, 
      # in sequential order, each separated by a _ sign.
      # e.g. C_A_E_F  for C then A then E then F
      # If the options are A B C D E F G
      # then the Path position vector would then be
      # 3 1 5 6
      
      PathName <- paste(PathName, 
                        as.character(Data$trt[next.position]), 
                        sep="_")
      # add the newly identified angle of this section of the frontier to the 
      # list of angles of options on the frontier
      # NOTE: the vector FrontierAngles should be 1 element shorter than the 
      # vector of options on the frontier (as it takes two options to make one angle)
      FrontierAngles <- c(FrontierAngles, 
                          AngleBlock[this.position, next.position])
      # rename next.position as this.position then repeat the while loop
      this.position <- next.position
    }
  }
  # remove extraneous attributes from the vector PathPosition
  names(PathPosition) <- NULL
  # return a list with (almost) everything you might be interested in
  return(list(Path=PathName, PathIndex=PathPosition, Angles=FrontierAngles, AngleBlock=AngleBlock, DistBlock=DistBlock))
}

# Short function for converting from angles to ICERs
# In order to avoid creating a character vector,
# dominating options are given the value negative infinity
# and dominated options are given the value positive infinity 
# (disinvestment icers are missing values)
Angles2Icers <- function(angles=NULL){
  icers <- rep(NA, length(angles))
  icers[angles < pi / 2] <- tan(angles[angles < pi / 2])
  icers[angles >= pi / 2 & angles < pi] <- Inf
  icers[angles > 3* pi / 2] <- -Inf
  return(icers)
}

# helper function for showing the ICERs of the options on the efficiency frontier
PathIcers <- function(Data=NULL, path=NULL){
  ######
  AngleBlock <- CalcAngleBlock(Data)
  Icer <- rep(NA, length(path)-1)
  for (i in 1:length(Icer)){
    this.angle <- AngleBlock[path[i], path[i+1]]
    Icer[i] <- Angles2Icers(this.angle)
  }
  return(Icer)
}

###########################################################################################################
######################### NET BENEFIT #####################################################################
###########################################################################################################

CalcNetBenefit <- function(Data, lambda=20){
  # Assuming cost unit is £1000 not £1
  nmb <- lambda * Data$qaly - Data$cost
  nhb <- (Data$cost - lambda * Data$qaly ) / lambda
  return(list(nmb=nmb, nhb=nhb))
}

