

#################################################################################################################
####################################FUNCTIONS####################################################################
#################################################################################################################

CalcAngleBlock <- function(Data=NULL){
  AngleBlock <- matrix(nrow=dim(Data)[1], ncol=dim(Data)[1])
  
  colnames(AngleBlock) <- rownames(AngleBlock) <- Data$trt
  for (i in 1:dim(Data)[1]){
    for (j in (1:dim(Data)[1])[-i]){
      c2 <- as.numeric(Data$cost[j])
      c1 <- as.numeric(Data$cost[i])
      
      q2 <- as.numeric(Data$qaly[j])
      q1 <- as.numeric(Data$qaly[i])
      
      this.angle <- atan2(c2 - c1, q2 - q1)
      if (this.angle < 0){ this.angle <- 2*pi + this.angle }
      AngleBlock[i,j] <- this.angle
    }  
  }
  return(AngleBlock)
}

CalcDistBlock <- function(Data=NULL){
  DistBlock <- matrix(nrow=dim(Data)[1], ncol=dim(Data)[1])
  colnames(DistBlock) <- rownames(DistBlock) <- Data$trt
  for (i in 1:dim(Data)[1]){
    for (j in (1:dim(Data)[1])[-i]){
      
      c2 <- as.numeric(Data$cost[j])
      c1 <- as.numeric(Data$cost[i])
      
      q2 <- as.numeric(Data$qaly[j])
      q1 <- as.numeric(Data$qaly[i])
      
      this.dist <- (  (c2 - c1)^2   + (q2 - q1)^2 )^0.5
      DistBlock[i,j] <- this.dist
    }  
  }
  return(DistBlock)
}

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
  
  continue <- T
  this.position <- start.position
  FrontierAngles <- c()
  
  while(continue==T){
    if(all(AngleBlock[this.position,]> pi / 2, na.rm=T) ){
      continue <- F}
    else {
      next.position <- which.min(AngleBlock[this.position,])
      PathPosition <- c(PathPosition, next.position)
      PathName <- paste(PathName, 
                        as.character(Data$trt[next.position]), 
                        sep="_")
      FrontierAngles <- c(FrontierAngles, 
                          AngleBlock[this.position, next.position])
      this.position <- next.position
    }
  }
  names(PathPosition) <- NULL
  return(list(Path=PathName, PathIndex=PathPosition, Angles=FrontierAngles, AngleBlock=AngleBlock, DistBlock=DistBlock))
}

Angles2Icers <- function(angles=NULL){
  icers <- rep(NA, length(angles))
  icers[angles < pi / 2] <- tan(angles[angles < pi / 2])
  icers[angles >= pi / 2 & angles < pi] <- Inf
  icers[angles > 3* pi / 2] <- -Inf
  return(icers)
}

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

