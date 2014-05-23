# 11/7/2012
# Path Robustness Simulations

# 1) Construct Fictitious Data

rm(list=ls())
setwd("X:/Path Robustness/")
source("PathRobustnessFunctions.R")

###########################################################################################################
###############################DATA #######################################################################
###########################################################################################################

# Using 10,000 draws from the the GERD model (Briggs et al)

Data <- read.csv("GERD_MC.csv")
# Tasks

# 1) Produce data frame for deterministic results only [done]
# 2) Plot deterministic results with letters as points [done]
# 3) Identify frontier for deterministic results using 
# algorithm [ done]

# 3) Format stochastic results as list of data.frames
# 4) plot stochastic results (perhaps for only 1000 items)
# 5) Calculate proportion of 1st part of MC results which match 
# first part of deterministic results
# 6) Same for other parts

# n.b. ROW 0 is deterministic; other rows are stochastic

Results.det <- data.frame(
  trt=LETTERS[1:6],
  qaly=NA,
  cost=NA
  )

Results.det$qaly <- c(
  Data[1,
       c("A.util", 
         "B.util", 
         "C.util",
         "D.util",
         "E.util",
         "F.util")
       ]
  )

Results.det$cost <- c(
  Data[1,
       c("A.cost",
         "B.cost",
         "C.cost", 
         "D.cost",
         "E.cost",
         "F.cost")
       ]
  )

# Plot this:
tiff("Points.tiff", 1000, 800)

plot(NA, type="n", xlim=c(35,50), ylim=c(600,1200), xlab="Health Effect", ylab="Cost")

for (i in 1:6){
  points(Results.det[i,"qaly"], Results.det[i,"cost"], pch=as.character(Results.det[i,"trt"]))
  }

dev.off()

########################################################################################
FrontierMeans <- CalcFrontier(Results.det)

# find ICERs of these elements
Angles2Icers(FrontierMeans$Angles)

# Test robustness

OptimalPath <- FrontierMeans$PathIndex
OptimalPath.label <- FrontierMeans$Path



# Now to do this for each of the other PSA runs

PSA.paths <- vector("list", 10000)
PSA.paths.labels <- vector("character", 10000)
det.opt.angles <- matrix(nrow=10000, ncol=3)
det.opt.icers <- matrix(nrow=10000, ncol=3)
  

for (i in 2:10001){
  tmp.data <-  data.frame(
      trt=LETTERS[1:6],
      qaly=NA,
      cost=NA
    )
  
  tmp.data$qaly <- c(
    Data[i,
         c("A.util", 
           "B.util", 
           "C.util",
           "D.util",
           "E.util",
           "F.util")
         ]
  )
  
  tmp.data$cost <- c(
    Data[i,
         c("A.cost",
           "B.cost",
           "C.cost", 
           "D.cost",
           "E.cost",
           "F.cost")
         ]
  )
  
  tmp.out <- CalcFrontier(tmp.data)
  
  det.opt.angles[i-1,] <- c(
    tmp.out$AngleBlock[3,1],
    tmp.out$AngleBlock[1,5],
    tmp.out$AngleBlock[5,2]
  )
  
  det.opt.icers[i-1,] <- Angles2Icers(det.opt.angles[i-1,])
  
  PSA.paths[[i-1]] <- tmp.out$PathIndex
  PSA.paths.labels[i-1] <- tmp.out$Path
}



table(PSA.paths.labels)[rev(order(table(PSA.paths.labels)))]


apply(det.opt.icers, 2, function(x) quantile(x, c(0,0.025, 0.10, 0.50, 0.90, 0.975,1)))

# How does median ICER compare with mean ICER?


tiff("Contours.tiff", 1000, 800)

# Want to know: if the det solution were followed, 
# what is the range of ICERs from option T to option T+1?
require(MASS)
contour(kde2d(Data$A.util, Data$A.cost), col="red", levels=seq(0.01, 0.07, by=0.01), 
        xlab="Health Effect (Weeks free of GERD)",
        ylab="Cost", main="Contour plot representations of PSA for each GERD option", 
        xlim=c(38,48), ylim=c(600,1200), drawlabels=F)
contour(kde2d(Data$B.util, Data$B.cost), col="green", levels=seq(0.01, 0.07, by=0.01), drawlabels=F, add=T)        
contour(kde2d(Data$C.util, Data$C.cost), col="blue", levels=seq(0.01, 0.07, by=0.01), drawlabels=F, add=T)        
contour(kde2d(Data$D.util, Data$D.cost), col="red", lty="dashed", levels=seq(0.01, 0.07, by=0.01), drawlabels=F, add=T)        
contour(kde2d(Data$E.util, Data$E.cost), col="green", lty="dashed", levels=seq(0.01, 0.07, by=0.01), drawlabels=F, add=T)        
contour(kde2d(Data$F.util, Data$F.cost), col="blue", levels=seq(0.01, 0.07, by=0.01), lty="dashed", drawlabels=F, add=T)        
legend("topleft", lwd=c(1,1,1,1,1,1), col=c("red", "green", "blue", "red", "green", "blue"), lty=c("solid", "solid", "solid", "dashed", "dashed", "dashed"), legend=c("A", "B", "C","D", "E", "F"))

# by using the same levels argument, we can compare between options in terms 
# of uncertainty. The more concentric circles, the less the uncertainty.

# Show Frontier

lines(as.numeric(Results.det[c(3,1,5,2),"qaly"]), as.numeric(Results.det[c(3,1,5,2), "cost"]), lwd=2, lty="dashed", col="grey")

for (i in 1:6){
  points(Results.det[i,"qaly"], Results.det[i,"cost"], pch=as.character(Results.det[i,"trt"]))
}

dev.off()


# Looking at net benefit at £20k

CalcNetBenefit <- function(Data, lambda=20){
  # Assuming cost unit is £1000 not £1
  nmb <- lambda * Data$qaly - Data$cost
  nhb <- (lambda * Data$qaly - Data$cost) / lambda
  return(list(nmb=nmb, nhb=nhb))
}

attach(Data)
dat.A <- data.frame(qaly=A.util, cost=A.cost)
dat.B <- data.frame(qaly=B.util, cost=B.cost)
dat.C <- data.frame(qaly=C.util, cost=C.cost)
dat.D <- data.frame(qaly=D.util, cost=D.cost)
dat.E <- data.frame(qaly=E.util, cost=E.cost)
dat.F <- data.frame(qaly=F.util, cost=F.cost)
detach(Data)


nb20.A <- CalcNetBenefit(dat.A, 20000)
nb20.B <- CalcNetBenefit(dat.B, 20000)
nb20.C <- CalcNetBenefit(dat.C, 20000)
nb20.D <- CalcNetBenefit(dat.D, 20000)
nb20.E <- CalcNetBenefit(dat.E, 20000)
nb20.F <- CalcNetBenefit(dat.F, 20000)


nb30.A <- CalcNetBenefit(dat.A, 30000)
nb30.B <- CalcNetBenefit(dat.B, 30000)
nb30.C <- CalcNetBenefit(dat.C, 30000)
nb30.D <- CalcNetBenefit(dat.D, 30000)
nb30.E <- CalcNetBenefit(dat.E, 30000)
nb30.F <- CalcNetBenefit(dat.F, 30000)


plot(density(nb20.B$nmb), xlim=c(700000, 1000000), main="Density plot of NMB @ £20k")
lines(density(nb20.A$nmb), lwd=2)
lines(density(nb20.C$nmb), lty=2)
lines(density(nb20.D$nmb), lwd=2, lty=2)
lines(density(nb20.E$nmb), lty=3)
lines(density(nb20.F$nmb), lwd=2, lty=3)
legend("topleft", lwd=c(2,1,1,2,1,2), lty=c(1,1,2,2,3,3), legend=c("A", "B", "C", "D", "E", "F"))

# Show uncertainty in ICER along path C A E B

plot(density(nb30.B$nmb), xlim=c(1100000, 1500000), main="Density plot of NMB @ £30k")
lines(density(nb30.A$nmb), lwd=2)
lines(density(nb30.C$nmb), lty=2)
lines(density(nb30.D$nmb), lwd=2, lty=2)
lines(density(nb30.E$nmb), lty=3)
lines(density(nb30.F$nmb), lwd=2, lty=3)
legend("topleft", lwd=c(2,1,1,2,1,2), lty=c(1,1,2,2,3,3), legend=c("A", "B", "C", "D", "E", "F"))


## Now want to know:
# 1) what proportion of PSA runs start with the optimal starting point?
# 2) what proportion of those that start with the correct start then 
# go onto the correct second position
# 3) etc

tmp.1 <- sapply(PSA.paths, function(x) x[1])

correct.1 <- length(which(tmp.1==OptimalPath[1]))/length(tmp.1)

tmp.2 <- sapply(PSA.paths, function(x) x[2])

correct.2 <- length(which(tmp.2==OptimalPath[2]))/length(tmp.2)

correct.12 <- length(which(tmp.1==OptimalPath[1] & tmp.2==OptimalPath[2]))/length(tmp.1)


############################## PLOTTING DATA ###########################################################
########################################################################################################

min(sapply(Sims, function(X) min(X[,1])))
max(sapply(Sims, function(X) max(X[,1])))

min(sapply(Sims, function(X) min(X[,2])))
max(sapply(Sims, function(X) max(X[,2])))


plot(NA, type="n", xlim=c(5,10), ylim=c(8, 15), 
     xlab="QALYs", ylab="Cost (£1000)", 
     main="PSA from nine hypothetical treatments")

points(Sims$A[,2] ~ Sims$A[,1], pch=".", col=1)
points(Sims$B[,2] ~ Sims$B[,1], pch=".", col=2)
points(Sims$C[,2] ~ Sims$C[,1], pch=".", col=3)
points(Sims$D[,2] ~ Sims$D[,1], pch=".", col=4)
points(Sims$E[,2] ~ Sims$E[,1], pch=".", col=5)
points(Sims$F[,2] ~ Sims$F[,1], pch=".", col=6)
points(Sims$G[,2] ~ Sims$G[,1], pch=".", col=7)
points(Sims$H[,2] ~ Sims$H[,1], pch=".", col=8)
points(Sims$I[,2] ~ Sims$I[,1], pch=".", col=9)

legend("topleft", pch=16, col=1:9, legend=LETTERS[1:9])


# now want means

Trt <- LETTERS[1:9]
TrtMeans <- data.frame(trt=Trt, cost=NA, qaly=NA)
TrtMeans[1,2] <- mean(Sims$A[,2])
TrtMeans[1,3] <- mean(Sims$A[,1])

TrtMeans[2,2] <- mean(Sims$B[,2])
TrtMeans[2,3] <- mean(Sims$B[,1])

TrtMeans[3,2] <- mean(Sims$C[,2])
TrtMeans[3,3] <- mean(Sims$C[,1])

TrtMeans[4,2] <- mean(Sims$D[,2])
TrtMeans[4,3] <- mean(Sims$D[,1])

TrtMeans[5,2] <- mean(Sims$E[,2])
TrtMeans[5,3] <- mean(Sims$E[,1])

TrtMeans[6,2] <- mean(Sims$F[,2])
TrtMeans[6,3] <- mean(Sims$F[,1])

TrtMeans[7,2] <- mean(Sims$G[,2])
TrtMeans[7,3] <- mean(Sims$G[,1])

TrtMeans[8,2] <- mean(Sims$H[,2])
TrtMeans[8,3] <- mean(Sims$H[,1])

TrtMeans[9,2] <- mean(Sims$I[,2])
TrtMeans[9,3] <- mean(Sims$I[,1])

for (i in 1:9){
  points(cost ~ qaly, data=TrtMeans[i,], pch=LETTERS[i], cex=2)
}

###########################################################################################
###########################################################################################

# Approach 1: Cluster analysis

# What do the data look like without identifying information?

AllPsa <- rbind(Sims$A, Sims$B, Sims$C, Sims$D, Sims$E, Sims$F, Sims$G, Sims$H, Sims$I)

fitk <- kmeans(AllPsa, 9)


# model based

require(mclust)
fit <- Mclust(AllPsa) # n.b. takes AGES

#plot(fit, AllPsa, xlim=c(5,10), ylim=c(8, 15))
plot(fit, AllPsa)




d2 <- data.frame(cluster=fitk$cluster, AllPsa)
d3 <- data.frame(cluster=fit$classification, AllPsa)
tmp <- rep(LETTERS[1:9], each=1000)
d2 <- data.frame(d2, trt=tmp)
d3 <- data.frame(d3, trt=tmp)
table(d2$cluster, d2$trt)
table(d3$cluster, d3$trt)

############################################################################################
############################################################################################

# What is the Frontier of the means?

FrontierMeans <- CalcFrontier(TrtMeans)

# find ICERs of these elements
Angles2Icers(FrontierMeans$Angles)

# Test robustness

OptimalPath <- FrontierMeans$PathIndex

RobustnessMatrix <- matrix(NA, nrow=1000, ncol=length(OptimalPath)-1)

for (i in 1:1000){
  Trt <- LETTERS[1:9]
  temp <- data.frame(trt=Trt, cost=NA, qaly=NA)
  for (j in 1:length(Trt)){
    temp[j,2] <- Sims[[j]][i,1]
    temp[j,3] <- Sims[[j]][i,2]
  }
  
  ptemp <- PathIcers(temp, FrontierMeans$PathIndex)
  RobustnessMatrix[i,] <- ptemp
}

# To Summarise:

SummaryRobustnessMatrix <- matrix(NA, nrow=5, ncol=length(OptimalPath)-1)

for (i in 1:dim(RobustnessMatrix)[2]){
  SummaryRobustnessMatrix[,i] <- quantile(RobustnessMatrix[,i], c(0.025, 0.05, 0.5, 0.95, 0.975), na.rm=T)
}


# Net benefit framework

NB <- list(
  NHB=matrix(NA, nrow=1000, ncol=9),
  NMB=matrix(NA, nrow=1000, ncol=9)
  )


for (i in 1:9){
  tmp <- CalcNetBenefit(Sims[[i]])  
  NB$NHB[,i] <- tmp$nhb
  NB$NMB[,i] <- tmp$nmb
}
  

plot(density(NB$NHB[,1]), xlim=c(5, 9), ylab="")
for (i in 2:9){
  lines(density(NB$NHB[,i]), col=i)
}

plot(density(NB$NMB[,1]), xlim=c(-200, 0), ylab="")
for (i in 2:9){
  lines(density(NB$NMB[,i]), col=i)
}
####################
###################
NB2 <- list(
  NHB=matrix(NA, nrow=1000, ncol=9),
  NMB=matrix(NA, nrow=1000, ncol=9)
  )


for (i in 1:9){
  tmp <- CalcNetBenefit(Sims[[i]], 30)  
  NB2$NHB[,i] <- tmp$nhb
  NB2$NMB[,i] <- tmp$nmb
}


plot(density(NB2$NHB[,1]), xlim=c(5, 9), ylab="")
for (i in 2:9){
  lines(density(NB2$NHB[,i]), col=i)
}

plot(density(NB2$NMB[,1]), xlim=c(-200, 0), ylab="")
for (i in 2:9){
  lines(density(NB2$NMB[,i]), col=i)
}

#####################################
# Testing the functions #############
#####################################

# Some data
Data <- data.frame(
  trt=LETTERS[1:6], 
  qaly=c(
   3.566773819,
   1.261496795,
   0.675316675,
   3.635356702,
   6.187548614,
   5.119131228),
  cost=c(
   3.900525701,
   0.987502088,
   7.484546366,
   1.174522543,
   2.630448942,
   5.908440068)                     
  )
  
plot(NA, type="n", ylim=c(0,10), xlim=c(0,10), main="", xlab="QALY", ylab="Cost (£10,000)")
for ( i in 1:dim(Data)[1]){
  points(cost[i] ~ qaly[i], data=Data, pch=as.character(Data$trt[i]))
}


CalcFrontier(Data)



##############################################################################
##############################################################################

X <- CalcFrontier(TrtMeans)

# Should atan2(y, x)
# if the values are negative subtract them from 2*pi

PathIcers(TrtMeans, c(3,1,4))
Data

          
#######################################################################################################
################### NET BENEFIT FRAMEWORK #############################################################
#######################################################################################################



############# OLDER CODE


#require(MASS)

# SimParams <- list(
#   A=list(
#     mu=c(6, 10),
#     Sigma=matrix(c(0.010412, 0.002603, 0.002603, 0.260308), nrow=2, byrow=T)
#     ),
#   B=list(
#     mu=c(7.5, 9.9),
#     Sigma=matrix(c(0.023428, 0.000625, 0.000625, 0.01412), nrow=2, byrow=T)
#     ),
#   C=list(
#     mu=c(7.5, 9.9),
#     Sigma=matrix(c(0.023425, 0.000625, 0.000625, 0.01412), nrow=2, byrow=T)
#     ),
#   D=list(
#     mu=c(6.9, 13.25),
#     Sigma=matrix(c(0.023428, 0.04217, 0.04217, 0.146423), nrow=2, byrow=T)
#     ),
#   E=list(
#     mu=c(8, 11.000),
#     Sigma=matrix(c(0.065077, -0.16373, -0.16373, 0.752291), nrow=2, byrow=T)
#     ),
#   F=list(
#     mu=c(8.3, 12.000),
#     Sigma=matrix(c(0.002603, 0.00026, 0.00026, 0.065077), nrow=2, byrow=T)
#     ),
#   G=list(
#     mu=c(9, 13.750),
#     Sigma=matrix(c(0.041649, 0.00859, 0.00859, 0.146423), nrow=2, byrow=T)
#     ),
#   H=list(
#     mu=c(8.7, 13.800),
#     Sigma=matrix(c(0.01412, -0.01921, -0.01921, 0.21085), nrow=2, byrow=T)
#     ),
#   I=list(
#     mu=c(9.2, 14.300),
#     Sigma=matrix(c(0.010412, -0.02952, -0.02952, 0.127551), nrow=2, byrow=T)
#     )
#   )
# 
# Sims <- list(
#   A=mvrnorm(1000, mu=SimParams$A$mu, Sigma=SimParams$A$Sigma),
#   B=mvrnorm(1000, mu=SimParams$B$mu, Sigma=SimParams$B$Sigma),
#   C=mvrnorm(1000, mu=SimParams$C$mu, Sigma=SimParams$C$Sigma),
#   D=mvrnorm(1000, mu=SimParams$D$mu, Sigma=SimParams$D$Sigma),
#   E=mvrnorm(1000, mu=SimParams$E$mu, Sigma=SimParams$E$Sigma),
#   F=mvrnorm(1000, mu=SimParams$F$mu, Sigma=SimParams$F$Sigma),
#   G=mvrnorm(1000, mu=SimParams$G$mu, Sigma=SimParams$G$Sigma),
#   H=mvrnorm(1000, mu=SimParams$H$mu, Sigma=SimParams$H$Sigma),
#   I=mvrnorm(1000, mu=SimParams$I$mu, Sigma=SimParams$I$Sigma)
#   )
# 
# 
# 
# colnames(Sims$A) <- c("qaly", "cost")
# colnames(Sims$B) <- c("qaly", "cost")
# colnames(Sims$C) <- c("qaly", "cost")
# colnames(Sims$D) <- c("qaly", "cost")
# colnames(Sims$E) <- c("qaly", "cost")
# colnames(Sims$F) <- c("qaly", "cost")
# colnames(Sims$G) <- c("qaly", "cost")
# colnames(Sims$H) <- c("qaly", "cost")
# colnames(Sims$I) <- c("qaly", "cost")
# 
# Sims$A <- data.frame(Sims$A)
# Sims$B <- data.frame(Sims$B)
# Sims$C <- data.frame(Sims$C)
# Sims$D <- data.frame(Sims$D)
# Sims$E <- data.frame(Sims$E)
# Sims$F <- data.frame(Sims$F)
# Sims$G <- data.frame(Sims$G)
# Sims$H <- data.frame(Sims$H)
# Sims$I <- data.frame(Sims$I)
# ########################################################################################################
