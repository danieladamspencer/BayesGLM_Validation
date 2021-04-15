setwd('/Users/Shared/Lab_Data/HCP_Motor_Task_Dan/visit1_data/137128/MNINonLinear/Results/tfMRI_MOTOR_LR')
motion <- read.table('Movement_Regressors.txt', header = FALSE)
colors <- rainbow(11)
plot(1:284, motion[,1], type='l', ylim=c(-0.6,1.15))
for(k in 2:12){
  lines(1:284, motion[,k], col=colors[k-1])
}

design <- readRDS('~/Downloads/sub10_visit1_sessionLR_design.rds')

library(clever)
myFD <- FD(X = motion[,1:6])
plot(1:284, myFD$measure, type='l')

drift1 <- (1:284)/284
drift <- cbind(drift1, drift1^2)
design2 <- nuisance_regression(design, cbind(motion, drift))

colors <- c('purple','blue','red', 'black')
plot(1:284, design[,1], col=colors[1], type='l', ylim=c(-0.5,1))
for(k in 2:4){
  lines(1:284, design[,k], col=colors[k])
}

lines(1:284, design2[,1], col=colors[1], type='l', lty=2)
for(k in 2:4){
  lines(1:284, design2[,k], col=colors[k], lty=2)
}

for(k in 1:4){
  print(summary(lm(design[,k] ~ as.matrix(cbind(motion, drift))))$r.squared)
}
