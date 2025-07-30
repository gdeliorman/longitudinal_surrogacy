#Setting 2: PANSS as a Surrogate for CGI
###Galecki's model
library(Surrogate)

##SAS output:
var_joint<- matrix(c(0.06345, NA,0.04101, NA,
                     NA, 0.06771, NA,0.04746 ,
                     0.04101, NA,0.07667, NA,
                     NA, 0.04746, NA, 0.0964), nrow=4, ncol=4)

Sigma<- var_joint
corr<- cov2cor(var_joint)

ICA<- ICA.ContCont(T0S0=corr[1,3], T1S1=corr[2,4], T0T0=Sigma[1,1], T1T1=Sigma[2,2], S0S0=Sigma[3,3], S1S1=Sigma[4,4], T0T1=seq(-1, 1, by=.1), 
                   T0S1=seq(-1, 1, by=.1), T1S0=seq(-1, 1, by=.1), S0S1=seq(-1, 1, by=.1))

rho<- ICA$ICA
p<-6
RH_square<- 1-(1-rho^2)^p
summary(RH_square)
min(RH_square)
max(RH_square)
quantile(RH_square)
mean(RH_square)
median(RH_square)
hist(RH_square)
length(RH_square)


hist(RH_square, main="", xlab=bquote("R"[Lambda]^2),  breaks = 15, labels = FALSE, ylim = c(0,15000), xlim = c(0,1))
grid()
box()


margin<- qt(0.95,df=length(RH_square)-1)*sd(RH_square)/sqrt(length(RH_square))
lowerinterval <- mean(RH_square) - margin
upperinterval <- mean(RH_square) + margin
