#Setting 2: PANSS as a Surrogate for CGI
##Linear model with random intercept

library(Surrogate)
library(matrixcalc)

###cgi as true endpoint in setting 4
Sigma<- matrix(c(0.03685 , NA,0.02377 ,NA,
                 NA,0.04968 , NA, 0.03410,
                 0.02377, NA,0.04581 , NA,
                 NA, 0.03410, NA, 0.06830), nrow=4, ncol=4)

##estimable corr
corr<- cov2cor(Sigma)
D_matrix<- matrix(c( 0.03301, NA, 0.02131,NA,
                     NA, 0.02588, NA, 0.02289,
                     0.02131, NA,0.02801 , NA,
                     NA, 0.02289, NA, 0.02428), nrow=4, ncol=4)

d_corr<- cov2cor(D_matrix)

ICA_1<- ICA.ContCont(T0S0=corr[1,3], T1S1=corr[2,4], T0T0=Sigma[1,1], T1T1=Sigma[2,2], S0S0=Sigma[3,3], S1S1=Sigma[4,4], T0T1=seq(-1, 1, by=.2), 
                     T0S1=seq(-1, 1, by=.2), T1S0=seq(-1, 1, by=.2), S0S1=seq(-1, 1, by=.2))

rho_delta_1<- ICA_1$ICA
pos_def_1<-ICA_1$Pos.Def

ICA_2<- ICA.ContCont(T0S0=d_corr[1,3], T1S1=d_corr[2,4], T0T0=D_matrix[1,1], T1T1=D_matrix[2,2], S0S0=D_matrix[3,3], S1S1=D_matrix[4,4], T0T1=seq(-1, 1, by=.1), 
                     T0S1=seq(-1, 1, by=.1), T1S0=seq(-1, 1, by=.1), S0S1=seq(-1, 1, by=.1))


pos_def_2<-ICA_2$Pos.Def
rho_delta_2<- ICA_2$ICA

p<-6
AR<- 0.6890
alpha<- (p-(p*AR)+2*AR)/(1+AR)
Q<- matrix(c(-1, 0,1,0,
             0, -1,0,1), nrow=2, ncol=4)
I6<- diag(6)

ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}
r<-as.matrix(ar1_cor(6, AR), nrow=6, ncol=6)
p1 <- matrix(c(1,1,1,1,1,1), nrow=1, ncol=6)



#RH_square<- data.frame()
for (i in 1:length(rho_delta_1)) {
  for (j in 1:length(rho_delta_2)) {
    
    ##covariances
    cov_sigma_T0T1<-  pos_def_1[i,1]* sqrt(Sigma[1,1]*Sigma[2,2])
    cov_sigma_T0S1<-  pos_def_1[i,3]* sqrt(Sigma[1,1]*Sigma[4,4])
    cov_sigma_T1S0<-  pos_def_1[i,4]* sqrt(Sigma[2,2]*Sigma[3,3])
    cov_sigma_S0S1<-  pos_def_1[i,6]* sqrt(Sigma[3,3]*Sigma[4,4])
    
    
    cov_d_T0T1<-  pos_def_2[j,1]* sqrt(D_matrix[1,1]*D_matrix[2,2])
    cov_d_T0S1<-  pos_def_2[j,3]* sqrt(D_matrix[1,1]*D_matrix[4,4])
    cov_d_T1S0<-  pos_def_2[j,4]* sqrt(D_matrix[2,2]*D_matrix[3,3])
    cov_d_S0S1<-  pos_def_2[j,6]* sqrt(D_matrix[3,3]*D_matrix[4,4])
    
    
    V_1<- sigma_u <- matrix (c((Sigma[1,1])+(Sigma[2,2])-2*cov_sigma_T0T1,
                               (Sigma[1,3])+(Sigma[2,4])- cov_sigma_T1S0-cov_sigma_T0S1,
                               (Sigma[1,3])+(Sigma[2,4])- cov_sigma_T1S0-cov_sigma_T0S1,
                               (Sigma[3,3])+(Sigma[4,4])-2*cov_sigma_S0S1), nrow=2, ncol=2)
    
    
    D_matrix[1,2]<- D_matrix[2,1]<- cov_d_T0T1
    D_matrix[1,4]<- D_matrix[4,1]<- cov_d_T0S1
    D_matrix[3,2]<- D_matrix[2,3]<- cov_d_T1S0
    D_matrix[3,4]<- D_matrix[4,3]<- cov_d_S0S1
    
    D_1<- Q%*% D_matrix %*% t(Q)
    D_1<- round(D_1, digits = 6)
    
    #D_1<- sigma_d <- matrix (c((D_matrix[1,1])+(D_matrix[2,2]) - 2*cov_d_T0T1,
    # (D_matrix[1,3])+(D_matrix[2,4])- cov_d_T1S0-cov_d_T0S1,
    # (D_matrix[1,3])+(D_matrix[2,4])- cov_d_T1S0-cov_d_T0S1,
    # (D_matrix[3,3])+(D_matrix[4,4])- 2*cov_d_S0S1), nrow=2, ncol=2)
    
    if (is.positive.definite(V_1)==TRUE & is.positive.definite(D_1)==TRUE ) {
      
      sigma_delta<- kronecker(D_1, t(p1)%*%(p1))+  kronecker(V_1, r)
      rh_square_check<- 1- (det(sigma_delta)/(det(sigma_delta[1:6, 1:6])*det(sigma_delta[7:12,7:12])))
      write(as.numeric(rh_square_check),"~/Downloads/RH_square_case_2_model_2.txt",  append=TRUE, sep = "\t")
      #RH_square<- c(RH_square, RH_square_comb)
      
    }
  }
}


RH_square_check<-read.delim(file= '~/Downloads/RH_square_case_2_model_2.txt', header = FALSE, sep = "\t", dec = ".")
summary(RH_square_check$V1)
min(RH_square_check$V1)
quantile(RH_square_check$V1)
mean(RH_square_check$V1)
median(RH_square_check$V1)
length(RH_square_check$V1)

hist( as.numeric(RH_square_check$V1), main="", xlab=bquote("R"[Lambda]^2),  breaks = 15, labels = FALSE,ylim = c(0,6000000))
grid()
box()



margin<- qt(0.95,df=length(as.numeric(RH_square_check$V1))-1)*sd(as.numeric(RH_square_check$V1))/sqrt(length(as.numeric(RH_square_check$V1)))
lowerinterval <- mean(as.numeric(RH_square_check$V1)) - margin
upperinterval <- mean(as.numeric(RH_square_check$V1)) + margin