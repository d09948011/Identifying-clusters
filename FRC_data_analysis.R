## load functions ##
rm(list=ls())
source("Identifying_Clusters_AUX.R")
source("Identifying_Clusters_MAIN.R")

# the matrix used to store needed informations
playoff_D_BIC = playoff_D = playoff_P_BIC = playoff_P = playoff_Y_BIC = playoff_Y = matrix(0,12,8)
colnames(playoff_D_BIC) = colnames(playoff_D) = colnames(playoff_P_BIC) = colnames(playoff_P) = colnames(playoff_Y_BIC) = colnames(playoff_Y) = c("Division","cluster", "PCP", "MSPE_P", "MSPE_Y", "quaterfinal", "semifinal", "final")

d = 1 # the threshold value for reestimate robots strengths 
#set d equal to the number of the robots for LCT method

for(div in 1:12){ # div in 1:12 represent "Carver", "Galileo", "Hopper", "Newton", "Roebling", "Turing", "Archimedes", "Carson", "Curie", "Daly", "Darwin", "Tesla" division 
 
  load(paste("2019_FRC_Data_" ,div,".RData", sep = ""))
  # Import RData file of a specific division data of the qualification and playoff stage       
  # The input files are named 2019_FRC_Data_div.RData (or 2018_FRC_Data_div.RData) with the following format:
  # the 1st column: match numbers
  # the 2nd to 4th columns: robot IDs in the red alliance
  # the 5th to 7th columns: robot IDs in the blue alliance
  # the 8th column: a score of the red alliace
  # the 9th column: a score of the blue alliance
  # the 10th to 15th column: NA
  # the remaining columns: covariate matrix X with X_ij being equal to 1 if the jth robot at the ith match 
  # assigned to the red alliance; -1 if the jth robot at the ith match assigned to the blue alliance; and 
  # 0 otherwise. 
  
  Data = data.frame(Qualification)  
  Data_po = data.frame(Playoff) # If the match is unnecessary, NA's are stored in that row

  # number of matches (M), number of robotics teams (K), and IDs in the qualification stage 
  M = length(Data[,1]) 
  K = length(Data[1,]) - 15
  
  # the corresponding response y0 and covariates x0, which are derived from qualification data
  yr = as.matrix(Data[1:M,8])
  yb = as.matrix(Data[1:M,9])
  y0 = as.matrix(yr-yb)
  ny = length(y0)
  x0 = as.matrix(Data[1:M,1:K+15])
  
  # remove matches with extremely high difference in scores (100 for 2019; 800 for 2018) to ensure better quality of data
  if (max(abs(y0)) > 100){
    outlier = which(abs(y0) >100 )
    M = M - length(outlier)
    x0 = x0[-outlier,]
    xn = x0[ , 1: ( K -1 )] - x0[ ,K] 
    y0 = y0[abs(y0)<100]
  }
  
  beta_Est = Classification.fn(x0,y0,d) # perform the classififcation 
  
  #######################################################
  ### prediction of the playoff based on six criteria ###
  #######################################################
  
  # data of the playoff
  M_po = length(t(Data_po[,1])) 
  K_po = length(Data_po[1,]) - 15
  
  yr_po = as.matrix(Data_po[1:M_po,8])
  yb_po = as.matrix(Data_po[1:M_po,9])
  y0_po = as.matrix(yr_po-yb_po)
  x0_po = as.matrix(Data_po[1:M_po,1:K_po+15])
  po = colSums((abs(x0_po)))
  
  # MSPE_D BIC
  b_inf_penality = log(1 - beta_Est[,2])+beta_Est[,1]*log(M)/M
  clu_opt = max(which(b_inf_penality == min(b_inf_penality), arr.ind = TRUE))
  cluster_info = beta_Est[clu_opt,5:(K+4)]
  b_opt = as.matrix(beta_esti.fn((outer(cluster_info,sort(unique(cluster_info)),FUN="==")*1)))

  playoff_D_BIC[div,]  = pred.fn(b_opt)
  
  # MSPE_D
  clu_opt = max(which(beta_Est[,2] == max(beta_Est[,2]), arr.ind = TRUE))
  cluster_info = beta_Est[clu_opt,5:(K+4)]
  b_opt = as.matrix(beta_esti.fn((outer(cluster_info,sort(unique(cluster_info)),FUN="==")*1)))
  
  playoff_D[div,]  = pred.fn(b_opt)
  
  # MSPE_P BIC
  b_inf_penality = log(beta_Est[,3])+beta_Est[,1]*log(M)/M
  clu_opt = max(which(b_inf_penality == min(b_inf_penality), arr.ind = TRUE))
  cluster_info = beta_Est[clu_opt,5:(K+4)]
  b_opt = as.matrix(beta_esti.fn((outer(cluster_info,sort(unique(cluster_info)),FUN="==")*1)))
  
  playoff_P_BIC[div,]  = pred.fn(b_opt)
  
  # MSPE_P
  clu_opt = max(which(beta_Est[,3] == min(beta_Est[,3]), arr.ind = TRUE))
  cluster_info = beta_Est[clu_opt,5:(K+4)]
  b_opt = as.matrix(beta_esti.fn((outer(cluster_info,sort(unique(cluster_info)),FUN="==")*1)))
  
  playoff_P[div,]  = pred.fn(b_opt)
  
  # MSPE_Y BIC
  b_inf_penality = log(beta_Est[,4])+beta_Est[,1]*log(M)/M
  clu_opt = max(which(b_inf_penality == min(b_inf_penality), arr.ind = TRUE))
  cluster_info = beta_Est[clu_opt,5:(K+4)]
  b_opt = as.matrix(beta_esti.fn((outer(cluster_info,sort(unique(cluster_info)),FUN="==")*1)))

  playoff_Y_BIC[div,]  = pred.fn(b_opt)
  
  # MSPE_Y
  clu_opt = max(which(beta_Est[,4] == min(beta_Est[,4]), arr.ind = TRUE))
  cluster_info = beta_Est[clu_opt,5:(K+4)]
  b_opt = as.matrix(beta_esti.fn((outer(cluster_info,sort(unique(cluster_info)),FUN="==")*1)))
  
  playoff_Y[div,]  = pred.fn(b_opt)

}

## save the result ##
filename="2019_data_analysis"
save.image(paste(filename,"_", d, ".RData", sep="")) 
  
  