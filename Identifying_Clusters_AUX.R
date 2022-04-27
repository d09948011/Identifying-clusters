# ========================================================================== #
## ------------------------------------------------------------------------ ##
## Code for An Effective Method for Identifying Clusters of Robot Strengths ##
## ---------- Jen-Chieh Teng, Chin-Tsang Chinag, and Alvin Lim ------------ ##
## ------------------------------------------------------------------------ ##
# ========================================================================== #


## ------------------------------------------------------------------------ ##
## ----------------------- Various Helper Functions ----------------------- ##
## ------------------------------------------------------------------------ ##


# Estimate beta for different grouping
beta_esti.fn = function(grp_info){
  num_grp = dim(grp_info)[2] # number of clusters
  clu_size = colSums(grp_info) # each clusters' size
  
  w_clu = clu_size[1:(num_grp-1)]/clu_size[num_grp]  # weight of 1, ..., c-1 clusters
  dim(clu_size) = c(num_grp,1) 
  
  xc = x0%*%grp_info 
  xcn = xc[ ,1:(num_grp-1)]-xc[ ,num_grp]%*%t(w_clu) # constraint sum beta = 0
  
  bhatc = solve(t(xcn)%*%xcn)%*%t(xcn)%*%y0
  binf = c(bhatc,-sum(bhatc*w_clu))
  dim(binf) = c(num_grp,1)
  binf = grp_info%*%binf
  binf
  }


# estimated three types of MSPE by removing one out cross-validation
esti.fn = function(grp_info){
  num_grp = dim(grp_info)[2]
  clu_size = colSums(grp_info)
  
  w_clu = clu_size[1:num_grp-1]/clu_size[num_grp]
  
  y_cv = c(y0)
  xc_cv = x0%*%grp_info
  xcn_cv = xc_cv[ ,1:(num_grp-1)]-xc_cv[ ,num_grp]%*%t(w_clu) 
  H0 = solve(t(xcn_cv)%*%xcn_cv)%*%t(xcn_cv)
  H_cv = xcn_cv%*%H0
  
  pr_cv0 = (H_cv%*%y0)[1:M]
  resid_cv0 = y_cv-pr_cv0
  
  w_0 = resid_cv0/(1-diag(H_cv))
  pr_cv = pr_cv0-t(t(H_cv)*w_0)
  resid_cv = y_cv-pr_cv
  Fhat_cv0 = (t(t(resid_cv)+diag(pr_cv))<=0)*1
  diag(Fhat_cv0) = 0
  Fhat_cv = colSums(Fhat_cv0)/(M-1)
  
  c(sum((y_cv*(0.5-Fhat_cv)>0)+0.5*(y_cv*(0.5-Fhat_cv)==0))/M, sum( ( (y_cv> 0)*1 - (1-Fhat_cv))^2 )/M, sum((diag(resid_cv)^2))/M)
}


# prediction of the playoff
pred.fn = function(beta_opt){
  
  ################ Pred. of Playoff based on emp. distribution #################
  yhatpo = x0_po %*% b_opt
  Fhatpo = colSums ( outer( y0 - x0 %*% b_opt , -( yhatpo ) , FUN = "<=")*1  ) /  M
  na_seq = which(is.na(Fhatpo))
  Fhatpo[c(na_seq)] = Fhatpo[c(na_seq)-1]
  pred_po = ((0.5 - Fhatpo) > 0) * 1 + ((0.5 - Fhatpo) < 0) * (-1) +  ((0.5 - Fhatpo) == 0) * (1/2)
  
  # pred. which teams will advance to semifinials
  semi_y = y0_po[1:12,]
  semi_x = x0_po[1:12,]
  semi_win = (semi_y >0)*1 + (semi_y <0)*(-1) 
  semi_win[is.na(semi_win)] = 0
  semi_pro = cbind(sum(semi_win[1:3]),sum(semi_win[4:6]),sum(semi_win[7:9]),sum(semi_win[10:12]))
  
  semi_pred_pro = cbind(sum(pred_po[1:3]),sum(pred_po[4:6]),sum(pred_po[7:9]),sum(pred_po[10:12]))
  semi = sum(semi_pro*semi_pred_pro>0)*1 + sum(semi_pro*semi_pred_pro==0)*(1/2)
  
  # pred. which teams will advance to the finial
  final_y = y0_po[13:18,]
  final_x = x0_po[13:18,]
  final_win = (final_y >0)*1 + (final_y <0)*(-1) 
  final_win[is.na(final_win)] = 0
  final_pro = cbind(sum(final_win[1:3]),sum(final_win[4:6]))
  
  final_pred_pro = cbind(sum(pred_po[13:15]),sum(pred_po[16:18]))
  final=sum(final_pro*final_pred_pro>0)*1+sum(final_pro*final_pred_pro==0)*(1/2)
  
  # pred. which team will win the championship
  champ_y = y0_po[19:21,]
  champ_x = x0_po[19:21,]
  champ_win = (champ_y >0)*1 + (champ_y <0)*(-1) 
  champ_win[is.na(champ_win)] = 0
  champ_pro = cbind(sum(champ_win))
  
  champ_pred_pro = cbind(sum(pred_po[19:21]))
  champ=sum(champ_pro*champ_pred_pro>0)*1+sum(champ_pro*champ_pred_pro==0)*(1/2)
  
  pred_inf = cbind(div,K-clu_opt+1,beta_Est[clu_opt,2],beta_Est[clu_opt,3],beta_Est[clu_opt,4],semi,final,champ)
  pred_inf
}



