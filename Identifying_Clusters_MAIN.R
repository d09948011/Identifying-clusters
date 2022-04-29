# ========================================================================== #
## ------------------------------------------------------------------------ ##
## Code for An Effective Method for Identifying Clusters of Robot Strengths ##
## ---------- Jen-Chieh Teng, Chin-Tsang Chinag, and Alvin Lim ------------ ##
## ------------------------------------------------------------------------ ##
# ========================================================================== #


# load auxillary functions & libraries for analysis 
source("Identifying_Clusters_AUX.R")

# LCT and TCL method
Classification.fn = function(x,y,m_0){
  
  # when the number of cluster is K
  b_info = matrix(0,K,K-1)
  b_info[1:K,K-1] = bhat0 = beta_esti.fn(diag(1,K,K))
  fit = hclust(dist(bhat0),method="centroid")
  Est = matrix(0,K-1,K+4)

  Est[1,1:4] = c(K,esti.fn(diag(1,K,K)))
  Est[1,5:(K+4)] = cutree(fit,K)

  for(w in (K-1):2){ # compute for the number of cluster is K-1, ..., 2

    fit = hclust(dist(bhat0),method="centroid")
    grp0 = grp = cutree(fit,k=w)
    INF = outer(grp,sort(unique(grp)),FUN="==") *1 
    b_info[1:K,w-1] = bhat0 = beta_esti.fn(INF)
    Est[K-w+1,] = c(w,esti.fn(INF),grp)
    c_list = c(which(colSums(INF)>m_0,arr.ind=TRUE))    # collect the number of the group whose size is greater than m_0
    c_list_l = length(c_list)
    
    diff = cri = 2 
    re_beta = c(bhat0)  
    while(diff > exp(-10) && c_list_l > 0 && cri != 1){ # the criteria for iterating non-hierarchical classification
      new_beta = old_beta = c(bhat0)
      for (s1 in c_list){    # the s1-group
        re_1 = which(c(grp0)==s1,arr.ind=TRUE)   # position
        
        for (s2 in re_1){     # the s2 position
          grp = c(grp0)
          grp[s2] = 0      # let one of robot in s1-group become group 0 
          INF2 = (outer(grp,sort(unique(grp)),FUN="==")*1)
          new_beta[s2] = beta_esti.fn(INF2)[s2]
        }
      }
      
      fit = hclust(dist(c(new_beta)),method="centroid") # perform the centroid linkage clustering method on reestimates
      grp0 = grp = cutree(fit,k=w)
      
      INF = outer(grp,sort(unique(grp)),FUN="==")*1
      bhat0 = beta_esti.fn(INF)
      re_beta = rbind(re_beta,c(bhat0))
      cri = sum(duplicated(re_beta))   
      diff = sum((bhat0 - old_beta)^2)/M
      
      c_list = c(which(colSums(INF)>m_0, arr.ind=TRUE))   
      c_list_l = length(c_list)
      b_info[1:K,w-1] = bhat0
    }
    
    Est[K-w+1,] = c(w,esti.fn(INF),grp)
  }
  Est 
}



