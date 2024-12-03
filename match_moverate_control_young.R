match_moverate_control_young <- function(dir) {

# packages needed
library(dplyr)
library(stringr)

# initialize variables
pvalue_m1=numeric(1000)
pvalue_pmc=numeric(1000)
tvalue_m1=numeric(1000)
tvalue_pmc=numeric(1000)
df_m1=numeric(1000)
df_pmc=numeric(1000)
ntrials=numeric()
plot_y_m1=matrix(nrow=1000, ncol=29)
plot_c_m1=matrix(nrow=1000, ncol=15)
plot_y_pmc=matrix(nrow=1000, ncol=29)
plot_c_pmc=matrix(nrow=1000, ncol=15)

# function that inserts NA in a vector at specified positions
insert_NA <- function(vec, diffindx) {
  # Calculate the size of the resulting vector
  result_size <- length(vec) + length(diffindx)
  
  # Create a vector of NAs with the calculated size
  modified_vec <- rep(NA, result_size)
  
  # Loop through the original vector and insert values
  k <- 1
  for (i in 1:result_size) {
    if (i %in% diffindx) {
      modified_vec[i] <- NA
    } else {
      modified_vec[i] <- vec[k]
      k <- k + 1
    }
  }
  
  return(modified_vec)
}

# loop trough 1000 iterations of the matching process
for (j in 1:length(pvalue_m1)) {
  
  data=read.csv(paste0(dir, "data_triallevel_publication.csv"))
  data$group <- factor(data$group,levels = c('young','control','stroke'),ordered = TRUE)
  # get ID numbers of stroke and control participants
  allsubjects_y=as.numeric(str_sub(unique(data$subj_id[data$group=='young']),2,3))
  allsubjects_c=as.numeric(str_sub(unique(data$subj_id[data$group=='control']),2,3))
  
  # get movement speed of participants, seperated for different groups, only 
  # for performance with left hand as all young participants performed the task with 
  # the left hand
  moverate=data[data$group=='young' & data$perfhand=='l', c('moverate')]
  trialidx=1:length(moverate)
  data_y=data.frame(moverate, trialidx)
  data_y=sample_n(data_y,nrow(data_y))
  
  moverate=data[data$group=='control' & data$perfhand=='l', c('moverate')]
  trialidx=1:length(moverate)
  data_c=data.frame(moverate, trialidx)
  data_c=sample_n(data_c,nrow(data_c))
  
  # initialize variables
  indx=1;
  indx2=1;
  nomatch=numeric();
  indices=data.frame(matrix(ncol=2, nrow=0))
  colnames(indices)=c('young','controls')
  
  # loops through trials of young participants finds all trials within +-1 % 
  # of movement speed from controls with left performing hand, selects a random 
  # trial of the matching trials
  for (i in 1:nrow(data_y)) {
    dum=data_y$moverate[i];
    dum2=data_c$moverate>dum-0.01*dum & data_c$moverate<dum+0.01*dum;
    dum3=which(dum2==TRUE);
    if (length(dum3)>0) {
      x=sample(1:length(dum3),1)
      newrow=c(data_y$trialidx[i],data_c$trialidx[dum3[x]]);
      indices[indx,]=newrow # saves trialnumbers of matching trials
      data_c[dum3[x],]=NA; # matching trial that was selected cannot be found again
      indx=indx+1;
    }
    else {
      nomatch[indx2]=data_y$trialidx[i]
      indx2=indx2+1;
    }
  }
  
  # print and save number of matching trials in this iteration
  print(paste('Iteration:', as.character(j), 'Number of matching trials:', as.character(nrow(indices))))
  ntrials[j]=nrow(indices)
  
  ## compute group difference in high gamma power in M1 with matching trials only 
  region='hgp_m1'
  
  # select matching trials out of all trials
  data_y=data[data$group=='young' & data$perfhand=='l', c(region,'group','subj_id')]
  data_c=data[data$group=='control' & data$perfhand=='l', c(region,'group','subj_id')]
  data_y=data_y[indices$young,]
  data_c=data_c[indices$controls,]
  data_all=rbind(data_y, data_c)
  
  # get subjectmean of high gamma power in M1 in selected trials for each group
  subjmean_y=tapply(data_all$hgp_m1[data_all$group=='young'], data_all$subj_id[data_all$group=='young'],'mean')
  subjmean_c=tapply(data_all$hgp_m1[data_all$group=='control'], data_all$subj_id[data_all$group=='control'],'mean')
  
  # t-test between groups with the alternative hypothesis that control participants 
  # have less high gamma power than young participants
  p=t.test(subjmean_c, subjmean_y, alternative='less')
  
  pvalue_m1[j]=p$p.value
  tvalue_m1[j]=p$statistic[[1]]
  df_m1[j]=p$parameter[[1]]
  
  # save movement speed values for plotting
  string=dimnames(subjmean_y)
  string=string[[1]]
  subjects=as.numeric(str_sub(string,2,3))
  diffindx=which(!(allsubjects_y %in% subjects))
  plot_y=as.numeric(subjmean_y)
  if (length(diffindx)>0) {
    plot_y=insert_NA(plot_y, diffindx) # # insert NA for participants with no selected matching trials
  }
  
  string=dimnames(subjmean_c)
  string=string[[1]]
  subjects=as.numeric(str_sub(string,2,3))
  diffindx=which(!(allsubjects_c %in% subjects))
  plot_c=as.numeric(subjmean_c)
  if (length(diffindx)>0) {
    plot_c=insert_NA(plot_c, diffindx)
  }
  
  plot_y_m1[j,]=plot_y
  plot_c_m1[j,]=plot_c
  
  ## compute group difference in high gamma power in PMC with matching trials only
  region='hgp_pmc'
  
  # select matching trials out of all trials
  data_y=data[data$group=='young' & data$perfhand=='l', c(region,'group','subj_id')]
  data_c=data[data$group=='control' & data$perfhand=='l', c(region,'group','subj_id')]
  data_y=data_y[indices$young,]
  data_c=data_c[indices$controls,]
  data_all=rbind(data_y, data_c)
  
  # get subjectmean of high gamma power in M1 in selected trials for each group
  subjmean_y=tapply(data_all$hgp_pmc[data_all$group=='young'], data_all$subj_id[data_all$group=='young'],'mean')
  subjmean_c=tapply(data_all$hgp_pmc[data_all$group=='control'], data_all$subj_id[data_all$group=='control'],'mean')
  
  # t-test between groups with the alternative hypothesis that control participants 
  # have less high gamma power than stroke patients
  p=t.test(subjmean_c, subjmean_y, alternative='less')
  
  pvalue_pmc[j]=p$p.value
  tvalue_pmc[j]=p$statistic[[1]]
  df_pmc[j]=p$parameter[[1]]
  
  # save movement speed values for plotting
  string=dimnames(subjmean_y)
  string=string[[1]]
  subjects=as.numeric(str_sub(string,2,3))
  diffindx=which(!(allsubjects_y %in% subjects))
  plot_y=as.numeric(subjmean_y)
  if (length(diffindx)>0) {
    plot_y=insert_NA(plot_y, diffindx)
  }
  string=dimnames(subjmean_c)
  string=string[[1]]
  subjects=as.numeric(str_sub(string,2,3))
  diffindx=which(!(allsubjects_c %in% subjects))
  plot_c=as.numeric(subjmean_c)
  if (length(diffindx)>0) {
    plot_c=insert_NA(plot_c, diffindx)
  }
  
  plot_y_pmc[j,]=plot_y
  plot_c_pmc[j,]=plot_c
  
} # end of loop through 1000 iterations
result=list("pvalue_m1" = pvalue_m1, 
            "pvalue_pmc" = pvalue_pmc, 
            "plot_y_m1" = plot_y_m1,
            "plot_y_pmc" = plot_y_pmc,
            "plot_c_m1" = plot_c_m1,
            "plot_c_pmc" = plot_c_pmc,
            "ntrials" = ntrials,
            "tvalue_m1" = tvalue_m1,
            "tvalue_pmc" = tvalue_pmc,
            "df_m1" = df_m1,
            "df_pmc" = df_pmc)
return(result)
}
