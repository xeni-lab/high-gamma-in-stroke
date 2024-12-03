match_moverate_stroke_control <- function(dir) {

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
plot_s_m1=matrix(nrow=1000, ncol=14)
plot_c_m1=matrix(nrow=1000, ncol=15)
plot_s_pmc=matrix(nrow=1000, ncol=14)
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
  allsubjects_s=as.numeric(str_sub(unique(data$subj_id[data$group=='stroke']),2,3))
  allsubjects_c=as.numeric(str_sub(unique(data$subj_id[data$group=='control']),2,3))
  
  # get movement speed of participants, seperated for different groups and performing hand within the groups
  moverate=data[data$group=='stroke' & data$perfhand=='r', 'moverate']
  trialidx=1:length(moverate) # number of trial in the original order
  data_sr=data.frame(moverate, trialidx)
  data_sr=sample_n(data_sr,nrow(data_sr)) # randomize order of trials to ensure random selection later on
  
  moverate=data[data$group=='stroke' & data$perfhand=='l', c('moverate')]
  trialidx=1:length(moverate)
  data_sl=data.frame(moverate, trialidx)
  data_sl=sample_n(data_sl,nrow(data_sl))
  
  moverate=data[data$group=='control' & data$perfhand=='r', c('moverate')]
  trialidx=1:length(moverate)
  data_cr=data.frame(moverate, trialidx)
  data_cr=sample_n(data_cr,nrow(data_cr))
  
  moverate=data[data$group=='control' & data$perfhand=='l', c('moverate')]
  trialidx=1:length(moverate)
  data_cl=data.frame(moverate, trialidx)
  data_cl=sample_n(data_cl,nrow(data_cl))
  
  # initialize variables
  indx=1;
  indx2=1;
  nomatchr=numeric();
  indicesr=data.frame(matrix(ncol=2, nrow=0))
  colnames(indicesr)=c('strokes','controls')
  
  # loops through trials of strokes with right performing hand, finds all trials within 
  # +-1 % of movement speed from controls with right performing hand, selects a random 
  # trial of the matching trials
  for (i in 1:nrow(data_sr)) {
    dum=data_sr$moverate[i];
    dum2=data_cr$moverate>dum-0.01*dum & data_cr$moverate<dum+0.01*dum;
    dum3=which(dum2==TRUE);
    if (length(dum3)>0) {
      x=sample(1:length(dum3),1)
      newrow=c(data_sr$trialidx[i],data_cr$trialidx[dum3[x]]);
      indicesr[indx,]=newrow # saves trialnumbers of matching trials
      data_cr[dum3[x],]=NA; # matching trial that was selected cannot be found again
      indx=indx+1;
    }
    else {
      nomatchr[indx2]=data_sr$trialidx[i]
      indx2=indx2+1;
    }
  }
  
  # initialize variables
  indx=1;
  indx2=1;
  nomatchl=numeric();
  indicesl=data.frame(matrix(ncol=2, nrow=0))
  colnames(indicesl)=c('strokes','controls')
  
  # loops through trials of strokes with left performing hand, finds all trials within 
  # +-1 % of movement speed from controls with left performing hand, selects a random 
  # trial of the matching trials
  for (i in 1:nrow(data_sl)) {
    dum=data_sl$moverate[i];
    dum2=data_cl$moverate>dum-0.01*dum & data_cl$moverate<dum+0.01*dum;
    dum3=which(dum2==TRUE);
    if (length(dum3)>0) {
      x=sample(1:length(dum3),1)
      newrow=c(data_sl$trialidx[i],data_cl$trialidx[dum3[x]]);
      indicesl[indx,]=newrow # saves trialnumbers of matching trials
      data_cl[dum3[x],]=NA; # matching trial that was selected cannot be found again
      indx=indx+1;
    }
    else {
      nomatchl[indx2]=data_sl$trialidx[i]
      indx2=indx2+1;
    }
  }
  
  # print and save number of matching trials in this iteration
  print(paste('Iteration:', as.character(j), 'Number of matching trials:', as.character(nrow(indicesr)+nrow(indicesl))))
  ntrials[j]=nrow(indicesr)+nrow(indicesl)
  
  ## compute group difference in high gamma power in M1 with matching trials only 
  region='hgp_m1'
  
  # select matching trials out of all trials
  data_sr=data[data$group=='stroke' & data$perfhand=='r', c(region,'group','subj_id')]
  data_sl=data[data$group=='stroke' & data$perfhand=='l', c(region,'group','subj_id')]
  data_cr=data[data$group=='control' & data$perfhand=='r', c(region,'group','subj_id')]
  data_cl=data[data$group=='control' & data$perfhand=='l', c(region,'group','subj_id')]
  
  data_sr=data_sr[indicesr$strokes,]
  data_cr=data_cr[indicesr$controls,]
  data_sl=data_sl[indicesl$strokes,]
  data_cl=data_cl[indicesl$controls,]
  
  data_all=rbind(data_sr,data_cr,data_sl, data_cl)
  
  # get subjectmean of high gamma power in M1 in selected trials for each group
  subjmean_s=tapply(data_all$hgp_m1[data_all$group=='stroke'], data_all$subj_id[data_all$group=='stroke'],'mean')
  subjmean_c=tapply(data_all$hgp_m1[data_all$group=='control'], data_all$subj_id[data_all$group=='control'],'mean')
  
  # t-test between groups with the alternative hypothesis that control participants 
  # have less high gamma power than stroke patients
  p=t.test(subjmean_s, subjmean_c, alternative='less')
  
  pvalue_m1[j]=p$p.value
  tvalue_m1[j]=p$statistic[[1]]
  df_m1[j]=p$parameter[[1]]
  
  # save movement speed values for plotting
  string=dimnames(subjmean_s) 
  string=string[[1]]
  subjects=as.numeric(str_sub(string,2,3))
  diffindx=which(!(allsubjects_s %in% subjects))
  plot_s=as.numeric(subjmean_s)
  if (length(diffindx)>0) {
    plot_s=insert_NA(plot_s, diffindx) # insert NA for participants with no selected matching trials
  }
  
  string=dimnames(subjmean_c)
  string=string[[1]]
  subjects=as.numeric(str_sub(string,2,3))
  diffindx=which(!(allsubjects_c %in% subjects))
  plot_c=as.numeric(subjmean_c)
  if (length(diffindx)>0) {
    plot_c=insert_NA(plot_c, diffindx)
  }
  
  plot_s_m1[j,]=plot_s
  plot_c_m1[j,]=plot_c
  
  ## compute group difference in high gamma power in PMC with matching trials only 
  region='hgp_pmc'
  
  # select matching trials out of all trials
  data_sr=data[data$group=='stroke' & data$perfhand=='r', c(region,'group','subj_id')]
  data_sl=data[data$group=='stroke' & data$perfhand=='l', c(region,'group','subj_id')]
  data_cr=data[data$group=='control' & data$perfhand=='r', c(region,'group','subj_id')]
  data_cl=data[data$group=='control' & data$perfhand=='l', c(region,'group','subj_id')]
  
  data_sr=data_sr[indicesr$strokes,]
  data_cr=data_cr[indicesr$controls,]
  data_sl=data_sl[indicesl$strokes,]
  data_cl=data_cl[indicesl$controls,]
  data_all=rbind(data_sr,data_cr,data_sl, data_cl)
  
  # get subjectmean of high gamma power in PMC in selected trials for each group
  subjmean_s=tapply(data_all$hgp_pmc[data_all$group=='stroke'], data_all$subj_id[data_all$group=='stroke'],'mean')
  subjmean_c=tapply(data_all$hgp_pmc[data_all$group=='control'], data_all$subj_id[data_all$group=='control'],'mean')
  
  # t-test between groups with the alternative hypothesis that stroke patients 
  # have less high gamma power than control participants
  p=t.test(subjmean_s, subjmean_c, alternative='less')
  
  pvalue_pmc[j]=p$p.value
  tvalue_pmc[j]=p$statistic[[1]]
  df_pmc[j]=p$parameter[[1]]
  
  # save movement speed values for plotting
  string=dimnames(subjmean_s)
  string=string[[1]]
  subjects=as.numeric(str_sub(string,2,3))
  diffindx=which(!(allsubjects_s %in% subjects))
  plot_s=as.numeric(subjmean_s)
  if (length(diffindx)>0) {
    plot_s=insert_NA(plot_s, diffindx)
  }
  
  string=dimnames(subjmean_c)
  string=string[[1]]
  subjects=as.numeric(str_sub(string,2,3))
  diffindx=which(!(allsubjects_c %in% subjects))
  plot_c=as.numeric(subjmean_c)
  if (length(diffindx)>0) {
    plot_c=insert_NA(plot_c, diffindx)
  }
  
  plot_s_pmc[j,]=plot_s
  plot_c_pmc[j,]=plot_c
  
} # end of loop through 1000 iterations
result=list("pvalue_m1" = pvalue_m1, 
            "pvalue_pmc" = pvalue_pmc, 
            "plot_s_m1" = plot_s_m1,
            "plot_s_pmc" = plot_s_pmc,
            "plot_c_m1" = plot_c_m1,
            "plot_c_pmc" = plot_c_pmc,
            "ntrials" = ntrials,
            "tvalue_m1" = tvalue_m1,
            "tvalue_pmc" = tvalue_pmc,
            "df_m1" = df_m1,
            "df_pmc" = df_pmc)
return(result)
}
