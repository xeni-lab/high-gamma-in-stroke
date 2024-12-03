#### necessary packages ####
install.packages('rstatix')
install.packages('lme4')
install.packages('stringr')
install.packages('dplyr')
install.packages('ggpubr')
install.packages('ggh4x')
install.packages('tidyr')
install.packages('ggsignif')
install.packages('ordinal')
install.packages('effects')
install.packages('MASS')

#### set path to folder with data ####
dir="path to data"

#### Results ####
#### participants: mean and sd of participant age and mmst ####
data=read.csv(paste0(dir, "data_subjectlevel.csv"))

mean(data$age[data$group=='stroke'])
sd(data$age[data$group=='stroke'])
mean(data$age[data$group=='control'])
sd(data$age[data$group=='control'])
mean(data$age[data$group=='young'])
sd(data$age[data$group=='young'])

mean(data$mmst[data$group=='stroke'])
mean(data$mmst[data$group=='control'])

#### clinical scores and structural imaging  ####
library(rstatix)
data=read.csv(paste0(dir, "data_subjectlevel.csv"))
data$group <- factor(data$group,levels = c('stroke','control','young'),ordered = TRUE)
my_comparisons=list(c('stroke', 'control'))

# mean uefm of stroke patients
mean(data$uefm[data$group=='stroke'])

# group differences in clincal scores
wilcox_test(data, uefm ~ group, comparisons=my_comparisons)
wilcox_test(data, arat ~ group, comparisons=my_comparisons)
t_test(data, nhpt ~ group, comparisons=my_comparisons)
t_test(data, bbt ~ group, comparisons=my_comparisons)
wilcox_test(data, gripforce_ratio ~ group, comparisons=my_comparisons)
t_test(data, keygripforce_ratio ~ group, comparisons=my_comparisons)
t_test(data, cst_ratio ~ group, comparisons=my_comparisons)

#### improvement in movement rate ####
library(rstatix)
data=read.csv(paste0(dir, "data_subjectlevel.csv"))
my_comparisons=list(c('stroke', 'control'), c('control','young'))
data$group <- factor(data$group,levels = c('stroke','control','young'),ordered = TRUE)

# group differences in movement rate
t_test(data, moverate ~ group, comparisons=my_comparisons, p.adjust.method='fdr')

# differences in movement rate between block 1 and block 6 in each group
dum=t.test(data$moverate_b1[data$group=='stroke'], data$moverate_b6[data$group=='stroke'], paired=TRUE)
p=dum$p.value
dum=t.test(data$moverate_b1[data$group=='control'], data$moverate_b6[data$group=='control'], paired=TRUE)
p[2]=dum$p.value
dum=t.test(data$moverate_b1[data$group=='young'], data$moverate_b6[data$group=='young'], paired=TRUE)
p[3]=dum$p.value

p.adjust(p, method='fdr')

# mean and sd of blockwise improvement in each group
mean(data$blockwise_improvement[data$group=='stroke'])
sd(data$blockwise_improvement[data$group=='stroke'])
mean(data$blockwise_improvement[data$group=='control'])
sd(data$blockwise_improvement[data$group=='control'])
mean(data$blockwise_improvement[data$group=='young'])
sd(data$blockwise_improvement[data$group=='young'])

# group difference in blockwise improvement
t_test(data, blockwise_improvement ~ group, comparisons=my_comparisons, p.adjust.method='fdr')

#### movement related high gamma power ####
library(rstatix)

data=read.csv(paste0(dir, "data_subjectlevel.csv"))
data$group <- factor(data$group,levels = c('stroke','control','young'),ordered = TRUE)
my_comparisons=list(c('stroke', 'control'), c('control','young'))

# group differences in high gamma power in M1 and PMC
dum=t_test(data, hgp_m1 ~ group, comparisons=my_comparisons, p.adjust.method='fdr')
p=dum$p
dum=t_test(data, hgp_pmc ~ group, comparisons=my_comparisons, p.adjust.method='fdr')
p[3:4]=dum$p

p.adjust(p, method='fdr')

# group difference in high gamma power between stroke patients
# and control participants matched for movement rate

# results of analysis used for publication: load(paste0(dir, "match_sc.Rdata"))

# function for the matching process, takes a few minutes
match_sc=match_moverate_stroke_control(dir)

# get median p-value of 1000 iterations for M1 and PMC, adjust for 2 comparisons
p=median(match_sc$pvalue_m1)
p[2]=median(match_sc$pvalue_pmc)

p.adjust(p, method='fdr')

# group difference in high gamma power between young participants
# and control participants matched for movement rate

# results of analysis used for publication: load(paste0(dir, "match_cy.Rdata"))

# function for the matching process, takes a few minutes
match_cy=match_moverate_control_young(dir)

# get median p-value of 1000 iterations for M1 and PMC, adjust for 2 comparisons
p=median(match_cy$pvalue_m1)
p[2]=median(match_cy$pvalue_pmc)

p

# effect of group on high-gamma power in a linear mixed-effects model including all groups
library(lme4)
data=read.csv(paste0(dir, "data_triallevel.csv"))

mdl=lmer(log(100+hgp_m1) ~ moverate + group + perfhand + (1+moverate|subj_id), data=data, REML=FALSE)
mdlnull=lmer(log(100+hgp_m1) ~ moverate + perfhand + (1+moverate|subj_id), data=data, REML=FALSE)
anova(mdl, mdlnull)

mdl=lmer(log(100+hgp_pmc) ~ moverate + group + perfhand + (1+moverate|subj_id), data=data, REML=FALSE)
mdlnull=lmer(log(100+hgp_pmc) ~ moverate + perfhand + (1+moverate|subj_id), data=data, REML=FALSE)
anova(mdl, mdlnull)

#### relation of high gamma power, structure and clinical scores ####
data=read.csv(paste0(dir, "data_subjectlevel.csv"))

# models of high gamma power and the integrity of the cst
mdl=lm(hgp_m1 ~ cst_ratio + perfhand, data=data[data$group=='stroke' & !is.na(data$cst_ratio),])
mdlnull=lm(hgp_m1 ~ perfhand, data=data[data$group=='stroke' & !is.na(data$cst_ratio),])
anova(mdl, mdlnull)
mdl=lm(hgp_pmc ~ cst_ratio + perfhand, data=data[data$group=='stroke' & !is.na(data$cst_ratio),])
mdlnull=lm(hgp_pmc ~ perfhand, data=data[data$group=='stroke' & !is.na(data$cst_ratio),])
anova(mdl, mdlnull)

# models of high gamma power and clinical scores
mdl=lm(hgp_m1 ~ uefm + perfhand, data=data[data$group=='stroke',])
mdlnull=lm(hgp_m1 ~ perfhand, data=data[data$group=='stroke',])
anova(mdl, mdlnull)
mdl=lm(hgp_pmc ~ uefm + perfhand, data=data[data$group=='stroke',])
mdlnull=lm(hgp_pmc ~ perfhand, data=data[data$group=='stroke',])
anova(mdl, mdlnull)

mdl=lm(hgp_m1 ~ arat + perfhand, data=data[data$group=='stroke',])
mdlnull=lm(hgp_m1 ~ perfhand, data=data[data$group=='stroke',])
anova(mdl, mdlnull)
mdl=lm(hgp_pmc ~ arat + perfhand, data=data[data$group=='stroke',])
mdlnull=lm(hgp_pmc ~ perfhand, data=data[data$group=='stroke',])
anova(mdl, mdlnull)

mdl=lm(hgp_m1 ~ nhpt + perfhand, data=data[data$group=='stroke',])
mdlnull=lm(hgp_m1 ~ perfhand, data=data[data$group=='stroke',])
anova(mdl, mdlnull)
mdl=lm(hgp_pmc ~ nhpt + perfhand, data=data[data$group=='stroke',])
mdlnull=lm(hgp_pmc ~ perfhand, data=data[data$group=='stroke',])
anova(mdl, mdlnull)

mdl=lm(hgp_m1 ~ bbt + perfhand, data=data[data$group=='stroke',])
mdlnull=lm(hgp_m1 ~ perfhand, data=data[data$group=='stroke',])
anova(mdl, mdlnull)
mdl=lm(hgp_pmc ~ bbt + perfhand, data=data[data$group=='stroke',])
mdlnull=lm(hgp_pmc ~ perfhand, data=data[data$group=='stroke',])
anova(mdl, mdlnull)

mdl=lm(hgp_m1 ~ gripforce_ratio + perfhand, data=data[data$group=='stroke',])
mdlnull=lm(hgp_m1 ~ perfhand, data=data[data$group=='stroke',])
anova(mdl, mdlnull)
mdl=lm(hgp_pmc ~ gripforce_ratio + perfhand, data=data[data$group=='stroke',])
mdlnull=lm(hgp_pmc ~ perfhand, data=data[data$group=='stroke',])
anova(mdl, mdlnull)

mdl=lm(hgp_m1 ~ keygripforce_ratio + perfhand, data=data[data$group=='stroke',])
mdlnull=lm(hgp_m1 ~ perfhand, data=data[data$group=='stroke',])
anova(mdl, mdlnull)
mdl=lm(hgp_pmc ~ keygripforce_ratio + perfhand, data=data[data$group=='stroke',])
mdlnull=lm(hgp_pmc ~ perfhand, data=data[data$group=='stroke',])
anova(mdl, mdlnull)

#### relation of high gamma power and movement rate ####
library(lme4)

data=read.csv(paste0(dir, "data_triallevel.csv"))

# linear mixed-effects models with non-significant interaction between movement rate and group
mdl=lmer(log(100+hgp_m1) ~ moverate * group + perfhand + (1+moverate|subj_id), data=data, REML=FALSE)
mdlnull=lmer(log(100+hgp_m1) ~ moverate + group + perfhand + (1+moverate|subj_id), data=data, REML=FALSE)
anova(mdl, mdlnull)
mdl=lmer(log(100+hgp_pmc) ~ moverate * group + perfhand + (1+moverate|subj_id), data=data, REML=FALSE)
mdlnull=lmer(log(100+hgp_pmc) ~ moverate + group + perfhand + (1+moverate|subj_id), data=data, REML=FALSE)
anova(mdl, mdlnull)

# linear mixed-effects model for the effect of movement rate on high-gamma power
mdl=lmer(log(100+hgp_m1) ~ moverate + group + perfhand + (1+moverate|subj_id), data=data, REML=FALSE)
mdlnull=lmer(log(100+hgp_m1) ~ group + perfhand + (1+moverate|subj_id), data=data, REML=FALSE)
dum=anova(mdl, mdlnull)
p=dum$`Pr(>Chisq)`[2]
mdl=lmer(log(100+hgp_pmc) ~ moverate + group + perfhand + (1+moverate|subj_id), data=data, REML=FALSE)
mdlnull=lmer(log(100+hgp_pmc) ~ group + perfhand + (1+moverate|subj_id), data=data, REML=FALSE)
dum=anova(mdl, mdlnull)
p[2]=dum$`Pr(>Chisq)`[2]

p.adjust(p, method='fdr')

# participant-level linear models for the effect of movement rate on high-gamma power
data=read.csv(paste0(dir, "data_subjectlevel.csv"))

mdl=lm(hgp_m1 ~ moverate + group + perfhand, data=data)
mdlnull=lm(hgp_m1 ~ group + perfhand, data=data)
dum=anova(mdl, mdlnull)
p=dum$`Pr(>F)`[2]
mdl=lm(hgp_pmc ~ moverate + group + perfhand, data=data)
mdlnull=lm(hgp_pmc ~ group + perfhand, data=data)
dum=anova(mdl, mdlnull)
p[2]=dum$`Pr(>F)`[2]

p.adjust(p, method='fdr')

#### high gamma peak frequency ####
library(rstatix)
data=read.csv(paste0(dir, "data_subjectlevel.csv"))
data$group <- factor(data$group,levels = c('stroke','control','young'),ordered = TRUE)
my_comparisons=list(c('stroke', 'control'), c('control','young'))

# mean peak frequency in the different groups
mean(data$peakfreq_m1[data$group=='young'])
mean(data$peakfreq_m1[data$group=='control'])
mean(data$peakfreq_m1[data$group=='stroke'])

# group differences in peak frequency
dum=wilcox_test(data, peakfreq_m1 ~ group, comparisons=my_comparisons, p.adjust.method='fdr')
p=dum$p
dum=wilcox_test(data, peakfreq_pmc ~ group, comparisons=my_comparisons, p.adjust.method='fdr')
p[3:4]=dum$p

p.adjust(p, method='fdr')

# linear models for the effect of movement rate on high gamma peak frequency
mdl=lm(peakfreq_m1 ~ moverate + group + perfhand, data=data)
mdlnull=lm(peakfreq_m1 ~ group + perfhand, data=data)
anova(mdl, mdlnull)
mdl=lm(peakfreq_pmc ~ moverate + group + perfhand, data=data)
mdlnull=lm(peakfreq_pmc ~ group + perfhand, data=data)
anova(mdl, mdlnull)

#### number of high gamma peaks ####
library(ordinal)
library(MASS)
library(rstatix)
library(effects)
data=read.csv(paste0(dir, "data_subjectlevel.csv"))
data$group <- factor(data$group,levels = c('stroke','control','young'),ordered = TRUE)
my_comparisons=list(c('stroke', 'control'), c('control','young'))

# group differences in number of high gamma peaks
dum=wilcox_test(data, npeaks_m1 ~ group, comparisons=my_comparisons, p.adjust.method='fdr')
p=dum$p
dum=wilcox_test(data, npeaks_pmc ~ group, comparisons=my_comparisons, p.adjust.method='fdr')
p[3:4]=dum$p

p.adjust(p, method='fdr')

## models between number of high gamma peaks, movement rate and group
data$npeaks_m1=as.factor(data$npeaks_m1)
data$npeaks_pmc=as.factor(data$npeaks_pmc)
data$group=as.factor(data$group)

# effect of movement rate
mdl=clm(npeaks_m1 ~ moverate + group + perfhand, data=data)
mdlnull=clm(npeaks_m1 ~ group + perfhand, data=data)
dum=anova(mdl, mdlnull)
p=dum$`Pr(>Chisq)`[2]
mdl=clm(npeaks_pmc ~ moverate + group + perfhand, data=data)
mdlnull=clm(npeaks_pmc ~ group + perfhand, data=data)
dum=anova(mdl, mdlnull)
p[2]=dum$`Pr(>Chisq)`[2]

p.adjust(p, method='fdr')

# effect of group
mdl=clm(npeaks_m1 ~ moverate + group + perfhand, data=data)
mdlnull=clm(npeaks_m1 ~ moverate + perfhand, data=data)
dum=anova(mdl, mdlnull)
p=dum$`Pr(>Chisq)`[2]
mdl=clm(npeaks_pmc ~ moverate + group + perfhand, data=data)
mdlnull=clm(npeaks_pmc ~ moverate + perfhand, data=data)
dum=anova(mdl, mdlnull)
p[2]=dum$`Pr(>Chisq)`[2]

p.adjust(p, method='fdr')

# cutoff movement rate for number of peaks in M1
mdl=clm(npeaks_m1 ~ moverate + group + perfhand, data=data)
eff=effect('moverate', mdl, xlevels=1000)
prob1=eff$prob[,1]
prob2=eff$prob[,2]

peak1=prob1[eff$x>1.25]
peak2=prob2[eff$x>1.25]
perfdum=eff$x[eff$x>1.25]
comp=peak1>peak2

i=1
while (comp[i]==FALSE) {
  i=i+1
}

perfdum[i] # cutoff movement rate for 1 peak highest probability

prob2=eff$prob[,2]
prob4=eff$prob[,4]

peak2=prob2[eff$x<1.5]
peak4=prob4[eff$x<1.5]
perfdum=eff$x[eff$x<1.5]
comp=peak4>peak2

i=1
while (comp[i]==TRUE) {
  i=i+1
}

perfdum[i] # cutoff movement rate for 4 peaks highest probability

#### relation of high gamma power and improvement ####
library(lme4)

# models of high gamma power and trial-wise improvement
data=read.csv(paste0(dir, "data_triallevel.csv"))

mdl=lmer(log(100+hgp_m1) ~ trialwise_improvement + moverate + group + perfhand + (1+moverate|subj_id),data=data[!is.na(data$trialwise_improvement),], REML=FALSE)
mdlnull=lmer(log(100+hgp_m1) ~ group + moverate + perfhand + (1+moverate|subj_id),data=data[!is.na(data$trialwise_improvement),], REML=FALSE)
anova(mdl, mdlnull)
mdl=lmer(log(100+hgp_pmc) ~ trialwise_improvement + moverate + group + perfhand + (1|subj_id),data=data[!is.na(data$trialwise_improvement),], REML=FALSE)
mdlnull=lmer(log(100+hgp_pmc) ~ group + moverate + perfhand + (1|subj_id),data=data[!is.na(data$trialwise_improvement),], REML=FALSE)
anova(mdl, mdlnull)

# models of high gamma power and and block-wise improvement
data=read.csv(paste0(dir, "data_subjectlevel.csv"))

mdl=lm(hgp_m1 ~ blockwise_improvement + moverate_b1 + group + perfhand, data=data)
mdlnull=lm(hgp_m1 ~ moverate_b1 + group + perfhand, data=data)
anova(mdl,mdlnull)
mdl=lm(hgp_pmc ~ blockwise_improvement + moverate_b1 + group + perfhand, data=data)
mdlnull=lm(hgp_pmc ~ moverate_b1 + group + perfhand, data=data)
anova(mdl,mdlnull)

#### Figures: initialize variables and functions ####
library(ggpubr)
library(rstatix)
library(ggh4x)
library(tidyr)
library(ggsignif)
library(lme4)
library(ordinal)
library(MASS)
library(ggeffects)

my_theme = theme(
  panel.background = element_rect(fill = "white"),
  panel.border=element_rect(linewidth=0.5, fill=NA))

colstroke=rgb(0, 114, 199, maxColorValue=255)
colcontrol=rgb(243, 159, 24, maxColorValue=255)
colyoung='seagreen'
colall='black'

pvalues_stars <- function(p) {
  output=character(length(p))
  for (i in 1:length(p)) {
    if (p[i]>0.05) {
      output[i]='ns';
    } else if (p[i]<=0.05 & p[i]>=0.01) {
      output[i]='*';
    } else if (p[i]<0.01 & p[i]>=0.001) {
      output[i]='**';
    } else if (p[i]<0.001) {
      output[i]='***';
    } else {
      output[i]='Error'
    }
  }
  return(output)
}

ggploteffects <- function(mdl, terms, col,xlab, ylab, xlim=NULL, ylim=NULL,linewidth=1, rug=FALSE) {
  if (missing(mdl)) {
    print('no input model')
  }
  if (missing(terms)) {
    print('you need to specify the term(s) you are interested in')
  }
  library(ggeffects)
  library(ggpubr)
  eff=ggeffect(mdl, terms = terms)
  p=ggplot(eff, aes(x=x, y=predicted)) +
    geom_line(color=col,linewidth=linewidth)+
    geom_ribbon(aes(ymin=conf.low, ymax=conf.high),fill=col,alpha=0.2)+
    theme(panel.background = element_rect(fill = "white", color = "black", size = .5))+
    theme(panel.border=element_rect(size=.5,fill=NA))+
    coord_cartesian(xlim=xlim, ylim=ylim)
  
  if (!missing(xlab)) {
    p=p+labs(x=xlab)
  }
  if (!missing(ylab)) {
    p=p+labs(y=ylab)
  }
  
  if (rug==TRUE) {
    p=p+geom_rug(sides="b")
  }
  return(p)
}

#### figure 1B ####
data=read.csv(paste0(dir, "data_triallevel.csv"))

stroke=tapply(data$moverate[data$group=='stroke'], data$trial[data$group=='stroke'], mean)
control=tapply(data$moverate[data$group=='control'], data$trial[data$group=='control'], mean)
young=tapply(data$moverate[data$group=='young'], data$trial[data$group=='young'], mean)
ci_stroke=tapply(data$moverate[data$group=='stroke'], data$trial[data$group=='stroke'], sd)
ci_control=tapply(data$moverate[data$group=='control'], data$trial[data$group=='control'], sd)
ci_young=tapply(data$moverate[data$group=='young'], data$trial[data$group=='young'], sd)
df=data.frame(stroke, control, young, ci_stroke, ci_control, ci_young)

ggplot(df, aes(1:240, stroke)) +
  geom_point(color = colstroke, size=1, shape=15) +
  geom_point(aes(1:240, control), color = colcontrol, size=1, shape=17) +
  geom_point(aes(1:240, young), color = colyoung, size=1, shape=16) +
  geom_ribbon(aes(1:240, ymin = stroke - ci_stroke, ymax = stroke + ci_stroke),
              fill = colstroke, alpha = 0.2) +
  geom_ribbon(aes(1:240, ymin = control - ci_control, ymax = control + ci_control),
              fill = colcontrol, alpha = 0.2) +
  geom_ribbon(aes(1:240, ymin = young - ci_young, ymax = young + ci_young),
              fill = colyoung, alpha = 0.2) +
  my_theme+
  scale_x_continuous(breaks = seq(0, 240, by = 80), limits = c(-15, 255),
                     expand = c(0, 0))+ 
  coord_cartesian(ylim=c(0.6,2.2))+
  xlab('Trial')+
  ylab('Movement rate (1/s)')


#### figure 1C ####
data=read.csv(paste0(dir, "data_subjectlevel.csv"))
data$group <- factor(data$group,levels = c('young','control','stroke'),ordered = TRUE)

my_comparisons <- list( c("young", "control"), c("control", "stroke") )

stat.test <- data %>%
  t_test(moverate ~ group, comparisons=my_comparisons, p.adjust.method='fdr')
stat.test <- stat.test %>% mutate(y.position = c(2.7, 2.7))
stat.test$p.adj.signif=pvalues_stars(stat.test$p.adj)

ggplot(data, aes(x = group, y = moverate, color = group)) +
  geom_boxplot(outlier.shape=NA, size=0.5) +
  geom_jitter(position = position_jitterdodge(jitter.width = 2),size=1, shape=16) +
  scale_color_manual(values = c(colyoung, colcontrol, colstroke)) +
  stat_pvalue_manual(stat.test, label = 'p.adj.signif', size=6)+
  my_theme+
  coord_cartesian(ylim=c(0.6,2.8))+ 
  scale_y_continuous(breaks=c(0.5,1,1.5,2,2.5))+
  xlab('Group')+
  ylab('Movement rate (1/s)')

#### figure 1D ####
data=read.csv(paste0(dir, "data_subjectlevel.csv"))
data$group <- factor(data$group,levels = c('young','control','stroke'),ordered = TRUE)

my_comparisons <- list( c("young", "control"), c("control", "stroke") )

ggplot(data, aes(x = group, y = blockwise_improvement, color = group)) +
  geom_boxplot(outlier.shape=NA, size=0.5) +
  geom_jitter(position = position_jitterdodge(jitter.width = 2),size=1, shape=16) +
  scale_color_manual(values = c(colyoung, colcontrol, colstroke)) +
  my_theme+
  scale_y_continuous(breaks = seq(0, 40, by = 20))+
  xlab('Group')+
  ylab('Block-wise improvement (%)')

#### figure 2B ####
data=read.csv(paste0(dir, "data_subjectlevel.csv"))
data=data[data$group %in% c('stroke','control'),]
data$group <- factor(data$group,levels = c('control','stroke'),ordered = TRUE)

ggplot(data, aes(x = group, y = cst_ratio, color = group)) +
  geom_boxplot(outlier.shape=NA, size=0.5,width=0.5) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5),size=1, shape=16) +
  scale_color_manual(values = c(colcontrol, colstroke))+
  my_theme+
  xlab('Group')+
  ylab('CST-FA ratio')

#### figure 3B ####
data=read.csv(paste0(dir, "data_subjectlevel.csv"))
regions <- c("hgp_m1", "hgp_pmc")
df=data[,c(regions, "group", 'subj_id')]
df$group <- factor(df$group,levels = c('young','control','stroke'),ordered = TRUE)

# Reshape the data into long format
df_long <- df %>%
  tidyr::pivot_longer(cols = starts_with(c("hgp")), names_to = "region")

df_long$region <- factor(df_long$region,levels = regions,ordered = TRUE)

df_long <- df_long %>%
  mutate(region_group = interaction(group, region, sep = "_"))

p_values= df_long %>% 
  group_by(region) %>% 
  t_test(value ~ group,p.adjust.method='fdr',comparisons=list(c('young','control'),c('stroke','control')))
p_values$p.adj=p.adjust(p_values$p, method='fdr')
p_values$p.adj.signif=c('','','','') # stars according to p-values get inserted later on manually

y_values=rep(150,4)

# Create boxplot using ggplot with faceting
ggplot(df_long, aes(x = region_group, y = value, color = group)) +
  geom_boxplot(outlier.shape=NA, size=0.5) +
  geom_jitter(position = position_jitterdodge(jitter.width = 2),size=1, shape=16) +
  scale_color_manual(values = c(colyoung, colcontrol, colstroke)) +
  scale_x_manual(as.double(c(1,2.15,3.3,5,6.15,7.3)),labels=c('Young','Control','Stroke','Young','Control','Stroke'))+
  my_theme+
  geom_signif(y_position = y_values, 
              xmin = c(1,2.15,5,6.15), 
              xmax = c(2.15,3.3,6.15,7.3), 
              annotations = p_values$p.adj.signif, 
              color='black', 
              size=0.3, )+
  annotate('text',label='*',x=1.575,y=151.5,size=6)+
  annotate('text',label='*',x=2.65,y=151.5,size=6)+
  annotate('text',label='*',x=5.575,y=151.5,size=6)+
  annotate('text',label='*',x=6.65,y=151.5,size=6)+
  xlab('Group')+
  ylab('High gamma power (%)')

#### figure 3C ####
load(paste0(dir, "match_sc.Rdata"))

mean_m1=c(colMeans(match_sc$plot_s_m1, na.rm=TRUE), colMeans(match_sc$plot_c_m1, na.rm=TRUE))
group=character(29)
group[1:14]='stroke'
group[15:29]='control'

df=data.frame(group, mean_m1)

stat.test <- df %>%
  t_test(mean_m1 ~ group, p.adjust.method='fdr')
stat.test <- stat.test %>% mutate(y.position = c(91))
stat.test$p.adj.signif='*'

ggplot(df, aes(x=group, y=mean_m1, color = group), na.rm=TRUE) +
  geom_boxplot(outlier.shape=NA, size=0.5, width=0.4)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.4),size=1, shape=16) +
  scale_color_manual(values = c(colcontrol, colstroke)) +
  stat_pvalue_manual(stat.test, label = 'p.adj.signif', size=6)+
  my_theme+
  scale_y_continuous(limits=c(NA,97))+
  xlab('Group')+
  ylab('High gamma power (%)')

#### figure 4A ####
data=read.csv(paste0(dir, "data_triallevel.csv"))

mdl=lmer(log(100+hgp_m1) ~ moverate + group + perfhand + (1+moverate|subj_id), data, REML=FALSE)

ggploteffects(mdl=mdl,terms='moverate [all]',col=colall, linewidth=0.5)+
  my_theme+
  xlab('Movement rate (1/s)')+
  ylab('Estimated log(100+HGP) (%)')+
  annotate('text',label='***',x=1.6,y=5.2,size=6)

#### figure 4B ####
data=read.csv(paste0(dir, "data_triallevel.csv"))

mdl=lmer(log(100+hgp_m1) ~ moverate * group + perfhand + (1+moverate|subj_id), data, REML=FALSE)

# young
eff=ggeffect(mdl, terms=c('moverate [all]','group'), length.out=100)
eff=eff[eff$group=='young',]

ggplot(eff, aes(x=x, y=predicted)) +
  geom_line(color=colyoung,linewidth=0.5)+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high),fill=colyoung,alpha=0.2)+
  my_theme+
  coord_cartesian(xlim=c(1, 3), ylim=c(4.7,5.4))+
  scale_y_continuous(breaks=c(4.8,5.2))+
  scale_x_continuous(breaks=c(1,2,3))+
  xlab('Movement rate (1/s)')+
  ylab('Estimated log(100+HGP) (%)')

# control
eff=ggeffect(mdl, terms=c('moverate [all]','group'))
eff=eff[eff$group=='control',]

ggplot(eff, aes(x=x, y=predicted)) +
  geom_line(color=colcontrol,linewidth=0.5)+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high),fill=colcontrol,alpha=0.2)+
  my_theme+
  coord_cartesian(xlim=c(0.8, 2.1), ylim=c(4.6,5.3))+
  scale_y_continuous(breaks=c(4.8,5.2))+
  scale_x_continuous(breaks=c(1,1.5,2))+
  xlab('Movement rate (1/s)')+
  ylab('Estimated log(100+HGP) (%)')

# stroke
eff=ggeffect(mdl, terms=c('moverate [all]','group'))
eff=eff[eff$group=='stroke',]

ggplot(eff, aes(x=x, y=predicted)) +
  geom_line(color=colstroke,linewidth=0.5)+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high),fill=colstroke,alpha=0.2)+
  my_theme+
  coord_cartesian(xlim=c(0.48, 1.7), ylim=c(4.6,5))+
  scale_y_continuous(breaks=c(4.7,4.9))+
  scale_x_continuous(breaks=c(0.5,1.0,1.5))+
  xlab('Movement rate (1/s)')+
  ylab('Estimated log(100+HGP) (%)')

#### figure 4C ####
data=read.csv(paste0(dir, "data_subjectlevel.csv"))

mdl=lm(hgp_m1 ~ moverate + group + perfhand, data)

eff=ggeffect(mdl, terms = 'moverate [all]')
ggplot(data,aes(moverate, hgp_m1)) + 
  geom_point(aes(col=group, shape=group),size=2)+
  scale_color_manual(values = c(colcontrol, colstroke, colyoung)) +
  scale_shape_manual(values=c(17, 15, 16))+
  geom_line(aes(x,predicted), eff)+
  geom_ribbon(aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high),data=eff,fill=colall,alpha=0.2)+
  my_theme+
  xlab('Movement rate (1/s)')+
  ylab('Estimated HGP (%)')+
  annotate('text',label='**',x=1.55,y=110,size=6)

#### figure 5 ####
data=read.csv(paste0(dir, "data_subjectlevel.csv"))
regions <- c("peakfreq_m1", "peakfreq_pmc")
df=data[,c(regions, "group", 'subj_id')]
df$group <- factor(df$group,levels = c('young','control','stroke'),ordered = TRUE)

# Reshape the data into long format
df_long <- df %>%
  tidyr::pivot_longer(cols = starts_with(c("peak")), names_to = "region")

df_long$region <- factor(df_long$region,levels = regions,ordered = TRUE)

df_long <- df_long %>%
  mutate(region_group = interaction(group, region, sep = "_"))

p_values= df_long %>% 
  group_by(region) %>% 
  wilcox_test(value ~ group,p.adjust.method='fdr',comparisons=list(c('young','control'),c('stroke','control')))
p_values$p.adj=p.adjust(p_values$p, method='fdr')
p_values$p.adj.signif=c('','','','')
p_values=p_values[-c(2,4),]

y_values=rep(93,2)

# Create boxplot using ggplot with faceting
ggplot(df_long, aes(x = region_group, y = value, color = group)) +
  geom_boxplot(outlier.shape=NA, size=0.5) +
  geom_jitter(position = position_jitterdodge(jitter.width = 2),size=1, shape=16) +
  scale_color_manual(values = c(colyoung, colcontrol, colstroke)) +
  scale_x_manual(as.double(c(1,2.15,3.3,5,6.15,7.3)),labels=c('Young','Control','Stroke','Young','Control','Stroke'))+
  my_theme+
  theme(legend.position='none')+ 
  geom_signif(y_position = y_values, 
              xmin = c(1,5), 
              xmax = c(2.15,6.15), 
              annotations = p_values$p.adj.signif, 
              color='black', 
              size=0.3)+
  annotate('text',label='***',x=1.575,y=93.35,size=6)+
  annotate('text',label='***',x=5.575,y=93.35,size=6)+
  coord_cartesian(ylim=c(NaN,95))+
  xlab('Group')+
  ylab('Peak frequency (Hz)')

#### figure 6B ####
data=read.csv(paste0(dir, "data_timecourses.csv"))
timepoints=read.csv(paste0(dir, "timepoints_buttonpresses.csv"))
data$pow=data$hgp_m1

# mean average time points of button presses over subjects
second=tapply(timepoints$second, timepoints$group, mean)
third=tapply(timepoints$third, timepoints$group, mean)
fourth=tapply(timepoints$fourth, timepoints$group, mean)

# calculate mean and median absolute deviation of high gamma power over subjects for each time bin
stroke=tapply(data$pow[data$group=='stroke'], data$time[data$group=='stroke'], mean)
control=tapply(data$pow[data$group=='control'], data$time[data$group=='control'], mean)
young=tapply(data$pow[data$group=='young'], data$time[data$group=='young'], mean)
ci_stroke=tapply(data$pow[data$group=='stroke'], data$time[data$group=='stroke'], mad)
ci_control=tapply(data$pow[data$group=='control'], data$time[data$group=='control'], mad)
ci_young=tapply(data$pow[data$group=='young'], data$time[data$group=='young'], mad)
time=round(data$time[data$subj_id=='y01'],3)
df <- data.frame(stroke, control, young, time)

# time course of young participant group
ggplot(df, aes(time, young)) +
  geom_line(color = colyoung) +
  geom_ribbon(aes(time, ymin = young - ci_young, ymax = young + ci_young),
              fill = colyoung, alpha = 0.2)+
  geom_vline(xintercept = 0, color = "black", linetype = 'dashed')+
  geom_vline(xintercept = second[3], color = "black", linetype = 'dashed')+
  geom_vline(xintercept = third[3], color = "black", linetype = 'dashed')+
  geom_vline(xintercept = fourth[3], color = "black", linetype = 'dashed')+
  xlim(-0.5, 1.5)+
  my_theme+
  scale_y_continuous(breaks=c(0,40,80))+
  xlab('Time (s)')+
  ylab('High gamma power (%)')

# time course of control participant group
ggplot(df, aes(time, control)) +
  geom_line(color = colcontrol) +
  geom_ribbon(aes(time, ymin = control - ci_control, ymax = control + ci_control),
              fill = colcontrol, alpha = 0.2)+
  geom_vline(xintercept = 0, color = "black", linetype = 'dashed')+
  geom_vline(xintercept = second[1], color = "black", linetype = 'dashed')+
  geom_vline(xintercept = third[1], color = "black", linetype = 'dashed')+
  geom_vline(xintercept = fourth[1], color = "black", linetype = 'dashed')+
  xlim(-0.5, 1.5)+
  my_theme+
  scale_y_continuous(breaks=c(0,30,60))+
  xlab('Time (s)')+
  ylab('High gamma power (%)')

# time course of stroke group
ggplot(df, aes(time, stroke)) +
  geom_line(color = colstroke) +
  geom_ribbon(aes(time, ymin = stroke - ci_stroke, ymax = stroke + ci_stroke),
              fill = colstroke, alpha = 0.2)+
  geom_vline(xintercept = 0, color = "black", linetype = 'dashed')+
  geom_vline(xintercept = second[2], color = "black", linetype = 'dashed')+
  geom_vline(xintercept = third[2], color = "black", linetype = 'dashed')+
  geom_vline(xintercept = fourth[2], color = "black", linetype = 'dashed')+
  xlim(-0.5, 1.5)+
  my_theme+
  xlab('Time (s)')+
  ylab('High gamma power (%)')

#### figure 6C ####
data=read.csv(paste0(dir, "data_subjectlevel.csv"))
data$npeaks_m1=as.factor(data$npeaks_m1)

mdl=clm(npeaks_m1 ~ group + moverate + perfhand, data=data)

eff=ggeffect(mdl, terms = 'moverate [all]')
eff$response.level=factor(eff$response.level)

ggplot(eff, aes(x=x,y=predicted, group=response.level,color=response.level)) +
  geom_line(aes(linetype=response.level))+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=response.level),color=NA,alpha=0.2)+
  scale_color_manual(values=c('#FF055F','#00008B','#FFC107','#000000'))+
  scale_fill_manual(values=c('#FF055F','#00008B','#FFC107','#000000'))+
  scale_linetype_manual(values=c(3,1,2,4))+
  annotate('text',label='*',x=1.5,y=0.9,size=6)+
  my_theme+
  scale_y_continuous(breaks=c(0,0.5,1))+
  scale_x_continuous(breaks=c(1,2))+
  xlab('Movement rate (1/s)')+
  ylab('Probability')