# Bird et al 2020 R code
# R version 3.6.1 (2019-07-05) -- "Action of the Toes"
# R studio 1.2.1335

#Load the required libraries:
library(metafor) # models the meta-regression
library(MuMIn) # for the model fitting
library(multcomp) # for the multiple comparisons
eval(metafor:::.MuMIn) # this is necessary for model fitting

# function to calculate I^2 (based on http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate)
I2 <- function(mod, df){
  W <- diag(1/df$vi)
  X <- model.matrix(mod)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  return(c((100 * sum(mod$sigma2) / (sum(mod$sigma2) + (mod$k-mod$p)/sum(diag(P)))),
           (100 * mod$sigma2 / (sum(mod$sigma2) + (mod$k-mod$p)/sum(diag(P))))))
}
#read in data (available at https://github.com/birdoptera/Honey_bee_antagonism)
d <- read.csv('bee_data_10.csv') 
#check out data
table(d$Study)
summary(d)

## Calculating additive and multiplicative treatments, log risk ratio and variance.

# calculate # died and survived
d$Experimental.died <- d$Experimental.proportion*d$Experimental.sample.size
d$Experimental.survived <- d$Experimental.sample.size - d$Experimental.died
d$Control.died <- d$Control.proportion*d$Control.sample.size
d$Control.survived <- d$Control.sample.size - d$Control.died

#this for loop calculates predicted additive treatments:
allt <- subset(d, Treatment == 0) # make a new data frame with same row names
trials <- unique(d$Trial) # make a list of unique trials
t <- trials[1]
for (t in trials){
  lug <- subset(d, Trial==t) #divides data set into groups by trial
  if (nrow(lug) == 3){ # removes trial groups where fewer than three, just in case any exist
    pestD <- subset(lug, Treatment=='pesticide') # identify rows for each treatment
    paraD <- subset(lug, Treatment=='parasite')
    combD <- subset(lug, Treatment=='combined')
    addD <- pestD # replicate pesticide row to keep other moderators (e.g. accommodation, life.stage) the same
    addD$Experimental.proportion <- (pestD$Experimental.proportion + paraD$Experimental.proportion)
    # additive treatment is sum of stressors/ control- divide by control later
    addD$Experimental.sample.size <- (pestD$Experimental.sample.size + paraD$Experimental.sample.size)/2
    #for trials where pests and paras had different sample sizes, we're taking an average
    addD$Experimental.died <- addD$Experimental.proportion*addD$Experimental.sample.size
    addD$Experimental.survived <- addD$Experimental.sample.size - addD$Experimental.died
    addD$Treatment <- 'additive'
    multiD <- addD
    # multiplicative treatment- from Aufauvre 2012, but found in many other places in the literature
    # I am using the term 'multiplicative' because it's used that way in Darling and Cote 2008
    # but the term seems to be used differently in the epidemiological literature
    # A + B - A*B and A + B/(1-A) (Aufauvre 2012) uses A + B/((1-A)/100) because they're using proportions*100 for some reason
    multiD$Experimental.proportion <- (pestD$Experimental.proportion + paraD$Experimental.proportion*(1-pestD$Experimental.proportion))
    multiD$Experimental.died <- multiD$Experimental.proportion*multiD$Experimental.sample.size
    multiD$Experimental.survived <- multiD$Experimental.sample.size - multiD$Experimental.died
    multiD$Treatment <- 'multiplicative'
    allt <- rbind(allt, pestD, paraD, addD, multiD, combD) # combine into one data frame
  }}
head(allt) #check results
summary(allt)

## Now calculating lRR and variance
#calculaing the lRR and variance:
d <- escalc("RR", ai = Experimental.died, bi = Experimental.survived, ci = Control.died,
            di = Control.survived, data = allt)

#now make separate dfs by additive, multiplicative, and all measured treatments
combo <- subset(d, Treatment == "combined")
ad <- subset(d, Treatment == "additive")
mult <- subset(d, Treatment == "multiplicative")
pest <- subset(d, Treatment == "pesticide")
para <- subset(d, Treatment == "parasite")

##combine subsets to compare by treatment
# for the main models, to compare the combined treatment to the predicted additive
# and multiplicative:
add <- droplevels(rbind(ad, combo))
multi <- droplevels(rbind(mult, combo))
# for the model comparing the measured treatments to determine the effects:
com <- droplevels(rbind(combo, para, pest))

## model fitting:
# main model: additive vs combined:

# I've categorized the pesticides according to 'neonicotinoid' and 'non-neonicotinoid', but one study
# mixed both together. I don't want this 'mixed' category to be the intercept because I don't care about it,
# so I'm making 'neonicotinoid' the reference and I'm going to ignore the 'mixed' result
# I'm also releveling 'virus'
add$Pesticide.group <- relevel(add$Pesticide.group, ref = "Neonicotinoid")
add$Parasite.kind <- relevel(add$Parasite.kind, ref = "Virus")
# here's the top model:
(add_mod <- rma.mv(yi, vi, mods = ~  Parasite.kind + Treatment +Accommodation + Day + Parasite.kind + Pesticide.group + Life.stage,
# random effects are added to account for within study and within trial (separate experiments in the study) heterogeneity.
# I also tried the random model structure Trial|Study but the results were similar, the model took longer to fit, and it was more
# difficult to extract I^2 scores.
                     random = ~1|Study/Trial, 
# 'test = 't'' adjusts confidence intervals and p-values "which slightly mimics the Knapp and Hartung (2003) method" according to the man pages for metafor.
# Viechtbauer (the author of metafor) has reported that the Knapp and Hartung method is one of the best ways to do this adjustment, on par with
# bootstrapping- and bootstrapping is a little complicated with mixed-effects models. 
                     data = add, test = "t"))
# test for colliniarity
vif(add_mod)
# doesn't seem to be any colliniarity present

# pairwise comparison. The 1s and 0s correspond to each predictor in the model, with the first being the intercept.
# the first three just give the same results as the model summary- the last three compare those three levels
# against each other. P-values and CIs are adjusted for multiple testing.
a_gl <- glht(add_mod, linfct=rbind(c(0,1,0,0,0,0,0,0,0,0), c(0,0,1,0,0,0,0,0,0,0),
                                   c(0,0,0,1,0,0,0,0,0,0), c(0,-1,1,0,0,0,0,0,0,0),
                                   c(0,-1,0,1,0,0,0,0,0,0), c(0,0,-1,1,0,0,0,0,0,0)
                                    ), test=adjusted("holm"))
#summary
(a_gls <- summary(glht(a_gl)))
# this gives us the CI
(a_glc <- confint(glht(a_gl))$confint)

# here's the function I included above testing for the I^2, which is a measure of how much the model explains the heterogeneity:
I2(add_mod, add)
#95.78458 79.56242 16.22216
# so there is 95.8% heterogeneity, which is pretty huge. The study random effect accounts for 79.6% of it and the trial random
# effect accounts for 16.2%. These ecology mixed-effects models tend to have
# really high I^2s because we're including disparate studies- but that's kind of the whole point.

# this test tells us if any of the studies have a major influence on the data:
(cd_add <- cooks.distance(add_mod, progbar=T, cluster = add$Study))
#sort by weight
sort(cd_add)
# these Cook's distances aren't very big, so I'm not worried about influence
# Kohler 2012   Doublet 2015 Aufauvre 2014 
# 2.14042305     2.43721659    3.14048181

## same but with multiplicative model: 

multi$Pesticide.group <- relevel(multi$Pesticide.group, ref = "Neonicotinoid")
multi$Parasite.kind <- relevel(multi$Parasite.kind, ref = "Virus")

(multi_mod <- rma.mv(yi, vi, mods = ~  Parasite.kind + Treatment +Accommodation + Day + Pesticide.group + Life.stage,
                   random = ~1|Study/Trial,
                   data = multi, test = "t"))
# test for colliniarity
vif(multi_mod)


# pairwise testing
m_gl <- glht(multi_mod, linfct=rbind(c(0,1,0,0,0,0,0,0,0,0), c(0,0,1,0,0,0,0,0,0,0),
                                   c(0,0,0,1,0,0,0,0,0,0), c(0,-1,1,0,0,0,0,0,0,0),
                                   c(0,-1,0,1,0,0,0,0,0,0), c(0,0,-1,1,0,0,0,0,0,0)), test=adjusted("holm"))

# I^2
I2(multi_mod, multi)
#95.75837 83.18497 12.57340

# cook's distance
(cd_multi <- cooks.distance(multi_mod, progbar=T, cluster = multi$Study))

sort(cd_multi)
# Gregoric 2016   Doublet 2015  Aufauvre 2014  
# 2.26244767     3.46195973     3.49052831 

## Model comparing pesticide, parasite and combined treatments:

# now we're going to do the same thing for the second model: comparing pesticides, parasites, and combined treatments:
# notice that pesticides and parasites aren't included here: because the parasite treatments don't have pesticide applications, vice versa,
# including those would be including false information

#running a no-intercept model to see how the treatments compare to control- we'll be doing pairwise testing to see
# how they compare to each other
(com_mod <- rma.mv(yi, vi, mods = ~ Treatment-1 + Accommodation
                     + Day + Life.stage,
                      random = ~ 1|Study/Trial,
                      data = com, test = "t"))
com_mod

#pairwise testing
c_gl <- glht(com_mod, linfct=rbind(c(-1,1,0,0,0,0), c(-1,0,1,0,0,0),
                                     c(0,-1,1,0,0,0)), test=adjusted("holm"))

(c_gls <- summary(glht(c_gl)))

(c_glc <- confint(glht(c_gl))$confint)

#I^2
I2(com_mod, com)
#95.03490 79.59377 15.44114

#Eggar's test, modified for multivariate models, as suggested by Viechtbaur,
# https://stats.stackexchange.com/questions/155693/metafor-package-bias-and-sensitivity-diagnostics
# with sample size as moderator - the usual test sees if there is a correlation
# between effect size and variance- Viechtbaur suggests that it's suspect when the variance is calculated
# from effect size, but for risk ratios, variance is calculated from both effect size and sample size,
# so I chose to test for correlation between effect size and sample size

# take the mean of control and experimental sample size
com$mss <- (com$Control.sample.size + com$Experimental.sample.size)/2

(eggar <- rma.mv(yi, vi, mods = ~ Treatment + Accommodation
                   + Day + Life.stage + mss,
                   random = ~ 1|Study/Trial,
                   data = com, test = "t"))
#mss                   -0.0015  0.0010   -1.4946  0.1368  -0.0035   0.0005 
# no correlation between sample size and lRR suggests no bias

(cd_com <- cooks.distance(com_mod, progbar=T, cluster = com$Study))
sort(cd_com)
#whoa- pretty big influence from Retschnig 2015 there. Enough reason to do a leave-one-out analysis
# there are no built in functions for rma.vi in metafor, so we have to hand-code it. ugh.

com_cd <- data.frame("Rmv_Study" = as.character(), "combined" = as.numeric(), 
                     "com_p" = as.numeric(), "pesticide" = as.numeric(), 
                     "pes_p" = as.numeric(), "parasite" = as.numeric(), 
                     "para_p" = as.numeric())
studies <- unique(com$Study)
i <- studies[1]
for (i in studies){
  com_t <- subset(com, Study != i)
  com_t_mod <- rma.mv(yi, vi, mods = ~  Treatment-1 +Accommodation + Day + Life.stage,
                      random = ~1|Study/Trial, 
                      data = com_t, test = "t")
  com_t_mod$pval[2]
  com_ct <- data.frame("Rmv_Study" = i, "combined" = com_t_mod$b[1], 
                       "com_p" = com_t_mod$pval[1], "pesticide" = com_t_mod$b[3], 
                       "pes_p" = com_t_mod$pval[3], "parasite" = com_t_mod$b[2], 
                       "para_p" = com_t_mod$pval[2])
  com_cd <- rbind(com_cd, com_ct)
}

com_cd

# plot leave one out analysis
library(tidyverse)
com_cd <- gather(com_cd, "Treatment", "Coef", c(2,4,6))
ggplot(com_cd, aes(x = Rmv_Study, y = Coef, color = Treatment, group = Treatment)) +
  geom_point() +
  geom_line() +
  labs(x = "Study removed", y = "log Risk Ratio") +
  scale_color_manual(values = c("#487d68","#E67E22", "#D35400")) +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90),
        panel.background = element_blank(),
        legend.key = element_blank())

# Aufauvre 2014 seems to reduce the effects of individual treatments and Gregorc 2016 increases them pretty strongly
# but I don't think its too serious

