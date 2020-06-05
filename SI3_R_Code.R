library(metafor) # models the meta-regression
library(MuMIn) # for the model fitting
library(multcomp) # for the multiple comparisons
eval(metafor:::.MuMIn) # this is necessary for model fitting

# function to calculate I^2 (based on http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate)
I2 <- function(mod, df){
  W <- diag(1/df$vv)
  X <- model.matrix(mod)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  return(c((100 * sum(mod$sigma2) / (sum(mod$sigma2) + (mod$k-mod$p)/sum(diag(P)))),
           (100 * mod$sigma2 / (sum(mod$sigma2) + (mod$k-mod$p)/sum(diag(P))))))
}

d <- read.csv('SI4_data.csv') #load data
## Calculating additive treatment, log risk ratio and variance.

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
    addD$Experimental.proportion <- pestD$Experimental.proportion + paraD$Experimental.proportion # additive treatment is sum of stressors/ control- divide by control later
    addD$Experimental.sample.size <- (pestD$Experimental.sample.size + paraD$Experimental.sample.size)/2 #for trials where pests and paras had different sample sizes I average them
    addD$Experimental.died <- addD$Experimental.proportion*addD$Experimental.sample.size
    addD$Experimental.survived <- addD$Experimental.sample.size - addD$Experimental.died
    addD$Treatment <- 'additive'
    allt <- rbind(allt, pestD, paraD, addD, combD) # add them all together
  }}
head(allt) #check results
summary(allt)
## Now calculating lRR and variance

#calculaing the lRR:
allt$lRR <- log(allt$Experimental.proportion/allt$Control.proportion)
#and variance. Because a few additive proportions are above 1, a few variances are slightly below zero, so I'm taking the absolute value
# these numbers are very small, so it shouldn't make a big difference:
allt$vv <- abs((1/allt$Experimental.died) - (1/allt$Experimental.sample.size) + (1/allt$Control.died) - (1/allt$Control.sample.size))
d <- allt # save to d

#now make separate dfs by additive, and all measured treatments
combo <- subset(d, Treatment == "combined")
ad <- subset(d, Treatment == "additive")
pest <- subset(d, Treatment == "pesticide")
para <- subset(d, Treatment == "parasite")

##combine subsets to compare by treatment
# for the main model, to compare the combined treatment to the predicted additive:
add <- droplevels(rbind(ad, combo))
# for the model comparing the measured treatments to determine the effects:
com <- droplevels(rbind(combo, para, pest))

## model fitting:
# main model: additive vs combined:

# construct a model with all considered treatments:
# (AIC scores aren't comparable for mixed-effects models, so model fitting is carried out with fixed-effects model)
add_dredge <- rma.mv(lRR, vv, mods = ~ Treatment*Accommodation + 
                        Treatment*Day + 
                        Treatment*Parasite.kind + 
                        Treatment*Pesticide.family +
                        Treatment*Life.stage+
                        Accommodation*Day,
                      data = add, method = "FE")
# dredge fits models of all variations of the above moderators and interactions orders them by increasing AICc scores.
# this takes a minute to run
(res_add <- dredge(add_dredge))
# I choose the top model because it has the lowest AICc score, and also because it's the most parsimonious (most simple) of all 
# models within two AICc's of each other. If models have similar scores, the law of parsimony suggests you go with the most simple model.

# I want to be able to compare combined to additive, so I'm making 'additive' the intercept. (It should be by default- it's supposed
# to be alphabetical, but it isn't, and I'm not sure why)
add$Treatment <- relevel(add$Treatment, ref = "additive")
# here's the top model:
(add_mod <- rma.mv(lRR, vv, mods = ~ Treatment*Accommodation + 
                     Parasite.kind + 
                     Treatment*Pesticide.family +
                     Accommodation*Day,
# random effects are added to account for within study and within trial (separate experiments in the study) heterogeneity.
# I also tried the random model structure Trial|Study but the results were similar, the model took longer to fit, and it was more
# difficult to extract I^2 scores.
                     random = ~1|Study/Trial, 
# 'test = 't'' adjusts confidence intervals and p-values "which slightly mimics the Knapp and Hartung (2003) method" according to the man pages
# Viechtbauer (the author of the package) has reported that the Knapp and Hartung method is one of the best ways to do this adjustment, on par with
# bootstrapping- and bootstrapping is a little complicated with mixed-effects models. 
                     data = add, test = "t"))

# Save the results in a data frame: I'm including the raw effects, transforming the effects, and also transforming the absolute value of the effects,
# since I think reporting all the results as greater than 1 makes more sense, given the way the anti-logged data gets transformed
# (logging makes the data linear enough to analyze, anti-logging makes it log scale again, so it gets a little weird to interpret)
add_data <- as.data.frame(add_mod$b)
names(add_data) <- "raw_effect"
add_data$trans_effect <- exp(add_mod$b)
add_data$abs_trans_effect <- exp(abs(add_mod$b))
add_data$cilb <- add_mod$ci.lb
add_data$ciub <- add_mod$ci.ub
add_data$t <- add_mod$zval
add_data$abs_trans_cilb <- exp(abs(add_mod$ci.lb))
add_data$abs_trans_ciub <- exp(abs(add_mod$ci.ub))
add_data$pval <- add_mod$pval
# I'm going to graph it later, so I don't want the mod names as row names, but they're really long, so I don't want them repeated
add_data$type <- rownames(add_data)
rownames(add_data) <- 1:nrow(add_data)
add_data <- add_data[,c(10,1:9)]
add_data
# write to a file (row.names = F because otherwise when I load it again, I'll have an annoying extra column)
write.csv(add_data, "Main_model_res.csv", row.names = F)

# some tests of the model:
# We can test for asymmetry with a rank correlation test. There are other tests
# available, but they don't work on multivariate models
ranktest(add_mod)

# Rank Correlation Test for Funnel Plot Asymmetry
# 
# Kendall's tau = -0.0296, p = 0.7051

#where the null hypothesis is that there is no bias- so we're failing to reject the null hypothesis of no bias

# here's the function I included above testing for the I^2, which is a measure of how much the model explains the heterogeneity:
I2(add_mod, add)
#98.369880 92.772347  5.597533
# so there is 98.4% heterogeneity, which is pretty huge. The study random effect accounts for 92.8% of it and the trial random
# effect accounts for 5.6%. These ecology mixed-effects models tend to have
# really high I^2s because we're including disparate studies- but that's kind of the whole point.

# this test tells us if any of the studies have a major influence on the data:
(cd_add <- cooks.distance(add_mod, progbar=T, cluster = add$Study))
#sort by weight
sort(cd_add)
# highest influence: Gregoric 15.9, Doublet 13.4 not too bad. For some reason, Retschnig 2015 appears here as NA, which I have
# not been able to figure out.

## Model comparing pesticide, parasite and combined treatments:

# now we're going to do the same thing for the second model: comparing pesticides, parasites, and combined treatments:
# notice that pesticides and parasites aren't included here: because the parasite treatments don't have pesticide applications, vice versa,
# including those would be including false information
com_dredge <- rma.mv(lRR, vv, mods = ~ Treatment*Accommodation + 
                         Treatment*Day + 
                         Treatment*Life.stage+
                         Accommodation*Day,
                        data = com, method = "FE")
(res_com <- dredge(com_dredge))
# here my best model is less parsimonious than the one after it, so I'm going to go with the second best fit, since the AICc is only 0.82 higher and
# it's more simple. The results are pretty similar, so it's not really that important.

# because I want to compare between the effects of the three treatments, I'm including the -1; this makes this an intercept-free model and gives us estimates
# of the total mean effects of the three treatments, rather than the difference between the reference treatment and the other treatments
(com_mod <- rma.mv(lRR, vv, mods = ~ Treatment*Accommodation -1 + 
                     Treatment*Day +
                     Accommodation*Day,
                      random = ~ 1|Study/Trial,
                      data = com, test = "t"))
com_mod
# again, I'm making a data frame with anti-logged and anti-logged absolute values:
com_data <- as.data.frame(com_mod$b)
names(com_data) <- "coef"
com_data$trans_effect <- exp(com_mod$b)
com_data$abs_trans_effect <- exp(abs(com_mod$b))
com_data$t <- com_mod$zval
com_data$cilb <- com_mod$ci.lb
com_data$ciub <- com_mod$ci.ub
com_data$abs_trans_cilb <- exp(abs(com_mod$ci.lb))
com_data$abs_trans_ciub <- exp(abs(com_mod$ci.ub))
com_data$pval <- com_mod$pval
com_data$type <- row.names(com_data)
row.names(com_data) <- 1:nrow(com_data)
com_data <- com_data[,c(10,1:9)]
com_data
write.csv(com_data, "Com_model_res.csv", row.names = F)
# this linear hypotheses test compares all three treatments and corrects the p-values because I'm running the test multiple times:
# all those zeros are telling it to ignore the other coefficients
(gl <- summary(glht(com_mod, linfct=cbind(contrMat(c("Combined"=1,"Parasite"=1,"Pesticide"=1), 
                                          type="Tukey"), 0, 0, 0, 0, 0, 0, 0)), test=adjusted("holm")))
# this gives me the confidence intervals for the above
(glc <- confint(glht(com_mod, linfct=cbind(contrMat(c("Combined"=1,"Parasite"=1,"Pesticide"=1), 
                                                   type="Tukey"), 0, 0, 0, 0, 0, 0, 0)), test=adjusted("holm")))
# now I'm combining those in a data frame:
glt <- as.data.frame(gl$test$coefficients)
glt$abs_exp_eff <- exp(abs(gl$test$coefficients))
glt$abs_exp <- exp(abs(glc$confint))
glt
write.csv(glt, "Pair_wise_test.csv")

I2(com_mod, com)


## Test to see if antagonism exists in the individual data sets in our studies
## of course, most of the researchers used different models- and had good reasons for doing so.
## and their data was better than ours, as extracting data from graphs isn't the most
## accurate. We're just trying to get a general idea of what antagonism might have been
## missed if the researchers weren't testing for it.

# load library for generalized linear models
library(logisticRR)

d <- read.csv('SI3_data.csv') # our data, again
head(d)
trials <- unique(d$Trial) # make a unique list of the trials
res <- data.frame(trials) # set up a new data frame
names(res) <- "Trial"
res$Interaction <- NA
res$pval <- NA
res$Type <- NA
res$Study <- NA
res$Accomodation <- NA
res$length <- NA
i <- 27
for (i in 1:length(trials)){
  nu <- subset(d, Trial == trials[i]) # chops the data to just the trial
  pest <- subset(nu, Treatment == "pesticide") # subset into treatments
  para <- subset(nu, Treatment == "parasite")
  combo <- subset(nu, Treatment == "combined")
  # make new data frames, with a row for each bee in study, 0 for survived, 1 for died
  dpe <- data.frame("died" = rep(c(0,1), c(pest$Experimental.survived, pest$Experimental.died)))
  dpe$Pesticide <- 1
  dpe$Parasite <- 0
  dc <- data.frame("died" = rep(c(0,1), c(pest$Control.survived, pest$Control.died)))
  dc$Pesticide <- 0
  dc$Parasite <- 0
  dpa <- data.frame("died" = rep(c(0,1), c(para$Experimental.survived, para$Experimental.died)))
  dpa$Pesticide <- 0
  dpa$Parasite <- 1
  dco <- data.frame("died" = rep(c(0,1), c(combo$Experimental.survived, combo$Experimental.died)))
  dco$Pesticide <- 1
  dco$Parasite <- 1
  dn <- rbind(dpe, dc, dpa, dco)
  gm <- glm(died ~ Pesticide*Parasite, family = binomial(link = "logit"), data = dn) # generalized linear model
  gms <- summary(gm)
  pest <- droplevels(pest)
  res$Study[i] <- as.character(pest$Study[1])
  res$Accomodation[i] <- as.character(pest$Accommodation[1])
  res$length[i] <- pest$Day[1]
  res$Interaction[i] <- gm$coefficients[[4]]
  res$pval[i] <- coef(gms)[4,4]
  if (res$pval[i] < 0.05){
    if (res$Interaction[i] > 0){
      res$Type[i] <- "Synergistic"
    }
    if (res$Interaction[i] < 0){
      res$Type[i] <- "Antagonistic"
    }}
  if (res$pval[i] >= 0.05){
    res$Type[i] <- "Additive"
  }    
}
# look at the number of synergistic, additive and antagonistic results
table(res$Type)

#Additive Antagonistic  Synergistic 
#26            3           10 
# lots more synergistic than antagonistic, but still some antagonistic results
summary(res$Interaction)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -2.6716 -0.3467  0.3758  1.5263  0.9681 35.4077 
hist(res$Interaction)

