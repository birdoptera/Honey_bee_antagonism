# Bird et al 2020 R code- for individual non-parametric bootstraps
# R version 3.6.1 (2019-07-05) -- "Action of the Toes"
# R studio 1.2.1335

d <- read.csv('bee_data_10.csv') # our data, again
head(d)

trials <- sort(unique(d$Trial)) # make a unique list of the trials
res <- data.frame(trials) # set up a new data frame
names(res) <- "Trial"
res$Study <- NA # study name
res$add_res <- NA #interaction result for additive model
res$add_i <- NA # coef dif for additive model
res$add_p <- NA # p-value for additive model
res$multi_res <- NA #interaction result for multiplicative model
res$multi_i <- NA # coef dif for multiplicative model
res$multi_p <- NA # p-value for multiplicaite model



# function for calculating log risk ratios then subtracting the combined
# from the predictor
ustat <- function(x,y){ 
  pred <- mean(x)
  combod <- mean(y)
  pred - combod
}

set.seed(36849) # set seed for reproducibility
B <- 10000 # number of bootstraps

i <- trials[1]
for (i in 1:length(trials)){
  print(i/(length(trials))) #shows loop progress
  nu <- subset(d, Trial == trials[i]) # chops the data to just one trial
  pest <- subset(nu, Treatment == "pesticide") # subset into treatments
  para <- subset(nu, Treatment == "parasite")
  combo <- subset(nu, Treatment == "combined")
  res$Study[i] <- as.character(pest$Study[1])
  # first we have to try to recreate the data frames of the studies:
  #cnn <- combo$Control.sample.size
  #cnp <- combo$Control.proportion
  if (cnp == 0){cnp <- 0.001} #risk ratios can't be calculated with 0s, so add a tiny bit to them
  pn <- (pest$Experimental.sample.size + pest$Experimental.sample.size)/2 #mean sample size for predictor
  cmn <- combo$Experimental.sample.size
  add <- pest$Experimental.proportion + para$Experimental.proportion #calculate additive
  multi <- (add - pest$Experimental.proportion*para$Experimental.proportion)
  add <- add
  cnd <- combo$Experimental.proportion
  aad <- round(add*pn,0)
  mud <- round(multi*pn,0)
  cmd <- round(cnd*cmn,0)
  addd <- rep(c(1,0), c(aad, (pn-aad))) # make a list of dead (1) and living (0) bees for each variable
  multd <- rep(c(1,0), c(mud, (pn-mud)))
  comd <- rep(c(1,0), c(cmd, (cmn-cmd)))
  dataA <- ustat(addd, comd) #calculate observed statistic
  NullD <- c(addd, comd) #create null set
  P <- matrix(sample(NullD, pn*B, replace=T), B, pn)  #bootstrap null for predictor
  smP <- split(P, row(P))
  C <- matrix(sample(NullD, cmn*B, replace=T), B, cmn) #bootstrap null for combined
  smC <- split(C, row(C))
  NullU <- mapply(ustat, smP, smC) # calculate statistic 10000 times with null data
  da <- mean(NullU > dataA) # calculate the proportion of null bootstraps that are greater than the observed
  res$add_i[i] <- dataA # save the difference in risks of mortality between the predicted and the combined
  res$add_p[i] <- da
  if (da > .975){res$add_res[i] <- "Synergism" # if the null is greater than the observed 0.975% of the time
  # it's a significant synergism
  }else if (da < .025){res$add_res[i] <- "Antagonism" # if it's less than 0.025% of the time, it's antagonism
  } else {res$add_res[i] <- "No_interaction"}   # else, its not significant
# why 0.975 and 0.025 rather than 0.95 and 0.05? Because it's a two-tailed test:
# the alpha is set to 0.05- >0.95 + < 0.05 would equal 0.1%, not 0.05%.
  
# now the same for the multiplicaive model
  dataM <- ustat(multd, comd)
  NullD <- c(multd, comd)
  P <- matrix(sample(NullD, pn*B, replace=T), B, pn)
  smP <- split(P, row(P))
  C <- matrix(sample(NullD, cmn*B, replace=T), B, cmn) 
  smC <- split(C, row(C))
  NullG <- mapply(ustat, smP, smC)    # 3. 10000 Test statistics
  dm <- mean(NullG > dataM)
  res$multi_i[i] <- dataM
  res$multi_p[i] <- dm
  if (dm > .975){res$multi_res[i] <- "Synergism"
    }else if (dm < .025){res$multi_res[i] <- "Antagonism"
    } else {res$multi_res[i] <- "No_interaction"}
}

#check out results
table(res$add_res)
table(res$multi_res)


