#### Hydrologic Data Processessing ####
library(xts)
library(prob)
library(TTR)
library(bigmemory)

#### Set working directory ####
setwd("/Users/Morrison/Documents/Research Projects/Riparian Bayesian Network/climate_bayesian")

#### Source common functions ####
source('~/Documents/Research Projects/Riparian Bayesian Network/climate_bayesian/r_code/1-functions.R')

hist_base <- climate.hydrology("data/HAD_hist_base.csv")

#### Add stage values (feet and cm) to data frame based on a regression equation ####
# The regression equation was developed using historical gage data, and changes depending on the site location

# SITE 1
# head_conv <- function(Q) 0.0011*Q

# SITE 2
# head_conv <- function(Q) 0.0011*Q

# SITE 3
# head_conv <- function(Q) 0.056*(Q^0.5612)

# SITE 4
# head_conv <- function(Q) 0.0014*(Q^0.99)

# SITE 5
head_conv <- function(Q) 0.1123*(Q^0.4963)

hist_base_stage <- stage(hist_base, head_conv)
hist_base <- data.frame(hist_base, hist_base_stage)

#### Calculate recession rates ####
# The 1-, 2-, ... , 7-day, 14-day recession rates were computed
recess_1d <- recession.rate.forward(hist_base, 1)
recess_2d <- recession.rate.forward(hist_base, 2)
recess_3d <- recession.rate.forward(hist_base, 3)
recess_4d <- recession.rate.forward(hist_base, 4)
recess_5d <- recession.rate.forward(hist_base, 5)
recess_6d <- recession.rate.forward(hist_base, 6)
recess_7d <- recession.rate.forward(hist_base, 7)
recess_14d <- recession.rate.forward(hist_base, 14)
recess_30d <- recession.rate.forward(hist_base, 30)
recess_45d <- recession.rate.forward(hist_base, 45)

hist_base <- cbind(hist_base, recess_1d, recess_2d, recess_3d, recess_4d, recess_5d, recess_6d, recess_7d, recess_14d, recess_30d, recess_45d)
colnames(hist_base)[5:14] <- c("recess_1d", "recess_2d", "recess_3d", "recess_4d", "recess_5d", "recess_6d", "recess_7d", "recess_14d", "recess_30d", "recess_45d")

#### Process Q-model bins ####
# Model Q-bins were selected so that inundation areas increased by less than 15% between the discharge associated with each bin
# q_bin <- c(719, 1209, 2100, 4000, 6296, 11300)
# q_bin <- c(1000, 1200, 1400, 1600, 1800, 2000, 2200)
q_bin <- c(1000, 1500, 2000, 2500, 3000, 3500, 4000)

# Calculate probability of flooding in each bin for the entire year
q_bin_prob <- q.prob(hist_base, q_bin)

# Calculate probabilities in each bin according to discrete timing states
# The discrete timing states are April-May, June-July, August-September
hist_base_monthly <- subset.month(hist_base)

apr_may <- list.to.df(hist_base_monthly[[1]])
jun_jul <- list.to.df(hist_base_monthly[[2]])
aug_sep <- list.to.df(hist_base_monthly[[3]])

q_prob_apr_may <- q.prob(apr_may, q_bin)
q_prob_jun_jul <- q.prob(jun_jul, q_bin)
q_prob_aug_sep <- q.prob(aug_sep, q_bin)

HB_q_prob_all <- data.frame(q_prob_apr_may, q_prob_jun_jul, q_prob_aug_sep)
# write.table(q_prob_all, "output/q_prob_all.txt", sep="\t")

# Test calculations
# R <- probspace(apr_may)
# Y <- prob(R, cfs>=719)
#
# C <- count(apr_may$cfs>=719)

#### Calculate probabilities for recession rates ####
# Test calculations
# S <- probspace(apr_may)
# P <- prob(S, recess_6d>=6)
# X <- S[["recess_1d"]]

# Set discrete states for recession rates (in cm/day)
recess_rates <- c(0, 1, 3, 6)

# Calculate recession rate probabilities for April-May
recess_1d_apr_may <- r.prob(apr_may, "recess_1d", recess_rates)
recess_2d_apr_may <- r.prob(apr_may, "recess_2d", recess_rates)
recess_3d_apr_may <- r.prob(apr_may, "recess_3d", recess_rates)
recess_4d_apr_may <- r.prob(apr_may, "recess_4d", recess_rates)
recess_5d_apr_may <- r.prob(apr_may, "recess_5d", recess_rates)
recess_6d_apr_may <- r.prob(apr_may, "recess_6d", recess_rates)
recess_7d_apr_may <- r.prob(apr_may, "recess_7d", recess_rates)
recess_14d_apr_may <- r.prob(apr_may, "recess_14d", recess_rates)
recess_30d_apr_may <- r.prob(apr_may, "recess_30d", recess_rates)
recess_45d_apr_may <- r.prob(apr_may, "recess_45d", recess_rates)

recess_all_apr_may <- cbind(recess_1d_apr_may, recess_2d_apr_may, recess_3d_apr_may, recess_4d_apr_may, recess_5d_apr_may, recess_6d_apr_may, recess_7d_apr_may, recess_14d_apr_may, recess_30d_apr_may, recess_45d_apr_may)
colnames(recess_all_apr_may) <- c("recess_1d", "recess_2d", "recess_3d", "recess_4d", "recess_5d", "recess_6d", "recess_7d", "recess_14d", "recess_30d", "recess_45d")
rownames(recess_all_apr_may) <- c("<0", "0-1", "1-3", "3-6", ">6")
recess_all_apr_may <- t(recess_all_apr_may)

# Calculate recession rate probabilities for June-July
recess_1d_jun_jul <- r.prob(jun_jul, "recess_1d", recess_rates)
recess_2d_jun_jul <- r.prob(jun_jul, "recess_2d", recess_rates)
recess_3d_jun_jul <- r.prob(jun_jul, "recess_3d", recess_rates)
recess_4d_jun_jul <- r.prob(jun_jul, "recess_4d", recess_rates)
recess_5d_jun_jul <- r.prob(jun_jul, "recess_5d", recess_rates)
recess_6d_jun_jul <- r.prob(jun_jul, "recess_6d", recess_rates)
recess_7d_jun_jul <- r.prob(jun_jul, "recess_7d", recess_rates)
recess_14d_jun_jul <- r.prob(jun_jul, "recess_14d", recess_rates)
recess_30d_jun_jul <- r.prob(jun_jul, "recess_30d", recess_rates)
recess_45d_jun_jul <- r.prob(jun_jul, "recess_45d", recess_rates)

recess_all_jun_jul <- cbind(recess_1d_jun_jul, recess_2d_jun_jul, recess_3d_jun_jul, recess_4d_jun_jul, recess_5d_jun_jul, recess_6d_jun_jul, recess_7d_jun_jul, recess_14d_jun_jul, recess_30d_jun_jul, recess_45d_jun_jul)
colnames(recess_all_jun_jul) <- c("recess_1d", "recess_2d", "recess_3d", "recess_4d", "recess_5d", "recess_6d", "recess_7d", "recess_14d", "recess_30d", "recess_45d")
rownames(recess_all_jun_jul) <- c("<0", "0-1", "1-3", "3-6", ">6")
recess_all_jun_jul <- t(recess_all_jun_jul)

# Calculate recession rate probabilities for August-September
recess_1d_aug_sep <- r.prob(aug_sep, "recess_1d", recess_rates)
recess_2d_aug_sep <- r.prob(aug_sep, "recess_2d", recess_rates)
recess_3d_aug_sep <- r.prob(aug_sep, "recess_3d", recess_rates)
recess_4d_aug_sep <- r.prob(aug_sep, "recess_4d", recess_rates)
recess_5d_aug_sep <- r.prob(aug_sep, "recess_5d", recess_rates)
recess_6d_aug_sep <- r.prob(aug_sep, "recess_6d", recess_rates)
recess_7d_aug_sep <- r.prob(aug_sep, "recess_7d", recess_rates)
recess_14d_aug_sep <- r.prob(aug_sep, "recess_14d", recess_rates)
recess_30d_aug_sep <- r.prob(aug_sep, "recess_30d", recess_rates)
recess_45d_aug_sep <- r.prob(aug_sep, "recess_45d", recess_rates)

recess_all_aug_sep <- cbind(recess_1d_aug_sep, recess_2d_aug_sep, recess_3d_aug_sep, recess_4d_aug_sep, recess_5d_aug_sep, recess_6d_aug_sep, recess_7d_aug_sep, recess_14d_aug_sep, recess_30d_aug_sep, recess_45d_aug_sep)
colnames(recess_all_aug_sep) <- c("recess_1d", "recess_2d", "recess_3d", "recess_4d", "recess_5d", "recess_6d", "recess_7d", "recess_14d", "recess_30d", "recess_45d")
rownames(recess_all_aug_sep) <- c("<0", "0-1", "1-3", "3-6", ">6")
recess_all_aug_sep <- t(recess_all_aug_sep)

HB_recess_prob_all <- cbind(recess_all_apr_may, recess_all_jun_jul, recess_all_aug_sep)

#### Populate network states based on scenarios ####
# Discrete states (1, 2, 3) based on the timing
HB_timing_state <- timing.states(hist_base)

# Discrete states (Y, N) based on inundation
HB_q1_state <- q.states(hist_base, q_bin, 1)
HB_q2_state <- q.states(hist_base, q_bin, 2)
HB_q3_state <- q.states(hist_base, q_bin, 3)
HB_q4_state <- q.states(hist_base, q_bin, 4)
HB_q5_state <- q.states(hist_base, q_bin, 5)
HB_q6_state <- q.states(hist_base, q_bin, 6)
HB_q7_state <- q.states(hist_base, q_bin, 7)

# Discrete states based on recession rates
HB_recess_state <- recess.states(hist_base, recess_rates, "recess_14d")