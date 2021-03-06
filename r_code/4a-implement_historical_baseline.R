#### Read results from SRH-2D model ####

set.seed(140)

# Change input table depending on the site number
HB_alldata <- read.table("data/site5_1000_to_4000.txt", header=FALSE, skip=0)

colnames(HB_alldata) <- c("cell", "flood_cfs")

#### Separate grid data according to the Q-bins previously defined ####
HB_q1_cells <- subset(HB_alldata, flood_cfs>q_bin[1] & flood_cfs<q_bin[2])
HB_q2_cells <- subset(HB_alldata, flood_cfs>=q_bin[2] & flood_cfs<q_bin[3])
HB_q3_cells <- subset(HB_alldata, flood_cfs>=q_bin[3] & flood_cfs<q_bin[4])
HB_q4_cells <- subset(HB_alldata, flood_cfs>=q_bin[4] & flood_cfs<q_bin[5])
HB_q5_cells <- subset(HB_alldata, flood_cfs>=q_bin[5] & flood_cfs<q_bin[6])
HB_q6_cells <- subset(HB_alldata, flood_cfs>=q_bin[6] & flood_cfs<q_bin[7])
HB_q7_cells <- subset(HB_alldata, flood_cfs>=q_bin[7])

#### Combine hydrology time series, timing states, and recession rate states ####
HB_states <- cbind(hist_base, HB_timing_state, HB_recess_state)
#### Trim data frame based on timing and recession rates ####
HB_states <- subset(HB_states, TIMING != "NA" & RECESSION != 1)
HB_timing_subset <- HB_states[,15]
HB_recess_subset <- HB_states[,16]

save(HB_states, file="output/site5/HB_states.Rdata")

#### Allow dimension naming using the "bigmemory" package ####
options(bigmemory.allow.dimnames=TRUE)

#### Create empty matrix to fill with data based on the Q-bin ####
HB_emptymatrix1 <- big.matrix(nrow=nrow(HB_states), ncol=(nrow(HB_q1_cells)+2), dimnames=list(NULL, c("timing", "recess", as.character(HB_q1_cells[,1]))))

HB_emptymatrix2 <- big.matrix(nrow=nrow(HB_states), ncol=(nrow(HB_q2_cells)+2), dimnames=list(NULL, c("timing", "recess", as.character(HB_q2_cells[,1]))))

HB_emptymatrix3 <- big.matrix(nrow=nrow(HB_states), ncol=(nrow(HB_q3_cells)+2), dimnames=list(NULL, c("timing", "recess", as.character(HB_q3_cells[,1]))))

HB_emptymatrix4 <- big.matrix(nrow=nrow(HB_states), ncol=(nrow(HB_q4_cells)+2), dimnames=list(NULL, c("timing", "recess", as.character(HB_q4_cells[,1]))))

HB_emptymatrix5 <- big.matrix(nrow=nrow(HB_states), ncol=(nrow(HB_q5_cells)+2), dimnames=list(NULL, c("timing", "recess", as.character(HB_q5_cells[,1]))))

HB_emptymatrix6 <- big.matrix(nrow=nrow(HB_states), ncol=(nrow(HB_q6_cells)+2), dimnames=list(NULL, c("timing", "recess", as.character(HB_q6_cells[,1]))))

HB_emptymatrix7 <- big.matrix(nrow=nrow(HB_states), ncol=(nrow(HB_q7_cells)+2), dimnames=list(NULL, c("timing", "recess", as.character(HB_q7_cells[,1]))))

#### Populate empty matrices
# The warnings that appear when the empty matrices are populated are due to the differing lengths of the Boolean comparisons.
HB_emptymatrix1[,1] <- as.vector(as.matrix(HB_timing_subset))
HB_emptymatrix1[,2] <- as.vector(as.matrix(HB_recess_subset))
for (i in 1:nrow(HB_q1_cells)) {
		HB_emptymatrix1[,i+2] <- HB_q1_cells[i,2]<=HB_states[,2] & HB_q1_cells[i,2]>HB_states[-1,2]
}

HB_emptymatrix2[,1] <- as.vector(as.matrix(HB_timing_subset))
HB_emptymatrix2[,2] <- as.vector(as.matrix(HB_recess_subset))
for (i in 1:nrow(HB_q2_cells)) {
  HB_emptymatrix2[,i+2] <- HB_q2_cells[i,2]<=HB_states[,2] & HB_q2_cells[i,2]>HB_states[-1,2]
}

HB_emptymatrix3[,1] <- as.vector(as.matrix(HB_timing_subset))
HB_emptymatrix3[,2] <- as.vector(as.matrix(HB_recess_subset))
for (i in 1:nrow(HB_q3_cells)) {
  HB_emptymatrix3[,i+2] <- HB_q3_cells[i,2]<=HB_states[,2] & HB_q3_cells[i,2]>HB_states[-1,2]
}

HB_emptymatrix4[,1] <- as.vector(as.matrix(HB_timing_subset))
HB_emptymatrix4[,2] <- as.vector(as.matrix(HB_recess_subset))
for (i in 1:nrow(HB_q4_cells)) {
  HB_emptymatrix4[,i+2] <- HB_q4_cells[i,2]<=HB_states[,2] & HB_q4_cells[i,2]>HB_states[-1,2]
}

HB_emptymatrix5[,1] <- as.vector(as.matrix(HB_timing_subset))
HB_emptymatrix5[,2] <- as.vector(as.matrix(HB_recess_subset))
for (i in 1:nrow(HB_q5_cells)) {
  HB_emptymatrix5[,i+2] <- HB_q5_cells[i,2]<=HB_states[,2] & HB_q5_cells[i,2]>HB_states[-1,2]
}

HB_emptymatrix6[,1] <- as.vector(as.matrix(HB_timing_subset))
HB_emptymatrix6[,2] <- as.vector(as.matrix(HB_recess_subset))
for (i in 1:nrow(HB_q6_cells)) {
  HB_emptymatrix6[,i+2] <- HB_q6_cells[i,2]<=HB_states[,2] & HB_q6_cells[i,2]>HB_states[-1,2]
}

HB_emptymatrix7[,1] <- as.vector(as.matrix(HB_timing_subset))
HB_emptymatrix7[,2] <- as.vector(as.matrix(HB_recess_subset))
for (i in 1:nrow(HB_q7_cells)) {
  HB_emptymatrix7[,i+2] <- HB_q7_cells[i,2]<=HB_states[,2] & HB_q7_cells[i,2]>HB_states[-1,2]
}

# a <- vector("list", nrow(q1_cells))
# names(a) <- (q1_cells[,1])
# a[[1]] <- data.frame(matrix(test1[,1:3][test1[,3]==1], ncol=3))
# a[[1]] <- lapply(a[[1]], factor)
# names(a[[1]]) <- c("TIMING", "RECESSION", "FLOOD")
# w <- a[[1]][1]
# x <- a[[1]][2]
# v <- a[[1]][3]
# a[[1]] <- data.frame(w,x,v)

#### Populate empty matrices with evidence ####
HB_evidence1 <- populate.evidence(HB_q1_cells, HB_emptymatrix1)
HB_evidence2 <- populate.evidence(HB_q2_cells, HB_emptymatrix2)
HB_evidence3 <- populate.evidence(HB_q3_cells, HB_emptymatrix3)
HB_evidence4 <- populate.evidence(HB_q4_cells, HB_emptymatrix4)
HB_evidence5 <- populate.evidence(HB_q5_cells, HB_emptymatrix5)
HB_evidence6 <- populate.evidence(HB_q6_cells, HB_emptymatrix6)
HB_evidence7 <- populate.evidence(HB_q7_cells, HB_emptymatrix7)
save(HB_evidence1, file="output/site5/HB_evidence1.Rdata")
save(HB_evidence2, file="output/site5/HB_evidence2.Rdata")
save(HB_evidence3, file="output/site5/HB_evidence3.Rdata")
save(HB_evidence4, file="output/site5/HB_evidence4.Rdata")
save(HB_evidence5, file="output/site5/HB_evidence5.Rdata")
save(HB_evidence6, file="output/site5/HB_evidence6.Rdata")
save(HB_evidence7, file="output/site5/HB_evidence7.Rdata")

#### Calculate recruitment probability by instantiating the BN with evidence ####
system.time(HB_q1_cell_probs <- lapply(HB_evidence1, apply, 1, function(x) cpquery(riparian.fit1, (POTENTIAL=="Y"), evidence=as.list(x), method="lw")))
save(HB_q1_cell_probs, file="output/site5/HB_q1_cell_probs.Rdata")
system.time(HB_q2_cell_probs <- lapply(HB_evidence2, apply, 1, function(x) cpquery(riparian.fit2, (POTENTIAL=="Y"), evidence=as.list(x), method="lw")))
save(HB_q2_cell_probs, file="output/site5/HB_q2_cell_probs.Rdata")
system.time(HB_q3_cell_probs <- lapply(HB_evidence3, apply, 1, function(x) cpquery(riparian.fit3, (POTENTIAL=="Y"), evidence=as.list(x), method="lw")))
save(HB_q3_cell_probs, file="output/site5/HB_q3_cell_probs.Rdata")
system.time(HB_q4_cell_probs <- lapply(HB_evidence4, apply, 1, function(x) cpquery(riparian.fit4, (POTENTIAL=="Y"), evidence=as.list(x), method="lw")))
save(HB_q4_cell_probs, file="output/site5/HB_q4_cell_probs.Rdata")
system.time(HB_q5_cell_probs <- lapply(HB_evidence5, apply, 1, function(x) cpquery(riparian.fit5, (POTENTIAL=="Y"), evidence=as.list(x), method="lw")))
save(HB_q5_cell_probs, file="output/site5/HB_q5_cell_probs.Rdata")
system.time(HB_q6_cell_probs <- lapply(HB_evidence6, apply, 1, function(x) cpquery(riparian.fit6, (POTENTIAL=="Y"), evidence=as.list(x), method="lw")))
save(HB_q6_cell_probs, file="output/site5/HB_q6_cell_probs.Rdata")
system.time(HB_q7_cell_probs <- lapply(HB_evidence7, apply, 1, function(x) cpquery(riparian.fit7, (POTENTIAL=="Y"), evidence=as.list(x), method="lw")))
save(HB_q7_cell_probs, file="output/site5/HB_q7_cell_probs.Rdata")

HB_q1_cell_prob_mn <- t(data.frame(lapply(HB_q1_cell_probs, mean, na.rm=TRUE), check.names=FALSE))
HB_q2_cell_prob_mn <- t(data.frame(lapply(HB_q2_cell_probs, mean, na.rm=TRUE), check.names=FALSE))
HB_q3_cell_prob_mn <- t(data.frame(lapply(HB_q3_cell_probs, mean, na.rm=TRUE), check.names=FALSE))
HB_q4_cell_prob_mn <- t(data.frame(lapply(HB_q4_cell_probs, mean, na.rm=TRUE), check.names=FALSE))
HB_q5_cell_prob_mn <- t(data.frame(lapply(HB_q5_cell_probs, mean, na.rm=TRUE), check.names=FALSE))
HB_q6_cell_prob_mn <- t(data.frame(lapply(HB_q6_cell_probs, mean, na.rm=TRUE), check.names=FALSE))
HB_q7_cell_prob_mn <- t(data.frame(lapply(HB_q7_cell_probs, mean, na.rm=TRUE), check.names=FALSE))
if (exists("HB_q7_cell_prob_mn")) HB_q_all_prob_mn <- rbind(HB_q1_cell_prob_mn, HB_q2_cell_prob_mn, HB_q3_cell_prob_mn, HB_q4_cell_prob_mn, HB_q5_cell_prob_mn, HB_q6_cell_prob_mn, HB_q7_cell_prob_mn)
if (!exists("HB_q7_cell_prob_mn")) HB_q_all_prob_mn <- rbind(HB_q1_cell_prob_mn, HB_q2_cell_prob_mn, HB_q3_cell_prob_mn, HB_q4_cell_prob_mn, HB_q5_cell_prob_mn, HB_q6_cell_prob_mn)

# q1_norm <- normalize(hydro_states, q1_cell_probs)
# q2_norm <- normalize(hydro_states, q2_cell_probs)
# q3_norm <- normalize(hydro_states, q3_cell_probs)
# q4_norm <- normalize(hydro_states, q4_cell_probs)
# q5_norm <- normalize(hydro_states, q5_cell_probs)
# q6_norm <- normalize(hydro_states, q6_cell_probs)

# q1_prob_mn_norm <- q1_cell_prob_mn * q1_norm
# q2_prob_mn_norm <- q2_cell_prob_mn * q2_norm
# q3_prob_mn_norm <- q3_cell_prob_mn * q3_norm
# q4_prob_mn_norm <- q4_cell_prob_mn * q4_norm
# q5_prob_mn_norm <- q5_cell_prob_mn * q5_norm
# q6_prob_mn_norm <- q6_cell_prob_mn * q6_norm

# row.names(q1_prob_mn_norm) <- row.names(q1_cell_prob_mn)
# row.names(q2_prob_mn_norm) <- row.names(q2_cell_prob_mn)
# row.names(q3_prob_mn_norm) <- row.names(q3_cell_prob_mn)
# row.names(q4_prob_mn_norm) <- row.names(q4_cell_prob_mn)
# row.names(q5_prob_mn_norm) <- row.names(q5_cell_prob_mn)
# row.names(q6_prob_mn_norm) <- row.names(q6_cell_prob_mn)

# if (exists("q6_prob_mn_norm")) q_all_mn_norm_e <- rbind(q1_prob_mn_norm, q2_prob_mn_norm, q3_prob_mn_norm, q4_prob_mn_norm, q5_prob_mn_norm, q6_prob_mn_norm)
# if (!exists("q6_prob_mn_norm")) q_all_mn_norm_e <- rbind(q1_prob_mn_norm, q2_prob_mn_norm, q3_prob_mn_norm, q4_prob_mn_norm, q5_prob_mn_norm)

#### Write csv file output ####
# write.csv(q1_prob_mn_norm, file="output/site5/q1_prob_mn_norm_e.csv")
# write.csv(q2_prob_mn_norm, file="output/site5/q2_prob_mn_norm_e.csv")
# write.csv(q3_prob_mn_norm, file="output/site5/q3_prob_mn_norm_e.csv")
# write.csv(q4_prob_mn_norm, file="output/site5/q4_prob_mn_norm_e.csv")
# write.csv(q5_prob_mn_norm, file="output/site5/q5_prob_mn_norm_e.csv")
# write.csv(q6_prob_mn_norm, file="output/site5/q6_prob_mn_norm_e.csv")

# Export File
write.csv(HB_q_all_prob_mn, file="output/site5/HB_q_all_prob_mn.csv")
HB_q_all_prob_mn2 <- read.csv("output/site5/HB_q_all_prob_mn.csv")
colnames(HB_q_all_prob_mn2) <- c("cell","mean_prob")
HB_q_all_prob_mn2 <- HB_q_all_prob_mn2[order(HB_q_all_prob_mn2$cell),]
write.csv(HB_q_all_prob_mn2, file="output/site5/HB_q_all_prob_mn.csv", row.names=FALSE)

# Site 2
# write.csv(q_all_mn_norm_e, file="output/site5/q_all_mn_norm_site5_e.csv")
# q_all_mn_norm_e2 <- read.csv("output/site5/q_all_mn_norm_site5_e.csv")
# colnames(q_all_mn_norm_e2) <- c("cell","mean_prob")
# q_all_mn_norm_e2 <- q_all_mn_norm_e2[order(q_all_mn_norm_e2$cell),]
# write.csv(q_all_mn_norm_e2, file="output/site5/q_all_mn_norm_site5_e.csv", row.names=FALSE)

# Site 3
# write.csv(q_all_mn_norm_e, file="output/site5/q_all_mn_norm_site5_e.csv")
# q_all_mn_norm_e2 <- read.csv("output/site5/q_all_mn_norm_site5_e.csv")
# colnames(q_all_mn_norm_e2) <- c("cell","mean_prob")
# q_all_mn_norm_e2 <- q_all_mn_norm_e2[order(q_all_mn_norm_e2$cell),]
# write.csv(q_all_mn_norm_e2, file="output/site5/q_all_mn_norm_site5_e.csv", row.names=FALSE)

# # Site 4
# write.csv(q_all_mn_norm_e, file="output/site5/q_all_mn_norm_site5_e.csv")
# q_all_mn_norm_e2 <- read.csv("output/site5/q_all_mn_norm_site5_e.csv")
# colnames(q_all_mn_norm_e2) <- c("cell","mean_prob")
# q_all_mn_norm_e2 <- q_all_mn_norm_e2[order(q_all_mn_norm_e2$cell),]
# write.csv(q_all_mn_norm_e2, file="output/site5/q_all_mn_norm_site5_e.csv", row.names=FALSE)

# Site 5
# write.csv(q_all_mn_norm_e, file="output/site5/q_all_mn_norm_site5_e.csv")
# q_all_mn_norm_e2 <- read.csv("output/site5/q_all_mn_norm_site5_e.csv")
# colnames(q_all_mn_norm_e2) <- c("cell","mean_prob")
# q_all_mn_norm_e2 <- q_all_mn_norm_e2[order(q_all_mn_norm_e2$cell),]
# write.csv(q_all_mn_norm_e2, file="output/site5/q_all_mn_norm_site5_e.csv", row.names=FALSE)