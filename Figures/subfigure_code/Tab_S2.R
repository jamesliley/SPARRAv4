###############################################################
#### Code used to extract the numbers provided in Table S2 ####
###############################################################

# Input EHR tables used for SPARRAv4
inputs1 <- c("AE2","PIS", "SMR00","SMR01","SMR01E","SMR04", "SystemWatch")
inputs2 <- c("SPARRALTC", "deaths")
inputs <- c(inputs1, inputs2)
# DSH extract with relevant information
fnum <- readLines("Analysis/full_model/Description/misc.txt")

Rs <- rep(0,length(inputs)) # Number of records per source
Is <- rep(0,length(inputs)) # Number of individuals per source
for (i in 1:length(inputs)) {
  w <- which(fnum == inputs[i])
  Rs[i] <- as.numeric(fnum[w+1])
  Is[i] <- as.numeric(fnum[w+2])
}

# Raw EHR
df <- data.frame("Input" = inputs1,
                 "Records" = Rs[inputs %in% inputs1],
                 "Individuals" = Is[inputs %in% inputs1])
df[order(df$Records, decreasing = TRUE),]

# Total number of records across all Raw EHR databases
sum(df$Records)

# Other
df <- data.frame("Input" = inputs2,
                 "Records" = Rs[inputs %in% inputs2],
                 "Individuals" = Is[inputs %in% inputs2])
df[order(df$Records, decreasing = TRUE),]

# Total number of individuals across all databases
# Need to verify that this number did not include individuals
# that were not in the tables above, but that were part of
# SPARRALTC or deaths
as.numeric(fnum[which(fnum == "Total patients")+1])