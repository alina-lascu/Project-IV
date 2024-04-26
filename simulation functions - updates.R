######## FREQUENCIES SIMULATIONS ########

freq <- function(snp, # number of random SNPs
                 k, # number of discrete subpopulations
                 Fst,# vector of Fst for each subpopulation 
                 ancestral = FALSE) # if true it prints ancestral frequencies
  
{
  anc <- runif(snp, min = 0.1, max = 0.9) # generate ancestral SNP frequencies
  
  if (ancestral) {
    cat('ancestral frequencies')
    print(anc)
    
  }
  subpop <- rbeta(snp, anc*(1-Fst[1])/Fst[1], 
                  (1-anc)*(1-Fst[1])/Fst[1]) # for first subpop
  
  for (i in 2:k){
    new <- rbeta(snp, anc*(1-Fst[i])/Fst[i], 
                 (1-anc)*(1-Fst[i])/Fst[i])
    subpop <- rbind(subpop, new)}
  
  rownames(subpop) <- NULL
  
  return(subpop)
}



######## NEED A NEW PCA DETECT FUNCTION ########

# first the auxiliary functions

normalisation <- function (x) # assuming no missing data!
{
  x <- as.vector(x)
  
  col.mean <- mean(x)
  n <- length(x)  
  p <- (sum(x) + 1) / (2 * (n + 1))  # estimate of the underlying allele freq
  
  y <- (x - col.mean) / sqrt (p * (1 - p)) # normalise the column
  
  return(y)
} 


# now a function that returns PCs

pca.detect <- function(data) # data matrix
{
  n <- nrow(data) # number of samples
  p <- ncol(data) - 1 # number of snps (remove the outcome column)
  x <- data[, -(p + 1)] # remove the outcome column
  Sn <- matrix(data = 0, nrow = n, ncol = n) # initialise matrix
  
  z <- apply(x, 2, normalisation) # normalised matrix
  
  rm(x) # we don't need to remember this
  
  Sn <- Sn + (z %*% t(z)) # compute unscaled matrix
  
  rm(z) # we don't need to remember this
  
  Sn <- Sn / p # finish computing Sn
  
  # now we do pca
  
  eig <- eigen(Sn)  # do eigenanalysis
  
  # return the first 10 eigenvectors here
  
  vec <- eig$vectors[,1:10]
  
  return(vec) 
  
}



#### INCLUDING ADMIXTURE ####

admix_update <- function(snp, # number of snps
                  subpop, # matrix of allele frequencies for all discrete subpop
                  parents, # index of the 2 populations mixed
                  alpha, # alpha parameter of beta distribution
                  beta) # beta parameter of beta distribution
{
  # first we create a vector of admixture proportions p
  # where p ~ Beta (alpha, beta). and we need snp of them
  
  p <- rbeta(snp, alpha, beta)
  
  # now create the frequencies of the admixed population 
  # weighted average between parents
  # for each snp one by one
  
  # get all allele frequencies of the parents
  
  #parent1 <- as.numeric(subpop[subpop[, snp + 1] == parents[1], 1:snp])
  #parent2 <- as.numeric(subpop[subpop[, snp + 1] == parents[2], 1:snp])
  
  parent1 <- subpop[parents[1], 1:snp]
  parent2 <- subpop[parents[2], 1:snp]
  
  #new <- rep(NA, snp + 1)
  
  # new <- rep(NA, snp + 1)
  new <- rep(NA, snp)
  
  for (i in 1:snp){
    new[i] <- p[i] * parent1[i] + (1 - p[i]) * parent2[i]
  }
  
  # new[snp + 1] <- nrow(subpop) + 1
  
  # adding the new row to the matrix of frequencies
  
  mixture <- rbind(subpop, new)
  
  # mixture$label[nrow(mixture)] <- paste("mix", parents[1], parents[2])
  
  return(mixture)
  
}

#### UPDATED VERSION ####

pop.struct <- function(snp, # how many snps
                       k, # how many discrete subpopulations
                       Fst, # vector of Fst for each subpop
                       t, # how many samples in each group
                       case_prop, # prop of each subpop in case
                       control_prop) # prop of each subpoop in control
{
  # first generate the snp frequencies
  
  subpop <- freq(snp, k, Fst)
  
  # include admixture
  
  subpop <- admix_update(snp, subpop, parents, alpha, beta)
  
  # then generate the genotypes using random_snps
  
  data.pc <- random_snps(snp, t, k + 1, case_prop, control_prop, subpop)
  
  # apply pca on data without last column
  
  # i think i already did this tho? when i set x in the pca.detect fct?
  pc <- pca.detect(data.pc)
  
  # this returns the first 10 PCs
  # save them separately
  
  PC1 <- pc[, 1]
  PC2 <- pc[, 2]
  PC3 <- pc[, 3]
  PC4 <- pc[, 4]
  PC5 <- pc[, 5]
  PC6 <- pc[, 6]
  PC7 <- pc[, 7]
  PC8 <- pc[, 8]
  PC9 <- pc[, 9]
  PC10 <- pc[, 10]
  
  rm(pc) # for memory purposes
  
  return(list(PC1 = PC1, PC2 = PC2, PC3 = PC3, PC4 = PC4, PC5 = PC5, 
              PC6 = PC6, PC7 = PC7, PC8 = PC8, PC9 = PC9, PC10 = PC10))
}


######## MODIFY LOGISTIC MODELS TO TRY THE DIFFERENT NUMBER OF AXES ########

#### Random SNPs ####

glm.random.pc <- function(snp, t, k, case_prop, control_prop){
  
  p.values <- c()
  
  subpop <- freq(snp, k, Fst)
  
  # include admixture
  
  subpop <- admix_update(snp, subpop, parents, alpha, beta)
  
  data <- as.data.frame(random_snps(snp, t, k + 1, case_prop, control_prop,
                                    subpop))
  
  # first we do regular regression
  
  for (i in 1:snp){
    
    model <- glm(disease_status ~ data[, i], data = data, family = binomial)
    d <- dim(coef(summary(model)))
    p <- coef(summary(model))[d[1], d[2]]
    p.values <- cbind(p.values, p)
  }
  
  p.values <- as.numeric(p.values)
  
  # regress out the first PC
  
  p.values.pca1 <- c()
  for (i in 1:snp){
    
    model <- glm(disease_status ~ PC1 + data[, i], 
                 data = data, family = binomial)
    d <- dim(coef(summary(model)))
    p <- coef(summary(model))[d[1], d[2]]
    p.values.pca1 <- cbind(p.values.pca1, p)
  }
  
  p.values.pca1 <- as.numeric(p.values.pca1)
  
  
  # regress out the first 2 PCs
  
  p.values.pca2 <- c()
  for (i in 1:snp){
    
    model <- glm(disease_status ~ PC1 + PC2 + data[, i], 
                 data = data, family = binomial)
    d <- dim(coef(summary(model)))
    p <- coef(summary(model))[d[1], d[2]]
    p.values.pca2 <- cbind(p.values.pca2, p)
  }
  
  p.values.pca2 <- as.numeric(p.values.pca2)
  
  # regress out the first 10 PCs
  
  p.values.pca10 <- c()
  for (i in 1:snp){
    
    model <- glm(disease_status ~ PC1 + PC2 + PC3 + PC4 + 
                   PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                   data[, i], 
                 data = data, family = binomial)
    d <- dim(coef(summary(model)))
    p <- coef(summary(model))[d[1], d[2]]
    p.values.pca10 <- cbind(p.values.pca10, p)
  }
  
  p.values.pca10 <- as.numeric(p.values.pca10)
  
  return(list(p.values = p.values, 
              p.values.pca1 = p.values.pca1,
              p.values.pca2 = p.values.pca2,
              p.values.pca10 = p.values.pca10))
}  


#### Diffferentiated SNPs ####

glm.diff.pc <- function(snp, t, k, case_prop, control_prop, diff.freq){
  
  p.values <- c()
  data <- as.data.frame(diff_snps(snp, t, k, case_prop, control_prop,
                                  diff.freq))
  # first do regular
  
  for (i in 1:snp){
    
    model <- glm(disease_status ~ data[, i], data = data, family = binomial)
    d <- dim(coef(summary(model)))
    p <- coef(summary(model))[d[1], d[2]]
    p.values <- cbind(p.values, p)
  }
  
  p.values <- as.numeric(p.values)
  
  # regress out the first PC
  
  p.values.pca1 <- c()
  for (i in 1:snp){
    
    model <- glm(disease_status ~ PC1 + data[, i], 
                 data = data, family = binomial)
    d <- dim(coef(summary(model)))
    p <- coef(summary(model))[d[1], d[2]]
    p.values.pca1 <- cbind(p.values.pca1, p)
  }
  
  p.values.pca1 <- as.numeric(p.values.pca1)
  
  
  # regress out the first 2 PCs
  
  p.values.pca2 <- c()
  for (i in 1:snp){
    
    model <- glm(disease_status ~ PC1 + PC2 + data[, i], 
                 data = data, family = binomial)
    d <- dim(coef(summary(model)))
    p <- coef(summary(model))[d[1], d[2]]
    p.values.pca2 <- cbind(p.values.pca2, p)
  }
  
  p.values.pca2 <- as.numeric(p.values.pca2)
  
  # regress out the first 10 PCs
  
  p.values.pca10 <- c()
  for (i in 1:snp){
    
    model <- glm(disease_status ~ PC1 + PC2 + PC3 + PC4 + 
                   PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                   data[, i], 
                 data = data, family = binomial)
    d <- dim(coef(summary(model)))
    p <- coef(summary(model))[d[1], d[2]]
    p.values.pca10 <- cbind(p.values.pca10, p)
  }
  
  p.values.pca10 <- as.numeric(p.values.pca10)
  
  return(list(p.values = p.values, 
              p.values.pca1 = p.values.pca1,
              p.values.pca2 = p.values.pca2,
              p.values.pca10 = p.values.pca10))
  
  
}  

#### Causal SNPs ####

glm.causal.pc <- function(snp, t, k, case_prop, control_prop, R){
  
  p.values <- c()
  data <- as.data.frame(causal_snps(snp, t, k, case_prop, control_prop,
                                    R))
  
  for (i in 1:snp){
    
    model <- glm(disease_status ~ data[, i], data = data, family = binomial)
    d <- dim(coef(summary(model)))
    p <- coef(summary(model))[d[1], d[2]]
    p.values <- cbind(p.values, p)
  }
  
  p.values <- as.numeric(p.values)
  
  # regress out the first PC
  
  p.values.pca1 <- c()
  for (i in 1:snp){
    
    model <- glm(disease_status ~ PC1 + data[, i], 
                 data = data, family = binomial)
    d <- dim(coef(summary(model)))
    p <- coef(summary(model))[d[1], d[2]]
    p.values.pca1 <- cbind(p.values.pca1, p)
  }
  
  p.values.pca1 <- as.numeric(p.values.pca1)
  
  
  # regress out the first 2 PCs
  
  p.values.pca2 <- c()
  for (i in 1:snp){
    
    model <- glm(disease_status ~ PC1 + PC2 + data[, i], 
                 data = data, family = binomial)
    d <- dim(coef(summary(model)))
    p <- coef(summary(model))[d[1], d[2]]
    p.values.pca2 <- cbind(p.values.pca2, p)
  }
  
  p.values.pca2 <- as.numeric(p.values.pca2)
  
  # regress out the first 10 PCs
  
  p.values.pca10 <- c()
  for (i in 1:snp){
    
    model <- glm(disease_status ~ PC1 + PC2 + PC3 + PC4 + 
                   PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                   data[, i], 
                 data = data, family = binomial)
    d <- dim(coef(summary(model)))
    p <- coef(summary(model))[d[1], d[2]]
    p.values.pca10 <- cbind(p.values.pca10, p)
  }
  
  p.values.pca10 <- as.numeric(p.values.pca10)
  
  return(list(p.values = p.values, 
              p.values.pca1 = p.values.pca1,
              p.values.pca2 = p.values.pca2,
              p.values.pca10 = p.values.pca10))
  
}  
