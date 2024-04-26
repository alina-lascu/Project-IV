######## CASE / CONTROL SIMULATIONS INCLUDING ADMIXTURE ########

# given the frequencies of the parents
# and given the ancestry risk 
# simulate admixed individuals

# okay so first, i assume to have the subpop already

# this gives me random snps for admixture


freq <- function(snp, # number of random SNPs
                 K, # number of discrete subpopulations
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
  
  for (i in 2:K){
    new <- rbeta(snp, anc*(1-Fst[i])/Fst[i], 
                 (1-anc)*(1-Fst[i])/Fst[i])
    subpop <- rbind(subpop, new)}
  
  rownames(subpop) <- NULL
  
  return(subpop)
}

# this is not for here
admixed <- function(snp, # number of SNPs
                    n, # number of genotypes created
                    alpha, # alpha parameter of beta distribution
                    beta, # beta parameter of beta distribution
                    parents, # vector of indexes of parent populations 
                    subpop) # matrix of frequencies of discrete populations
{
  # first we initialise the things we need
  
  # first we store the parent frequencies separately
  
  parent1 <- subpop[parents[1], ] 
  parent2 <- subpop[parents[2], ]
  
  # this might actually be better bc then i know the proportion for each one so i can use it to show the inconsistencies (like in yu)
  p <- rbeta(n, alpha, beta) # vector of admixture proportion for each individual
  
  out <- c() # initialise out vector
  
  for (i in 1:n){ # so for all individuals
    # first i generate a uniform thing for all the snps
    u <- runif(snp, min = 0, max = 1) # generate vector of unif for each snp
    
    gen <- rep(NA, snp) # initialise genotype of 1 individual
    
    for (j in 1:snp){
      if (u[j] < p[i]) 
        gen[j] <- rbinom(1, size = 1, prob = as.numeric(parent1[j])) +
          rbinom(1, size = 1, prob = as.numeric(parent1[j]))
      else
        gen[j] <- rbinom(1, size = 1, prob = as.numeric(parent2[j])) +
          rbinom(1, size = 1, prob = as.numeric(parent2[j]))
    }  
    
    out <- rbind(out, gen) # adding the individual to the data matrix
    
  }
  
  row.names(out) <- NULL
  
  return(out)
}



### so... this function here only returns to me 2 vectors
# one with the ancestry proportion and one with the disease status.
# then i do separate functions for sample creation

admixed_status <- function(n,   # how many individuals (total so cases + controls combined)
                         r)   # ancestry risk
{
  # find admixture weight for each individual
  
  p <- runif(n, min = 0, max = 1)  
  
  # decide whether they are case or control
  # for each individual, create a vector with the probabilities of disease
  
  prob <- 0.5 * log(r) * r^p / (r - 1)
  
  # initialise vector of disease status for every individual
  disease <- rep(NA, n)
  
  # now create another uniform to compare with
  
  u <- runif(n, min = 0, max = 1)
  
  for (i in 1:n){
    if (u[i] < prob[i]) 
      disease[i] <- 1
    else
      disease[i] <- 0
  }  
  return(list(p = p, disease = disease))
    
}  


#### now we simulate the SNPs normally using the results above

# first generate random SNPs for admixture

random_admixed <- function(snp,  # how many snps
                           n,    # how many individuals
                           subpop,   # parent frequencies
                           p,   # ancestry proportion
                           disease,   # disease status
                           parents)   # parent indexes

{
  # first get the parent frequencies
  parent1 <- subpop[parents[1], ] 
  parent2 <- subpop[parents[2], ]
  
  out <- c() # initialise out vector
  
  for (i in 1:n){ # so for all individuals
    # first i generate a uniform thing for all the snps
    u <- runif(snp, min = 0, max = 1) # generate vector of unif for each snp
    
    gen <- rep(NA, snp) # initialise genotype of 1 individual
    
    for (j in 1:snp){
      if (u[j] < p[i]) 
        gen[j] <- rbinom(1, size = 1, prob = as.numeric(parent1[j])) +
          rbinom(1, size = 1, prob = as.numeric(parent1[j]))
      else
        gen[j] <- rbinom(1, size = 1, prob = as.numeric(parent2[j])) +
          rbinom(1, size = 1, prob = as.numeric(parent2[j]))
    }  
    
    out <- rbind(out, gen) # adding the individual to the data matrix
    
  }
  
  row.names(out) <- NULL
  out <- cbind(out, disease)  # add a column for disease status
  
  return(out)
}

# then differentiated SNPs

diff_admixed <- function(snp,  # how many snps
                           n,    # how many individuals
                           diff.freq,   # vector of different parent frequencies frequencies
                           p,   # ancestry proportion
                           disease,   # disease status
                           parents)   # parent indexes
  
{
  # first get the parent frequencies
  parent1 <- diff.freq[parents[1]]  # this is just a number
  parent2 <- diff.freq[parents[2]]  # this is just a number
  
  out <- c() # initialise out vector
  
  for (i in 1:n){ # so for all individuals
    # first i generate a uniform thing for all the snps
    u <- runif(snp, min = 0, max = 1) # generate vector of unif for each snp
    
    gen <- rep(NA, snp) # initialise genotype of 1 individual
    
    for (j in 1:snp){
      if (u[j] < p[i]) # here we use the one of parent 1
        gen[j] <- rbinom(1, size = 1, prob = as.numeric(parent1)) +
          rbinom(1, size = 1, prob = as.numeric(parent1))
      else
        gen[j] <- rbinom(1, size = 1, prob = as.numeric(parent2)) +
          rbinom(1, size = 1, prob = as.numeric(parent2))
    }  
    
    out <- rbind(out, gen) # adding the individual to the data matrix
    
  }
  
  row.names(out) <- NULL
  out <- cbind(out, disease)  # add a column for disease status
  
  return(out)
}


# and now for causal SNPs - this time we have to pay attention to the disease status


causal_admixed <- function(snp,  # how many snps
                           n,    # how many individuals
                           subpop,   # parent frequencies
                           p,   # ancestry proportion
                           disease,   # disease status
                           parents,   # parent indexes
                           R)   # relative risk
  
{
  # first get the parent frequencies
  
  parent1 <- subpop[parents[1], ] 
  parent2 <- subpop[parents[2], ]
  
  out <- c() # initialise out vector
  
  # internal function for risk
  
  prob_case <- function(x){
    vec <- c( (1-x)^2, 2*R*x*(1-x), R^2 * x^2)/
      ((1-x)^2 + 2*R*x*(1-x) + R^2 * x^2)
    return(vec)
  }
  
  prob_control <- function(x){
    vec <- c( (1-x)^2, 2*x*(1-x), x^2)
    return(vec)
  }
  
  for (i in 1:n){ # so for all individuals
    
    gen <- rep(NA, snp) # initialise genotype of 1 individual
    
    # get their mixed frequency for all SNPs
    
    mix.freq <- rep(NA, snp)  # initialise frequencies
    
    for (j in 1:snp){
      mix.freq[j] <- p[i] * parent1[j] + (1 - p[i]) * parent2[j]
    }
    
    # check if they are cases or controls
    
    if (disease[i] == 1) { # they are cases so we use multiplicative risk model
      for (j in 1:snp) # do it for all snps independently 
        gen[j] <- sample(0:2, size = 1, prob = prob_case(mix.freq[j]))
    }
    else # they are controls so we use normal probabilities
      for (j in 1:snp)
        gen[j] <- sample(0:2, size = 1, prob = prob_control(mix.freq[j]))
    
    out <- rbind(out, gen) # adding the individual to the data matrix
    
  }
  
  row.names(out) <- NULL
  out <- cbind(out, disease)  # add a column for disease status
  
  return(out)
}


# function to infer population structure

pop.struct <- function(snp, # how many snps
                       K, # how many discrete subpopulations
                       Fst, # vector of Fst for each subpop
                       n, # how many individuals in total
                       parents, # parent indexes
                       p, # admixture proportion
                       disease) # disease status (for completeness)
{
  # first generate the snp frequencies
  
  subpop <- freq(snp, K, Fst)
  
  # then save the parents separately
  
  parent1 <- subpop[parents[1], ] 
  parent2 <- subpop[parents[2], ]
  
  # then generate random snps for population detection
  
  data.pc <- random_admixed(snp, n, subpop, p, disease, parents)
  
  # apply PCA
  
  pc <- pca.detect(data.pc)
  
  # this returns the first 10 PCs
  # save them separately
  
  PC1 <- pc[, 1]
  PC2 <- pc[, 2]
  #PC3 <- pc[, 3]
  #PC4 <- pc[, 4]
  #PC5 <- pc[, 5]
  #PC6 <- pc[, 6]
  #PC7 <- pc[, 7]
  #PC8 <- pc[, 8]
  #PC9 <- pc[, 9]
  #PC10 <- pc[, 10]
  
  rm(pc) # for memory purposes
  
  return(list(PC1 = PC1, PC2 = PC2))#, 
              #PC3 = PC3, PC4 = PC4, PC5 = PC5, 
              #PC6 = PC6, PC7 = PC7, PC8 = PC8, PC9 = PC9, PC10 = PC10))
}

# then to run the glms

##### GLM FOR RANDOM SNPS #####

glm.random.pc.admixed <- function(snp, K, Fst, n, p, disease, parents){
  
  p.values <- c()
  
  subpop <- freq(snp, K, Fst)   # first generate the parent frequencies
  
  # then generate random snps
  
  #parent1 <- subpop[parents[1], ] 
  #parent2 <- subpop[parents[2], ]
  
  # then generate random snps for population detection
  
  data <- as.data.frame(random_admixed(snp, n, subpop, p, disease, parents))
  
  # first we do regular regression
  
  for (i in 1:snp){
    
    model <- glm(disease ~ data[, i], data = data, family = binomial)
    d <- dim(coef(summary(model)))
    p <- coef(summary(model))[d[1], d[2]]
    p.values <- cbind(p.values, p)
  }
  
  p.values <- as.numeric(p.values)
  
  # regress out the first PC
  
  p.values.pca1 <- c()
  for (i in 1:snp){
    
    model <- glm(disease ~ PC1 + data[, i], 
                 data = data, family = binomial)
    d <- dim(coef(summary(model)))
    p <- coef(summary(model))[d[1], d[2]]
    p.values.pca1 <- cbind(p.values.pca1, p)
  }
  
  p.values.pca1 <- as.numeric(p.values.pca1)
  
  
  # regress out the first 2 PCs
  
  #p.values.pca2 <- c()
  #for (i in 1:snp){
    
   # model <- glm(disease ~ PC1 + PC2 + data[, i], 
    #             data = data, family = binomial)
    #d <- dim(coef(summary(model)))
    #p <- coef(summary(model))[d[1], d[2]]
    #p.values.pca2 <- cbind(p.values.pca2, p)
  #}
  
  #p.values.pca2 <- as.numeric(p.values.pca2)
  
  # regress out the first 10 PCs
  
  #p.values.pca10 <- c()
  #for (i in 1:snp){
    
   # model <- glm(disease ~ PC1 + PC2 + PC3 + PC4 + 
    #               PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
     #              data[, i], 
      #           data = data, family = binomial)
    #d <- dim(coef(summary(model)))
    #p <- coef(summary(model))[d[1], d[2]]
    #p.values.pca10 <- cbind(p.values.pca10, p)
#  }
  
  #p.values.pca10 <- as.numeric(p.values.pca10)
  
  return(list(p.values = p.values, 
              p.values.pca1 = p.values.pca1))#,
              #p.values.pca2 = p.values.pca2))#,
            #  p.values.pca10 = p.values.pca10))
}  


#### GLM FOR DIFFERENTIATED SNPS ####

glm.diff.pc.admixed <- function(snp, n, p, disease, parents, diff.freq){
  
  p.values <- c()

  # then generate random snps for population detection
  
  data <- as.data.frame(diff_admixed(snp, n, diff.freq, p, disease, parents))
  
  # first we do regular regression
  
  for (i in 1:snp){
    
    model <- glm(disease ~ data[, i], data = data, family = binomial)
    d <- dim(coef(summary(model)))
    p <- coef(summary(model))[d[1], d[2]]
    p.values <- cbind(p.values, p)
  }
  
  p.values <- as.numeric(p.values)
  
  # regress out the first PC
  
  p.values.pca1 <- c()
  for (i in 1:snp){
    
    model <- glm(disease ~ PC1 + data[, i], 
                 data = data, family = binomial)
    d <- dim(coef(summary(model)))
    p <- coef(summary(model))[d[1], d[2]]
    p.values.pca1 <- cbind(p.values.pca1, p)
  }
  
  p.values.pca1 <- as.numeric(p.values.pca1)
  
  
  # regress out the first 2 PCs
  
  #p.values.pca2 <- c()
  #for (i in 1:snp){
    
   # model <- glm(disease ~ PC1 + PC2 + data[, i], 
    #             data = data, family = binomial)
    #d <- dim(coef(summary(model)))
    #p <- coef(summary(model))[d[1], d[2]]
    #p.values.pca2 <- cbind(p.values.pca2, p)
  #}
  
  #p.values.pca2 <- as.numeric(p.values.pca2)
  
  # regress out the first 10 PCs
  
  #p.values.pca10 <- c()
  #for (i in 1:snp){
    
   # model <- glm(disease ~ PC1 + PC2 + PC3 + PC4 + 
    #               PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
     #              data[, i], 
      #           data = data, family = binomial)
    #d <- dim(coef(summary(model)))
    #p <- coef(summary(model))[d[1], d[2]]
    #p.values.pca10 <- cbind(p.values.pca10, p)
  #}
  
  #p.values.pca10 <- as.numeric(p.values.pca10)
  
  return(list(p.values = p.values, 
              p.values.pca1 = p.values.pca1))#,
             # p.values.pca2 = p.values.pca2))#,
              #p.values.pca10 = p.values.pca10))
}  

#### GLM FOR CAUSAL SNPS ####

glm.causal.pc.admixed <- function(snp, K, Fst, n, p, disease, parents, R){
  
  p.values <- c()
  
  subpop <- freq(snp, K, Fst)   # first generate the parent frequencies
  
  # then generate random snps
  
  #parent1 <- subpop[parents[1], ] 
  #parent2 <- subpop[parents[2], ]
  
  # then generate random snps for population detection
  
  data <- as.data.frame(causal_admixed(snp, n, subpop, p, disease, parents, R))
  
  # first we do regular regression
  
  for (i in 1:snp){
    
    model <- glm(disease ~ data[, i], data = data, family = binomial)
    d <- dim(coef(summary(model)))
    p <- coef(summary(model))[d[1], d[2]]
    p.values <- cbind(p.values, p)
  }
  
  p.values <- as.numeric(p.values)
  
  # regress out the first PC
  
  p.values.pca1 <- c()
  for (i in 1:snp){
    
    model <- glm(disease ~ PC1 + data[, i], 
                 data = data, family = binomial)
    d <- dim(coef(summary(model)))
    p <- coef(summary(model))[d[1], d[2]]
    p.values.pca1 <- cbind(p.values.pca1, p)
  }
  
  p.values.pca1 <- as.numeric(p.values.pca1)
  
  
  # regress out the first 2 PCs
  
  #p.values.pca2 <- c()
  #for (i in 1:snp){
    
   # model <- glm(disease ~ PC1 + PC2 + data[, i], 
    #             data = data, family = binomial)
    #d <- dim(coef(summary(model)))
    #p <- coef(summary(model))[d[1], d[2]]
    #p.values.pca2 <- cbind(p.values.pca2, p)
  #}
  
  #p.values.pca2 <- as.numeric(p.values.pca2)
  
  # regress out the first 10 PCs
  
  #p.values.pca10 <- c()
  #for (i in 1:snp){
    
   # model <- glm(disease ~ PC1 + PC2 + PC3 + PC4 + 
    #               PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
     #              data[, i], 
      #           data = data, family = binomial)
    #d <- dim(coef(summary(model)))
    #p <- coef(summary(model))[d[1], d[2]]
    #p.values.pca10 <- cbind(p.values.pca10, p)
  #}
  
  #p.values.pca10 <- as.numeric(p.values.pca10)
  
  return(list(p.values = p.values, 
              p.values.pca1 = p.values.pca1))#,
              #p.values.pca2 = p.values.pca2))#,
        #      p.values.pca10 = p.values.pca10))
}  
