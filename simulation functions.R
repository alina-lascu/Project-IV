# attempt at simulation

######## RANDOM SNPs - base functions ########

#### Subpop allele frequencies ####
# function for generating subpop allele frequencies
# returns a matrix where each row represents a subpopulation
# and each column a snp

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

#### Data simulation ####

## function to simulate n individuals from subpopulation k

subjects <- function(snp, # number of SNPs 
                     n, # number of genotypes generated
                     k,
                     subpop) # subpopulation index
                     
{
  out <- c()
 
  for (i in 1:n){
    gen <- rbinom(snp, size = 1, prob = subpop[k, ]) +
      rbinom(snp, size = 1, prob = subpop[k, ])
    
    out <- rbind(out, gen) 
  }
  
  row.names(out) <- NULL
  return(out) # genotypes of n individuals from population k
}


## create data matrix with all samples


random_snps <- function(snp, # how many snps we're looking at
                        t, # how many samples in each group
                        k, # how many discrete subpopulations
                        case_prop, # prop of subpops in cases
                        control_prop,
                        subpop) # prop of subpops in controls
{
  
# for cases first
cases = c()
for (i in 1:k){
  x <- subjects(snp = snp, n = t * case_prop[i], k = i, subpop = subpop)
  cases <- rbind(cases, x)
} # this gives the cases 
disease_status <- rep(1, t)
cases <- cbind(cases, disease_status)

# now for controls
controls = c()
for (i in 1:k){
  x <- subjects(snp = snp, n = t * control_prop[i], k = i,
                subpop = subpop)
  controls <- rbind(controls, x)
} # this gives the controls 
disease_status <- rep(0, t)
controls <- cbind(controls, disease_status)

data <- rbind(cases, controls)
return(data)
}



######## DIFFERENTIATED SNPS ########

# these snps have very different frequencies between the subpops
# e.g. (0.8, 0.2) for subpopulations k = (1,2)

# function to generate n individuals from subpopulation k
# for set frequencies of snps

diff <- function(snp, # number of SNPs 
                 n, # number of genotypes generated
                 k, # subpopulation index
                 diff.freq) # vector of set allele freq for the subpops
  
{
  out <- c()
  
  for (i in 1:n){
    gen <- rbinom(snp, size = 1, prob = diff.freq[k]) +
      rbinom(snp, size = 1, prob = diff.freq[k])
    
    out <- rbind(out, gen) 
  }
  
  row.names(out) <- NULL
  return(out) # genotypes of n individuals from population k
}


## function to get data matrix for diff snps

diff_snps <- function(snp, # how many snps we're looking at
                        t, # how many samples in each group
                        k, # how many discrete subpopulations
                        case_prop, # prop of subpops in cases
                        control_prop, # prop of subpops in controls
                      diff.freq) # set allele freq in subpops 
{
  
  # for cases first
  cases = c()
  for (i in 1:k){
    x <- diff(snp = snp, n = t * case_prop[i], k = i, 
              diff.freq = diff.freq) # e.g. diff.freq = c(0.2, 0.7)
    cases <- rbind(cases, x)
  } # this gives the cases 
  disease_status <- rep(1, t)
  cases <- cbind(cases, disease_status)
  
  # now for controls
  controls = c()
  for (i in 1:k){
    x <- diff(snp = snp, n = t * control_prop[i], k = i,
              diff.freq = diff.freq)
    controls <- rbind(controls, x)
  } # this gives the controls 
  disease_status <- rep(0, t)
  controls <- cbind(controls, disease_status)
  
  data <- rbind(cases, controls)
  return(data)
}


######## CAUSAL SNPS ########
# given the subpop allele frequency 
# using multiplicative risk model

# function to simulate n case individuals from subpop k
# for causal alleles

causal_case <- function(snp, # number of SNPs 
                 n, # number of genotypes generated
                 k, # subpopulation index
                 R,
                 subpop) # relative risk for the causal allele
  
{
  out <- c()
  
  # internal function to return the genotype prob for each snp
  prob <- function(x){
    vec <- c( (1-x)^2, 2*R*x*(1-x), R^2 * x^2)/
      ((1-x)^2 + 2*R*x*(1-x) + R^2 * x^2)
    return(vec)
  }

  
  for (i in 1:n){ # do it for n individuals
    ind <- c()
    for (j in 1:snp) { # do it for all snps
      s <- sample(0:2, size = 1, prob = prob(subpop[k, j]))
      ind <- cbind(ind, s)
    }
    
    out <- rbind(out, ind) 
  }
  
  colnames(out) <- NULL
  row.names(out) <- NULL
  return(out) # genotypes of n individuals from population k
}

# and for controls we just simulate random snps
# so we don't need another function


causal_snps <- function(snp, # how many snps we're looking at
                        t, # how many samples in each group
                        k, # how many discrete subpopulations
                        case_prop, # prop of subpops in cases
                        control_prop, # prop of subpops in controls
                        R) # relative risk for the causal allele 
{
  # first we simulate underlying allele frequencies
  
  subpop <- freq(snp, k, Fst)
  
  # for cases first - using multiplicative risk model
  cases = c()
  for (i in 1:k){
    x <- causal_case(snp = snp, n = t * case_prop[i], k = i,
                     R = R, subpop = subpop)
    cases <- rbind(cases, x)
  } # this gives the cases 
  disease_status <- rep(1, t)
  cases <- cbind(cases, disease_status)
  
  # now for controls - just random snps like before
  controls = c()
  for (i in 1:k){
    x <- subjects(snp = snp, n = t * control_prop[i], k = i,
                  subpop = subpop)
    controls <- rbind(controls, x)
  } # this gives the controls 
  disease_status <- rep(0, t)
  controls <- cbind(controls, disease_status)
  
  data <- rbind(cases, controls)
  return(data)
}


# causal-differentiated SNPs
causal_diff_snps <- function(snp, # how many snps we're looking at
                        t, # how many samples in each group
                        k, # how many discrete subpopulations
                        case_prop, # prop of subpops in cases
                        control_prop, # prop of subpops in controls
                        R, # relative risk for the causal allele 
                        diff.cause) # diff frequencies for causal SNPs
{
  # make a "subpop" with the different frequencies
  
  subpop <- rep(diff.cause[1], snp)
  subpop <- rbind(subpop, rep(diff.cause[2], snp))
  
  # for cases first - using multiplicative risk model
  cases = c()
  for (i in 1:k){
    x <- causal_case(snp = snp, n = t * case_prop[i], k = i,
                     R = R, subpop = subpop)
    cases <- rbind(cases, x)
  } # this gives the cases 
  disease_status <- rep(1, t)
  cases <- cbind(cases, disease_status)
  
  # now for controls - just random snps like before
  controls = c()
  for (i in 1:k){
    x <- subjects(snp = snp, n = t * control_prop[i], k = i,
                  subpop = subpop)
    controls <- rbind(controls, x)
  } # this gives the controls 
  disease_status <- rep(0, t)
  controls <- cbind(controls, disease_status)
  
  data <- rbind(cases, controls)
  return(data)
}


######## AUXILIAR FUNCTIONS ########

# function to get proportion of significant p-values
# this is the type 1 error / power 
significant <- function(p.values, # vector of all p-values
                        alpha) # significance level
{
  significant <- sum(p.values <= alpha) / length(p.values)
  return(significant)
}

# function for qq-plot

qq <- function(x, ...) {
  plot((1:length(x))/(1 + length(x)), sort(x),
       type = 'l', ...)
  abline(0, 1, col = "red")
}



######## PCA DETECTION FUNCTIONS ########
#### Auxiliary Functions ####

# first the function that normalises each column separately

normalisation <- function (x) # assuming no missing data!
{
  x <- as.vector(x)
  
  col.mean <- mean(x)
  n <- length(x)  
  p <- (sum(x) + 1) / (2 * (n + 1))  # estimate of the underlying allele freq
  
  y <- (x - col.mean) / sqrt (p * (1 - p)) # normalise the column
  
  return(y)
} 

# now a function that returns PC's

pca.detect <- function(data) # data matrix
{
  n <- nrow(data) # number of samples
  p <- ncol(data) - 1 # number of snps
  x <- data[, -(p + 1)]
  Sn <- matrix(data = 0, nrow = n, ncol = n) # initialise matrix
  
  z <- apply(x, 2, normalisation) # normalised matrix
  
  rm(x) # we don't need to remember this
  
  Sn <- Sn + (z %*% t(z)) # compute unscaled matrix
  
  rm(z) # we don't need to remember this
  
  Sn <- Sn / p # finish computing Sn
  
  # now we do pca
  
  eig <- eigen(Sn)  # do eigenanalysis
  vec <- eig$vectors[,1:2] # get the first 2 eigenvectors
  
  return(vec)
  
}

#### Create PS data and detect PCA ####

# function to generate snps for pop structure and apply pca on it

pop.struct <- function(snp, # how many snps
                       k, # how many discrete subpopulations
                       Fst, # vector of Fst for each subpop
                       t, # how many samples in each group
                       case_prop, # prop of each subpop in case
                       control_prop) # prop of each subpoop in control
{
  # first generate the snp frequencies
  
  subpop <- freq(snp, k, Fst)
  
  # then generate the genotypes using random_snps
  
  data.pc <- random_snps(snp, t, k, case_prop, control_prop, subpop)
  
  # apply pca on data without last column
  
  # i think i already did this tho? when i set x in the pca.detect fct?
  pc <- pca.detect(data.pc)
  
  # save the first 2 pc's
  
  PC1 <- pc[,1]
  PC2 <- pc[,2]
  
  rm(pc) # for memory purposes
  
  #p <- plot(PC2 ~ PC1, col = data.pc[,'disease_status'] + 3)
  
  return(list(PC1 = PC1, PC2 = PC2))
}


######## LOGISTIC MODELS ########

#### Random SNPs ####

glm.random.pc <- function(snp, t, k, case_prop, control_prop){
  
  p.values <- c()
  
  subpop <- freq(snp, k, Fst)
  
  data <- as.data.frame(random_snps(snp, t, k, case_prop, control_prop,
                                    subpop))
  
  # first we do regular regression
  
  for (i in 1:snp){
    
    model <- glm(disease_status ~ data[, i], data = data, family = binomial)
    d <- dim(coef(summary(model)))
    p <- coef(summary(model))[d[1], d[2]]
    p.values <- cbind(p.values, p)
  }
  
  p.values <- as.numeric(p.values)
  
  # regress out the PC's as found above
  
  p.values.pca <- c()
  for (i in 1:snp){
    
    model <- glm(disease_status ~ PC1 + PC2 + data[, i], 
                 data = data, family = binomial)
    d <- dim(coef(summary(model)))
    p <- coef(summary(model))[d[1], d[2]]
    p.values.pca <- cbind(p.values.pca, p)
  }
  
  p.values.pca <- as.numeric(p.values.pca)
  
  
  return(list(p.values = p.values, 
              p.values.pca = p.values.pca))
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
  
  # now we regress out the PC's from above
  
  p.values.pca <- c()
  for (i in 1:snp){
    
    model <- glm(disease_status ~ PC1 + PC2 + data[, i], 
                 data = data, family = binomial)
    d <- dim(coef(summary(model)))
    p <- coef(summary(model))[d[1], d[2]]
    p.values.pca <- cbind(p.values.pca, p)
  }
  
  p.values.pca <- as.numeric(p.values.pca)
  
  
  return(list(p.values = p.values, 
              p.values.pca = p.values.pca))
  
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
  
  
  # regress out the PC's found earlier
  
  p.values.pca <- c()
  for (i in 1:snp){
    
    model <- glm(disease_status ~ PC1 + PC2 + data[, i], 
                 data = data, family = binomial)
    d <- dim(coef(summary(model)))
    p <- coef(summary(model))[d[1], d[2]]
    p.values.pca <- cbind(p.values.pca, p)
  }
  
  p.values.pca <- as.numeric(p.values.pca)
  
  
  return(list(p.values = p.values, 
              p.values.pca = p.values.pca))
  
}  


#### Causal Differentiated SNPs GLM ####

glm.causal.diff <- function(snp, t, k, case_prop, control_prop, R, diff.cause){
  
  p.values <- c()
  data <- as.data.frame(causal_diff_snps(snp, t, k, case_prop, control_prop,
                                    R, diff.cause))
  
  for (i in 1:snp){
    
    model <- glm(disease_status ~ data[, i], data = data, family = binomial)
    d <- dim(coef(summary(model)))
    p <- coef(summary(model))[d[1], d[2]]
    p.values <- cbind(p.values, p)
  }
  
  p.values <- as.numeric(p.values)
  
  
  # regress out the PC's found earlier
  
  p.values.pca <- c()
  for (i in 1:snp){
    
    model <- glm(disease_status ~ PC1 + PC2 + data[, i], 
                 data = data, family = binomial)
    d <- dim(coef(summary(model)))
    p <- coef(summary(model))[d[1], d[2]]
    p.values.pca <- cbind(p.values.pca, p)
  }
  
  p.values.pca <- as.numeric(p.values.pca)
  
  
  return(list(p.values = p.values, 
              p.values.pca = p.values.pca))
  
}  
