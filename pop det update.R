# population detection functions. but correct this time
# pls pls pls

# okay so first we generate the discrete frequencies
# this is the same as before

library(ggplot2)

freq <- function(snp, # number of random SNPs
                 K, # number of discrete subpopulations. uppercase
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
  
 # label <- c()
#  for (i in 1:k){
 #   label <- rbind(label, i)
  #}
  
  #subpop <- cbind(subpop, label) 
  
  # we don't need labels
  
  rownames(subpop) <- NULL
  
  return(subpop)
}


# so... i still have a function that creates individuals from discrete populations

# R is case sensitive so use K as number of populations and k as pop index

individuals <- function(snp, # number of SNPs 
                        n, # number of genotypes generated 
                        k, # subpopulation index (lowercase)
                        subpop) # data frame of allele frequencies
  
{
  out <- c()
  
  for (i in 1:n){
    gen <- rbinom(snp, size = 1, prob = as.numeric(subpop[k, ])) +
      rbinom(snp, size = 1, prob = as.numeric(subpop[k, ]))
    
    out <- rbind(out, gen) 
  }
  
  #out <- cbind(out, rep(k, n)) # ohhh we are adding the label here
  
  row.names(out) <- NULL
  return(out) # genotypes of n individuals from population k
} # this is just for discrete populations!


# so now if we have admixture, it's about how we generate the individuals
# so we need a function that generates n individuals from an admixture 
# it generates random SNPs for the individuals

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


# so now when i generate my samples i need to be more careful

sample_generator <- function(snp,  # how many snps 
                             discrete, # vector of how many from each discrete
                             subpop, # frequencies of discrete
                             mix, #  how many admixed
                             parents, # indexes of admixed
                             alpha,
                             beta)

{
  # how many discrete populations
  K <- nrow(subpop)
  
  data <- c()  # initialise data matrix
  
  # first add the discrete people
  
  for (i in 1:K){  # for each population create label[i] genotypes
    gen <- individuals(snp = snp, n = discrete[i], k = i, subpop = subpop)
    gen <- cbind(gen, rep(i, discrete[i])) # add a label
    data <- rbind(data, gen)
    }
  
  # now add the admixed people
  # in this function we only have 1 admixture but obv it can be extended to more
  
  gen <- admixed(snp, n = mix, alpha, beta, parents, subpop)
  gen <- cbind(gen, rep(73, mix))  # add label 73 to admixed people
  data <- rbind(data, gen)
  
  return(data)
}




