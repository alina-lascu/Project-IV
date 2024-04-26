# so i need
# 1. function for ancestral allele frequency
# 2. allele freq for the diverged populations
# 3. admixture proportion
# 4. admixture allele frequency
# 5. generate individuals from allele frequency
# 6. generate the sample (also with population labels!!)
# 7. apply pca (and return the first 4 pcs)
# 8. anova test 
# 9. plot function

library(ggplot2)

######## BASIC SIMULATING FUNCTIONS ########

#### Discrete Subpopulations Allele Frequencies ####

# taken from other file. but make sure to add labels!

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
  
  label <- c()
  for (i in 1:k){
    label <- rbind(label, i)
  }
  
  subpop <- cbind(subpop, label)
  
  rownames(subpop) <- NULL
  
  
  
  
  
 # subpop <- as.data.frame(subpop)
  
 #  adding subpopulation label
  
#  name <- c('A', 'B', 'C', 'D', 'E')
  
 # subpop <- cbind(subpop, name[1:k])
  
#  colnames(subpop)[ncol(subpop)] <- "label"
  
 
  
  return(subpop)
}

#### ADMIXTURE FUNCTIONS ####

# so we have a mix between population A and B...
# i think i will add it as another row to the subpop matrix

admix <- function(snp, # number of snps
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
  
  new <- rep(NA, snp + 1)
  
  for (i in 1:snp){
    new[i] <- p[i] * parent1[i] + (1 - p[i]) * parent2[i]
  }
  
  new[snp + 1] <- nrow(subpop) + 1
               
  # adding the new row to the matrix of frequencies
  
  mixture <- rbind(subpop, new)
  
 # mixture$label[nrow(mixture)] <- paste("mix", parents[1], parents[2])
  
  return(mixture)
  
}


######## GENERATING INDIVIDUALS ########

# this will create n individuals from subpopulation k

individuals <- function(snp, # number of SNPs 
                     n, # number of genotypes generated 
                     k, # subpopulation index
                     mixture) # data frame of allele frequencies
  
{
  out <- c()
  
  for (i in 1:n){
    gen <- rbinom(snp, size = 1, prob = as.numeric(mixture[k, ])) +
      rbinom(snp, size = 1, prob = as.numeric(mixture[k, ]))
    
    out <- rbind(out, gen) 
  }
  
  out <- cbind(out, rep(k, n))
  
  row.names(out) <- NULL
  return(out) # genotypes of n individuals from population k
}


# and now we have to generate the whole sample
# and make sure they are labeled - or actually just use a separate
# vector of labels 

sample_generator <- function(snp, # how many snps we're looking at
                        label, # vector of how many from each subpop
                        subpop) # matrix of underlying frequencies
{
  
  # count how many subpopulations we have
  k <- nrow(subpop)
  
  data <- c()
  
  # go through all subpopulations in order
  
  for (i in 1:k){ # for each subpopulation k, create label[k] genotypes
    x <- individuals(snp = snp, n = label[i], k = i, mixture = subpop)
    if (label[i] != 0)
    data <- rbind(data, x)
    }

  return(data)
}


######## PCA ########

# now that we have the data it's time to analyse it 

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
  
  # return the whole eigen object for now
  
  val <- eig$values[1:4]
  vec <- eig$vectors[,1:4]
  
  return(list(val = val, vec = vec)) 
  
}

# plotting function

plot_pc <- function(data, # the genotype data
                    vec) # the eigenvectors we got from pca
{
  axis1 <- vec[, 1]
  axis2 <- vec[, 2]
  
  #data <- as.data.frame(data)
  #data$label <- c(1,1,2,2,3,3)
  
  # try to label them here...
  
  #plot(axis2 ~ axis1, col = data[, ncol(data)])

  ggplot(as.data.frame(data), aes(x = axis1, y = axis2, col = as.factor(data[, ncol(data)]),
                                  shape = as.factor(data[, ncol(data)]))) +
    geom_point(size = 2) +
    labs(x = "First Eigenvector", y = "Second Eigenvector", col = "Population") +
    scale_color_manual(values = c("#970973", "#F30C74", "#EF6125"), labels = c("A", "B", "C")) +
    scale_shape_manual(values = c(15, 17, 19), guide = "none")+
    guides(color = guide_legend(override.aes = list(shape = c(15, 17, 19))))+
    theme_classic()
  
}


######## ANOVA ########

anova_test <- function(label,
                       e){ # the eigen object
  
  # first separate the eigenvectors
  
  vec <- e$vec
  
  # then add the labels
  
  lab <- c()
  
  # okay and now we fill them in one by one
  
  for (i in 1:length(label))
    lab <- c(lab, rep(i, label[i]))
    
  # now make a dataframe with all 
  
  data_anova <- as.data.frame(cbind(vec, lab))
  colnames(data_anova) <- c("eig1","eig2","eig3",
                            "eig4","lab")
  
  # now do the anova test for each vector
  # and store the p-values separately
  p.val <- c()
  
  for (i in 1:4){
    model <- aov(data_anova[, i] ~ lab, data = data_anova)
    p <- summary(model)[[1]]$"Pr(>F)"[1]
    p.val <- c(p.val, p)
  }
 
  return(p.val)
}



