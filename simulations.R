#### DEFINING VARIABLES ####

# things i need to define before running any function

case_prop <- c(0.4, 0.6)
control_prop <- c(0.6, 0.4)
snp_pca <- 1000 # number of snps used for pca
snp <- 1000 # number of snps used in testing
t <- 800 # number of people in each group
Fst <- c(0.05, 0.1) # vector of population differentiation
k <- 2 # number of discrete subpops
R <- 1.2 # relative risk
diff.freq <- c(0.2, 0.7) # diff allele frequencies
diff.cause <- c(0.1, 0.2) # diff causal allele frequencies

# and now we run some functions 

#### Population Structure detection ####

# first we run the PS detection bit

x <- pop.struct(snp_pca, k, Fst, t, case_prop, control_prop)

# save the pc's as pc1 and pc2
PC1 <- x$PC1
PC2 <- x$PC2

#### Simulations for the 3 categories of SNPs ####

# Random SNPs

random <- glm.random.pc(snp, t, k, case_prop, control_prop)
significant(random$p.values, 0.01)
significant(random$p.values.pca, 0.01)

# Differentiated SNPs

differentiated <- glm.diff.pc(snp, t, k, case_prop, control_prop, diff.freq)
significant(differentiated$p.values, 0.01)
significant(differentiated$p.values.pca, 0.01)

# Causal SNPs (pls work <3)

causal <- glm.causal.pc(snp, t, k, case_prop, control_prop, R)
significant(causal$p.values, 0.01)
significant(causal$p.values.pca, 0.01)

# Causal Differentiated SNPs (cmon)

causal.cause <- glm.causal.diff(snp, t, k, case_prop, control_prop, R, diff.cause)
significant(causal.cause$p.values, 0.01)
significant(causal.cause$p.values.pca, 0.01)


#### DOING A LOOP TO ASSESS POWER LOSS ####
# these stay the same
case_prop <- c(0.4, 0.6)
control_prop <- c(0.6, 0.4)
snp_pca <- 1000 # number of snps used for pca
snp <- 1000 # number of snps used in testing
t <- 800 # number of people in each group
Fst <- c(0.01, 0.02) # vector of population differentiation
k <- 2 # number of discrete subpops

# these change with the cases
R <- 1.2 # relative risk
diff.cause <- c(0.4, 0.5) # diff causal allele frequencies

# initialise vectors of results
cd <- c()
cd.pc <- c()

# do the loop

for (i in 1:10){
  
  # first find the eigenvectors
  x <- pop.struct(snp_pca, k, Fst, t, case_prop, control_prop)
  
  PC1 <- x$PC1
  PC2 <- x$PC2
  
  # now run the glm
  causal.cause <- glm.causal.diff(snp, t, k, case_prop, control_prop, R, diff.cause)
  cd <- cbind(cd, significant(causal.cause$p.values, 0.01))
  cd.pc <- cbind(cd.pc, significant(causal.cause$p.values.pca, 0.01))
  
}

# also compare them to random SNPs with different relative risks R

c <- c()
c.pc <- c()
R <- 1.5

for (i in 1:10){
  x <- pop.struct(snp_pca, k, Fst, t, case_prop, control_prop)
  
  PC1 <- x$PC1
  PC2 <- x$PC2
  
  causal <- glm.causal.pc(snp, t, k, case_prop, control_prop, R)
  c <- cbind(c, significant(causal$p.values, 0.01))
  c.pc <- cbind(c.pc, significant(causal$p.values.pca, 0.01))
}

# Let's do a for loop and average the results over 10 simulations
c13 <- c(0.682, 0.653, 0.659, 0.691, 0.647, 0.667, 0.679, 0.692, 0.679, 0.679)
c.pc13 <- c(0.722, 0.691, 0.697, 0.724, 0.703, 0.709, 0.705, 0.71, 0.705, 0.714)

c12 <- c(0.451, 0.39, 0.419, 0.414, 0.419, 0.395, 0.406, 0.402, 0.378, 0.429)
c.pc12 <- c(0.396, 0.343, 0.36, 0.372, 0.344, 0.331, 0.339, 0.353, 0.334, 0.364)

# initialise vectors for the p-values
r <- c()
r.pc <- c()
d <- c()
d.pc <- c()
c <- c()
c.pc <- c()


for (i in 1:3){
  
  # pop structure
  
  x <- pop.struct(snp_pca, k, Fst, t, case_prop, control_prop)
  
  PC1 <- x$PC1
  PC2 <- x$PC2
  
  # random
  
  random <- glm.random.pc(snp, t, k, case_prop, control_prop)
  r <- cbind(r, significant(random$p.values, 0.01))
  r.pc <- cbind(r.pc, significant(random$p.values.pca, 0.01))
  
  # differentiated
  
  differentiated <- glm.diff.pc(snp, t, k, case_prop, control_prop, diff.freq)
  d <- cbind(d, significant(differentiated$p.values, 0.01))
  d.pc <- cbind(d.pc, significant(differentiated$p.values.pca, 0.01))
  
  # causal
  
  causal <- glm.causal.pc(snp, t, k, case_prop, control_prop, R)
  c <- cbind(c, significant(causal$p.values, 0.01))
  c.pc <- cbind(c.pc, significant(causal$p.values.pca, 0.01))
  
}

# visualise the p.values

ran <- data.frame(
  x = c("Uncorrected", "Eigenstrat"),
  y = c(mean(r), mean(r.pc))
)

ggplot(ran, aes(x = x, y = y, fill = x)) +
  geom_bar(stat = "identity") +
  theme_minimal() +  # Apply a minimal theme
  labs(
    title = "Comparison of Two Values",
    x = "Category",
    y = "Value"
  )


# for example
data <- random_snps(100, 2, 2, c(0.5, 0.5), c(0.5, 0.5), subpop)
subpop <- freq(100, 2, Fst)
pca.detect(data)




######## EXTRA SIMULATIONS #######

# doing it slightly slower 

# things i need to define before running any function

case_prop <- c(0.5, 0.2, 0.3)
control_prop <- c(0.2, 0.5, 0.3)
snp_pca <- 1000 # number of snps used for pca
snp <- 1000 # number of snps used in testing
t <- 480 # number of people in each group
Fst <- c(0.01, 0.02) # vector of population differentiation
k <- 2 # number of discrete subpops
R <- 1.3 # relative risk
diff.freq <- c(0.2, 0.7) # diff allele frequencies

parents <- c(1, 2)
alpha <- 19
beta <- 13

# first run population structure inference
x <- pop.struct(snp_pca, k, Fst, t, case_prop, control_prop)

PC1 <- x$PC1
PC2 <- x$PC2
PC3 <- x$PC3
PC4 <- x$PC4
PC5 <- x$PC5
PC6 <- x$PC6
PC7 <- x$PC7
PC8 <- x$PC8
PC9 <- x$PC9
PC10 <- x$PC10

# Random SNPs

random <- glm.random.pc(snp, t, k, case_prop, control_prop)
significant(random$p.values, 0.001)
significant(random$p.values.pca1, 0.001)
significant(random$p.values.pca2, 0.001)
significant(random$p.values.pca10, 0.001)

# Differentiated SNPs

differentiated <- glm.diff.pc(snp, t, k, case_prop, control_prop, diff.freq)
significant(differentiated$p.values, 0.01)
significant(differentiated$p.values.pca1, 0.01)
significant(differentiated$p.values.pca2, 0.01)
significant(differentiated$p.values.pca10, 0.01)

# Causal SNPs (pls work <3)

causal <- glm.causal.pc(snp, t, k, case_prop, control_prop, R)
significant(causal$p.values, 0.01)
significant(causal$p.values.pca1, 0.01)
significant(causal$p.values.pca2, 0.01)
significant(causal$p.values.pca10, 0.01)

