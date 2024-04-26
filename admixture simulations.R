snp <- 10
Fst <- c(0.01, 0.02)
n <- 10
r <- 2
R <- 1.5
K <- 2

subpop <- freq(snp, K, Fst)
subpop

boop <- admixed_status(n, r)
disease <- boop$disease
p <- boop$p

parents <- c(1,2)
diff.freq <- c(0.2, 0.8)

random_admixed(snp, n, subpop, p, disease, parents)
diff_admixed(snp, n, diff.freq, p, disease, parents)
causal_admixed(snp, n, subpop, p, disease, parents, R)



#### FIRST GET THE EIGENVECTORS ####

snp.pc <- 1000
snp <- 2000
Fst <- c(0.08, 0.12)
n <- 800
r <- 3
K <- 2
parents <- c(1,2)
R <- 1.3
diff.freq <- c(0.2, 0.7)


# actually, first i get the proportions so

boop <- admixed_status(n, r)
disease <- boop$disease
p <- boop$p

# now i apply the population structure

x <- pop.struct(snp.pc, K, Fst, n, parents, p, disease)
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

# now we run the random snps glm

random <- glm.random.pc.admixed(snp, K, Fst, n, p, disease, parents)
significant(random$p.values, 0.01)
significant(random$p.values.pca1, 0.01)
significant(random$p.values.pca2, 0.01)
significant(random$p.values.pca10, 0.01)

# now run differentiated snps glm

diff.freq <- c(0.2, 0.7)
diff <- glm.diff.pc.admixed(snp, n, p, disease, parents, diff.freq)
significant(diff$p.values, 0.01)
significant(diff$p.values.pca1, 0.01)
significant(diff$p.values.pca2, 0.01)
significant(diff$p.values.pca10, 0.01)

# now run causal snps glm

R <- 1.5
causal <- glm.causal.pc.admixed(snp, K, Fst, n, p, disease, parents, R)
significant(causal$p.values, 0.01)
significant(causal$p.values.pca1, 0.01)
significant(causal$p.values.pca2, 0.01)
significant(causal$p.values.pca10, 0.01)


#### OKAY NOW WRITE A LOOP FOR ALL ####

# initialise vectors for the p-values
ra <- c()
r.pc <- c()

d <- c()
d.pc <- c()

c <- c()
c.pc <- c()



for (i in 1:3){
  # proportions
  boop <- admixed_status(n, r)
  disease <- boop$disease
  p <- boop$p
  
  # pop structure
  x <- pop.struct(snp.pc, K, Fst, n, parents, p, disease)
  PC1 <- x$PC1
  #PC2 <- x$PC2
  
  # random SNPs
  
  random <- glm.random.pc.admixed(snp, K, Fst, n, p, disease, parents)
  ra <- cbind(ra, significant(random$p.values, 0.01))
  r.pc <- cbind(r.pc, significant(random$p.values.pca1, 0.01))
  #r.pc2 <- cbind(r.pc2, significant(random$p.values.pca2, 0.01))
  
  # differentiated SNPs
  
  diff <- glm.diff.pc.admixed(snp, n, p, disease, parents, diff.freq)
  d <- cbind(d, significant(diff$p.values, 0.01))
  d.pc <- cbind(d.pc, significant(diff$p.values.pca1, 0.01))
  #d.pc2 <- cbind(d.pc2, significant(diff$p.values.pca2, 0.01))
  
  # causal SNPs
  
  causal <- glm.causal.pc.admixed(snp, K, Fst, n, p, disease, parents, R)
  c <- cbind(c, significant(causal$p.values, 0.01))
  c.pc <- cbind(c.pc, significant(causal$p.values.pca1, 0.01))
  #c.pc2 <- cbind(c.pc2, significant(causal$p.values.pca2, 0.01))
}
