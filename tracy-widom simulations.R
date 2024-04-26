snp <- 3000
k <- 2 # number of discrete subpopulations
Fst <- c(0.2, 0.3) # how different they are

subpop <- freq(snp, k, Fst)
subpop

parents = c(1, 2)
alpha = 19
beta = 13

mixture <- admix(snp, subpop, parents, alpha, beta)

label <- c(300, 300, 300)
data <- sample_generator(snp, label, mixture)

eval <- pca.detect_tw(data)

tracy_widom(eval)

################

snp <- 3000
k <- 3 # number of discrete subpopulations
Fst <- c(0.05, 0.07, 0.02) # how different they are

subpop <- freq(snp, k, Fst)
subpop

parents <- c(1, 2)
alpha <- 19
beta <- 13

parents <- c(2, 3)
al <- 3.5
be <- 1.5

mixture <- admix(snp, subpop, parents, alpha, beta)

mixture <- admix(snp, mixture, parents, al, be)

label <- c(100, 124, 314, 820, 120)
data <- sample_generator(snp, label, mixture)

eval <- pca.detect_tw(data)

tracy_widom(eval)


#### 2 DISCRETE SUBPOPULATIONS ####

# 2 admixed actually 

snp <- 5000
K <- 2
Fst <- c(0.3, 0.4)

discrete <- c(1, 1) # how many people from the discrete pops
mix <- 300 # how many people from the admixture

parents = c(1, 2)
alpha = 5
beta = 2

res <- c()

for (i in 1:20){
  subpop <- freq(snp, K, Fst)
  data <- sample_generator(snp, discrete, subpop, mix, parents, alpha, beta)
  eval <- pca.detect_tw(data)
  res <- cbind(res, tracy_widom(eval))
}
res


#### 2-ADMIXTURE ####

snp <- 5000
k <- 2
Fst <- c(0.03, 0.04)
label <- c(0, 150, 150)
parents <- c(1,2)
alpha <- 5
beta <- 2

res <- c()

for (i in 1:20){
  subpop <- freq(snp, k, Fst)
  mixture <- admix(snp, subpop, parents, alpha, beta)
  data <- sample_generator(snp, label, mixture)
  eval <- pca.detect_tw(data)
  res <- cbind(res, tracy_widom(eval))
}

res



#### 3-DISCRETE-SUBPOPULATIONS ####

snp <- 3000
k <- 3
Fst <- c(0.2, 0.4, 0.3)
label <- c(150, 150, 150) 
res <- c()
for (i in 1:20){
  subpop <- freq(snp, k, Fst)
  data <- sample_generator(snp, label, subpop)
  eval <- pca.detect_tw(data)
  res <- cbind(res, tracy_widom(eval))
}
res



#### 3 with 1 mix ####


snp <- 3000
K <- 3
Fst <- c(0.02, 0.03, 0.06)
discrete <- c(150, 150, 150)
mix <- 150
parents <- c(1,2)
alpha <- 5
beta <- 2

res <- c()

for (i in 1:10){
  subpop <- freq(snp, k, Fst)
  mixture <- admix(snp, subpop, parents, alpha, beta)
  data <- sample_generator(snp, label, mixture)
  eval <- pca.detect_tw(data)
  res <- cbind(res, tracy_widom(eval))
}

res


#### some more admixture ####

snp <- 3000
k <- 3
Fst <- c(0.02, 0.03, 0.06)
label <- c(150, 150, 0, 150, 150)
parents <- c(1,2)
parents_extra <- c(2,3)
alpha <- 5
beta <- 2
alpha_extra <- 10
beta_extra <- 8

res <- c()

for (i in 1:40){
  subpop <- freq(snp, k, Fst)
  mixture <- admix(snp, subpop, parents, alpha, beta)
  mixture <- admix(snp, mixture, parents_extra, alpha_extra, beta_extra)
  data <- sample_generator(snp, label, mixture)
  eval <- pca.detect_tw(data)
  res <- cbind(res, tracy_widom(eval))
}

res

# admixture with beta(5,2)
# so i first need to simulate the population
# and then i apply pca and tw

snp <- 5000
K <- 2
Fst <- c(0.3, 0.4)
alpha <- 1   # beta(1,1) is like a uniform distribution
beta <- 1
parents <- c(1,2)
n <- 700 # number of admixed individuals

# so i just want a sample of 300 admixed people

subpop <- freq(snp, K, Fst)
data <- admixed(snp, n, alpha, beta, parents, subpop)
eval <- pca.detect_tw(data)
tracy_widom(eval)



###
snp <- 3000
K <- 3
Fst <- c(0.2, 0.4, 0.3)

# include some admixed individuals (300)

label <- c(150, 150, 150)
n <- 300
alpha <- 5
beta <- 2
parents <- c(1,2)
res <- c()
for (i in 1:20){
subpop <- freq(snp, K, Fst)
data <- sample_generator(snp, label, subpop)
data <- data[, -ncol(data)]
mix <- admixed(snp, n, alpha, beta, parents, subpop)
data <- rbind(data, mix)
eval <- pca.detect_tw(data)
res <- cbind(res, tracy_widom(eval))
}

# try just with the admixed (without parent populations)

res <- c()
n <- 300
parents <- c(1,2)
alpha <- 5
beta <- 2
Fst <- c(0.3, 0.4)

for (i in 1:20){
  subpop <- freq(snp, K, Fst)
  data <- admixed(snp, n, alpha, beta, parents, subpop)
  eval <- pca.detect_tw(data)
  res <- cbind(res, tracy_widom(eval))
  }
