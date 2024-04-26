snp <- 5000 # we have 10 snps
k <- 2 # number of discrete subpopulations
Fst <- c(0.01, 0.02) # how different they are

subpop <- freq(snp, k, Fst)
subpop

parents = c(1, 2)
alpha = 3.5
beta = 1.5

mixture <- admix(snp, subpop, parents, alpha, beta)

n <- 10
out <- individuals(snp, n, 1, mixture)

label <- c(200, 200) # will give us 20 from each population
data <- sample_generator(snp, label, mixture)

e <- pca.detect(data)


plot_pc(data, e$vec)


######## Scenarios ######## 
# maybe later i can look at how changing the parameters
# would change this... i dont know if it's too 
# many things

#### 2-discrete subpopulations ####

snp <- 5000
k <- 2 # number of discrete subpopulations
Fst <- c(0.01, 0.02) # how different they are

subpop <- freq(snp, k, Fst)
subpop

parents = c(1, 2)
alpha = 19
beta = 13

mixture <- admix(snp, subpop, parents, alpha, beta)


label <- c(400, 300, 0) # will give us 20 from each population
data <- sample_generator(snp, label, mixture)

e <- pca.detect(data)

vec <- e$vec

axis1 <- vec[, 1]
axis2 <- vec[, 2]

anova_test(label, e)

ggplot(as.data.frame(data), aes(x = axis1, y = axis2, fill = as.factor(data[, ncol(data)]),
                                shape = as.factor(data[, ncol(data)]))) +
  geom_point(size = 2) +
  labs(x = "First Eigenvector", y = "Second Eigenvector", fill = "Population") +
  scale_fill_manual(values = c("#970973", "#FF6FF6"), labels = c("A", "B")) +
  scale_shape_manual(values = c(22, 24), guide = "none")+
  guides(fill = guide_legend(override.aes = list(shape = c(22, 24))))+
  theme_classic()

#### and then also some admixture between those 2 ####

label <- c(400, 300, 200) # will give us 20 from each population
data <- sample_generator(snp, label, mixture)

e <- pca.detect(data)

vec <- e$vec

axis1 <- vec[, 1]
axis2 <- vec[, 2]

anova_test(label, e)

ggplot(as.data.frame(data), aes(x = axis1, y = axis2, fill = as.factor(data[, ncol(data)]),
                                shape = as.factor(data[, ncol(data)]))) +
  geom_point(size = 2) +
  labs(x = "First Eigenvector", y = "Second Eigenvector", fill = "Population") +
  scale_fill_manual(values = c("#970973", "#FF6FF6", "#aa59ff"), labels = c("A", "B", "Q")) +
  scale_shape_manual(values = c(22, 24, 23), guide = "none")+
  guides(fill = guide_legend(override.aes = list(shape = c(22, 24, 23))))+
  theme_classic()



#### 3-discrete subpopulations ####

snp <- 5000
k <- 3 # number of discrete subpopulations
Fst <- c(0.01, 0.02, 0.04) # how different they are

subpop <- freq(snp, k, Fst)
subpop

parents = c(1, 2)
alpha = 19
beta = 13

mixture <- admix(snp, subpop, parents, alpha, beta)


label <- c(400, 300, 275, 0) # will give us 20 from each population
data <- sample_generator(snp, label, mixture)

e <- pca.detect(data)

library(ggplot2)

vec <- e$vec

axis1 <- vec[, 1]
axis2 <- vec[, 2]

anova_test(label, e)

# do the plotting separately

ggplot(as.data.frame(data), aes(x = axis1, y = axis2, fill = as.factor(data[, ncol(data)]),
                                shape = as.factor(data[, ncol(data)]))) +
  geom_point(size = 2) +
  labs(x = "First Eigenvector", y = "Second Eigenvector", fill = "Population") +
  scale_fill_manual(values = c("#970973", "#FF6FF6", "#FF005D"), labels = c("A", "B", "C")) +
  scale_shape_manual(values = c(22, 24, 21), guide = "none")+
  guides(fill = guide_legend(override.aes = list(shape = c(22, 24, 21))))+
  theme_classic()



#### admixture ####


label <- c(400, 300, 275, 200) # will give us 20 from each population
data <- sample_generator(snp, label, mixture)

e <- pca.detect(data)

vec <- e$vec

axis1 <- vec[, 1]
axis2 <- vec[, 2]

anova_test(label, e)

# do the plotting separately

ggplot(as.data.frame(data), aes(x = -axis1, y = axis2, fill = as.factor(data[, ncol(data)]),
                                shape = as.factor(data[, ncol(data)]))) +
  geom_point(size = 2) +
  labs(x = "First Eigenvector", y = "Second Eigenvector", fill = "Population") +
  scale_fill_manual(values = c("#970973", "#FF6FF6", "#FF005D","#aa59ff"), labels = c("A", "B", "C", "Q")) +
  scale_shape_manual(values = c(22, 24, 21, 23), guide = "none")+
  guides(fill = guide_legend(override.aes = list(shape = c(22, 24, 21, 23))))+
  theme_classic()

#### admixed but not all ####

label <- c(0, 300, 275, 200) # will give us 20 from each population
data <- sample_generator(snp, label, mixture)

e <- pca.detect(data)

vec <- e$vec

axis1 <- vec[, 1]
axis2 <- vec[, 2]

anova_test(label, e)

# do the plotting separately

ggplot(as.data.frame(data), aes(x = axis1, y = axis2, fill = as.factor(data[, ncol(data)]),
                                shape = as.factor(data[, ncol(data)]))) +
  geom_point(size = 2) +
  labs(x = "First Eigenvector", y = "Second Eigenvector", fill = "Population") +
  scale_fill_manual(values = c("#FF6FF6", "#FF005D","#aa59ff"), labels = c("B", "C", "Q")) +
  scale_shape_manual(values = c(24, 21, 23), guide = "none")+
  guides(fill = guide_legend(override.aes = list(shape = c(24, 21, 23))))+
  theme_classic()

