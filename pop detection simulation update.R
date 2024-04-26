### try the population detection simulations again...

# okay so. 

# try the scenario with 2 discrete and their admixture

snp <- 5000
K <- 2
Fst <- c(0.01, 0.02)

discrete <- c(400, 300) # how many people from the discrete pops
mix <- 1 # how many people from the admixture

parents = c(1, 2)
alpha = 19
beta = 13

subpop <- freq(snp, K, Fst)

data <- sample_generator(snp, discrete, subpop, mix, parents, alpha, beta)
data <- data[-nrow(data),]

e <- pca.detect(data)

vec <- e$vec

axis1 <- vec[, 1]
axis2 <- vec[, 2]

label <- c(discrete, mix)
anova_test(label, e)

ggplot(as.data.frame(data), aes(x = -axis1, y = axis2, fill = as.factor(data[, ncol(data)]),
                                shape = as.factor(data[, ncol(data)]))) +
  geom_point(size = 2) +
  labs(x = "First Eigenvector", y = "Second Eigenvector", fill = "Population") +
  scale_fill_manual(values = c("#970973", "#FF6FF6", "#aa59ff"), labels = c("A", "B", "Q")) +
  scale_shape_manual(values = c(22, 24, 23), guide = "none")+
  guides(fill = guide_legend(override.aes = list(shape = c(22, 24, 23))))+
  theme_classic()


## 2 discrete plot

ggplot(as.data.frame(data), aes(x = axis1, y = axis2, fill = as.factor(data[, ncol(data)]),
                                shape = as.factor(data[, ncol(data)]))) +
  geom_point(size = 2) +
  labs(x = "First Eigenvector", y = "Second Eigenvector", fill = "Population") +
  scale_fill_manual(values = c("#970973", "#FF6FF6"), labels = c("A", "B")) +
  scale_shape_manual(values = c(22, 24), guide = "none")+
  guides(fill = guide_legend(override.aes = list(shape = c(22, 24))))+
  theme_classic()


# 3 discrete plot

ggplot(as.data.frame(data), aes(x = axis1, y = axis2, fill = as.factor(data[, ncol(data)]),
                                shape = as.factor(data[, ncol(data)]))) +
  geom_point(size = 2) +
  labs(x = "First Eigenvector", y = "Second Eigenvector", fill = "Population") +
  scale_fill_manual(values = c("#970973", "#FF6FF6", "#FF005D"), labels = c("A", "B", "C")) +
  scale_shape_manual(values = c(22, 24, 21), guide = "none")+
  guides(fill = guide_legend(override.aes = list(shape = c(22, 24, 21))))+
  theme_classic()




# 3 discrete + admixed

snp <- 5000
K <- 3
Fst <- c(0.01, 0.02, 0.04)

discrete <- c(400, 300, 275) # how many people from the discrete pops
mix <- 1 # how many people from the admixture

parents = c(1, 2)
alpha = 19
beta = 13

subpop <- freq(snp, K, Fst)

data <- sample_generator(snp, discrete, subpop, mix, parents, alpha, beta)
data <- data[-nrow(data),]

e <- pca.detect(data)

vec <- e$vec

axis1 <- vec[, 1]
axis2 <- vec[, 2]

label <- c(discrete, mix)
anova_test(label, e)

ggplot(as.data.frame(data), aes(x = axis1, y = -axis2, fill = as.factor(data[, ncol(data)]),
                                shape = as.factor(data[, ncol(data)]))) +
  geom_point(size = 2) +
  labs(x = "First Eigenvector", y = "Second Eigenvector", fill = "Population") +
  scale_fill_manual(values = c("#970973", "#FF6FF6", "#FF005D","#aa59ff"), labels = c("A", "B", "C", "Q")) +
  scale_shape_manual(values = c(22, 24, 21, 23), guide = "none")+
  guides(fill = guide_legend(override.aes = list(shape = c(22, 24, 21, 23))))+
  theme_classic()


## admixed but not all

snp <- 5000
K <- 3
Fst <- c(0.01, 0.02, 0.04)

discrete <- c(0, 300, 275) # how many people from the discrete pops
mix <- 200 # how many people from the admixture

parents = c(1, 2)
alpha = 19
beta = 13

subpop <- freq(snp, K, Fst)

data <- sample_generator(snp, discrete, subpop, mix, parents, alpha, beta)
data <- data[-1, ]

e <- pca.detect(data)

vec <- e$vec

axis1 <- vec[, 1]
axis2 <- vec[, 2]

label <- c(discrete, mix)
anova_test(label, e)

ggplot(as.data.frame(data), aes(x = axis1, y = -axis2, fill = as.factor(data[, ncol(data)]),
                                shape = as.factor(data[, ncol(data)]))) +
  geom_point(size = 2) +
  labs(x = "First Eigenvector", y = "Second Eigenvector", fill = "Population") +
  scale_fill_manual(values = c("#FF6FF6", "#FF005D","#aa59ff"), labels = c("B", "C", "Q")) +
  scale_shape_manual(values = c(24, 21, 23), guide = "none")+
  guides(fill = guide_legend(override.aes = list(shape = c(24, 21, 23))))+
  theme_classic()

