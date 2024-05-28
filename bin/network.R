# Script used to create network representations of the distance matrices
setwd("//wsl$/Ubuntu/home/cisnerazerx/Documents/paper3/classification/VirTaK")

c1 <- read.csv("CLUSTER22_arfiviricetes/matrix_network.csv", check.names = F, row.names = 1)
C1 <- 1-c1
C1_mat <- as.matrix(C1)

library(RColorBrewer)
metadata <- read.delim("CLUSTER22_arfiviricetes/metadata.txt")
coul <- brewer.pal(nlevels(as.factor(metadata$Source)), "Set2")
my_color <- coul[as.numeric(as.factor(metadata$Source))]

library(igraph)
net <- graph_from_adjacency_matrix(C1_mat, mode='undirected', weighted = T, diag = F)
par(mar=c(0,0,0,0))
set.seed(4)
plot(net, layout=layout.fruchterman.reingold, main="", vertex.size=6, vertex.label.cex=0.5, vertex.label.color="black", vertex.frame.color="black", vertex.color=my_color)

library(reshape2)
C1_mat_melt <- melt(C1_mat)
C1_mat_melt <- subset(C1_mat_melt, Var1 != Var2)
Q1 <- unname(quantile(C1_mat_melt$value, prob=c(.25))) # 75%
MED <- median(C1_mat_melt$value) # 50%
Q3 <- unname(quantile(C1_mat_melt$value, prob=c(.75))) # 25%

C1_mat[C1_mat<Q1] <- 0
net75 <- graph_from_adjacency_matrix(C1_mat, mode='undirected', weighted = T, diag = F)
par(mar=c(0,0,0,0))
set.seed(4)
plot(net75, layout=layout.fruchterman.reingold, main="", vertex.size=6, vertex.label.cex=0.5, vertex.label.color="black", vertex.frame.color="black", vertex.color=my_color)

C1_mat[C1_mat<MED] <- 0
net50 <- graph_from_adjacency_matrix(C1_mat, mode='undirected', weighted = T, diag = F)
par(mar=c(0,0,0,0))
set.seed(4)
plot(net50, layout=layout.fruchterman.reingold, main="", vertex.size=6, vertex.label.cex=0.5, vertex.label.color="black", vertex.frame.color="black", vertex.color=my_color)

C1_mat[C1_mat<Q3] <- 0
net25 <- graph_from_adjacency_matrix(C1_mat, mode='undirected', weighted = T, diag = F)
par(mar=c(0,0,0,0))
set.seed(4)
plot(net25, layout=layout.fruchterman.reingold, main="", vertex.size=6, vertex.label.cex=0.5, vertex.label.color="black", vertex.frame.color="black", vertex.color=my_color)

C1_mat[C1_mat<0.5] <- 0
net05 <- graph_from_adjacency_matrix(C1_mat, mode='undirected', weighted = T, diag = F)
par(mar=c(0,0,0,0))
set.seed(4)
plot(net05, layout=layout.fruchterman.reingold, main="", vertex.size=6, vertex.label.cex=0.5, vertex.label.color="black", vertex.frame.color="black", vertex.color=my_color)

####################################################################
c1 <- read.csv("clusters.csv", check.names = F, row.names = 1)
C1 <- 1-c1
C1_mat <- as.matrix(C1)
C1_mat[C1_mat<0.055] <- 0
net <- graph_from_adjacency_matrix(C1_mat, mode='undirected', weighted = T, diag = F)
par(mar=c(0,0,0,0))
set.seed(4)
plot(net, layout=layout.fruchterman.reingold, main="", vertex.size=6, vertex.label.cex=0.5, vertex.label.color="black", vertex.frame.color="black")
