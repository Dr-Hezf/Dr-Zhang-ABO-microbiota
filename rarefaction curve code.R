rm(list = ls())
getwd()
setwd("C:/Users/hejia/Desktop")

install.packages("doBy")
install.packages("vegan")
install.packages("ggplot2")
install.packages("ggalt")
install.packages("ragg")
library(vegan)   #compute the shannon、simpson、chao1、ACE
library(ggplot2)  
library(doBy) 
library(ggalt)  
library(ragg)


alpha_index <- function(x, method = 'richness', tree = NULL, base = exp(1)) {
  if (method == 'richness') result <- rowSums(x > 0)    
  else if (method == 'chao1') result <- estimateR(x)[2, ]    
  else if (method == 'ace') result <- estimateR(x)[4, ]   
  else if (method == 'shannon') result <- diversity(x, index = 'shannon', base = base)    
  else if (method == 'simpson') result <- diversity(x, index = 'simpson')    #Gini-Simpson 
  else if (method == 'pielou') result <- diversity(x, index = 'shannon', base = base) / log(estimateR(x)[1, ], base)    #Pielou 
  else if (method == 'gc') result <- 1 - rowSums(x == 1) / rowSums(x)    #goods_coverage
  else if (method == 'pd' & !is.null(tree)) {    #PD_whole_tree
    pd <- pd(x, tree, include.root = FALSE)
    result <- pd[ ,1]
    names(result) <- rownames(pd)
  }
  result
}

#based on step，to compute Alpha diversity
alpha_curves <- function(x, step, method = 'richness', rare = NULL, tree = NULL, base = exp(1)) {
  x_nrow <- nrow(x)
  if (is.null(rare)) rare <- rowSums(x) else rare <- rep(rare, x_nrow)
  alpha_rare <- list()
  
  for (i in 1:x_nrow) {
    step_num <- seq(0, rare[i], step)
    if (max(step_num) < rare[i]) step_num <- c(step_num, rare[i])
    
    alpha_rare_i <- NULL
    for (step_num_n in step_num) alpha_rare_i <- c(alpha_rare_i, alpha_index(x = rrarefy(x[i, ], step_num_n), method = method, tree = tree, base = base))
    names(alpha_rare_i) <- step_num
    alpha_rare <- c(alpha_rare, list(alpha_rare_i))
  }
  
  names(alpha_rare) <- rownames(x)
  alpha_rare
}



otu<-read.csv("张-otu-质控图.csv",header = T,row.names = 1,stringsAsFactors = F,check.names = F)   #OTU data


##chao1
chao1_curves <- alpha_curves(otu,step = 2000,method = 'chao1') 

# ggplot2
plot_chao1 <- data.frame()
for (i in names(chao1_curves)) {
  chao1_curves_i <- (chao1_curves[[i]])
  chao1_curves_i <- data.frame(rare = names(chao1_curves_i), alpha = chao1_curves_i, sample = i, stringsAsFactors = FALSE)
  plot_chao1 <- rbind(plot_chao1, chao1_curves_i)
}

rownames(plot_chao1) <- NULL
plot_chao1$rare <- as.numeric(plot_chao1$rare)
plot_chao1$alpha <- as.numeric(plot_chao1$alpha)

#ggplot2 
yt_chao1 <- ggplot(plot_chao1, aes(rare, alpha, color = sample)) +
  geom_line() +
  labs(x = 'Number of sequences', y = 'chao1', color = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill = 'transparent')) +
  theme(axis.title = element_text(size = 150),axis.text = element_text(size=100,angle=0)) +            
  geom_vline(xintercept = min(rowSums(otu)), linetype = 2) +
  scale_x_continuous(breaks = seq(0, 40000, 5000), labels = as.character(seq(0, 40000, 5000)))

yt_chao1  

ggsave("yt_chao1-4.pdf",yt_chao1,width = 200,dpi = 30,height = 200,units = "cm",limitsize = FALSE)    




##richnes
richness_curves <- alpha_curves(otu,step = 2000,method = 'richness') 

plot_richness <- data.frame()
for (i in names(richness_curves)) {
  richness_curves_i <- (richness_curves[[i]])
  richness_curves_i <- data.frame(rare = names(richness_curves_i), alpha = richness_curves_i, sample = i, stringsAsFactors = FALSE)
  plot_richness <- rbind(plot_richness, richness_curves_i)
}

rownames(plot_richness) <- NULL
plot_richness$rare <- as.numeric(plot_richness$rare)
plot_richness$alpha <- as.numeric(plot_richness$alpha)

#ggplot2 
yt_richness <- ggplot(plot_richness, aes(rare, alpha, color = sample)) +
  geom_line() +
  labs(x = 'Number of sequences', y = 'richness', color = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill = 'transparent')) +
  theme(axis.title = element_text(size = 150),axis.text = element_text(size=100,angle=0)) +            
  geom_vline(xintercept = min(rowSums(otu)), linetype = 2) +
  scale_x_continuous(breaks = seq(0, 40000, 5000), labels = as.character(seq(0, 40000, 5000)))

yt_richness  

ggsave("yt_richness.pdf",yt_richness,width = 200,dpi = 30,height = 200,units = "cm",limitsize = FALSE)    



##shannon

shannon_curves <- alpha_curves(otu,step = 2000,method = 'shannon') 


plot_shannon <- data.frame()
for (i in names(shannon_curves)) {
  shannon_curves_i <- (shannon_curves[[i]])
  shannon_curves_i <- data.frame(rare = names(shannon_curves_i), alpha = shannon_curves_i, sample = i, stringsAsFactors = FALSE)
  plot_shannon <- rbind(plot_shannon, shannon_curves_i)
}

rownames(plot_shannon) <- NULL
plot_shannon$rare <- as.numeric(plot_shannon$rare)
plot_shannon$alpha <- as.numeric(plot_shannon$alpha)

#ggplot2 
yt_shannon <- ggplot(plot_shannon, aes(rare, alpha, color = sample)) +
  geom_line() +
  labs(x = 'Number of sequences', y = 'shannon', color = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill = 'transparent')) +
  theme(axis.title = element_text(size = 150),axis.text = element_text(size=100,angle=0)) +             
  geom_vline(xintercept = min(rowSums(otu)), linetype = 2) +
  scale_x_continuous(breaks = seq(0, 40000, 5000), labels = as.character(seq(0, 40000, 5000)))

yt_shannon  

ggsave("yt_shannon.pdf",yt_shannon,width = 200,dpi = 30,height = 200,units = "cm",limitsize = FALSE)    
