rm(list = ls())
getwd()
setwd("C:/Users/hejia/Desktop")
#安装、加载必须的R包
install.packages("doBy")
install.packages("vegan")
install.packages("ggplot2")
install.packages("ggalt")
install.packages("ragg")
library(vegan)   #用于计算shannon、simpson、chao1、ACE指数等
library(ggplot2)  #用于作图
library(doBy)  #用于分组统计
library(ggalt)  #用于绘制拟合曲线
library(ragg)

##定义函数
#计算多种Alpha多样性指数、结果返回至向量
alpha_index <- function(x, method = 'richness', tree = NULL, base = exp(1)) {
  if (method == 'richness') result <- rowSums(x > 0)    #丰富度指数
  else if (method == 'chao1') result <- estimateR(x)[2, ]    #Chao1 指数
  else if (method == 'ace') result <- estimateR(x)[4, ]    #ACE 指数
  else if (method == 'shannon') result <- diversity(x, index = 'shannon', base = base)    #Shannon 指数
  else if (method == 'simpson') result <- diversity(x, index = 'simpson')    #Gini-Simpson 指数
  else if (method == 'pielou') result <- diversity(x, index = 'shannon', base = base) / log(estimateR(x)[1, ], base)    #Pielou 均匀度
  else if (method == 'gc') result <- 1 - rowSums(x == 1) / rowSums(x)    #goods_coverage
  else if (method == 'pd' & !is.null(tree)) {    #PD_whole_tree
    pd <- pd(x, tree, include.root = FALSE)
    result <- pd[ ,1]
    names(result) <- rownames(pd)
  }
  result
}

#根据抽样步长（step），统计每个稀释梯度下的 Alpha 多样性指数，结果返回至列表
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


#导入数据
otu<-read.csv("张-otu-质控图.csv",header = T,row.names = 1,stringsAsFactors = F,check.names = F)   #读取OTU原始文件


##chao1指数
#以2000的步长计算chao1
chao1_curves <- alpha_curves(otu,step = 2000,method = 'chao1') #步长根据需要调整

#获得 ggplot2 作图文件
plot_chao1 <- data.frame()
for (i in names(chao1_curves)) {
  chao1_curves_i <- (chao1_curves[[i]])
  chao1_curves_i <- data.frame(rare = names(chao1_curves_i), alpha = chao1_curves_i, sample = i, stringsAsFactors = FALSE)
  plot_chao1 <- rbind(plot_chao1, chao1_curves_i)
}

rownames(plot_chao1) <- NULL
plot_chao1$rare <- as.numeric(plot_chao1$rare)
plot_chao1$alpha <- as.numeric(plot_chao1$alpha)

#ggplot2 作图
yt_chao1 <- ggplot(plot_chao1, aes(rare, alpha, color = sample)) +
  geom_line() +
  labs(x = 'Number of sequences', y = 'chao1', color = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill = 'transparent')) +
  theme(axis.title = element_text(size = 150),axis.text = element_text(size=100,angle=0)) +              #axis.text调节刻度标签大小；axis.title调整轴标题字体
  geom_vline(xintercept = min(rowSums(otu)), linetype = 2) +
  scale_x_continuous(breaks = seq(0, 40000, 5000), labels = as.character(seq(0, 40000, 5000)))

yt_chao1   #这一步是为了了解电脑已经完成任务

ggsave("yt_chao1-4.pdf",yt_chao1,width = 200,dpi = 30,height = 200,units = "cm",limitsize = FALSE)     #当样本量很大，单位需要改成cm，大小要足够大，才能够显示全图




##richness指数
#以2000的步长计算richness指数
richness_curves <- alpha_curves(otu,step = 2000,method = 'richness') #步长根据需要调整

#获得 ggplot2 作图文件
plot_richness <- data.frame()
for (i in names(richness_curves)) {
  richness_curves_i <- (richness_curves[[i]])
  richness_curves_i <- data.frame(rare = names(richness_curves_i), alpha = richness_curves_i, sample = i, stringsAsFactors = FALSE)
  plot_richness <- rbind(plot_richness, richness_curves_i)
}

rownames(plot_richness) <- NULL
plot_richness$rare <- as.numeric(plot_richness$rare)
plot_richness$alpha <- as.numeric(plot_richness$alpha)

#ggplot2 作图
yt_richness <- ggplot(plot_richness, aes(rare, alpha, color = sample)) +
  geom_line() +
  labs(x = 'Number of sequences', y = 'richness', color = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill = 'transparent')) +
  theme(axis.title = element_text(size = 150),axis.text = element_text(size=100,angle=0)) +              #axis.text调节刻度标签大小；axis.title调整轴标题字体
  geom_vline(xintercept = min(rowSums(otu)), linetype = 2) +
  scale_x_continuous(breaks = seq(0, 40000, 5000), labels = as.character(seq(0, 40000, 5000)))

yt_richness   #这一步是为了了解电脑已经完成任务

ggsave("yt_richness.pdf",yt_richness,width = 200,dpi = 30,height = 200,units = "cm",limitsize = FALSE)     #当样本量很大，单位需要改成cm，大小要足够大，才能够显示全图



##shannon
#以2000的步长计算shannon
shannon_curves <- alpha_curves(otu,step = 2000,method = 'shannon') #步长根据需要调整

#获得 ggplot2 作图文件
plot_shannon <- data.frame()
for (i in names(shannon_curves)) {
  shannon_curves_i <- (shannon_curves[[i]])
  shannon_curves_i <- data.frame(rare = names(shannon_curves_i), alpha = shannon_curves_i, sample = i, stringsAsFactors = FALSE)
  plot_shannon <- rbind(plot_shannon, shannon_curves_i)
}

rownames(plot_shannon) <- NULL
plot_shannon$rare <- as.numeric(plot_shannon$rare)
plot_shannon$alpha <- as.numeric(plot_shannon$alpha)

#ggplot2 作图
yt_shannon <- ggplot(plot_shannon, aes(rare, alpha, color = sample)) +
  geom_line() +
  labs(x = 'Number of sequences', y = 'shannon', color = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill = 'transparent')) +
  theme(axis.title = element_text(size = 150),axis.text = element_text(size=100,angle=0)) +              #axis.text调节刻度标签大小；axis.title调整轴标题字体
  geom_vline(xintercept = min(rowSums(otu)), linetype = 2) +
  scale_x_continuous(breaks = seq(0, 40000, 5000), labels = as.character(seq(0, 40000, 5000)))

yt_shannon   #这一步是为了了解电脑已经完成任务

ggsave("yt_shannon.pdf",yt_shannon,width = 200,dpi = 30,height = 200,units = "cm",limitsize = FALSE)     #当样本量很大，单位需要改成cm，大小要足够大，才能够显示全图

