#!/usr/bin/env Rscript

# 检查并安装缺失的包
load_or_install <- function(package, from_bioconductor=FALSE) {
  if (!require(package, character.only = TRUE)) {
    if (from_bioconductor) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(package)
    } else {
      install.packages(package, dependencies = TRUE)
    }
    library(package, character.only = TRUE)
  }
}

# 加载所需的包
load_or_install("argparser")
load_or_install("ggplot2", from_bioconductor=TRUE)
load_or_install("tibble")
load_or_install("dplyr")
load_or_install("minpack.lm")


p <- arg_parser("core and pan gene analysis from gene/family PAV table.",
                hide.opts = T )

# Add command line arguments
p <- add_argument(p, "pav", help="input: genePAV.table", type="character" )
p <- add_argument(p, "sim", help="Number of simulations,  0 for all combinations", type="numeric" ,default = 0 )
p <- add_argument(p, "outpre", help="output prefix", type="character", default = "Result" )

argv <- parse_args(p)

input <- argv$pav
sim <- argv$sim
pre <- argv$outpre

## 数据格式
#ID  ind1 ind2 ind3 ind4
#gene1 0 1 0 0
#gene2 0 0 3 0
#gene3 1 0 5 1
#gene4 1 7 1 1

pav <- read.delim(input, header = T, row.names = 1)
pav <- pav[rowSums(pav)>0,] ## 去除0行 
samples <- colnames(pav)

## 0 1 矩阵
pav[pav > 0] = 1

data <- tibble(samples =  character(), 
               sampleNum =numeric(),
               Pan =numeric(),
               Core = numeric(),
               )

if (sim == 0) {
  for (i in 1:ncol(pav)) {
    ## 列出所有组合
    comb <- combn(samples,i)
    # print(comb)
    maxcomb <- 1000
    if (ncol(comb) > maxcomb ) {
      comb <- comb[,1:maxcomb]
    }
    for ( j in 1:ncol(comb)) {
      subPAV <- pav[, comb[,j]]
      if (i == 1) {
        Npan <- sum(subPAV)
        Ncore <- Npan
      } else {
        subPAV <- subPAV[rowSums(subPAV)>0,]
        Npan <- nrow(subPAV)
        Ncore <- nrow(subPAV[rowSums(subPAV) == i,])
      }
      name <- paste(comb[,j], sep="", collapse = "/")
    
      #print(paste(name,":", i, Npan, Ncore, Nvar))
      data <- add_row(data, 
                     samples = name,
                     sampleNum = i,
                     Pan = Npan,
                     Core = Ncore,
                     )
    }
  }
} 

if (sim > 0) {
  for (i in 1:sim) {
    sim_samples <- sample(samples, replace = F) #不重复抽样
    sim_pav <- pav[,sim_samples]
    for (j in 1:ncol(pav)) {
      subsim <- sim_pav[,1:j]
      
      if (j == 1) {
        Npan <- sum(subsim)
        Ncore <- Npan
      }else{
        sum = rowSums(subsim)
        Ncore <- nrow(subsim[sum == j,])
        Npan <- nrow(subsim[sum > 0,])
      }
      data <- add_row(data, 
                      samples = as.character(i),
                      sampleNum = j,
                      Pan = Npan,
                      Core = Ncore)
    }
  }
}


## 汇总
dsum <- dplyr::group_by(data, sampleNum) %>%
  dplyr::summarise(Pan.mean = mean(Pan),
            Pan.sd = sd(Pan),
            Core.mean = mean(Core),
            Core.sd = sd(Core))
dsum <- round(dsum, 2)

#write.table(data, file=paste(pre, "CorePan.txt", sep="."),sep = "\t",quote = F, row.names = F)
#write.table(dsum, file=paste(pre, "CorePan_summary.txt", sep="."),sep = "\t",quote = F, row.names = F)




## nls拟合
# core: 指数回归, y=Ae^(Bx)+C
# start估计
c.0 <- min(dsum$Core.mean) * 0.5
model.0 <- lm(log(Core.mean - c.0) ~ sampleNum, data=dsum)
start.0 <- list(a=as.numeric(exp(coef(model.0)[1])), b=as.numeric(coef(model.0)[2]), c=c.0)

# 拟合
model.c <- nlsLM(Core.mean ~ a * exp(b * sampleNum) + c, 
                 data = dsum, 
                 start = start.0)

#sink(file = paste(pre, "CorePan_core_model.txt", sep="."),append = FALSE)
print(summary(model.c), digits = 10)
#sink()

# pan: 幂律回归, y=Ax^B+C = A * e^(B*log(x))+C,  0<B<1 means open pan-genome 
# start估计

# 假设C < Pan.mean
# c.1 <- min(dsum$Pan.mean)*0.9
# model.1 <- lm(log(Pan.mean - c.1) ~ log(sampleNum), data=dsum)
# start.1 <- list(a = as.numeric(exp(coef(model.1)[1])), 
#                b = as.numeric(coef(model.1)[2]), 
#                c = c.1)


# 假设C > Pan.mean
c.1 <- min(dsum$Pan.mean)*1.2
model.1 <- lm(log(c.1 - Pan.mean) ~ log(sampleNum), data=dsum)
start.1 <- list(a = -as.numeric(exp(coef(model.1)[1])), 
                b = as.numeric(coef(model.1)[2]), 
                c = c.1)

model.p <- nlsLM(Pan.mean ~ a * (sampleNum^b) +c ,
                 data = dsum,
                 start = start.1)

#sink(file = paste(pre, "CorePan_pan_model.txt", sep="."),append = FALSE)
print(summary(model.p), digits = 10)
#sink()

#summary(model.p)
#coef(summary(model.p))

## 颜色设置
colors <- c("darkgreen","red")
## 点图+拟合曲线
p1 <- ggplot() +
  geom_point(aes(x=data$sampleNum, y=data$Core, color = "Core"), alpha=0.2) +  # core 
  geom_point(aes(x=data$sampleNum, y=data$Pan,  color = "Pan") , alpha=0.2) +  # pan  
  geom_point(aes(x = dsum$sampleNum, y = dsum$Core.mean),color = "black" ) +  # core mean
  geom_point(aes(x = dsum$sampleNum, y = dsum$Pan.mean),color = "black" ) +  # pan mean
  # core拟合曲线
  geom_smooth(method = "nlsLM", 
              method.args = list( start = start.0, control = nls.lm.control(maxiter = 1000)),
              formula =  y ~ a * exp(b * x) + c ,
              se = FALSE,
              data = dsum,
              aes(x = sampleNum, y= Core.mean),
              color = "red",
              ) +
  # pan拟合曲线
  geom_smooth(method = "nlsLM", 
              method.args = list( start = start.1,control = nls.lm.control(maxiter = 1000)),
              formula =  y ~ a * x ^ b + c,
              se = FALSE,
              data = dsum,
              aes(x = sampleNum, y= Pan.mean),
              color = "darkgreen",
  ) +
  scale_x_continuous(limits = c(0,NA)) +  # 设置x显示范围
  ylab("Number of TE family") +
  xlab("Number of Genomes") +
  # 设置图例
  scale_color_manual(values=colors, 
                     breaks=c("Pan", "Core"),
                     name = NULL
  ) + 
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.position = c(0.9, 0.6))

ggsave(paste(pre, "CorePan_fitmodel.pdf", sep = "."), plot = p1 , device = "pdf", width = 7,height = 7 )




p2 <- ggplot() +
  geom_point(aes(x=data$sampleNum, y=data$Core, color = "Core"), alpha=0.2) +  # core 
  geom_point(aes(x=data$sampleNum, y=data$Pan,  color = "Pan") , alpha=0.2) +  # pan  
  geom_point(aes(x = dsum$sampleNum, y = dsum$Core.mean),color = "black" ) +  # core mean
  geom_point(aes(x = dsum$sampleNum, y = dsum$Pan.mean),color = "black" ) +  # pan mean
  # core拟合曲线
  geom_smooth(
              se = FALSE,
              data = dsum,
              aes(x = sampleNum, y= Core.mean),
              color = "red",
  ) +
  # pan拟合曲线
  geom_smooth(
              se = FALSE,
              data = dsum,
              aes(x = sampleNum, y= Pan.mean),
              color = "darkgreen",
  ) +
  scale_x_continuous(limits = c(0,NA)) +  # 设置x显示范围
  ylab("Number of TE family") +
  xlab("Number of Genomes") +
  # 设置图例
  scale_color_manual(values=colors, 
                     breaks=c("Pan", "Core"),
                     name = NULL
  ) + 
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.position = c(0.9, 0.6))
ggsave(paste(pre, "CorePan_fitsmooth.pdf", sep = "."), plot = p2 , device = "pdf", width = 7,height = 7 )

