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
load_or_install("readr")
load_or_install("stringr")
load_or_install("tidyr")
load_or_install("locfit")
load_or_install("tibble")
load_or_install("dplyr")
load_or_install("minpack.lm")
load_or_install("Matrix")
load_or_install("ggplot2")
load_or_install("limma", from_bioconductor=TRUE)     # 从 Bioconductor 安装
load_or_install("edgeR", from_bioconductor=TRUE)     # 从 Bioconductor 安装
