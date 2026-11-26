#!/usr/bin/env Rscript
# Reproduce the revision's bacterial alpha/beta diversity statistics.
library(vegan)
set.seed(42)
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
out <- args[1]
dest <- file.path(out,"beta_diversity")
meta <- read.delim(file.path(dest,"sample_metadata.tsv"))
rownames(meta) <- meta$Sample
alpha <- read.delim(file.path(out,"rarefaction","sample_rarefaction.tsv"))
groups <- c("Carbonate","Plastic","Shark_egg","ARMS")
alpha_results <- list()
for (depth in c(1000,2000)) {
    for (hits in c(1,2,5)) {
        a <- subset(alpha,Depth==depth & Minimum_hits==hits & Substrate %in% groups)
        stopifnot(nrow(a)==17)
        kw <- kruskal.test(Median_Shannon ~ Substrate,data=a)
        cp <- subset(a,Substrate %in% c("Carbonate","Plastic"))
        wt <- wilcox.test(Median_Shannon ~ Substrate,data=cp,exact=TRUE)
        alpha_results[[length(alpha_results)+1]] <- data.frame(
            Depth=depth,Minimum_hits=hits,N=nrow(a),KW_statistic=unname(kw$statistic),
            KW_p=kw$p.value,Wilcoxon_W=unname(wt$statistic),Wilcoxon_p=wt$p.value)
    }
}
write.table(do.call(rbind,alpha_results),file.path(dest,"alpha_tests.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
beta_results <- list()
axes <- list()
for (hits in c(1,2,5)) {
    m <- as.matrix(read.delim(file.path(dest,paste0("bray_mean_min",hits,".tsv")),row.names=1,check.names=FALSE))
    stopifnot(identical(rownames(m),colnames(m)),setequal(rownames(m),meta$Sample))
    ord <- wcmdscale(as.dist(sqrt(m)),k=2,eig=TRUE)
    coords <- data.frame(Sample=rownames(ord$points),PC1=ord$points[,1],PC2=ord$points[,2])
    coords$Substrate <- meta[coords$Sample,"Substrate"]
    write.table(coords,file.path(dest,paste0("PCoA_min",hits,".tsv")),sep="\t",row.names=FALSE,quote=FALSE)
    pct <- 100*ord$eig[1:2]/sum(ord$eig[ord$eig>0])
    axes[[length(axes)+1]] <- data.frame(Minimum_hits=hits,PC1_pct=pct[1],PC2_pct=pct[2])
    comparisons <- list(Replicated_substrates=groups,Carbonate_vs_plastic=c("Carbonate","Plastic"))
    for (comparison in names(comparisons)) {
        md <- meta[meta$Substrate %in% comparisons[[comparison]],,drop=FALSE]
        md$Substrate <- factor(md$Substrate)
        ds <- as.dist(m[md$Sample,md$Sample])
        pa <- adonis2(ds ~ Substrate,data=md,permutations=9999,sqrt.dist=TRUE)
        bd <- betadisper(ds,md$Substrate,type="median",bias.adjust=TRUE,sqrt.dist=TRUE)
        pd <- permutest(bd,permutations=9999)
        print(pa); print(pd)
        beta_results[[length(beta_results)+1]] <- data.frame(
            Minimum_hits=hits,Comparison=comparison,N=nrow(md),
            R2=pa[1,"R2"],Pseudo_F=pa[1,"F"],PERMANOVA_p=pa[1,"Pr(>F)"],
            PERMDISP_F=pd$tab[1,"F"],PERMDISP_p=pd$tab[1,"Pr(>F)"])
    }
}
write.table(do.call(rbind,beta_results),file.path(dest,"beta_tests.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,axes),file.path(dest,"PCoA_axis_variance.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
capture.output(sessionInfo(),file=file.path(dest,"R_sessionInfo.txt"))
