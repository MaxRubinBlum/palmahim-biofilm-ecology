#!/usr/bin/env Rscript
# PERMANOVA, Mantel and Procrustes tests for Supplementary Table 5.
# Data required:
#   abundance.csv            samples x MAGs
#   guild_abundance.csv      samples x strict guilds OR construct upstream
#   primary_abundance.csv    samples x primary-production traits
#   metadata.csv             columns sample, habitat
#
# The singleton crab sample should be excluded from substrate PERMANOVA;
# use only replicated substrate classes, as stated in the Supplementary Methods.

suppressPackageStartupMessages(library(vegan))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) stop("Usage: 02b...R abundance.csv guild.csv primary.csv metadata.csv outdir")
read_mat <- function(f) as.matrix(read.csv(f, row.names=1, check.names=FALSE))
A <- read_mat(args[1]); G <- read_mat(args[2]); P <- read_mat(args[3])
meta <- read.csv(args[4], stringsAsFactors=FALSE); rownames(meta) <- meta$sample
outdir <- args[5]; dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
common <- Reduce(intersect, list(rownames(A), rownames(G), rownames(P), rownames(meta)))
A <- A[common,,drop=FALSE]; G <- G[common,,drop=FALSE]; P <- P[common,,drop=FALSE]; meta <- meta[common,,drop=FALSE]

hbc <- function(X) vegdist(decostand(X, method="hellinger"), method="bray")
Dtax <- hbc(A); Dguild <- hbc(G); Dprimary <- hbc(P)

# Replicated substrate classes only.
tab <- table(meta$habitat)
keep_h <- names(tab)[tab >= 2]
keep <- meta$habitat %in% keep_h
m2 <- droplevels(transform(meta[keep,,drop=FALSE], habitat=factor(habitat)))

do_permanova <- function(X, label) {
  d <- hbc(X[keep,,drop=FALSE])
  fit <- adonis2(d ~ habitat, data=m2, permutations=9999)
  data.frame(representation=label, F=fit$F[1], R2=fit$R2[1], P=fit$`Pr(>F)`[1])
}
perm <- rbind(do_permanova(A,"taxonomic"), do_permanova(G,"guild"), do_permanova(P,"primary"))
write.csv(perm, file.path(outdir,"permanova_replicated_substrates.csv"), row.names=FALSE)

mantel_one <- function(d2, label) {
  z <- mantel(Dtax, d2, permutations=9999, method="pearson")
  data.frame(comparison=paste("taxonomic vs", label), r=unname(z$statistic), P=z$signif)
}
write.csv(rbind(mantel_one(Dguild,"guild"), mantel_one(Dprimary,"primary")),
          file.path(outdir,"mantel.csv"), row.names=FALSE)

# PCoA coordinates used for symmetric Procrustes / PROTEST.
ord <- function(d) cmdscale(d, eig=TRUE, k=min(5, attr(d,"Size")-1))$points
Ot <- ord(Dtax); Og <- ord(Dguild); Op <- ord(Dprimary)
pro_one <- function(o2, label) {
  z <- protest(Ot, o2, permutations=9999, symmetric=TRUE)
  data.frame(comparison=paste("taxonomic vs", label), correlation=sqrt(1-z$ss), P=z$signif)
}
write.csv(rbind(pro_one(Og,"guild"), pro_one(Op,"primary")),
          file.path(outdir,"procrustes.csv"), row.names=FALSE)
