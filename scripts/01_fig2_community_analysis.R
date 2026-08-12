# Fig. 2 community analysis
# Data required:
#   comm    : samples x MAG relative-abundance matrix
#   metadata: data frame with row names matching comm and column 'habitat'

library(vegan)

# Shannon diversity
shannon <- diversity(comm, index = "shannon")
kruskal.test(shannon ~ metadata$habitat)
pairwise.wilcox.test(shannon, metadata$habitat, p.adjust.method = "BH")

# Hellinger transformation and Bray-Curtis dissimilarity
comm_hel <- decostand(comm, method = "hellinger")
D <- vegdist(comm_hel, method = "bray")

# PCoA
pcoa <- cmdscale(D, eig = TRUE, k = 2)

# PERMANOVA
adonis2(D ~ habitat, data = metadata, permutations = 999)

# PERMDISP: diagnostic for within-group dispersion
bd <- betadisper(D, metadata$habitat)
permutest(bd, permutations = 999)

# Because the sampling design is unbalanced, interpret PERMANOVA cautiously
# whenever dispersion differs significantly among groups.
