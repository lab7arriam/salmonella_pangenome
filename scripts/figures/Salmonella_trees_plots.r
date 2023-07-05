library(ggtree)
library(treeio)
library(ggplot2)
library(ape)
library(ggnewscale)

library(rcartocolor)
library(RColorBrewer)
library(dplyr)

# Create color palette
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'seq',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]


#Read files
beast_tree <- read.tree("./core_snp_big.fasta.raxml.bestTree")
genotype <- read.table("./sm.csv", sep=",", stringsAsFactor=F)



#beast_tree$edge.length <- (log(beast_tree$edge.length)+(min(log(beast_tree$edge.length))*-1))e
# or this gives the same
beast_tree$edge.length <- log(beast_tree$edge.length*(1/min(beast_tree$edge.length)))

ggtree(beast_tree)
# Make index as names of assemblies, it is necessary for using gheatmap
rownames(genotype) <- genotype$V1
genotype <- genotype[2:4]
colnames(genotype)[1] ="Host"
colnames(genotype)[2] ="Serovar"
colnames(genotype)[3] ="Source"

# Create names for annotation (the values that you can see on a legend)
names_serovar <- (c(genotype$Serovar))
names_set_serovar <- names_serovar[!duplicated(names_serovar)]
colors_serovar <-  sample(color, length(names_set_serovar) + 1)

names_hosts <- (c(genotype$Host))
names_set_hosts <- names_hosts[!duplicated(names_hosts)]
colors_hosts <-  sample(color, length(names_set_hosts) + 1)

names_source <- (c(genotype$Source))
names_set_source <- names_source[!duplicated(names_source)]
colors_source <-  sample(color, length(names_set_source) + 1)


colors_manua <- c('#f3f3f3',
                  '#f3f3f3',
                  '#8dd3c7',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#bebada',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#fdb462',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#57b65f',
                  '#80b1d3',
                  '#e6ff9a',
                  '#ffd0c9',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#e5e5e5',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#b3de69',
                  '#7fc97f',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#fb8072',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#bc80bd',
                  '#fccde5',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3',
                  '#f3f3f3')

length(names_set_serovar)


#Building a tree, branch.length='none' to align leaves plot.phylo()  branch.length='none',
p <- ggtree(beast_tree, ladderize=FALSE)


h1 <- gheatmap(p, genotype[1], offset = 30, width = 0.091,
               colnames_position = "top", colnames_offset_y = 3) + scale_fill_manual(breaks=names_set_hosts,
                                                          values=c('#e5e5e5',
                                                                   '#1f78b4',
                                                                   '#b2df8a',
                                                                   '#0F8887',
                                                                   '#ffffb3',
                                                                   '#fbb4ae',
                                                                   '#fdbf6f',
                                                                   '#ff7f00',
                                                                   '#cab2d6'), name=c("Host")) + coord_cartesian(clip = "off")


genotype[2]

h2 <- h1 + new_scale_fill()

h3 <- gheatmap(h2, genotype[2], offset = 5, width = 0.091, colnames_position = "top", colnames_offset_y = 3) +
  scale_fill_manual(breaks=names_set_serovar,
                     values=colors_manua, name=c("Serovar")) + coord_cartesian(clip = "off") +
  vexpand(.009)
