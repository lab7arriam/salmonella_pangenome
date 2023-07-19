library(dplyr)
library(ggplot2)
library(ggfortify)
library(factoextra)


draw_autoplot <- function(clust_kmeans,sc_kmeans){
  p<-autoplot(clust_kmeans, 
              data=sc_kmeans, label = TRUE,frame = TRUE,
              frame.type = 'convex',label.size = 3.3,label.col='black', alpha=0.8)+
    scale_colour_manual(values = c('grey','grey','grey', 'grey','grey', 'grey', 'grey'))+
    theme_bw()+
    theme(axis.text.y = element_text(color='black', 
                                     size=14),
          axis.title.y=element_text(color="black", 
                                    size=16),
          panel.background = element_blank(), 
          axis.line = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(color='black', 
                                     size=14),
          axis.title.x = element_text(color="black", 
                                      size=16),
          legend.position = 'none'
    ) + 
    scale_x_continuous(expand = c(.1, .1)) 
  scale_y_continuous(expand = c(.1, .1))
  return(p)
}


distance_matrix <- read.table('/home/anton/Salmonella_pangenome_2023/All_jac.csv', row.names=1,header = TRUE,  sep = "\t")

sc_dist<-data.frame(t(distance_matrix))
shil_stat <- fviz_nbclust(sc_dist, kmeans, method = 'silhouette', k.max=4)
optimal_num <- which.max(shil_stat$data[, 2]) 

dist_kmeans <- kmeans(sc_dist,centers=optimal_num)
dist_kmeans$cluster
sc_dist$name=rownames(sc_dist)
sc_dist$clust= as.factor(dist_kmeans$cluster)

draw_autoplot(dist_kmeans, sc_dist)


