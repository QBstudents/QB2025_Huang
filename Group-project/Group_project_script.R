rm(list = ls())
getwd()
setwd("/cloud/project/QB2025_Huang/Group-project")

load("/cloud/project/QB2025_Huang/Group-project/longdataBac_objects2_datadryad.rda")
Bacteria <- longdataBac_datadryad #rename
rm(longdataBac_datadryad)

#Matrix based on diff habitat type 
bac_by_site <- with(Bacteria, tapply(Counts, list(PlotID, Sender), sum, default = 0)) 
write.table(bac_by_site, file = "bacteria_div.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

#drop some rows 

filtered_bac <- bac_by_site[!grepl("_3|_2", rownames(bac_by_site)), ]

bac <- filtered_bac
View(bac)

# Make a resemblance matrix based on Bray-Curtis
bac.db <- vegdist(bac, method = "bray", upper = TRUE, diag = TRUE)
View(bac.db)

#Heatmap
order_b <- rev(attr(bac, "Labels")) #Define Order of Sites
levelplot(as.matrix(bac.db)[, order_b], aspect = "iso", col.regions = mako, 
          xlab = "Bacteria Site", ylab = "Bacteria Site", scales = list(cex = 0.5),
          main = "Bray-Curtis Distance")

#Cluster analysis 
bac.ward <- hclust(bac.db, method = "ward.D2")

#plot with `hclust`:
par(mar = c(1, 5, 2, 2) + 0.1)
plot(bac.ward, main = "Bacteria: Ward's Clustering", 
     ylab = "Squared Bray-Curtis Distance")

#plot with `heatmap.2`: 
gplots::heatmap.2(as.matrix(bac_by_site), 
                  distfun = function(x) vegdist(x, method = "bray"),
                  hclustfun = function(x) hclust(x, method = "ward.D2"),
                  col = viridis, trace = "none", density.info = "none")

#PCoA 
bac.pcoa <- cmdscale(bac.db, eig = TRUE, k = 3) #performed a PCoA

# Variation explained by the first three axes: 
explainvar1_b <- round(bac.pcoa$eig[1]/sum(bac.pcoa$eig), 3)*100
explainvar2_b <- round(bac.pcoa$eig[2]/sum(bac.pcoa$eig), 3)*100
explainvar3_b <- round(bac.pcoa$eig[3]/sum(bac.pcoa$eig), 3)*100
sum.eig <- sum(explainvar1_b, explainvar2_b, explainvar3_b)

#identify influential species
bacREL <- bac_by_site
for(i in 1:nrow(bac_by_site)){
  bacREL[i, ] = bac_by_site[i, ]/sum(bac_by_site[i, ])
}

bac.pcoa <- add.spec.scores.class(bac.pcoa, bacREL, method = "pcoa.scores")

# Plot the PCoA ordination:
par(mar = c(5, 5, 1, 2) + 0.1)
plot(bac.pcoa$points[, 1], bac.pcoa$points[, 2], ylim = c(-0.2, 0.7),
     xlab = paste("PCoA 1 (", explainvar1_b, "%)", sep = ""),
     ylab = paste("PCoA 2 (", explainvar2_b, "%)", sep = ""),
     pch = 16, cex = 2.0, type = "n", cex.lab = 1.5, 
     cex.axis = 1.2, axes = FALSE)
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)
points(bac.pcoa$points[, 1], bac.pcoa$points[, 2],
       pch = 19, cex = 3, bg = "gray", col = "gray")
text(bac.pcoa$points[, 1], bac.pcoa$points[, 2],
     labels = row.names(bac.pcoa$points))
text(bac.pcoa$cproj[, 1], bac.pcoa$cproj[, 2],
     labels = row.names(bac.pcoa$cproj), col = "black")   #add species scores
