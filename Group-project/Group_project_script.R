rm(list = ls())
getwd()
setwd("/cloud/project/QB2025_Huang/Group-project")

# Load data ----
load("/cloud/project/QB2025_Huang/Group-project/longdataBac_objects2_datadryad.rda")
Bacteria <- longdataBac_datadryad #rename
rm(longdataBac_datadryad)
load("Bac_wide_plot_final2_datadryad.rda")
Bac_env <- Bac_wide_plot #rename
rm(Bac_wide_plot)
bac <- read.table("/cloud/project/QB2025_Huang/Group-project/Cleaned_data/bacteria_div.txt", 
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

## Data cleaning ----
#Matrix based on diff habitat type 
#write.table(bac_by_site, file = "bacteria_div.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

#drop some rows 
bac.reduced <- bac[!grepl("_3|_2|_1", rownames(bac)), ]
env.reduced <- Bac_env[!grepl("_3|_2|_1", Bac_env$PlotID), ]
rownames(env.reduced) <- env.reduced$PlotID
env.reduced <- env.reduced[, -c(1, 2, 3, 9)]
#write.table(env.reduced, file = "bac_env.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
Bacteria.reduced <- Bacteria[!grepl("_3|_2|_1", Bacteria$PlotID), ]
xy <- aggregate(cbind(POINT_X, POINT_Y) ~ PlotID, data = Bacteria.reduced, FUN = mean)
#write.table(xy, file = "bac_xy.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)


## Distance matrix
# Bray-Curtis: 
bac.bc <- vegdist(bac.reduced, method = "bray", upper = TRUE, diag = TRUE)
# Jaccard 
bac.jac <- vegdist(bac.reduced, method = "jaccard", upper = TRUE, diag = TRUE)

### Explore a bit ----
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


## Try hypothesis testing ----

### PERMANOVA ----
# None of these are significant, no matter using jaccard or bray
land_type <- env.reduced$Landscape
habitat <- env.reduced$Habitat
#adonis2(bac.reduced ~ land_type, method = "jaccard", permutation = 999)
adonis2(bac.reduced ~ land_type, method = "bray", permutation = 999)
#adonis2(bac.reduced ~ habitat, method = "jaccard", permutation = 999)
adonis2(bac.reduced ~ habitat, method = "bray", permutation = 999)

### Mantel test ----
env.bc <- vegdist(scale(env.reduced[3:5]), method = "euclid") # env matrix
mantel(bac.bc, env.bc)

### Constrained Ordination ----
env <- as.matrix(env.reduced[3:5]) #Continous env conditions
bac.dbrda <- dbrda(bac.reduced ~ ., as.data.frame(env)) # using abundance based distance
bac.dbrda_j <- dbrda(bac.jac ~ ., as.data.frame(env)) # using incidence based distance
ordiplot(bac.dbrda)
#ordiplot(bac.dbrda_j) #This does not have ANOVA result
bc.explainvar1 <-  round(bac.dbrda$CCA$eig[1] /
                           sum(c(bac.dbrda$CCA$eig, bac.dbrda$CA$eig)), 3) * 100
bc.explainvar2 <- round(bac.dbrda$CCA$eig[2] /
                          sum(c(bac.dbrda$CCA$eig, bac.dbrda$CA$eig)), 3) * 100
# Plot the ordination plot:
par(mar = c(5,5,4,4) + 0.1)
plot(scores(bac.dbrda, display = "wa"), xlim = c(-1.3, 1.1), ylim = c(-1.1, 2.7), 
     xlab = paste("dbRDA 1 (", bc.explainvar1, "%)",sep = ""), 
     ylab = paste("dbRDA 2 (", bc.explainvar2, "%)", sep = ""),
     pch = 16, cex = 2.0, type = "n", cex.lab = 1.5,
     cex.axis = 1.2, axes = FALSE)
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)
points(scores(bac.dbrda, display = "wa"), pch = 19, cex = 1, bg = "gray", 
       col = "gray")
text(scores(bac.dbrda, display = "wa"), 
     labels = substr(row.names(scores(bac.dbrda, display = "wa")), 1, 4))
bc.vectors <- scores(bac.dbrda, display = "bp")
arrows(0, 0, bc.vectors[, 1], bc.vectors[, 2],
       lwd = 2, lty = 1, length = 0.2, col = "blue")
text(bc.vectors[, 1], bc.vectors[, 2], pos = 3, 
     labels = row.names(bc.vectors), col = "blue")
axis(side = 3, lwd.ticks = 2, cex.axis = 1.2, las = 1, col = "red", 
     lwd = 2.2, at = pretty(range(bc.vectors[, 1])) * 2, 
     labels = pretty(range(bc.vectors[, 1])))
axis(side = 4, lwd.ticks = 2, cex.axis = 1.2, las = 1, col = "red", 
     lwd = 2.2, at = pretty(range(bc.vectors[, 2])) * 2, 
     labels = pretty(range(bc.vectors[, 2])))

### Variation Partitioning ----
bac.dbrda$call
bac.env.mod <- model.matrix( ~ MAP + MAT + average_Temp_DL, data = as.data.frame(env))[,-1]
bac.rs <- rowSums(bac.reduced)/sum(bac.reduced) # relative abundance of each site
bac.pcnmw <- pcnm(dist(xy), w = bac.rs, dist.ret = T)
bac.pcnmw$values > 0 # filter for the non-negatice eigenvalues as these are the only that are meaningful

# Model selection that drops the redundant eigenvalues: 
bac.space <- as.data.frame(scores(bac.pcnmw))
bac.pcnm.mod0 <- dbrda(bac.bc ~ 1, bac.space)
bac.pcnm.mod1 <- dbrda(bac.bc ~ ., bac.space)
step.bac.pcnm <- ordiR2step(bac.pcnm.mod0, bac.pcnm.mod1, perm.max = 200)

plot(step.bac.pcnm) # Visualize the biplot for spatial factors' influence on sites composition
step.bac.pcnm$anova
space.bac.mod <- model.matrix(~ PCNM2 + PCNM3 + PCNM5 + PCNM1 + 
                            PCNM13 + PCNM16 + PCNM6, bac.space)[,-1]

# Use partial constrained ordination, which requires two explanatory matrix instead of one (2 layers of parameter data)
bac.total.env <- dbrda(bac.bc ~ bac.env.mod)
bac.total.space <- dbrda(bac.bc ~ space.bac.mod)
bac.env.cond.space <- dbrda(bac.bc ~ bac.env.mod + Condition(space.bac.mod))
bac.space.cond.env <- dbrda(bac.bc ~ space.bac.mod + Condition(bac.env.mod))

permutest(bac.env.cond.space, permutation = 999)
permutest(bac.space.cond.env, permutation = 999)
permutest(bac.total.env, permutation = 999)
permutest(bac.total.space, permutation = 999)

bac.varpart <- varpart(bac.bc, bac.env.mod, space.bac.mod)
bac.varpart

par(mar = c(2,2,2,2))
plot(bac.varpart)
text(1, 0.25, "Space")
text(0, 0.25, "Env")
mtext("Variation Partitioning of \nDoubs Fish Diversity", side = 3, line = -3)
