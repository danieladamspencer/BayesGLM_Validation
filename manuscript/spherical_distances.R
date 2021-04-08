library(gifti)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/Applications/workbench/')

setwd('~/Google Drive/My Drive/RESEARCH/BayesGLM_Validation/')

fname_sphere <- 'Q1-Q6_R440.L.sphere.32k_fs_LR.surf.gii'
fname_midthick <- 'Q1-Q6_R440.L.midthickness.32k_fs_LR.surf.gii'

resamplegii(fname_sphere, target_fname = 'sphere_6k.surf.gii', hemisphere = 'left', resamp_res = 6000)
resamplegii(fname_midthick, target_fname = 'midthick_6k.surf.gii', hemisphere = 'left', resamp_res = 6000)

surf_sphere <- read_gifti('sphere_6k.surf.gii')
surf_midthick <- read_gifti('midthick_6k.surf.gii')

V <- nrow(surf_sphere$data$pointset)
locs_sphere <- surf_sphere$data$pointset
dist_sphere <- dist(locs_sphere)
locs_midthick <- surf_midthick$data$pointset
dist_midthick <- dist(locs_midthick)

#establish direct neighbors
triangle <- (surf_sphere$data$triangle + 1)
neighbors <- matrix(0, V, V)
for(ii in 1:nrow(triangle)){
  tri_ii <- triangle[ii,]
  loc1 <- tri_ii[1]
  loc2 <- tri_ii[2]
  loc3 <- tri_ii[3]
  neighbors[loc1, loc2] <- neighbors[loc2, loc1] <- 1
  neighbors[loc1, loc3] <- neighbors[loc3, loc1] <- 1
  neighbors[loc3, loc2] <- neighbors[loc2, loc3] <- 1
}
neighbors[upper.tri(neighbors)] <- 0

dist_sphere_mat <- dist_midthick_mat <- matrix(0, V, V)
dist_sphere_mat[lower.tri(dist_sphere_mat)] <- dist_sphere
dist_midthick_mat[lower.tri(dist_midthick_mat)] <- dist_midthick

dist_sphere_neighbors <- dist_sphere_mat[neighbors==1]
dist_midthick_neighbors <- dist_midthick_mat[neighbors==1]

pdf('spherical_distances.pdf', width=7, height=4)
par(mfrow=c(1,2))
set.seed(5829743)
samp <- sort(sample(1:length(dist_sphere_neighbors), 1000, replace=FALSE))
plot(dist_sphere_neighbors[samp], dist_midthick_neighbors[samp], pch=19, col=rgb(0,0,0,0.25), xlim=c(3,7), ylim=c(1,8),
     xlab='Spherical Distance', ylab='Midthickness Distance')
abline(a=0, b=1, col='red')
hist(dist_sphere_neighbors/dist_midthick_neighbors, breaks=50, xlab='Spherical Distance Distortion', main='')
abline(v=1, col='red')
dev.off()
