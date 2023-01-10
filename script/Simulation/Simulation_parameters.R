# Simulation - goal is to test true highly variable genes
# Scdd different number of hvgs and different level of sparsity 
library(Seurat)
library(splatter)
library(scater)
# Both methods based on overall gene expression variability, was not able to generate the true highly variable genes
# use TMS HSC facs data 
sce <- readRDS("HSC_droplet.rds")

# Estimate parameters from mock data
params <- splatEstimate(as.matrix(sce@assays$RNA@counts))

sim <- splatSimulate(params, nGenes = 10000)
sim2 <- splatSimulateGroups(params,
                            batchCells=400,
                            nGenes = 10000,
                            group.prob = c(0.2, 0.8),
                            de.prob = 0,
                            de.facScale = c(0.6,0.4),
                            verbose = FALSE)
sim2 <- logNormCounts(sim2)
sim2 <- runPCA(sim2)
plotPCA(sim2, colour_by = "Group") 

sim3 <- splatSimulateGroups(params,
                            batchCells=400,
                            nGenes = 10000,
                            group.prob = c(0.2, 0.8),
                            de.prob = 0,
                            de.facScale = c(0.8,0.4),
                            verbose = FALSE)

sim3 <- logNormCounts(sim3)
sim3 <- runPCA(sim3)
plotPCA(sim3, colour_by = "Group") 

sim2.para <- splatEstimate(as.matrix(sim2@assays@data@listData[["TrueCounts"]]))
sim4 <- splatter:::splatSimDropout(sim2, setParam(params, "dropout.mid", 2))
sim4 <- logNormCounts(sim4)
sim4 <- runPCA(sim4)
plotPCA(sim4, colour_by = "Group") 


sim5 <- splatter:::splatSimDropout(sim3, setParam(params, "dropout.mid", 2))
sim5 <- logNormCounts(sim5)
sim5 <- runPCA(sim5)
plotPCA(sim5, colour_by = "Group") 

saveRDS(sim2, "Simulation1.rds")
saveRDS(sim3, "Simulation2.rds")
saveRDS(sim4, "Simulation3.rds")
saveRDS(sim5, "Simulation4.rds")
