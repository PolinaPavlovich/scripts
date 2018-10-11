# install all needed packages 
# #cran
# install.packages(c("ggplot2", "Matrix", "dbscan", "pheatmap", "fpc", "zoo",
#"dynamicTreeCut", "factoextra", "digest", "RColorBrewer", "doParallel", 
#"dplyr", "matrixStats"))
# #bioconductor
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("BiocParallel", "scran", "scater", "monocle", 
#"KEGGREST", "org.Mm.eg.db", "biomaRt", "AnnotationDbi"))
# biocLite("SingleCellExperiment", dependencies=TRUE, lib="~/R/library")

# setting necessary parameters
dataDirectory <- "/g/lancrin/People/Pauline/data/Shiny_Conclus/Example_full_workflow"
experimentName <- "Bergiers"

# loading the countMatrix
filename <- file.path(dataDirectory, "GSE96982_countMatrix.txt")
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE96982&format=file&file=GSE96982%5FcountMatrix%2Etxt%2Egz",
              destfile = paste0(filename, ".gz"))
system(paste0("gunzip ", paste0(filename, ".gz")))

countMatrix <- read.delim(filename, stringsAsFactors = FALSE)
rm(filename)

# loading the colData downloaded from the website
# https://www.ncbi.nlm.nih.gov/geo/browse/?view=samples&series=96982
# by clicking on Export -> Amount: All search results and Format: Tab
# and moving this file to the dataDirectory

colData <- read.delim(file.path(dataDirectory, "sample.tsv"))

# preparing the data for CONCLUS
rownames(countMatrix) <- countMatrix[,1]
countMatrix <- countMatrix[,-1]

colData$state <- gsub(":.*", "", colData$Title)
colData$cellBarcode <- gsub(".*:", "", colData$Title)

# replacing cellBarcodes with c1, c2 because it is easier to work with
colData$cellName <- paste0("c", 1:nrow(colData))
rownames(colData) <- colData$cellName
countMatrix <- countMatrix[,colData$cellBarcode]
stopifnot(all(colData$cellBarcode == colnames(countMatrix)))

colnames(countMatrix) <- colData$cellName

colData <- colData[,13:15]

# deleting controls
colData <- colData[grepl("E_minus", colData$state) |
                      grepl("E_plus", colData$state) |
                      grepl("i8TFs_minus", colData$state) |
                      grepl("i8TFs_plus", colData$state),]
countMatrix <- countMatrix[,colData$cellName]

# setting seed for reproducibility
set.seed(42)
# loading functions
source("/g/lancrin/People/Pauline/data/Shiny_Conclus/conclus_bioconductor/visualisation_and_clustering_functions.R")

### CONCLUS workflow ###
# 1. Normalisation
sceObject <- normaliseCountMatrix(countMatrix, species = "mmu", 
                                  method = "default",
                                  colData = colData)

# checking what changed after the normalisation
dim(sceObject)
coldataSCE <- as.data.frame(colData(sceObject))
rowdataSCE <- as.data.frame(rowData(sceObject))

# substituting ENSEMBL IDs with SYMBOLs for well-annotated genes
rownames(sceObject)[!is.na(rowData(sceObject)$SYMBOL)] = 
    rowData(sceObject)$SYMBOL[!is.na(rowData(sceObject)$SYMBOL)]

# deleting hemoblobin genes (recommended only hor hematopoietic systems
# where these genes usually are detected at a very high level)
sceObject <- sceObject[!grepl("Hba", rownames(sceObject)) &
                           !grepl("Hbb", rownames(sceObject)),]
dim(sceObject)

# 2. Test step (optional)
testClustering(sceObject, dataDirectory, experimentName)

# 3. Running the analysis and saving the consensus clustering solution
sceObjectCONCLUS <- runCONCLUS(sceObject, dataDirectory, experimentName, 
                               plotPDFcellSim = TRUE, # FALSE for > 2500 cells
                               deepSplit = 2,
                               cores = 14, # 14 for servers, 1 for PC
                               statePalette = c("bisque", "cadetblue2", 
                                                "coral1", "cornflowerblue"),
                               deleteOutliers = FALSE  # TRUE takes more time
                               )
exportClusteringResults(sceObjectCONCLUS, dataDirectory, experimentName, 
                        "clusters_table.tsv")

# 4. Correcting clustering manually (optional)
sceObjectCONCLUS <- addClusteringManually(fileName = "clusters_table.tsv", 
    dataDirectory = dataDirectory, 
    experimentName = experimentName,
    sceObject = sceObjectCONCLUS, 
    columnName = "clusters")

# 4.1 Redo the analysis with manual clustering
sceObjectCONCLUS <- runCONCLUS(sceObjectCONCLUS, dataDirectory, experimentName, 
                        preClustered = T,
                        cores = 14, # 14 for servers, 1 for PC
                        statePalette = c("bisque", "cadetblue2", 
                                         "coral1", "cornflowerblue"))

# 5. Plotting heatmaps
markersClusters <- getMarkerGenes(dataDirectory, sceObjectCONCLUS, 
                                  experimentName = experimentName,
                                  genesNumber = 30)
orderClusters <- T
orderGenes <- T
meanCentered <- T
plotCellHeatmap(markersClusters, sceObjectCONCLUS, dataDirectory, 
                experimentName, 
                paste0("clusters",
                       length(levels(colData(sceObjectCONCLUS)$clusters)),
                       "_meanCentered",meanCentered,
                       "_orderClusters",orderClusters,
                       "_orderGenes",orderGenes), 
                meanCentered = meanCentered, 
                colorPalette = brewer.pal(10, "Paired"),
                orderClusters = orderClusters,
                orderGenes = orderGenes,
                fontsize_row = 3,
                statePalette = c("bisque", "cadetblue2", 
                                 "coral1", "cornflowerblue"),
                color = colorRampPalette(c("#023b84","#4b97fc", 
                                           "#FEE395", 
                                           "#F4794E", "#D73027",
                                           "#a31008","#7a0f09"))(100))

# 6. Export data (optional)
exportData(sceObjectCONCLUS, dataDirectory, experimentName)
