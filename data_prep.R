# SOUND Meeting September 2018
# Tutorial: Single-cell visualization
# Simon Anders, Svetlana Ovchinnikova

# Script 1: Data preparation

# As example data, we use the "CiteSeq" data set, which is described in the following publication:
#   M. Stoeckius et al.: "Simultaneous epitope and transcriptome measurement in single cells"
#   Nature Methods, 14, 865-8 (2017)  -  https://doi.org/10.1038/nmeth.4380
# Stoeckius et al. sequenced from human cord blood mononuclear cells (CBMCs), to demonstrate a technique
# to simultaneously assess transcriptome and epitopes on a single-cell levels.
# The raw data is available on GEO, here:
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz

# We first load the matrix of UMI counts for the RNA data. (We will not use ADT data here.)

counts <- as.matrix( read.csv( "~/Downloads/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", header=TRUE, row.names=1 ) )

# Stoeckius et al. have spiked in about 4% mouse fibroblast cells to check for dublets. 
# We first calculate for each cell the ratio of human to mouse genes.

# Which cellsrows are mouse and which human genes?
gene_species <- sapply( strsplit( rownames(counts), "_" ), `[`, 1 ) 
table( gene_species )

# Get for each cell the sum of UMI counts from either species
cell_counts_by_species <- sapply( unique( gene_species ), function(x) 
   colSums( counts[ gene_species == x, , drop=FALSE ] ) )
head(cell_counts_by_species)
hist( log2( cell_counts_by_species[,"HUMAN"] / cell_counts_by_species[,"MOUSE"] ), 30 )

# Subset the counts matrix columns to only those cells that at least 10x as many human as mouse UMIs,
# and the rows to only the human genes, then remove the "HUMAN_" prefix from the row names
counts <- counts[ 
  gene_species == "HUMAN",
  cell_counts_by_species[,"HUMAN"] / cell_counts_by_species[,"MOUSE"] > 10 ]
rownames(counts) <- sub( "HUMAN_", "", rownames(counts) )
counts[ 1:5, 1:5 ]

# For some reasons, the cells in this matrix are sorted (nearly but not quite) by total UMI count
# For plotting, it will later be better to have the cells in random order to avoid spurious
# patterns when overplotting. Hence, let's shuffle the columns.
counts <- counts[ , sample(seq_len(ncol(counts))) ]

# Now get the sum of UMIs per cell for normalization
countsums <- colSums( counts )
hist( log10( countsums ) )

# Let's center these around one to get "size factors", by which we can then divide 
# the counts to normalize across cells
sf <- countsums / exp(mean(log(countsums)))
hist( log10( sf ) )

# Now, calculate for each gene mean and variance of the normalized counts
means <- apply( counts, 1, function(x) mean( x / sf ) )
vars <- apply( counts, 1, function(x) var( x / sf ) )

#Plot this:
plot( vars ~ means, pch=".", log="xy" )

# The expected variance due to Poisson noise is equal to the mean, multiplied by
# a correction factor for the fact that each cell has a different size factor,
# which turns out to be 'mean(1/sf)'.
# Mark this expected variance with a line
xg <- 10^seq( -6, 3, length.out=100 )
lines( xg, xg * mean(1/sf), col="red" )

# The plot looks nicer if we use the variance-to-mean ratio as y axis:
plot( vars / means ~ means, pch=".", log="xy" )
abline( h = mean(1/sf), col="red")



# We will also later need a t-SNE plot of the cells

# For this, we pick only those gene that show some above-Poisson variation and 
# have not too low mean
variableGenes <- vars / means > 1.5  &  means > 1e-3
plot( vars / means ~ means, pch=".", log="xy", col = 1+variableGenes )

# Perform a square root transformation on the normalized counts
# I prefer ths over the logarithm used more commonly.
exprs <- t( apply( counts, 1, function(x) sqrt( x / sf ) ) )

# Calculate a t-SNE plot on these
library( Rtsne )
tsne <- Rtsne( t( exprs[ variableGenes, ] ) )
rownames( tsne$Y ) <- colnames( exprs )

# Here is the t-SNE plot
plot( tsne$Y, asp=1, pch="." )

# Save the data
save( counts, countsums, sf, means, vars, tsne, variableGenes, 
         file = "citeseq.rda" )

