# This script demonstrates the use of "rlc"

# Install rlc, and the supporting library, JsRCom, from our github page:
devtools::install_github("anders-biostat/JsRCom" )
devtools::install_github("anders-biostat/rlc")

# Load the rlc library
library( rlc )

# Load the example data, which has been prepared with the script "data_prep.R"
load( "citeseq.rda" )

# Open a page in the browser, with space for 2x2 plots, labelled A1, A2, B1, B2
openPage( layout = "table2x2", useViewer=FALSE )

# The gene to show in the second plot
gene <- "CD79A"

# First plot: variance vs means;
# clicking on a point changes 'gene' and redraws the other two plots
lc_scatter(
  dat( 
    x = means, 
    y = vars / means,
    logScaleX=10, 
    logScaleY=10, 
    size=2.5,
    on_click = function(i) {
      gene <<- names(means)[i]
      updateChart( "A2" )
      updateChart( "B1" )
    } ), 
  "A1"
)

# Second plot: raw UMI counts for the selected gene versus total count sums 
lc_scatter(
  dat( 
    x = countsums, 
    y = counts[ gene, ], 
    title = gene, 
    size = 1 ),
  "A2"
)

# Third plot: a t-SNE plot (coordiantes pre-calculted with Rtsne), 
# coloured accoding to expression of the selcted gene.
lc_scatter(
  dat( 
    x = tsne$Y[,1], 
    y = tsne$Y[,2], 
    colourValue = sqrt( counts[ gene, ] / sf ),
    palette = RColorBrewer::brewer.pal( 9, "YlOrRd" ),
    size = 1 ),
  "B1"
)


# If one marks points in the t-SNE plot (by dragging the mouse while 
# holding the Shift key), one can get a list of the marked points
# like this:

getMarked("B1")

# Once we have marked some cells, we can get a list of genes that distinguish
# these cells from all the others

# For this, we work with Anscombe (i.e., square root) transformed, then, 
# L2-normalized counts. 
# To speed things up a bit, we use only the highly variable genes
# chosen earlier, when the t-SNE was calculated in 'data_prep.R'
exprs <- t( apply( counts[ variableGenes, ], 1, function(x) sqrt( x / sf ) ) )

# Now, we calculate mean and standard deviation for all genes
# across the marked and across the other cells.

library( genefilter )

exprsMarked <- exprs[ , getMarked("B1") ]
exprsUnmarked <- exprs[ , -getMarked("B1") ] 
df <- data.frame(
  meanMarked = rowMeans( exprsMarked ),
  meanUnmarked = rowMeans( exprsUnmarked ),
  sdMarked = rowSds( exprsMarked ),
  sdUnmarked = rowSds( exprsUnmarked )
)

# We calculate a "separation score", simple the difference in mean devided by the 
# sum of the standard deviations, to find genes which separate the marked cells 
# well from the others
df$sepScore <- ( df$meanMarked - df$meanUnmarked ) / pmax( df$sdMarked + df$sdUnmarked, 0.002 )

# Show the genes with highest separation score
head( df[ order(-df$sepScore), ], 20 )

# The 'hwriter' package allows to turn this table into HTML code
hwriter::hwrite( head( df[ order(-df$sepScore), ], 20 ) )

# In 'rlc_example_2.R', we show how this table can be sent top appear on the browser