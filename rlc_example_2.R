# Repeat all the steps explained in 'rlc_example.R'

library( rlc )

load( "citeseq.rda" )

openPage( layout = "table2x2", useViewer=FALSE )

gene <- "CD79A"

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

lc_scatter(
  dat( 
    x = countsums, 
    y = counts[ gene, ], 
    title = gene, 
    size = 1 ),
  "A2"
)

# For the third chart, we now add a new property, a function 'markedUpdated',
# which will be called whenever the user marks cells. It gets the indexes of
# the marked points and then calls the 'showHighGenes' function that we will
# define below.

lc_scatter(
  dat( 
    x = tsne$Y[,1], 
    y = tsne$Y[,2], 
    colourValue = sqrt( counts[ gene, ] / sf ),
    palette = RColorBrewer::brewer.pal( 9, "YlOrRd" ),
    size = 1,
    markedUpdated = function() {
       showHighGenes( getMarked( "B1" ) )
    }),
  "B1"
)

exprs <- t( apply( counts[ variableGenes, ], 1, function(x) sqrt( x / sf ) ) )

showHighGenes <- function( marked ){
  
  # If no genes are marked, clear output, do nothing else
  if( length(marked) == 0 ) {
    lc_html( "", "B2" )
    return()
  }

  exprsMarked <- exprs[ , marked ]
  exprsUnmarked <- exprs[ , -marked ] 
  df <- data.frame(
    meanMarked = rowMeans( exprsMarked ),
    meanUnmarked = rowMeans( exprsUnmarked ),
    sdMarked = rowSds( exprsMarked ),
    sdUnmarked = rowSds( exprsUnmarked )
  )

  df$sepScore <- ( df$meanMarked - df$meanUnmarked ) / pmax( df$sdMarked + df$sdUnmarked, 0.002 )

  html_table <- hwriter::hwrite( head( df[ order(-df$sepScore), ], 20 ) )
 
  # Now, send the html table to appear in slot B2
  lc_html( html_table, "B2" )
}  