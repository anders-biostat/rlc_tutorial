# Exploring single-cell data with R/LinkedCharts

This tutorial demonstrates how our R/LinkedCharts package can be used to created linked scatter charts. 

## Example data

As example data, we use single-cell transcriptomics data from the "CiteSeq" paper (Stoeckius et al., *Simultaneous epitope and transcriptome measurement in single cells*, [Nature Methods, 14:865](https://doi.org/10.1038/nmeth.4380), 2017). This paper describes a new method to simultaneously measure the transcriptome as well as the "epitome" (the abundance of selected surface markers) of single cells, using droplet sequencing. They demonstrate it with a 
sample of cord blood mononuclear cells.

Both the RNA-Seq and the epitome data are available from GEO ([GSE100866](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100866). We will use only the RNA-Seq data, which is available as file [`GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz`](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100866&format=file&file=GSE100866%5FCBMC%5F8K%5F13AB%5F10X%2DRNA%5Fumi%2Ecsv%2Egz) from the GEO page. 

The file contains a table with one row per gene and one column per cell, giving for each gene and each cell the number of transcript molecules (UMIs -- Unique Molecular Identifiers) detected in the cDNA of the cell that mapped to the gene. In the following, we will explore this data matrix.


## Data preparation

We first perform some very simple data preparation, using only base R commands. If you have worked with expression
data before, this will look completely familiar.

First, we load the data matrix

```r
# The path to the data file downloaded to GEO
fileName <- "~/Downloads/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz"

# Load the file. (This takes a while as it is a large file.)
countMatrix <- as.matrix( read.csv( gzfile( fileName ), row.names = 1) )
```

Stoecklin et al.\ have spiked in a few mouse cells to their human cord blood sample for quality
control purposes, and have therefore mapped everything against both the human and the mouse genome.
We keep only the cells with nearly only human transcripts and also remove the mouse gene counts from
the data matrix.

```r
# Calculate for each cell ratio of molecules mapped to human genes
# versus mapped to mouse genes
human_mouse_ratio <- 
  colSums( countMatrix[ grepl( "HUMAN" , rownames(countMatrix) ), ] ) / 
  colSums( countMatrix[ grepl( "MOUSE" , rownames(countMatrix) ), ] )

# Keep only the cells with at least 10 times more human than mouse genes
# and keep only the counts for the human genes.
countMatrix <- countMatrix[ 
  grepl( "HUMAN" , rownames(countMatrix) ), 
  human_mouse_ratio > 10 ]

# Remove the "HUMAN_" prefix from the gene names
rownames(countMatrix) <- sub( "HUMAN_", "", rownames(countMatrix) )
```

Now, we have a matrix of counts that looks like this

```r
> # size of the matrix
> dim(countMatrix)
[1] 20400  8005

> # top left corner of the matrix
> countMatrix[ 1:5, 1:5 ]
         TACAGTGTCTCGGACG GTTTCTACATCATCCC ATCATGGAGTAGGCCA GTACGTATCCCATTTA ATGTGTGGTCGCCATG
A1BG                    0                0                0                0                0
A1BG-AS1                0                0                0                0                0
A1CF                    0                0                0                0                0
A2M                     0                0                0                0                0
A2M-AS1                 0                0                0                0                0
```

Each column corresponds to a cell, identified by its DNA barcode, and each row to a gene. The 
numbers are the number of transcript molecules (UMIs -- unique molecular identifiers) from a
given cell that were mapped to the given gene.

Next, we calculate the sums of the counts for each cell, and divide them by 1000. We will call these
the "size factors" and use them for normalization: when comparing the counts from two cells, we will
consider them as comparable after dividing them by their respective size factors.

```r
sf <- colSums(countMatrix) / 1000
```

For our first example, we will ask whcih genes actually differ in an inmteresting way from cell to cell.
To this end, we calculate for each gene the mean and variance of the normalized counts

```r
means <- apply( countMatrix, 1, function(x) mean( x / sf ) )
vars <- apply( countMatrix, 1, function(x) var( x / sf ) )
```

We save these variable into a file, to be used later

```r
save( countMatrix, sf, means, vars, file = "citeseq_data.rda" )
```

## Analysis the old-fashioned way

Here is a plot of the variance-to-mean ratios versus the means

```r
plot( means, vars / means, 
   pch = ".", log = "xy", col = adjustcolor( "black", 0.4 ) )

abline( h = mean( 1/sf ), 
	col = adjustcolor( "red", 0.5 ) )
```

![](var_mean_simple.png)

All the genes that are scattering close to the red line have no or only little signal:
the Poisson noise dominated over any biological signal of interest that may there be.
Only the genes that are clearly above the red line can convey useful information.

(To understand why the red line marks the expected strength of Poisson noise, and why it
is has to be placed at `mean( 1/sf )`, please see our write-up at ... [Sorry, still under preparation at the moment.])

We can see groups of genes emerging from the horizontal clouds at different points. Which genes are these?
Are they active in only few cells or in many cells? 

It would be helpful if, for each gene, we could get a plot like this one

```r
gene <- "CD79A"

plot( sf, countMatrix[ gene, ], 
   log = "x", col = adjustcolor( "black", 0.4 ), main = gene )
```

![](CD79A_simple.png)

But where is gene "CD79A" on the previos plot? And, given a specific point on the overview plot (the first one), how can
we quickly get the detail plot (the second one)?

Here we can click on a point in the left chart, to change what is displayed in the right plot:

![](lc_upper_mockup.png)

**[[Todo: Replace this with a JS/LinkedChart app that mimics the R/LinkedChart one made below]]**

Now, we show how to create such linked charts with R/LinkedCharts.

## Using LinkedCharts

If you haven't done so yet, install the R/LinkedChart package with
```r
devtools::install_github( "anders-biostat/rlc" )
```

Load the R/LinkedCharts library
```r
library( rlc )
```

To produce the first scatter plot, use

```r
openPage( useViewer=FALSE )

lc_scatter(
  dat( 
    x = means, 
    y = vars / means ),
    logScaleX=10, 
    logScaleY=10,
  ), 
  "A1"
)
```
The `openPage` function opens a new page, either in your web browser (for `useViewer=FALSE`) or
in the "Viewer" pane of RStudio (for `useViewer=TRUE`, the default). We have chosen to use a browser 
window, as the viewer pane is a bit small to show two plots side by side.

The `lc_scatter` function inserts a scatter plot. The first argument sepcifies values for the plots parameters. Unsurpisingly, we have to pass two vectors with x and y coordinates of the points, and furthermore, we specify
a few more optional parameters: Both axes should be logarithmic, with tic marks labelling multiples
of 10. Note how all the parameters are wrapped in a call to a function named `dat`. This will be explained later.

The second parameter ("`place`") specifies where to put the plot on the page: here, at position `A1`, i.e., top left.

We will explain further details below; before, we add the second plot:

```r
gene <- "CD79A"

lc_scatter(
  dat( 
    x = sf, 
    y = countMatrix[ gene, ], 
    logScaleX = 10, 
    title = gene,
    opacity = 0.4 ),
  "A2"
)
```

As before, we have used `lc_scatter` to insert a scatter plot, now at position `A2` (to the right of `A1`). As before,
we use the vector `sf` for the `x` axis, and one row of `countMatrix` for the `y` axis. The variable `gene` contains
the rowname of the row to display. We also use two further optional parameters: `title` sets the plot title and maked
the gene name (as stored in the variable `gene`) appear above the chart, and `opacity = 0.2` renders the points semi-transparent, so that we can better see whether several points are sitting on top of each other.

Now, we need to *link* the two charts. When the user clicks on a point in the left chart (chart A1), then the value of the variable `gene` should be changed to the name of the gene onto which the user has clicked. If we then cause the second chart (A2) to be redrawn, it will show the data for that gene. To explain to R/LinkedChart that thsi should happen
when the user clicks on a data point, we change the code for the firts chart to the following:

```r
lc_scatter(
  dat( 
    x = means, 
    y = vars / means ),
    logScaleX=10, 
    logScaleY=10,
    on_click = function( i ) {
       gene <<- rownames(countMatrix)[ i ]
       updateChart( "A2" )
    }
  ), 
  "A1"
)
```

As you can see, we have merely added one more optional parameter to the `dat` block. It is called `on_click` and set
to a function. This function will now be called whenever the user clicks on a point and it gets passed, as its argument
`i`, the index of this point. We want ``gene`` to contain not the number of the gene's row but its name. Hence, we
simply look up the name of the gene as the `i`-th rowname of count matrix and write this into the variable. 

Note that
we assign to the variable `gene` using the `<<-` operator. This variant of the usual `<` (or `=`) operator assigns to
a global rather than a local variable. Had we used the normal `<-`, R would have created a local variable called `gene`
within the function and forgot it immediatly afterwards, rather than changing the global variable that we had defined before.

After changing the variable, we simply instruct R/LinkedCharts to update chart A2: `updateChart( "A2" )`. 

*This is all we need to link the two charts.*

Here, again, is the complete code. It is not much longer than the R code one would have written for normal static plots: Adding interactivity is easy with R/LinkedCharts.

```r
library( rlc )

load( "citeseq_data.rda" )

openPage( useViewer=FALSE )

gene <- "CD79A"

lc_scatter(
  dat( 
    x = means, 
    y = vars / means ),
    logScaleX=10, 
    logScaleY=10,
    on_click = function( i ) {
       gene <<- rownames(countMatrix)[ i ]
       updateChart( "A2" )
    }
  ), 
  "A1"
)

lc_scatter(
  dat( 
    x = sf, 
    y = countMatrix[ gene, ], 
    logScaleX = 10, 
    title = gene,
    opacity = 0.4 ),
  "A2"
)
```


How does it work? The call to `updateChart` causes
all the code within the `dat` function of chart A2 to be re-evaluated, and then, the new data to be sent to the browser. As `gene` has changed, the re-evaluation of the code in the `dat` call will now give different values to `y` and to
`title`, and the plot changes. This is why all the plot parameters have to be placed inside the `dat` call: The `dat` function is a little trick to keep R from evaluating this code straight away when you call lc_scatter but rather to
keep the code as is in order to be able to reevaluate it whenever needed. (In other word, `dat` just defines an anonymous
function without arguments, but does so in a manner that allows for a ligh, unobstructive syntax.)


## Adding a t-SNE plot

So far, our app demonstrated well the principle of using R/LinkedChart, but practitioners of single-cell transcriptomics might still find it unconvincing. Hence, let's add a bit more to make it actually useful.

In single-cell RNA-Seq, one often employs dimension reduction methods, such as [t-SNE](https://lvdmaaten.github.io/tsne/). These produce a plot with points representing the cells such that cells with similar transcriptome are closeby.

Let's use the `Rtsne` package to calculate a t-SNE plot for our expression data.

```r
install.packages( "Rtsne" )   # if needed

library( Rtsne )

# Normalize the counts by dividing by the size factors and take the square root
# (The square root stabilizes the variance; please see our write-up on 
# this [under preparation] for details, if you are interested.)
expr <- apply( countMatrix, 2, function( x ) sqrt( x / sf ) )

# Calculate the t-SNE embedding. This may take a while.
tsne <- Rtsne( expr )

# Delete the 'expr' matrix again, as it takes up a lot of memory.
rm( expr )
```

The embedding is in the `Y` field of the object returned by `Rtsne`.
We display it with a third scatter plot:

```r
lc_scatter(
  dat( 
    x = tsne$Y[,1], 
    y = tsne$Y[,2], 
    colourValue = sqrt( countMatrix[ gene, ] / sf ),
    palette = RColorBrewer::brewer.pal( 9, "YlOrRd" ),
    size = 1,
  "B1"
)
```

This plot is now placed a position "B1", i.e., under A1. We use the two columns of the
2x*N* matrix in `tsne$Y` for x and y. Each data point represents a cell, and we use
the points' colours to indicate the expression strength of the currently selected gene.
To this end, we read the corresponding row from the countMatrix the same way as we
did to get `y` in chart B1, but now also normalize by dividing by the size factors. Furthermore,
we take teh square root; this is not really necessary but results in a more even use of the 
colour palette. 

For the colour palette, we set the parameter `palette`. We have to provide a vector of colours, and
LinkedCharts will interpolate between these. Here, we have used the "YlOrRd" (for "Yellow-Orange-Red")
palette provided by the RColorBrewer package.

Note how it easy it is again to link this chart with the existing ones. If the use clicks on a point
on chart A1 to select a gene, the `gene` variable will be changed and this will affect what data is
shown now in both the charts A2 and B1. The only thing that is missing is that we need to change
the even handler in A1 to also update B2. So, we simply add this to the call to `updateChart` above.

```r
...
    on_click = function( i ) {
       gene <<- rownames(countMatrix)[ i ]
       updateChart( c( "A2", "B1" ) )
    }
...
```

Have a look [here](LINK_MISSING) for the final code, if you are unsure how this fits in.

Play a bit with the app. Now you can easily find out for any gene with high variability whether it is
expressed throughout, or whether is is specific to one of the clusters visible in the t-SNE embedding.
This will help us to understand the t-SNE plot.

## 


## What else is there?

R/LinkedCharts has many more possibilities than this simple example could demonstrate. Here are a few:

- There are many **more options** that you can put in `dat`. We have seen the mandatory parameters `x` and `y`, and
the optional parameters `logScaleX` and `logScaleY` for logarithmic axes, `size` and opacity to control the 
appearance of the points, and `on_click` to specify how to react to user interaction.

- **Point appearance** can be specified point by point, by passing a vector of values instead of a single value. Of 
course, you can also change the **shape** and **colour** of the points (parameters **todo: fill in**)

  See here [put a link to the t-SNE tutorial] for some examples how to use this.

- The scatter charts themselves have functionality for **zooming, panning, saving plots and marking points**. Press the  
arrow button at the chart's top right corner to see.

- There are not only scatter charts: You can also add **lines** to your plots, or use a **beeswarm** chart.

- There are **heatmaps** (with and without **clustering**), with many possibilities, described in their own tutorial [hopefully written soon]

- You can create your own HTML page and embed the LinkedCharts plots inside. All you have to do is to hand the HTML file
to `openPage` and to mark the places where the charts should go with `id` tags (e.g. by putting an empty `div` tag with an ID as placeholder). Then, you can use the IDs for the `place` argument (instead of `A1`, `A2`, etc). This allows you 
to make beautiful apps with little work.

- If you want a web app that runs independent of R, you can use plain LinkedCharts instead of R/LinkedCharts. You have to write in JavaScript instead of in R, but you only have know enough to load the data and write code as simple as what we have written in R above.


## Where to find more information

**[TO DO: Links to further documentation go here.]**


## Some of the ideas behind R/LinkedCharts

### A minimalistic model-view-controller pattern

If you have experience with GUI programming, you might be familiar with the [model-view-controller](https://en.wikipedia.org/wiki/Model%E2%80%93view%E2%80%93controller) (MVC) design pattern. What we have just built 
is a most minimalistic version of this.

- First, we have a number of global variables: `countMatrix`, `means` and `vars` contain the data that we want
to show in the two charts, and `gene` specifies for which gene we want to see details in the second plot. In the 
MVC concept, this is the **model**. For complex applications, one would typically implement a class to hold all the
data that makes up the model. As R/LinkedChart is designed to allow you quickly and simply get a job done, we suggest
to simply use a few global variables, as we have done here. If you want to, you can also use a more sophisticated and
organized model. After all, you can put whatever code you like into the `dat` call. You just explain to R there how and where to find your data.

- Then, we have the calls to R/LinkedCharts functions like `lc_scatter`, with its core part, the `dat` call. They
describe how each chart should look like and how it finds its data in the model. These are the **views**.

- In the traditional MVC pattern, the **controller** is the code that handles user interactions. We have not made it a separate part but rather folded it into each chart, specifically, into its event handler parameters such as the `on_click` function above. To keep things simple, we suggest that these controller codes directly manipulate some subsets of the global variables (the model) and then take care of telling the appropriate charts to update themselves. For complex GUIs, this might seem messy, but for the quick-and-not-too-dirty programming style that one wants to use during exploratory data analysis, it is suitable -- and very straight-forward.


### LinkedCharts and JsRCom

The work horse under the hood of R/LinkedCharts is our LinkedCharts library, which is written in JavaScript, not in R. It can also used directly, in order to write stand-alone client-side data exploration apps, and offers all the functionality of R/LinkedCharts and some more. R/LinkedCharts uses our JsRCom library [documentation/tutorial to be written] to establish a bidirectional [WebSocket](https://en.wikipedia.org/wiki/WebSocket) communication between the R session on one end and the web browser session with its JavaScript interpreter on the other end. The WebSocket link is used to send data and commands (e.g., to add charts to the page) from R to JavaScript and to send event notifications about user actions (such as a mouse click) from JavaScript to R. 