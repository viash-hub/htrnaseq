
#' Displays the annotation of the wells in a plateLayout
#' @param plateData a data.table object containing the information of the plate. This must contain a "WellID", "Compound" and "Dose" column in which "Dose" must be a factor.
#' @param plateName The plate name
#' @param studyType A character string indicating wethere the data comes from a "CompoundProfiling" or "MiniTox" background. For both, the Compound - Dose combination will be displayed in a plate layut with a colouring scheme according to the dose levels. This will decide the colour scheme. For any other type, specify the "valueVariable", "textVariable" and "colours"  parameters. For the latter, the names should corresponding with the levels of the valueVariable to be coloured.
#' @param valueVariable The name of the variable in 'plateData' to be visualized in a plate layout. Ignored if study type is specified.
#' @param textVariable The name of the variable in 'plateData' to be shown in the wells of the plate layout. If NULL, the valueVariable is shown. Ignored if study type is specified.
#' @param colours A named character vector containing the colours for the different levels of the valuevariable. The names should correspond to the dose levels. if not specified, a scheme of blues will be provided.
#' @param breaks Numeric vector indicating breaks for plot coloring.
#' @param colourWellText Colour to display the text in the wells.
#' @param layout Integer vector of length two with number of rows and colums in a plate, e.g. \code{c(16,24)}
#' @param legend.title A title for the legend
#' @param plot.title A title for the plot, will be contracted with the plate name
#' @param ... additional arguments for \code{plateLayout.default} function
#' @import data.table
#' @importFrom platetools fill_plate
#' @export
plateLayout.annotation <- function(plateData, plateName = character(), studyType = c("CompoundProfiling", "MiniTox", "OligoProfiling", "NTCPlate"), 
                                   valueVariable = "Dose", textVariable = NULL, breaks = NULL, colours = NULL, colourWellText = "black", 
                                   layout = c(16,24), legend.title = "Dose", plot.title ="Plate Annotation - ", textFontSize = 9, ...){
  
  WellID <- Label <- Dose <- NULL
  
  if(!(all(c("WellID", "SampleName") %in% colnames(plateData)))){
    stop(" 'WellID' and 'SampleName' column required in plateData object")
  }
  
  
  plateData[, WellID := paste0(sub(".*([[:alpha:]]).+", "\\1", plateData$WellID), sprintf("%02d",as.numeric(sub(".*[[:alpha:]](.+)", "\\1", plateData$WellID))))]
  
  plateData <- platetools::fill_plate(plateData, "WellID", plate = 384)
  
  plateData$column <- factor(sprintf("%02d",as.numeric(sub(".*[[:alpha:]](.+)", "\\1", plateData$WellID))),
                             levels = sprintf("%02d", seq(1,layout[2])))
  plateData$row <- factor(sub(".*([[:alpha:]]).+", "\\1", plateData$WellID),
                          levels = LETTERS[seq(1,layout[1])])
  
  
  
  if(!is.null(valueVariable)){
    plateData[, values := as.character(plateData[,..valueVariable][[1]])]
    valueVar <- "values"
  }else{
    plateData[, values := 'grey']
    valueVar <- "values"
    colours <- setNames('grey', 'grey')
  }
  
  
  if(is.null(colours)){
    
    blues <- colorRampPalette(c("#d6e0ff","#2171B5"))
    greens <- colorRampPalette(c("light green","dark green"))
    greys <- colorRampPalette(c("light grey","dark grey"))
    
    numLevels <- sort(as.numeric(as.character(unique(plateData[,values])[grepl("^[[:digit:]]+([.][[:digit:]]+)?$",trimws(unique(plateData[,values])))])))
    otherLevels <- sort(as.character(unique(plateData[,values])[!grepl("^[[:digit:]]+([.][[:digit:]]+)?$",trimws(unique(plateData[,values])))]))
  
    colours <- c(blues(length(numLevels)), greens(length(otherLevels)) ,"red")
    names(colours) <- c(numLevels, otherLevels, "failed")
    
  }  
  
  if(!is.null(textVariable)){
    plateData[, Label :=  do.call(paste, c(.SD, sep = "\n ")), .SDcols = textVariable]
    plateData[, Label :=  gsub("-","-\n", Label)]
    plateData[, Label :=  gsub("_","_\n", Label)]
    textVar <- "Label"
  }else{
    textVar <- NULL
  }
  
  
  if(is.null(breaks)){
    breaks <- seq_len(length(colours))
  }

  
  plateLayout(plateData = plateData, valueVariable = valueVar, textVariable = textVar, plateName = plateName, breaks = breaks, 
                      colourWellText = colourWellText, legend.title = legend.title, layout = layout, colours = colours, plot.title = 
                      plot.title, textFontSize = textFontSize, ...)
  
  

}



#' Create a heatmap of values in a plateLayout view. The values can be library sizes, number of genes, qcScore (0/1) or a factor.
#' @param plateData A data.table of the values to be visualized with at least the column of interest (specified in 'varOfInterest') and a 'WellID' column indicating the wells in the plate. The WellID is a combination of a letter (row in the plate) and an integer (column in the plate).
#' @param valueVariable The name of the variable in 'plateData' to be visualized in a plate layout
#' @param textVariable The name of the variable in 'plateData' to be shown in the wells of the plate layout. Defaults to the valueVariable and if NULL, no text will be displayed.
#' @param breaks Numeric vector indicating breaks for plot coloring.
#' @param colours Colours to be used for levels specified by the breaks. If NULL, a colour scheme of purples is shown.
#' @param colourWellText Colour to display the text in the wells.
#' @param layout Integer vector of length two with number of rows and colums in a plate, e.g. \code{c(16,24)}
#' @param makeContourColours Logical, whether or not the plate layout will contain a contour colours for the wells based on the parameters in 'contourColours' and 'categories'
#' @param contourVariable The variable used for the contour colouring
#' @param contourColours Character vector specifying a colour for each range in 'categories'
#' @param labelsCategories Character vector specifying the names (labels) for each range in 'categories'
#' @param categories if contour Variable is not a factor, a numeric vector specifying the categories to divide the 'varOfInterest', including the lower and upper limits.
#' @param plateName The plate name
#' @param plot.title A title for the plot, will be contracted with the plate name
#' @param legend.title A title for the legend
#' @param displayHeatmap Logical, whether to display the plateLayout heatmap
#' @param saveHeatmap Logical, whether to save the plateLayout heatmap
#' @param outputDir The directory where the plateLayout heatmap should be saved
#' @param prefix The prefix to the file name of the saved plateLayout heatmap
#' @param ... additional arguments for \code{ComplexHeatmap::Heatmap} function
#' @importFrom platetools fill_plate
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#' @importFrom grid grid.text grid.rect gpar legendGrob	gpar
#' @importFrom grDevices dev.off png
#' @importFrom graphics title
#' @export
plateLayout <- function(plateData, valueVariable, textVariable = valueVariable,
                        breaks = NULL, colours = NULL, colourWellText = "white", textFontSize = 6,
                        layout = c(16,24), makeContourColours = FALSE,contourVariable = character(), 
                        contourColours = c("red", "orange", "seagreen3"),lwdContours = c(1,1,1), 
                        labelsCategories = c('1', '2', '3'), categories = NULL,plateName = character(), 
                        plot.title = character(), legend.title = NULL, legendFontSize = 15, row_split = rep("A", 16), col_split = rep("A", 24),
                        legendFontSizeTitle = 15, displayHeatmap = TRUE, saveHeatmap = FALSE, outputDir = ".", prefix = "") {
  
  WellID <- NULL
  
  if(!(all(c("WellID", "SampleName") %in% colnames(plateData)))){
    stop(" 'WellID' and 'SampleName' column required in plateData object")
  }
  
  
  plateData[, WellID := paste0(sub(".*([[:alpha:]]).+", "\\1", plateData$WellID), sprintf("%02d",as.numeric(sub(".*[[:alpha:]](.+)", "\\1", plateData$WellID))))]
  
  plateData <- platetools::fill_plate(plateData, "WellID", plate = 384)
  
  plateData$column <- factor(sprintf("%02d",as.numeric(sub(".*[[:alpha:]](.+)", "\\1", plateData$WellID))),
                             levels = sprintf("%02d", seq(1,layout[2])))
  plateData$row <- factor(sub(".*([[:alpha:]]).+", "\\1", plateData$WellID),
                          levels = LETTERS[seq(1,layout[1])])
  
    
  plateValues <- plateLayoutFormat(plateData, varOfInterest = valueVariable, rows = layout[1], cols = layout[2])
  if(!is.null(textVariable)){
    plateText <- plateLayoutFormat(plateData, varOfInterest = textVariable, rows = layout[1], cols = layout[2])
  }
  plot.title <- gsub("^([a-z])","\\U\\1",gsub("([A-Z])"," \\1", plot.title, perl=TRUE), perl=TRUE)
  mainTitle <- paste0(plot.title, plateName)
  plateContourColours <- matrix("", nrow=layout[1], ncol=layout[2]) 
  
  if(makeContourColours){
    
    contourData <- plateData[WellType == "nonEmpty" | WellType == "Treated Wells" ,]
    
    if(is.numeric(contourData[, ..contourVariable][[1]])){
      contourData$contours <- cut(contourData[, ..contourVariable][[1]], categories, left = TRUE, right = TRUE, labels = labelsCategories)
    }
    else{
      contourData$contours <- contourData[, ..contourVariable][[1]]
    }
    names(contourColours) <- labelsCategories
    names(lwdContours) <- labelsCategories
    for (i in seq_len(layout[1])) {
      for (j in seq_len(layout[2])) {
        
        tryCatch({
          
          sampleHit <- which(as.character(contourData$WellID)==paste0(LETTERS[i], sprintf("%02d",j)))  
          if(length(sampleHit)==1){
            plateContourColours[i,j] <- as.character(contourData[sampleHit,'contours'][[1]])
          }
          
        },
        error = function(e) {
          print(paste0(LETTERS[i], sprintf("%02d", j), " is missing."))
        }
        )
      }
    }
  }
  
  plateValues$contours <- plateContourColours
  colnames(plateValues$values) <- seq_len(ncol(plateValues$values))
  
  
  if(is.null(breaks)) {
    breakValues <- plateValues$values
    breakValues[which(is.na(breakValues))] <- 0 
    if(all(breakValues>=0)){
      breaks <- computeBreaks(7, max(plateValues$values, na.rm=TRUE))
    }else{
      breaks <- quantile(plateValues$values,  probs = seq(0, 1, 0.125))
    }
  }
    
  if(is.null(colours)){
      colours <- tryCatch({
        
        colorRamp2(breaks = breaks,colors = brewer.pal(length(breaks), "Purples"))
      },
      error = function(cond){
        return(c("#9370DB","white"))
      }
      )
    }
    
    ht <- Heatmap(plateValues$values, 
                  column_title=mainTitle, column_title_side="top",
                  rect_gp=gpar(lwd=0.4),
                  cluster_rows=FALSE, cluster_columns=FALSE,    
                  col = colours, row_title = NULL,
                  row_split = row_split, column_split = col_split,
                  row_names_side="left",
                  cluster_row_slices = FALSE,
                  cluster_column_slices = FALSE,
                  show_heatmap_legend = TRUE,
                  heatmap_legend_param = list(title = ifelse(is.null(legend.title),paste0(valueVariable,"\n"), paste0(legend.title,"\n")),
                                              grid_height = unit(9, "mm"), border = "black",
                                              labels_gp = gpar(fontsize = legendFontSize), title_gp = gpar(fontsize = legendFontSizeTitle)),
                  
                  cell_fun = function(j,i,x,y,width,height,fill) {
                    
                    tol = 1e-5
                    if(is.na(plateValues$values[i,j])){
                      grid.rect(x,y,width,height,gp=gpar(fill="white", alpha=0.7, lwd=0.7,col="white"))
                    }
                    else if(!is.null(textVariable)){
                     grid.text(plateText$values[i,j],x,y,just="centre", gp=gpar(fontsize=textFontSize, col=colourWellText))
                    }
                    if(makeContourColours){
                      if (!is.na(plateValues$contours[i,j])) {
                        grid.rect(x,y,width,height,gp=gpar(col=contourColours[as.character(plateValues$contours[i,j])],
                                                           fill=NA,lwd=lwdContours[as.character(plateValues$contours[i,j])]))
                      }
                    }
                  })
    
    if(displayHeatmap){
      print(ht)
    }
    
    if(saveHeatmap){
      
      png(file.path(outputDir, paste0(prefix,gsub(" |-","",plot.title),"_",plateName,".png")), width = 30, height = 10, units ="cm", res =1200)
      print(ht)
      dev.off()
    }
    
    return(ht)
  }
  


#' Return numerical matrix with number of reads that corresponds to the plate layout
#' @param data A data.frame of the values to be visualized with at least the column of interest (specified in 'varOfInterest') and a 'WellID' column indicating the wells in the plate. The WellID is a combination of a letter (row in the plate) and an integer (column in the plate).
#' @param varOfInterest The name of the variable in 'data' to be visualized in a plate layout
#' @param rows number of rows in a plate layout
#' @param cols number of columns in a plate layout
#' @param verbose if \code{TRUE}, samples missing from the plate will be reported
#' @export
plateLayoutFormat <- function(data, varOfInterest, rows=16, cols=24, verbose=FALSE) {
  
  #sampleNames <- rownames(data)
  plateValues <- matrix(NA, nrow=rows, ncol=cols)
  for (i in seq_len(rows)) {
    for (j in seq_len(cols)) {
      
      tryCatch({
        
        sampleHit <- which(as.character(data$WellID)==paste0(LETTERS[i], sprintf("%02d",j)))  
        if(length(sampleHit)==1){
          plateValues[i,j] <- data[sampleHit, ..varOfInterest][[1]]
        }
        
      },
      error = function(e) {
        if (verbose == TRUE) {
          print(paste0(LETTERS[i], sprintf("%02d", j), " is missing."))
        }
      }
      )
    }
  }
  
  row.names(plateValues) <- LETTERS[1:rows]
  return(list("values"=plateValues))
}



#' Helper function to automate break selection for raw count data
#'
#' This function creates an exponentially increasing vector for given number
#' breaks between zero and some element of choice. It is particularly useful for
#' raw counts or raw counts per million.
#' 
#' @param nBreaks Number of breaks to be generated
#' @param maxElement Maximum value of data entries
#' @export
computeBreaks <- function(nBreaks, variable){
  
  maxElement <- max(variable, na.rm = TRUE)
  if(length(unique(variable))==1){
    
    breaks <-  c(0,0.5,ifelse( maxElement < 1, 1,  maxElement))
    
  }else{
    
    coefSystem <- solve(rbind(c(1,1), c(1, (nBreaks-1)))) %*% c(0, log(maxElement))
    coefExp <- c(exp(coefSystem[1]), coefSystem[2])
    breaks <- coefExp[1]*exp((1:(nBreaks-1))*coefExp[2])
    
  }
  

  return(c(0, breaks))
  
}

