#'convertMtxToH5: Converts the triple .mtx files from single cell experiments to a .h5 counts file
#' 
#' @description 
#' 
#' @details 
#' 
#' @param 
#' 
#' @return 

library(Seurat)
library(DropletUtils)

args = commandArgs(trailingOnly=T)

sampleName = args[1]
mtx_file = args[2]
features_file=args[3]
barcodes_file=args[4]

counts = Seurat::ReadMtx(mtx = mtx_file, cells = barcodes_file, features=features_file)
outfile = paste0(sampleName,".h5")
DropletUtils::write10xCounts(x= counts,path=outfile)
