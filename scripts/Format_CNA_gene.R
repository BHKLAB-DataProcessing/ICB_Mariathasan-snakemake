library(data.table)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

case = read.csv( file.path(output_dir, "cased_sequenced.csv"), stringsAsFactors=FALSE , sep=";" )

cna = as.data.frame( read.csv( file.path(input_dir, "CNA_gene.txt") , stringsAsFactors=FALSE , header=TRUE , sep="\t" ) )
# colnames(cna) =  sapply( colnames(cna) , function( x ) paste( "P" , unlist( strsplit( x , "X" , fixed=TRUE ) )[2] , sep="" ) ) 

cna = cna[ , colnames(cna) %in% case[ case$cna %in% 1 , ]$patient ]

write.table( cna , file= file.path(output_dir, "CNA_gene.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=TRUE )
