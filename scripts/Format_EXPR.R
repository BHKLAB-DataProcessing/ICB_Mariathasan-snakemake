library(data.table)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

expr = as.data.frame( fread( file.path(input_dir, "EXPR.txt.gz") , stringsAsFactors=FALSE  , sep="\t"))
rownames(expr) = expr[ , 1 ]
expr = expr[ , -1 ]

# colnames(expr)  = paste( "P" , colnames(expr)  , sep="" )

case = read.csv( file.path(output_dir, "cased_sequenced.csv"), stringsAsFactors=FALSE , sep=";" )
expr = expr[ , colnames(expr) %in% case[ case$expr %in% 1 , ]$patient ]

write.table( expr , file= file.path(output_dir, "EXPR.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=TRUE )
