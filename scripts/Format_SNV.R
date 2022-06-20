library(data.table)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

snv = as.data.frame( fread( file.path(input_dir, "SNV.txt.gz") , stringsAsFactors=FALSE , sep=";" ))

data = cbind( snv[ , c("patient" , "gene" ) ] , "Missense_Mutation" , NA , NA , NA , NA , NA )

colnames(data) = c( "Sample" , "Gene" , "Effect" , "Pos" , "Chr", "Ref" , "Alt", "MutType" )
data$Sample  = paste( "P" , data$Sample  , sep="" )

for( i in 1:nrow(snv)){
	mut = snv$mutation[i]

	if( length( grep("ins" , mut) ) ){
		if( grep("ins" , mut) ){
		 	data$Ref[i] = "_"
		 	data$Alt[i] = unlist( strsplit( unlist( strsplit( mut ,"ins" ))[2], "_" ))[1]
	 	}
	} else{
		if( length( grep("del" , mut) ) ){
			if( grep("del" , mut) ){
			 	data$Ref[i] = unlist( strsplit( unlist( strsplit( mut ,"del" ))[2], "_" ))[1]
			 	data$Alt[i] = "_"
			 }
		} else{
		 	data$Ref[i] = unlist( strsplit( unlist( strsplit( mut ,">" ))[1], "[0-9]+" , perl=TRUE ))[2]
		 	data$Alt[i] = unlist( strsplit( unlist( strsplit( mut ,">" ))[2], "_" ))[1]
		} 		
	}
}

data$Ref = ifelse( data$Ref %in% "-" , "" , data$Ref )
data$Alt = ifelse( data$Alt %in% "-" , "" , data$Alt )

data$MutType =  apply( data[ , c( "Ref", "Alt" ) ] , 1 , function(x){ ifelse( nchar(x[1]) != nchar(x[2]) , "INDEL", "SNV") } )


case = read.csv( file.path(output_dir, "cased_sequenced.csv"), stringsAsFactors=FALSE , sep=";" )
data = data[ data$Sample %in% case[ case$snv %in% 1 , ]$patient , c( "Sample" , "Gene" , "Chr" , "Pos" , "Ref" , "Alt" , "Effect" , "MutType" ) ]

write.table( data , file= file.path(output_dir, "SNV.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )
