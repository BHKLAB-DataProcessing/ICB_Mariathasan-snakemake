args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/Get_Response.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/format_clin_data.R")

clin_original = read.csv( file.path(input_dir, "CLIN.txt"), stringsAsFactors=FALSE , sep="\t")

selected_cols <- c( "Best.Confirmed.Overall.Response","Sex","ANONPT_ID","os","censOS","Tissue" )
clin = cbind( clin_original[ , selected_cols  ] , "PD-1/PD-L1" , NA , NA , NA , NA , NA , NA , NA , NA , NA )
colnames(clin) = c( "recist" , "sex" , "patient" , "t.os"  ,"os" , "primary" , "drug_type" , "age" , "t.pfs", "pfs"  , "histo" , "stage" , "dna" , "rna" , "response.other.info" , "response" )

clin$recist[ clin$recist %in% "NE"] = NA
clin$response = Get_Response( data=clin )
clin$patient = paste( "P" , clin$patient , sep="" )
clin_original$ANONPT_ID <- paste( "P" , clin_original$ANONPT_ID , sep="" )

case = read.csv( file.path(output_dir, "cased_sequenced.csv"), stringsAsFactors=FALSE , sep=";" )
clin$rna[ clin$patient %in% case[ case$expr %in% 1 , ]$patient ] = "tpm"
clin$dna[ clin$patient %in% case[ case$cna %in% 1 , ]$patient ] = "tgs"

clin$primary = ifelse( clin$primary %in% "bladder" , "Bladder" , 
				ifelse( clin$primary %in% "kidney" , "Kidney" ,
					ifelse( clin$primary %in% "liver" , "Liver" ,
					ifelse( clin$primary %in% "lung" , "Lung" ,
					ifelse( clin$primary %in% "lymph node" , "Lymph_node" ,
					ifelse( clin$primary %in% "ureter" , "Ureter" , 
					ifelse( clin$primary %in% "other" , "Unknown" ,  NA )))))))

clin$primary[ is.na(clin$primary) ] = "Unknown"

clin = clin[ , c("patient" , "sex" , "age" , "primary" , "histo" , "stage" , "response.other.info" , "recist" , "response" , "drug_type" , "dna" , "rna" , "t.pfs" , "pfs" , "t.os" , "os" ) ]
clin$primary[ clin$primary %in% "Ureter" ] = "Ureteral"

clin <- format_clin_data(clin_original, 'ANONPT_ID', selected_cols, clin)

write.table( clin , file=file.path(output_dir, "CLIN.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )

