library("IMvigor210CoreBiologies")
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]

##################################
## To load a Genomic alterations as assessed by FMOne panel (Foundation Medicine, Inc.)
data(fmone)

##################################
## Get CNA

amp = assayDataElement(fmone, "amplification")
gain = assayDataElement(fmone, "gain")
deletion = assayDataElement(fmone, "deletion")

cna = matrix(0,nrow(amp),ncol(amp))
colnames(cna) = colnames(amp)
rownames(cna) = rownames(amp)

for(i in 1:nrow(cna)){
	for(j in 1:ncol(cna)){
		if( ! amp[i,j] %in% "" ){ cna[i,j] = 2 }
		if( ! gain[i,j] %in% "" ){ cna[i,j] = 1 }
		if( ! deletion[i,j] %in% "" ){ cna[i,j] = -1 }
	}
}
colnames(cna) = paste0('P', pData(fmone)[colnames(cna),]$ANONPT_ID)

write.table( cna , file= file.path(work_dir, "CNA_gene.txt") , quote=FALSE , sep="\t" , col.names=TRUE , row.names=TRUE )

##################################
## Get SNV

known = assayDataElement(fmone, "known_short")
likely = assayDataElement(fmone, "likely_short")
colnames(known) = paste0('P', pData(fmone)[colnames(known),]$ANONPT_ID)
colnames(likely) = paste0('P', pData(fmone)[colnames(likely),]$ANONPT_ID)

snv = NULL
for(i in 1:ncol(known)){
	for(j in 1:nrow(known)){
		if( ! known[j,i] %in% "" ){ 
			patient = colnames(known)[i]
			gene = rownames(known)[j]
			mut = unlist(strsplit( known[j,i], '(' , fixed= TRUE))[1]
			vaf = unlist(strsplit( unlist(strsplit( known[j,i], '(' , fixed= TRUE))[2], ",", fixed= TRUE ))[1]
			cov = unlist(strsplit( unlist(strsplit( unlist(strsplit( known[j,i], '(' , fixed= TRUE))[2], ",", fixed= TRUE ))[2], ")", fixed= TRUE ))[1]

			snv = rbind(snv,c(patient,gene,mut,"known_short",vaf,cov))
		}
	}
}

for(i in 1:ncol(likely)){
	for(j in 1:nrow(likely)){
		if( ! likely[j,i] %in% "" ){ 
			patient = colnames(likely)[i]
			gene = rownames(likely)[j]
			mut = unlist(strsplit( likely[j,i], '(' , fixed= TRUE))[1]
			vaf = unlist(strsplit( unlist(strsplit( likely[j,i], '(' , fixed= TRUE))[2], ",", fixed= TRUE ))[1]
			cov = unlist(strsplit( unlist(strsplit( unlist(strsplit( likely[j,i], '(' , fixed= TRUE))[2], ",", fixed= TRUE ))[2], ")", fixed= TRUE ))[1]

			snv = rbind(snv,c(patient,gene,mut,"likely_short",vaf,cov))
		}
	}
}

colnames(snv) = c("patient","gene","mutation","type","vaf","coverage")
snv = snv[order(snv[,1],snv[,2]),]

gz <- gzfile(file.path(work_dir, 'SNV.txt.gz'), "w")
write.table( snv , file=gz , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )
close(gz)
# write.table( snv , file= "/data/Mariathasan/SNV.txt" , quote=FALSE , sep=";" , col.names=TRUE , row.names=TRUE )
# system("gzip --force /data/Mariathasan/SNV.txt > /data/Mariathasan/SNV.txt.gz")


