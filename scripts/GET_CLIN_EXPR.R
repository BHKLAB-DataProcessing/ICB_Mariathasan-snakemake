library("IMvigor210CoreBiologies")
library(data.table)
library(biomaRt)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]

##################################
## To load a CountDataSet
data(cds)

##################################
## Get Clinical Data
clin <- pData(cds)
dup <- data.frame(table(clin$ANONPT_ID))
dup <- dup[dup$Freq > 1, ]

for(dup_id in dup$Var1){
  clin[clin$ANONPT_ID == dup_id, ]$ANONPT_ID <- unlist(
    lapply(c(1:dup[dup$Var1 == dup_id, ]$Freq), function(num){
      if(num - 1 > 0){
        return(paste0(dup_id, '.', (num - 1)))
      }else{
        return(dup_id)
      }
      
    })
  )
}

clin$ANONPT_ID <- paste0('P', clin$ANONPT_ID)

write.table(clin,file=file.path(work_dir, "CLIN.txt"),sep="\t",quote=F)

##################################
## Get Clinical Data

data = counts(cds)
annot = fData(cds)

rownames(data) = annot[rownames(data),"Symbol"]
data = data[!rownames(data)%in%"",]

##################################
##Remove Duplicated Genes
t_uniq <- data[!(rownames(data)%in%rownames(data[duplicated(rownames(data)),])),]
t_dup <- data[(rownames(data)%in%rownames(data[duplicated(rownames(data)),])),]

t_dup <- t_dup[order(rownames(t_dup)),]
id <- unique(rownames(t_dup))

t.dup.rm <- NULL
for(j in 1:length(id)){
	tmp <- t_dup[which(rownames(t_dup)%in%id[j]),]
	tmp = apply(tmp,2,mean, na.rm=T)
	t.dup.rm <- rbind(t.dup.rm,tmp)			
}
data <- rbind(t_uniq,t.dup.rm)
rownames(data) <- c(rownames(t_uniq),id)


##################################
## Get Gene Length
genes <- rownames(data)
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene_coords=getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "start_position","end_position"), filters="hgnc_symbol",values=genes, mart=human)
size=gene_coords$end_position - gene_coords$start_position
names(size) = gene_coords[,"hgnc_symbol"]

##################################
##Remove Duplicated Genes
t_uniq <- size[ !( names(size) %in% names(size)[ duplicated(names(size)) ]) ]
t_dup <- size[ ( names(size) %in% names(size)[ duplicated(names(size)) ]) ]

t_dup <- t_dup[order(names(t_dup))]
id <- unique(names(t_dup))

t.dup.rm <- NULL
for(j in 1:length(id)){
	tmp <- t_dup[which(names(t_dup)%in%id[j])]
	tmp = mean(tmp, na.rm=T)
	t.dup.rm <- c(t.dup.rm,tmp)			
}
size <- c(t_uniq,t.dup.rm)
names(size) <- c(names(t_uniq),id)

data = data[names(size),]

##################################
## Compute TPM data

GetTPM <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

TPM = log2( GetTPM(data,size) + 1 )
colnames(TPM) = paste0(clin[colnames(TPM),]$ANONPT_ID)
TPM <- data.frame(TPM)
TPM <- cbind(gene_id=rownames(TPM), TPM)
gz <- gzfile(file.path(work_dir, 'EXPR.txt.gz'), "w")
write.table( TPM , file=gz , quote=FALSE , sep="\t" , col.names=TRUE , row.names=FALSE )
close(gz)

# write.table( TPM , file= "/data/Mariathasan/EXPR.txt" , quote=FALSE , sep=";" , col.names=TRUE , row.names=TRUE )
# system("gzip --force /data/Mariathasan/EXPR.txt > /data/Mariathasan/EXPR.txt.gz")
