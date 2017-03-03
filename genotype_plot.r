library(ggplot2)
library(reshape2)
genotypes=read.table("/home/ealves/workfiles/driver_genotypes.tsv",header=TRUE,sep='\t')
gen_filter=sapply(genotypes,function(value) value!="./.") 
 gen_filter=gen_filter[,7:81] 
row.names(gen_filter)=paste(genotypes[,2],genotypes[,3],genotypes[,4],genotypes[,5],sep='.')
melted=melt(as.matrix(gen_filter))
ggplot(melted, aes(x = Var2, y = Var1, fill = value, h=10)) + geom_tile() +
scale_fill_manual(values = c("white", "black")) +theme(legend.position = "none") +xlab("sample")+ylab("mutation")+theme(axis.text.x = element_text(angle = 90, hjust = 1,size=5),axis.text.y=element_text(size=5))
  
gen_agg=aggregate(. ~ gene+ chrom + start ,data=genotypes,sum) 
  
  egfr_table=cbind(rowSums(gen_filter[,1:18]==TRUE),rowSums(gen_filter[,29:161]))
  
  
  g1 <- ggplot(data = egfr_table, aes(x = row.names(egfr_table), y = EGFR_pos)) +
  geom_bar(stat = "identity") + ggtitle("EGFR pos") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(1,0,1,0), "mm")) +
  scale_y_reverse() + coord_flip()

g2 <- ggplot(data = DATA, aes(x = row.names(egfr_table), y = EGFR_neg)) +
  geom_bar(stat = "identity") + ggtitle("EGFR neg") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(1,5,1,0), "mm")) + coord_flip()

multiplot(g1, g2, cols = 2)
  coverage=read.csv("/home/ealves/datafiles/all_samples_coverage.csv",header=TRUE)
  plot(coverage$avg,log="y",pch=".",xlab="Amplicon",ylab="Coverage")
   plot(coverage$pct,pch=".",xlab="Amplicon",ylab="Percentage")
  gene_coverage=read.csv("/home/ealves/workfiles/gene_coverage.csv",header=TRUE)
  plot(gene_coverage$gt500/gene_coverage$tot,pch=".",xlab="Gene",ylab="Percent>500")
  coverage$gene=sapply(strsplit(as.character(coverage$amplicon),"_"),"[[",1)
  gen_agg=aggregate(. ~ gene ,data=coverage,sum) 
 
  genotypes=read.table("/home/ealves/workfiles/gefitinib_genotypes.tsv",header=TRUE,sep='\t')
gen_filter=sapply(genotypes,function(value) value!="NOCALL"&value!="HOMREF")
  gen_filter=gen_filter[,-(22:23)]  #removed 184ST AND RIT
  gen_filter=gen_filter[,-64] #189 ST
  gen_filter=gen_filter[,-24] #209 
 gen_filter=gen_filter[,7:81] 
row.names(gen_filter)=paste(genotypes[,2],genotypes[,3],genotypes[,4],genotypes[,5],sep='.')
  top_mutation=gen_filter[rowSums(gen_filter[,1:75]==TRUE)>3,]
  melted=melt(as.matrix(top_mutation))
ggplot(melted, aes(x = Var2, y = Var1, fill = value, h=10)) + geom_tile() +
scale_fill_manual(values = c("white", "black")) +theme(legend.position = "none") +xlab("sample")+ylab("mutation")+theme(axis.text.x = element_text(angle = 90, hjust = 1,size=5),axis.text.y=element_text(size=5))
  
   gef_table=cbind(rowSums(top_mutation[,1:28]==TRUE),rowSums(top_mutation[,19:75]))
  
  colnames(gef_table)=c("DisCon","NoDisCon")
  g1 <- ggplot(data = as.data.frame(gef_table), aes(x = row.names(top_mutation), y = DisCon)) +
  geom_bar(stat = "identity") + ggtitle("Dis Control") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(1,0,1,0), "mm")) +
  scale_y_reverse() + coord_flip()

g2 <- ggplot(data = as.data.frame(gef_table), aes(x = row.names(top_mutation), y = NoDisCon)) +
  geom_bar(stat = "identity") + ggtitle("No Dis Control") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y=element_text(size=3),plot.margin = unit(c(1,5,1,0), "mm")) + coord_flip()

multiplot(g1, g2, cols = 2)

  insert.after.column.number <- function (mat, cnum, what) {
  if (length(what) != nrow(mat)) what <- c(what, rep_len(NA, nrow(mat) - length(what)))
  cbind(mat[, 1:cnum], what, mat[,(cnum + 1):ncol(mat) ])
}

  
   genotypes=read.table("/home/eduardoalves/RDP_genotypes_gefit_r.tsv",header=TRUE,sep='\t')
genotypes$TC_184ST=NULL  #removed 184ST 
genotypes$TC_184RIT=NULL  #removed 184RIT 
genotypes$TC_189ST=NULL  #removed 189ST 
genotypes$TC209DT=NULL  #removed 209DT
gen_filter=sapply(genotypes,function(value) value!="NOCALL"&value!="HOMREF")
 gen_filter=gen_filter[,6:80] 
row.names(gen_filter)=paste(genotypes[,2],genotypes[,3],genotypes[,4],genotypes[,5],sep='.')
  gen_filter <- insert.after.column.number(gen_filter, 23, NA)
  colnames(gen_filter)=gsub("gts.006",'',colnames(gen_filter))
  colnames(gen_filter)[24]=''
  melted=melt(as.matrix(gen_filter))
ggplot(melted, aes(x = Var2, y = Var1, fill = value, h=5)) + geom_tile(aes(width=0.8, height=0.8)) +
scale_fill_manual(values = c("grey", "black")) +theme(legend.position = "none", panel.background = element_blank()) +xlab("sample")+ylab("mutation")+theme(axis.text.x = element_text(angle = 90, vjust = 1,size=5),axis.text.y=element_text(size=5))
  gen_order=gen_filter[order(rowSums(gen_filter[,-24])),]
  melted=melt(as.matrix(gen_order))
ggplot(melted, aes(x = Var2, y = Var1, fill = value, h=5)) + geom_tile(aes(width=0.8, height=0.8)) + scale_fill_manual(values = c("grey", "black")) +theme(legend.position = "none", panel.background = element_blank()) +xlab("disease control                          no disease control")+ylab("mutation")+theme(axis.title.x = element_text(hjust=0.2))+theme(axis.text.x = element_text(angle = 90, vjust = 1,size=5),axis.text.y=element_text(size=5))
  driver_by_gene=read.table("/home/ealves/workfiles/driver_summary_by_gene.tsv",header=TRUE,sep='\t')
ggsave("RDP_genotype.pdf") 
  
 #color by impact_so
  genotypes=read.table("/home/eduardoalves/gefitinib_genotypes_impact_r.tsv",header=TRUE,sep='\t') 
genotypes=genotypes[order(genotypes$gene,genotypes$start),]
colnames(genotypes)=gsub("gts.006",'',colnames(genotypes))
genotypes$TC_184ST=NULL  #removed 184ST 
genotypes$TC_184RIT=NULL  #removed 184RIT 
genotypes$TC_189ST=NULL  #removed 189ST 
genotypes$TC209DT=NULL  #removed 209DT
genotypes <- insert.after.column.number(genotypes, 29, NA) 
colnames(genotypes)[30]='   '
gen_filter=sapply(genotypes,function(value) value!="NOCALL"&value!="HOMREF")
gen_filter=gen_filter[,7:81] 
gen_filter=cbind(as.matrix(genotypes[,1]),paste(genotypes[,2],genotypes[,3],genotypes[,4],genotypes[,5],sep='.'),as.matrix(genotypes[,6]),gen_filter)
colnames(gen_filter)[1]="gene"
colnames(gen_filter)[2]="id"
colnames(gen_filter)[3]="impact"
gen_df=as.data.frame(gen_filter) 
gen_df$id <- factor(gen_df$id, levels = gen_df$id)
  
gen_df=gen_df[101:220,]


melted=melt(gen_df,id=c("id","gene","impact"))
melted$impact=factor(ifelse(melted$value,as.character(melted$impact),"reference"))
melted$impact=factor(ifelse(melted$impact=="non_coding_transcript_exon_variant","non_coding_exon",as.character(melted$impact)))
melted$impact=factor(ifelse(melted$impact=="intron_variant","intron",as.character(melted$impact)))
melted$impact=factor(ifelse(melted$impact=="missense_variant","missense",as.character(melted$impact)))
melted$impact=factor(ifelse(melted$impact=="frameshift_variant","frameshift",as.character(melted$impact)))
melted$impact=factor(ifelse(melted$impact=="splice_donor_variant","splice",as.character(melted$impact)))
melted$impact=factor(ifelse(melted$impact=="splice_acceptor_variant","splice",as.character(melted$impact)))
melted$impact=factor(ifelse(melted$impact=="splice_region_variant","splice",as.character(melted$impact)))
melted$impact=factor(ifelse(melted$impact=="inframe_deletion","inframe_indel",as.character(melted$impact)))
melted$impact=factor(ifelse(melted$impact=="inframe_insertion","inframe_indel",as.character(melted$impact)))
melted$impact=factor(ifelse(melted$impact=="5_prime_UTR_variant","5_prime",as.character(melted$impact)))

  ggplot(melted, aes(x = variable, y = id, fill = impact,h=5)) + geom_tile(aes(width=0.8, height=0.8)) +theme(legend.position = "right", panel.background = element_blank()) + scale_fill_manual(values = c("frameshift"="red", "inframe_indel"="yellow","intron"="blue","missense"="green","non_coding_exon"="purple","reference"="grey","splice"="orange","start_lost"="pink","stop_gained"="brown","5_prime"="black")) +xlab("disease control                          no disease control")+ylab("gene")+theme(axis.title.x = element_text(hjust=0.2))+theme(axis.text.x = element_text(angle = 90, vjust = 1,size=5),axis.text.y=element_blank())+facet_grid(gene ~ .,scale="free",space="free",switch="y")+theme(panel.margin=unit(0.1,"mm"),strip.text.y = element_text(angle=180,size=4),legend.key.size=unit(2,"mm"),legend.text=element_text(size=4),legend.title=element_text(size=5),axis.ticks=element_blank())
ggsave("genotype_gene2.pdf")

egfr pos

 genotypes=read.table("/home/eduardoalves/gefitinib_genotypes_pos_r.tsv",header=TRUE,sep='\t') 
genotypes=genotypes[order(genotypes$gene,genotypes$start),]
colnames(genotypes)=gsub("gts.006",'',colnames(genotypes))
genotypes <- insert.after.column.number(genotypes, 13, NA) 
colnames(genotypes)[14]='   '
gen_filter=sapply(genotypes,function(value) value!="NOCALL"&value!="HOMREF")
gen_filter=gen_filter[,7:25] 
gen_filter=cbind(as.matrix(genotypes[,1]),paste(genotypes[,2],genotypes[,3],genotypes[,4],genotypes[,5],sep='.'),as.matrix(genotypes[,6]),gen_filter)
colnames(gen_filter)[1]="gene"
colnames(gen_filter)[2]="id"
colnames(gen_filter)[3]="impact"
gen_df=as.data.frame(gen_filter) 
gen_df$id <- factor(gen_df$id, levels = gen_df$id)
melted=melt(gen_df,id=c("id","gene","impact"))
melted$impact=factor(ifelse(melted$value,as.character(melted$impact),"reference"))
melted$impact=factor(ifelse(melted$impact=="non_coding_transcript_exon_variant","non_coding_exon",as.character(melted$impact)))
melted$impact=factor(ifelse(melted$impact=="intron_variant","intron",as.character(melted$impact)))
melted$impact=factor(ifelse(melted$impact=="missense_variant","missense",as.character(melted$impact)))
melted$impact=factor(ifelse(melted$impact=="frameshift_variant","frameshift",as.character(melted$impact)))
melted$impact=factor(ifelse(melted$impact=="splice_donor_variant","splice",as.character(melted$impact)))
melted$impact=factor(ifelse(melted$impact=="splice_acceptor_variant","splice",as.character(melted$impact)))
melted$impact=factor(ifelse(melted$impact=="splice_region_variant","splice",as.character(melted$impact)))
melted$impact=factor(ifelse(melted$impact=="inframe_deletion","inframe_indel",as.character(melted$impact)))
melted$impact=factor(ifelse(melted$impact=="inframe_insertion","inframe_indel",as.character(melted$impact)))
melted$impact=factor(ifelse(melted$impact=="5_prime_UTR_variant","5_prime",as.character(melted$impact)))

  ggplot(melted, aes(x = variable, y = id, fill = impact,h=5)) + geom_tile(aes(width=0.8, height=0.8)) +theme(legend.position = "right", panel.background = element_blank()) + scale_fill_manual(values = c("frameshift"="red", "inframe_indel"="yellow","intron"="blue","missense"="green","non_coding_exon"="purple","reference"="grey","splice"="orange","start_lost"="pink","stop_gained"="brown","5_prime"="black")) +xlab("disease control                          no disease control")+ylab("gene")+theme(axis.title.x = element_text(hjust=0.2))+theme(axis.text.x = element_text(angle = 90, vjust = 1,size=5),axis.text.y=element_blank())+facet_grid(gene ~ .,scale="free",space="free",switch="y")+theme(panel.margin=unit(0.1,"mm"),strip.text.y = element_text(angle=180,size=4),legend.key.size=unit(2,"mm"),legend.text=element_text(size=4),legend.title=element_text(size=5),axis.ticks=element_blank())
ggsave("genotype_gene_pos.pdf")

4 categories

 genotypes=read.table("/home/eduardoalves/gefitinib_genotypes_cat_r.tsv",header=TRUE,sep='\t') 
genotypes=genotypes[order(genotypes$gene,genotypes$start),]
colnames(genotypes)=gsub("gts.006",'',colnames(genotypes))
genotypes$TC_184ST=NULL  #removed 184ST 
genotypes$TC_184RIT=NULL  #removed 184RIT 
genotypes$TC_189ST=NULL  #removed 189ST 
genotypes$TC209DT=NULL  #removed 209DT
genotypes <- insert.after.column.number(genotypes, 40, NA)
genotypes <- insert.after.column.number(genotypes, 24, NA) 
genotypes <- insert.after.column.number(genotypes, 13, NA)  
colnames(genotypes)[43]=' '
colnames(genotypes)[26]='  '
colnames(genotypes)[14]='   '
gen_filter=sapply(genotypes,function(value) value!="NOCALL"&value!="HOMREF")
gen_filter=gen_filter[,7:81] 
gen_filter=cbind(as.matrix(genotypes[,1]),paste(genotypes[,2],genotypes[,3],genotypes[,4],genotypes[,5],sep='.'),as.matrix(genotypes[,6]),gen_filter)
colnames(gen_filter)[1]="gene"
colnames(gen_filter)[2]="id"
colnames(gen_filter)[3]="impact"
gen_df=as.data.frame(gen_filter) 
gen_df$id <- factor(gen_df$id, levels = gen_df$id)
  
gen_df=gen_df[101:220,]


melted=melt(gen_df,id=c("id","gene","impact"))
melted$impact=factor(ifelse(melted$value,as.character(melted$impact),"reference"))
melted$impact=factor(ifelse(melted$impact=="non_coding_transcript_exon_variant","non_coding_exon",as.character(melted$impact)))
melted$impact=factor(ifelse(melted$impact=="intron_variant","intron",as.character(melted$impact)))
melted$impact=factor(ifelse(melted$impact=="missense_variant","missense",as.character(melted$impact)))
melted$impact=factor(ifelse(melted$impact=="frameshift_variant","frameshift",as.character(melted$impact)))
melted$impact=factor(ifelse(melted$impact=="splice_donor_variant","splice",as.character(melted$impact)))
melted$impact=factor(ifelse(melted$impact=="splice_acceptor_variant","splice",as.character(melted$impact)))
melted$impact=factor(ifelse(melted$impact=="splice_region_variant","splice",as.character(melted$impact)))
melted$impact=factor(ifelse(melted$impact=="inframe_deletion","inframe_indel",as.character(melted$impact)))
melted$impact=factor(ifelse(melted$impact=="inframe_insertion","inframe_indel",as.character(melted$impact)))
melted$impact=factor(ifelse(melted$impact=="5_prime_UTR_variant","5_prime",as.character(melted$impact)))

  ggplot(melted, aes(x = variable, y = id, fill = impact,h=5)) + geom_tile(aes(width=0.8, height=0.8)) +theme(legend.position = "right", panel.background = element_blank()) + scale_fill_manual(values = c("frameshift"="red", "inframe_indel"="yellow","intron"="blue","missense"="green","non_coding_exon"="purple","reference"="grey","splice"="orange","start_lost"="pink","stop_gained"="brown","5_prime"="black")) +xlab("dc    no dc              dc                           no dc     \nEGFR+                        EGFR-")+ylab("gene")+theme(axis.title.x = element_text(hjust=0.2))+theme(axis.text.x = element_text(angle = 90, vjust = 1,size=5),axis.text.y=element_blank())+facet_grid(gene ~ .,scale="free",space="free",switch="y")+theme(panel.margin=unit(0.1,"mm"),strip.text.y = element_text(angle=180,size=4),legend.key.size=unit(2,"mm"),legend.text=element_text(size=4),legend.title=element_text(size=5),axis.ticks=element_blank())
ggsave("genotype_cat.pdf")



frequencies per gene ..
genotypes=read.table("/home/eduardoalves/driver_genotypes_r.tsv",header=TRUE,sep='\t') 
colnames(genotypes)=gsub("gts.006",'',colnames(genotypes))
genotypes$TC_184ST=NULL  #removed 184ST 
genotypes$TC_184RIT=NULL  #removed 184RIT 
genotypes$TC_189ST=NULL  #removed 189ST 
genotypes$TC209DT=NULL  #removed 209DT
genotypes$TC280ACTM=NULL 
genotypes$TC_185ST=NULL 
genotypes$TC_185RIM=NULL 
genotypes$TC_186ST=NULL
genotypes$TC_182_ACT=NULL
gen_filter=sapply(genotypes,function(value) value!="NOCALL"&value!="HOMREF")
gen_filter=cbind(as.matrix(genotypes[,1]),genfilter[,6:157])
colnames(gen_filter)[1]="gene"
gen_agg=aggregate(. ~ gene,gen_filter,sum) 
gen_freq=data.frame(as.matrix(gen_agg[,1]),rowSums(gen_agg[,2:153]>0))
gen_freq=gen_freq[order(gen_freq[,2]),]
gen_agg=gen_agg[order(gen_freq[,2]),]
gen_agg=as.data.frame(cbind(as.matrix(gen_agg[,1]),gen_agg[,2:153]>0))
gen_agg$gene <- factor(gen_agg$gene, levels = gen_freq[,1])
colnames(gen_agg)[1]="gene"
melted=melt(gen_agg[-(1,39),],id="gene")
ggplot(melted, aes(x = variable, y = gene, fill = value,h=5)) + geom_tile(aes(width=0.8, height=0.8)) +scale_fill_manual(values = c("grey", "black")) +theme(legend.position = "none", panel.background = element_blank()) +xlab("samples")+ylab("gene")+theme(axis.title.x = element_text(hjust=0.2))+theme(axis.text.x = element_text(angle = 90, vjust = 1,size=3),axis.text.y=element_text(size=3))
ggsave("driver_cohort.pdf")
gen_freq=as.data.frame(gen_freq)
gen_freq$gene <- factor(gen_freq$gene, levels = gen_freq$gene)
colnames(gen_freq)[1]="gene"
ggplot(data = gen_freq[-(1:39),], aes(x = gene, y = freq)) +geom_bar(stat = "identity") + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y=element_text(size=3),plot.margin = unit(c(1,5,1,0), "mm")) + coord_flip()

calculate frequecies for the 4 groups

 genotypes=read.table("/home/eduardoalves/gefitinib_genotypes_cat_r.tsv",header=TRUE,sep='\t') 
genotypes=genotypes[order(genotypes$gene,genotypes$start),]
colnames(genotypes)=gsub("gts.006",'',colnames(genotypes))
genotypes$TC_184ST=NULL  #removed 184ST 
genotypes$TC_184RIT=NULL  #removed 184RIT 
genotypes$TC_189ST=NULL  #removed 189ST 
genotypes$TC209DT=NULL  #removed 209DT
gen_filter=as.data.matrix(sapply(genotypes,function(value) value!="NOCALL"&value!="HOMREF"))
gen_filter=data.frame(genotypes[,1],gen_filter[,7:81])
colnames(gen_filter)[1]="gene"
gen_agg=aggregate(. ~ gene,gen_filter,sum) 
gen_freq=data.frame(as.matrix(gen_agg[,1]),rowSums(gen_agg[,2:8]>0),rowSums(gen_agg[,9:19]>0),rowSums(gen_agg[,20:35]>0),rowSums(gen_agg[,36:76]>0))

