# modified on 15 Mar 2017
# this script will produce a tiled grid colored by impact for a variant annotation saved in Gemini database with two phenotypes
# usage: Rscript genotype_plot.r database tsv_output phenotype_lable1 phenotype_label2 pdf_output
library(ggplot2)
library(reshape2)

args=commandArgs()
db=args[1]
tsv=args[2]
pheno1=args[3]
pheno2=args[4]
pdf=args[5]

# query gemini database -- gemini must be in $PATH
cmd=paste("gemini query -q 'select  v.gene,v.chrom,v.start,v.ref,v.alt,v.vep_hgvsp,v.impact_so, (gts).(phenotype==1),(gts).(phenotype==2)   from variants v where   --header ",db," > ", tsv, sep='.')
system(cmd)

# get the number of samples for the first phenotype to introduce a brek in plot
pheno1_count=as.integer(system(paste("gemini query -q 'select count(1) from samples where phenotype=1 ",db," intern=TRUE")))

# load the results, sort by gene and start, clean sample names from gemini-prefix gts
genotypes=read.table(tsv,header=TRUE,sep='\t') 
genotypes=genotypes[order(genotypes$gene,genotypes$start),]
colnames(genotypes)=gsub("gts",'',colnames(genotypes))

# introduce break between the phenotypes in plot by creating a empty column
genotypes <- insert.after.column.number(genotypes, pheno1_count, NA) 
colnames(genotypes)[pheno1_count+1]='   '

# filter out homozygous reference and NOCALL genotypes
gen_filter=sapply(genotypes,function(value) value!="NOCALL"&value!="HOMREF")

# reformat for column names
gen_filter=gen_filter[,7:] 
gen_filter=cbind(as.matrix(genotypes[,1]),paste(genotypes[,2],genotypes[,3],genotypes[,4],genotypes[,5],sep='.'),as.matrix(genotypes[,6]),gen_filter)
colnames(gen_filter)[1]="gene"
colnames(gen_filter)[2]="id"
colnames(gen_filter)[3]="impact"

# convert into data frame in order to melt (2 column format input for ggplot)
gen_df=as.data.frame(gen_filter) 
gen_df$id <- factor(gen_df$id, levels = gen_df$id)
melted=melt(gen_df,id=c("id","gene","impact"))

# rename no impact columns to reference for coloring gray
melted$impact=factor(ifelse(melted$value,as.character(melted$impact),"reference"))

#plot using colors based on impact
ggplot(melted, aes(x = variable, y = id, fill = impact,h=5)) + geom_tile(aes(width=0.8, height=0.8)) +theme(legend.position = "right", panel.background = element_blank()) + scale_fill_manual(values = c("frameshift"="red", "inframe_indel"="yellow","intron"="blue","missense"="green","non_coding_exon"="purple","reference"="grey","splice"="orange","start_lost"="pink","stop_gained"="brown","5_prime"="black")) +xlab(paste(pheno1,pheno2,sep="\t")+ylab("gene")+theme(axis.title.x = element_text(hjust=0.2))+theme(axis.text.x = element_text(angle = 90, vjust = 1,size=5),axis.text.y=element_blank())+facet_grid(gene ~ .,scale="free",space="free",switch="y")+theme(panel.margin=unit(0.1,"mm"),strip.text.y = element_text(angle=180,size=4),legend.key.size=unit(2,"mm"),legend.text=element_text(size=4),legend.title=element_text(size=5),axis.ticks=element_blank())

#save as pdf
ggsave(pdf)
