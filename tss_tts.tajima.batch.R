#!/usr/bin/env Rscript
## Script to read Tajima.D files generated from SLiM simulations/vcftools
## and produces tables for SNP counts and D values along binned intervals

## Install from github: "greymonroe/polymorphology"
library(data.table)
library(polymorphology)
library(R.utils)
library(dplyr)

dirpath = "/home/fli21/cassava-mutation-bias/"
# dirpath = "~/Documents/Monroe/Cassava Mutation Bias/"

## Read GFF file
cas_gff<-fread(paste0(dirpath, "Mesculenta_305_v6.1.gene_exons.gff3"), skip=3)
setnames(cas_gff,c("V1","V3","V4","V5"),c("chr","type","start","stop"))

## Read all Tajima D files in directory
#args = commandArgs(trailingOnly=TRUE)
# args = paste0(dirpath, "cassava_dm0.5_gds0.3_ids0.1_sf0.3_i10.Tajima.D")
args = Sys.glob(paste0(dirpath,"/tajima/*.Tajima.D"))

taj_tss_all<-data.table()

# aggregate all Tajima's D files
for(file in args){
  data<-fread(file)
  data$CHROM<-as.character(data$CHROM)
  ## Change Tajima file to have chromosome column consistent with GFF file
  data$CHROM <- "Chromosome01"
  taj_tss<-tss_tts.tajima(gff=cas_gff[chr=="Chromosome01"], tajima=data)
  taj_tss$dm=as.numeric(gsub(".+dm(.+)_gds.+","\\1", file))
  taj_tss$gds=as.numeric(gsub(".+gds(.+)_ids.+_i.+","\\1", file))
  taj_tss$ids=as.numeric(gsub(".+ids(.+)_sf.+","\\1", file))
  taj_tss$sf=as.numeric(gsub(".+sf(.+)_i.+","\\1", file))
  taj_tss$i=gsub(".+i(.+).Tajima.D","\\1", file)
  taj_tss$file=gsub(dirpath,"", args[1])
  taj_tss$trt<-paste0("dm", taj_tss$dm, "_gds", taj_tss$gds, "_ids", taj_tss$ids, "_sf", taj_tss$sf )
  ## Rowbind replicate to aggregate dataframe, bottleneck line
  taj_tss_all<-rbind(taj_tss, taj_tss_all)
}


## data.table that filters NA values, aggregates mean on D and N_SNPS, bins of 10 (6000 / 600)
trt_dt<-taj_tss_all[, .(D=mean(D, na.rm=T), NSNPS=mean(N_SNPS, na.rm=T)), by=.(bins=as.numeric(cut(pos, 600)), loc, trt, dm, gds, ids, sf)]
## Write this datatable to file so don't have to aggregate in the future
fwrite(trt_dt, paste0(dirpath,args[1],".trt",sep=""))

# create plots
pdf(paste0(dirpath, "/tajimasD.pdf"),width=4, height=2)
ggplot(trt_dt %>% filter(dm == 0.5, gds == 0, ids == 0), aes(x=bins, y=D, group=loc, col=loc))+
  geom_point(size=0.2)+
  geom_line(size = 0.1)+
  facet_grid(cols = vars(loc), rows = vars(sf))+
  theme_bw(base_size = 6)+
  # theme(legend.position = c(0.95,0.1),
  #       legend.key.size = unit(1, 'mm'))+
  scale_color_manual(values=c("orange","blue"))+
  ggtitle('Mutation rate = 0.5 w/o selective bias')
ggplot(trt_dt %>% filter(dm == 1.0, gds == 0, ids == 0), aes(x=bins, y=D, group=loc, col=loc))+
  geom_point(size=0.2)+
  geom_line(size = 0.1)+
  facet_grid(cols = vars(loc), rows = vars(sf))+
  theme_bw(base_size = 6)+
  scale_color_manual(values=c("orange","blue"))+
  ggtitle('Mutation rate = 1.0 w/o selective bias')
ggplot(trt_dt %>% filter(dm == 0.5, gds != 0, ids != 0), aes(x=bins, y=D, group=loc, col=loc))+
  geom_point(size=0.2)+
  geom_line(size = 0.1)+
  facet_grid(cols = vars(loc), rows = vars(sf))+
  theme_bw(base_size = 6)+
  scale_color_manual(values=c("orange","blue"))+
  ggtitle('Mutation rate = 0.5 w/ genic selection = 0.3 & intergenic selection = 0.1')
ggplot(trt_dt %>% filter(dm == 1.0, gds != 0, ids != 0), aes(x=bins, y=D, group=loc, col=loc))+
  geom_point(size=0.2)+
  geom_line(size = 0.1)+
  facet_grid(cols = vars(loc), rows = vars(sf))+
  theme_bw(base_size = 6)+
  scale_color_manual(values=c("orange","blue"))+
  ggtitle('Mutation rate = 1.0 w/ genic selection = 0.3 & intergenic selection = 0.1')
dev.off()

