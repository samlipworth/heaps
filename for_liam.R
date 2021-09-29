library(tidyverse)
library(data.table)
library(micropan)

guuids<-read_tsv('./data/guuids') 
cg<-filter(guuids,grepl('chromo',guuid)) # cg = chromosomal guuids
pg<-filter(guuids,!grepl('chromo',guuid)) # pg = plasmid guuids
guuids$guuid<-str_replace_all(guuids$guuid,'_.*','')

spec<-read_tsv('./data/kleborate.tsv') %>% select(strain,species)

guuids<-left_join(guuids,spec,by=c("guuid"="strain")) 
guuids<-filter(guuids,grepl('Kleb',species) | grepl('Esch',species))

c<-fread('~/data/gene_presence_absence.Rtab') # this is just the Panaroo matrix
x<-names(c)
chromo<-x[grepl('chromo',x)]
chromo_pangenome<-dplyr::select(c,chromo)

#nb these are concatenated, i.e. all plasmids in one fasta file
plasmid_pangenome <- fread('~/data/gene_presence_absence_concat_plas.Rtab')


dates<-read.delim('~/data/guuids_dates',sep=' ',header=F)
names(dates)<-c('acc','cd','t','guuid')
dates$year<-year(dates$cd)
d2009<-filter(dates, year==2009) %>% filter(guuid %in% guuids$guuid) 
d2018<-filter(dates, year==2018) %>% filter(guuid %in% guuids$guuid)

x<-names(chromo_pangenome)
x<-str_replace_all(x,'_.*','')
names(chromo_pangenome)<-x

cpg2009<-dplyr::select(chromo_pangenome,d2009$guuid)
cpg2009<-cpg2009[rowSums(cpg2009)>0]
cpg2018<-dplyr::select(chromo_pangenome,d2018$guuid)
cpg2018<-cpg2018[rowSums(cpg2018)>0]
out_cpg_2018=NULL
for(i in 1:100){
  print(i)
  t<-micropan::heaps(t(cpg2018),n.perm = 100)
  out_cpg_2018=rbind(out_cpg_2018,t[2])
}

out_cpg_2009=NULL
for(i in 1:100){
  print(i)
  t<-micropan::heaps(t(cpg2009),n.perm = 100)
  out_cpg_2009=rbind(out_cpg_2009,t[2])
}


d2009<-filter(d2009,guuid %in% names(plasmid_pangenome)) # not all isolates have plasmids
d2018<-filter(d2018,guuid %in% names(plasmid_pangenome))
ppg2009<-dplyr::select(plasmid_pangenome,d2009$guuid)
ppg2009<-ppg2009[rowSums(ppg2009)>0]

ppg2018<-dplyr::select(plasmid_pangenome,d2018$guuid)
ppg2018<-ppg2018[rowSums(ppg2018)>0]


out_ppg_2018=NULL
for(i in 1:100){
  print(i)
  t<-micropan::heaps(t(ppg2018),n.perm = 100)
  out_ppg_2018=rbind(out_ppg_2018,t[2])
}

out_ppg_2009=NULL
for(i in 1:100){
  print(i)
  t<-micropan::heaps(t(ppg2009),n.perm = 100)
  out_ppg_2009=rbind(out_ppg_2009,t[2])
}
