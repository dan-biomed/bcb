# R assignment
library(dplyr)
library(tidyverse)
library(stringr)
library(reshape2)
##this is delimited by tab
snp<-read.table("/Users/dan/Desktop/BCB546X-Spring2017/UNIX_Assignment/snp_position.txt", header = T,
sep="\t")
genotype<-read.table("/Users/dan/Desktop/BCB546X-Spring2017/UNIX_Assignment/fang_et_al_genotypes.txt",
sep="\t")
dim(snp)
dim(genotype)
##take the columns of interest, ie the SNP_ID, Chromosome, and position columns.
##same output, comma is to be clear (snp_cut_col<-snp[c(1,3,4)] snp_cut_col<-snp[,c(1,3,4)]
##alternatively, I can use:
snp_cut_col<-as.data.frame(snp[c(1,3,4)])
##taking out all the maize varieties and teosinte varieties into separate objects.
maize<-genotype %>% filter(V3=="ZMMIL"|V3=="ZMMLR"|V3=="ZMMMR")
teosinte<-genotype %>% filter(V3=="ZMPBA"|V3=="ZMPIL"|V3=="ZMPJA")
##since the new varieties object has no header, add
#header (SNP_ID column from genotype) back in as header so merging is possible later using this as common column.
id<-genotype[1,]
maize_header<-rbind(id,maize)
teosinte_header<-rbind(id, teosinte)
##then transpose maize and teosinte's rows for columns, into objects t_maize and t_teosinte.
t_maize_header<-as.data.frame(t(maize_header))
t_teosinte_header <- as.data.frame(t(teosinte_header))
##check if there are repeating row names
length(unique(t_maize_header))
length(t_maize_header)
length(unique(t_teosinte_header))
length(t_teosinte_header)
##next join the snp_cut_col to transposed maize or teosinte genotypes using SNP_ID as common column
colnames(t_maize_header)[1]<-"SNP_ID"
maize_join<-full_join(snp_cut_col,t_maize_header,by="SNP_ID")
colnames(t_teosinte_header)[1]<-"SNP_ID"
teosinte_join <- full_join(snp_cut_col,t_teosinte_header,by="SNP_ID")
##replace ?/? with ?? in the joined maize file, into new object maize_join_sub
maize_join_sub<-as.data.frame(lapply(maize_join, function(x) {
gsub("\\?\\/\\?", "\\?\\?",x)
}))

##for maize need 20 files, 10 for increasing order of chromosomes with ??
maize_join_sub_in<-maize_join_sub %>% arrange(Chromosome)

##make the first 10 files of maize with each chromosome per file.
maize_chrom01_incre<-subset(maize_join_sub_in, Chromosome=="1")
maize_chrom02_incre<-subset(maize_join_sub_in, Chromosome=="2")
maize_chrom03_incre<-subset(maize_join_sub_in, Chromosome=="3")
maize_chrom04_incre<-subset(maize_join_sub_in, Chromosome=="4")
maize_chrom05_incre<-subset(maize_join_sub_in, Chromosome=="5")
maize_chrom06_incre<-subset(maize_join_sub_in, Chromosome=="6")
maize_chrom07_incre<-subset(maize_join_sub_in, Chromosome=="7")
maize_chrom08_incre<-subset(maize_join_sub_in, Chromosome=="8")
maize_chrom09_incre<-subset(maize_join_sub_in, Chromosome=="9")
maize_chrom10_incre<-subset(maize_join_sub_in, Chromosome=="10")
##then replace the ?? in joined maize file with --,
maize_join_sub_de<-as.data.frame(lapply(maize_join, function(x) {
gsub("\\?\\/\\?", "\\-\\-",x)
}))
##rearrange the same object, in decreasing chromosome order, into new object maize_merge_sub, and replace ?? with --
maize_join_sub_de<-maize_join_sub %>% arrange(desc(Chromosome)) %>%
lapply(function(x) {
gsub("\\?\\?", "\\-\\-",x)
})
##the above will make a list, so then make this into dataframe.
maize_join_sub_de<-as.data.frame(maize_join_sub_de)
##make the second 10 files of maize with decreasing chromosome order, 1 chromosome per file
maize_chrom01_decre<-subset(maize_join_sub_de, Chromosome=="1")
maize_chrom02_decre<-subset(maize_join_sub_de, Chromosome=="2")
maize_chrom03_decre<-subset(maize_join_sub_de, Chromosome=="3")
maize_chrom04_decre<-subset(maize_join_sub_de, Chromosome=="4")
maize_chrom05_decre<-subset(maize_join_sub_de, Chromosome=="5")
maize_chrom06_decre<-subset(maize_join_sub_de, Chromosome=="6")
maize_chrom07_decre<-subset(maize_join_sub_de, Chromosome=="7")
maize_chrom08_decre<-subset(maize_join_sub_de, Chromosome=="8")
maize_chrom09_decre<-subset(maize_join_sub_de, Chromosome=="9")
maize_chrom10_decre<-subset(maize_join_sub_de, Chromosome=="10")
##for teosinte part
##replace ?/? with ?? in the joined teosinte file, into new object teosinte_join_sub
teosinte_join_sub<-as.data.frame(lapply(teosinte_join, function(x) {
gsub("\\?\\/\\?", "\\?\\?",x)
}))
##for teosinte need 20 files, 10 for increasing order of chromosomes
teosinte_join_sub_in<-teosinte_join_sub %>% arrange(Chromosome)
##make the first 10 files of teosinte with each chromosome per file.
teosinte_chrom01_incre<-subset(teosinte_join_sub_in, Chromosome=="1")
teosinte_chrom02_incre<-subset(teosinte_join_sub_in, Chromosome=="2")
teosinte_chrom03_incre<-subset(teosinte_join_sub_in, Chromosome=="3")
teosinte_chrom04_incre<-subset(teosinte_join_sub_in, Chromosome=="4")
teosinte_chrom05_incre<-subset(teosinte_join_sub_in, Chromosome=="5")
teosinte_chrom06_incre<-subset(teosinte_join_sub_in, Chromosome=="6")
teosinte_chrom07_incre<-subset(teosinte_join_sub_in, Chromosome=="7")
teosinte_chrom08_incre<-subset(teosinte_join_sub_in, Chromosome=="8")
teosinte_chrom09_incre<-subset(teosinte_join_sub_in, Chromosome=="9")
teosinte_chrom10_incre<-subset(teosinte_join_sub_in, Chromosome=="10")
##then replace the ?? in joined teosinte file with --,
teosinte_join_sub_de<-as.data.frame(lapply(teosinte_join, function(x) {
gsub("\\?\\/\\?", "\\-\\-",x)
}))
##rearrange the same object, in decreasing chromosome order, into new object teosinte_merge_sub, and replace ?? with --
teosinte_join_sub_de<-teosinte_join_sub %>% arrange(desc(Chromosome)) %>%
lapply(function(x) {
gsub("\\?\\?", "\\-\\-",x)
})
teosinte_join_sub_de<-as.data.frame(teosinte_join_sub_de)
##make the second 10 files of maize with decreasing chromosome order, 1 chromosome per file
teosinte_chrom01_decre<-subset(teosinte_join_sub_de, Chromosome=="1")
teosinte_chrom02_decre<-subset(teosinte_join_sub_de, Chromosome=="2")
teosinte_chrom03_decre<-subset(teosinte_join_sub_de, Chromosome=="3")
teosinte_chrom04_decre<-subset(teosinte_join_sub_de, Chromosome=="4")
teosinte_chrom05_decre<-subset(teosinte_join_sub_de, Chromosome=="5")
teosinte_chrom06_decre<-subset(teosinte_join_sub_de, Chromosome=="6")
teosinte_chrom07_decre<-subset(teosinte_join_sub_de, Chromosome=="7")
teosinte_chrom08_decre<-subset(teosinte_join_sub_de, Chromosome=="8")
teosinte_chrom09_decre<-subset(teosinte_join_sub_de, Chromosome=="9")
teosinte_chrom10_decre<-subset(teosinte_join_sub_de, Chromosome=="10")





##below is to make data tidy
##need to change headers and add in the SNP, chromosome, and position headers
colnames(maize_join_sub_in)<-as.character(unlist(maize_join_sub_in[984,]))
colnames(maize_join_sub_in)[1]<-"SNP_ID"
colnames(maize_join_sub_in)[2]<-"Chromosome"
colnames(maize_join_sub_in)[3]<-"Position"
##cut out the troublesome rows
##maize_join_sub_in<-maize_join_sub_in[-(984:986),]

maize_group<-maize_join_sub_in
colnames(maize_group)<-as.character(unlist(maize_group[986,]))
colnames(maize_group)[1]<-"SNP_ID"
colnames(maize_group)[2]<-"Chromosome"
colnames(maize_group)[3]<-"Position"

dim(maize_join_sub_in)
tidied_maize<-maize_join_sub_in %>% gather(key=Sample_ID, value=genotypes, 4:1576)
dim(maize_group)
tidied_group<-melt(maize_group, id=c("SNP_ID","Chromosome","Position"))
colnames(tidied_group)[4]<-"Group"
colnames(tidied_group)[5]<-"genotypes"

add_group<-full_join(tidied_maize, tidied_group, by=c("SNP_ID","Position","Chromosome","genotypes"))
##drop the last 3 rows
add_group<-add_group[-(3050451:3050453),]

ggplot(add_group)+geom_bar(aes(x=Chromosome))


##then if homozygous, then assign "homozygous" to the column homohetero
add_group<-within(add_group, {
  homohetero<-NA
  homohetero[genotypes=="A/A"|genotypes=="T/T"|genotypes=="C/C"|genotypes=="G/G"]<-"homozygous"
  homohetero[genotypes=="A/T"|genotypes=="A/C"|genotypes=="A/G"|genotypes=="C/T"|genotypes=="C/G"|genotypes=="C/A"|
               genotypes=="T/A"|genotypes=="T/C"|genotypes=="T/G"|genotypes=="G/A"|genotypes=="G/C"|genotypes=="G/T"]<-"heterozygous"
})

add_group<-as.data.frame(lapply(add_group, function(x) {
  gsub("\\?\\?", NA,x)
}))

add_group<-add_group %>% arrange(Sample_ID, Group)


ggplot(add_group,aes(x=Sample_ID, fill=homohetero))+
  geom_bar(position="fill")+
  facet_wrap(~Group)

ggplot(add_group)+geom_bar(aes(x=homohetero))+facet_wrap(~Group)








# UNIX ASSIGNMENT
## Data Inspection

##First check for structures and dimensions of the 2 files.

$ awk -F "\t" '{print NF; exit}' fang_et_al_genotypes.txt				### inspect number of columns in fang_et_al_genotypes.txt, which there are 986

$ awk -F "\t" '{print NF; exit}' snp_position.txt					### inspect number of columns in snp_position.txt, which there are 15



## Data Processing


$ grep -v "^#" snp_position.txt | cut -f 1,3,4 > snp_positions_cut_columns.txt		### since only interested in SNP_ID, Chromosome, and position, cut these columns to new file to be joined

$ grep -v "^#" fang_et_al_genotypes.txt | cut -f 3-986 > fang_cut_columns.txt		### From the assignment it wasn't clear if Sample ID needed to be kept, but this line will allow the files (after trasposition) to line up on column 1, with matching SNP_ID
 

### I could transpose at this point, but it seems to me that once I transpose the fang_cut_columns.txt I'd have a hard time extracting the teosinte and maize, so I'll extract these two first into two files.

$ grep Group fang_cut_columns.txt > header_only.txt					### I know that "Group" is in the header, and I want to keep the header, so I grep out this line to combine with other files later. Alternatively, I could have used $ cut -f 2 fang_cut_columns.txt | head -n 2

$ grep ZMMIL fang_cut_columns.txt > ZMMIL.txt						### Since the assignment asked for different orders of the varieties, and I don't know where the row started for each variety, I grep the three kinds of maize (and also for teosinte) and output them to new files, then cat them together
$ grep ZMMLR fang_cut_columns.txt > ZMMLR.txt
$ grep ZMMMR fang_cut_columns.txt > ZMMMR.txt

$ cat header_only.txt ZMMIL.txt ZMMLR.txt ZMMMR.txt >maize.txt

$ grep ZMPBA fang_cut_columns.txt > ZMPBA.txt
$ grep ZMPIL fang_cut_columns.txt > ZMPIL.txt
$ grep ZMPJA fang_cut_columns.txt > ZMPJA.txt

$ cat header_only.txt zmpba.txt zmpil.txt zmpja.txt > teosinte.txt


$ awk -f transpose.awk maize.txt > transposed_maize.txt					### The next step is to transpose the maize and teosinte files
$ awk -f transpose.awk teosinte.txt > transposed_teosinte.txt

$ sort -k1,1 snp_positions_cut_columns.txt >sorted_snp_positions_cut_columns.txt	### sort the files to be joined
$ sort -k1,1 transposed_maize.txt >sorted_transposed_maize.txt

											### Next step is to join the transposed maize or teosinte to the snp position file	
$ join --header -t $'\t' -1 1 -2 1 sorted_snp_positions_cut_columns.txt sorted_transposed_maize.txt > joined_maize.txt
$ head -n1 joined_maize.txt								
$ join --header -t $'\t' -1 1 -2 1 sorted_snp_positions_cut_columns.txt sorted_transposed_teosinte.txt > joined_teosinte.txt
$ head -n1 joined_maize.txt								### checked the order of a couple data points, seems correct. 


$ sort -k2n,2 -k1,1 joined_maize.txt > sorted_maize.txt					### sorted the joined_maize.txt by alphanumeric order of chromosome (which is second column) first (-k2n,2) and second column is both start and end key. Then this is sorted by first column alphabetically.
$ sort -k2n,2 -k1,1 joined_teosinte.txt > sorted_teosinte.txt				### sorted joined_teosinte.txt

$ sort -k2nr,2 -k1,1 joined_maize.txt > reverse_sorted_maize.txt			### reverse sorted joined_maize.txt
$ sort -k2nr,2 -k1,1 joined_teosinte.txt > reverse_sorted_teosinte.txt			### reverse sorted joined_teosinte.txt


$ sed "s^?/?^??^g" sorted_maize.txt > sed_maize.txt					### replace all ?/? with ?? in sorted_maize.txt
$ sed "s^?/?^??^g" sorted_teosinte.txt > sed_teosinte.txt				### replace all ?/? with ?? in sorted_teosinte.txt

$ sed "s^?/?^--^g" reverse_sorted_maize.txt > sed_reverse_maize.txt			### replace all ?/? with -- in reverse_sorted_maize.txt
$ sed "s^?/?^--^g" reverse_sorted_teosinte.txt > sed_reverse_teosinte.txt		### replace all ?/? with -- in reverse_sorted_teosinte.txt


####grep out each chromosomes in the four sed files. 

awk '$2 == 1' sed_maize.txt > 01maize_chromosome01.txt			
awk '$2 == 2' sed_maize.txt > 02maize_chromosome02.txt
awk '$2 == 3' sed_maize.txt > 03maize_chromosome03.txt
awk '$2 == 4' sed_maize.txt > 04maize_chromosome04.txt
awk '$2 == 5' sed_maize.txt > 05maize_chromosome05.txt
awk '$2 == 6' sed_maize.txt > 06maize_chromosome06.txt
awk '$2 == 7' sed_maize.txt > 07maize_chromosome07.txt
awk '$2 == 8' sed_maize.txt > 08maize_chromosome08.txt
awk '$2 == 9' sed_maize.txt > 09maize_chromosome09.txt
awk '$2 == 10' sed_maize.txt > 10maize_chromosome10.txt
awk '$2 == 1' sed_reverse_maize.txt > 11maize_chromosome01decreasing.txt			
awk '$2 == 2' sed_reverse_maize.txt > 12maize_chromosome02decreasing.txt	
awk '$2 == 3' sed_reverse_maize.txt > 13maize_chromosome03decreasing.txt	
awk '$2 == 4' sed_reverse_maize.txt > 14maize_chromosome04decreasing.txt	
awk '$2 == 5' sed_reverse_maize.txt > 15maize_chromosome05decreasing.txt	
awk '$2 == 6' sed_reverse_maize.txt > 16maize_chromosome06decreasing.txt	
awk '$2 == 7' sed_reverse_maize.txt > 17maize_chromosome07decreasing.txt	
awk '$2 == 8' sed_reverse_maize.txt > 18maize_chromosome08decreasing.txt	
awk '$2 == 9' sed_reverse_maize.txt > 19maize_chromosome09decreasing.txt	
awk '$2 == 10' sed_reverse_maize.txt > 20maize_chromosome10decreasing.txt	
awk '$2 == 1' sed_teosinte.txt > 21teosinte_chromosome01.txt			
awk '$2 == 2' sed_teosinte.txt > 22teosinte_chromosome02.txt
awk '$2 == 3' sed_teosinte.txt > 23teosinte_chromosome03.txt
awk '$2 == 4' sed_teosinte.txt > 24teosinte_chromosome04.txt
awk '$2 == 5' sed_teosinte.txt > 25teosinte_chromosome05.txt
awk '$2 == 6' sed_teosinte.txt > 26teosinte_chromosome06.txt
awk '$2 == 7' sed_teosinte.txt > 27teosinte_chromosome07.txt
awk '$2 == 8' sed_teosinte.txt > 28teosinte_chromosome08.txt
awk '$2 == 9' sed_teosinte.txt > 29teosinte_chromosome09.txt
awk '$2 == 10' sed_teosinte.txt > 30teosinte_chromosome10.txt
awk '$2 == 1' sed_reverse_teosinte.txt > 31maize_chromosome01decreasing.txt			
awk '$2 == 2' sed_reverse_teosinte.txt > 32maize_chromosome02decreasing.txt	
awk '$2 == 3' sed_reverse_teosinte.txt > 33maize_chromosome03decreasing.txt	
awk '$2 == 4' sed_reverse_teosinte.txt > 34maize_chromosome04decreasing.txt	
awk '$2 == 5' sed_reverse_teosinte.txt > 35maize_chromosome05decreasing.txt	
awk '$2 == 6' sed_reverse_teosinte.txt > 36maize_chromosome06decreasing.txt	
awk '$2 == 7' sed_reverse_teosinte.txt > 37maize_chromosome07decreasing.txt	
awk '$2 == 8' sed_reverse_teosinte.txt > 38maize_chromosome08decreasing.txt	
awk '$2 == 9' sed_reverse_teosinte.txt > 39maize_chromosome09decreasing.txt	
awk '$2 == 10' sed_reverse_teosinte.txt > 40maize_chromosome10decreasing.txt


