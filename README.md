# BCB546X UNIX ASSIGNMENT
# Data Inspection

##First check for structures and dimensions of the 2 files.

$ awk -F "\t" '{print NF; exit}' fang_et_al_genotypes.txt				#### inspect number of columns in fang_et_al_genotypes.txt, which there are 986
986

$ awk -F "\t" '{print NF; exit}' snp_position.txt					#### inspect number of columns in snp_position.txt, which there are 15
15



# Data Processing


$ grep -v "^#" snp_position.txt | cut -f 1,3,4 > snp_positions_cut_columns.txt		#### since only interested in SNP_ID, Chromosome, and position, cut these columns to new file to be joined

$ grep -v "^#" fang_et_al_genotypes.txt | cut -f 3-986 > fang_cut_columns.txt		#### From the assignment it wasn't clear if Sample ID needed to be kept, but this line will allow the files (after trasposition) to line up on column 1, with matching SNP_ID
 

											#### I could transpose at this point, but it seems to me that once I transpose the fang_cut_columns.txt I'd have a hard time extracting the teosinte and maize, so I'll extract these two first into two files.

$ grep Group fang_cut_columns.txt > header_only.txt					#### I know that "Group" is in the header, and I want to keep the header, so I grep out this line to combine with other files later. Alternatively, I could have used $ cut -f 2 fang_cut_columns.txt | head -n 2

$ grep ZMMIL fang_cut_columns.txt > ZMMIL.txt						#### Since the assignment asked for different orders of the varieties, and I don't know where the row started for each variety, I grep the three kinds of maize (and also for teosinte) and output them to new files, then cat them together
$ grep ZMMLR fang_cut_columns.txt > ZMMLR.txt
$ grep ZMMMR fang_cut_columns.txt > ZMMMR.txt

$ cat header_only.txt ZMMIL.txt ZMMLR.txt ZMMMR.txt >maize.txt

$ grep ZMPBA fang_cut_columns.txt > ZMPBA.txt
$ grep ZMPIL fang_cut_columns.txt > ZMPIL.txt
$ grep ZMPJA fang_cut_columns.txt > ZMPJA.txt

$ cat header_only.txt zmpba.txt zmpil.txt zmpja.txt > teosinte.txt


$ awk -f transpose.awk maize.txt > transposed_maize.txt					#### The next step is to transpose the maize and teosinte files
$ awk -f transpose.awk teosinte.txt > transposed_teosinte.txt

$ sort -k1,1 snp_positions_cut_columns.txt >sorted_snp_positions_cut_columns.txt	#### sort the files to be joined
$ sort -k1,1 transposed_maize.txt >sorted_transposed_maize.txt

											#### Next step is to join the transposed maize or teosinte to the snp position file	
$ join --header -t $'\t' -1 1 -2 1 sorted_snp_positions_cut_columns.txt sorted_transposed_maize.txt > joined_maize.txt
$ head -n1 joined_maize.txt								
$ join --header -t $'\t' -1 1 -2 1 sorted_snp_positions_cut_columns.txt sorted_transposed_teosinte.txt > joined_teosinte.txt
$ head -n1 joined_maize.txt								#### checked the order of a couple data points, seems correct. 


$ sort -k2n,2 -k1,1 joined_maize.txt > sorted_maize.txt					#### sorted the joined_maize.txt by alphanumeric order of chromosome (which is second column) first (-k2n,2) and second column is both start and end key. Then this is sorted by first column alphabetically.
$ sort -k2n,2 -k1,1 joined_teosinte.txt > sorted_teosinte.txt				#### sorted joined_teosinte.txt

$ sort -k2nr,2 -k1,1 joined_maize.txt > reverse_sorted_maize.txt			#### reverse sorted joined_maize.txt
$ sort -k2nr,2 -k1,1 joined_teosinte.txt > reverse_sorted_teosinte.txt			#### reverse sorted joined_teosinte.txt


$ sed "s^?/?^??^g" sorted_maize.txt > sed_maize.txt					#### replace all ?/? with ?? in sorted_maize.txt
$ sed "s^?/?^??^g" sorted_teosinte.txt > sed_teosinte.txt				#### replace all ?/? with ?? in sorted_teosinte.txt

$ sed "s^?/?^--^g" reverse_sorted_maize.txt > sed_reverse_maize.txt			#### replace all ?/? with -- in reverse_sorted_maize.txt
$ sed "s^?/?^--^g" reverse_sorted_teosinte.txt > sed_reverse_teosinte.txt		#### replace all ?/? with -- in reverse_sorted_teosinte.txt


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


