Make directories with names 10, 20, 30, 40, 50, 60, 70, 80, 90
Copy .fa and usearch binary to each dir named 10 to 90
Run this inside each folder
./usearch11.0.667_i86linux32 -cluster_fast pET_full_without_his_tag.fa  -id XX -msaout msaout.fa -threads 4 

Replace -id XX by -id 0.1 for dir 10, 0.2 for dir 20, 0.3 for dir 30 etc..
Remove the fasta and usearch from the folder after clustering
