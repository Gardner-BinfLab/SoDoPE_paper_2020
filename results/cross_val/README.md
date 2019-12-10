## Clustering how to
- Run the notebook Data_cleanup first.
- Use pET\_full_without\_his\_tag.fa from the results folder using the following command from USEARCH (needs to be downloaded separately from https://www.drive5.com/usearch/download.html).
```./usearch11.0.667_i86linux32 -cluster_fast pET_full_without_his_tag.fa  -id 0.1 -msaout msaout.fa -threads 4```
- Copy all msaout.fa* files in clusters folder. 
- Run Cross_\validation\_dataset

