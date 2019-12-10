## Clustering how to
- Run the notebook Data_cleanup first.
- Use pET\_full_without\_his\_tag.fa from results folder with following command from usearch (needs to be downloaded)
```./usearch11.0.667_i86linux32 -cluster_fast pET_full_without_his_tag.fa  -id 0.1 -msaout msaout.fa -threads 4```
- copy all msaout.fa* files in clusters folder. 
- run Cross_\validation\_dataset

