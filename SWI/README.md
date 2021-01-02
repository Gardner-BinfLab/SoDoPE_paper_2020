Compile and run the program using the following commands. The input file must be a fasta file with an amino acid sequence(s).
```console
gcc swi.c -lm -o swi -std=gnu99
chmod +x swi
./swi <input.faa>
```
Please note that the input fasta should be a single line fasta. For eg:
```
>sp|Q99TG3|ACCA_STAAM Acetyl-coenzyme A carboxylase carboxyl transferase subunit alpha OS=Staphylococcus aureus (strain Mu50 / ATCC 700699) OX=158878 GN=accA PE=3 SV=1
MLDFEKPLFEIRNKIESLKESQDKNDVDLQEEIDMLEASLERETKKIYTNLKPWDRVQIARLQERPTTLDYIPYIFDSFMELHGDRNFRDDPAMIGGIGFLNGRAVTVIGQQRGKDTKDNIYRNFGMAHPEGYRKALRLMKQAEKFNRPIFTFIDTKGAYPGKAAEERGQSESIATNLIEMASLKVPVIAIVIGEGGSGGALGIGIANKVLMLENSTYSVISPEGAAALLWKDSNLAKIAAETMKITAHDIKQLGIIDDVISEPLGGAHKDVEQQALAIKSAFVAQLDSLESLSRDEIANDRFEKFRNIGSYIE
>sp|P49145|5HT1D_RABIT 5-hydroxytryptamine receptor 1D OS=Oryctolagus cuniculus OX=9986 GN=HTR1D PE=3 SV=1
MSPSNQSAEGLPQEAANRSLNATGTPEAWDPGTLQALKISLAVVLSIITVATVLSNTFVLTTILLTRKLHTPANYLIGSLATTDLLVSILVMPISIAYTITHTWNFGQVLCDIWVSSDIT
```
The results can be saved by using `>`. For eg: `./swi sequeneces.fa > results.csv`
