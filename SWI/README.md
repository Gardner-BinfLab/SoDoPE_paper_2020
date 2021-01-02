### Solubility prediction using SWI
Either the C version or Python version can be used. C version is a bit faster but rather limited. Python version is flexible as compared to C but maybe slightly slower. For comparision: C version took 5 secs whereas Python version took ~30 secs on i7 to calculate solubility of all reviewed sequences from UniProt (N=560,434).

#### C
Compile and run the program using the following commands. The input file must be a fasta file with an amino acid sequence(s).
```console
gcc swi.c -lm -o swi -std=gnu99
chmod +x swi
./swi <input.faa>
```
Please note that for the C version, the input fasta must be a single line fasta. For eg:
```
>sp|Q99TG3|ACCA_STAAM Acetyl-coenzyme A carboxylase carboxyl transferase subunit alpha OS=Staphylococcus aureus (strain Mu50 / ATCC 700699) OX=158878 GN=accA PE=3 SV=1
MLDFEKPLFEIRNKIESLKESQDKNDVDLQEEIDMLEASLERETKKIYTNLKPWDRVQIARLQERPTTLDYIPYIFDSFMELHGDRNFRDDPAMIGGIGFLNGRAVTVIGQQRGKDTKDNIYRNFGMAHPEGYRKALRLMKQAEKFNRPIFTFIDTKGAYPGKAAEERGQSESIATNLIEMASLKVPVIAIVIGEGGSGGALGIGIANKVLMLENSTYSVISPEGAAALLWKDSNLAKIAAETMKITAHDIKQLGIIDDVISEPLGGAHKDVEQQALAIKSAFVAQLDSLESLSRDEIANDRFEKFRNIGSYIE
>sp|P49145|5HT1D_RABIT 5-hydroxytryptamine receptor 1D OS=Oryctolagus cuniculus OX=9986 GN=HTR1D PE=3 SV=1
MSPSNQSAEGLPQEAANRSLNATGTPEAWDPGTLQALKISLAVVLSIITVATVLSNTFVLTTILLTRKLHTPANYLIGSLATTDLLVSILVMPISIAYTITHTWNFGQVLCDIWVSSDIT
```
The results can be saved by using `>`. For eg: `./swi sequences.fa > results.csv`

#### Python
Requirements

- Python 3.6+
- Pandas (Installation: `python3 -m pip install --user pandas`)
- Numpy (Installation: `python3 -m pip install --user numpy`)

Usage:

```python3 swi.py -f sequences.fa```

The results will be automatically saved in the current directory.

#### Cite:
- Bikash K Bhandari, Paul P Gardner, Chun Shen Lim. (2020). Solubility-Weighted Index: fast and accurate prediction of protein solubility. Bioinformatics. DOI:[10.1093/bioinformatics/btaa578](https://dx.doi.org/10.1093/bioinformatics/btaa578)
