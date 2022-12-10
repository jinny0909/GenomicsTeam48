# GenomicsTeam48
JHU Fall 22 Computational Genomics 
### Hee Yun Suh, David Pu, Kyoungjin Lim 
## 🧬 Topic: Protein Sequence Alignment Tool for E. Coli Identification

Alignment at the amino acid level has multiple advantages over alignment at the nucleotide level, primarily because of the larger amino acid alphabet and consequently lower signal-to-noise ratio within amino acid sequences compared to DNA [6]. Based on the methods employed by Bowtie[3], BWA[4], and PALADIN[7], we aimed to implement a simple protein sequence aligner capable of efficiently storing and querying a database of proteins using generated reads and different alignment algorithms.

##  ▶️ Running our code

All our BWA alignment, FM index, DP local alignment are called from ```testing_new.py``` file. 
After cloning our github repository, please run the following command:

``` 
python3 testing_new.py
```

If you would like to run your own alignment using a different database and different reads, you can use the align_tool function.
Please keep database in a fasta file, have the read in string format, decided on maximum_mismatch, use a=0 for exact matching, and (import align_tool.py) call the function:

``` 
align_tool(fasta_file, read, max_mismatch, a)
```

## 📂 Files and repository structure

|   FileName     |    Description   |  Link   |
| :-------------: |:-------------:| :-----:|
| align_tool.py | Alignment Tool for one Read |  [🔗](https://github.com/jinny0909/GenomicsTeam48/blob/main/align_tool.py) |
| blosum.py | Local Alignment using Dynamic Programming and BLOSUM 62 |  [🔗](https://github.com/jinny0909/GenomicsTeam48/blob/main/blosum.py) |
| bowtie.py | Bowtie Alignment implementation | [🔗](https://github.com/jinny0909/GenomicsTeam48/blob/main/bowtie.py) |  |
| bwa.py | BWA Alignment implementation | [🔗](https://github.com/jinny0909/GenomicsTeam48/blob/main/bwa.py)  |
| driver.py | Driver function for FM index | [🔗](https://github.com/jinny0909/GenomicsTeam48/blob/main/driver.py) |
| finalproject.py | Code used for midpoint presentation |  [🔗](https://github.com/jinny0909/GenomicsTeam48/blob/main/finalproject.py) |
| openRF.py | Helper functions to find Open Reading Frame | [🔗](https://github.com/jinny0909/GenomicsTeam48/blob/main/openRF.py) |
| proteinIndex.py | FM index and query helper functions |  [🔗](https://github.com/jinny0909/GenomicsTeam48/blob/main/proteinIndex.py)|
| proteinToDNA.py | Helper functions to convert protein to DNA  |  [🔗](https://github.com/jinny0909/GenomicsTeam48/blob/main/proteinToDNA.py)|
| testing_new.py | Main testing function for alignment tools  |  [🔗](https://github.com/jinny0909/GenomicsTeam48/blob/main/testing_new.py)|
| uniprot-compressed | Protein DB used for testing |  [💿](https://github.com/jinny0909/GenomicsTeam48/blob/main/uniprot-compressed_true_download_true_format_fasta_query_accession_3-2022.12.08-22.52.25.99.fasta.gz)|





## 🔖 External code referenced 

### David
- Notebooks from class; A lot of the information on FM Indexes was also pulled from the class notes, and I was not sure whether to cite those in the report.
- Initial FM index implementation was written purely based on lecture notes, but while debugging and searching for answers I came across an excellent implementation on GitHub: https://github.com/egonelbre/fm-index/tree/master/src. 
- I tried not to copy individual functions directly and mainly focused on restructuring my FMIndex class to isolate the actions of querying, resolving counts, and resolving offsets as the above implementation does. I did ultimately use their method of tracking the search bounds during a query as that was where my main difficulties were.

### Kyoungjin 

Although bowtie is not fully implemented, I have referred to the existing public repositories for bowtie implentation. 
- official bowtie repository:  https://github.com/BenLangmead/bowtie
- https://github.com/DonQwerty/Genome-Aligners/blob/master/src/bowtie.cpp A lot of the class and function structures are referred from this respository. 

### Hee Yun
- Because the BWA paper by Li and Durbin had very specific pseudo-code, that was followed.  This also means my code is similar to others (like the one I found below).
- https://github.com/Jwomers/burrows_wheeler_alignment/blob/012f2d97ab609c5d4404faab0dfac07d85896b7a/BWA.py
