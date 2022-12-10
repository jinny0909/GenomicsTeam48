# GenomicsTeam48
JHU Fall 22 Computational Genomics 
### Hee Yun Suh, David Pu, Kyoungjin Lim 
## 🧬 Topic: Protein Sequence Alignment Tool for E. Coli Identification

Alignment at the amino acid level has multiple advantages over alignment at the nucleotide level, primarily because of the larger amino acid alphabet and consequently lower signal-to-noise ratio within amino acid sequences compared to DNA [6]. Based on the methods employed by Bowtie[3], BWA[4], and PALADIN[7], we aimed to implement a simple protein sequence aligner capable of efficiently storing and querying a database of proteins using generated reads and different alignment algorithms.

## 📂 Files and repository structure




## 🔖 External code referenced 

### David
-Notebooks from class
-Initial FM index implementation was written purely based on lecture notes, but while debugging and searching for answers I came across an excellent implementation on GitHub: https://github.com/egonelbre/fm-index/tree/master/src. 
-I tried not to copy individual functions directly and mainly focused on restructuring my FMIndex class to isolate the actions of querying, resolving counts, and resolving offsets as the above implementation does.
