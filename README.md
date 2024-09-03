# FM-index-prj

This is a script that do the FM index referring to a file FASTA (DNA or RNA).

To have a complete view I also define the command line arguments here:

- sys.arg[0] --> path/of/FMindexCode
- sys.arg[1] --> type of read function (the difference depends on if you have or not the file)
  - 1 if it's your file
  - 2 if you have the NCBI code
- sys.arg[2] --> path/of/FASTAsequence or the NCBI code (depends on the previously argument)
- sys.arg[3] --> method full or personalized 
  - full, if you want the full sequence
  - personalized, if you want a specific region
- sys.arg[4] --> the start of the region (only if you choose personalized)
- sys.arg[5] --> the stop of the region (only if you choose personalized)
- sys.arg[6:] --> the patterns

Arguments from 0 to 3 are mandatory (if you choose the personalized method, 4 and 5 are also mandatory)

The program must be run from the terminal (here some example of how run it):

Windows:
```bash
python path/of/FMindexCode 1 path/of/FASTAsequence full atg ggg cgtg

python path/of/FMindexCode 2 NG_005346.1 personalized 0 250 ggggg
```

macOS:

```bash
python3 path/of/FMindexCode 1 path/of/FASTAsequence full atg ggg cgtg

python3 path/of/FMindexCode 2 NG_005346.1 personalized 0 250 ggggg
```

Below there are a short descriptions of all the function created.

This is the class that contain all the functions
```py
class FM_index():
  ...
```

The function 
```py
def dir_creation():
  ...
```
it provides to create the source directory where inside it are saved all the sequence file obtained as output by all the other function

The function 
```py
def readFile():
  ...
```
it provides to the reading part of the input file. The function splitted in two part: the first one corresponds to the fact that the executor has the FASTA file on his computer; the second one corresponds to the fact that the executor doesn't have the FASTA file, but he has the NCBI code of the sequence (so thanks to the API he can use the code to take the FASTA file)

The function 
```py
def single_file_directory():
  ...
```
it provides to create a directory named with the header of the FASTA file and save it inside the the source directory previously created.

```py
def suffix_array():
  ...
```
The suffix is a sorted array (SA) of all suffixes of a string. It is a simple space efficient alternative to suffix trees.

The function
```py
def BWTransform():
  ...
```
it provides to do the borrows wheeler transform of the FASTA file. Rearranges a character string into runs of similar characters. This is useful for compression! More importantly, the transformation is reversible, without needing to store any additional data except the position of the first original character. The BWT is thus a "free" method of improving the efficiency of text compression algorithms, costing only some extra computation.
The transform is done by sorting all the circular shifts of a text in lexicographic order and by extracting the last column and the index of the original string in the set of sorted permutations of the start text.

The function
```py
def C_c():
  ...
```
it provides to create a dataframe, is a table that, for each character c in the alphabet, contains the number of occurrences of lexically smaller characters in the text.

The function [Occ(c, k)]
```py
def Occ_matrix():
  ...
```
it is the number of occurrences of character c in the prefix L[1..k]. Ferragina and Manzini (the creator of this method) showed that it is possible to compute Occ(c, k) in constant time.

The function 
```py
def pattern_input():
  ...
```
it provides to take as input patterns that the executor wants search inside the FASTA sequence.

The function
```py
def backward_search_for_pattern():
  ...
```
using the structures previously created it provides to calculate the number of occurrences of the pattern inside the sequence, allows to find the first index (position) of the pattern inside the sequence and also can allows to find all the indexes (positions) of the pattern.

The function
```py
def final_results():
  ...
```
it provides to take the output of the function backward search and stores them in a dataframe.

The function
```py
def DNA_plot():
  ...
```
provides for tracing the patterns found along the entire length of the sequence.

The function
```py
def save():
  ...
```
it provide to save all the results in the correct folder.