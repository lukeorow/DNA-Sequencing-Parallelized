# DNA GC Content Calculation

## About
This python-based tool was created to analyze DNA sequences and calculate their GC content. The GC content is a percentage of the DNA sequence that contains Guanine and Cytosine.
GC content is a metric used in bioinformatics that give insigts of DNA stability and can be used to find evidence of evolutionary relationships between organisms.

We implemented this tool both sequentially and with parallelism to show the runtime differences betwen the two. The sequential calculation uses a single thread to go through
the entire sequence. The multiprocessing parallel calculation uses python's ProcessPoolExecutor to utilize multiple CPU cores for parallelization. When tested on very large datasets,
this parallel approach was significantly quicker than the sequential calculation.

## How to Run
1. Have Python Downloaded
2. Have an IDE Downloaded (VSCode)
3. Clone this repo via Git or download as ZIP
4. Open in IDE
5. Prepare your DNA sequence file and name it dna_sequence.txt and put it in the root directory
   (file should only contain these four characters: 'A', 'T', 'G', and 'C')
7. Run the dna.py code
