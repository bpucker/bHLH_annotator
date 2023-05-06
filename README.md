# Automatic annotation of the bHLH gene family in plants

The bHLH_annotator allows the automatic identification and functional annotation of the bHLH transcription factor family in novel plant sequence data sets.  Coding sequences or peptide sequences derived from a de novo genome and transcriptome assembly can be analyzed with this pipeline. 

A phylogenetic approach is performed for the annotation of the candidates, based on a [bait collection](https://github.com/bpucker/bHLH_annotator/blob/main/supplements/) of bHLHs and outgroup sequences (non-bHLHs with a high sequence similarity to bHLHs):

<img alt="bHLH_annotator steps" src="https://github.com/bpucker/bHLH_annotator/blob/main/bHLH_annotator.png?raw=true" width="50%" height="70%">

For the identification of initial bHLH candidates (**step 1**), two search options are available:

 - **BLAST option** (default): Candidates are identified based on sequence similarity to the bait collection. This option is recommended if also bHLHs with a lost domain should be identified.
 - **HMMER option**: Candidates which harbour the HMM motif of the bait collection are identified. This includes candidates with a high specification, that are not represented by the bait collection.
 
The initial candidates are sorted out based on their phylogenetic relationship to the bHLH and outgroup baits (**step 2**). The functional annotation of the candidates is assigned by identifying ortholog reference sequences (**step 4**).  As default references, annotated *A. thaliana* bHLHs are used. Further, bHLH-specific characteristics are analyzed: Presence of the bHLH domain (**step 5**), DNA-binding properties (**step 5**), and the identification of subfamily specific motifs (**step 6**). A phylogenetic tree is constructed with *A. thaliana* bHLHs to allow a detailed investigation on the foundation of a well-studied species (**step 7**). 

For large datasets like *de novo* transcriptome assemblies, the collapse option is recommended (**step 8 and 9**) which collapses paralogous groups by defining a representative candidate. The parallel option is also recommenced to reduce the pipeline runtime and consumption of memory resources.

The data files used in each step can be customised by the user to allow an investigation suiting the own research purpose (described below). A more detailed description of the pipeline and the bait collection can be found [here](https://doi.org/10.1101/2023.05.02.539087). 

## Setup

### Installation in a conda environment 
The easiest way for installation is by creating an conda environment with the dedicated ```environment.yml``` file. This automatically installs all dependencies.
```batch
git clone https://github.com/bpucker/bHLH_annotator
cd bHLH_annotator
conda env create -f environment.yml
conda activate bHLH_annotator
```
### Manual installation of the dependencies
The following dependencies are necessary for the execution of the pipeline:
 - [Python3](https://www.python.org/): ```sudo apt install python3.11``` (other versions are also compatible)
	  - [dendropy]( https://dendropy.readthedocs.io/en/main/): ```sudo apt
   install python3-pip && python3 -m pip install -U dendropy```
	   - [pandas](https://pandas.pydata.org/docs/index.html): ```pip install pandas```
	   - [numpy](https://numpy.org/doc/stable/#): ```pip install numpy```
	   - [matplotlib](https://matplotlib.org/stable/index.html): ```pip  install  matplotlib```
 - [BLAST](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/): ```
sudo apt install ncbi-blast+```
- [HMMER](http://hmmer.org/documentation.html): ```conda install -c
  bioconda hmmer``` 
 - [MAFFT](https://mafft.cbrc.jp/alignment/software/linuxportable.html) : ```sudo apt install mafft```
 - [Muscle5](https://drive5.com/muscle5/) (precompiled binaries recommended)
- [FastTree2](http://www.microbesonline.org/fasttree/#Install): ```sudo
   apt-get install -y fasttree```
 - [RAxML-NG](https://github.com/amkozlov/raxml-ng) (precompiled binaries recommended)

The bHLH_annotator can be cloned from github:
```batch
git clone https://github.com/bpucker/bHLH_annotator
cd bHLH_annotator
```

## Usage 
The pipeline is executed through the following command:
```
cd <PATH>/bHLH_annotator
python3 bHLH_annotator.py --subject <PATH> --out <OUTPUT> --info <DEFINITION_FILE>
```
The ```--subject``` file defines the path to the input FASTA file containing coding or peptide sequences. The output directory is defined with the ```--out``` command. In the output directory, a RESULT folder is created containing the output files created in the pipeline steps. The ```--info```file represents the bHLH_annotator.csv. This file is necessary as it defines the input data files utilized in the pipeline.  

#### Optional arguments regarding subject file
|Command|Description|Default
|--|--|--
|```--name STR```    |Prefix of output file names|--
|```--cdsinput ```    |Changes expected input to CDS|--
|```--keepnames ```    |Prevents splitting of sequence names at first space|--
|```--collapse```    |Reduces paralogs to one representative|--
|```--parallel```    |Parallel option for classification|--
|```--paralogdist <FLOAT>```    |X*average nearest neighbour distance is used as cutoff to identify paralogs| ```10.0```
|```--numprocesscandidates <INT>```    |Number of candidates processed at the same time in the parallel option|```200``` 


#### Optional arguments for tool adjustments
|Command|Description|Default
|--|--|--
|```--search <STR>```    |Search option for the inistial search (```blast|hmmer```)|```blast```
|```--mode_aln <STR>```    |Alignment tool (```muscle|mafft```)|```muscle```
|```--mode_tree <STR>```    |Tool for tree construction (```fasttree|raxml```)|```fasttree```
|```--blastp <STR>```    |Path to blastp|```blastp```
|```--makeblastdb <STR>```    |Path to makeblastdb |```makeblastdb```
|```--hmmsearch <STR>```    |Path to hmmsearch|```hmmsearch```
|```--mafft <STR>```   |Path to MAFFT|```mafft```
|```--muscle <STR>```    |Path to muscle|```muscle```
|```--fasttree <STR>```    |Path to FastTree|```fasttree```
|```--raxml <STR>```    |Path to RAxML|```raxml```

  #### Optional arguments to adjust the candidate search and classification
|Command|Description|Default
|--|--|--
|```--bitcutp <INT>```    |BLASTp bitscore cutoff|```60```
|```--bitcutp	 <INT>```    | BLASTp similarity cutoff |```40.0``` 
|```--poscutp <INT>```    | Max number of BLASTp hits per bait|```100```
|```--lencutp <INT>```    |Min BLASTp alignment length|```80``` 
|```--filterdomain```    |Filter candidates not matching the HMM motif of the bait collection|-- 
|```--minscore <FLOAT>```    |Minimal score to be considered as ingroup|```0.5``` 
|```--numneighbours <INT>```    | Neighbours to consider for classification|```10``` 
|```--neighbourdist <FLOAT>```    | X*average nearest neighbour distance is used as minimal distance cutoff to be considered as a neighbour|```5``` 
|```--minneighbours <INT>```    |Minimal number of bait neighbours to be considered as ingroup|```0``` 

  #### Optional arguments regarding performance
|Command|Description|Default
|--|--|--
|```--cpu <INT>```    |Number of threads|```4```
|```--cpumax <INT>```    | Maximal number of threads for classification (step 2) |value of ```--cpu``` 
|```--cpub <INT>```    |Number of threads for BLASTp search (step 1)|value of ```--cpu```
|```--cpur <INT>```    |Number of threads for alignment/tree construction|value of ```--cpu```       


## Adjustment of input data files
The data input files required as resources for the steps of the pipeline are defined in the ```bHLH_annotator.csv``` file. The default files are stored in the ```data```folder. The following files are required as resources:

|File|Description|Argument
|--|--|--
Baits | Contains the bHLH and outgroup sequences of the bait collection (mandatory) |```--baits <PATH>```
BaitsInfo| Info file defining each bait as bHLH or outgroup sequence (mandatory) | ```--baitsinfo <PATH>``` 
OptimisedBaits | Optimised bait collection containing only phylogenetic distinct baits that are used for tree construction (recommended) | ```--optimisedbaits <PATH>```
HMM | HMM motif representing the bHLH domain (step 2 for HMMER search, step 5, and ```--filterdomain``` option)| ```--hmm <PATH>```
Reference | References with alternative name, functional annotation and subfamily (step 4 and 9) | ```--reference <PATH>```
Motifs | HMM motifs of subfamily specific motifs (step 5) | ```--motifs <PATH>```
Ath |*A. thaliana* sequences used for tree construction (step 7 and 9) | ```--ath <PATH>```

If the optional files are not defined in the ```bHLH_annotator.csv``` file or via argument, the dependent pipeline steps are skipped. If the files are defined both in the csv file and via argument, the argument is prioritized. The files refered to need to be placed in the 'data' folder or defined using complete paths.

 
## Requirements
Python, dendropy, pandas, numpy, matplotlib, BLAST+, HMMER, MAFFT or MUSCLE5, FastTree2 or RAxML 

## Reference
Thoben C. and Pucker B. (2023). Automatic annotation of the bHLH gene family in plants. bioRxiv 2023.05.02.539087; doi: [10.1101/2023.05.02.539087](https://doi.org/10.1101/2023.05.02.539087)




