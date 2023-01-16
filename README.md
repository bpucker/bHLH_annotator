# bHLH_annotator
automatic functional annotation of the bHLH transcription factor family in plants

## Usage 

```
cd <PATH>/bHLH_annotator
python3 bHLH_annotator.py --subject <PATH> --out <OUTPUT> --info <DEFINITION_FILE>
```

```
Mandatory:
--out           STR     Output folder
--subject       STR     Subject file   | --subjectdir    STR     Subject dir
--info          STR     CSV definition file | --baits   STR    Baits file  --baitsinfo     STR   Baits info file  
					
Optional:
--family		LIST	Families the search is executed on (from family_info.csv)	[bHLH]
--search        STR     Tool for initial search (blast|hmmer)						[blast]	-> for option hmmer the HMM motif has to be defined in the info file
--mode_aln		STR		Tool for alignment											[muscle]
--mode_tree          STR     Tool for tree construction (fasttree|raxml)			[fasttree]

--name          STR     Prefix of output file names
--cdsinput      NONE    Changes expected input to CDS
--keepnames     NONE    Prevents splitting of sequence names at first space
--collapse      NONE    Reduces paralogs to one representative
--parameter_graphs      NONE    Create graphs showing candidate distribution for variation of blast and classification parameters

--cpu           INT     Number of threads 												[4]
--cpumax         INT     Number of threads for in ingroup/outgroup classification      [cpu]
--cpub          INT     Number of threads for BLASTp									[cpu]
--cpur          INT     Number of threads for alignment/tree construction				[cpu]
					
--blastp        STR     Path to blastp 													[blastp]
--makeblastdb   STR     Path to makeblastdb 											[makeblastdb]
--hmmsearch     STR     Path to hmmsearch 												[hmmsearch]
--mafft         STR     Path to MAFFT 													[mafft]
--muscle		STR		Path to muscle													[muscle]
--raxml         STR     Path to RAxML 													[raxml]				
--fasttree      STR     Path to FastTree 												[fasttree]

	Candidate search settings: 
--bitcutp		INT		BLASTp bitscore cutoff											[60] 					
--simcutp       FLOAT   BLASTp similarity cutoff 										[40.0]
--poscutp       INT     Max number of BLASTp hits per bait 								[100]
--lencutp		INT     Min BLASTp alignment length 									[80]

	Candidate filtering settings:
--filterdomain  NONE    Filters candidates by hmm motif -> only if HMM motif is defined in the info file
--minscore		FLOAT	Minimal score to be considered as ingroup						[0.5]
--numneighbours INT     Neighbours to consider for classification 						[10]
--neighbourdist FLOAT   Cutoff in neighbour identification 							[5.0]
--minneighbours INT     Minimal number of bait neighbours for classification 			[0]
--parallel    NONE      Parallelization at ingroup/outgroup classification
--numprocesscandidates  INT Maximal number of candidates per parallel ingroup/outgroup classification  [200]
--paralogdist   FLOAT   Distance cutoff in paralog maksing 							[10]
```

## Family definition
To use the annotator for a family, specific files have to be defined in bHLH_annotator.csv which is given via the "--info" argument. 
The files refered to need to be placed in the 'data' folder or defined using complete paths.  

```
Mandatory:
Family			Family name	
Baits			File name of baits file
BaitsInfo		File name of info file

Optional: 
ThinnedBaits	File name of reduced bait set file for tree construction, otherwise baits file is used
HMM				File name of HMM motif for hmmer search
Reference		File name of reference file
Ath				File name of Ath family members
Motifs			File name of motifs file for motif search
```

## Requirements
dendropy, BLAST, HMMER, MAFFT, MUSCLE, RAxML or FastTree2

## Reference
This repository.



