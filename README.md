# bHLH_annotator
automatic annotation of bHLH gene family in plants



table about all genome-wide bHLH studies (table)


bHLH land mark sequences (done)

bHLH sequences of all species (subfolder/archive)



------
 TOOL
------

cd <PATH>/bHLH_annotator
python3 family_annotator.py --subject <PATH> --out <OUTPUT>  

Mandatory:
--out           STR     Output folder
--subject       STR     Subject file   | --subjectdir    STR     Subject dir
					
Optional:
--family		LIST	Families the search is executed on (from family_info.csv)	[bHLH]
--search        STR     Tool for initial search (blast|hmmer)						[blast]	-> for option hmmer the HMM motif has to be defined in the info file
--mode_aln		STR		Tool for alignment											[muscle]
--mode_tree          STR     Tool for tree construction (fasttree|raxml)			[fasttree]

--name          STR     Prefix of output file names
--cdsinput      NONE    Changes expected input to CDS
--keepnames     NONE    Prevents splitting of sequence names at first space
--collapse      NONE    Reduces paralogs to one representative

--cpu           INT     Number of threads 												[4]
--cpub          INT     Number of threads for BLASTp									[cpu]
--cpur          INT     Number of threads for RAxML										[cpu]
					
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
--lencutp		INT     Min BLASTp alignment length 									[60]

	Candidate filtering settings:
--minscore		FLOAT	Minimal score to be considered as ingroup						[0.8]
--numneighbours INT     Neighbours to consider for classification 						[30]
--neighbourdist FLOAT   Cutoff in neighbour identification 							[3]
--minneighbours INT     Minimal number of bait neighbours for classification 			[0]
--paralogdist   FLOAT   Distance cutoff in paralog maksing 							[10]


-------------------
 FAMILY DEFINITION
-------------------
To use the annotator for a family, specific files have to be defined in family_info.csv. 
The files refered to need to be placed in the 'data' folder.  

Mandatory:
Family			Family name	
Baits			File name of baits file
BaitsInfo		File name of info file

Optional: 
ThinnedBaits	File name of reduced bait set file for tree construction, otherwise baits file is used
HMM				File name of HMM motif for hmmer search
Reference		File name of reference file
Ath				File name of Ath family members
RegEx			File name of RegEx file for domain check
Motifs			File name of motifs file for motif search