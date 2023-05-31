### most functions were copied from MYB_annotator: http://dx.doi.org/10.1186/s12864-022-08452-5
### some functions were copied from KIPEs: https://doi.org/10.3390/plants9091103 and MaMYB: https://doi.org/10.1371/journal.pone.0239275 ###

### WARNING: do not use underscores in the bait IDs ###

### definition of bHLH family in data/bHLH_annotator.csv:
# Baits	BaitsInfo	ThinnedBaits	HMM	Reference	Ath	Motifs
# bHLH_baits.fasta	bHLH_baits.txt	bHLH_baits_thinned.fasta	bHLH_baits_thinned.hmm	AthRefbHLHs.txt	bHLH_at.fasta	bHLH_motifs.hmm


import os, glob, sys, subprocess, dendropy, datetime
from operator import itemgetter
import math, time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
try:
	import hashlib
except ImportError:
	pass

__version__ = "v1.03"
__usage__ = """
                    python3 bHLH_annotator.py
                    --out <OUTPUT_DIR>
                    --subject <SUBJECT_FILE (peptide,transcript,genomic sequences)> | --subjectdir <SUBJECT_FOLDER_WITH_SEQ_FILES>
                    --info <PATH_TO_CSV_DEFINITION_FILE> | --baits <PATH_TO_BAIT_FILE> --baitsinfo <PATH_TO_BAITS_INFO_FILE>
                    
                    
                    optional:
                    --search <INITIAL_SEARCH_TOOL>(blast|hmmer)[blast]
                    --mode_aln <ALIGNMENT_TOOL>(muscle|maddt)[muscle]
                    --mode_tree <TREE_TOOL>(fasttree|raxml)[fasttree]
                    --name <STRING_USED_AS_PREFIX_IN_FILENAMES>                   
                    --cdsinput <CHANGES_EXPECTED_INPUT_TO_CDS>
                    --keepnames <PREVENTS_CUTTING_OF_NAMES_AT_FIRST_SPACE>                    
                    --collapse <REDUCES IN-PARALOGS_TO_ONE_REPRESENTATIVE>
                    --parameter_graphs<CREATES_GRAPHS_WITH_PARAMETER> 
                    
                    --cpu <NUMBER_OF_THREADS>[4]
                    --cpumax <MAX_CPUs_FOR_IN_OUT_CLASSIFICATION>[cpu]
                    --cpub <CPUs_TO_USE_FOR_BLASTp>[cpu]
                    --cpur <CPUs_TO_USE_FOR_RAxML>[cpu]
                                                                                
                    --blastp <PATH_TO_AND_INCLUDING_BINARY>[blastp]
                    --hmmsearch <PATH_TO_HMMSEARCH>[hmmsearch]
                    --makeblastdb <PATH_TO_AND_INCLUDING_BINARY>[makeblastdb]
                    --mafft <PATH_TO_MAFFT>[raxml]              
                    --muscle <PATH_TO_MUSCLE>[raxml]                
                    --fasttree <PATH_TO_FASTTREE>[fasttree]
                    --raxml <PATH_TO_RAXML>[raxml-ng]  
		                
                    --simcutp <BLASTP_SIMILARITY_CUTOFF>[40.0]
                    --poscutp <BLASTP_POSSIBLE_HIT_NUMBER_PER_BAIT_CUTOFF>[100]
                    --lencutp <BLASTP_MIN_LENGTH_CUTOFF>[80]
                    --bitcutp <BLASTP_BITSCORE_CUTOFF>[60]
                    
                    --filterdomain <DOMAIN_FILTER_FOR_CLASSIFICATION> [False]
                    --minscore <MINIMAL_SCORE> [0.5]
                    --numneighbours <NUMBER_OF_NEIGHBOURS_FOR_CLASSIFICATION> [10]
                    --neighbourdist <NEIGHBOUR_DISTANCE> [5]
                    --minneighbours <MINIMAL_NUMBER_OF_NEIGHBOURS> [0]
                    --parallel <PARALLEL_CLASSIFICATION_MODE> [False]
                    --numprocesscandidates <NUMBER_OF_CANDIDATES_PER_PARALLEL_CLASSIFICATION> [200]
                    --paralogdist <DISTANCE_OF_PARALOGS_IN_MASKING_STEP> [10]
                           
                          
                    bug reports and feature requests: b.pucker@tu-braunschweig.de
                    """
#%%

def translate( seqs ):
    """! @brief translates the given nucleotide sequence into peptide and splits at each star (stop codon) """
    
    genetic_code = {    'CTT': 'L', 'ATG': 'M', 'AAG': 'K', 'AAA': 'K', 'ATC': 'I',
                                'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T',
                                'AGC': 'S', 'ACA': 'T', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N',
                                'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H',
                                'ACG': 'T', 'CCG': 'P', 'AGT': 'S', 'CAG': 'Q', 'CAA': 'Q',
                                'CCC': 'P', 'TAG': '*', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C',
                                'CGA': 'R', 'CCA': 'P', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R',
                                'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G',
                                'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'GAG': 'E', 'TCG': 'S',
                                'TTA': 'L', 'GAC': 'D', 'TCC': 'S', 'GAA': 'E', 'TCA': 'S',
                                'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A',
                                'GTG': 'V', 'TTC': 'F', 'GTT': 'V', 'GCT': 'A', 'ACC': 'T',
                                'TTG': 'L', 'CGT': 'R', 'TAA': '*', 'CGC': 'R'
                            }
    
    final_peptide_seqs = {}
    for key in seqs.keys():
        seq = seqs[ key ].upper()
        peptide = []
        for i in range( int( len( seq ) / 3.0 ) ):
            codon = seq[i*3:i*3+3]
            try:
                peptide.append( genetic_code[ codon ] )
            except:
                peptide.append( "*" )
        final_peptide_seqs.update( { key: "".join( peptide ) } )
    return final_peptide_seqs


def clean_input_FASTA_file( raw_subject_file, subject_file, mapping_table, cds_input, trim_names ):
    """! @brief clean input FASTA file """
    
    forbidden_characters = [ ";", ":", "(", ")", "_", "=" ]
    
    with open( mapping_table, "w" ) as out:
        out.write( "InitialID\tCleanID\n" )
        sequences = {}
        with open( raw_subject_file ) as f:
            header = f.readline()[1:].strip()
            if trim_names:
                if " " in header:
                    header = header.split(' ')[0]
                    if "\t" in header:
                        header = header.split('\t')[0]
            out.write( header + "\t" )
            if " " in header:
                header = header.split(' ')[0]
            if "\t" in header:
                header = header.split('\t')[0]
            for each in forbidden_characters:
                header = header.replace( each, "-" )
            header = header.encode("ascii", "ignore").decode()    #removal of non-ASCII characters
            out.write( header + "\n" )
            seq = []
            line = f.readline()
            while line:
                if line[0] == '>':
                        sequences.update( { header: "".join( seq ) } )
                        header = line.strip()[1:]
                        if trim_names:
                            if " " in header:
                                header = header.split(' ')[0]
                                if "\t" in header:
                                    header = header.split('\t')[0]
                        out.write( header + "\t" )
                        if " " in header:
                            header = header.split(' ')[0]
                        if "\t" in header:
                            header = header.split('\t')[0]
                        for each in forbidden_characters:
                            header = header.replace( each, "-" )
                        header = header.encode("ascii", "ignore").decode()
                        out.write( header + "\n" )
                        seq = []
                else:
                    seq.append( line.strip() )
                line = f.readline()
            sequences.update( { header: "".join( seq ) } )
    
    if cds_input:
        cds_file = subject_file.replace( ".pep.fasta", ".cds.fasta" )
        pep_sequences = translate( sequences )
        with open( cds_file, "w" ) as out:    #construct file with clean CDS
            for key in sequences.keys():
                out.write( '>' + key + "\n" + sequences[ key ] + "\n" )
    
        with open( subject_file, "w" ) as out:    #construct file with clean PEPs
            for key in pep_sequences.keys():
                out.write( '>' + key + "\n" + pep_sequences[ key ] + "\n" )
            
    else:
        with open( subject_file, "w" ) as out:    #construct file with clean PEPs
            for key in sequences.keys():
                out.write( '>' + key + "\n" + sequences[ key ] + "\n" )


def load_sequences( fasta_file ):
    """! @brief load candidate gene IDs from file """
    
    sequences = {}
    with open( fasta_file ) as f:
        header = f.readline()[1:].strip()
        seq = []
        line = f.readline()
        while line:
            if line[0] == '>':
                    sequences.update( { header: "".join( seq ) } )
                    header = line.strip()[1:]
                    seq = []
            else:
                seq.append( line.strip() )
            line = f.readline()
        sequences.update( { header: "".join( seq ) } )    
    return sequences


def load_bait_fam_anno( info_file, seq_file ):
    """! @brief load IDs into two list (in and out) """
    
    ins, outs = [], []
    in_list, out_list = [], []
    sequences = load_sequences(seq_file)
    
    with open( info_file, "r" ) as f:
        line = f.readline()
        while line:
            parts = line.strip().split('\t')
            if parts[1] == "in":
                ins.append( parts[0] )
            else:
                outs.append( parts[0] )
            line = f.readline()
            
    for s in sequences:
        if s in ins:
            in_list.append(s)
        elif s in outs:       
            out_list.append(s)
             
    return in_list, out_list


def load_fam_classification_from_file( tmp_result_table ):
    """! @brief load family member classification from file """
    
    fam_classification = {}
    with open( tmp_result_table, "r" ) as f:
        f.readline()    #remove header
        line = f.readline()
        while line:
            parts = line.strip().split('\t')
            fam_classification.update( { parts[1]: float( parts[4] ) } )
            line = f.readline()
    return fam_classification


def check_fam_IDs_across_files( fam_bait_seq_file, fam_bait_seq_file_thinned,fam_info_file, fam_ath_file,ref_file, fam ):
    """! @brief check IDs across the different files """
    
    forbidden_characters = [ ";", ":", "(", ")", "_" ]
    fam_status = True
    # --- check FASTA file for forbidden characters --- #
    seqs_all = load_sequences( fam_bait_seq_file )
    header_string = "".join( seqs_all.keys() )
    for each in forbidden_characters:
        if each in header_string:
            sys.stderr.write( "Forbidden character detected in " + fam + " IDs (bait FASTA file): " + each + " (occurrences:" + str( header_string.count( each ) ) + ")\n" )
            sys.stderr.flush()
            fam_status = False
    
    # --- check structure of info file --- #
    info_IDs = {}
    with open( fam_info_file, "r" ) as f:
        line = f.readline()
        while line:
            parts = line.strip().split('\t')
            if len( parts ) < 2:
                sys.stderr.write( "Issue in " + fam + " info file (number of columns not 2): " + line + "\n" )
                sys.stderr.flush()
            else:
                info_IDs.update( { parts[0]: None } )
                if parts[1] not in [ "in", "out" ]:
                    sys.stderr.write( "Issue in " + fam + " info file (unexpected status; only 'in' and 'out' permitted): " + line.replace("\n","") + "\n" )
                    sys.stderr.flush()
            line = f.readline()

    # --- check thinned FASTA file for forbidden characters --- #
    seqs = load_sequences( fam_bait_seq_file_thinned )
    header_string = "".join( seqs.keys() )
    for each in forbidden_characters:
        if each in header_string:
            sys.stderr.write( "Forbidden character detected in " + fam + " IDs (thinned baits FASTA file): " + each + " (occurrences:" + str( header_string.count( each ) ) + ")\n" )
            sys.stderr.flush()
            fam_status = False

    # --- compare IDs between thinned bait and info file --- #
    missing_in_fasta = []
    for seq_ID in info_IDs.keys():
        try:
            seqs_all[ seq_ID ]
        except KeyError:
            missing_in_fasta.append( seq_ID )
    if len( missing_in_fasta ) > 0:
        sys.stderr.write( "Unmatched " + fam + " IDs (missing in bait FASTA file): " + ";".join( missing_in_fasta ) + "\n" )
        sys.stderr.flush()
        fam_status = False
    
    missing_in_info = []
    for seq_ID in seqs.keys():
        try:
            info_IDs[ seq_ID ]
        except KeyError:
            missing_in_info.append(seq_ID )
    if len( missing_in_info ) > 0:
        sys.stderr.write( "Unmatched " + fam + " IDs (missing in info file): " + ";".join( missing_in_info ) + "\n" )
        sys.stderr.flush()
        fam_status = False
    
    # --- check FASTA file for forbidden characters --- #
    seqs_ath = []
    if len(fam_ath_file) > 0:
        seqs_ath = load_sequences( fam_ath_file )
        header_string = "".join( seqs.keys() )
        for each in forbidden_characters:
            if each in header_string:
                sys.stderr.write( "Forbidden character detected in Ath" + fam + " IDs (Ath %s members FASTA file): " % (fam) + each + " (occurrences:" + str( header_string.count( each ) ) + ")\n" )
                sys.stderr.flush()
                fam_status = False
    
    # --- compare IDs between ref and ath/all baits fasta file --- #        
    if len( ref_file ) > 0:    #only try to check if file actually exists
        missing_ref_ids = []
        with open( ref_file, "r" ) as f:
            line = f.readline()
            while line:
                x = line.strip()
                if "\t" in x:
                    x = line.split('\t')[0]
                try:
                    seqs_all[ x ]
                except KeyError:
                    try:
                        seqs_ath[ x ]
                    except KeyError:
                        missing_ref_ids.append( x )
                line = f.readline()
        if len( missing_ref_ids ) > 0:
            sys.stderr.write( "Reference " + fam + " IDs missing in bait collection or Ath members FASTA file): " + ";".join( missing_ref_ids ) + "\n" )
            sys.stderr.flush()
            fam_status = False
        
    return fam_status


def load_subject_name_mapping_table( mapping_table_file ):
    """! @brief load subject name mapping table """
    
    mapping_table = {}
    with open( mapping_table_file, "r" ) as f:
        line = f.readline()
        while line:
            parts = line.strip().split('\t')
            mapping_table.update( { parts[1]: parts[0] } )    #clean ID back to original one
            line = f.readline()
    return mapping_table


def load_ref_members( ref_members_file ):
    """! @brief load IDs from given file """
    
    refmem = {}
    with open( ref_members_file, "r" ) as f:
        line = f.readline()
        while line:
            if "\t" in line:
                parts = line.strip().split("\t")
                if len( parts ) == 4:
                    refmem.update( { parts[0]: { 'id': parts[0], 'name': parts[2], 'function': parts[3], 'group': "" } } )
                elif len( parts ) == 5:
                    refmem.update( { parts[0]: { 'id': parts[0], 'name': parts[2], 'function': parts[3], 'group': parts[4] } } )                
                if len(refmem[parts[0]]['name']) < 1:
                    refmem[parts[0]]['name'] = parts[0]                            
            else:
                refmem.update( { line.strip(): { 'id': parts[0], 'name': line.strip(), 'function': "n/a" , 'group': "" } } )
            line = f.readline()
    return refmem


def load_motifs_from_file( motifs_file ):
    """! @brief load motifs from file """
    
    motifs = {}
    with open( motifs_file, "r" ) as f:
        line = f.readline()
        while line:
            if "\t" in line:
                parts = line.strip().split()
                motifs.update( { parts[0]: parts[1] } )
            else:
                motifs.update( { "motif-" + str( len( motifs.keys() ).zfill(3) +1 ): line.strip() } )
            line = f.readline()
    return motifs


def md5_calculator( input_file ):
    """! @brief calculate md5sum of given file """
    
    with open( input_file, "rb" ) as f:
        content = f.read()
    try:
        return hashlib.md5( content ).hexdigest()
    except NameError:
        return "n/a"


def generate_documentation_file(doc_file,fam_bait_seq_file_all ,fam_bait_seq_file, fam_info_file, fam_hmm_motif, fam_ath_file, motifs_file ,output_folder, raw_subject_file,
                                search, mode_aln, mode_tree, blastp, makeblastdb, hmmsearch, cpumax ,cpub, cpur, mafft, muscle, raxml, fasttree, ref_file,
                                bitscore_cutoff_p, similarity_cutoff_p, possibility_cutoff_p, length_cutoff_p, cds_input, 
                                min_score_cutoff, neighbour_cutoff, mean_factor_cutoff, min_neighbour_cutoff, dist_cutoff_factorB, fam,
                                filter_domain, parallel_mode, num_process_candidates, name, keepnames, collapse,parameter_graphs):
    """! @brief write documentation file with specified inputs and parameters """
    
    with open( doc_file, "w" ) as out:

        out.write( "bHLH_annotator.py version: " + __version__ + "\n" )
        fam_bait_seq_file_md5 = md5_calculator( fam_bait_seq_file_all )
        out.write( fam + " bait file: " + fam_bait_seq_file_all + "\t" + fam_bait_seq_file_md5 + "\n" )
        fam_info_file_md5 = md5_calculator( fam_info_file )
        out.write( fam + " bait info file: " + fam_info_file + "\t" + fam_info_file_md5 + "\n" )
        fam_reduced_bait_file_md5 = md5_calculator( fam_bait_seq_file )
        out.write(" thinned bait file: " + fam_bait_seq_file + "\t" + fam_reduced_bait_file_md5 + "\n" )
        raw_subject_file_md5 = md5_calculator( raw_subject_file )
        out.write( "Subject FASTA file: " + raw_subject_file + "\t" + raw_subject_file_md5 + "\n" )
        out.write( "Output folder: " + output_folder + "\n" )
        
        #--- optional --- #
        out.write( "Tool for initial candidate selection: " + search + "\n" )
        out.write( "Tool for alignment: " + mode_aln + "\n" )
        out.write( "Tool for tree construction: " + mode_tree + "\n" )
        out.write( "Maximal CPUs for ingroup/outgroup classification: " + str( cpumax ) + "\n" )
        out.write( "CPUs for BLASTp: " + str( cpub ) + "\n" )
        out.write( "CPUs for RAxML: " + str( cpur ) + "\n" )
        if len(fam_hmm_motif) > 0:
            fam_bait_motif_file_md5 = md5_calculator(fam_hmm_motif )
            out.write("HMM bait motif file: " + fam_hmm_motif + "\t" + fam_bait_motif_file_md5  + "\n" )
        else:
            out.write( "HMM bait motif file: " + "n/a\n" )
        if len( ref_file ) > 0:
            ref_file_md5 = md5_calculator( ref_file )
            out.write( "Reference " + fam + " file: " + ref_file + "\t" + ref_file_md5 + "\n" )
        else:
            out.write( "Reference " + fam + " file: n/a\n" )
        if len(fam_ath_file) > 0:
            fam_at_file_md5 = md5_calculator(fam_ath_file )
            out.write("Ath family members file: " + fam_ath_file + "\t" + fam_at_file_md5  + "\n" )
        else:
            out.write( "Ath family members file: " + "n/a\n" )
        if len(motifs_file) > 0:
            fam_motifs_file_md5 = md5_calculator(motifs_file)
            out.write("Motif search file: " + motifs_file + "\t" + fam_motifs_file_md5  + "\n" )
        else:
            out.write( "Motif search file: " + "n/a\n" )
        if cds_input:
            out.write( "Type of input: CDS\n" )
        else:
            out.write( "Type of input: PEP\n" )
        out.write( "Prefix of output file names: " + name + "\n" )
        out.write( "Preventing of splitting sequence names at first space (--keepnames): " + str(keepnames) + "\n" )
        out.write( "Reducing paralogs to one representative (--collapse): " + str(collapse) + "\n" )
        out.write( "Create graphs for parameter values (--parametergraphs): " + str(parameter_graphs) + "\n" )
        
        # --- paths to tools --- #
        out.write( "blastp path: " + blastp + "\n" )
        out.write( "makeblastdb path: " + makeblastdb + "\n" )
        out.write( "hmmsearch path: " + hmmsearch + "\n" )
        out.write( "mafft tool path: " + mafft + "\n" )
        out.write( "muscle tool path: " + muscle + "\n" )
        out.write( "raxml path: " + raxml + "\n" )
        out.write( "fasttree path: " + fasttree + "\n" )
        
        # ---- BLAST filter criteria --- #
        out.write( "Minimal bitscore cutoff: " + str( bitscore_cutoff_p ) + "\n" )
        out.write( "Minimal BLASTp hit similarity cutoff: " + str( similarity_cutoff_p ) + "\n" )
        out.write( "Maximal number of BLASTp hits per bait: " + str( possibility_cutoff_p ) + "\n" )
        out.write( "Minimal BLASTp hit alignment length: " + str( length_cutoff_p ) + "\n" )
        
        # --- tree analysis settings --- #
        out.write( "Minimal score for candidate to be considered as ingroup: " + str( min_score_cutoff ) + "\n" )
        out.write( "Number of neighbours to consider in classification: " + str( neighbour_cutoff ) + "\n" )
        out.write( "Factor for branch length cutoff in identification of neighbours: " + str( mean_factor_cutoff ) + "\n" )
        out.write( "Minimal number of neighbours required for classification as " + fam + ": " + str( min_neighbour_cutoff ) + "\n" )
        out.write( "Filter candidates by hmm motif: " + str(filter_domain) + "\n")
        out.write( "Distance cutoff for paralog clade masking: " + str( dist_cutoff_factorB ) + "\n" )
          
        out.write( "Parallelization at ingroup/outgroup classification: " + str(parallel_mode) + "\n")
        out.write( "Maximal number of candidates per parallel ingroup/outgroup classification: " + str(num_process_candidates) + "\n" )
           
        # --- add tool versions --- #
        try:
            muscle_version_raw = subprocess.Popen( args= muscle + " -help", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True )
            muscle_version =  muscle_version_raw.stdout.read()
            out.write ( "Muscle version: " + str(  muscle_version )[4:20] + "\n" )	#remove characters introduced through binary
        except:
            out.write ( "Muscle version detection failed.\n" )	#if no MAFFT installation was detected
        try:
            mafft_version_raw = subprocess.Popen( args=mafft + " --version", stderr=subprocess.PIPE, shell=True )
            mafft_version = mafft_version_raw.stderr.read()
            out.write ( "MAFFT version: " + str( mafft_version )[2:-3] + "\n" )	#remove characters introduced through binary
        except:
            out.write ( "MAFFT version detection failed.\n" )	#if no MAFFT installation was detected
        out.write ( "FastTree version: PLEASE_ADD_MANUALLY\n"  )	#version not available via command
        try:
            raxml_version_raw = subprocess.Popen( args=raxml + " --version", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True )
            raxml_version = str( raxml_version_raw.stdout.read() ).strip()
            out.write ( "RAxML version: " + ( raxml_version[4:65]) + "...\n" )    #remove characters introduced through binary
        except:
            out.write ( "RAxML version detection failed.\n" )    #if no RAxML installation was detected
        try:
            hmmsearch_version_raw = subprocess.Popen( args=hmmsearch + " -h", stdout=subprocess.PIPE, shell=True )
            hmmsearch_version = str( hmmsearch_version_raw.stdout.read() ).strip().split("#")[2]
            out.write ( "hmmsearch version: " + ( hmmsearch_version ) + "...\n" )    #remove characters introduced through binary
        except:
            out.write ( "hmmsearch version detection failed.\n" )    #if no hmmsearch installation was detected


def load_BLAST_results( blast_result_file, similarity_cutoff, possibility_cutoff, length_cutoff,bitscore ):
	"""! @brief load BLAST results """
	
	valid_blast_hits = {}
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if float( parts[2] ) > similarity_cutoff:	#similarity is sufficient
				if float( parts[3] ) > length_cutoff:	#substantial part of query is matched
					if float(parts[-1]) >= bitscore:
    					 try:
    						 valid_blast_hits[ parts[1] ].append( { 'gene': parts[0], 'score': float( parts[-1] ) } )
    					 except KeyError:
    						 valid_blast_hits.update( { parts[1]: [ { 'gene': parts[0], 'score': float( parts[-1] ) } ] } )
			line = f.readline()
	
	# --- reduce BLAST hit number to given number of candidate possibilities ---- #
	final_valid_blast_hits = {}
	for key in list(valid_blast_hits.keys()):
		hits = sorted( valid_blast_hits[ key ], key=itemgetter( 'score' ) )[::-1]
		genes = []
		for hit in hits:
			if hit['gene'] not in genes:
				if len( genes ) < possibility_cutoff:
					genes.append( hit['gene'] )
		final_valid_blast_hits.update( { key: genes } )
	
	return final_valid_blast_hits


def analyze_blast_result(result_folder, blast_result):
    """! @brief analyzes the number of identified candidates for a variation of blast parameters """
    if not os.path.isdir(result_folder):
        os.mkdir(result_folder)
    
    data = pd.read_csv( blast_result, delimiter = "\t", header=None)    
    data["id"] = data[1] #new candidate id
    data["bit_score"] = data[11]
    data["length"] = data[3]    
    data["similarity"] = data[2] 
    data["hit"] = data[0]
    data = data.reindex(columns=["id","bit_score","length","similarity","hit"])
      
    param_figures = result_folder + "blast_parameter_"
    for param in ["bit_score","length","similarity"]:        
        x, y = [], []     
        for val in range(0,100,10):
            x.append(val)
            y.append( len(data[data[param] > val].drop_duplicates(subset="id")))
        fig = plt.figure()
        plt.plot(x,y,label=param)
        plt.xlabel(param)
        plt.ylabel("number of candidates")
        plt.tight_layout()
        plt.savefig(param_figures + param +".jpg",dpi=300)
        plt.close(fig)

    combined_figures = result_folder + "1_blast_combined_simcutp_"
    for similarity_cutoff in range(10,100,10):
        x,y,z = [],[],[] 

        for bitscore in range(0,100,10):
            for length in range(0,100,10):
                x.append(bitscore)
                y.append(length)               
                z.append(len(data[(data["bit_score"] >= bitscore) & (data["length"] > length) & (data["similarity"]>similarity_cutoff)].drop_duplicates(subset="id")))   
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.plot_trisurf(np.array(x), np.array(y), np.array(z), alpha = 0.8, cmap = "viridis")
        ax.set_xlabel("bit_score")
        ax.set_ylabel("length")
        ax.set_zlabel("number of candidates")     
        plt.title("--simcutp: " + str(similarity_cutoff))  
        plt.tight_layout()
        plt.savefig(combined_figures + str(int(similarity_cutoff)) +".jpg",dpi=300)
        plt.close(fig)

def load_hmmsearch_results( seq_search_result_file ):
    """! @brief load all hmmsearch hits into a dictionary """
    
    hmm_search_results = {}
    with open( seq_search_result_file, "r" ) as f:
        line = f.readline()
        while line:
            if line[0] != "#":
                if "\t" in line:
                    hmm_search_results.update( { line.split('\t')[0]: None } )
                elif " " in line:
                    hmm_search_results.update( { line.split(' ')[0]: None } )
                else:
                    print( line )
            line = f.readline()
    return hmm_search_results


def load_alignment( aln_file, tmp_mapping ):
    """! @brief load alignment and replace query IDs by real sequence names """
    
    sequences = {}
    
    with open( aln_file ) as f:
        header = f.readline()[1:].strip()
        try:
            header = tmp_mapping[ header ]
        except KeyError:
            pass
        seq = []
        line = f.readline()
        while line:
            if line[0] == '>':
                    sequences.update( { header: "".join( seq ) } )
                    header = line.strip()[1:]
                    try:
                        header = tmp_mapping[ header ]
                    except KeyError:
                        pass
                    seq = []
            else:
                seq.append( line.strip() )
            line = f.readline()
        sequences.update( { header: "".join( seq ) } )    
    return sequences


def alignment_trimming( aln_file, cln_aln_file, occupancy ):
    """! @brief remove all alignment columns with insufficient occupancy """
    
    alignment = load_alignment( aln_file, {} )
    # --- if there is an alignment (expected case) 
    if len( list(alignment.keys()) ) > 0:
        # --- identify valid residues in aligned sequences (columns with sufficient occupancy) --- #
        valid_index = []
        for idx, aa in enumerate( list(alignment.values())[0] ):
            counter = 0
            for key in list(alignment.keys()):
                if alignment[ key ][ idx ] != "-":
                    counter += 1
            if counter / float( len( list(alignment.keys()) ) ) > occupancy:
                valid_index.append( idx )
        
        # --- generate new sequences --- #
        with open( cln_aln_file, "w" ) as out:
            for key in list(alignment.keys()):
                seq = alignment[ key ]
                new_seq = []
                for idx in valid_index:
                    new_seq.append( seq[ idx ] )
                out.write( ">" + key + '\n' + "".join( new_seq ) + '\n' )
    # --- just in case the alignment file is empyt (is this possible?) ---#
    else:
        with open( cln_aln_file, "w" ) as out:
            out.write( "" )


def tree_constructor( X_aln_input_file, X_aln_file, X_cln_aln_file, X_bait_seq_file, X_cand_file, X_output_folder, Xname, Xnumber, mode_aln, mode_tree, mafft, muscle, raxml, fasttree, cpur ):
    """! @brief handles the construction of alignments and phylogenetic tree
            @note second FASTA file can be an empty string to run this function just based on one FASTA file
    """
    
    if not os.path.isfile( X_aln_input_file ):
        if len( X_cand_file ) > 0:
            p = subprocess.Popen( args= "cat " + X_bait_seq_file + " " + X_cand_file + " > " + X_aln_input_file, shell=True )
            p.communicate()
        else:
            p = subprocess.Popen( args= "cp " + X_bait_seq_file + " " + X_aln_input_file, shell=True )
            p.communicate()
      
    if not os.path.isfile( X_aln_file ):       
        if mode_aln == "mafft":
            p = subprocess.Popen( args= mafft + " --quiet --thread " + str(cpur) + " " + X_aln_input_file + " > " + X_aln_file, shell=True )     
        else:
            p = subprocess.Popen(args = muscle + " -align " + X_aln_input_file + " -output " + X_aln_file  + " -threads " + str(cpur), shell=True)
        p.communicate()
               
    if not os.path.isfile( X_cln_aln_file ):
        alignment_trimming( X_aln_file, X_cln_aln_file, occupancy=0.1 )
    
    if mode_tree == "raxml":    #RAxML
        prefix = X_output_folder + Xname + Xnumber + "RAxML_tree"
        tree_file = prefix + ".raxml.bestTree"
        if not os.path.isfile( tree_file ):
            p = subprocess.Popen( args= " ".join( [ raxml, "--all --threads " + str( cpur ) + " --model LG+G8+F --msa", X_cln_aln_file, "--prefix", prefix ] ), shell=True )
            p.communicate()
    else:    #FastTree2
        tree_file = X_output_folder  + Xname + Xnumber + "FastTree_tree.tre"
        if not os.path.isfile( tree_file ):
            p = subprocess.Popen( args= " ".join( [ fasttree, "-wag -nopr <", X_cln_aln_file, ">", tree_file ] ), shell=True )
            p.communicate()
    return tree_file


def parallel_tree_constructor( num_process_candidates, tree_output_folder,range_candidate_file ,aln_input_file, aln_file, cln_aln_file, X_bait_seq_file, X_cand_file, Xname, Xnumber, mode_aln, mode_tree, mafft, muscle, raxml, fasttree, cpu_max,cpur,parallel_mode):   
    """! @brief handles the construction of alignments and phylogenetic tree; the number of candidates per tree can be limited, so that several trees with a given number of candidates are constructed in parallel 
            @note second FASTA file can be an empty string to run this function just based on one FASTA file
    """
    trees = []
    candidates = load_sequences(X_cand_file)
    candidates_number = min(num_process_candidates,len(candidates))
    num_par =   math.ceil(len(candidates)/candidates_number) 
    cpu_par = math.floor(cpu_max / num_par) 
    
    if not parallel_mode:
        num_par = 1
        candidates_number = len(candidates)
        cpu_par = cpu_max
    elif cpu_par < 8:
        cpu_par = cpur
    
    if not os.path.isdir(tree_output_folder):
        os.mkdir(tree_output_folder)
    
    #extract range
    for i in range(num_par):                                 
        start = i*candidates_number
        end = min((i+1)*candidates_number, len(candidates))
        cand = tree_output_folder + Xname + Xnumber + str(i) + "_" + range_candidate_file
        aln_input = tree_output_folder + Xname + Xnumber + str(i) + "_" +  aln_input_file
        
        if not os.path.isfile( aln_input ):  
            with open(cand,"w") as out:
                for c in list(candidates.keys())[start:end]:
                    out.write( '>' + c + "\n" + candidates[ c ] + "\n" )       
            p = subprocess.Popen( args= "cat " + X_bait_seq_file + " " + cand + " > " + aln_input, shell=True )
            p.communicate()              

    #align 
    processes = []
    cpu_use = 0     
    for i in range(num_par):                                 
        aln_input = tree_output_folder + Xname + Xnumber  + str(i) + "_" +  aln_input_file
        aln = tree_output_folder + Xname + Xnumber + str(i) + "_" +  aln_file
        
        #wait for other processes to finish
        if cpu_use+cpu_par > cpu_max:         
            while cpu_use+cpu_par > cpu_max:
                time.sleep(10)
                cpu_use = len([p for p in processes if p.poll() is None])*cpu_par

        if not os.path.isfile( aln ):      
            if mode_aln == "mafft":
                p = subprocess.Popen( args= mafft + " --quiet --thread " + str(cpu_par) + " " + aln_input + " > " + aln, shell=True )     
                processes.append(p)
            else:
                p = subprocess.Popen(args = muscle + " -align " + aln_input + " -output " + aln  + " -threads " + str(cpu_par), shell=True)
                processes.append(p)                            
            cpu_use += cpu_par 
    [p.wait() for p in processes]# -> all processes parallel

    #tree
    processes = []
    cpu_use = 0 
    for i in range(num_par):                                 
        cln_aln = tree_output_folder + Xname + Xnumber + str(i) + "_" +  cln_aln_file
        aln = tree_output_folder + Xname + Xnumber + str(i) + "_" +  aln_file      
               
        #wait for other processes to finish
        if cpu_use+cpu_par > cpu_max:
            while cpu_use+cpu_par > cpu_max:
                time.sleep(10)
                cpu_use = len([p for p in processes if p.poll() is None])*cpu_par

        if not os.path.isfile( cln_aln ):
            alignment_trimming( aln, cln_aln, occupancy=0.1 )
        
        if mode_tree == "raxml":    #RAxML
            prefix = tree_output_folder + Xname + Xnumber + str(i) + "_" + "RAxML_tree"
            tree_file = prefix + ".raxml.bestTree"
            trees.append(tree_file)
            if not os.path.isfile( tree_file ):
                p = subprocess.Popen( args= " ".join( [ raxml, "--all --threads " + str( cpu_par ) + " --model LG+G8+F --msa",cln_aln, "--prefix", prefix ] ), shell=True )
                processes.append(p)
        else:    #FastTree2
            tree_file = tree_output_folder + Xname + Xnumber + str(i) + "_"  + "FastTree_tree.tre"
            trees.append(tree_file)
            if not os.path.isfile( tree_file ):
                p = subprocess.Popen( args= " ".join( [ fasttree, "-wag -nopr <",cln_aln, ">", tree_file ] ), shell=True )
                processes.append(p)
        cpu_use += cpu_par 
    [p.wait() for p in processes]# -> all processes parallel  
    
    return trees
        

def split_into_ingroup_and_outgroup( tree_file, in_list, out_list, neighbour_cutoff, mean_factor_cutoff, min_neighbour_cutoff, hmmsearch_results ):
    """! @brief split subject sequences into intgroup and outgroup based on reference baits and bait like non-family sequences """

    # --- preparation of data structure --- #
    groups_around_ref_gene = {}
    for gene in ( in_list+out_list ):
        groups_around_ref_gene.update( { gene: [] } )
    
    # --- find node objects of reference genes --- #
    tree = dendropy.Tree.get_from_path( tree_file, "newick" )
    pdm = dendropy.PhylogeneticDistanceMatrix.from_tree( tree )
    my_mean_nearest_taxon_distance = pdm.mean_nearest_taxon_distance()
    
    ref_node_objects = {}
    for node in tree.taxon_namespace:
        try:
            groups_around_ref_gene[ node.label ]
            ref_node_objects.update( { node.label: node } )
        except KeyError:
            pass
    
    ref_gene_nodes = []
    ref_gene_nodes_dict_to_check = {}
    for gene in ( in_list+out_list ):
        ref_gene_nodes.append( ref_node_objects[ gene ] )
        ref_gene_nodes_dict_to_check.update( { ref_node_objects[ gene ]: None } )
    
    results = {}
    for i, t1 in enumerate( tree.taxon_namespace ):
        try:
            ref_gene_nodes_dict_to_check[ t1 ]
        except KeyError:    #only run analysis for non-reference sequences
            path_distances = []
            patristic_distances = {}
            for t2 in tree.taxon_namespace:    #calculate distance to all other sequences in tree
                path_distance = pdm.path_edge_count( t1, t2)
                patr_distance = pdm.patristic_distance( t1, t2 )
                path_distances.append( { "key": t2.label, "val": path_distance } )
                patristic_distances.update( { t2.label: patr_distance } )
            in_counter = 0
            out_counter = 0
            cand_counter = 0
            sorted_distances = sorted( path_distances, key=itemgetter("val") )
            for each in sorted_distances[ : min( [ len( path_distances ), neighbour_cutoff ] ) ]:
                patr = patristic_distances[ each["key"] ]
                if patr < mean_factor_cutoff*my_mean_nearest_taxon_distance:    #exclude outliers on extremely long branches
                    if each["key"] in in_list:    #check if smalles path_distances are to in- or outgroup baits
                        in_counter += 1
                    if each["key"] in out_list:
                        out_counter += 1
                    else:
                        cand_counter += 1
            hmm = '-' if len(hmmsearch_results) == 0 else "yes" if t1.label in hmmsearch_results else "no"                  
            if in_counter+out_counter > min_neighbour_cutoff:
                results.update( { t1.label: { "score": float( in_counter ) / ( in_counter + out_counter ), "in": in_counter, "out": out_counter, 'hmm':hmm, "tree":os.path.basename(tree_file)} } )
            else:
                results.update( { t1.label: { "score": 0.0, "in": in_counter, "out": out_counter, "hmm":hmm, "tree":os.path.basename(tree_file)} } )
            #score ranges from 0 (non-family member) to 1 (family member)
    return results


def analyze_ingroup_outgroup_result(result_folder, tree_files, in_list, out_list, hmmresult, filter_domains):    
    """! @brief analyzes the number of classified ingroup candidates for a variation of parameters """
    
    neighbours =  range(10,100,10)
    scores = [r/10 for r in range(1,11)]
    
    if not os.path.isdir(result_folder):
        os.mkdir(result_folder)
    
    for dist_factor in [5,10]:  
        for min_neighbour_cutoff in range(0,11):
        
            x, y, z = [], [], []     
            neighbour_params = result_folder + "2_classification_neighbourdist_%i_minneighbours_%i.jpg" % (int(dist_factor), int(min_neighbour_cutoff))
            for n in neighbours: 
                
                fam_classification = {}
                for tree_file in tree_files: 
                    classification = split_into_ingroup_and_outgroup( tree_file, in_list, out_list, n, dist_factor, min_neighbour_cutoff, hmmresult)
                    fam_classification.update(classification)
                         
                for s in scores:    
                    results = [r for r in fam_classification if fam_classification[r]["score"]>s  and ( r in hmmresult or not filter_domains)] 
                    x.append(s)    
                    y.append(n)
                    z.append(len(results))
            
            fig = plt.figure()
            ax = plt.axes(projection='3d')
            ax.plot_trisurf(np.array(x), np.array(y),  np.array(z), alpha = 0.8, cmap = 'viridis')
            ax.set_xlabel('score') 
            ax.set_ylabel('neighbours')
            ax.set_zlabel('number of candidates') 
            plt.title("--neighbourdist: %i --minneighbours: % i" % (int(dist_factor), int(min_neighbour_cutoff)))
            plt.tight_layout()
            plt.savefig(neighbour_params,dpi=300)
            plt.close(fig)


def write_tree_symbol(df, field_labels, field_shapes, field_colors, output_file):
    """! @brief creates an iTol dataset """
    with open( output_file ,"w") as out: 
        out.write( "\n".join(["DATASET_BINARY","SEPARATOR COMMA",
         "DATASET_LABEL,label1",
         "COLOR,#ff0000",
         "\n"]))    
        out.write("FIELD_SHAPES,%s\n" % (",".join(field_shapes)))
        out.write("FIELD_LABELS,%s\n" % (",".join(field_labels)))
        out.write("FIELD_COLORS,%s\n" % (",".join(field_colors)))
        out.write("DATA\n" )    
        out.write(df.to_csv(index=False,header=False,sep=","))


def fam_classification_dataset_file( tree_dataset_file, aln_hmmsearch_results, in_list, out_list, fam_classification , min_score_cutoff, candidates_domain_filter):                    
    """! @brief creates an iTol dataset to visualize ingroup and outgroup classification """
    row_list=[]
    
    for bait in in_list:
        missing_domain = "-1" if (bait in aln_hmmsearch_results or len(aln_hmmsearch_results)==0) else "1"
        row_list.append( {"CleanID":bait,"missing_hmm_domain":missing_domain, "bhlh":1})
        
    for bait in out_list:
        missing_domain = "-1" if (bait in aln_hmmsearch_results or len(aln_hmmsearch_results)==0) else "1"
        row_list.append( {"CleanID":bait,"bait_out":"1", "missing_hmm_domain":missing_domain, "bhlh":-1})
        
    for candidate in fam_classification:
        missing_domain = "-1" if (candidate in aln_hmmsearch_results or len(aln_hmmsearch_results)==0) else "1"
        is_ingroup = fam_classification[ candidate ]["score"] > min_score_cutoff and ( candidate in aln_hmmsearch_results or not candidates_domain_filter)
        cnd_in = "1" if is_ingroup else "-1"
        
        if is_ingroup:
            row_list.append({"CleanID":candidate, "cnd_in":cnd_in ,"missing_hmm_domain":missing_domain,  "bhlh":-1})                              

    data = pd.DataFrame(row_list).fillna(value="-1").reindex(columns=["CleanID","missing_hmm_domain","bait_out","cnd_in","bhlh"])
    write_tree_symbol(data,data.columns.values[1:], ["3","1","2","2"],["#3d85c6","#bcbcbc","#f44336","	#00ff00"],tree_dataset_file)
 
    return tree_dataset_file                    


def member_group_assignment( ref_members, tree_file, member_candidates ):
    """! @brief assign new candidates to reference members e.g. the A.thaliana family """
    
    new2ref_mapping_table, new_per_ref_mem = {}, {}
    
    #new2ref_mapping_table = { candiate1: ref1, candidate2: ref1, candiate3: ref14, ... }
    #new_per_ref_mem = { ref1: [ candidate1, candidate2 ], ref2: [candidate15], ref3: [], ... }
    
    # --- preparation of data structure --- #
    my_ref_members = list( sorted( ref_members.keys() ) )
    for gene in my_ref_members:    #reference members
        new_per_ref_mem.update( { gene: [] } )
    
    for gene in member_candidates:    #candidate genes of new species
        new2ref_mapping_table.update( { gene: None } )
    
    # --- find node objects of reference genes --- #
    tree = dendropy.Tree.get_from_path( tree_file, "newick" )
    pdm = dendropy.PhylogeneticDistanceMatrix.from_tree( tree )
    
    ref_node_objects = {}
    new_node_objects = {}
    for node in tree.taxon_namespace:
        try:
            new_per_ref_mem[ node.label ]
            ref_node_objects.update( { node.label: node } )
        except KeyError:
            try:
                new2ref_mapping_table[ node.label ]
                new_node_objects.update( { node.label: node } )
            except KeyError:
                pass
    
    ref_gene_nodes = []
    ref_gene_nodes_dict_to_check = {}
    candidate_gene_nodes = []
    canidate_gene_nodes_dict_to_check = {}
    for gene in my_ref_members:
        ref_gene_nodes.append( ref_node_objects[ gene ] )
        ref_gene_nodes_dict_to_check.update( { ref_node_objects[ gene ]: None } )
    for gene in member_candidates:
        candidate_gene_nodes.append( new_node_objects[ gene ] )
        canidate_gene_nodes_dict_to_check.update( { new_node_objects[ gene ]: None } )
    
    for i, t1 in enumerate( candidate_gene_nodes ):
        edge_distances = []
        patr_distances = []
        for t2 in ref_gene_nodes:    #calculate distance to all other sequences in tree
            edge_distances.append( pdm.path_edge_count( t1, t2) )
            patr_distances.append( pdm.patristic_distance( t1, t2 ) )
        ref_member = my_ref_members[ edge_distances.index( min( edge_distances ) ) ]
        new2ref_mapping_table[ t1.label ] = { 'label': ref_member, 'edges': min( edge_distances ), 'patr': patr_distances[ edge_distances.index( min( edge_distances ) ) ] }
        new_per_ref_mem[ ref_member ].append( t1.label )

    return new2ref_mapping_table, new_per_ref_mem


def motif_check(  sequences_file, motifs_file, hmmresult, hmmoutput, hmmsearch ):
    """! @brief screen sequences for motifs """
    
    # get hmmcheck results
    p = subprocess.Popen( args=  hmmsearch + " --domtblout " + hmmresult + " -o " + hmmoutput + " " + motifs_file + " " + sequences_file , shell=True )
    p.communicate()
    hmmcheck = pd.read_csv(hmmresult, delim_whitespace=True,comment="#",header=None)
    
    # highest c-Evalue for each sequence/domain pair 
    hmmcheck = hmmcheck[hmmcheck[11]<0.00001].sort_values(11).drop_duplicates([0,3])
    
    motifs = list(set(hmmcheck[3]))
    seqs = load_sequences(sequences_file)    
    results = {}
    for key in sorted( seqs.keys() ):       
        results.update( { key: {} } )
        hits = hmmcheck[hmmcheck[0]==key]
        for ID in motifs:
            try:
                coords = hits[hits[3]==ID].values[0]
                match = seqs[key][coords[17]-1:coords[18]]
                results[ key ].update( { str(ID): match } )          
            except IndexError:
                results[ key ].update( { str(ID): "" } )
    
    motifs = [str(m) for m in motifs]    
    return results, motifs


def establish_paralog_groups( tree_file, member_candidates, dist_cutoff_factorB ):
    """! @brief construct paralog groups """
    
    candidate_mapping_table = {}
    for gene in member_candidates:    #candidate genes of new species
        candidate_mapping_table.update( { gene: None } )
    
    # --- find node objects of reference genes --- #
    tree = dendropy.Tree.get_from_path( tree_file, "newick" )
    pdm = dendropy.PhylogeneticDistanceMatrix.from_tree( tree )
    my_mean_nearest_taxon_distance = pdm.mean_nearest_taxon_distance()
    
    new_node_objects = {}    #get new family candidate node objects
    for node in tree.taxon_namespace:
        try:
            candidate_mapping_table[ node.label ]
            new_node_objects.update( { node.label: node } )
        except KeyError:
            pass
    
    candidate_gene_nodes = [] # nodes of new candidates
    for gene in member_candidates:
        candidate_gene_nodes.append( new_node_objects[ gene ] )
    
    black_list = {}
    paralog_collection = []
    for i, t1 in enumerate( candidate_gene_nodes ):
        try:
            black_list[ t1.label ]
        except KeyError:
            paralogs = [ t1.label ]
            edge_distances = []
            patr_distances = {}
            for t2 in tree.taxon_namespace:    #calculate distance to all other sequences in tree
                try:
                    black_list[ t2.label ]
                except KeyError:
                    if t1.label != t2.label:
                        edge_distances.append( { 'id': t2.label, 'dist': pdm.path_edge_count( t1, t2) } )
                        patr_distances.update( { t2.label: pdm.patristic_distance( t1, t2 ) } )
            for each in list( sorted( edge_distances, key=itemgetter('dist') ) ):
                try:
                    candidate_mapping_table[ each['id'] ]
                    if patr_distances[ each['id'] ] < ( my_mean_nearest_taxon_distance*dist_cutoff_factorB ):
                        paralogs.append( each['id'] )
                        black_list.update( { each['id']: None } )
                except KeyError:
                    break    #next neighbour is not a new candidate => break extension of paralog group
            paralog_collection.append( paralogs )
            black_list.update( { t1.label: None } )

    return paralog_collection


def get_represenative_paralog_per_group( paralog_groups, clean_members, repr_clean_file ):
    """! @brief select longest sequence as representative per paralog group """
    
    paralog_representatives = {}
    with open( repr_clean_file, "w" ) as out:
        for group in paralog_groups:
            if len( group ) == 1:
                out.write( '>' + group[0] + "\n" + clean_members[ group[0] ] + "\n" )
                paralog_representatives.update( { group[0]:  clean_members[ group[0] ] } )
            else:
                seqs_len_sorting = []
                for each in group:
                    seqs_len_sorting.append( { 'id': each, 'len': len( clean_members[ each ] ), 'seq': clean_members[ each ] } )
                representative = list( sorted( seqs_len_sorting, key=itemgetter('len', 'id') ) )[-1]
                out.write( '>' + representative['id'] + "\n" + representative['seq'] + "\n" )
                paralog_representatives.update( { representative['id']:  representative['seq'] } )
    return paralog_representatives


def load_candidate_mem_to_mem_mapping_table( new_2_ref_mem_mapping_file ):
    """! @brief load mapping table """
    
    new2ref_mapping_table = {}
    with open( new_2_ref_mem_mapping_file, "r" ) as f:
        f.readline()    #remove header
        line = f.readline()
        while line:
            parts = line.strip().split('\t')
            new2ref_mapping_table.update( { parts[0]: parts[1] } )
            line = f.readline()
    return new2ref_mapping_table


def construct_fasta_file_w_repr_and_aths( ref_mems, new2ref_mapping_table, repr_and_ath_for_tree, repr_and_ath_fasta_file ):
    """! @brief rename sequences with group for final tree construction """
    
    mem2group = {}
    mem_id2name = {}
    mem_name2id = {}
    for each in list( ref_mems.values() ):
        mem2group.update( { each['id']: each['group'] } )
        mem_id2name.update( { each['id']: each['name'] } )
        mem_name2id.update( { each['name']: each['id'] } )
    
    with open( repr_and_ath_fasta_file, "w" ) as out:
        for key in list( repr_and_ath_for_tree.keys() ):
            try:
                group = mem2group[ key ]
            except KeyError:
                try:
                    group = mem2group[ new2ref_mapping_table[ key ] ]
                except KeyError:
                    try:
                        group = mem2group[ mem_name2id[ new2ref_mapping_table[ key ] ] ]
                    except KeyError:
                        group = ""
            try:
                out.write( '>' + mem_id2name[ key ] + "-" + group + "\n" + repr_and_ath_for_tree[ key ] + "\n" )
            except KeyError:
                out.write( '>' + key + "-" + group + "\n" + repr_and_ath_for_tree[ key ] + "\n" )


def analyse_bHLH_domain(binding_properties_file, domain_file,fin_aln_file, member,subject_name_mapping_table):
    """! @brief analyses dna binding properties and extracts bHLH domain based on reference bHLH AtbHLH008 """
    
    ref_name = "AtbHLH008"
    #ref_name = "AT1G09530"
    residues = [347, 351, 354,355,339,398]  # H-E-R motif + first & last domain residue (indexes)
    alignment = load_sequences(fin_aln_file)
    
    # get alignment positions of residues
    refs = [a for a in alignment if ref_name in a]
    residues_al=[-1,-1,-1,-1,-1,-1]
    if len(refs)>0:
        ref_name = refs[0]
        seq_idx = -1  
        for idx, aa in enumerate( alignment[ ref_name ] ):
            if aa != "-":
                seq_idx += 1           
                if seq_idx in residues:               
                    residues_al[residues.index(seq_idx)] = idx         
            if not -1 in residues_al:
                break
        
        # analyze dna binding properties
        basic = ["K","R","H"] 
        h = residues_al[0] # pos9 
        e = residues_al[1] # pos13
        r1 = residues_al[2] # pos16
        r2 = residues_al[3] # pos17      
        start = residues_al[4] # first position: pos1
        end = residues_al[5] # last position: pos60 
        
        with open(binding_properties_file, "w") as out:
            with open(domain_file, "w") as out_domain:          
                out.write("OriginalID\tCleanID\tnum basic residues\tdna binding\te-box binding\tg-box binding\tH-E-R motif\n") 
                for mem in member:               
                    basic_region = alignment[mem][start:r2+1].replace("-","")               
                    num_basic = len([aa for aa in basic_region if aa in basic])                    
                    out.write( "\t".join( list( map( str, [subject_name_mapping_table[ mem ],
                        mem,
                        num_basic,
                        int(num_basic >= 5),
                        int(alignment[mem][e] == 'E' and alignment[mem][r1]=='R'),
                        int(alignment[mem][h] == 'H' and alignment[mem][e]=='E' and alignment[mem][r2] == 'R'),
                        alignment[mem][h]+alignment[mem][e] + alignment[mem][r2]
                        ] ) ) ) + "\n" )
                    out_domain.write(">" + mem + "\n" + alignment[mem][start:end+1].replace("-","") + "\n")
                             
    return binding_properties_file


def get_from_definition_file(clmn,family_definition, fam_definition_file):
    """! @brief gets value from column in definition file; warning, if the file path doesnt exist """
    file = ""
    
    if family_definition is not None:   
        if clmn in family_definition.columns:
            file = str(family_definition.loc[0,clmn])     
            if not os.path.isfile(file):
    
               if os.path.isfile("data/" + file):
                    file = "data/" + file
               else:
                    if len(file) > 1 :
                        sys.stdout.write( 'WARNING: "' + clmn + '" file defined in ' + fam_definition_file  + " does not exist: "  + file + " \n")         
                    file = "" 
    return file

def get_from_arguments_file(arguments,argument_name):
    """! @brief gets file path defined in arguments; warning, if the file path doesnt exist """
    file = arguments[ arguments.index(argument_name)+1 ]  
    if not os.path.isfile(file):
       if os.path.isfile("data/" + file):
            file = "data/" + file
       else:
            if len(file) > 1 :
                sys.stdout.write( "WARNING: File defined in " + argument_name  + " argument does not exist: "  + file + " \n")         
            file = "" 
    return file

#%%


def main( arguments ):
    """! @brief run everything """
    fam = "bHLH"    
#%%    
    
    output_folder = arguments[ arguments.index("--out")+1 ]
    if output_folder[-1] != "/":
        output_folder += "/"
    if not os.path.exists( output_folder ):
        os.makedirs( output_folder )
                     
    if "--subject" in arguments:
        raw_subject_files = [ arguments[ arguments.index("--subject")+1 ] ]
    else:
        subject_file_dir = arguments[ arguments.index("--subjectdir")+1 ]
        if not subject_file_dir[-1] == "/":
            subject_file_dir + "/"
        extensions = [ ".fasta", ".fa", ".fas", ".FASTA", ".FA", ".FAS", ".fna", ".FNA", ".cds", ".CDS", ".pep", ".PEP" ]
        raw_subject_files = [ ]
        for each in extensions:
            raw_subject_files += glob.glob( subject_file_dir + "*" + each )
        raw_subject_files = list( sorted( raw_subject_files ) )
   
    # --- Search options --- # 
    if "--search" in arguments:
        search = arguments[ arguments.index("--search")+1 ]
        if search not in [ "blast", "hmmer" ]:
            search = "blast"
    else:
        search = "blast"   
    if "--mode_aln" in arguments:
        mode_aln = arguments[ arguments.index("--mode_aln")+1 ]
        if mode_aln not in [ "mafft", "muscle" ]:
            mode_aln = "muscle"
    else:
        mode_aln = "muscle"    
    if "--mode_tree" in arguments:
        mode_tree = arguments[ arguments.index("--mode_tree")+1 ]
        if mode_tree not in [ "fasttree", "raxml" ]:
            mode_tree = "fasttree"
    else:
        mode_tree = "fasttree"
       
    # --- Execution Options --- #
    if "--name" in arguments:
        name = arguments[ arguments.index("--name")+1 ]
    else:
        name = ""
    if "--cdsinput" in arguments:
        cds_input = True
    else:
        cds_input = False   
    if "--keepnames" in arguments:
        trim_names = False
    else:
        trim_names = True
    if "--collapse" in arguments:
        collapse_mode = True
    else:
        collapse_mode = False
    if "--parameter_graphs" in arguments:
        parameter_graphs = True
    else:
        parameter_graphs = False

    # --- CPU usage --- #
    if "--cpu" in arguments:
        cpu = int( arguments[ arguments.index("--cpu")+1 ] )
    else:
        cpu = 4 
    if "--cpumax" in arguments:
        cpu_max = int( arguments[ arguments.index("--cpumax")+1 ] )
    else:
        cpu_max = cpu   
    if "--cpub" in arguments:
        cpub = int( arguments[ arguments.index("--cpub")+1 ] )
    else:
        cpub = cpu + 0
    if "--cpur" in arguments:
        cpur = int( arguments[ arguments.index("--cpur")+1 ] )
    else:
        cpur = cpu + 0
    
    # --- Tool paths --- #
    if "--blastp" in arguments:
        blastp = arguments[ arguments.index("--blastp")+1 ]
    else:
        blastp = "blastp"
    if "--makeblastdb" in arguments:
        makeblastdb = arguments[ arguments.index("--makeblastdb")+1 ]
    else:
        makeblastdb = "makeblastdb"
    if "--hmmsearch" in arguments:
        hmmsearch = arguments[ arguments.index("--hmmsearch")+1 ]
    else:
        hmmsearch = "hmmsearch"
    if "--mafft" in arguments:
        mafft = arguments[ arguments.index("--mafft")+1 ]
    else:
        mafft = "mafft"
    if "--muscle" in arguments:
        muscle = arguments[ arguments.index("--muscle")+1 ]
    else:
        muscle = "muscle"
    if "--fasttree" in arguments:
        fasttree = arguments[ arguments.index("--fasttree")+1 ]
    else:
        fasttree = "fasttree"
    if "--raxml" in arguments:
        raxml = arguments[ arguments.index("--raxml")+1 ]
    else:
        raxml = "raxml-ng"
       
    # --- BLAST hit cutoffs --- #
    if "--simcutp" in arguments:
        similarity_cutoff_p = float( arguments[ arguments.index("--simcutp")+1 ] )
    else:
        similarity_cutoff_p=40.00
    if "--poscutp" in arguments:
        possibility_cutoff_p = int( arguments[ arguments.index("--poscutp")+1 ] )
    else:
        possibility_cutoff_p=100
    if "--lencutp" in arguments:
        length_cutoff_p = int( arguments[ arguments.index("--lencutp")+1 ] )
    else:
        length_cutoff_p= 80      
    if "--bitcutp" in arguments:
        bitscore_p = int( arguments[ arguments.index("--bitcutp")+1 ] )
    else:
        bitscore_p=60
    
    # --- Candidate cutoffs --- #
    if "--filterdomain" in arguments:
        candidates_domain_filter = True
    else:
        candidates_domain_filter = False 
    if "--minscore" in arguments:
        min_score_cutoff = float( arguments[ arguments.index("--minscore")+1 ] )
    else:
        min_score_cutoff=0.5    #minimal score for candidate to be considered as ingroup      

    if "--numneighbours" in arguments:
        neighbour_cutoff = int( arguments[ arguments.index("--numneighbours")+1 ] )
    else:
        neighbour_cutoff=10    #numbers of closest neighbour that is considered in ingroup/outgroup classification   
    if "--neighbourdist" in arguments:
        mean_factor_cutoff = float( arguments[ arguments.index("--neighbourdist")+1 ] )
    else:
        mean_factor_cutoff=5.0    #X*average nearest neighbor distance  
    if "--minneighbours" in arguments:
        min_neighbour_cutoff = int( arguments[ arguments.index("--minneighbours")+1 ] )
    else:
        min_neighbour_cutoff = 0    #minimal number of valid bait sequences (ingroup+outgroup) in range - 1       
    if "--parallel" in arguments:
        parallel_mode = True
    else:
        parallel_mode = False
    if "--numprocesscandidates" in arguments:
        num_process_candidates = int( arguments[ arguments.index("--numprocesscandidates")+1 ] )
    else:
        num_process_candidates = 200    #max number of candidates per ingroup/outgroup classification  
    if "--paralogdist" in arguments:
        dist_cutoff_factorB = float( arguments[ arguments.index("--paralogdist")+1 ] )
    else:
        dist_cutoff_factorB=10.0    #X*average nearest neighbour distance used as cutoff in the monophyletic tip masking


    # --- get inputs from definition file --- # 
    if "--info" in arguments:
        fam_definition_file = arguments[ arguments.index("--info")+1 ] 
        family_definition = pd.read_csv(fam_definition_file,sep="\t")
        family_definition = family_definition.fillna('')
        
        if len(family_definition.index) < 1: 
            sys.exit( "ERROR: No files defined for " + fam + " family in " + fam_definition_file + "\n" )
            
        if not (("Baits" in family_definition.columns or "--baits" in arguments) and ("BaitsInfo" in family_definition.columns or "--baitsinfo" in arguments)):    
            sys.exit( 'ERROR: Mandatory columns "Baits" and/or "BaitsInfo" are missing in ' + fam_definition_file + "\n" )        
    else:
        fam_definition_file = ""
        family_definition = None
    
        
    #Baits
    if "--baits" in arguments:
        fam_bait_seq_file_all = get_from_arguments_file(arguments, "--baits")
    else: 
        fam_bait_seq_file_all = get_from_definition_file("Baits", family_definition, fam_definition_file)
    if fam_bait_seq_file_all == "":
        sys.exit( 'ERROR: "Baits" file must be defined for execution of annotator \n')  
        
    #BaitsInfo
    if "--baitsinfo" in arguments:
        fam_info_file = get_from_arguments_file(arguments, "--baitsinfo")
    else:  
        fam_info_file = get_from_definition_file("BaitsInfo", family_definition, fam_definition_file)
    if fam_info_file == "":
        sys.exit( 'ERROR: "BaitsInfo" file must be defined for execution of annotator \n') 
    
    #ThinnedBaits
    if "--optimisedbaits" in arguments:
        fam_bait_seq_file  = get_from_arguments_file(arguments, "--optimisedbaits")
    else: 
        fam_bait_seq_file  = get_from_definition_file("OptimisedBaits", family_definition, fam_definition_file)   
    if fam_bait_seq_file  == "":
        fam_bait_seq_file = fam_bait_seq_file_all
    
    #HMM
    if "--hmm" in arguments:
        fam_bait_hmm = get_from_arguments_file(arguments, "--hmm")
    else:
        fam_bait_hmm = get_from_definition_file("HMM", family_definition, fam_definition_file) #HMM   
    
    #Reference
    if "--reference" in arguments:
        ref_file = get_from_arguments_file(arguments, "--reference")
    else:
        ref_file = get_from_definition_file("Reference", family_definition, fam_definition_file) #Reference
    
    #Ath
    if "--ath" in arguments:
        fam_ath_file = get_from_arguments_file(arguments, "--ath")
    else:
        fam_ath_file = get_from_definition_file("Ath", family_definition, fam_definition_file) #Ath  
    
    #Motif
    if "--motifs" in arguments:
        motif_file = get_from_arguments_file(arguments, "--motifs")
    else:    
        motif_file = get_from_definition_file("Motifs", family_definition, fam_definition_file) #Motifs
        
    #Landmark
    if "--landmark" in arguments:
        raw_landmark_file = get_from_arguments_file(arguments, "--landmark")
    else:    
        raw_landmark_file = get_from_definition_file("Landmark", family_definition, fam_definition_file) #Landmark
    
    
    # Check hmm file for domain filter and hmmer search       
    if (search != "blast" or candidates_domain_filter) and not os.path.isfile( fam_bait_hmm ):            
        if search != "blast":
            sys.stdout.write( 'ERROR: Option "hmmer" for --search only possible if HMM motif defined for ' + fam + " family" + "\n")                   
        if candidates_domain_filter :
            sys.stdout.write( "ERROR: Option --filterdomain only possible if HMM motif defined for " + fam + " family" + "\n")                                    
        sys.exit()

#%% 
  
    # --- separated analyses for each subject --- #             
    for jidx, raw_subject_file in enumerate( raw_subject_files ):    #use jidx to generate unique IDs for all jobs
        
        # --- prepare output folder for each job if there are multiple --- #   
        if len( raw_subject_files ) == 1:
            job_output_folder = output_folder 
        else:
            job_ID = raw_subject_file.split("/")[-1].split(".")[0]
            job_output_folder = output_folder +"/" + str( jidx ).zfill(5) + "_" + job_ID + "/"
        
        if not job_output_folder[-1] +"/":
            job_output_folder = job_output_folder + "/" 
        
        if not os.path.exists( job_output_folder ):
            os.makedirs( job_output_folder )
          
                   
        # --- 00 validation of inputs --- #
        #check if all baits are listed in the info file
        subject_file = job_output_folder + "clean_subject_sequences.pep.fasta"
        mapping_table_file = job_output_folder + "raw_subject_to_clean_subject_mapping_table.txt"
        if not os.path.isfile( subject_file ):
            clean_input_FASTA_file( raw_subject_file, subject_file, mapping_table_file, cds_input, trim_names )    #remove illegal characters from subject sequence headers
        landmark_file = job_output_folder + "clean_landmark_flavonoid_bHLHs.pep.fasta"
        landmark_mapping_table_file = job_output_folder +  "raw_landmark_to_clean_landmark_mapping_table.txt"
        if os.path.isfile(raw_landmark_file) and not os.path.isfile(landmark_file):
            clean_input_FASTA_file( raw_landmark_file, landmark_file,landmark_mapping_table_file, cds_input, trim_names )
               
        fam_check_status = check_fam_IDs_across_files( fam_bait_seq_file_all, fam_bait_seq_file ,fam_info_file, fam_ath_file,ref_file, fam )
        if not fam_check_status:
            sys.exit( "ERROR: analysis is stopped due to inconstistency of " + fam +" IDs between files" )
        subject_name_mapping_table = load_subject_name_mapping_table( mapping_table_file )
                    
        result_folder = job_output_folder + "RESULTS/"
        if not os.path.exists( result_folder ):
            os.makedirs( result_folder )           
        doc_file = result_folder + name + "00_documentation.txt"           
        generate_documentation_file(doc_file, fam_bait_seq_file_all, fam_bait_seq_file, fam_info_file, fam_bait_hmm, fam_ath_file, motif_file, job_output_folder, raw_subject_file,
                                    search, mode_aln ,mode_tree, blastp, makeblastdb,hmmsearch, cpu_max ,cpub, cpur, mafft, muscle,raxml, fasttree, ref_file,
                                    bitscore_p, similarity_cutoff_p, possibility_cutoff_p, length_cutoff_p, cds_input,
                                    min_score_cutoff,neighbour_cutoff, mean_factor_cutoff, min_neighbour_cutoff, dist_cutoff_factorB, fam,
                                    candidates_domain_filter, parallel_mode, num_process_candidates, name, trim_names, collapse_mode, parameter_graphs)        
        start_time = datetime.datetime.now()   
        
        
        # --- 01 find initial candidates --- #
        seq_search_result_file = job_output_folder +"01_seq_search_results.txt"
        blast_analyze_folder = job_output_folder +"01_blast_seq_search_analysis/"
        if not os.path.isfile( seq_search_result_file ):
            if search == "blast":
                blast_db = job_output_folder + "01_blastdb"
                p = subprocess.Popen( args= makeblastdb + " -in " + subject_file + " -out " + blast_db + " -dbtype prot", shell=True )
                p.communicate()
                
                p = subprocess.Popen( args= blastp + " -query " + fam_bait_seq_file_all + " -db " + blast_db + " -out " + seq_search_result_file + " -outfmt 6 -evalue 0.001 -num_threads " + str( cpub ), shell=True )
                p.communicate()
            else:                          
                p = subprocess.Popen( args= hmmsearch + " --tblout " + seq_search_result_file + " " + fam_bait_hmm + " " + subject_file + " > " + job_output_folder+"01_hmmsearch_waste.txt", shell=True )
                p.communicate()
        
        if search == "blast":
            seq_search_results = load_BLAST_results( seq_search_result_file, similarity_cutoff_p, possibility_cutoff_p, length_cutoff_p ,bitscore_p)    #load valid BLASTp results
            if parameter_graphs:
                analyze_blast_result(blast_analyze_folder, seq_search_result_file)            
        else:
            seq_search_results = load_hmmsearch_results( seq_search_result_file )    #load valid hmmsearch results
        
        subject_sequences = load_sequences( subject_file )
        if cds_input:
            cds_subject_sequences = load_sequences( subject_file.replace( ".pep.fasta", ".cds.fasta" ) )
        
        candidate_file = result_folder + name + "01_initial_candidates.pep.fasta"
        with open( candidate_file, "w" ) as out:
            for each in seq_search_results.keys():      
                out.write( '>' + each + "\n" + subject_sequences[ each ] + "\n" )
                   
        if cds_input:
            cds_candidate_file = candidate_file.replace( ".pep.fasta", ".cds.fasta" )
            with open( cds_candidate_file, "w" ) as out:
                for each in seq_search_results.keys():
                    out.write( '>' + each + "\n" + cds_subject_sequences[ each ] + "\n" )
    
        
        # --- 02 construct phylogenetic tree and analyze tree file --- #                
        tree_output_folder = job_output_folder + "02_in_out_%s_analysis_trees/" % (fam)
        in_list, out_list = load_bait_fam_anno( fam_info_file, fam_bait_seq_file )
        range_candidate_file =  "alignment_candidates.fasta"              
        aln_input_file =  "alignment_input.fasta" 
        aln_file =  "alignment_input.fasta.aln" 
        cln_aln_file =  "alignment_input.fasta.aln.cln"  
        
        # hmm motif            
        hmmsearch_seq_file = job_output_folder +"02_hmmsearch_sequences.fasta"
        aln_hmmsearch_results_file = job_output_folder +"02_hmmsearch_results.txt"
        aln_hmmsearch_results = {}
        if os.path.isfile( fam_bait_hmm ) and not os.path.isfile( aln_hmmsearch_results_file  ) :                  
            p = subprocess.Popen( args= "cat " + fam_bait_seq_file + " " +  candidate_file + " > " + hmmsearch_seq_file, shell=True )
            p.communicate()                 
            p = subprocess.Popen( args= hmmsearch + " --domtblout " + aln_hmmsearch_results_file + " " + fam_bait_hmm + " " + hmmsearch_seq_file + " > " + job_output_folder + "02_hmmsearch_waste.txt", shell=True )
            p.communicate()
        if os.path.isfile( aln_hmmsearch_results_file  ) :
            aln_hmmsearch_results = load_hmmsearch_results( aln_hmmsearch_results_file )
   
        # --- 02 first classification: --- #
        clean_members_file_f =  tree_output_folder + name + "first_clean_%ss.pep.fasta" % (fam)
        tmp_result_table_f =  tree_output_folder + name + "first_in_out_%s_analysis_results.txt" % (fam)
        if not os.path.isfile( tmp_result_table_f ):              
            
            fam_classification = {}
            sys.stdout.write( "Number of ingroup %s baits: " % (fam) + str( len( in_list ) ) + "\n" )
            sys.stdout.write( "Number of outgroup %s baits: " % (fam) + str( len( out_list ) ) + "\n" )
            sys.stdout.flush()                            
                                             
            # construct phylogenetic trees       
            tree_files = parallel_tree_constructor( num_process_candidates, tree_output_folder,range_candidate_file ,aln_input_file, aln_file, cln_aln_file, fam_bait_seq_file, candidate_file,"first_", "", mode_aln, mode_tree, mafft, muscle, raxml, fasttree, cpu_max,cpur, parallel_mode )
            for tree_file in tree_files: 
                classification = split_into_ingroup_and_outgroup( tree_file, in_list, out_list, neighbour_cutoff, mean_factor_cutoff, min_neighbour_cutoff, aln_hmmsearch_results)
                fam_classification.update(classification)
                           
            with open( clean_members_file_f, "w" ) as out:
                with open( tmp_result_table_f, "w" ) as out2:
                    out2.write( "ClassificationTree\tOriginalID\tCleanID\tHMMmotif\tScore\tIngroupMatches\tOutgroupMatches\n" )
                    candidate_order = list( sorted( fam_classification.keys() ) )
                    for candidate in candidate_order:
                        out2.write( "\t".join( list( map( str,  [fam_classification[ candidate ]["tree"],
                                                                    subject_name_mapping_table[ candidate ],
                                                                    candidate,
                                                                    fam_classification[ candidate ]["hmm"],
                                                                    fam_classification[ candidate ]["score"],
                                                                    fam_classification[ candidate ]["in"],
                                                                    fam_classification[ candidate ]["out"],
                                                                    ] ) ) ) + "\n" )
                        if fam_classification[ candidate ]["score"] > min_score_cutoff and ( candidate in aln_hmmsearch_results or not candidates_domain_filter):                            
                            out.write( ">" + candidate + "\n" + subject_sequences[ candidate ] + "\n" )
            
            # create analysis of candidate number for parameter variation
            classification_analyze_folder = job_output_folder +"02_first_in_out_analysis/"                
            if parameter_graphs:
                analyze_ingroup_outgroup_result(classification_analyze_folder, tree_files, in_list, out_list, aln_hmmsearch_results, candidates_domain_filter)                 
            
            # create tree dataset file in iTOL format
            tree_dataset_file = tree_output_folder + "first_tree_info_dataset.txt"               
            fam_classification_dataset_file( tree_dataset_file, aln_hmmsearch_results, in_list, out_list, fam_classification , min_score_cutoff, candidates_domain_filter)
                    
        # --- 02 second classification: 
        clean_members_file_s =  tree_output_folder + name + "second_clean_%ss.pep.fasta" % (fam)
        tmp_result_table_s =  tree_output_folder + name + "second_in_out_%s_analysis_results.txt" % (fam)
        clean_members_file = result_folder + name + "02_clean_%ss.pep.fasta" % (fam)
        tmp_result_table = result_folder + name + "02_in_out_%s_analysis_results.txt" % (fam)
        if not os.path.isfile( tmp_result_table ):
        
            fam_classification = {}
            sys.stdout.write( "Number of ingroup %s baits: " % (fam) + str( len( in_list ) ) + "\n" )
            sys.stdout.write( "Number of outgroup %s baits: " % (fam) + str( len( out_list ) ) + "\n" )
            sys.stdout.flush()                            
                          
            # construct phylogenetic trees       
            tree_files = parallel_tree_constructor( num_process_candidates, tree_output_folder,range_candidate_file ,aln_input_file, aln_file, cln_aln_file, fam_bait_seq_file, clean_members_file_f,"second_", "", mode_aln, mode_tree, mafft, muscle, raxml, fasttree, cpu_max,cpur, parallel_mode )
            for tree_file in tree_files: 
                classification = split_into_ingroup_and_outgroup( tree_file, in_list, out_list, neighbour_cutoff, mean_factor_cutoff, min_neighbour_cutoff, aln_hmmsearch_results)
                fam_classification.update(classification)
                           
            with open( clean_members_file_s, "w" ) as out:
                with open( tmp_result_table_s, "w" ) as out2:
                    out2.write( "ClassificationTree\tOriginalID\tCleanID\tHMMmotif\tScore\tIngroupMatches\tOutgroupMatches\n" )
                    candidate_order = list( sorted( fam_classification.keys() ) )
                    for candidate in candidate_order:
                        out2.write( "\t".join( list( map( str,  [fam_classification[ candidate ]["tree"],
                                                                    subject_name_mapping_table[ candidate ],
                                                                    candidate,
                                                                    fam_classification[ candidate ]["hmm"],
                                                                    fam_classification[ candidate ]["score"],
                                                                    fam_classification[ candidate ]["in"],
                                                                    fam_classification[ candidate ]["out"],
                                                                    ] ) ) ) + "\n" )
                        if fam_classification[ candidate ]["score"] > min_score_cutoff and ( candidate in aln_hmmsearch_results or not candidates_domain_filter):                            
                            out.write( ">" + candidate + "\n" + subject_sequences[ candidate ] + "\n" )

            # create analysis of candidate number for parameter variation
            classification_analyze_folder = job_output_folder +"02_second_in_out_analysis/"                
            if parameter_graphs:
                analyze_ingroup_outgroup_result(classification_analyze_folder, tree_files, in_list, out_list, aln_hmmsearch_results, candidates_domain_filter) 

            # create tree dataset file in iTOL format
            tree_dataset_file = tree_output_folder + "second_tree_info_dataset.txt"               
            fam_classification_dataset_file( tree_dataset_file, aln_hmmsearch_results, in_list, out_list, fam_classification , min_score_cutoff, candidates_domain_filter)
            
            # final candidate file
            p = subprocess.Popen( args= "cp " +  clean_members_file_s + " " +  clean_members_file, shell=True )         
            p = subprocess.Popen( args= "cp " +  tmp_result_table_s + " " +  tmp_result_table, shell=True )
            p.communicate()   
            if cds_input:    #generate corresponding CDS clean member output file
                cds_clean_members_file = clean_members_file.replace( ".pep.fasta", ".cds.fasta" )
                with open( cds_clean_members_file, "w" ) as out:
                    for candidate in candidate_order:
                        if fam_classification[ candidate ]["score"] >= min_score_cutoff and ( candidate in aln_hmmsearch_results or not candidates_domain_filter):
                            out.write( ">" + candidate + "\n" + cds_subject_sequences[ candidate ] + "\n" )
                           
        else:
            fam_classification = load_fam_classification_from_file( tmp_result_table )
        clean_members = load_sequences( clean_members_file )
        if len( list( clean_members.keys() ) ) < 1:
            sys.exit( "ERROR: no " + fam + "s detected." )            
        
        
        # --- 03 construct a final tree --- #
        fin_aln_input_file = job_output_folder + "03_fin_alignment_input.fasta"
        fin_aln_file = job_output_folder + "03_fin_alignment_input.fasta.aln"
        fin_cln_aln_file = job_output_folder + "03_fin_alignment_input.fasta.aln.cln"
        tree_file = tree_constructor( fin_aln_input_file, fin_aln_file, fin_cln_aln_file, fam_bait_seq_file, clean_members_file, result_folder, name, "03_final_", mode_aln, mode_tree, mafft, muscle, raxml, fasttree, cpur)
        
        
        # --- 04 find closest reference for baits --- #            
        group_around_ref_file = result_folder + name + "04a_candidates_group_around_%s_bait.txt"  % (fam)   #produce table sorted by baits
        new_2_ref_mapping_file = result_folder + name + "04a_candidate_2_%s_bait_mapping_file.txt"  % (fam)   #produce table sorted by subject sequences
        if not os.path.isfile(group_around_ref_file ):
            fam_bait_seqs = load_sequences( fam_bait_seq_file)       
            if not os.path.isfile( new_2_ref_mapping_file ):
                
                new2ref_mapping_table, new_per_ref_members = member_group_assignment( fam_bait_seqs, tree_file, clean_members.keys() )
        
                with open( group_around_ref_file, "w" ) as out:
                    out.write( "RefMember\tNewMembers\n" )
                    gene_order = list( sorted( new_per_ref_members.keys() ) )
                    for gene in gene_order:
                        original_gene_names = []
                        for x in new_per_ref_members[ gene ]:
                            original_gene_names.append( subject_name_mapping_table[ x ] )
                        out.write( gene + "\t" + ";".join( original_gene_names )  + "\n" )        
                
                with open( new_2_ref_mapping_file, "w" ) as out:
                    out.write( "NewMember\tRefMember\tEdgeDistance\tPatristicDistance\n" )
                    gene_order = list( sorted( new2ref_mapping_table.keys() ) )
                    for gene in gene_order:
                        out.write( "\t".join( list( map( str, [gene,
                                                               new2ref_mapping_table[ gene ]['label'],    #map back to member name
                                                               new2ref_mapping_table[ gene ]['edges'],
                                                               new2ref_mapping_table[ gene ]['patr']
                                                               ] ) ) ) + "\n" )    #label, edges, patr
        
        
        # --- 04 find closest reference in reference file --- #            
        group_around_ref_file = result_folder + name + "04b_candidates_group_around_reference.txt"   #produce table sorted by reference members
        new_2_ref_mapping_file = result_folder + name + "04b_candidate_2_reference_mapping_file.txt"     #produce table sorted by subject sequences
        if len( ref_file ) > 0 and not os.path.isfile(group_around_ref_file ):    #only performed if reference file is provided
            ref_input_file = job_output_folder + "04_reference_sequences.fasta"
            ref_ath_fin_aln_input_file = job_output_folder + "04_reference_ath_fin_alignment_input.fasta"
            ref_ath_fin_aln_file = job_output_folder + "04_reference_ath_fin_alignment_input.fasta.aln"
            ref_ath_fin_cln_aln_file = job_output_folder + "04_reference_ath_fin_alignment_input.fasta.aln.cln"
            
            ref_members = load_ref_members( ref_file )    #load IDs and trivial name from additional text file (dictionary)
            fam_bait_seqs = load_sequences( fam_bait_seq_file_all)
            fam_ath_seqs = load_sequences (fam_ath_file) if len(fam_ath_file) > 0 else {}               
            with open(ref_input_file, "w") as out:
                for ref in ref_members.keys():
                    try:
                        ref_seq = fam_bait_seqs[ref]
                    except KeyError:
                        ref_seq = fam_ath_seqs[ref]
                    out.write(">" + ref + "\n" + ref_seq + "\n" )                
            ref_tree_file = tree_constructor( ref_ath_fin_aln_input_file, ref_ath_fin_aln_file, ref_ath_fin_cln_aln_file, ref_input_file, clean_members_file,job_output_folder, name, "04_reference_", mode_aln, mode_tree, mafft, muscle, raxml, fasttree, cpur )
        
            if not os.path.isfile( new_2_ref_mapping_file ):              
                new2ref_mapping_table, new_per_ref_members = member_group_assignment( ref_members, ref_tree_file, clean_members.keys() )     
                with open( group_around_ref_file, "w" ) as out:
                    out.write( "RefMember\tFunction\tNewMembers\n" )
                    gene_order = list( sorted( new_per_ref_members.keys() ) )
                    for gene in gene_order:
                        original_gene_names = []
                        for x in new_per_ref_members[ gene ]:
                            original_gene_names.append( subject_name_mapping_table[ x ] )
                        out.write( ref_members[ gene ]['name'] + "\t" + ref_members[ gene ]['function'] + "\t" + ";".join( original_gene_names )  + "\n" )                     
                with open( new_2_ref_mapping_file, "w" ) as out:
                    out.write( "NewMember\tRefMember\tEdgeDistance\tPatristicDistance\tAnnotation\n" )
                    gene_order = list( sorted( new2ref_mapping_table.keys() ) )
                    for gene in gene_order:
                        out.write( "\t".join( list( map( str, [gene,
                                                               ref_members[ new2ref_mapping_table[ gene ]['label'] ]['name'],    #map back to member name
                                                               new2ref_mapping_table[ gene ]['edges'],
                                                               new2ref_mapping_table[ gene ]['patr'],
                                                               ref_members[new2ref_mapping_table[ gene ]['label'] ]['function']    #map back to member name
                                                               ] ) ) ) + "\n" )    #label, edges, patr
        

        #--- 05 check for domain and dna binding properties --- #
        binding_properties_file = result_folder + name + "05_bhlh_dna_binding_group.txt" #produce table with analysis of dna binding properties
        domain_check_file = result_folder + name + "05_bhlh_domain_check.txt" #FASTA file containing family domain sequences
        domain_check_file_pep = result_folder + name + "05_bhlh_domains.pep.fasta" #FASTA file containing family domain sequences
        if not os.path.isfile( binding_properties_file ):
            analyse_bHLH_domain(binding_properties_file, domain_check_file_pep, fin_aln_file,clean_members,subject_name_mapping_table)
        
        if (not os.path.isfile(domain_check_file)) and os.path.isfile(aln_hmmsearch_results_file) and len(aln_hmmsearch_results) > 0:
            with open(  domain_check_file_pep, "w" ) as out:
                hmmcheck = pd.read_csv(aln_hmmsearch_results_file, delim_whitespace=True,comment="#",header=None)
                hmmcheck = hmmcheck[hmmcheck[11]<0.00001].sort_values(11).drop_duplicates([0,3])  # highest c-Evalue for each
                for candidate in sorted(clean_members.keys() ):
                    try:
                        best_hit = hmmcheck[hmmcheck[0]==candidate].values[0]
                        match = clean_members[candidate][best_hit[17]-1:best_hit[18]] 
                        out.write(">" + subject_name_mapping_table[candidate] + "\n" + match + "\n" )  
                    except IndexError:
                        continue
            with open( domain_check_file, "w" ) as out:
                  out.write( "OriginalGeneID\tCleanGeneID\tbHLH domain status\n" )
                  candidates = list( sorted( clean_members.keys() ) )
                  for candidate in candidates:
                      domain = "yes" if candidate in aln_hmmsearch_results else "no"
                      out.write(subject_name_mapping_table[candidate] + "\t" + candidate + "\t" + domain + "\n" )
                  

        # --- 06 check for different motifs --- #
        motif_check_file_summary = result_folder + name + "06_motif_check.txt"    #produce table with motif (0/1)
        motif_check_file_seqs = result_folder + name + "06_motif_check_details.txt"    #produce table with motif sequences
        motif_check_hmm_result = job_output_folder + "06_motif_hmmsearch.txt"
        motif_check_hmm_output = job_output_folder + "06_motif_hmmsearch_output.txt"
        if len( motif_file ) > 0 and not os.path.isfile(motif_check_file_summary):
            motif_check_results, motif_names = motif_check( clean_members_file, motif_file, motif_check_hmm_result, motif_check_hmm_output ,hmmsearch )    #prim key = seqID, secondary key = motifs
            with open( motif_check_file_summary, "w" ) as out1:
                out1.write( "RefMember\t"+"\t".join(motif_names)+"\n" )
                with open( motif_check_file_seqs, "w" ) as out2:
                    out2.write( "RefMember\t"+"\t".join(motif_names)+"\n" )
                    candidates = list( sorted( clean_members.keys() ) )
                    for candidate in candidates:
                        new_line_details = [ subject_name_mapping_table[ candidate ] ]
                        new_line_summary = [ subject_name_mapping_table[ candidate ] ]
                        for mot in motif_names:
                            new_line_details.append( motif_check_results[ candidate ][ mot ] )
                            if len( motif_check_results[ candidate ][ mot ] ) > 0:
                                new_line_summary.append( "1" )
                            else:
                                new_line_summary.append( "0" )
                        out1.write("\t".join(new_line_summary) + "\n")
                        out2.write("\t".join(new_line_details) + "\n")
                    
                                
        # --- 07 construct a final tree with Ath members --- #
        if len( fam_ath_file) > 0:
            ath_fin_aln_input_file = job_output_folder + "07_ath_fin_alignment_input.fasta"
            ath_fin_aln_file = job_output_folder + "07_ath_fin_alignment_input.fasta.aln"
            ath_fin_cln_aln_file = job_output_folder + "07_ath_fin_alignment_input.fasta.aln.cln"
            tree_constructor( ath_fin_aln_input_file, ath_fin_aln_file, ath_fin_cln_aln_file, fam_ath_file, clean_members_file, result_folder, name, "07_ath_", mode_aln, mode_tree, mafft, muscle, raxml, fasttree, cpur)
            
            
        # --- 08 find in species-specific paralogs (in-paralogs) --- #
        if collapse_mode:
            paralog_groups = establish_paralog_groups( tree_file, clean_members.keys(), dist_cutoff_factorB )    #get list of sublist; each sublist represents one paralog group
            paralog_group_file = result_folder + name + "08_repr_%ss.txt" % (fam)
            repr_clean_file = result_folder + name + "08_repr_%ss.pep.fasta" % (fam)
            if not os.path.isfile( repr_clean_file ):
                rep_per_group = get_represenative_paralog_per_group( paralog_groups, clean_members, repr_clean_file )
                with open( paralog_group_file, "w" ) as out:
                    out.write( "RepresentativeSeqID\tMembersOfParalogGroup\n" )
                    for gene in list( rep_per_group.keys() ):
                        for group in paralog_groups:
                            if gene in group:
                                out.write( gene + "\t" + ";".join( group ) + "\n" )
                if cds_input:
                    cds_repr_clean_file = repr_clean_file.replace( ".pep.fasta", ".cds.fasta" )
                    with open( cds_repr_clean_file, "w" ) as out:
                        for key in list( rep_per_group.keys() ):
                            out.write( '>' + key + "\n" + cds_subject_sequences[ key ] + "\n" )
                            
            repr_fin_aln_input_file = job_output_folder + "08_repr_fin_alignment_input.fasta"
            repr_fin_aln_file = job_output_folder + "08_repr_fin_alignment_input.fasta.aln"
            repr_fin_cln_aln_file = job_output_folder + "08_repr_fin_alignment_input.fasta.aln.cln"
            tree_constructor( repr_fin_aln_input_file, repr_fin_aln_file, repr_fin_cln_aln_file, fam_bait_seq_file, repr_clean_file, result_folder, name, "08_representatives_baits_", mode_aln, mode_tree, mafft, muscle, raxml, fasttree, cpur)
                        
            # --- 09 represent cluster only by longest sequence --- #
            if len( fam_ath_file ) > 0:
                repr_ath_fin_aln_input_file = job_output_folder + "09a_repr_ath_fin_alignment_input.fasta"
                repr_ath_fin_aln_file = job_output_folder + "09a_repr_ath_fin_alignment_input.fasta.aln"
                repr_ath_fin_cln_aln_file = job_output_folder + "09a_repr_ath_fin_alignment_input.fasta.aln.cln"
                tree_file = tree_constructor( repr_ath_fin_aln_input_file, repr_ath_fin_aln_file, repr_ath_fin_cln_aln_file, fam_ath_file, repr_clean_file,result_folder, name, "09a_representatives_ath_", mode_aln, mode_tree, mafft, muscle, raxml, fasttree, cpur )
                
                # --- define groups for Ath members and include these group names in tip labels of a phylogenetic tree --- #
                if len( ref_file ) > 0:    #only performed if reference member file is provided
                    ref_members = load_ref_members( ref_file ) # id : id name function group
                    new2ref_mapping_table = load_candidate_mem_to_mem_mapping_table( new_2_ref_mapping_file ) # cnd : ref
                    repr_and_ath_for_tree = load_sequences( repr_ath_fin_aln_input_file ) # id : seq of ath and cnds
                    repr_and_ath_fasta_file = job_output_folder + name + "09_repr_ath_%s.fasta" % (fam)
                    construct_fasta_file_w_repr_and_aths( ref_members, new2ref_mapping_table, repr_and_ath_for_tree, repr_and_ath_fasta_file )
                    group_aln_input_file = job_output_folder + "09b_group_alignment_input.fasta"
                    group_aln_file = job_output_folder + "09b_group_alignment_input.fasta.aln"
                    group_cln_aln_file = job_output_folder + "09b_group_alignment_input.fasta.aln.cln"
                    tree_file = tree_constructor( group_aln_input_file, group_aln_file, group_cln_aln_file, repr_and_ath_fasta_file, "", result_folder, name, "09b_representatives_ath_detailed_", mode_aln, mode_tree, mafft, muscle, raxml, fasttree, cpur)
        
        # --- 09 construct a tree with landmark sequences --- #            
        land_aln_input_file = job_output_folder + "10_landmark_alignment_input.fasta"
        land_aln_file = job_output_folder + "10_landmark_alignment_input.fasta.aln"
        land_cln_aln_file = job_output_folder + "10_landmark_alignment_input.fasta.aln.cln"
        tree_constructor( land_aln_input_file, land_aln_file, land_cln_aln_file, landmark_file, clean_members_file, result_folder, name, "10_landmark_", mode_aln, mode_tree, mafft, muscle, raxml, fasttree, cpur)


        # --- documentation of execution time --- #
        end_time = datetime.datetime.now()
        duration = end_time - start_time
        seconds = duration.total_seconds()
        with open(doc_file, "a") as out:
            out.write("\n" + "Execution time: " + f"{int(seconds // 3600)}h:{int((seconds % 3600) // 60)}min:{int(seconds % 60)}sec")               

        sys.stdout.write("\n" + "Successfully finished execution. " + "\n")                
  
  
    
if '--out' in sys.argv and '--subject' in sys.argv and "--info" in sys.argv:
    main( sys.argv )
elif '--out' in sys.argv and '--subjectdir' in sys.argv and "--info" in sys.argv:
    main( sys.argv )
elif '--out' in sys.argv and '--subject' in sys.argv and "--baits" in sys.argv and "--baitsinfo" in sys.argv:
    main( sys.argv )
elif '--out' in sys.argv and '--subjectdir' in sys.argv and "--baits" in sys.argv and "--baitsinfo" in sys.argv:
    main( sys.argv )
else:
    sys.exit( __usage__ ) 