### most functions were copied from MYB_annotator: http://dx.doi.org/10.1186/s12864-022-08452-5
### some functions were copied from KIPEs: https://doi.org/10.3390/plants9091103 and MaMYB: https://doi.org/10.1371/journal.pone.0239275 ###

### WARNING: do not use underscores in the bait IDs ###

#TODO: complete
### definition of familys in data/family_annotator.csv ###
###     Family	Baits	Information	Reference	HMM	Ath	RegEx ###
###     bHLH	bHLH_baits.fasta	bHLH_info.txt	-	-	-   - ###

import os, glob, sys, re, subprocess, dendropy
from operator import itemgetter
import pandas as pd
try:
	import hashlib
except ImportError:
	pass

__version__ = "v0.000"

#TODO: usage
__usage__ = """
                    python3 family_annotator.py
                    --baits <MYB_SEQ_FILE>
                    --info <MYB_CLASSIFICATION_FILE>
                    --out <OUTPUT_DIR>
                    --subject <SUBJECT_FILE (peptide,transcript,genomic sequences)> | --subjectdir <SUBJECT_FOLDER_WITH_SEQ_FILES>
                    
                    optional:
                    --search <INITIAL_SEARCH_TOOL>(blast|hmmer)[blast]
                    --mode <TREE_BUILDER>(fasttree|raxml)[fasttree]
                    --hmm <MYB_BAIT_HMM_FILE>
                    --refmybs <REF_MYB_FILE>
                    --ath <ATH_MYB_FILE_FOR_FINAL_TREE>
                    --name <STRING_USED_AS_PREFIX_IN_FILENAMES>
                    --collapse <REDUCES IN-PARALOGS_TO_ONE_REPRESENTATIVE>
                    --motifs <MOTIFS_TO_CHECK_FILE>
                    --cpu <NUMBER_OF_THREADS>[4]
                    --cpub <CPUs_TO_USE_FOR_BLASTp>[cpu]
                    --cpur <CPUs_TO_USE_FOR_RAxML>[cpu]
                    --cdsinput <CHANGES_EXPECTED_INPUT_TO_CDS>
                    --keepnames <PREVENTS_CUTTING_OF_NAMES_AT_FIRST_SPACE>
                    
                    --aln_mode <PATH_TO_aln_mode>[aln_mode]
                    --blastp <PATH_TO_AND_INCLUDING_BINARY>[blastp]
                    --hmmsearch <PATH_TO_HMMSEARCH>[hmmsearch]
                    --makeblastdb <PATH_TO_AND_INCLUDING_BINARY>[makeblastdb]
                    
                    --fasttree <PATH_TO_FASTTREE>[fasttree]
                    --raxml <PATH_TO_RAXML>[raxml]                
                    
                    --simcutp <BLASTP_SIMILARITY_CUTOFF>[0.8]
                    --poscutp <BLASTP_POSSIBLE_HIT_NUMBER_PER_BAIT_CUTOFF>[100]
                    --lencutp    <BLASTP_MIN_LENGTH_CUTOFF>[50]
                    
                    --numneighbours <NUMBER_OF_NEIGHBOURS_FOR_CLASSIFICATION> [10]
                    --neighbourdist <NEIGHBOUR_DISTANCE> [3]
                    --minneighbours <MINIMAL_NUMBER_OF_NEIGHBOURS> [0]
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
            fam_classification.update( { parts[1]: float( parts[2] ) } )
            line = f.readline()
    return fam_classification


def check_fam_IDs_across_files( fam_bait_seq_file, fam_info_file, ref_file, fam ):
    """! @brief check IDs across the different files """
    
    fam_status = True
    # --- check FASTA file for forbidden characters --- #
    forbidden_characters = [ ";", ":", "(", ")", "_" ]
    seqs = load_sequences( fam_bait_seq_file )
    header_string = "".join( seqs.keys() )
    for each in forbidden_characters:
        if each in header_string:
            sys.stderr.write( "Forbidden character detected in " + fam + " IDs (bait FASTA file): " + each + "(occurrences:" + str( header_string.count( each ) ) + ")\n" )
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
                    sys.stderr.write( "Issue in " + fam + " info file (unexpected status; only 'in' and 'out' permitted): " + line + "\n" )
                    sys.stderr.flush()
            line = f.readline()
    
    # --- compare IDs between bait and info file --- #
    missing_in_fasta = []
    for seq_ID in info_IDs.keys():
        try:
            seqs[ seq_ID ]
        except KeyError:
            missing_in_fasta.append( seq_ID )
    if len( missing_in_fasta ) > 0:
        sys.stderr.write( "Unmatched " + fam + " IDs (missing in bait FASTA file) :" + ";".join( missing_in_fasta ) + "\n" )
        sys.stderr.flush()
        fam_status = False
    
    missing_in_info = []
    for seq_ID in seqs.keys():
        try:
            info_IDs[ seq_ID ]
        except KeyError:
            missing_in_info.append(seq_ID )
    if len( missing_in_info ) > 0:
        sys.stderr.write( "Unmatched " + fam + " IDs (missing in info file) :" + ";".join( missing_in_info ) + "\n" )
        sys.stderr.flush()
        fam_status = False
    
    if len( ref_file ) > 0:    #only try to check if file actually exists
        missing_ref_ids = []
        with open( ref_file, "r" ) as f:
            line = f.readline()
            while line:
                x = line.strip()
                if "\t" in x:
                    x = line.split('\t')[0]
                try:
                    seqs[ x ]
                except KeyError:
                    missing_ref_ids.append( x )
                line = f.readline()
        if len( missing_ref_ids ) > 0:
            sys.stderr.write( "Reference " + fam + " IDs missing in FASTA and/or info file) :" + ";".join( missing_ref_ids ) + "\n" )
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
                if len( parts ) == 2:
                    refmem.update( { parts[0]: { 'id': parts[0], 'name': parts[1], 'function': "n/a", 'group': "" } } )
                elif len( parts ) == 3:
                    refmem.update( { parts[0]: { 'id': parts[0], 'name': parts[1], 'function': parts[2], 'group': "" } } )
                elif len( parts ) == 4:
                    refmem.update( { parts[0]: { 'id': parts[0], 'name': parts[1], 'function': parts[2], 'group': parts[3] } } )
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


def generate_documentation_file(doc_file,fam_bait_seq_file_all ,fam_bait_seq_file, fam_info_file, output_folder, raw_subject_file,
                                search, mode_aln, mode_tree, blastp, makeblastdb, hmmsearch ,cpub, cpur, mafft, muscle, raxml, fasttree, ref_file,
                                bitscore_cutoff_p, similarity_cutoff_p, possibility_cutoff_p, length_cutoff_p, cds_input, 
                                min_score_cutoff, neighbour_cutoff, mean_factor_cutoff, min_neighbour_cutoff, dist_cutoff_factorB, fam):
    """! @brief write documentation file with specified inputs and parameters """
    
    with open( doc_file, "w" ) as out:
        #TODO
        #out.write( "Please cite 'Pucker B (2021). Automatic identification and annotation of MYB gene family members in plants. doi:10.1101/2021.10.16.464636' when using MYB_annotator.py.\n\n" )
        out.write( "family_annotator.py version: " + __version__ + "\n" )
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
        out.write( "CPUs for BLASTp: " + str( cpub ) + "\n" )
        out.write( "CPUs for RAxML: " + str( cpur ) + "\n" )
        if len( ref_file ) > 0:
            ref_file_md5 = md5_calculator( ref_file )
            out.write( "Reference " + fam + " file: " + ref_file + "\t" + ref_file_md5 + "\n" )
        else:
            out.write( "Reference " + fam + " ffile: n/a\n" )
        if cds_input:
            out.write( "Type of input: CDS\n" )
        else:
            out.write( "Type of input: PEP\n" )
        
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
        out.write( "Distance cutoff for paralog clade masking: " + str( dist_cutoff_factorB ) + "\n" )
        
        
        # --- add tool versions --- #
        #TODO aln tool version?
        # try:
        #     aln_tool_version_raw = subprocess.Popen( args=aln_tool + " --version", stderr=subprocess.PIPE, shell=True )
        #     aln_tool_version = aln_tool_version_raw.stderr.read()
        #     out.write ( "aln_tool version: " + str( aln_tool_version )[2:-3] + "\n" )    #remove characters introduced through binary
        # except:
        #     out.write ( "aln_tool version detection failed.\n" )    #if no aln_tool installation was detected
        # out.write ( "FastTree version: PLEASE_ADD_MANUALLY\n"  )    #version not available via command
        try:
            raxml_version_raw = subprocess.Popen( args=raxml + " --version", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True )
            raxml_version = str( raxml_version_raw.stdout.read() ).strip()
            out.write ( "RAxML version: " + ( raxml_version[4:65]) + "...\n" )    #remove characters introduced through binary
        except:
            out.write ( "RAxML version detection failed.\n" )    #if no RAxML installation was detected
        try:
            hmmsearch_version_raw = subprocess.Popen( args=hmmsearch + " -h", stdout=subprocess.PIPE, shell=True )
            hmmsearch_version = str( hmmsearch_version_raw.stdout.read() ).strip().split('\n')[1]
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
        #TODO test methods         
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


def split_into_ingroup_and_outgroup( tree_file, in_list, out_list, neighbour_cutoff, mean_factor_cutoff, min_neighbour_cutoff ):
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
            for t2 in ref_gene_nodes:    #calculate distance to all other sequences in tree
                path_distance = pdm.path_edge_count( t1, t2)
                patr_distance = pdm.patristic_distance( t1, t2 )
                path_distances.append( { 'key': t2.label, 'val': path_distance } )
                patristic_distances.update( { t2.label: patr_distance } )
            in_counter = 0
            out_counter = 0
            sorted_distances = sorted( path_distances, key=itemgetter('val') )
            for each in sorted_distances[ : min( [ len( path_distances ), neighbour_cutoff ] ) ]:
                patr = patristic_distances[ each['key'] ]
                if patr < mean_factor_cutoff*my_mean_nearest_taxon_distance:    #exclude outliers on extremely long branches
                    if each['key'] in in_list:    #check if smalles path_distances are to in- or outgroup baits
                        in_counter += 1
                    else:
                        out_counter += 1
            if in_counter+out_counter > min_neighbour_cutoff:
                results.update( { t1.label: { 'score': float( in_counter ) / ( in_counter + out_counter ), 'in': in_counter, 'out': out_counter } } )
            else:
                results.update( { t1.label: { 'score': 0.0, 'in': in_counter, 'out': out_counter } } )
            #score ranges from 0 (non-family member) to 1 (family member)
    return results


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


def motif_check( seqs, motifs ):
    """! @brief screen sequences for motifs """
    
    results = {}
    for key in sorted( seqs.keys() ):
        results.update( { key: {} } )
        seq = seqs[ key ]
        for ID in motifs.keys():
            try:
                match = re.findall( motifs[ ID ], seq )[0]
                results[ key ].update( { key: match } )
            except:
                results[ key ].update( { key: "" } )
    return results


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

#%%


def main( arguments ):
    """! @brief run everything """
      
    output_folder = arguments[ arguments.index('--out')+1 ]
    
    fam_definition_file = "family_info.csv"
    family_definition = pd.read_csv(fam_definition_file,sep="\t",index_col=0)

#%%
    if '--family' in arguments:
        if ',' in arguments[ arguments.index('--family')+1 ]:
            family = arguments[ arguments.index('--family')+1 ].split(',')
        else:
            family = [arguments[ arguments.index('--family')+1 ]]
    else: 
        family = ['bHLH']
                   
    if '--subject' in arguments:
        raw_subject_files = [ arguments[ arguments.index('--subject')+1 ] ]
    else:
        subject_file_dir = arguments[ arguments.index('--subjectdir')+1 ]
        if not subject_file_dir[-1] == "/":
            subject_file_dir + "/"
        extensions = [ ".fasta", ".fa", ".fas", ".FASTA", ".FA", ".FAS", ".fna", ".FNA", ".cds", ".CDS", ".pep", ".PEP" ]
        raw_subject_files = [ ]
        for each in extensions:
            raw_subject_files += glob.glob( subject_file_dir + "*" + each )
        raw_subject_files = list( sorted( raw_subject_files ) )
   
    
    # --- Search options --- # 
    if '--search' in arguments:
        search = arguments[ arguments.index('--search')+1 ]
        if search not in [ "blast", "hmmer" ]:
            search = "blast"
    else:
        search = "blast"   
    if '--mode_aln' in arguments:
        mode_aln = arguments[ arguments.index('--mode_aln')+1 ]
        if mode_aln not in [ "mafft", "muscle" ]:
            mode_aln = "muscle"
    else:
        mode_aln = "muscle"    
    if '--mode_tree' in arguments:
        mode_tree = arguments[ arguments.index('--mode_tree')+1 ]
        if mode_tree not in [ "fasttree", "raxml" ]:
            mode_tree = "fasttree"
    else:
        mode_tree = "fasttree"
    
    
    # --- Execution Options --- #
    if "--name" in arguments:
        name = arguments[ arguments.index('--name')+1 ]
    else:
        name = ""
    if '--cdsinput' in arguments:
        cds_input = True
    else:
        cds_input = False   
    if "--keepnames" in arguments:
        trim_names = False
    else:
        trim_names = True
    if '--collapse' in arguments:
        collapse_mode = True
    else:
        collapse_mode = False


    # --- CPU usage --- #
    if '--cpu' in arguments:
        cpu = int( arguments[ arguments.index('--cpu')+1 ] )
    else:
        cpu = 4   
    if '--cpub' in arguments:
        cpub = int( arguments[ arguments.index('--cpub')+1 ] )
    else:
        cpub = cpu + 0
    if '--cpur' in arguments:
        cpur = int( arguments[ arguments.index('--cpur')+1 ] )
    else:
        cpur = cpu + 0
    

    # --- Tool paths --- #
    if '--blastp' in arguments:
        blastp = arguments[ arguments.index('--blastp')+1 ]
    else:
        blastp = "blastp"
    if '--makeblastdb' in arguments:
        makeblastdb = arguments[ arguments.index('--makeblastdb')+1 ]
    else:
        makeblastdb = "makeblastdb"
    if '--hmmsearch' in arguments:
        hmmsearch = arguments[ arguments.index('--hmmsearch')+1 ]
    else:
        hmmsearch = "hmmsearch"
    if "--mafft" in arguments:
        mafft = arguments[ arguments.index('--mafft')+1 ]
    else:
        mafft = "mafft"
    if "--muscle" in arguments:
        muscle = arguments[ arguments.index('--muscle')+1 ]
    else:
        muscle = "muscle"
    if "--fasttree" in arguments:
        fasttree = arguments[ arguments.index('--fasttree')+1 ]
    else:
        fasttree = "fasttree"
    if '--raxml' in arguments:
        raxml = arguments[ arguments.index('--raxml')+1 ]
    else:
        raxml = "raxml"
    
    
    # --- BLAST hit cutoffs --- #
    if "--simcutp" in arguments:
        similarity_cutoff_p = float( arguments[ arguments.index('--simcutp')+1 ] )
    else:
        similarity_cutoff_p=40.00
    if '--poscutp' in arguments:
        possibility_cutoff_p = int( arguments[ arguments.index('--poscutp')+1 ] )
    else:
        possibility_cutoff_p=100
    if '--lencutp' in arguments:
        length_cutoff_p = int( arguments[ arguments.index('--lencutp')+1 ] )
    else:
        length_cutoff_p=60      
    if '--bitcutp' in arguments:
        bitscore_p = int( arguments[ arguments.index('--bitcutp')+1 ] )
    else:
        bitscore_p=60
    

    # --- Candidate cutoffs --- #
    if '--minscore' in arguments:
        min_score_cutoff = int( arguments[ arguments.index('--minscore')+1 ] )
    else:
        min_score_cutoff=0.8    #minimal score for candidate to be considered as ingroup      
    if '--numneighbours' in arguments:
        neighbour_cutoff = int( arguments[ arguments.index('--numneighbours')+1 ] )
    else:
        neighbour_cutoff=30    #numbers of closest neighbour that is considered in ingroup/outgroup classification   
    if '--neighbourdist' in arguments:
        mean_factor_cutoff = float( arguments[ arguments.index('--neighbourdist')+1 ] )
    else:
        mean_factor_cutoff=10.0    #X*average nearest neighbor distance  
    if "--minneighbours" in arguments:
        min_neighbour_cutoff = int( arguments[ arguments.index('--minneighbours')+1 ] )
    else:
        min_neighbour_cutoff = 0    #minimal number of valid bait sequences (ingroup+outgroup) in range - 1    
    if '--paralogdist' in arguments:
        dist_cutoff_factorB = float( arguments[ arguments.index('--paralogdist')+1 ] )
    else:
        dist_cutoff_factorB=10.0    #X*average nearest neighbour distance used as cutoff in the monophyletic tip masking
    
    
    if output_folder[-1] != "/":
        output_folder += "/"
    if not os.path.exists( output_folder ):
        os.makedirs( output_folder )
 
    if not ('Baits' in family_definition.columns and 'BaitsInfo' in family_definition.columns):    
        sys.exit( "ERROR: Mandatory columns 'Baits' and/or 'BaitsInfo' are missing in information file: " + fam_definition_file )
#%% 
        
    # --- separated analyses for each family --- #
    for fam in family:
       
        # --- get inputs from definition file --- #        
        if not fam in family_definition.index: 
            sys.stdout.write( "ERROR: Family not defined in definition file:" + fam + "\n" )
            continue        
        
        fam_bait_seq_file_all = "data/" + family_definition.loc[fam,"Baits"]
        fam_info_file = "data/" + family_definition.loc[fam,"BaitsInfo"] 


        if 'ThinnedBaits' in family_definition.columns:
            fam_bait_seq_file = family_definition.loc[fam,'ThinnedBaits']
            fam_bait_seq_file = "data/" + fam_bait_seq_file  if len(fam_bait_seq_file ) > 1 else fam_bait_seq_file_all # to remove -          
        else:
            fam_bait_seq_file = fam_bait_seq_file_all
 
        if search == "hmmer" and 'HMM' in family_definition.columns:
            fam_bait_hmm = family_definition.loc[fam,'HMM']
            fam_bait_hmm  = "data/" + fam_bait_hmm  if len(fam_bait_hmm ) > 1 else "" # to remove -
        else:
            fam_bait_hmm = ''    
 
        if 'Reference' in family_definition.columns:
            ref_file = family_definition.loc[fam,'Reference']
            ref_file = "data/" + ref_file if len(ref_file) > 1 else "" # to remove -
        
        if 'Ath' in family_definition.columns:      
            fam_ath_file = family_definition.loc[fam,'Ath']
            fam_ath_file  = "data/" + fam_ath_file  if len(fam_ath_file ) > 1 else "" # to remove -
        else:
            fam_ath_file = ''
        
        if 'Motifs' in family_definition.columns:      
            motif_file = family_definition.loc[fam,'Motifs']
            motif_file  = "data/" + motif_file if len(motif_file ) > 1 else "" # to remove -
        else:
            motif_file = ''
        
    
        # --- separated analyses for each subject --- #           
        fam_ID = '' if len(family)==1 else family   
        for jidx, raw_subject_file in enumerate( raw_subject_files ):    #use jidx to generate unique IDs for all jobs
            
            # --- prepare output folder for each job if there are multiple --- #   
            if len( raw_subject_files ) == 1:
                job_output_folder = output_folder + fam_ID 
            else:
                job_ID = raw_subject_file.split('/')[-1].split('.')[0]
                job_output_folder = output_folder + fam_ID +'/' + str( jidx ).zfill(5) + "_" + job_ID + "/"
            
            if not job_output_folder[-1] +'/':
                job_output_folder = job_output_folder + '/' 
            
            if not os.path.exists( job_output_folder ):
                os.makedirs( job_output_folder )
              
                
            # --- 00 validation of inputs --- #
            #check if all baits are listed in the info file
            subject_file = job_output_folder + "clean_subject_sequences.pep.fasta"
            mapping_table_file = job_output_folder + "raw_subject_to_clean_subject_mapping_table.txt"
            if not os.path.isfile( subject_file ):
                clean_input_FASTA_file( raw_subject_file, subject_file, mapping_table_file, cds_input, trim_names )    #remove illegal characters from subject sequence headers
            #TODO reduced? for ref file check if in at??
            fam_check_status = check_fam_IDs_across_files( fam_bait_seq_file_all, fam_info_file, ref_file, fam )
            if not fam_check_status:
                sys.exit( "ERROR: analysis is stopped due to inconstistency of " + fam +" IDs between files" )
            subject_name_mapping_table = load_subject_name_mapping_table( mapping_table_file )
                        
            result_folder = job_output_folder + "RESULTS/"
            if not os.path.exists( result_folder ):
                os.makedirs( result_folder )           
            doc_file = result_folder + name + "00_documentation.txt"           
            generate_documentation_file(doc_file, fam_bait_seq_file_all, fam_bait_seq_file, fam_info_file, job_output_folder, raw_subject_file,
                                        search, mode_aln ,mode_tree, blastp, makeblastdb,hmmsearch, cpub, cpur, mafft, muscle,raxml, fasttree, ref_file,
                                        bitscore_p, similarity_cutoff_p, possibility_cutoff_p, length_cutoff_p, cds_input,
                                        min_score_cutoff,neighbour_cutoff, mean_factor_cutoff, min_neighbour_cutoff, dist_cutoff_factorB, fam)        
                
            
            # --- 01 find initial candidates --- #
            seq_search_result_file = job_output_folder +"01_seq_search_results.txt"
            if not os.path.isfile( seq_search_result_file ):
                if search == "blast":
                    blast_db = job_output_folder + "01_blastdb"
                    p = subprocess.Popen( args= makeblastdb + " -in " + subject_file + " -out " + blast_db + " -dbtype prot", shell=True )
                    p.communicate()
                    
                    p = subprocess.Popen( args= blastp + " -query " + fam_bait_seq_file_all + " -db " + blast_db + " -out " + seq_search_result_file + " -outfmt 6 -evalue 0.0001 -num_threads " + str( cpub ), shell=True )
                    p.communicate()
                else:      
                    #TODO: test
                    if not len(fam_bait_hmm) > 1:
                        sys.stdout.write( "ERROR: Missing information in " + fam_definition_file + " for " + fam + " family:" + "\n" + "Hmmer used as search option without defined motif in HMM column" + "\n" )  
                        continue                       
                    elif not os.path.isfile( fam_bait_hmm ):
                        sys.stdout.write( "ERROR: No file at path of HMM motif defined in " + fam_definition_file + " for " + fam + " family:" + "\n" + "File " + fam_bait_hmm + " does not exist")  
                        continue
                        
                    p = subprocess.Popen( args= hmmsearch + " --tblout " + seq_search_result_file + " " + fam_bait_hmm + " " + subject_file + " > " + job_output_folder+"01_hmmsearch_waste.txt", shell=True )
                    p.communicate()
            
            if search == "blast":
                seq_search_results = load_BLAST_results( seq_search_result_file, similarity_cutoff_p, possibility_cutoff_p, length_cutoff_p ,bitscore_p)    #load valid BLASTp results
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
        
        
            # --- 02 construct phylogenetic tree --- #
            aln_input_file = job_output_folder + "02_alignment_input.fasta"
            aln_file = job_output_folder + "02_alignment_input.fasta.aln"
            cln_aln_file = job_output_folder + "02_alignment_input.fasta.aln.cln"
            tree_file = tree_constructor( aln_input_file, aln_file, cln_aln_file, fam_bait_seq_file, candidate_file, job_output_folder, "02_", "", mode_aln, mode_tree, mafft, muscle, raxml, fasttree, cpur )
                 
            
            # --- 02 analyze tree file --- #
            clean_members_file = result_folder + name + "02a_clean_%ss.pep.fasta" % (fam)
            tmp_result_table = result_folder + name + "02b_in_out_%s_analysis_results.txt" % (fam)
            if not os.path.isfile( tmp_result_table ):
                in_list, out_list = load_bait_fam_anno( fam_info_file, fam_bait_seq_file )
                sys.stdout.write( "Number of ingroup %s baits: " % (fam) + str( len( in_list ) ) + "\n" )
                sys.stdout.write( "Number of outgroup %s baits: " % (fam) + str( len( out_list ) ) + "\n" )
                sys.stdout.flush()
                fam_classification = split_into_ingroup_and_outgroup( tree_file, in_list, out_list, neighbour_cutoff, mean_factor_cutoff, min_neighbour_cutoff )
                #dictionary with subject IDs: values are in the range of 0 (non-family member) to 1 (family member)
                
                with open( clean_members_file, "w" ) as out:
                    with open( tmp_result_table, "w" ) as out2:
                        out2.write( "OriginalID\tCleanID\tScore\tIngroupMatches\tOutgroupMatches\n" )
                        candidate_order = list( sorted( fam_classification.keys() ) )
                        for candidate in candidate_order:
                            if fam_classification[ candidate ]['score'] > min_score_cutoff:
                                out.write( '>' + candidate + "\n" + subject_sequences[ candidate ] + "\n" )
                            out2.write( "\t".join( list( map( str, [subject_name_mapping_table[ candidate ],
                                                                    candidate,
                                                                    fam_classification[ candidate ]['score'],
                                                                    fam_classification[ candidate ]['in'],
                                                                    fam_classification[ candidate ]['out']
                                                                    ] ) ) ) + "\n" )
                if cds_input:    #generate corresponding CDS clean member output file
                    cds_clean_members_file = clean_members_file.replace( ".pep.fasta", ".cds.fasta" )
                    with open( cds_clean_members_file, "w" ) as out:
                        for candidate in candidate_order:
                            if fam_classification[ candidate ]['score'] > 0.5:
                                out.write( '>' + candidate + "\n" + cds_subject_sequences[ candidate ] + "\n" )
            else:
                fam_classification = load_fam_classification_from_file( tmp_result_table )
            clean_members = load_sequences( clean_members_file )
            if len( list( clean_members.keys() ) ) < 1:
                sys.exit( "ERROR: no " + fam + "s detected." )
                
                
            # --- 03 find closest reference --- #
            if len( ref_file ) > 0:    #only performed if reference file is provided
                group_around_ref_file = result_folder + name + "03a_group_around_ref_%ss.txt"  % (fam)   #produce table sorted by reference members
                new_2_ref_mapping_file = result_folder + name + "03b_new_2_ref_%s_mapping_file.txt"  % (fam)   #produce table sorted by subject sequences
                if not os.path.isfile( new_2_ref_mapping_file ):
                    ref_members = load_ref_members( ref_file )    #load IDs and trivial name from additional text file (dictionary)
                    new2ref_mapping_table, new_per_ref_members = member_group_assignment( ref_members, tree_file, clean_members.keys() )
            
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
            
            #TODO
            # --- 04 check for presence of MYB domains --- #
            # myb_domain_check_file = result_folder + name + "04a_MYB_domain_check.txt"        #produce table with MYB domain status and sequence
            # myb_domain_fasta_file = result_folder + name + "04c_MYB_domain_check.pep.fasta"        #FASTA file containing MYB domain sequence
            # myb_domain_doc_file = result_folder + name + "04d_MYB_domain_check.doc.txt"        #text file summarizing MYB domain detection results
            # if not os.path.isfile( myb_domain_check_file ):
            #     MYB_domain_check_wrapper( clean_mybs_file, myb_domain_check_file, myb_domain_fasta_file, myb_domain_doc_file, subject_name_mapping_table )
        
        
            #TODO
            # --- 04 check for different motifs --- #
            # motif_check_file_summary = result_folder + name + "04b_motif_check.txt"    #produce table with motif (0/1)
            # motif_check_file_seqs = result_folder + name + "04b_motif_check_details.txt"    #produce table with motif sequences
            # if len( motifs_file ) > 0:
            #     motifs = load_motifs_from_file( motifs_file )
            #     motif_check_results = motif_check( clean_members, motifs )    #prim key = seqID, secondary key = motifs
            #     motif_names = list( sorted( list( motifs.keys() ) ) )
            #     with open( motif_check_file_summary, "w" ) as out1:
            #         out1.write( "RefMember\t"+"\t".join(motif_names)+"\n" )
            #         with open( motif_check_file_seqs, "w" ) as out2:
            #             out2.write( "RefMember\t"+"\t".join(motif_names)+"\n" )
            #             candidates = list( sorted( clean_members.keys() ) )
            #             for candidate in candidates:
            #                 new_line_details = [ subject_name_mapping_table[ candidate ] ]
            #                 new_line_summary = [ subject_name_mapping_table[ candidate ] ]
            #                 for mot in motif_names:
            #                     new_line_details.append( motif_check_results[ candidate ][ mot ] )
            #                     if len( motif_check_results[ candidate ][ mot ] ) > 0:
            #                         new_line_summary.append( "1" )
            #                     else:
            #                         new_line_summary.append( "0" )
            #             out1.write("\t".join(new_line_summary) + "\n")
            #             out2.write("\t".join(new_line_details) + "\n")
                        
                        
            # --- 05 construct a final tree --- #
            fin_aln_input_file = job_output_folder + "05_fin_alignment_input.fasta"
            fin_aln_file = job_output_folder + "05_fin_alignment_input.fasta.aln"
            fin_cln_aln_file = job_output_folder + "05_fin_alignment_input.fasta.aln.cln"
            tree_file = tree_constructor( fin_aln_input_file, fin_aln_file, fin_cln_aln_file, fam_bait_seq_file, clean_members_file, result_folder, name, "05_", mode_aln, mode_tree, mafft, muscle, raxml, fasttree, cpur)
            
            
            # --- 06 construct a final tree with Ath members --- #
            if len( fam_ath_file) > 0:
                ath_fin_aln_input_file = job_output_folder + "06_ath_fin_alignment_input.fasta"
                ath_fin_aln_file = job_output_folder + "06_ath_fin_alignment_input.fasta.aln"
                ath_fin_cln_aln_file = job_output_folder + "06_ath_fin_alignment_input.fasta.aln.cln"
                tree_file = tree_constructor( ath_fin_aln_input_file, ath_fin_aln_file, ath_fin_cln_aln_file, fam_ath_file, clean_members_file, result_folder, name, "06_", mode_aln, mode_tree, mafft, muscle, raxml, fasttree, cpur)
                
                
            # --- 07 find in species-specific paralogs (in-paralogs) --- #
            if collapse_mode:
                paralog_groups = establish_paralog_groups( tree_file, clean_members.keys(), dist_cutoff_factorB )    #get list of sublist; each sublist represents one paralog group
                paralog_group_file = result_folder + name + "07a_repr_%ss.txt" % (fam)
                repr_clean_file = result_folder + name + "07b_repr_%ss.pep.fasta" % (fam)
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
                            
            
                # --- represent cluster only by longest sequence --- #
                if len( fam_ath_file ) > 0:
                    repr_ath_fin_aln_input_file = job_output_folder + "07c_repr_ath_fin_alignment_input.fasta"
                    repr_ath_fin_aln_file = job_output_folder + "07c_repr_ath_fin_alignment_input.fasta.aln"
                    repr_ath_fin_cln_aln_file = job_output_folder + "07c_repr_ath_fin_alignment_input.fasta.aln.cln"
                    tree_file = tree_constructor( repr_ath_fin_aln_input_file, repr_ath_fin_aln_file, repr_ath_fin_cln_aln_file, fam_ath_file, repr_clean_file,result_folder, name, "07c_", mode_aln, mode_tree, mafft, muscle, raxml, fasttree, cpur )
                    
                    # --- define groups for Ath members and include these group names in tip labels of a phylogenetic tree --- #
                    if len( ref_file ) > 0:    #only performed if reference member file is provided
                        ref_members = load_ref_members( ref_file ) # id : id name function group
                        new2ref_mapping_table = load_candidate_mem_to_mem_mapping_table( new_2_ref_mapping_file ) # cnd : ref
                        repr_and_ath_for_tree = load_sequences( repr_ath_fin_aln_input_file ) # id : seq of ath and cnds
                        repr_and_ath_fasta_file = result_folder + name + "08a_repr_ath_%s.fasta" % (fam)
                        construct_fasta_file_w_repr_and_aths( ref_members, new2ref_mapping_table, repr_and_ath_for_tree, repr_and_ath_fasta_file )
                        group_aln_input_file = job_output_folder + "08b_group_alignment_input.fasta"
                        group_aln_file = job_output_folder + "08b_group_alignment_input.fasta.aln"
                        group_cln_aln_file = job_output_folder + "08b_group_alignment_input.fasta.aln.cln"
                        tree_file = tree_constructor( group_aln_input_file, group_aln_file, group_cln_aln_file, repr_and_ath_fasta_file, "", result_folder, name, "08b_", mode_aln, mode_tree, mafft, muscle, raxml, fasttree, cpur)
                    # TODO
                    # --- check for presence of MYB domains --- #
                    #repr_domain_check_file = result_folder + name + "08d_%s_domain_check.txt" % (fam)   #produce table with family domain status and sequence
                    #repr_domain_fasta_file = result_folder + name + "08e_%s_domain_check.fasta" % (fam)   #FASTA file containing family sequence
                    #repr_domain_doc_file = result_folder + name + "08f_%s_domain_check.doc.txt" % (fam)   #text file summarizing family domain detection results
                    #if not os.path.isfile( repr_domain_check_file ):
                    #    MYB_domain_check_wrapper( repr_clean_myb_file, repr_myb_domain_check_file, repr_myb_domain_fasta_file, repr_myb_domain_doc_file, subject_name_mapping_table )
                        
        # TODO                
        # --- summarize stats of all species --- #
        #if len( raw_subject_files ) > 1:    #only useful if there are more than one species
        #    summary_file4 = output_folder + "4_domain_detection_summary.txt"
        #    summarize_domain_counts( summary_file4, raw_subject_files, "04d", output_folder, name )
        #    if collapse_mode:
        #        summary_file8 = output_folder + "8_domain_detection_summary.txt"
        #        summarize_domain_counts( summary_file4, raw_subject_files, "08f", output_folder, name )
                            
 
                            
if "--search" in sys.argv:    #initial search could be based on hmmsearch
    if sys.argv[ sys.argv.index( '--search' )+1 ] == "hmmer":
        if '--out' in sys.argv and '--subject' in sys.argv: #!!!
            main( sys.argv )
        elif '--hmm' in sys.argv and '--out' in sys.argv and '--subjectdir' in sys.argv:
            main( sys.argv )
        else:
            sys.exit( __usage__ )
    else:
        if '--out' in sys.argv and '--subject' in sys.argv:
            main( sys.argv )
        elif '--out' in sys.argv and '--subjectdir' in sys.argv:
            main( sys.argv )
        else:
            sys.exit( __usage__ )
    
else:    #initial search is based on BLAST (default)
    if '--out' in sys.argv and '--subject' in sys.argv:
        main( sys.argv )
    elif '--out' in sys.argv and '--subjectdir' in sys.argv:
        main( sys.argv )
    else:
        sys.exit( __usage__ ) 