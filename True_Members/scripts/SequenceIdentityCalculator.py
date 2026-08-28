class SequenceIdentityCalculator:
    
    def __init__(self):
        pass

    def compute_sequence_identity(self,
                                  output_tsv_path: str,
                                  current_seed_fasta_path: str,
                                  new_members_fasta_path: str,
                                  members_fasta_path: str,
                                  previous_used_seeds_fasta_path: str ):
        
        # create databses for new seed
        # mmseqs createdb current_seed_fasta_path new_seedDB
        
        # create db for new members
        # mmseqs createdb new_members_fasta_path new_membersDB
        
        # create db for all members
        # mmseqs createdb members_fasta_path membersDB
        
        # create db for previous used seeds
        # mmseqs createdb previous_used_seeds_fasta_path old_seedsDB
        
        # computation 1: search new seed against db of all unique sequences to get sequence identity for all to new seed
        #mmseqs search new_seedDB membersDB res_new_seed_to_all tmp
        #mmseqs convertalis new_seedDB membersDB res_new_seed_to_all res_new_seed_to_all.tsv \
        #--format-output "query,target,fseqid,alnlen"
        
        # computation 2: search all old seeds against all new sequences
        #mmseqs search old_seedsDB new_membersDB res_old_seeds_to_new_members tmp
        #mmseqs convertalis old_seedsDB new_membersDB res_old_seeds_to_new_members res_old_seeds_to_new_members.tsv \
        #    --format-output "query,target,fseqid,alnlen"
    
        # add output of two tsv files into master tsv
        pass
    
    def _create_db(self, fasta_path: str, db_name: str):
        # create a mmseqs2 database from a fasta file
        pass
    
    def _search_db(self, query_db: str, target_db: str, result_db: str):
        # search query_db against target_db and store results in result_db
        pass
    
    def _convert_results_to_tsv(self, query_db: str, target_db: str, result_db: str, output_tsv_path: str):
        # convert mmseqs2 results to a tsv file with sequence identity information
        pass
    
    def _add_results_to_master_tsv(self, new_results_tsv_path: str, master_tsv_path: str):
        # add new results to the master tsv file
        pass