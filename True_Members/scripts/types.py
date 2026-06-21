from dataclasses import dataclass


@dataclass
class JobStats:
    annotation_files_jackhmmer: int = 0
    ids_jackhmmer: int = 0
    ids_foldseek_webserver: int = 0
    ids_combined: int = 0
    structure_downloads: int = 0
    uniprot_records: int = 0
    filtered_structures_foldseek: int = 0
    filtered_structures_usalign: int = 0
    filtered_structures_total: int = 0
    hmm_search_results: int = 0
    all_members: int = 0
    true_members: int = 0
    