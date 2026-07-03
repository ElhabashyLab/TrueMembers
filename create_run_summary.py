"""
the goal if this script is create a summary of a true members run.
It will summarize the results of all Jobs in a summary dict and save it to a json file which path is provided as an argument.

The input of this script will be a path in the runs dir to a pipeline run and the output path were to save the summary json file.

the jason should look something like this for each Job in the pipeline run:

{
    "Q5SK37": {
        iteration: 1, (can be found in: run_dir/final_true_members_results/Used_Seeds_i_*)
        new_true_members: 100, (can be found in: run_dir/final_true_members_results/Used_Seeds_i_*)
        annotation_files_jackhmmer: 4, (can be found in /runs/job_{seed_id}_iteration_{iteration_id}_2026-06-27_13-33-13/results/number_tracker.csv)
        ids_jackhmmer: 89882, (can be found in /runs/job_{seed_id}_iteration_{iteration_id}_2026-06-27_13-33-13/results/number_tracker.csv)
        ids_foldseek_webserver: 0, (can be found in /runs/job_{seed_id}_iteration_{iteration_id}_2026-06-27_13-33-13/results/number_tracker.csv)
        ids_combined: 89882, (can be found in /runs/job_{seed_id}_iteration_{iteration_id}_2026-06-27_13-33-13/results/number_tracker.csv)
        structure_downloads: 32325, (can be found in /runs/job_{seed_id}_iteration_{iteration_id}_2026-06-27_13-33-13/results/number_tracker.csv)
        uniprot_records: 32325, (can be found in /runs/job_{seed_id}_iteration_{iteration_id}_2026-06-27_13-33-13/results/number_tracker.csv)
        filtered_structures_foldseek: 0, (can be found in /runs/job_{seed_id}_iteration_{iteration_id}_2026-06-27_13-33-13/results/number_tracker.csv)
        filtered_structures_usalign: 5574, (can be found in /runs/job_{seed_id}_iteration_{iteration_id}_2026-06-27_13-33-13/results/number_tracker.csv)
        filtered_structures_total: 5573, (can be found in /runs/job_{seed_id}_iteration_{iteration_id}_2026-06-27_13-33-13/results/number_tracker.csv)
        hmm_search_results: 5531, (can be found in /runs/job_{seed_id}_iteration_{iteration_id}_2026-06-27_13-33-13/results/number_tracker.csv)
        all_members: 5573, (can be found in /runs/job_{seed_id}_iteration_{iteration_id}_2026-06-27_13-33-13/results/number_tracker.csv)
        true_members: 201 (can be found in /runs/job_{seed_id}_iteration_{iteration_id}_2026-06-27_13-33-13/results/number_tracker.csv)
    },
    ...
}

I put an example file for Used_Seeds_i_* and number tracker in the True_Members/test_data/test_results dir.
"""