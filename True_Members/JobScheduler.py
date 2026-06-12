import os
import pandas as pd


class JobScheduler:
    def __init__(self, config, current_iteration, seed_list_path, used_seeds_path, all_members_path):
        self.config = config      
        self.ITERATIONS = config["global"]["iterations"]
        self.CONVERGENCE_RUN = config["global"]["convergence_run"]
        self.CONVERGENCE_THRESHOLD = config["global"]["convergence_threshold"]
        
        self.seed_list_path = seed_list_path
        self.used_seeds_path = used_seeds_path
        self.all_members_path = all_members_path
        os.makedirs("true_members_results", exist_ok=True)
        
        self.current_iteration = current_iteration
    
    def start(self):
        while True:
            if not self.current_iteration is 1:
                seed_list = pd.read_csv(f"{self.seed_list_path}")
                if seed_list.empty:
                    print("Seed list is empty. Stopping pipeline.")
                    self.end_pipeline()
                    break
                if not self.CONVERGENCE_RUN and self.current_iteration > self.ITERATIONS:
                    print(f"Reached maximum number of iterations ({self.ITERATIONS}). Stopping pipeline.")
                    self.end_pipeline()
                    break
                else:
                    if self.CONVERGENCE_RUN:
                        converged = self.compute_convergence()
                        if converged:
                            print(f"Number of True Members converged with threshold {self.CONVERGENCE_THRESHOLD}. Stopping pipeline.")
                            self.end_pipeline()
                            break
                        else:
                            print(f"Number of True Members did not converge with threshold {self.CONVERGENCE_THRESHOLD}. Continuing pipeline.")
                            next_iteration = self.current_iteration + 1
                            print(f"Starting iteration {next_iteration} until convergence.")
                            self.start_job()
                    else:
                        next_iteration = self.current_iteration + 1
                        print(f"Continuing to next iteration ({next_iteration}).")
                        self.start_job()
            else:
                print(f"Starting iteration {self.current_iteration} of {self.ITERATIONS}")
                self.start_job()
                
                
                
    def end_pipeline(self):
        # copy seeds and all members to results directory
        pass


    def compute_convergence(self):
        return False


    def start_job(self):
        # Create Slurm script
        # Create true_members_config
        # Create root_dir
        # create Job 
        # start Job
        pass