# TrueMembers
Automated pipeline for the identification and curation of homologous proteins within a protein family. The workflow combines HMM-based searches with iterative sequence ans structural refinement and filtering steps to improve sensitivity while maintaining high specificity.



<p align="center"> <img src="TrueMembers.png" alt="Figure description" width=\linewidth> </p> <p align="center"><em>Figure 1. (A) Flowchart illustrating key steps of the TrueMembers pipeline" </em></p>


# External software tools
Please install the following software and add the paths to the respective binaries in the configuration file(s).

- jackhmmer (required)
Download and install HMMER from the HMMER webpage to a directory of your choice.
http://hmmer.org/download.html

- USalign (required)
https://github.com/pylelab/USalign

-Foldseek (required)
https://github.com/steineggerlab/foldseek


# Databases
The TrueMember tool requires the following minimal databases: UniProt (Swiss-Prot/TrEMBL), UniRef100, UniRef90, and the AlphaFold protein structure database. The following script automatically downloads and sets up all required databases, ensuring that the environment is correctly initialized for downstream analysis without manual configuration.



##Authors
Aron Wishtner
Hadeer Elhabashy

Contact
If you have any questions or inquiries, please feel free to contact Hadeer Elhabashy at (Elhabashylab [@] gmail.com)

License
TrueMember code in this repository is licensed under the MIT License.
⭐ If this tool helped your research, please consider starring the repository.
