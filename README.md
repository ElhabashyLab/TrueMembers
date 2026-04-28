# TrueMembers
Automated pipeline for the identification and curation of homologous proteins within a protein family. The workflow combines HMM-based searches with iterative sequence ans structural refinement and filtering steps to improve sensitivity while maintaining high specificity.



<p align="center"> <img src="TrueMembers.png" alt="Figure description" width=\linewidth> </p> <p align="center"><em>Figure 1. (A) Flowchart illustrating key steps of the TrueMembers pipeline" </em></p>


# External software tools
After installation and before running compute jobs, the paths to the respective binaries of the following external tools have to be set in your EVcouplings job configuration file(s).

- jackhmmer (required)
Download and install HMMER from the HMMER webpage to a directory of your choice.

- USalign (optional)
evcouplings uses the hhfilter tool to filter sequence alignments. Installation is only required if you need this functionality.


# Databases
Automatic database setup
The TrueMember minimally needs a sequence database for alignment generation, and structure mapping information for comparison of evolutionary couplings to 3D structures.

Sequence and structure mapping databases for EVcouplings can be automatically downloaded using the included command line tool evcouplings_dbupdate. This tool will fetch the UniProt (SwissProt/TrEMBL), UniRef100 and UniRef90 databases, and generate SIFTS-based structure mapping tables.



##Authors
Aron Wishtner
Hadeer Elhabashy

Contact
If you have any questions or inquiries, please feel free to contact Hadeer Elhabashy at (Elhabashylab [@] gmail.com)

License
TrueMember code in this repository is licensed under the MIT License.
⭐ If this tool helped your research, please consider starring the repository.
