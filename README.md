CLIP-seq CLUSTER DETECTION TOOL 
================

The software developed is available in this page. This software is used to get Clusters on experimental data of RNA-Binding Proteins (RBP), including profiles of clusters.

--
REQUIREMENTS

Java 1.5 or more

Minimum 512 Mb of RAM (Recommended 4 Gb of RAM)

Any processor (Recommended multi-core processor)

Enough space to save results
	

	
--
CONFIGURATION AND INPUTS

You must to create a set-up file (with extension .ini). In that file write values for parameters required to run the application.
Parameters are: 

DataSet=Data set (in SAM format) of reads obtained experimentally.

Filter_Minimum_Ocurrences=This is a filter to reduce no-significant clusters by have a little number of reads aligned. By default is 5.

Filter_Minimum_Read_length=This is a filter to reduce no-significant reads by have a little size of sequence. By default is 5.

Folder_Results=Path where the results will be saved.

E-mail=Mail to send a message when the experiment is finished.

FASTA_Folder=Path of FASTA files.

SNP_Folder=Path of SNP files.

Export_SAM=Flag to export cluster's list in SAM. By default is FALSE.

Export_BED=Flag to export cluster's list in BED. By default is FALSE.

Filter_Entrophy=Flag to use Shannon's Entropy to filter noise reads. By default is TRUE.

Mutations=Type of mutations to consider (Options:All, T2C, A2C, A2G, A2T, C2A, C2G, C2T, G2A, G2C, G2T, T2A, T2G).

Insertions=Type of insertions to consider (Options: All, A, C, G, T).

Deletions=Type of deletions to consider (Options: All, A, C, G, T).


Additional parameters in v1.5:

Background=Path of background file.

Background_Score=Filter of minimum read's score value in background to take into account in comparision.


If there's no FASTA files, the application show a message and not execute.
If folder results non exists, the application will create it.


The setup.ini file (in the same root of this file) is an example of a configuration setup.

--
EXECUTION

To run, just write on terminal:
java -jar ClusterDetection_v*.jar Setup_file_name.ini


--
OUTPUTS

1) Total_clusters: List of all clusters (Optionally can be in format .SAM and .BED)

2) Filtered_Clusters: List of clusters that pass the filters

3) Folders of details of clusters, organized by chromosome
