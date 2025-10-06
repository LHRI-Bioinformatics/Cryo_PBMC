Instruction on working with the PBMC data

1. Clone the code to local working directory,
 
   $git clone https://github.niaid.nih.gov/LHRI-Bioinformatics/Cryo_PBMC.git
   

2. The GEO accession# for our data is GSE304704, the list of Run accessions can be found and downloaded there. fastq files can be download by using these RUN accessions,

   $mkdir data
   
   $cd data

   $fasterq-dump --split-files SRR34882987
   

3. Run the cellranger pipeline and save the filtered count data into a folder "data/filtered_count"; it is recommended to make a directory "analysis" for running the code and saving the analysis result,

   $mkdir data/filtered_count

   $mkdir analysis/


