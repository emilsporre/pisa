PISA analysis pipeline

This file details how to run the PISA analysis pipeline on ubuntu, including maxquant and MSstats data anlysis in R. 
This tutorial uses MaxQuant 2.0.3.0 with dotnet 3.1.426 on ubuntu 22.04, but you should be able to run a different version of MaxQuant on a different OS and still use the MSstats analysis from this tutorial.

#### Contents ####
1. Preparation of .raw files
2. Preparation of mqpar.xml file
3. Run MaxQuant peptide search
4. Prepare annotation files
5. Run PISA analysis R script
6. Results

#### 1. Preparation of .raw files ####
When running a new MaxQuant analysis, move your .raw files to a new folder named "date_experiment_organism", e.g. "20230216_ascorbate_arabidopsis", in the /pisa/data folder.

If you are only interested in the R data analysis, just move the folder of your finished maxquant analysis to the same folder and proceed with step 4.

#### 2. Preparation of mqpar.xml file ####
To run a MaxQuant search you need a mqpar.xml file that specifies all the required settings. There is an example template: 
/pisa/files/mqpar_top15_tmt_linux_template.xml

When running a new MaxQuant analysis, this file is to be copied to the folder where the .raw files are and key elements of the file needs to be changed. By default, this file is for 6 fractions, with TMT labels and no normalization in MaxQuant. The elements that need be changed from run to run are the paths to the raw files, the path to the .fasta background file (organism specific) and the isobaric label adjustments (varies between batches of labels). The adjustment factors come as a sheet of paper with the package of isobaric labels.

At the top of the mqpar.xml file is an area that looks like so:

<fastaFiles>
      <FastaFileInfo>
         <fastaFilePath>/home/user/pisa/resources/arabidopsis/arabidopsis_20240624.fasta</fastaFilePath>
         <identifierParseRule>>(.*)</identifierParseRule>
         <descriptionParseRule>>(.*)</descriptionParseRule>
         <taxonomyParseRule></taxonomyParseRule>
         <variationParseRule></variationParseRule>
         <modificationParseRule></modificationParseRule>
         <taxonomyId></taxonomyId>
      </FastaFileInfo>
   </fastaFiles>
   
Here the <fastaFilePath> string needs to match the organism in question

Around line 150 there is an area that looks like so:

<filePaths>
      <string>/home/user/pisa/data/20230216_ascorbate_arabidopsis/HF_20230216_DDA_350to1500mz_120min_1to32_700nLmin_60C_4ul_TMT_frac_C1.raw</string>
      <string>/home/user/pisa/data/20230216_ascorbate_arabidopsis/HF_20230216_DDA_350to1500mz_120min_1to32_700nLmin_60C_4ul_TMT_frac_C2.raw</string>
      <string>/home/user/pisa/data/20230216_ascorbate_arabidopsis/HF_20230216_DDA_350to1500mz_120min_1to32_700nLmin_60C_4ul_TMT_frac_C3.raw</string>
      <string>/home/user/pisa/data/20230216_ascorbate_arabidopsis/HF_20230216_DDA_350to1500mz_120min_1to32_700nLmin_60C_4ul_TMT_frac_C4.raw</string>
      <string>/home/user/pisa/data/20230216_ascorbate_arabidopsis/HF_20230216_DDA_350to1500mz_120min_1to32_700nLmin_60C_4ul_TMT_frac_C5.raw</string>
      <string>/home/user/pisa/data/20230216_ascorbate_arabidopsis/HF_20230216_DDA_350to1500mz_120min_1to32_700nLmin_60C_4ul_TMT_frac_C6.raw</string>
   </filePaths>
   
These strings must be the paths to the relevant .raw files.

Around line 270 there is an area that looks like so, repeated for each isobaric label in the TMTPro 16-plex product:

<isobaricLabels>
         	<IsobaricLabelInfo>
               <internalLabel>TMTpro16plex-Lys126C</internalLabel>
               <terminalLabel>TMTpro16plex-Nter126C</terminalLabel>
               <correctionFactorM2>0</correctionFactorM2>
               <correctionFactorM1>0</correctionFactorM1>
               <correctionFactorP1>9.09</correctionFactorP1>
               <correctionFactorP2>0.32</correctionFactorP2>
               <tmtLike>True</tmtLike>
            </IsobaricLabelInfo>

The correction factors must be the same as on the paper sheet delivered with the labels for every label.


#### 3. Run MaxQuant peptide search ####

To run the analysis, simply run the following command in the terminal, where the path is substituted with the path to your mqpar.xml file:

maxquant path

e.g.:

maxquant /home/user/pisa/data/20230216_ascorbate_arabidopsis/mqpar_top15_tmt_linux.xml

Much of MaxQuant is single threaded and the search can take a long time even on powerful computers, especially if the search space is large (i.e. there are many proteins in the supplied .fasta file)

The terminal output should look like so:

Configuring 
Assemble run info 
Finish run info 
Testing fasta files 
Testing raw files 
Feature detection 
Deisotoping 
MS/MS preparation 
Calculating peak properties 
Combining apl files for first search 
Preparing searches 
MS/MS first search 
Read search results for recalibration 
Mass recalibration 
Calculating masses 
MS/MS preparation for main search 
Combining apl files for main search 
MS/MS main search 
Preparing combined folder  
Correcting errors 
Reading search engine results 
Preparing reverse hits 
Finish search engine results 
Filter identifications (MS/MS) 
Calculating PEP 
Copying identifications 
Applying FDR 
Assembling second peptide MS/MS 
Combining second peptide files 
Second peptide search 
Reading search engine results (SP) 
Finish search engine results (SP) 
Filtering identifications (SP) 
Applying FDR (SP) 
Re-quantification 
Reporter quantification 
Retention time alignment 
Matching between runs 1 
Matching between runs 2 
Matching between runs 3 
Matching between runs 4 
Prepare protein assembly 
Assembling proteins 
Assembling unidentified peptides 
Finish protein assembly 
Updating identifications 
Estimating complexity 
Prepare writing tables  
Writing tables 
Finish writing tables

#### 4. Prepare annotation files ####
Once the the MaxQuant analysis has finished, an annotation file is required before the R script that handles normalization, protein quantification and statistical comparison can be run. For MSStats, the annotation file is called annotation.mq3 with path /home/user/pisa/files/annotation.mq3 and looks like so:

Run									Fraction	TechRepMixture	Channel	Condition	Mixture	BioReplicate
HF_20230125_DDA_350to1500mz_120min_1to32_700nLmin_60C_4ul_TMT_frac_D1	1	1		channel1	0		Mixture1	0
HF_20230125_DDA_350to1500mz_120min_1to32_700nLmin_60C_4ul_TMT_frac_D1	1	1		channel2	1		Mixture1	1
HF_20230125_DDA_350to1500mz_120min_1to32_700nLmin_60C_4ul_TMT_frac_D1	1	1		channel3	0		Mixture1	0
...
HF_20230125_DDA_350to1500mz_120min_1to32_700nLmin_60C_4ul_TMT_frac_D6	6	1		channel14	3		Mixture1	3
HF_20230125_DDA_350to1500mz_120min_1to32_700nLmin_60C_4ul_TMT_frac_D6	6	1		channel15	2		Mixture1	2
HF_20230125_DDA_350to1500mz_120min_1to32_700nLmin_60C_4ul_TMT_frac_D6	6	1		channel16	3		Mixture1	3

The file presumes that there are 6 fractions and that TMTPro 16-plex is being used. There are different files depending on the number of concentrations used, annotation.mq1 for 1 concentration, annotation.mq2 for 2 concentrations etc.
OBSERVE! This pipeline has thus far only been run for 3 concentrations, errors may occur when changing that.
The most important part is that the correct conditions are coupled to the correct channel. The channels are the different isobaric labels in ascending order, i.e. channel1 is the 126C label, and channel 16 is the 134N label. The file presumes that condition 1 & 3 alternate the first 8 channels and that conditions 2 & 4 are alternating in channels 9-16. As long as this is the case, the file does not need to be updated. BioReplicate should be the same as condition. The R script automatically fills in the .raw file names in the Run column.

In short, unless the label - condition pattern has changed, this file need not be changed at all. If it has changed, it is the channel - condition pattern that need be adjusted.

#### 5. Run PISA analysis R script ####

The path to the script is /home/user/pisa/ms_analysis_pisa_msstats_20230130.R

To run it, simply execute the following command in the terminal:

Rscript /home/user/pisa/ms_analysis_pisa_msstats_20230130.R folder name nbr_of_concentrations q_value_cutoff outlier_channels

e.g.:

Rscript /home/user/pisa/ms_analysis_pisa_msstats_20230130.R 20230216_ascorbate_arabidopsis 3 0.05

would run the script for the 20230216_acorbate_arabidopsis data for 3 concentrations and with a qvalue cutoff of 0.05
If one wants to exclude certain replicates as outliers, simply add on the number of the channel(s) in question after the qvalue cutoff in the command with spaces in between, eg:

Rscript /home/user/pisa/ms_analysis_pisa_msstats_20230130.R 20230216_ascorbate_arabidopsis 3 0.05 3 7

would run the script identically to the above command, but would also exclude channels 3 and 7. Which channel is associated with which replicate is specified in the annotation file prepared above.

#### 6. Results ####

The results are automatically written into the /home/user/pisa/results folder with one table containing one row per protein detected and one .pdf file with various plots of interest.
