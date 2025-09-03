These are the scripts and data for our preprint on Argentinian beans. 
In accordance with both my values for transparency and the requirements of the journal, I'm submitting all my data and scripts for analysis onto this repository. 

No project is perfect, and thus I treat this as a stepping stone and learning experience for better science down the road. Feel free to reach out with comments and corrections if you feel so inclined. 

This paper is the product of several people; namely, W George and A Chaparro. I took on this project after their departure and drove it to completion with final assays and analyses. 
The research and data generated are the work of all of us together, including as all AK Schmidt. Original pilot data was obtained by AD Steinbrenner a number of years ago. 
I credit W George to teaching me how to code in R and getting me started on my coding journey when I needed it most. 

As such, the data and code presented here are not solely mine, but our collaborative effort to identify what is going on with these puzzling beans from Jujuy. 
The base code from Figure 1A and B is by ADS, but I optimized it and annotated it. 
Much of the code from Figure 2B and C is by WG, and I cleaned it up and annotated it.

The rest of the code from Figure 3 as well as from the supplementals is my work with debugging and package recommendations from various GPTs. 

The story behind this paper is told through the three main figures and can be summarized as such: 
INR is a protein that mediates legume defenses against herbivory. 
The incidental discovery of low expression INR beans from Argentina is due to a shared genetic basis. (Fig 1)
Key phenotypic assays of ethylene, expression, and marker gene expression suggested they would also be susceptible to herbivores. (Figs 2 & 3)
Contrary to our hypothesis, herbivores fared comparatively well on these low expression plants to WT. (Fig 3)
Thus, low INR expression does not translate to herbivore susceptibility. 

Additionally, one key philosophical lesson in all of this is that we get what we screen for. Low expression is just that - low expression. 
But it is perfectly valid to predict that the presence of a "watchman" molecule at all in the system will still activate downstream defenses. Maybe not as much, but it's not zero. 

For a good comparison, check out our tritrophic interactions paper. 

Each script that I provide has its own internal annotation and commentary on why we did what we did, who generated which dataset, as well as any other choices we made in presenting the data. 

The code for Figure 1A loads data from an ethylene screen by ADS, wrangles the data, and performs a T-test to categorize the accessions as responsive or unresponsive. The plot displays the two groups in two facets, and shows the raw data of wound + water treatment vs. wound + In11 on the beans. Colouration is done in Illustrator. 

The code for Figure 1B loads world data from rnaturalearth and coalesces our 21 accessions by their USDA GRIN geotagged longitude and latitude values into a data frame. To plot, I grouped them into a summary data frame and displayed them as bubbles on a map that shows context of where in Argentina they are from. I left extra space to annotate the accessions in Illustrator. 

The code for Figure 2A loads data collected by AKS, wrangles the data by decoding the initial raw data values into the assignments we ended up using in our project, and exports the summary data frame. It then tests the data for normality before performing a t-test and plotting the two groups. The notable thing about this data is that homozygous and heterozygous alleles are grouped together and shown for ease as one group against the homozygous low-expression allele as a ratio of In11/H2O. 

The code for Figure 2B binds two data sets curated by WGS, decodes the internal alias values we used for the accessions, and tests them for normality. It then selects the four promoters of interest that we selected based on earlier expression data and runs a Welch's ANOVA on their data before plotting. 

The code for Figure 2C binds six data sets curated by WGS, decodes the internal alias values we used for the accessions, and filters out data that I was not able to interpret or wasn't clear. Because of the six different days of testing, and the wildly different values the plate reader gives out, we opted to use Z scores as a normalization factor across all the data over the weeks this data was compiled. It then tests for normality and runs a Welch's ANOVA before plotting the values in the same order as 2B. Lastly, there are some stats I ran in order to make specific summary claims about the data. 

The two code sets for Figure 3A load separate data frames curated by WGS and BQB, decodes the data, tests for normality, and runs either an ANOVA or t-test on the data before plotting. 

The code for Figure 3B loads a data frame with parental and NIL ethylene peak data curated by BQB. It wrangles it before splitting it into four data frames to test for normality before running an ANOVA on the data (since we have 3 treatments per biological replicate) and plotting the two NILs together as two facets. Colouration was performed in Illustrator. 

The code for Figure 3C loads a data frame with MYB expression data from qPCR curated by BQB. Two biological replicaets from one accession died, but they were to help establish a baseline for the two NILs, so we omitted them from testing and analysis. The rest of the replicates were tested for normality and were submitted to an ANOVA before plotting. Colouration was performed in Illustrator. 

The code for Figure 3D loads a combined data frame of herbivory initial and finall masses curated by BQB and filters out the intial masses so we can compare the final masses of the caterpillars. It tests the data for normality before running a Welch's ANOVA and plotting. 

The bash script for creating the BAM files takes the paired-end read files and sorts them against the Phaseolus vulgaris G19833 (Pvulgaris_442_v2.0.fa; Phytozome) genome. 
The bash script for quality control metrics analyzes all of the BAM files and ouputs a summary table of statistics for the data generated from Whole Genome Sequencing (WGS).
The bash script for SNP calls exises the INR locus and calls SNPs in each of the 21 accessions against G19833 and loads them into a multi track VCF and multi fasta for viewing on third party viewers or IGV/IGB. 

Code for supplementals forthcoming. DOI forthcoming. 
