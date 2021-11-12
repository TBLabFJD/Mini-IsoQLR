# Mini-IsoQLR: Pipeline for Isoform Quantification using Long-Reads sequencing data
This pipeline was developed to detect and quantify isoforms from the expression of minigenes, whose cDNA was sequenced using Oxford Nanopore Technologies (ONT).
This protocol uses [GMAP](https://academic.oup.com/bioinformatics/article/21/9/1859/409207) aligner, which aligns cDNA sequences to a genome, using the parameter `--format=2` which generates a GFF3 file which contains the coordinates of the exons from all reads. Using this information, `Mini-IsoQLR.R` classify the mapped reads into isoforms. To do so, it:
1.  defines the consensus breakpoints (start and end of exons which are present in more than 5% of the reads by default), 
2.  sets a lower and higher threshold that is used to rescue reads whose breakpoints do not exactly coincide with the consensus breakpoints but are close to one,
3.  assigns the consensus breakpoints to the reads to establish the exons,
4.  classifies and quantifies reads into isoforms by concatenating the established exons.

[![Workflow](https://github.com/TBLabFJD/Mini-IsoQLR/blob/master/Workflow.png)](https://github.com/TBLabFJD/Mini-IsoQLR)
Complete bioinformatic protocol we followed using data sequended with a MinION sequencer(Oxford Nanopore Technology). Basecalling, demultiplexing and quality filtering is not included in the pipeline. This pipeline starts from the FASTQ file.

## License
Mini-IsoQLR source code is provided under the [**Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)**](https://creativecommons.org/licenses/by-nc-sa/4.0/). Mini-IsoQLR includes several third party packages provided under other open source licenses, please check them for additional details.

## Requirements
This pipeline was tested using the following program/library versions:
- GMAP (Genomic Mapping and Alignment Program) Package (version 2021-05-27)
- R (version 3.6.3)
- R libraries:
  - ggplot2 (version 3.3.5)
  - plyr(version 1.8.6)
  - optparse (version 1.7.1)


## Getting Started
The first step is to build the reference used by GMAP. To do so run the following command:
```sh
output_dir="/path/to/build/reference/genome/"
genomePrefix="genomePrefix"
fasta="/path/to/genome.fasta"

gmap_build -D ${output_dir} -d ${genomePrefix} ${fasta}
```

The second step is to perform the alignment using GMAP as follow:
```sh
treads=8 # number of threads
reference="/path/to/reference/genome/with/prefix"
genomePrefix="genomePrefix"
fastq="/path/to/sample.fastq"
gff3="/path/to/output/file.gff3"
error="/path/to/error/file.err" # it is used later on to know the number and ID of unmapped reads

gmap -n1 -t ${treads} --cross-species --gff3-add-separators=0 -f 2 -z auto -D ${reference} -d ${genomePrefix} ${fastq} > ${gff3} 2> ${error}
```
Note that the output of the alignment is printed through the stdout and the logfile through the stderr

The third and last step is to run the R script `Mini-IsoQLR.R`.
```sh
gff3="/path/to/input/file.gff3"
outputDir="/path/to/output/dir/"
error="/path/to/error/file.err"
runName="RunName" # This name will appear in the output file names and figures

mkdir ${outputDir}
Rscript Mini-IsoQLR.R -i ${gff3} -o ${outputDir} -l ${error} -r ${runName}
```
#### Arguments
| Argument | Default value | Summary |
| ------ | ------ | ------ |
| Required Arguments |
| -i/--input | null | input gff3 file (output from GMAP) |
| -o/--outputdir | null | output directory |
| -l/--logfile | null | error file from GMAP |
| -r/--runame | null | run name that will be used as a file prefix and in figure titles |
| Optional Tool Arguments |
| -b/--beginning | 0 | beginning position of the segment of study (trimming) |
| -f/--final | null | final position of the segment of study (trimming) |
| -t/--threshold | 5 | threshold (0-100%) used to filter the breakpoints present in mode than x % of the reads |
| -p/--padding | 5 | number of bases to each side from a defined break point to consider a read as part of that group |
| -a/--abundance | 5 | only isoforms with a percentage equal or higher will be displayed on the combined plot |


### Output
**Figures:**
- `RunName.all_isoform_information(.jpeg|.pdf)` - Figure showing all detected isoforms, their proportion and exon coordinates
- `RunName.breakpoints(.jpeg|.pdf)` - Figure representing the distribution of the breakpoints of all reads (it represents the coordinates that are in between the specify coordinates)
- `RunName.breakpoint_superplot(.jpeg|.pdf)` - Figure showing the percentage of reads with a breakpoint in each position of the whole genome (it does not take into account the specify coordinates, but all). Similar to the `RunName.breakpoints.jpeg` but it takes into account all positions and the percentage is calculated for each position (like an histogram where the number of bins = number of positions in the x-axis)
- `RunName.combined_plot_vertical(.jpeg|.pdf)` - Two aligned plots (by the nucleotide position) showing A) isoforms detected in more than 5% (by default) of the reads, their proportion and exon coordinates, and B) the distribution of the breakpoints of all reads (it represents the coordinates that are in between the specify coordinates)

**Tables:**
- `RunName.breakpoints_info.tsv` - Table containing the following information of all breakpoints: 
  - **breakpoint:** breakpoint position (coordinate), 
  - **type:** breakpoint type (start/end), 
  - **number_of_reads:** number of reads with that exact breakpoint, 
  - **percent_of_reads:** percentage of reads with that exact break point,
  - **left_point:** lower end of the threshold that is used to rescue reads whose breakpoints do not exactly coincide with the consensus breakpoints,
  - **right_point:** higher end of the threshold that is used to rescue reads whose breakpoints do not exactly coincide with the consensus breakpoints,
  - **final_number_of_reads:** number of reads with that exact breakpoint plus the rescue ones,
  - **increase_percentage:** percentage of increase when adding rescued reads in each breakpoint
- `RunName.isoform_freq.tsv` - Table containing the following information of all isoforms:
  - **Group:** isoform o category type (Unmapped, No_consensus_break_point_reads, Only_vector_reads)
  - **N_of_reads:** number of reads for each category or isoform 
  - **Prerc_partial:** percentage of abundance of each isoform (taking into account only classified reads into isoforms)
  - **Perc_total:** percentage of abundance of each isoform and category (taking into account all reads)
- `RunName.read_clasification.tsv` - Table containing the classification of all reads:
  - **read_ids:** read identifier,
  - **group:** group which the read is classified in (Unmapped, No_consensus_break_point_reads, Only_vector_reads, isoform_id)




## Run the example



```sh
repo_path="/local/path/to/the/Mini-IsoQLR/repository/"

example_dir="${repo_path}/example"


# Building the reference genome

output_dir="${example_dir}/references/gmap_index_pPSL3/"
genomePrefix="pPSL3"
fasta="${example_dir}/references/pSPL3_PAX6_Ex5-7.fasta"

gmap_build -D ${output_dir} -d ${genomePrefix} ${fasta}



# Mapping

treads=8 # number of threads
reference="${example_dir}/references/gmap_index_pPSL3/pPSL3"
genomePrefix="pPSL3"
fastq="${example_dir}/fastq/fastq_BARCODE01.fastq"
gff3="${example_dir}/mapped/cluster_cons_BARCODE01.gff3"
error="${example_dir}/mapped/log_BARCODE01.err" 

gmap -n1 -t ${treads} --cross-species --gff3-add-separators=0 -f 2 -z auto -D ${reference} -d ${genomePrefix} ${fastq} > ${gff3} 2> ${error}



# Running Mini-IsoQLR.R

outputDir="${example_dir}/results"
runName="Multiplex1_BARCODE01" # This name will appear in the output file names and figures

mkdir ${outputDir}
Rscript ${repo_path}/Mini-IsoQLR.R -i ${gff3} -o ${outputDir} -l ${error} -r ${runName}
```



### The figures should look like:

**`Multiplex1_BARCODE01.combined_plot_vertical.jpeg`**
[![Workflow](https://github.com/TBLabFJD/Mini-IsoQLR/blob/master/example/results/Multiplex1_BARCODE01.combined_plot_vertical.jpeg)](https://github.com/TBLabFJD/Mini-IsoQLR)


**`Multiplex1_BARCODE01.all_isoform_information.jpeg`**
[![Workflow](https://github.com/TBLabFJD/Mini-IsoQLR/blob/master/example/results/Multiplex1_BARCODE01.all_isoform_information.jpeg)](https://github.com/TBLabFJD/Mini-IsoQLR)


**`Multiplex1_BARCODE01.breakpoints.jpeg`** 
[![Workflow](https://github.com/TBLabFJD/Mini-IsoQLR/blob/master/example/results/Multiplex1_BARCODE01.breakpoints.jpeg)](https://github.com/TBLabFJD/Mini-IsoQLR)


**`Multiplex1_BARCODE01.breakpoint_superplot.jpeg`**
[![Workflow](https://github.com/TBLabFJD/Mini-IsoQLR/blob/master/example/results/Multiplex1_BARCODE01.breakpoint_superplot.jpeg)](https://github.com/TBLabFJD/Mini-IsoQLR)




