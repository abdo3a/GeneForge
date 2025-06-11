# GeneForge
A Nextflow Pipeline for Gene Prediction and Functional Annotation

**GeneForge** is a flexible Nextflow pipeline for gene prediction and optional functional annotation of eukaryotic genomes. It integrates BRAKER3 and funannotate for gene prediction, compares their quality using BUSCO, and supports functional annotation with tools like Phobius, InterProScan, and eggNOG-mapper. The pipeline offers three run modes: prediction only, annotation only, or both, with customizable parameters for species, input paths, and annotation preferences.

## Features

- **Dual Prediction**: Runs BRAKER3 and funannotate in parallel for robust gene predictions.
- **Quality Assessment**: Compares predictions using BUSCO completeness scores.
- **Functional Annotation**: Optional annotation with Phobius, InterProScan, eggNOGMapper, and funannotate.
- **Flexible Run Modes**:
   - *prediction_only*: Gene prediction and BUSCO comparison.
   - *annotation_only*: Functional annotation using existing predictions.
   - *both*: Prediction followed by annotation, using the best BUSCO-scored prediction or user-specified choice.
- **Containerized**: Uses Docker for reproducibility, with a custom image for most tools and BRAKER3’s official image.
- **Customizable**: Supports user-defined species, input paths, BUSCO lineages, and computational resources.

## Prerequisites

- **Nextflow**: Version ≥ 22.04 (install via wget -qO- https://get.nextflow.io | bash).
- **Docker**: Required for containerized execution (install via your package manager or Docker’s website).
- **System**: Linux/Unix-based system with sufficient RAM (≥16GB recommended) and CPU cores (≥8 recommended).
- **Disk Space**: Depends on genome size and RNA-Seq data; 100GB+ recommended.

## Installation

1. **Clone the Repository**:
````bash
git clone https://github.com/yourusername/GeneForge.git
cd GeneForge
````
2. **Build the Docker Image**: Save the provided ```Dockerfile``` and build the custom image:
````bash
docker build -t genome-annotation:latest .
````
The BRAKER3 image (```teambraker/braker3:latest```) is pulled automatically by Nextflow.

3. **Install Nextflow** (if not already installed):
````bash
curl -s https://get.nextflow.io | bash
chmod +x nextflow
mv nextflow /usr/local/bin/
````
## Input Files

Place input files in your working directory (e.g., ```/path/to/project```):

- **Required**:
  - Genome FASTA: Unmasked genome sequence (e.g., ```genome.fasta```).
  - Masked Genome FASTA: Repeat-masked genome (e.g., ```genome_masked.fasta```).
  - Protein Evidence FASTA: Homologous protein sequences (e.g., ```proteins.fa```).
  - Funannotate Database: Directory with funannotate databases (e.g., ```/opt/funannotate_db```).
- **Optional**:
  - **RNA-Seq Reads**: Paired-end FASTQ files (gzipped) in a directory (e.g., ```data/rnaseq/sample_{R1,R2}_001.fastq.gz```).
    If provided, requires:
    - **Adapter FASTA**: Sequencing adapters for trimming (e.g., ```adapters.fa```).
  - **For annotation_only mode**, ensure GFF and protein files exist in the output directory (```results/```).  

## Usage

Run the pipeline with ```nextflow run genome_annotation.nf```, specifying parameters and the configuration file (```custom.config```). Paths are relative to the Docker mount point (```--docker_mount_base```, default: ```/workspace```).

**Example Commands**
1. **Prediction Only (With RNA-Seq)**:
````bash
nextflow run genome_annotation.nf \
    -c custom.config \
    --rnaseq_dir "/workspace/data/rnaseq" \
    --genome_fasta "/workspace/genome.fasta" \
    --masked_genome "/workspace/genome_masked.fasta" \
    --protein_evidence "/workspace/proteins.fa" \
    --adapters "/workspace/adapters.fa" \
    --outdir "/workspace/results" \
    --species "YourSpeciesName" \
    --threads 16 \
    --memory "32G" \
    --busco_db "protists" \
    --busco_lineage "protists_odb10" \
    --run_mode prediction_only \
    -with-report -with-timeline -with-dag
````
2. **Prediction Only (Without RNA-Seq)**:
````bash
nextflow run genome_annotation.nf \
    -c custom.config \
    --genome_fasta "/workspace/genome.fasta" \
    --masked_genome "/workspace/genome_masked.fasta" \
    --protein_evidence "/workspace/proteins.fa" \
    --outdir "/workspace/results" \
    --species "YourSpeciesName" \
    --threads 16 \
    --memory "32G" \
    --busco_db "protists" \
    --busco_lineage "protists_odb10" \
    --run_mode prediction_only \
    -with-report
````
3. **Annotation Only (Using BRAKER3 Prediction)**:
````bash
nextflow run genome_annotation.nf \
    -c custom.config \
    --outdir "/workspace/results" \
    --genome_fasta "/workspace/genome.fasta" \
    --species "YourSpeciesName" \
    --funannotate_db "/opt/funannotate_db" \
    --busco_db "protists" \
    --threads 16 \
    --memory "32G" \
    --run_mode annotation_only \
    --enable_functional_annotation true \
    --annotation_prediction braker \
    -with-report
````
4. **Both Prediction and Annotation (Without RNA-Seq, Best BUSCO)**:
````bash
nextflow run genome_annotation.nf \
    -c custom.config \
    --genome_fasta "/workspace/genome.fasta" \
    --masked_genome "/workspace/genome_masked.fasta" \
    --protein_evidence "/workspace/proteins.fa" \
    --outdir "/workspace/results" \
    --species "YourSpeciesName" \
    --threads 16 \
    --memory "32G" \
    --busco_db "protists" \
    --busco_lineage "protists_odb10" \
    --run_mode both \
    --enable_functional_annotation true \
    --annotation_prediction best_busco \
    -with-report
````
5. **Custom Mount Point**: To use a different Docker mount point (e.g., ```/custom_mount```):
````bash
nextflow run genome_annotation.nf \
    -c custom.config \
    --genome_fasta "/custom_mount/genome.fasta" \
    --outdir "/custom_mount/results" \
    --docker_mount_base "/custom_mount" \
    [...]
````
**Optional Parameters**
- ```--gc_probability```: Set BRAKER3 GC probability (e.g., ```0.6377```). Default: ```false``` (omitted).
- ```--star_index_nbases```: STAR genome index parameter (default: ```10```). Only used with RNA-Seq.
- ```--bam_sort_memory```: Memory for BAM sorting (default: ```10G```). Only used with RNA-Seq.
- ```--min_protlen```: Minimum protein length for funannotate (default: ```100```).

## Outputs
Outputs are saved in ```--outdir``` (e.g., ```/workspace/results```):
- **Prediction**:
   - ```braker3/```: BRAKER3 predictions, including ```braker.gtf``` (if generated).
   - ```funannotate_update/```: Funannotate predictions, including GFF3 and proteins.
   - ```merged_annotations/fixed_merged.gff```: BRAKER3 GFF with tRNA annotations, fixed for overlaps (if BRAKER3 GTF exists).
   - ```busco/```: BUSCO results for BRAKER3 and funannotate, with ```busco_comparison.txt``` summarizing completeness scores.
   - ```cleaned_reads/```, ```mapped/```, ```stringtie/```, ```transcripts/``` (only if RNA-Seq provided).
- **Annotation** (if ```--enable_functional_annotation true```):
   - ```functional_annotation/anno/```: Funannotate annotations.
   - ```functional_annotation/phobius/```, ```interproscan/```, ```eggnog/```: Tool-specific results.
- **Logs**: ```report.html```, ```timeline.html```, ```dag.dot``` for pipeline execution details.

## Run Modes

- **prediction_only**: Runs BRAKER3 and funannotate predictions, compares BUSCO scores.
- **annotation_only**: Performs functional annotation on existing BRAKER3 or funannotate predictions (specified by ```--annotation_prediction```).
- **both**: Runs predictions and annotates the user-specified or best BUSCO-scored prediction.

## Notes

- **RNA-Seq Absence**: If ```--rnaseq_dir``` is not specified or the directory doesn’t exist, the pipeline skips RNA-Seq processing (```CleanReads```, ```IndexGenome```, ```MapReads```, ```MergeBAMs```, ```StringTie```, ```GTFtoFasta```, ```FunannotateTrain```) and runs predictions using only genome, masked genome, and protein evidence. BRAKER3 and funannotate will omit RNA-Seq-related options (e.g., ```--bam```, ```--stringtie```).
- **BRAKER3 GTF**: Without RNA-Seq, BRAKER3 may still produce a GTF using protein evidence, but UTR addition is skipped.
- **Funannotate**: Without RNA-Seq, funannotate relies on protein evidence and tRNA annotations, which may reduce prediction accuracy.
- **BUSCO Comparison**: Always performed to compare BRAKER3 and funannotate predictions, even without RNA-Seq.

## Troubleshooting

- **Docker Errors**: Ensure Docker is running and your user is in the ```docker``` group (```sudo usermod -aG docker $USER```).
- **File Not Found**: Verify input paths are correct relative to ```--docker_mount_base```. Check file permissions.
- **BRAKER3 Fails**: Ensure ```--busco_lineage``` matches your organism. Check ```results/braker3/``` logs. Without RNA-Seq, predictions may be less accurate.
- **Resource Issues**: Increase ```--threads``` or ```--memory``` for large datasets. Monitor ```report.html```.
- **Funannotate Database**: Ensure --funannotate_db points to a valid database directory.

## License

This project is licensed under the MIT License. See the LICENSE file for details.

## Acknowledgments
This work is supported through the Sequencing analysis (SequAna) core facility at the University of Konstanz [https://www.biologie.uni-konstanz.de/sequana/]

**Tools**: BRAKER3, funannotate, BUSCO, Trimmomatic, STAR, StringTie, tRNAscan-SE, Phobius, InterProScan, eggNOG-mapper.
**Nextflow**: For workflow orchestration.
**Docker**: For reproducible environments.


For more information or help, please contact [abdoallah.sharaf@uni-konstanz.de].

This is the complete README in markdown format, covering everything from installation to usage, including file formats and options. Feel free to let me know if you'd like to make any additional changes or if there's something else you'd like to add!









