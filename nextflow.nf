$HOSTNAME = ""
params.outdir = 'results'  

params.publishdict = [:]

def pathChecker(input, path, type){
	def cmd = "mkdir -p check && mv ${input} check/. "
	if (!input || input.empty()){
		input = file(path).getName().toString()
		cmd = "mkdir -p check && cd check && ln -s ${path} ${input} && cd .."
		if (path.indexOf('s3:') > -1 || path.indexOf('S3:') >-1){
			recursive = (type == "folder") ? "--recursive" : ""
			cmd = "mkdir -p check && cd check && aws s3 cp ${recursive} ${path} ${workDir}/${input} && ln -s ${workDir}/${input} . && cd .."
		} else if (path.indexOf('gs:') > -1 || path.indexOf('GS:') >-1){
			if (type == "folder"){
				cmd = "mkdir -p check ${workDir}/${input} && cd check && gsutil rsync -r ${path} ${workDir}/${input} && cp -R ${workDir}/${input} . && cd .."
			} else {
				cmd = "mkdir -p check && cd check && gsutil cp ${path} ${workDir}/${input} && cp -R ${workDir}/${input} . && cd .."
			}
		} else if (path.indexOf('/') == -1){
			cmd = ""
		}
}
	return [cmd,input]
}
if (!params.feature_reference){params.feature_reference = ""} 
if (!params.VDJ_reference){params.VDJ_reference = ""} 
if (!params.reads){params.reads = ""} 
if (!params.mate){params.mate = ""} 
if (!params.custom_additional_genome){params.custom_additional_genome = ""} 
if (!params.Metadata){params.Metadata = ""} 
if (!params.cmo_set){params.cmo_set = ""} 
if (!params.custom_additional_gtf){params.custom_additional_gtf = ""} 
if (!params.mask_gtf){params.mask_gtf = ""} 
if (!params.db_feather){params.db_feather = ""} 
if (!params.motif_db){params.motif_db = ""} 
if (!params.tf_lists){params.tf_lists = ""} 
// Stage empty file to be used as an optional input where required
ch_empty_file_1 = file("$baseDir/.emptyfiles/NO_FILE_1", hidden:true)
ch_empty_file_2 = file("$baseDir/.emptyfiles/NO_FILE_2", hidden:true)
ch_empty_file_3 = file("$baseDir/.emptyfiles/NO_FILE_3", hidden:true)
ch_empty_file_4 = file("$baseDir/.emptyfiles/NO_FILE_4", hidden:true)
ch_empty_file_5 = file("$baseDir/.emptyfiles/NO_FILE_5", hidden:true)
ch_empty_file_6 = file("$baseDir/.emptyfiles/NO_FILE_6", hidden:true)
ch_empty_file_7 = file("$baseDir/.emptyfiles/NO_FILE_7", hidden:true)

g_11_2_g_20 = params.feature_reference && file(params.feature_reference, type: 'any').exists() ? file(params.feature_reference, type: 'any') : ch_empty_file_2
g_12_3_g_20 = params.VDJ_reference && file(params.VDJ_reference, type: 'any').exists() ? file(params.VDJ_reference, type: 'any') : ch_empty_file_3
if (params.reads){
Channel
	.fromFilePairs( params.reads,checkExists:true , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 ) 
	.set{g_23_0_g_22}
  } else {  
	g_23_0_g_22 = Channel.empty()
 }

Channel.value(params.mate).set{g_24_1_g_22}
g_40_2_g50_58 = params.custom_additional_genome && file(params.custom_additional_genome, type: 'any').exists() ? file(params.custom_additional_genome, type: 'any') : ch_empty_file_1
g_43_1_g51_0 = params.Metadata && file(params.Metadata, type: 'any').exists() ? file(params.Metadata, type: 'any') : ch_empty_file_1
g_68_9_g_20 = params.cmo_set && file(params.cmo_set, type: 'any').exists() ? file(params.cmo_set, type: 'any') : ch_empty_file_5
g_69_3_g50_58 = params.custom_additional_gtf && file(params.custom_additional_gtf, type: 'any').exists() ? file(params.custom_additional_gtf, type: 'any') : ch_empty_file_2
g_71_1_g70_1 = params.mask_gtf && file(params.mask_gtf, type: 'any').exists() ? file(params.mask_gtf, type: 'any') : ch_empty_file_1
g_75_1_g82_8 = file(params.db_feather, type: 'any')
g_76_2_g82_8 = file(params.motif_db, type: 'any')
g_77_1_g82_1 = file(params.tf_lists, type: 'any')

//* @style @array:{bcl_directory,sampleSheet} @multicolumn:{bcl_directory,sampleSheet}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 1
}
//* platform
//* platform
//* autofill

process Demultiplexer_prep {


output:
 val bcl_and_sampleSheets  ,emit:g_4_bcl00_g_0 
 val bcl_directory  ,emit:g_4_bcl_directory16_g_20 

container 'quay.io/viascientific/pipeline_base_image:1.0'

when:
params.run_bclConvert == "yes"

exec:
bcl_directory = params.Demultiplexer_prep.bcl_directory
sampleSheet = params.Demultiplexer_prep.sampleSheet
bcl_directory2 = bcl_directory.collect{ '"' + it + '"'}
sampleSheet2 = sampleSheet.collect{ '"' + it + '"'}
bcl_and_sampleSheets = [ bcl_directory, sampleSheet ].transpose().flatten()


}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 8
    $MEMORY = 50
}
//* platform
//* platform
//* autofill

process bclConvert {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${bcl_files}_fastq$/) "bclConvert/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${bcl_files}_reports$/) "bclConvert_report/$filename"}
input:
 tuple file(bcl_files),file(samplesheet)

output:
 path "${bcl_files}_fastq"  ,emit:g_0_reads00_g_33 
 path "${bcl_files}_reports"  ,emit:g_0_outputDir11 

disk { 1500.GB * task.attempt }

container 'quay.io/viascientific/bclconvert:4.3.6'

when:
params.run_bclConvert == "yes"
script:
bclconvert_parameters = params.bclConvert.bclconvert_parameters
SampleSheetLocation = samplesheet.name.startsWith('input.') ? "${bcl_files}/SampleSheet.csv" : "${samplesheet}"
"""	
bcl-convert ${bclconvert_parameters} --bcl-input-directory ${bcl_files} --output-directory fastq --sample-sheet ${SampleSheetLocation}
mv fastq ${bcl_files}_fastq
mkdir ${bcl_files}_reports 
find .
mv ${bcl_files}_fastq/Reports ${bcl_files}_reports/.
mv ${bcl_files}_fastq/Logs ${bcl_files}_reports/.
"""

}


process cellranger_fastq_prep {

input:
 tuple val(name), file(reads)
 val mate

output:
 path "reads/*"  ,emit:g_22_reads00_g_25 

disk 300.GB 

script:
bcl_directory = ["/reads"]
sample = name
nameAll = reads.toString()
nameArray = nameAll.split(' ')
if (mate == "pair"){
    read1 = nameArray[0]
    read2 = nameArray[1]
    if (nameAll.contains('.gz')) {
        runGzip = ''
    } else {
        read1 = read1 + ".gz"
        read2 = read2 + ".gz"
        runGzip = "ls * | xargs -i echo gzip -f {} | sh"
    }
    mvReads = "cp " +read1 +" reads/"+ name + "_S1_L001_R1_001.fastq.gz && cp " +read2 +" reads/"+ name + "_S1_L001_R2_001.fastq.gz"  
} else {
    read1 = nameArray[0]
    if (nameAll.contains('.gz')) {
        runGzip = ''
    } else {
        read1 = read1 + ".gz"
        runGzip = "ls * | xargs -i echo gzip -f {} | sh"
    }
    mvReads = "cp " +read1 +" reads/"+ name + "_S1_L001_R1_001.fastq.gz"
}

"""
$runGzip
mkdir reads
$mvReads



"""
}


process cellranger_fastq_collect {

input:
 path reads

output:
 val bcl_directory  ,emit:g_25_bcl_directory05_g_20 
 path "reads_fastq"  ,emit:g_25_reads17_g_20 

disk { 1000.GB * task.attempt }
errorStrategy 'retry'
maxRetries 2

script:
bcl_directory = ["/reads"]
"""	
mkdir reads_fastq
mv $reads reads_fastq/.
"""
}


process flatten_cellranger_reads {

input:
 path reads

output:
 val "single"  ,emit:g_33_mate00_g_35 
 path "allreads/*/*_{R1,R2}_001.fastq.gz", includeInputs: true   ,emit:g_33_reads10_g_34 

stageInMode 'copy'
disk { 1000.GB * task.attempt }
errorStrategy 'retry'
maxRetries 2

when:
(params.run_FastQC && (params.run_FastQC == "yes"))

script:
"""
# rename Undetermined fastq's to prevent name clash in multiqc 
for r in */Undetermined*L00*_R*; do
    rname=\$(basename \$r)
    dirname=\${r%/*}
    prefix=\${r%_fastq/*}
    mv \$r \$dirname/\${prefix}_\${rname}
done

mkdir allreads
mv ${reads} allreads/.
"""

}


process file_to_set_conversion_for_reads {

input:
 path reads

output:
 tuple val(name),file(reads)  ,emit:g_34_reads01_g_35 

script:
name = reads.baseName
"""
echo "done"	
"""
}



process FastQC_after_mkfastq {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.(html|zip)$/) "fastqc/$filename"}
input:
 val mate
 tuple val(name), file(reads)

output:
 path '*.{html,zip}'  ,emit:g_35_FastQCout04_g_52 
 tuple val(name), file("reads/*")  ,emit:g_35_reads11 

when:
(params.run_FastQC && (params.run_FastQC == "yes"))

script:
"""
fastqc ${reads} 
mkdir reads
mv ${reads} reads/.
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 24
}
//* platform
//* platform
//* autofill

process MultiQC {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /multiqc_report.html$/) "multiqc/$filename"}
input:
 path "fastqc/*"

output:
 path "multiqc_report.html" ,optional:true  ,emit:g_52_outputHTML00 
 path "*" ,optional:true  ,emit:g_52_outputDir11 

container 'quay.io/biocontainers/multiqc:1.25.2--pyhdfd78af_0'

script:
multiqc_parameters = params.MultiQC.multiqc_parameters
plots_flat_numseries = params.MultiQC.plots_flat_numseries

"""
multiqc ${multiqc_parameters}  -d -dd 2 --cl-config "plots_flat_numseries: ${plots_flat_numseries}" .

"""


}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 3
}
//* autofill

process cellranger_multi_library_prep {


output:
 val librarySettings  ,emit:g_67_librarySettings00_g_5 

when:
params.run_cellranger_multi == "yes"

exec:
//libraries
def group;
fastq_id = params.cellranger_multi_library_prep.fastq_id
group = params.cellranger_multi_library_prep.group
feature_types = params.cellranger_multi_library_prep.feature_types
librarySettings = [:]
librarySettings["fastq_id"] = fastq_id.collect{ '"' + it + '"'}
librarySettings["group_libraries"] = group.collect{ '"' + it + '"'}
librarySettings["feature_types"] = feature_types.collect{ '"' + it + '"'}

//* @spreadsheet:{fastq_id,group,feature_types}

}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 3
}
//* autofill

process cellranger_multi_prep {

input:
 val librarySettings

output:
 path "*.csv"  ,emit:g_5_csvFile04_g_20 
 val settings  ,emit:g_5_settings18_g_20 

when:
(params.run_cellranger_multi && (params.run_cellranger_multi == "yes")) || !params.run_cellranger_multi

script:

settings = [:]
settings["fastq_id"] = librarySettings["fastq_id"]
settings["group_libraries"] = librarySettings["group_libraries"]
settings["feature_types"] = librarySettings["feature_types"]
def group;
//samples
sample_id = params.cellranger_multi_prep.sample_id
group = params.cellranger_multi_prep.group
cmo_ids = params.cellranger_multi_prep.cmo_ids
hashtag_ids = params.cellranger_multi_prep.hashtag_ids
ocm_barcode_ids = params.cellranger_multi_prep.ocm_barcode_ids
probe_barcode_ids = params.cellranger_multi_prep.probe_barcode_ids
description = params.cellranger_multi_prep.description

//gene-expression
create_bam = params.cellranger_multi_prep.create_bam
expect_cells = params.cellranger_multi_prep.expect_cells
force_cells = params.cellranger_multi_prep.force_cells
cellranger_multi_chemistry = params.cellranger_multi_prep.cellranger_multi_chemistry
r1_length = params.cellranger_multi_prep.r1_length
check_library_compatibility = params.cellranger_multi_prep.check_library_compatibility
include_introns = params.cellranger_multi_prep.include_introns
//feature
r1_length_feature = params.cellranger_multi_prep.r1_length_feature
//vdj
r1_length_vdj = params.cellranger_multi_prep.r1_length_vdj

settings["create_bam"] = create_bam
settings["expect_cells"] = expect_cells
settings["force_cells"] = force_cells
settings["cellranger_multi_chemistry"] = cellranger_multi_chemistry
settings["r1_length"] = r1_length
settings["check_library_compatibility"] = check_library_compatibility
settings["include_introns"] = include_introns
settings["r1_length_feature"] = r1_length_feature
settings["r1_length_vdj"] = r1_length_vdj

group_samples = []
if (group instanceof List) {
	group_samples = group.collect{ '"' + it + '"'}	
} 
description2 = description.collect{ '"' + it + '"'}

settings["sample_id"]=sample_id.collect{ '"' + it + '"'}
settings["cmo_ids"]=cmo_ids.collect{ '"' + it + '"'}
settings["hashtag_ids"]=hashtag_ids.collect{ '"' + it + '"'}
settings["ocm_barcode_ids"]=ocm_barcode_ids.collect{ '"' + it + '"'}
settings["probe_barcode_ids"]=probe_barcode_ids.collect{ '"' + it + '"'}
settings["group_samples"]=group_samples
settings["description"]=description2


//* @style @multicolumn:{create_bam,expect_cells,force_cells,cellranger_multi_chemistry,r1_length, check_library_compatibility,include_introns}  @spreadsheet:{sample_id,group,cmo_ids,hashtag_ids,ocm_barcode_ids,probe_barcode_ids,description} 
"""
#!/usr/bin/env python

import subprocess,sys
import csv,os  

def is_not_blank(s):
	return bool(s and not s.isspace() and s != "null" and not s.startswith('NO_FILE'))

lanes = ${settings["group_libraries"]}
lane_of_samples = ${settings["group_samples"]}

# Group config files by using lanes:  lanes + lane_of_samples
all_lanes = lanes + lane_of_samples
unique_lanes = list(set(all_lanes))
print(unique_lanes)
# remove "" from list if it has "all" in it
if '' in unique_lanes :
	unique_lanes.remove('')
	if 'all' not in unique_lanes :
		unique_lanes.append("all")

# if unique_lanes has "all" and it has more than 1 lanes then remove "all"
# other lanes will use "all"
if (len(unique_lanes) > 1) and ('all' in unique_lanes):
	unique_lanes.remove('all')
	
print(unique_lanes)

for l in range(len(unique_lanes)):
	f = open('config_'+unique_lanes[l]+'.csv', "w")
	

"""



}

//* params.gtf =  ""  //* @input
//* params.genome =  ""  //* @input
//* params.commondb =  ""  //* @input
//* params.genome_source =  ""  //* @input
//* params.gtf_source =  ""  //* @input
//* params.commondb_source =  ""  //* @input @optional

def downFile(path, task){
	println workDir
    if (path.take(1).indexOf("/") == 0){
      target=path
      if (task.executor == "awsbatch" || task.executor == "google-batch") {
      	a=file(path)
    	fname = a.getName().toString()
    	target = "${workDir}/${fname}"
    	if (!file(target).exists()){
    		a.copyTo(workDir)
    	}
      }
    } else {
      a=file(path)
      fname = a.getName().toString()
      target = "${workDir}/${fname}"
      if (!file(target).exists()){
    		a.copyTo(workDir)
      } 
    }
    return target
}

def getLastName (str){
	if (str.indexOf("/") > -1){
		return  str.substring(str.lastIndexOf('/')+1,str.length())
	} 
	return ""
}

process Check_and_Build_Module_Check_Genome_GTF {


output:
 path "${newNameFasta}"  ,emit:g50_21_genome00_g50_58 
 path "${newNameGtf}"  ,emit:g50_21_gtfFile10_g50_57 

container "${ params.IMAGE_BASE ? "${params.IMAGE_BASE}/pipeline_base_image:1.0" : "quay.io/viascientific/pipeline_base_image:1.0" }"


when:
params.run_Download_Genomic_Sources == "yes"

script:
genomeSource = !file("${params.genome}").exists() ? params.genome_source : params.genome
genomeName = getLastName(genomeSource)

gtfSource = !file("${params.gtf}").exists() ? params.gtf_source : params.gtf
gtfName = getLastName(gtfSource)


newNameGtf = gtfName
newNameFasta = genomeName
if (gtfName.contains('.gz')) { newNameGtf =  newNameGtf - '.gz'  } 
if (genomeName.contains('.gz')) { newNameFasta =  newNameFasta - '.gz'  } 

runGzip = ""
if (gtfName.contains('.gz') || genomeName.contains('.gz')) {
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} 

slashCountGenome = params.genome_source.count("/")
cutDirGenome = slashCountGenome - 3;

slashCountGtf = params.gtf_source.count("/")
cutDirGtf = slashCountGtf - 3;

"""
if [ ! -e "${params.genome_source}" ] ; then
    echo "${params.genome_source} not found"
	if [[ "${params.genome_source}" =~ "s3" ]]; then
		echo "Downloading s3 path from ${params.genome_source}"
		aws s3 cp ${params.genome_source} ${workDir}/${genomeName} && ln -s ${workDir}/${genomeName} ${genomeName}
	elif [[ "${params.genome_source}" =~ "gs" ]]; then
		echo "Downloading gs path from ${params.genome_source}"
		gsutil cp  ${params.genome_source} ${genomeName}
	else
		echo "Downloading genome with wget"
		wget --no-check-certificate --secure-protocol=TLSv1 -l inf -nc -nH --cut-dirs=$cutDirGenome -R 'index.html*' -r --no-parent  ${params.genome_source}
	fi

else 
	ln -s ${params.genome_source} ${genomeName}
fi

if [ ! -e "${params.gtf_source}" ] ; then
    echo "${params.gtf_source} not found"
	if [[ "${params.gtf_source}" =~ "s3" ]]; then
		echo "Downloading s3 path from ${params.gtf_source}"
		aws s3 cp  ${params.gtf_source} ${workDir}/${gtfName} && ln -s ${workDir}/${gtfName} ${gtfName}
	elif [[ "${params.gtf_source}" =~ "gs" ]]; then
		echo "Downloading gs path from ${params.gtf_source}"
		gsutil cp  ${params.gtf_source} ${gtfName}
	else
		echo "Downloading gtf with wget"
		wget --no-check-certificate --secure-protocol=TLSv1 -l inf -nc -nH --cut-dirs=$cutDirGtf -R 'index.html*' -r --no-parent  ${params.gtf_source}
	fi

else 
	ln -s ${params.gtf_source} ${gtfName}
fi

$runGzip

"""




}


process Check_and_Build_Module_convert_gtf_attributes {

input:
 path gtf

output:
 path "out/${gtf}"  ,emit:g50_57_gtfFile01_g50_58 

when:
params.replace_geneID_with_geneName == "yes"

script:
"""
#!/usr/bin/env perl 

## Replace gene_id column with gene_name column in the gtf file
## Also check if any transcript_id defined in multiple chromosomes.
system("mkdir out");

open(my \$out_valid, ">out/${gtf}");
open(my \$out_notvalid, ">notvalid_${gtf}");
my \$fileName = "${gtf}";
# To track transcript IDs and ensure each appears on only one chromosome.
my %transcript_chrom;

# Open input file.
open(my \$in, "<", \$fileName) or die "Cannot open \$fileName: \$!";

while (my \$line = <\$in>) {
    chomp \$line;
    my @fields = split(/\\t/, \$line);
	my \$feature = \$fields[2];
	
    if (@fields < 9) {
        print \$out_notvalid "\$line\n";
        next;
    }

    # Split the attribute column by semicolons.
    my @attrs = split(/;/, \$fields[8]);
    # Remove any empty elements (if present).
    @attrs = grep { /\\S/ } @attrs;

    # Create a hash for easy lookup of attributes.
    my %attr_hash;
    foreach my \$attr (@attrs) {
        \$attr =~ s/^\\s+|\\s+\$//g;  # trim leading and trailing whitespace
        if ( \$attr =~ /^(\\S+)\\s+"([^"]+)"/ ) {
            my (\$key, \$value) = (\$1, \$2);
            \$attr_hash{\$key} = \$value;
        }
    }

    # Determine the gene_id value: use gene_name if available.
    my \$geneId = "";
    if (exists \$attr_hash{"gene_name"}) {
        \$geneId = \$attr_hash{"gene_name"};
    } elsif (exists \$attr_hash{"gene_id"}) {
        \$geneId = \$attr_hash{"gene_id"};
    }

    # Determine transcript_id from available attributes.
    my \$transcript_id = "";
    if (exists \$attr_hash{"transcript_id"} && \$attr_hash{"transcript_id"} ne "") {
        \$transcript_id = \$attr_hash{"transcript_id"};
    } elsif (exists \$attr_hash{"transcript_name"}) {
        \$transcript_id = \$attr_hash{"transcript_name"};
    } elsif (exists \$attr_hash{"gene_id"}) {
        \$transcript_id = \$attr_hash{"gene_id"};
    } 
    
    # If both geneId and transcript_id were found, continue processing.
    if (\$geneId ne "" && \$transcript_id ne "") {
        # Check that the same transcript_id is not found on multiple chromosomes.
        if (exists \$transcript_chrom{\$transcript_id}) {
            if (\$transcript_chrom{\$transcript_id} ne \$fields[0]) {
                print \$out_notvalid "\$transcript_id: \$transcript_chrom{\$transcript_id} vs \$fields[0]\\n";
                next;
            }
        } else {
            \$transcript_chrom{\$transcript_id} = \$fields[0];
        }

        # Remove any existing gene_id and transcript_id attributes from the list.
        @attrs = grep { !/^\\s*(gene_id|transcript_id)\\b/ } @attrs;

        # Prepend the updated gene_id and transcript_id attributes.
        # if \$feature = "gene" no need to add transcript_id column (which will be gene_id)
        if (\$feature ne "gene"){
        	unshift @attrs, 'transcript_id "' . \$transcript_id . '"';
        }
        unshift @attrs, 'gene_id "' . \$geneId . '"';

        # Reassemble the attribute field (adding a semicolon after each item).
        \$fields[8] = join("; ", @attrs) . ";";
        # --------------------------------------------------------------------

        # Print the updated line to the valid output file.
        print \$out_valid join("\\t", @fields) . "\\n";
    } else {
        # If either geneId or transcript_id is missing, output to the not valid file.
        print \$out_notvalid "\$line\\n";
    }
}

close \$in;
close \$out_valid;
close \$out_notvalid;
"""
}


process Check_and_Build_Module_Add_custom_seq_to_genome_gtf {

input:
 path genome
 path gtf
 path custom_fasta
 path custom_gtf

output:
 path "${genomeName}_custom.fa"  ,emit:g50_58_genome00_g50_52 
 path "${gtfName}_custom_sorted.gtf"  ,emit:g50_58_gtfFile10_g50_53 

container "${ params.IMAGE_BASE ? "${params.IMAGE_BASE}/custom_sequence_to_genome_gtf:1.0" : "quay.io/viascientific/custom_sequence_to_genome_gtf:1.0" }"

when:
params.add_sequences_to_reference == "yes"

script:
genomeName = genome.baseName
gtfName = gtf.baseName
is_custom_genome_exists = custom_fasta.name.startsWith('NO_FILE') ? "False" : "True" 
is_custom_gtf_exists = custom_gtf.name.startsWith('NO_FILE') ? "False" : "True" 
"""
#!/usr/bin/env python 
import requests
import os
import pandas as pd
import re
import urllib
from Bio import SeqIO

def add_to_fasta(seq, sqid, out_name):
	new_line = '>' + sqid + '\\n' + seq + '\\n'
	with open(out_name + '.fa', 'a') as f:
		f.write(new_line)

def createCustomGtfFromFasta(fastaFile, outCustomGtfFile):

    fasta_sequences = SeqIO.parse(open(fastaFile),'fasta')
    with open(outCustomGtfFile, "w") as out_file:
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            last = len(sequence)
            line1 = "{gene}\\tKNOWN\\tgene\\t{first}\\t{last}\\t.\\t+\\t.\\tgene_id \\"{gene}\\"; gene_version \\"1\\"; gene_type \\"protein_coding\\"; gene_source \\"KNOWN\\"; gene_name \\"{gene}\\"; gene_biotype \\"protein_coding\\"; gene_status \\"KNOWN\\"; level 1;".format(gene=name, first="1", last=last)
            line2 = "{gene}\\tKNOWN\\ttranscript\\t{first}\\t{last}\\t.\\t+\\t.\\tgene_id \\"{gene}\\"; gene_version \\"1\\"; transcript_id \\"{gene}_trans\\"; transcript_version \\"1\\"; gene_type \\"protein_coding\\"; gene_source \\"KNOWN\\"; transcript_source \\"KNOWN\\"; gene_status \\"KNOWN\\"; gene_name \\"{gene}\\"; gene_biotype \\"protein_coding\\"; transcript_type \\"protein_coding\\"; transcript_status \\"KNOWN\\"; transcript_name \\"{gene}_1\\"; level 1; tag \\"basic\\"; transcript_biotype \\"protein_coding\\"; transcript_support_level \\"1\\";".format(gene=name, first="1", last=last)
            line3 = "{gene}\\tKNOWN\\texon\\t{first}\\t{last}\\t.\\t+\\t.\\tgene_id \\"{gene}\\"; gene_version \\"1\\"; transcript_id \\"{gene}_trans\\"; transcript_version \\"1\\"; exon_number 1; gene_type \\"protein_coding\\"; gene_source \\"KNOWN\\"; transcript_source \\"KNOWN\\"; gene_status \\"KNOWN\\"; gene_name \\"{gene}\\"; gene_biotype \\"protein_coding\\"; transcript_type \\"protein_coding\\"; transcript_status \\"KNOWN\\"; transcript_biotype \\"protein_coding\\"; transcript_name \\"{gene}_1\\"; exon_number 1; exon_id \\"{gene}.1\\"; level 1; tag \\"basic\\"; transcript_support_level \\"1\\";".format(gene=name, first="1", last=last)
            out_file.write("{}\\n{}\\n{}\\n".format(line1, line2, line3))

	
os.system('cp ${genomeName}.fa ${genomeName}_custom.fa')  
os.system('cp ${gtfName}.gtf ${gtfName}_custom.gtf')  

if ${is_custom_genome_exists}:
	os.system("tr -d '\\r' < ${custom_fasta} > ${custom_fasta}_tmp && rm ${custom_fasta} && mv ${custom_fasta}_tmp ${custom_fasta}")
	os.system('cat ${custom_fasta} >> ${genomeName}_custom.fa')
	if ${is_custom_gtf_exists}:
		os.system("tr -d '\\r' < ${custom_gtf} > ${custom_gtf}_tmp && rm ${custom_gtf} && mv ${custom_gtf}_tmp ${custom_gtf}")
		os.system("mv ${custom_gtf} ${custom_fasta}.gtf")
	else:
		createCustomGtfFromFasta("${custom_fasta}", "${custom_fasta}.gtf")
	os.system('cat ${custom_fasta}.gtf >> ${gtfName}_custom.gtf')

	
os.system('samtools faidx ${genomeName}_custom.fa')
os.system('igvtools sort ${gtfName}_custom.gtf ${gtfName}_custom_sorted.gtf')
os.system('igvtools index ${gtfName}_custom_sorted.gtf')

"""
}

//* params.gtf2bed_path =  ""  //* @input
//* params.bed =  ""  //* @input

process Check_and_Build_Module_Check_BED12 {

input:
 path gtf

output:
 path "${gtfName}.bed"  ,emit:g50_53_bed03_g50_54 

container "${ params.IMAGE_BASE ? "${params.IMAGE_BASE}/rnaseq:4.0" : "quay.io/viascientific/rnaseq:4.0" }"

when:
params.run_Download_Genomic_Sources == "yes"

script:
gtfName  = gtf.baseName
beddir = ""
if (params.bed.indexOf('/') > -1){
	beddir  = params.bed.substring(0, params.bed.lastIndexOf('/')) 
}
"""

if [ ! -e "${params.bed}" ] ; then
    echo "${params.bed} not found"
    perl ${params.gtf2bed_path} $gtf > ${gtfName}.bed
else 
	cp -n ${params.bed} ${gtfName}.bed
fi
if [ "${beddir}" != "" ] ; then
	mkdir -p ${beddir}
	cp -n ${gtfName}.bed ${params.bed} 
fi
"""




}

//* params.gtf2bed_path =  ""  //* @input
//* params.genome_sizes =  ""  //* @input

process Check_and_Build_Module_Check_chrom_sizes_and_index {

input:
 path genome

output:
 path "${genomeName}.chrom.sizes"  ,emit:g50_52_genomeSizes02_g50_54 

when:
params.run_Download_Genomic_Sources == "yes"

script:
genomeName  = genome.baseName
genome_sizes_dir = ""
if (params.genome_sizes.indexOf('/') > -1){
	genome_sizes_dir  = params.genome_sizes.substring(0, params.genome_sizes.lastIndexOf('/')) 
}

"""
if [ ! -e "${params.genome_sizes}" ] ; then
    echo "${params.genome_sizes} not found"
    cat ${genome} | awk '\$0 ~ ">" {print c; c=0;printf substr(\$1,2,100) "\\t"; } \$0 !~ ">" {c+=length(\$0);} END { print c; }' > ${genomeName}.chrom.sizes
    ##clean first empty line
    sed -i '1{/^\$/d}' ${genomeName}.chrom.sizes
    if [ "${genome_sizes_dir}" != "" ] ; then
    	mkdir -p ${genome_sizes_dir}
		cp -n ${genomeName}.chrom.sizes ${params.genome_sizes} 
	fi
else 
	cp ${params.genome_sizes} ${genomeName}.chrom.sizes
fi

"""




}


process Check_and_Build_Module_check_files {

input:
 path gtf
 path genome
 path genomeSizes
 path bed

output:
 path "*/${gtf2}" ,optional:true  ,emit:g50_54_gtfFile01_g_7 
 path "*/${genome2}" ,optional:true  ,emit:g50_54_genome10_g_7 
 path "*/${genomeSizes2}" ,optional:true  ,emit:g50_54_genomeSizes22 
 path "*/${bed2}" ,optional:true  ,emit:g50_54_bed33 

container "${ params.IMAGE_BASE ? "${params.IMAGE_BASE}/pipeline_base_image:1.0" : "quay.io/viascientific/pipeline_base_image:1.0" }"

stageInMode 'copy'

script:
(cmd1, gtf2) = pathChecker(gtf, params.gtf, "file")
(cmd2, genome2) = pathChecker(genome, params.genome, "file")
(cmd3, genomeSizes2) = pathChecker(genomeSizes, params.genome_sizes, "file")
(cmd4, bed2) = pathChecker(bed, params.bed, "file")
"""
$cmd1
$cmd2
$cmd3
$cmd4
"""
}

//* params.gtf_type = 'ensembl'  //*  @dropdown @options:"ensembl","gencode"

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 4
    $MEMORY = 64
}
//* platform
//* platform
//* autofill

process cellranger_mkref {

input:
 path genome
 path gtf

output:
 path "ref"  ,emit:g_7_reference00_g_18 

when:
(params.run_mkref && (params.run_mkref == "yes")) || !params.run_mkref

script:
optional_mkgtf_filtering_parameters = params.cellranger_mkref.optional_mkgtf_filtering_parameters
"""
if [ "${params.gtf_type}" == "gencode" ]; then
	# Define string patterns for GTF tags
# NOTES:
# - Since GENCODE release 31/M22 (Ensembl 97), the "lincRNA" and "antisense"
#   biotypes are part of a more generic "lncRNA" biotype.
# - These filters are relevant only to GTF files from GENCODE. The GTFs from
#   Ensembl release 98 have the following differences:
#   - The names "gene_biotype" and "transcript_biotype" are used instead of
#     "gene_type" and "transcript_type".
#   - Readthrough transcripts are present but are not marked with the
#     "readthrough_transcript" tag.
#   - Only the X chromosome versions of genes in the pseudoautosomal regions
#     are present, so there is no "PAR" tag.
	BIOTYPE_PATTERN="(protein_coding|lncRNA|IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|TR_V_pseudogene|TR_J_pseudogene)"
	GENE_PATTERN="gene_type \\"\${BIOTYPE_PATTERN}\\""
	TX_PATTERN="transcript_type \\"\${BIOTYPE_PATTERN}\\""
	READTHROUGH_PATTERN="tag \\"readthrough_transcript\\""
	PAR_PATTERN="tag \\"PAR\\""


# Construct the gene ID allowlist. We filter the list of all transcripts
# based on these criteria:
#   - allowable gene_type (biotype)
#   - allowable transcript_type (biotype)
#   - no "PAR" tag (only present for Y chromosome PAR)
#   - no "readthrough_transcript" tag
# We then collect the list of gene IDs that have at least one associated
# transcript passing the filters.
	cat "${gtf}" \
    | awk '\$3 == "transcript"' \
    | grep -E "\$GENE_PATTERN" \
    | grep -E "\$TX_PATTERN" \
    | grep -Ev "\$READTHROUGH_PATTERN" \
    | grep -Ev "\$PAR_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\\1/' | sort | uniq > "gene_allowlist"

	# Filter the GTF file based on the gene allowlist
	# Copy header lines beginning with "#"
	grep -E "^#" "${gtf}" > filtered_${gtf}
	# Filter to the gene allowlist
	grep -Ff "gene_allowlist" "${gtf}" >> filtered_${gtf}
	
	cellranger mkgtf filtered_${gtf} mkgtf_${gtf} ${optional_mkgtf_filtering_parameters}
else 
	cellranger mkgtf $gtf mkgtf_${gtf} --attribute=gene_biotype:protein_coding ${optional_mkgtf_filtering_parameters}
fi

cellranger mkref --genome=ref --fasta=${genome} --genes=mkgtf_${gtf}
"""

}

//* params.transcriptome =  ""  //* @input

process cellranger_ref_checker {

input:
 path transcriptome

output:
 path "*/${transcriptome}" ,optional:true  ,emit:g_18_reference01_g_20 

container 'quay.io/viascientific/pipeline_base_image:1.0'
stageInMode 'copy'

script:
(cmd1, transcriptome) = pathChecker(transcriptome, params.transcriptome, "folder")
"""
$cmd1
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 16
    $MEMORY = 128
}
//* platform
//* platform
//* autofill



process cellranger_multi {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${run_id}_outs$/) "multi/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_web_summary.html$/) "cellranger_report/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /final_${config}$/) "multi/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*filtered_contig_annotations.csv$/) "cellranger_vdj_reports/$filename"}
input:
 path reads
 path ref
 path feature_reference
 path VDJ_reference
 path config
 val bcl_directory
 val fastq_start_directory
 path fastq_start_reads
 val settings
 path cmo_set

output:
 path "${run_id}_outs"  ,emit:g_20_outputDir00_g_59 
 path "*_web_summary.html"  ,emit:g_20_outputHTML11 
 path "final_${config}"  ,emit:g_20_csvFile22 
 path "*filtered_contig_annotations.csv" ,optional:true  ,emit:g_20_csvFile33 

disk { 1500.GB * task.attempt }
stageInMode 'copy'

when:
(params.run_cellranger_multi && (params.run_cellranger_multi == "yes")) || !params.run_cellranger_multi

script:
cpu = task.cpus - 1
memory = task.memory.toGiga() - 8
config_id = config.toString().replace(".csv","").replace("config_","")
run_id = config_id
bcl_directory2 = ""
if (bcl_directory){
	bcl_directory2 = bcl_directory.collect{ '"' + it + '"'}
} else if (fastq_start_directory){
	bcl_directory2 = fastq_start_directory.collect{ '"' + it + '"'}
}


"""
#!/usr/bin/env python

import subprocess,sys
import csv,os  
import glob

def is_not_blank(s):
	return bool(s and not s.isspace() and s != "null" and not s.startswith('NO_FILE'))

# These variables are defined in cellranger_multi_prep -> advanced tab -> Header Script section
bcl_directory = ${bcl_directory2}
fastq_id = ${settings["fastq_id"]}
lanes10x = ${settings["group_libraries"]}
feature_types = ${settings["feature_types"]}
sample_id = ${settings["sample_id"]}
cmo_ids = ${settings["cmo_ids"]} 
hashtag_ids = ${settings["hashtag_ids"]} 
ocm_barcode_ids = ${settings["ocm_barcode_ids"]} 
probe_barcode_ids = ${settings["probe_barcode_ids"]} 
lane_of_samples = ${settings["group_samples"]}
description = ${settings["description"]}


def is_nonempty_string_or_array(variable):
    if isinstance(variable, str):
        # Check if the string is not empty after stripping whitespace
        return variable.strip() != ""
    elif isinstance(variable, (list, tuple)):
        # Check if it's not empty AND at least one element is nonempty
        return len(variable) > 0 and any(str(item).strip() != "" for item in variable)
    else:
        # If it's neither a string nor a list/tuple, decide how to handle
        # For this example, we'll return False
        return False

# We can only use one of them in [samples] section, so check which one has the value 
optional_columns = [
    ("cmo_ids", cmo_ids),
    ("hashtag_ids", hashtag_ids),
    ("ocm_barcode_ids", ocm_barcode_ids),
    ("probe_barcode_ids", probe_barcode_ids),
]

nonempty_columns = [(name, col_data) for name, col_data in optional_columns if is_nonempty_string_or_array(col_data)]

# Check that only one optional column is nonempty
if len(nonempty_columns) > 1:
    raise ValueError("You can only use one of the following columns in [samples] section: cmo_ids, hashtag_ids, ocm_barcode_ids, probe_barcode_ids")

unique_lane_id = "${config_id}"

with open(r'final_${config}', 'a') as f:
	w = csv.writer(f)
	w.writerow(["[gene-expression]"])
	ref_path = os.path.abspath("${ref}")
	w.writerow(["reference",ref_path])
	if is_not_blank("${cmo_set}") and ("${cmo_set}" != "NA"):
		cmo_path = os.path.abspath("${cmo_set}")
		w.writerow(["cmo-set",cmo_path])
	if is_not_blank("${settings["expect_cells"]}"):
		w.writerow(["expect-cells","${settings["expect_cells"]}"])
	if is_not_blank("${settings["force_cells"]}"):
		w.writerow(["force-cells","${settings["force_cells"]}"])
	if is_not_blank("${settings["cellranger_multi_chemistry"]}"):
		w.writerow(["chemistry","${settings["cellranger_multi_chemistry"]}"])
	if is_not_blank("${settings["r1_length"]}"):
		w.writerow(["r1-length","${settings["r1_length"]}"])
	if is_not_blank("${settings["include_introns"]}"):
		w.writerow(["include-introns","${settings["include_introns"]}"])
	if is_not_blank("${settings["check_library_compatibility"]}"):
		w.writerow(["check-library-compatibility","${settings["check_library_compatibility"]}"])
	if is_not_blank("${settings["create_bam"]}"):
		w.writerow(["create-bam","${settings["create_bam"]}"])
	if (is_not_blank("${feature_reference}") and ("${feature_reference}" != "NA")) or is_not_blank("${settings["r1_length_feature"]}"):
		w.writerow(["[feature]"])
		if is_not_blank("${feature_reference}") and ("${feature_reference}" != "NA"):
			feature_reference_path = os.path.abspath("${feature_reference}")
			w.writerow(["reference",feature_reference_path])
		if is_not_blank("${settings["r1_length_feature"]}"):
			w.writerow(["r1-length","${settings["r1_length_feature"]}"])
	if is_not_blank("${VDJ_reference}") and ("${VDJ_reference}" != "NA"):
		w.writerow(["[vdj]"])
		vdj_path = os.path.abspath("${VDJ_reference}")
		w.writerow(["reference",vdj_path])
		if is_not_blank("${settings["r1_length_vdj"]}"):
			w.writerow(["r1-length","${settings["r1_length_vdj"]}"])

	if fastq_id:
		w.writerow(["[libraries]"])
		w.writerow(["fastq_id","fastqs","physical_library_id","feature_types"])
		for i in range(len(fastq_id)):
			if (lanes10x[i] == "all" or lanes10x[i] == "") or (unique_lane_id == lanes10x[i]) :
				for b in range(len(bcl_directory)):
					bcl_name = bcl_directory[b].rstrip('/').rsplit('/', 1)[1]
					print(bcl_name)
					for root, dirs, files in os.walk(bcl_name + "_fastq"):
						if any(fname.startswith(fastq_id[i]) for fname in files):
							bcl_path = os.path.abspath(root)
							print("directory: "+bcl_path)
							print(glob.glob(os.path.join(root, fastq_id[i] + "*")))
							w.writerow([fastq_id[i],bcl_path,"",feature_types[i]])

	if is_nonempty_string_or_array(sample_id): 
		w.writerow(["[samples]"])
		w.writerow(["sample_id", nonempty_columns[0][0], "description"])
		for i in range(len(sample_id)):
			if (lane_of_samples[i] == "all" or lane_of_samples[i] == "") or (unique_lane_id == lane_of_samples[i]) :
				w.writerow([sample_id[i],nonempty_columns[0][1][i],description[i]])
				
print("## Config File:")
with open('final_${config}', 'r') as f:
	print(f.read())

find_process = subprocess.run(["find", "."], cwd=os.getcwd(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
print(find_process.stdout)

cmd = ["cellranger","multi","--id=${run_id}","--csv=final_${config}","--localcores=${cpu}","--localmem=${memory}"]
print(cmd)

p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = p.communicate()
print(stdout.decode())
print(stderr.decode())
if p.returncode != 0:
	sys.exit(p.returncode)


	
subprocess.run("mv ${run_id}/outs ${run_id}_outs && rm -rf ${run_id}", shell=True, check=True)
subprocess.run("for i in \$(ls ${run_id}_outs/per_sample_outs); do cp ${run_id}_outs/per_sample_outs/\${i}/web_summary.html ${config_id}_\${i}_web_summary.html; done", shell=True, check=True)
try:
	subprocess.run("for i in \$(ls ${run_id}_outs/per_sample_outs); do cp ${run_id}_outs/per_sample_outs/\${i}/vdj_b/filtered_contig_annotations.csv ${config_id}_\${i}_vdj_b_filtered_contig_annotations.csv; done", shell=True, check=True)
except subprocess.CalledProcessError:
	print("INFO: vdj_b/filtered_contig_annotations.csv was not found to publish it to report tab.")
try:
	subprocess.run("for i in \$(ls ${run_id}_outs/per_sample_outs); do cp ${run_id}_outs/per_sample_outs/\${i}/vdj_t/filtered_contig_annotations.csv ${config_id}_\${i}_vdj_t_filtered_contig_annotations.csv; done", shell=True, check=True)
except subprocess.CalledProcessError:
	print("INFO: vdj_t/filtered_contig_annotations.csv was not found to publish it to report tab.")
try:
	subprocess.run("for i in \$(ls ${run_id}_outs/per_sample_outs); do cp ${run_id}_outs/per_sample_outs/\${i}/vdj_t_gd/filtered_contig_annotations.csv ${config_id}_\${i}_vdj_t_gd_filtered_contig_annotations.csv; done", shell=True, check=True)
except subprocess.CalledProcessError:
	print("INFO: vdj_t_gd/filtered_contig_annotations.csv was not found to publish it to report tab.")


"""


}


process Multi_h5_explorer {

input:
 path output_dir

output:
 path "*.h5"  ,emit:g_59_h5_file00_g_61 

script:
"""
#!/bin/bash




outputdir=${output_dir}


if test -d \$outputdir/per_sample_outs/
		then
		for file in \$outputdir/per_sample_outs/* ;do 
			samplename=\${file/\$outputdir\\/per_sample_outs\\//""}
			echo \$samplename
			echo \$file
			cp \${file}/*_filtered_feature_bc_matrix.h5 ./\$samplename.h5
			if test -f \${file}/*filtered_feature_bc_matrix.h5 
				then
					echo "passed"
				else
					cp \${file}/count/*_filtered_feature_bc_matrix.h5 ./\$samplename.h5
			fi
		done
fi



"""

}


process file_to_set_conversion_for_h5 {

input:
 path h5

output:
 tuple val(name),file(h5)  ,emit:g_61_h5_file00_g51_0 

script:
name = h5.baseName
"""
echo "done"	
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 200
}
//* platform
//* platform
//* autofill

process scRNA_Analysis_Module_Quality_Control_and_Filtering {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_filtering_report.html$/) "QC_Reports/$filename"}
input:
 tuple val(name), file(input_file)
 path metadata

output:
 path "${name}.rds"  ,emit:g51_0_rdsFile00_g51_14 
 tuple val(name),file("${name}_filtering_report.html")  ,emit:g51_0_outputFileHTML11 
 path "${name}_filter_summary.tsv"  ,emit:g51_0_outFileTSV20_g51_34 

container "quay.io/viascientific/scrna_seurat:2.0"

when:
(params.run_scRNA_Analysis && (params.run_scRNA_Analysis == "yes")) || !params.run_scRNA_Analysis || params.run_pySCENIC == "yes"

script:

remove_mitochondiral_genes = params.scRNA_Analysis_Module_Quality_Control_and_Filtering.remove_mitochondiral_genes
remove_ribosomal_genes = params.scRNA_Analysis_Module_Quality_Control_and_Filtering.remove_ribosomal_genes
	
min_genes = params.scRNA_Analysis_Module_Quality_Control_and_Filtering.min_genes
max_genes = params.scRNA_Analysis_Module_Quality_Control_and_Filtering.max_genes
min_UMIs = params.scRNA_Analysis_Module_Quality_Control_and_Filtering.min_UMIs
max_UMIs = params.scRNA_Analysis_Module_Quality_Control_and_Filtering.max_UMIs
percent_mitochondrial = params.scRNA_Analysis_Module_Quality_Control_and_Filtering.percent_mitochondrial
percent_ribosomal = params.scRNA_Analysis_Module_Quality_Control_and_Filtering.percent_ribosomal

doublet_removal = params.scRNA_Analysis_Module_Quality_Control_and_Filtering.doublet_removal
doublet_percentage = params.scRNA_Analysis_Module_Quality_Control_and_Filtering.doublet_percentage

normalization_method = params.scRNA_Analysis_Module_Quality_Control_and_Filtering.normalization_method
variable_features = params.scRNA_Analysis_Module_Quality_Control_and_Filtering.variable_features

remove_mitochondiral_genes_arg = remove_mitochondiral_genes == 'true' ? '--remove-mitochondrial-genes' : ''
remove_ribosomal_genes_arg = remove_ribosomal_genes == 'true' ? '--remove_ribosomal_genes' : ''

//* @style @multicolumn:{remove_mitochondiral_genes, remove_ribosomal_genes}, {min_genes, max_genes}, {min_UMIs, max_UMIs}, {percent_mitochondrial, percent_ribosomal}, {doublet_removal, doublet_percentage}, {normalization_method, variable_features}

"""
build_QC_report.py --output-prefix ${name} --input-file ${input_file} --metadata-file ${metadata} \
--min-genes ${min_genes} --max-genes ${max_genes} \
--min-UMIs ${min_UMIs} --max-UMIs ${max_UMIs} \
--percent-mitochondrial-cutoff ${percent_mitochondrial} --percent-ribosomal-cutoff ${percent_ribosomal} ${remove_mitochondiral_genes_arg} ${remove_ribosomal_genes_arg} \
--variable-features ${variable_features} --normalization-method ${normalization_method} \
--doublet-removal ${doublet_removal} --doublet-percentage ${doublet_percentage}
"""



}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 50
}
//* platform
//* platform
//* autofill

process scRNA_Analysis_Module_Merge_Seurat_Objects {

input:
 path seurat_obj

output:
 path "merged_filtered_seurat.rds"  ,emit:g51_14_rdsFile00_g51_17 

container "quay.io/viascientific/scrna_seurat:2.0"

script:
"""
#!/usr/bin/env Rscript

library(Seurat)
list_of_samples <- list.files(pattern = "*.rds")

if (length(list_of_samples)==1) {
	list_of_seurat=list()
	obj = readRDS(list_of_samples[1])
	list_of_seurat[[list_of_samples[1]]]=obj
} else {
list_of_seurat <- list()
for(i in 1:length(list_of_samples)){
  # print name
  print(list_of_samples[i])
  list_of_seurat[[i]] <- readRDS(list_of_samples[i])
}

}
saveRDS(list_of_seurat, file="merged_filtered_seurat.rds")

"""


}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 140
}
//* platform
//* platform
//* autofill

process scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction {

input:
 path seurat_object

output:
 path "Reduced_and_Corrected.rds"  ,emit:g51_17_rdsFile00_g51_19 

container "quay.io/viascientific/scrna_seurat:2.0"

script:
varFeatures = params.scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction.varFeatures
selmethod = params.scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction.selmethod
Batch_Effect_Correction = params.scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction.Batch_Effect_Correction
WNN = params.scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction.WNN

//* @style @multicolumn:{varFeatures, selmethod},{Batch_Effect_Correction, WNN}

"""
#!/usr/bin/env Rscript

# libraries
library(Seurat)
library(dplyr)
#install.packages("harmony",repos = "http://cran.us.r-project.org")
library(harmony)

selmethod <- "${selmethod}"
varFeatures <- "${varFeatures}"

Data=readRDS("${seurat_object}")
Multi_sample=0
if (length(Data)==1) {
	Data=Data[[1]]
	if (DefaultAssay(Data)=="SCT"){
		Data=RunPCA(Data,npcs=100)
	} else {
		Data <- FindVariableFeatures(Data,selection.method=selmethod,nfeatures=as.numeric(varFeatures))
		if (all(Data[["percent.mt"]]==0)) {
		Data <- ScaleData(Data)
		} else {
		Data <- ScaleData(Data,vars.to.regress="percent.mt")
		}
		Data=RunPCA(Data,npcs=100)
	}
} else {
Multi_sample=1
	if (DefaultAssay(Data[[1]])=="SCT") {
		variable.features=SelectIntegrationFeatures(object.list = Data, nfeatures = as.numeric(varFeatures))
		Data <- merge(Data[[1]],Data[-1])
		VariableFeatures(Data) <- variable.features

		if (all(Data[["percent.mt"]]==0)) {
			Data <- ScaleData(Data)
		} else {
			Data <- ScaleData(Data,vars.to.regress="percent.mt")
			}
		Data=RunPCA(Data,npcs=100)
		if (as.logical("${Batch_Effect_Correction}")){
		Data=RunHarmony(Data,assay.use = DefaultAssay(Data),group.by.vars = "sample",max.iter.harmony = 10000,max.iter.cluster = 10000)
		}
	} else {
		Data <- merge(Data[[1]],Data[-1])
		Data <- FindVariableFeatures(Data,selection.method=selmethod,nfeatures=as.numeric(varFeatures))
		if (all(Data[["percent.mt"]]==0)) {
			Data <- ScaleData(Data)
		} else {
			Data <- ScaleData(Data,vars.to.regress="percent.mt")
		}
		Data=RunPCA(Data,npcs=100)
		if (as.logical("${Batch_Effect_Correction}")){
		Data=RunHarmony(Data,assay.use = DefaultAssay(Data),group.by.vars = "sample",max.iter.harmony = 10000,max.iter.cluster = 10000)
		}

	}
}

if ("${WNN}"!="") {
original.assay=DefaultAssay(Data)

DefaultAssay(Data)="${WNN}"

Data=NormalizeData(Data,normalization.method = "CLR",margin=2)

VariableFeatures(Data)=rownames(Data)

Data=ScaleData(Data)

Data=RunPCA(Data,reduction.name = "wpca")

if (Multi_sample==1) {
	Data=RunHarmony(Data,group.by.vars = "sample",assay.use = "${WNN}",reduction = "wpca",reduction.save = "wharmony")

}

DefaultAssay(Data)=original.assay

}

saveRDS(Data,"Reduced_and_Corrected.rds")

"""


}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 16
    $MEMORY = 140
}
//* platform
//* platform
//* autofill

process scRNA_Analysis_Module_Clustering_and_Find_Markers {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /final_report.html$/) "Final_Report/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tsv$/) "ClusterMarkers/$filename"}
input:
 path seurat_object

output:
 path "final_report.html"  ,emit:g51_19_outputHTML00 
 path "Final_Analysis.rds"  ,emit:g51_19_rdsFile10_g51_36 
 path "*.tsv"  ,emit:g51_19_outFileTSV22 

container "quay.io/viascientific/scrna_seurat:2.0"

script:
min_resolution = params.scRNA_Analysis_Module_Clustering_and_Find_Markers.min_resolution
max_resolution = params.scRNA_Analysis_Module_Clustering_and_Find_Markers.max_resolution
num_pc = params.scRNA_Analysis_Module_Clustering_and_Find_Markers.num_pc
find_markers_for_all_resolution = params.scRNA_Analysis_Module_Clustering_and_Find_Markers.find_markers_for_all_resolution

find_markers_for_all_resolution_arg = find_markers_for_all_resolution == 'true' ? '--all-resolution-cluster-markers' : ''

//* @style @multicolumn:{min_resolution, max_resolution}, {num_pc, find_markers_for_all_resolution}

"""
build_clustering_and_find_markers.py --sample-path ${seurat_object} --min-resolution ${min_resolution} --max-resolution ${max_resolution} --num-pc ${num_pc} ${find_markers_for_all_resolution_arg}
"""

}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 4
}
//* platform
//* platform
//* autofill

process scRNA_Analysis_Module_filter_summary {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*$/) "filter_summary/$filename"}
input:
 path input_files

output:
 path "*"  ,emit:g51_34_outputFileHTML00 

container "quay.io/viascientific/scrna_seurat:2.0"

script:
	
"""
build_filtration_report.py --input-dir .

mkdir output
mv by_criteria_summary.tsv output
mv filtration_summary_report.Rmd output
mv overall_filtration_summary.tsv output
"""
}


process scRNA_Analysis_Module_sc_annotation {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.rds$/) "scViewer/$filename"}
input:
 path input_rds

output:
 path "*.rds"  ,emit:g51_36_rdsFile00_g51_22 

container 'quay.io/mustafapir/sc_annotation:1.0.0'

script:

tissue_type = params.scRNA_Analysis_Module_sc_annotation.tissue_type

"""
if [ ${params.run_annotation} = "yes" ]; then
  run_sctype.R \
    --input ${input_rds} \
    --output Final_Analysis_annotated.rds \
    --organism ${params.species} \
    --tissue "${tissue_type}" && \
  rm -f ${input_rds}
fi

"""

}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 50
}
//* platform
//* platform
//* autofill

process scRNA_Analysis_Module_SCEtoLOOM {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /Data.loom$/) "LOOM/$filename"}
input:
 path seurat_object

output:
 path "Data.loom" ,optional:true  ,emit:g51_30_outputFileOut00_g82_1 

container "quay.io/mustafapir/scrna_seurat:2.0.2"

script:
Generate_loom_file = params.scRNA_Analysis_Module_SCEtoLOOM.Generate_loom_file

"""
#!/usr/bin/env Rscript

#library
library(Seurat)
library(biomaRt)
 
if (as.logical("${Generate_loom_file}") || as.logical("${params.run_pySCENIC == 'yes'}")) {

Data=readRDS("${seurat_object}")

create_annotation_table <- function(organism = c("human", "mouse", "d_melanogaster")) {
  organism <- match.arg(organism)
  
  # Select ENSEMBL dataset depending on organism
  dataset <- switch(
    organism,
    human = "hsapiens_gene_ensembl",
    mouse = "mmusculus_gene_ensembl",
    d_melanogaster = "dmelanogaster_gene_ensembl"
  )
  
  # Connect to Ensembl
  ensembl <- useEnsembl(biomart = "genes", dataset = dataset)
  
  # Retrieve Ensembl gene ID, gene name, gene type
  genes <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
    mart = ensembl
  )
  
  # Ensure column names match your exact requested format
  colnames(genes) <- c("ensembl_id", "gene_name", "gene_type")
  
  return(genes)
}

# annotation=read.csv("https://huggingface.co/datasets/ctheodoris/Genecorpus-30M/raw/main/example_input_files/gene_info_table.csv",header = T,row.names = 1)
annotation=create_annotation_table("${params.species}")

annotation=annotation[annotation\$gene_name%in%names(table(annotation\$gene_name))[table(annotation\$gene_name)==1],]
rownames(annotation)=annotation\$gene_name
metadata=Data@meta.data
matrix=Data@assays\$RNA@counts
matrix=matrix[rowSums(matrix)>0,]
matrix=matrix[intersect(rownames(matrix),rownames(annotation)),]
annotation=annotation[intersect(rownames(matrix),rownames(annotation)),]
matrix=matrix[rownames(matrix)[order(rownames(matrix),decreasing = F)],]
NewData=CreateSeuratObject(matrix,meta.data = metadata)
NewData[["RNA"]]@meta.features\$ensembl_id=annotation[rownames(NewData),"ensembl_id"]
NewData.loom <- SeuratDisk::as.loom(NewData, filename = "Data.loom", verbose = FALSE)

}
"""


}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 60
}
//* platform
//* platform
//* autofill

process scRNA_Analysis_Module_Create_h5ad {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.h5ad$/) "H5AD_file/$filename"}
input:
 path seurat_obj

output:
 path "*.h5ad"  ,emit:g51_22_h5ad_file01_g70_12 

container "quay.io/viascientific/scrna_seurat:2.0"

script:
"""
#!/usr/bin/env Rscript

# libraries
library(Seurat)
library(SeuratDisk)

# read data
seurat_obj <- readRDS("${seurat_obj}")

for (var in colnames(seurat_obj@meta.data)) {
  if (is.factor(seurat_obj@meta.data[,var])) {
    seurat_obj@meta.data[,var]=as.character(seurat_obj@meta.data[,var])
  }
}

# save h5ad file
seu_name <- gsub(".rds","","${seurat_obj}")
SaveH5Seurat(seurat_obj, filename = paste0(seu_name,".h5Seurat"))
Convert(paste0(seu_name,".h5Seurat"), dest = "h5ad")
"""


}


process RNA_Velocity_Module_prepare_input_velocyto {

input:
 path outs

output:
 path "output_files/*"  ,emit:g70_5_outputDir00_g70_1 

when:
params.run_velocity == "yes"

script:

try {
	myVariable = bam
} catch (MissingPropertyException e) {
	bam = ""
}

try {
	myVariable = bai
} catch (MissingPropertyException e) {
	bai = ""
}

try {
	myVariable = barcodes
} catch (MissingPropertyException e) {
	barcodes = ""
}

run_name = outs.toString().startsWith('NO_FILE')
    ? name
    : outs.toString().replaceFirst(/_outs$/, '')

"""
mkdir -p output_files

echo ${run_name}

if [[ ${outs} == NO_FILE* ]]; then
    # Create a folder using the task name so the structure matches the loop below
    mkdir -p output_files/${run_name}
    
    mv ${bam} output_files/${run_name}/input.bam
    mv ${bai} output_files/${run_name}/input.bam.bai
    mv ${barcodes} output_files/${run_name}/input_barcodes.tsv.gz


# --- SCENARIO B: Cell Ranger Output Directory ---
else
    # Check for Cell Ranger MULTI (Branching structure)
    if [ -d ${outs}/per_sample_outs ]; then
        echo "Detected Cell Ranger MULTI structure"
        
        # Loop through every directory inside per_sample_outs
        for sample_dir in ${outs}/per_sample_outs/*; do
            
            # Extract the sample name (e.g., "Sample_Alpha")
            s_name=\$(basename \$sample_dir)
            
            # Create a specific folder for this sample
            target_dir="output_files/\$s_name"
            mkdir -p "\$target_dir"
            
            # Move and Rename BAM
            # Note: In multi, it is named 'sample_alignments.bam'
            if [ -f "\$sample_dir/count/sample_alignments.bam" ]; then
                mv "\$sample_dir/count/sample_alignments.bam" "\$target_dir/input.bam"
                mv "\$sample_dir/count/sample_alignments.bam.bai" "\$target_dir/input.bam.bai"
                
                # Move Barcodes (Ensure we get the filtered matrix, not raw)
                mv "\$sample_dir/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz" "\$target_dir/input_barcodes.tsv.gz"
            else
                echo "Warning: No count data found for \$s_name in \$sample_dir"
                rm -rf "\$target_dir" # Clean up empty dir
            fi
        done

    # Check for Cell Ranger COUNT (Flat structure)
    elif [ -f ${outs}/possorted_genome_bam.bam ]; then
        echo "Detected Cell Ranger COUNT structure"
        
        # Use the run_name derived from the folder name
        mkdir -p output_files/${run_name}
        
        mv "${outs}/possorted_genome_bam.bam" output_files/${run_name}/input.bam
        mv "${outs}/possorted_genome_bam.bam.bai" output_files/${run_name}/input.bam.bai
        mv "${outs}/filtered_feature_bc_matrix/barcodes.tsv.gz" output_files/${run_name}/input_barcodes.tsv.gz

    else
        echo "ERROR: Unknown Cell Ranger structure in ${outs}"
        exit 1
    fi
fi
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 16
    $MEMORY = 50
}
//* platform
//* platform
//* autofill

process RNA_Velocity_Module_velocyto {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_output.loom$/) "loom_out/$filename"}
input:
 tuple val(name), file(bam), file(barcodes)
 path mask_gtf
 path gtf_velo

output:
 path "${name}_output.loom"  ,emit:g70_1_loom00_g70_12 

container "quay.io/biocontainers/velocyto.py:0.17.17--py310h581d4b6_7"

when:
params.run_velocity == "yes"

script:

mask_gtf_option = mask_gtf.name.startsWith('NO_FILE') ? "" : "-m ${mask_gtf}"

"""
echo ${name}
mkdir -p velocyto_out

cp ${bam} ./local_input.bam
velocyto run -b ${barcodes} ${mask_gtf_option} -o velocyto_out local_input.bam ${gtf_velo}

mv velocyto_out/*.loom ${name}_output.loom

"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 8
    $MEMORY = 50
}
//* platform
//* platform
//* autofill

process RNA_Velocity_Module_process_anndata {

input:
 path loom_file
 path h5ad_file

output:
 path "processed_adata.h5ad"  ,emit:g70_12_h5ad00_g70_15 

container 'quay.io/viascientific/scvelo_shiny:1.1.4'

when:
params.run_velocity == "yes"

script:

"""

preprocess_anndata.py \
    --h5ad ${h5ad_file} \
    --loom ${loom_file} \
    --output 'processed_adata.h5ad'

"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 16
    $MEMORY = 64
}
//* platform
//* platform
//* autofill

process RNA_Velocity_Module_process_scVelo {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /scvelo_out.h5ad$/) "scVelo_out/$filename"}
input:
 path input_adata

output:
 path "scvelo_out.h5ad"  ,emit:g70_15_h5ad00 

container "quay.io/viascientific/scvelo_shiny:1.1.4"

script:

// clusters = "seurat_clusters" //* @input @description:"Name of clusters column"
clusters = params.run_annotation == 'yes' ? 'sctype_classification' : 'seurat_clusters'

conditions = params.RNA_Velocity_Module_process_scVelo.conditions
k_steps = params.RNA_Velocity_Module_process_scVelo.k_steps
n_macrostates = params.RNA_Velocity_Module_process_scVelo.n_macrostates
target_clusters = params.RNA_Velocity_Module_process_scVelo.target_clusters

target_clusters_arg = (target_clusters == '') ? '' : "--target_clusters ${target_clusters}"

"""


precompute_analysis.py ${input_adata} \
    --group_col ${clusters} \
    --condition_col ${conditions} \
    --k_steps ${k_steps} \
    --n_macrostates ${n_macrostates} \
    ${target_clusters_arg} \
    --output_file scvelo_out.h5ad \
    --n_jobs 1

"""

}


process Trajectory_Module_slingshot {

input:
 path input_rds

output:
 path "sce.rds"  ,emit:g80_1_rdsFile01_g80_4 

container 'quay.io/viascientific/slingshot:1.0.1'

when:
params.run_slingshot == "yes"

script:

reduction = params.Trajectory_Module_slingshot.reduction

clusters = params.run_annotation == 'yes' ? 'sctype_classification' : 'seurat_clusters'
"""

slingshot_analysis.R ${input_rds} ${reduction} ${clusters}

"""

}


//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 16
    $MEMORY = 128
}
//* platform
//* platform
//* autofill

process Trajectory_Module_fitgam {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /(sce_fitgam.rds|sce.rds)$/) "trajectory_out/$filename"}
input:
 path obj
 path sce

output:
 tuple file("sce_fitgam.rds"), file("sce.rds")  ,emit:g80_4_rdsFile00 

container 'quay.io/viascientific/slingshot:1.0.1'

script:

threads = task.cpus
num_genes = params.Trajectory_Module_fitgam.num_genes
"""

fitgam.R ${obj} ${sce} ${threads} ${num_genes}
"""

}

//* params.db_feather =  ""  //* @input
//* params.motif_db =  ""  //* @input
//* params.tf_lists =  ""  //* @input


//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 16
    $MEMORY = 64
}
//* platform
//* platform
//* autofill

process pySCENIC_module_pySCENIC_GRN {

input:
 path loom_file
 path tf_lists

output:
 path "adjacencies.csv"  ,emit:g82_1_csvFile00_g82_8 

if (params.GRN_method == "regdiffusion_GPU"){
    container 'quay.io/viascientific/pyscenic:1.0.2_gpu'

    machineType 'g2-standard-(4|8|16)'
    containerOptions '-v /var/lib/nvidia/lib64:/usr/local/nvidia/lib64 --device /dev/nvidia0:/dev/nvidia0 --device /dev/nvidia-uvm:/dev/nvidia-uvm --device /dev/nvidiactl:/dev/nvidiactl -e LD_LIBRARY_PATH=/usr/local/nvidia/lib64:/usr/local/cuda/lib64:${LD_LIBRARY_PATH}'
} else if (params.GRN_method == "regdiffusion_CPU"){
    container 'quay.io/viascientific/pyscenic:1.0.2_gpu'
} else {
    container 'quay.io/viascientific/pyscenic:1.0.2'
}

when:
params.run_pySCENIC == "yes"

script:

threads = task.cpus
GRN_method = params.GRN_method
// GRN_method = "GRNBoost2" //* @dropdown @options:"GRNBoost2","regdiffusion_GPU","regdiffusion_CPU" @label:"Method to use in GRN inference" @description:"Algorithm to use for gene regulatory network inference. GRNBoost2 is the original method used in pySCENIC workflow but it is the slowest one. Using `regdiffusiton_GPU` will speed up the analysis significantly but will cost more."
num_of_hvg = params.pySCENIC_module_pySCENIC_GRN.num_of_hvg

arg_num_of_hvg = (num_of_hvg == "") ? "" : "--num_hvg ${num_of_hvg}"

"""

if [ ${GRN_method} = "regdiffusion_GPU" ]; then
    run_regdiffusion_gpu.py \
        ${loom_file} ${tf_lists} \
        --output adjacencies.csv \
        --num_workers ${threads} \
        ${arg_num_of_hvg}
elif [ ${GRN_method} = "regdiffusion_CPU" ]; then
    run_regdiffusion_gpu.py \
        ${loom_file} ${tf_lists} \
        --output adjacencies.csv \
        --num_workers ${threads} \
        ${arg_num_of_hvg} \
        --cpu
else
    arboreto_with_multiprocessing.py \
        ${loom_file} ${tf_lists} \
        --output adjacencies.csv \
        --method grnboost2 \
        --seed 29 \
        --num_workers ${threads} \
        --sparse
fi

"""

}


//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 16
    $MEMORY = 64
}
//* platform
//* platform
//* autofill

process pySCENIC_module_pySCENIC_ctx_auc {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /pyscenic_out.zip$/) "pySCENIC_out/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /scenic_integrated.loom$/) "pySCENIC_loom/$filename"}
input:
 path adj
 path db_feather
 path motif_db
 path loom_file

output:
 path "pyscenic_out.zip"  ,emit:g82_8_zip_file00 
 path "scenic_integrated.loom"  ,emit:g82_8_loom11 

container 'quay.io/viascientific/pyscenic:1.0.3'

when:
params.run_pySCENIC == "yes"

script:

threads = task.cpus
mask_dropouts = params.pySCENIC_module_pySCENIC_ctx_auc.mask_dropouts

mask_dropouts_option = mask_dropouts ? '--mask_dropouts' : ''
auc_threshold = params.pySCENIC_module_pySCENIC_ctx_auc.auc_threshold
"""

f_db_path=${db_feather}
f_db_names=\$(echo "\$f_db_path"/*.feather)

f_motif_path=${motif_db}

pyscenic ctx ${adj} \
    \$f_db_names \
    --annotations_fname ${motif_db} \
    --expression_mtx_fname ${loom_file} \
    --output regulons.csv \
    ${mask_dropouts_option} \
    --num_workers ${threads}

pyscenic aucell \
    ${loom_file} \
    regulons.csv \
    --output pyscenic_out.loom \
    --num_workers ${threads} \
    --auc_threshold ${auc_threshold}

integrate_pyscenic_output.py \
    -i ${loom_file} \
    -p pyscenic_out.loom \
    -o scenic_integrated.loom \
    --export_auc_csv aucell_matrix.csv \
    -t ${threads}

zip pyscenic_out.zip adjacencies.csv regulons.csv aucell_matrix.csv

"""

}


workflow {


Demultiplexer_prep()
g_4_bcl00_g_0 = Demultiplexer_prep.out.g_4_bcl00_g_0.flatten().collate(2).map({item ->  if (item[1]) { return tuple(file(item[0]),file(item[1])) } else { return  tuple(file(item[0]))   } })
g_4_bcl_directory16_g_20 = Demultiplexer_prep.out.g_4_bcl_directory16_g_20


bclConvert(g_4_bcl00_g_0)
g_0_reads00_g_33 = bclConvert.out.g_0_reads00_g_33
(g_0_reads00_g_20) = [g_0_reads00_g_33]
g_0_outputDir11 = bclConvert.out.g_0_outputDir11


cellranger_fastq_prep(g_23_0_g_22,g_24_1_g_22)
g_22_reads00_g_25 = cellranger_fastq_prep.out.g_22_reads00_g_25


cellranger_fastq_collect(g_22_reads00_g_25.collect())
g_25_bcl_directory05_g_20 = cellranger_fastq_collect.out.g_25_bcl_directory05_g_20
g_25_reads17_g_20 = cellranger_fastq_collect.out.g_25_reads17_g_20


if (!((params.run_FastQC && (params.run_FastQC == "yes")))){
g_0_reads00_g_33.set{g_33_reads10_g_34}
g_33_mate00_g_35 = Channel.empty()
} else {

flatten_cellranger_reads(g_0_reads00_g_33.collect())
g_33_mate00_g_35 = flatten_cellranger_reads.out.g_33_mate00_g_35
g_33_reads10_g_34 = flatten_cellranger_reads.out.g_33_reads10_g_34
}


file_to_set_conversion_for_reads(g_33_reads10_g_34.flatten())
g_34_reads01_g_35 = file_to_set_conversion_for_reads.out.g_34_reads01_g_35


if (!((params.run_FastQC && (params.run_FastQC == "yes")))){
g_34_reads01_g_35.set{g_35_reads11}
g_35_FastQCout04_g_52 = Channel.empty()
} else {

FastQC_after_mkfastq(g_33_mate00_g_35,g_34_reads01_g_35)
g_35_FastQCout04_g_52 = FastQC_after_mkfastq.out.g_35_FastQCout04_g_52
g_35_reads11 = FastQC_after_mkfastq.out.g_35_reads11
}


MultiQC(g_35_FastQCout04_g_52.flatten().toList())
g_52_outputHTML00 = MultiQC.out.g_52_outputHTML00
g_52_outputDir11 = MultiQC.out.g_52_outputDir11


cellranger_multi_library_prep()
g_67_librarySettings00_g_5 = cellranger_multi_library_prep.out.g_67_librarySettings00_g_5


cellranger_multi_prep(g_67_librarySettings00_g_5)
g_5_csvFile04_g_20 = cellranger_multi_prep.out.g_5_csvFile04_g_20
g_5_settings18_g_20 = cellranger_multi_prep.out.g_5_settings18_g_20


Check_and_Build_Module_Check_Genome_GTF()
g50_21_genome00_g50_58 = Check_and_Build_Module_Check_Genome_GTF.out.g50_21_genome00_g50_58
g50_21_gtfFile10_g50_57 = Check_and_Build_Module_Check_Genome_GTF.out.g50_21_gtfFile10_g50_57


if (!(params.replace_geneID_with_geneName == "yes")){
g50_21_gtfFile10_g50_57.set{g50_57_gtfFile01_g50_58}
} else {

Check_and_Build_Module_convert_gtf_attributes(g50_21_gtfFile10_g50_57)
g50_57_gtfFile01_g50_58 = Check_and_Build_Module_convert_gtf_attributes.out.g50_57_gtfFile01_g50_58
}



if (!(params.add_sequences_to_reference == "yes")){
g50_21_genome00_g50_58.set{g50_58_genome00_g50_52}
(g50_58_genome01_g50_54) = [g50_58_genome00_g50_52]
g50_57_gtfFile01_g50_58.set{g50_58_gtfFile10_g50_53}
(g50_58_gtfFile10_g50_54) = [g50_58_gtfFile10_g50_53]
} else {

Check_and_Build_Module_Add_custom_seq_to_genome_gtf(g50_21_genome00_g50_58,g50_57_gtfFile01_g50_58,g_40_2_g50_58,g_69_3_g50_58)
g50_58_genome00_g50_52 = Check_and_Build_Module_Add_custom_seq_to_genome_gtf.out.g50_58_genome00_g50_52
(g50_58_genome01_g50_54) = [g50_58_genome00_g50_52]
g50_58_gtfFile10_g50_53 = Check_and_Build_Module_Add_custom_seq_to_genome_gtf.out.g50_58_gtfFile10_g50_53
(g50_58_gtfFile10_g50_54) = [g50_58_gtfFile10_g50_53]
}


Check_and_Build_Module_Check_BED12(g50_58_gtfFile10_g50_53)
g50_53_bed03_g50_54 = Check_and_Build_Module_Check_BED12.out.g50_53_bed03_g50_54


Check_and_Build_Module_Check_chrom_sizes_and_index(g50_58_genome00_g50_52)
g50_52_genomeSizes02_g50_54 = Check_and_Build_Module_Check_chrom_sizes_and_index.out.g50_52_genomeSizes02_g50_54

g50_58_gtfFile10_g50_54= g50_58_gtfFile10_g50_54.ifEmpty(ch_empty_file_1) 
g50_58_genome01_g50_54= g50_58_genome01_g50_54.ifEmpty(ch_empty_file_2) 
g50_52_genomeSizes02_g50_54= g50_52_genomeSizes02_g50_54.ifEmpty(ch_empty_file_3) 
g50_53_bed03_g50_54= g50_53_bed03_g50_54.ifEmpty(ch_empty_file_4) 


Check_and_Build_Module_check_files(g50_58_gtfFile10_g50_54,g50_58_genome01_g50_54,g50_52_genomeSizes02_g50_54,g50_53_bed03_g50_54)
g50_54_gtfFile01_g_7 = Check_and_Build_Module_check_files.out.g50_54_gtfFile01_g_7
(g50_54_gtfFile02_g70_1) = [g50_54_gtfFile01_g_7]
g50_54_genome10_g_7 = Check_and_Build_Module_check_files.out.g50_54_genome10_g_7
g50_54_genomeSizes22 = Check_and_Build_Module_check_files.out.g50_54_genomeSizes22
g50_54_bed33 = Check_and_Build_Module_check_files.out.g50_54_bed33


cellranger_mkref(g50_54_genome10_g_7,g50_54_gtfFile01_g_7)
g_7_reference00_g_18 = cellranger_mkref.out.g_7_reference00_g_18

g_7_reference00_g_18= g_7_reference00_g_18.ifEmpty(ch_empty_file_1) 


cellranger_ref_checker(g_7_reference00_g_18)
g_18_reference01_g_20 = cellranger_ref_checker.out.g_18_reference01_g_20

g_0_reads00_g_20= g_0_reads00_g_20.ifEmpty(ch_empty_file_1) 
g_25_bcl_directory05_g_20= g_25_bcl_directory05_g_20.ifEmpty("") 
g_4_bcl_directory16_g_20= g_4_bcl_directory16_g_20.ifEmpty("") 
g_25_reads17_g_20= g_25_reads17_g_20.ifEmpty(ch_empty_file_6) 


if (!((params.run_cellranger_multi && (params.run_cellranger_multi == "yes")) || !params.run_cellranger_multi)){
g_5_csvFile04_g_20.set{g_20_csvFile22}
g_20_outputDir00_g_59 = Channel.empty()
g_20_outputDir00_g70_5 = Channel.empty()
g_20_outputHTML11 = Channel.empty()
} else {

cellranger_multi(g_0_reads00_g_20.collect(),g_18_reference01_g_20,g_11_2_g_20,g_12_3_g_20,g_5_csvFile04_g_20.flatten(),g_25_bcl_directory05_g_20,g_4_bcl_directory16_g_20,g_25_reads17_g_20.collect(),g_5_settings18_g_20,g_68_9_g_20)
g_20_outputDir00_g_59 = cellranger_multi.out.g_20_outputDir00_g_59
(g_20_outputDir00_g70_5) = [g_20_outputDir00_g_59]
g_20_outputHTML11 = cellranger_multi.out.g_20_outputHTML11
g_20_csvFile22 = cellranger_multi.out.g_20_csvFile22
g_20_csvFile33 = cellranger_multi.out.g_20_csvFile33
}


Multi_h5_explorer(g_20_outputDir00_g_59)
g_59_h5_file00_g_61 = Multi_h5_explorer.out.g_59_h5_file00_g_61


file_to_set_conversion_for_h5(g_59_h5_file00_g_61.flatten())
g_61_h5_file00_g51_0 = file_to_set_conversion_for_h5.out.g_61_h5_file00_g51_0



scRNA_Analysis_Module_Quality_Control_and_Filtering(g_61_h5_file00_g51_0,g_43_1_g51_0)
g51_0_rdsFile00_g51_14 = scRNA_Analysis_Module_Quality_Control_and_Filtering.out.g51_0_rdsFile00_g51_14
g51_0_outputFileHTML11 = scRNA_Analysis_Module_Quality_Control_and_Filtering.out.g51_0_outputFileHTML11
g51_0_outFileTSV20_g51_34 = scRNA_Analysis_Module_Quality_Control_and_Filtering.out.g51_0_outFileTSV20_g51_34


scRNA_Analysis_Module_Merge_Seurat_Objects(g51_0_rdsFile00_g51_14.collect())
g51_14_rdsFile00_g51_17 = scRNA_Analysis_Module_Merge_Seurat_Objects.out.g51_14_rdsFile00_g51_17


scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction(g51_14_rdsFile00_g51_17)
g51_17_rdsFile00_g51_19 = scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction.out.g51_17_rdsFile00_g51_19


scRNA_Analysis_Module_Clustering_and_Find_Markers(g51_17_rdsFile00_g51_19)
g51_19_outputHTML00 = scRNA_Analysis_Module_Clustering_and_Find_Markers.out.g51_19_outputHTML00
g51_19_rdsFile10_g51_36 = scRNA_Analysis_Module_Clustering_and_Find_Markers.out.g51_19_rdsFile10_g51_36
g51_19_outFileTSV22 = scRNA_Analysis_Module_Clustering_and_Find_Markers.out.g51_19_outFileTSV22


scRNA_Analysis_Module_filter_summary(g51_0_outFileTSV20_g51_34.collect())
g51_34_outputFileHTML00 = scRNA_Analysis_Module_filter_summary.out.g51_34_outputFileHTML00


scRNA_Analysis_Module_sc_annotation(g51_19_rdsFile10_g51_36)
g51_36_rdsFile00_g51_22 = scRNA_Analysis_Module_sc_annotation.out.g51_36_rdsFile00_g51_22
(g51_36_rdsFile00_g51_30,g51_36_rdsFile00_g80_1,g51_36_rdsFile00_g80_4) = [g51_36_rdsFile00_g51_22,g51_36_rdsFile00_g51_22,g51_36_rdsFile00_g51_22]


scRNA_Analysis_Module_SCEtoLOOM(g51_36_rdsFile00_g51_30)
g51_30_outputFileOut00_g82_1 = scRNA_Analysis_Module_SCEtoLOOM.out.g51_30_outputFileOut00_g82_1
(g51_30_outputFileOut03_g82_8) = [g51_30_outputFileOut00_g82_1]


scRNA_Analysis_Module_Create_h5ad(g51_36_rdsFile00_g51_22)
g51_22_h5ad_file01_g70_12 = scRNA_Analysis_Module_Create_h5ad.out.g51_22_h5ad_file01_g70_12

g_20_outputDir00_g70_5= g_20_outputDir00_g70_5.ifEmpty(ch_empty_file_1) 


RNA_Velocity_Module_prepare_input_velocyto(g_20_outputDir00_g70_5)
g70_5_outputDir00_g70_1 = RNA_Velocity_Module_prepare_input_velocyto.out.g70_5_outputDir00_g70_1.flatten().map { item -> return tuple( item.name, item / "input.bam", item / "input_barcodes.tsv.gz" ) }



RNA_Velocity_Module_velocyto(g70_5_outputDir00_g70_1,g_71_1_g70_1,g50_54_gtfFile02_g70_1.first())
g70_1_loom00_g70_12 = RNA_Velocity_Module_velocyto.out.g70_1_loom00_g70_12


RNA_Velocity_Module_process_anndata(g70_1_loom00_g70_12.collect(),g51_22_h5ad_file01_g70_12)
g70_12_h5ad00_g70_15 = RNA_Velocity_Module_process_anndata.out.g70_12_h5ad00_g70_15


RNA_Velocity_Module_process_scVelo(g70_12_h5ad00_g70_15)
g70_15_h5ad00 = RNA_Velocity_Module_process_scVelo.out.g70_15_h5ad00


Trajectory_Module_slingshot(g51_36_rdsFile00_g80_1)
g80_1_rdsFile01_g80_4 = Trajectory_Module_slingshot.out.g80_1_rdsFile01_g80_4


Trajectory_Module_fitgam(g51_36_rdsFile00_g80_4,g80_1_rdsFile01_g80_4)
g80_4_rdsFile00 = Trajectory_Module_fitgam.out.g80_4_rdsFile00


pySCENIC_module_pySCENIC_GRN(g51_30_outputFileOut00_g82_1,g_77_1_g82_1)
g82_1_csvFile00_g82_8 = pySCENIC_module_pySCENIC_GRN.out.g82_1_csvFile00_g82_8


pySCENIC_module_pySCENIC_ctx_auc(g82_1_csvFile00_g82_8,g_75_1_g82_8,g_76_2_g82_8,g51_30_outputFileOut03_g82_8)
g82_8_zip_file00 = pySCENIC_module_pySCENIC_ctx_auc.out.g82_8_zip_file00
g82_8_loom11 = pySCENIC_module_pySCENIC_ctx_auc.out.g82_8_loom11


}

workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
