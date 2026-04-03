$HOSTNAME = ""
params.outdir = 'results'  

params.publishdict = [:]

def pathChecker(input, path, type){
	def cmd = "mkdir -p check && mv ${input} check/. "
	if (!input || input.empty()){
		input = file(path).getName().toString()
		cmd = "mkdir -p check && cd check && ln -s ${path} ${input} && cd .."
		if (path.indexOf('s3:') > -1 || path.indexOf('S3:') >-1){
			def recursive = (type == "folder") ? "--recursive" : ""
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
if (!params.cmo_set){params.cmo_set = ""} 
if (!params.feature_reference){params.feature_reference = ""} 
if (!params.VDJ_reference){params.VDJ_reference = ""} 
if (!params.feature_ref){params.feature_ref = ""} 
if (!params.reads){params.reads = ""} 
if (!params.mate){params.mate = ""} 
if (!params.custom_additional_genome){params.custom_additional_genome = ""} 
if (!params.Metadata){params.Metadata = ""} 
// Stage empty file to be used as an optional input where required
ch_empty_file_1 = file("$baseDir/.emptyfiles/NO_FILE_1", hidden:true)
ch_empty_file_2 = file("$baseDir/.emptyfiles/NO_FILE_2", hidden:true)
ch_empty_file_3 = file("$baseDir/.emptyfiles/NO_FILE_3", hidden:true)
ch_empty_file_4 = file("$baseDir/.emptyfiles/NO_FILE_4", hidden:true)
ch_empty_file_5 = file("$baseDir/.emptyfiles/NO_FILE_5", hidden:true)
ch_empty_file_6 = file("$baseDir/.emptyfiles/NO_FILE_6", hidden:true)
ch_empty_file_7 = file("$baseDir/.emptyfiles/NO_FILE_7", hidden:true)

g_10_4_g_55 = params.cmo_set && file(params.cmo_set, type: 'any').exists() ? file(params.cmo_set, type: 'any') : ch_empty_file_4
g_10_4_g_20 = params.cmo_set && file(params.cmo_set, type: 'any').exists() ? file(params.cmo_set, type: 'any') : ch_empty_file_4
g_11_2_g_55 = params.feature_reference && file(params.feature_reference, type: 'any').exists() ? file(params.feature_reference, type: 'any') : ch_empty_file_2
g_11_2_g_20 = params.feature_reference && file(params.feature_reference, type: 'any').exists() ? file(params.feature_reference, type: 'any') : ch_empty_file_2
g_12_3_g_55 = params.VDJ_reference && file(params.VDJ_reference, type: 'any').exists() ? file(params.VDJ_reference, type: 'any') : ch_empty_file_3
g_12_3_g_20 = params.VDJ_reference && file(params.VDJ_reference, type: 'any').exists() ? file(params.VDJ_reference, type: 'any') : ch_empty_file_3
g_17_2_g_13 = params.feature_ref && file(params.feature_ref, type: 'any').exists() ? file(params.feature_ref, type: 'any') : ch_empty_file_2
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

//* @style @array:{bcl_directory,mkfastq_sampleSheet} @multicolumn:{bcl_directory,mkfastq_sampleSheet}


process mkfastq_prep {


output:
 path "*"  ,emit:g_4_bcl00_g_0 
 val bcl_directory  ,emit:g_4_bcl_directory13_g_13 

container 'quay.io/ummsbiocore/pipeline_base_image:1.0'

when:
(params.run_mkfastq && (params.run_mkfastq == "yes")) || !params.run_mkfastq

script:
bcl_directory = params.mkfastq_prep.bcl_directory
mkfastq_sampleSheet = params.mkfastq_prep.mkfastq_sampleSheet
bcl_directory2 = bcl_directory.collect{ '"' + it + '"'}
mkfastq_sampleSheet2 = mkfastq_sampleSheet.collect{ '"' + it + '"'}
"""
#!/usr/bin/env python

import subprocess

def is_not_blank(s):
    return bool(s and not s.isspace())

def copy_or_link(src, dest, is_directory=False):
    if src.startswith('s3://'):
        cmd = ["aws", "s3", "cp", src, dest]
        if is_directory:
            cmd.append("--recursive")
    elif src.startswith('gs://'):
        # For gsutil, the recursive flag must be before the source
        cmd = ["gsutil", "-m", "cp"]
        if is_directory:
            cmd.append("-r")
        cmd.append(src)
        if is_directory:
            cmd.append(".")
        else:
            cmd.append(dest)
    else:
        cmd = ["ln", "-s", src, dest]
    print(cmd)
    subprocess.run(cmd)

bcl_directory = ${bcl_directory2}
mkfastq_sampleSheet = ${mkfastq_sampleSheet2}

for i in range(len(bcl_directory)):
    bcl = bcl_directory[i]
    sampleSheet = mkfastq_sampleSheet[i]
    # Remove last slash if exist
    bcl = bcl.rstrip('/')
    # Get name of the bcl directory
    bcl_name = bcl.rsplit('/', 1)[1]
    
    # Handle BCL
    if is_not_blank(bcl_name):
        copy_or_link(bcl + "/", bcl_name, is_directory=True)

    # Handle SampleSheet
    if not is_not_blank(sampleSheet):
        sampleSheet = bcl + "/" + "SampleSheet.csv"
    copy_or_link(sampleSheet, bcl_name + "/SampleSheet.csv")

"""

}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 4
    $MEMORY = 50
    $TIME = 1000
}
//* platform
//* platform
//* autofill

process mkfastq {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${bcl_files}_fastq$/) "mkfastq/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${bcl_files}_reports$/) "mkfastq/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${bcl_files}_laneBarcode.html$/) "mkfastq_report/$filename"}
input:
 path bcl_files

output:
 path "${bcl_files}_fastq"  ,emit:g_0_reads00_g_13 
 path "${bcl_files}_reports"  ,emit:g_0_outputDir11 
 path "${bcl_files}_laneBarcode.html" ,optional:true  ,emit:g_0_outputHTML22 

when:
(params.run_mkfastq && (params.run_mkfastq == "yes")) || !params.run_mkfastq

script:
cellranger_mkfastq_parameters = params.mkfastq.cellranger_mkfastq_parameters
"""	
cellranger mkfastq --id=mkfastq --run=${bcl_files} --csv=${bcl_files}/SampleSheet.csv --output-dir=fastq --rc-i2-override=true ${cellranger_mkfastq_parameters}
mv fastq ${bcl_files}_fastq
mkdir ${bcl_files}_reports 

mv ${bcl_files}_fastq/Reports ${bcl_files}_reports/.
mv ${bcl_files}_fastq/Stats ${bcl_files}_reports/.
cp ${bcl_files}_reports/Reports/html/*/all/all/all/laneBarcode.html ${bcl_files}_laneBarcode.html
"""

}

//libraries
fastq_id = params.cellranger_multi_prep.fastq_id
lanes_10x = params.cellranger_multi_prep.lanes_10x
physical_library_id = params.cellranger_multi_prep.physical_library_id
feature_types = params.cellranger_multi_prep.feature_types
//gene-expression
create_bam = params.cellranger_multi_prep.create_bam
expect_cells = params.cellranger_multi_prep.expect_cells
force_cells = params.cellranger_multi_prep.force_cells
cellranger_multi_chemistry = params.cellranger_multi_prep.cellranger_multi_chemistry
r1_length = params.cellranger_multi_prep.r1_length
include_introns = params.cellranger_multi_prep.include_introns
//* params.cmo_set =  ""  //* @input @single_file @optional @description:"Optional. CMO set CSV file, declaring CMO constructs and associated barcodes."
//feature
//* params.feature_reference =  ""  //* @input  @optional @single_file @description:"Feature reference CSV file, declaring Feature Barcode constructs and associated barcodes. Required for Feature Barcode libraries, otherwise optional."
r1_length_feature = params.cellranger_multi_prep.r1_length_feature
//vdj
//* params.VDJ_reference =  ""  //* @input  @optional @single_file @description:"VDJ reference CSV file."
r1_length_vdj = params.cellranger_multi_prep.r1_length_vdj
//samples
sample_id = params.cellranger_multi_prep.sample_id
lane_of_sample_10x = params.cellranger_multi_prep.lane_of_sample_10x
cmo_ids = params.cellranger_multi_prep.cmo_ids
hashtag_ids = params.cellranger_multi_prep.hashtag_ids
ocm_barcode_ids = params.cellranger_multi_prep.ocm_barcode_ids
probe_barcode_ids = params.cellranger_multi_prep.probe_barcode_ids
description = params.cellranger_multi_prep.description

fastq_id2 = fastq_id.collect{ '"' + it + '"'}
lanes_10x2 = lanes_10x.collect{ '"' + it + '"'}
physical_library_id2 = physical_library_id.collect{ '"' + it + '"'}
feature_types2 = feature_types.collect{ '"' + it + '"'}
sample_id2 = sample_id.collect{ '"' + it + '"'}
cmo_ids2 = cmo_ids.collect{ '"' + it + '"'}
hashtag_ids2 = hashtag_ids.collect{ '"' + it + '"'}
ocm_barcode_ids2 = ocm_barcode_ids.collect{ '"' + it + '"'}
probe_barcode_ids2 = probe_barcode_ids.collect{ '"' + it + '"'}

lane_of_sample2 = []
if (lane_of_sample_10x instanceof List) {
	lane_of_sample2 = lane_of_sample_10x.collect{ '"' + it + '"'}	
} 
description2 = description.collect{ '"' + it + '"'}

//* @style @spreadsheet:{fastq_id,lanes_10x,physical_library_id,feature_types},{sample_id,lane_of_sample_10x,cmo_ids,hashtag_ids,ocm_barcode_ids,probe_barcode_ids,description} 

process cellranger_multi_prep {


output:
 path "*.csv"  ,emit:g_5_csvFile05_g_20 

when:
(params.run_cellranger_multi && (params.run_cellranger_multi == "yes")) || !params.run_cellranger_multi

script:

"""
#!/usr/bin/env python

import subprocess,sys
import csv,os  

def is_not_blank(s):
	return bool(s and not s.isspace() and s != "null" and not s.startswith('NO_FILE'))

lanes = ${lanes_10x2}
lane_of_samples = ${lane_of_sample2}

# Group config files by using lanes:  lanes + lane_of_samples
all_lanes = lanes + lane_of_samples
unique_lanes = list(set(all_lanes))
print(unique_lanes)
# remove "" from list if it has "any" in it
if '' in unique_lanes :
	unique_lanes.remove('')
	if 'any' not in unique_lanes :
		unique_lanes.append("any")

# if unique_lanes has "any" and it has more than 1 lanes then remove "any"
# other lanes will use "any"
if (len(unique_lanes) > 1) and ('any' in unique_lanes):
	unique_lanes.remove('any')
	
print(unique_lanes)

for l in range(len(unique_lanes)):
	f = open('config_'+unique_lanes[l]+'.csv', "w")
	

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
 val bcl_directory  ,emit:g_25_bcl_directory05_g_13 
 path "reads_fastq"  ,emit:g_25_reads16_g_13 

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

//* params.feature_ref =  ""  //* @input  @optional @single_file @description:"Feature reference CSV file, optional."

//libraries
sample = params.cellranger_count_prep.sample
lane_of_library = params.cellranger_count_prep.lane_of_library
library_types = params.cellranger_count_prep.library_types
//gene-expression
cellranger_count_create_bam = params.cellranger_count_prep.cellranger_count_create_bam
cellranger_count_expect_cells = params.cellranger_count_prep.cellranger_count_expect_cells
cellranger_count_r1_length = params.cellranger_count_prep.cellranger_count_r1_length
cellranger_count_chemistry = params.cellranger_count_prep.cellranger_count_chemistry
cellranger_count_parameters = params.cellranger_count_prep.cellranger_count_parameters

sample2 = sample.collect{ '"' + it + '"'}
lane_of_library2 = lane_of_library.collect{ '"' + it + '"'}
library_types2 = library_types.collect{ '"' + it + '"'}
//* @style @spreadsheet:{sample,lane_of_library,library_types}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 4
    $MEMORY = 100
}
//* platform
//* platform
//* autofill

process cellranger_count_prep {


output:
 path "*.csv"  ,emit:g_27_csvFile04_g_13 

when:
(params.run_cellranger_count && (params.run_cellranger_count == "yes")) || !params.run_cellranger_count

script:

"""
#!/usr/bin/env python

import subprocess
import csv  

def is_not_blank(s):
	return bool(s and not s.isspace() and s != "null" and not s.startswith('NO_FILE'))

# Group config files by using lanes:
all_lanes = ${lane_of_library2}
unique_lanes = list(set(all_lanes))
print(unique_lanes)
# remove "" from list if it has "any" in it
if '' in unique_lanes :
	unique_lanes.remove('')
	if 'any' not in unique_lanes :
		unique_lanes.append("any")

# if unique_lanes has "any" and it has more than 1 lanes then remove "any"
# other lanes will use "any"
if (len(unique_lanes) > 1) and ('any' in unique_lanes):
	unique_lanes.remove('any')
	
print(unique_lanes)

for l in range(len(unique_lanes)):
	f = open('config_'+unique_lanes[l]+'.csv', "w")


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

container 'quay.io/ummsbiocore/pipeline_base_image:1.0'

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

shell:
'''
#!/usr/bin/env perl 

## Replace gene_id column with gene_name column in the gtf file
## Also check if any transcript_id defined in multiple chromosomes.
system("mkdir out");

open(OUT1, ">out/!{gtf}");
open(OUT2, ">notvalid_!{gtf}");
my %transcipt;
my $file = "!{gtf}";
open IN, $file;
while( my $line = <IN>)  {
    chomp;
    @a=split("\\t",$line);
    @attr=split(";",$a[8]);
    my %h;
    for my $elem (@attr) {
        ($first, $rest) = split ' ', $elem, 2;
        $h{$first} = $rest.";";
    }
    my $geneId = "";
    my $transcript_id = "";
    if (exists $h{"gene_name"}){
        $geneId = $h{"gene_name"};
    } elsif (exists $h{"gene_id"}){
        $geneId = $h{"gene_id"};
    }
    if (exists $h{"transcript_id"}){
        $transcript_id = $h{"transcript_id"};
    } elsif (exists $h{"transcript_name"}){
        $transcript_id = $h{"transcript_name"};
    } elsif (exists $h{"gene_id"}){
        $transcript_id = $h{"gene_id"};
    }
    if ($geneId ne "" && $transcript_id ne ""){
        ## check if any transcript_id defined in multiple chromosomes.
        if (exists $transcipt{$transcript_id}){
             if ($transcipt{$transcript_id} ne $a[0]){
               print OUT2 "$transcript_id: $transcipt{$transcript_id} vs $a[0]\\n";
                next;
                }
        } else {
             $transcipt{$transcript_id} = $a[0];
        }
        $a[8]=join(" ",("gene_id",$geneId,"transcript_id",$transcript_id));
        print OUT1 join("\\t",@a), "\\n";
    }  else {
        print OUT2 "$line";
    }
}
close OUT1;
close OUT2;
close IN;
'''
}


process Check_and_Build_Module_Add_custom_seq_to_genome_gtf {

input:
 path genome
 path gtf
 path custom_fasta

output:
 path "${genomeName}_custom.fa"  ,emit:g50_58_genome00_g50_52 
 path "${gtfName}_custom_sorted.gtf"  ,emit:g50_58_gtfFile10_g50_53 

when:
params.add_sequences_to_reference == "yes"

script:
genomeName = genome.baseName
gtfName = gtf.baseName
is_custom_genome_exists = custom_fasta.name.startsWith('NO_FILE') ? "False" : "True" 

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

container "${ params.IMAGE_BASE ? "${params.IMAGE_BASE}/rnaseq:4.0" : "quay.io/ummsbiocore/rnaseq:4.0" }"

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

container 'quay.io/ummsbiocore/pipeline_base_image:1.0'
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
 path "*/${transcriptome}" ,optional:true  ,emit:g_18_reference01_g_13 

container 'quay.io/ummsbiocore/pipeline_base_image:1.0'
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
 path cmo_set
 path config
 val bcl_directory
 val fastq_start_directory
 path fastq_start_reads

output:
 path "${run_id}_outs"  ,emit:g_20_outputDir01_g_54 
 path "*_web_summary.html"  ,emit:g_20_outputHTML11 
 path "final_${config}"  ,emit:g_20_csvFile22 
 path "run_bamtofastq" ,optional:true  ,emit:g_20_runFile30_g_54 
 path "*filtered_contig_annotations.csv" ,optional:true  ,emit:g_20_csvFile44 

when:
(params.run_cellranger_multi && (params.run_cellranger_multi == "yes")) || !params.run_cellranger_multi

script:
config_id = config.toString().replace(".csv","").replace("config_","")
run_id = "run_"+config_id
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
fastq_id = ${fastq_id2}
lanes10x = ${lanes_10x2}
physical_library_id = ${physical_library_id2}	
feature_types = ${feature_types2}
sample_id = ${sample_id2}
cmo_ids = ${cmo_ids2} 
hashtag_ids = ${hashtag_ids2} 
ocm_barcode_ids = ${ocm_barcode_ids2} 
probe_barcode_ids = ${probe_barcode_ids2} 
lane_of_samples = ${lane_of_sample2}
description = ${description2}


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


# if feature_types both VDJ-T and Multiplexing Capture then first run demultiplexing
demultiplexing = False
if 'VDJ-T' in feature_types and 'Multiplexing Capture' in feature_types:
	demultiplexing = True

unique_lane_id = "${config_id}"

with open(r'final_${config}', 'a') as f:
	w = csv.writer(f)
	w.writerow(["[gene-expression]"])
	ref_path = os.path.abspath("${ref}")
	w.writerow(["reference",ref_path])
	if is_not_blank("${cmo_set}") and ("${cmo_set}" != "NA"):
		cmo_path = os.path.abspath("${cmo_set}")
		w.writerow(["cmo-set",cmo_path])
	if is_not_blank("${expect_cells}"):
		w.writerow(["expect-cells","${expect_cells}"])
	if is_not_blank("${cellranger_multi_chemistry}"):
		w.writerow(["chemistry","${cellranger_multi_chemistry}"])
	if is_not_blank("${r1_length}"):
		w.writerow(["r1-length","${r1_length}"])
	if is_not_blank("${include_introns}"):
		w.writerow(["include-introns","${include_introns}"])
	if is_not_blank("${create_bam}"):
		w.writerow(["create-bam","${create_bam}"])
	if (is_not_blank("${feature_reference}") and ("${feature_reference}" != "NA")) or is_not_blank("${r1_length_feature}"):
		w.writerow(["[feature]"])
		if is_not_blank("${feature_reference}") and ("${feature_reference}" != "NA"):
			feature_reference_path = os.path.abspath("${feature_reference}")
			w.writerow(["reference",feature_reference_path])
		if is_not_blank("${r1_length_feature}"):
			w.writerow(["r1-length","${r1_length_feature}"])
	if is_not_blank("${VDJ_reference}") and ("${VDJ_reference}" != "NA"):
		w.writerow(["[vdj]"])
		vdj_path = os.path.abspath("${VDJ_reference}")
		w.writerow(["reference",vdj_path])
		if is_not_blank("${r1_length_vdj}"):
			w.writerow(["r1-length","${r1_length_vdj}"])

	if fastq_id:
		w.writerow(["[libraries]"])
		w.writerow(["fastq_id","fastqs","physical_library_id","feature_types"])
		for i in range(len(fastq_id)):
			if (lanes10x[i] == "any" or lanes10x[i] == "") or (unique_lane_id == lanes10x[i]) :
				for b in range(len(bcl_directory)):
					bcl_name = bcl_directory[b].rstrip('/').rsplit('/', 1)[1]
					print(bcl_name)
					for root, dirs, files in os.walk(bcl_name + "_fastq"):
						if any(fname.startswith(fastq_id[i]) for fname in files):
							bcl_path = os.path.abspath(root)
							print("directory: "+bcl_path)
							print(glob.glob(os.path.join(root, fastq_id[i] + "*")))
							if (demultiplexing == False or (demultiplexing == True and feature_types[i] != "VDJ-T")):
								w.writerow([fastq_id[i],bcl_path,physical_library_id[i],feature_types[i]])

	if is_nonempty_string_or_array(sample_id): 
		w.writerow(["[samples]"])
		w.writerow(["sample_id", nonempty_columns[0][0], "description"])
		for i in range(len(sample_id)):
			if (lane_of_samples[i] == "any" or lane_of_samples[i] == "") or (unique_lane_id == lane_of_samples[i]) :
				w.writerow([sample_id[i],nonempty_columns[0][1][i],description[i]])
				
print("## Config File:")
with open('final_${config}', 'r') as f:
	print(f.read())
    
cmd = ["cellranger","multi","--id=${run_id}","--csv=final_${config}","--localcores=16","--localmem=120"]
print(cmd)

p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = p.communicate()
print(stdout.decode())
print(stderr.decode())
if p.returncode != 0:
	sys.exit(p.returncode)

prefix = ""
if demultiplexing == True:
	subprocess.run("touch run_bamtofastq", shell=True, check=True)
	prefix = "demultiplexing_"
	subprocess.run("mv final_${config} demultiplexing_final_${config}", shell=True, check=True)
	subprocess.run("cp ${config} ${run_id}/outs/.", shell=True, check=True)
	
subprocess.run("mv ${run_id}/outs "+prefix+"${run_id}_outs && rm -rf ${run_id}", shell=True, check=True)
subprocess.run("for i in \$(ls "+prefix+"${run_id}_outs/per_sample_outs); do cp "+prefix+"${run_id}_outs/per_sample_outs/\${i}/web_summary.html "+prefix+"lane10x_${config_id}_\${i}_web_summary.html; done", shell=True, check=True)
try:
	subprocess.run("for i in \$(ls "+prefix+"${run_id}_outs/per_sample_outs); do cp "+prefix+"${run_id}_outs/per_sample_outs/\${i}/vdj_b/filtered_contig_annotations.csv "+prefix+"lane10x_${config_id}_\${i}_vdj_b_filtered_contig_annotations.csv; done", shell=True, check=True)
except subprocess.CalledProcessError:
	print("INFO: vdj_b/filtered_contig_annotations.csv was not found to publish it to report tab.")
try:
	subprocess.run("for i in \$(ls "+prefix+"${run_id}_outs/per_sample_outs); do cp "+prefix+"${run_id}_outs/per_sample_outs/\${i}/vdj_t/filtered_contig_annotations.csv "+prefix+"lane10x_${config_id}_\${i}_vdj_t_filtered_contig_annotations.csv; done", shell=True, check=True)
except subprocess.CalledProcessError:
	print("INFO: vdj_t/filtered_contig_annotations.csv was not found to publish it to report tab.")
try:
	subprocess.run("for i in \$(ls "+prefix+"${run_id}_outs/per_sample_outs); do cp "+prefix+"${run_id}_outs/per_sample_outs/\${i}/vdj_t_gd/filtered_contig_annotations.csv "+prefix+"lane10x_${config_id}_\${i}_vdj_t_gd_filtered_contig_annotations.csv; done", shell=True, check=True)
except subprocess.CalledProcessError:
	print("INFO: vdj_t_gd/filtered_contig_annotations.csv was not found to publish it to report tab.")


"""


}


process Multi_h5_explorer {

input:
 path output_dir

output:
 path "*.h5"  ,emit:g_59_h5_file01_g_60 

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


process bamtofastq_10x {

input:
 path run_bamtofastq
 path demultiplexed_folders

output:
 tuple file("config_*"), file("bamtofastq_filtered/*")  ,emit:g_54_config_reads06_g_55 

script:
"""
for i in \$(ls ${demultiplexed_folders}/per_sample_outs); do mkdir -p bamtofastq && bamtofastq --reads-per-fastq=2200000000 ${demultiplexed_folders}/per_sample_outs/\${i}/count/sample_alignments.bam bamtofastq/\${i}; done
cp ${demultiplexed_folders}/config_* .
mkdir -p bamtofastq_filtered
for i in \$(ls ${demultiplexed_folders}/per_sample_outs); do
	string=\$(samtools view -H ${demultiplexed_folders}/per_sample_outs/\${i}/count/sample_alignments.bam | grep "Gene Expression")
	gexid=\$(echo \${string#*library_id\\":} |cut -d ','  -f1)
	# demultiplexed_folders-> demultiplexing_run_12_outs
	foldername="${demultiplexed_folders}"
	a=\$(echo \${foldername#*demultiplexing_}) ##run_12_outs
	runid=\$(echo \${a%*_outs})  ##run_12
	echo "runid: \${runid}"
	echo "gexid: \${gexid}"
	mkdir -p bamtofastq_filtered/\${i}
	mv bamtofastq/\${i}/\${runid}_\${gexid}_*  bamtofastq_filtered/\${i}
done

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



process cellranger_multi_after_demultiplexing {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*${run_id}_outs$/) "multi/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_web_summary.html$/) "cellranger_report/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*final_${config}$/) "multi/$filename"}
input:
 path reads
 path ref
 path feature_reference
 path VDJ_reference
 path cmo_set
 val bcl_directory
 tuple file(config), file(bamtofastq)
 val fastq_start_directory

output:
 path "*${run_id}_outs"  ,emit:g_55_outputDir00 
 path "*_web_summary.html"  ,emit:g_55_outputHTML11 
 path "*final_${config}"  ,emit:g_55_csvFile22 

when:
(params.run_cellranger_multi && (params.run_cellranger_multi == "yes")) || !params.run_cellranger_multi

script:
config_id = config.toString().replace(".csv","").replace("config_","")
run_id = bamtofastq+"_run_"+config_id
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

bcl_directory = ${bcl_directory2}
fastq_id = ${fastq_id2}
lanes10x = ${lanes_10x2}
physical_library_id = ${physical_library_id2}	
feature_types = ${feature_types2}
sample_id = ${sample_id2}
cmo_ids = ${cmo_ids2} 
lane_of_samples = ${lane_of_sample2}
description = ${description2}


unique_lane_id = "${config_id}"

with open(r'${bamtofastq}_final_${config}', 'a') as f:
	w = csv.writer(f)
	w.writerow(["[gene-expression]"])
	ref_path = os.path.abspath("${ref}")
	w.writerow(["reference",ref_path])
	w.writerow(["check-library-compatibility","false"])
	if is_not_blank("${cmo_set}") and ("${cmo_set}" != "NA"):
		cmo_path = os.path.abspath("${cmo_set}")
		w.writerow(["cmo-set",cmo_path])
	if is_not_blank("${expect_cells}"):
		w.writerow(["expect-cells","${expect_cells}"])
	if is_not_blank("${cellranger_multi_chemistry}"):
		w.writerow(["chemistry","${cellranger_multi_chemistry}"])
	if is_not_blank("${r1_length}"):
		w.writerow(["r1-length","${r1_length}"])
	if is_not_blank("${include_introns}"):
		w.writerow(["include-introns","${include_introns}"])
	if is_not_blank("${create_bam}"):
		w.writerow(["create-bam","${create_bam}"])
	if (is_not_blank("${feature_reference}") and ("${feature_reference}" != "NA")) or is_not_blank("${r1_length_feature}"):
		w.writerow(["[feature]"])
		if is_not_blank("${feature_reference}") and ("${feature_reference}" != "NA"):
			feature_reference_path = os.path.abspath("${feature_reference}")
			w.writerow(["reference",feature_reference_path])
		if is_not_blank("${r1_length_feature}"):
			w.writerow(["r1-length","${r1_length_feature}"])
	if is_not_blank("${VDJ_reference}") and ("${VDJ_reference}" != "NA"):
		w.writerow(["[vdj]"])
		vdj_path = os.path.abspath("${VDJ_reference}")
		w.writerow(["reference",vdj_path])
		if is_not_blank("${r1_length_vdj}"):
			w.writerow(["r1-length","${r1_length_vdj}"])

	if fastq_id:
		w.writerow(["[libraries]"])
		w.writerow(["fastq_id","fastqs","physical_library_id","feature_types"])
		
		## To Support multiple recursive directories
		for path, currentDirectory, files in os.walk("${bamtofastq}"):
			if any(file.startswith("bamtofastq") for file in files):
				cur_path = os.path.abspath(".")
				w.writerow(["bamtofastq",cur_path+"/"+path,"","Gene Expression"])
		for i in range(len(fastq_id)):
			if (lanes10x[i] == "any" or lanes10x[i] == "") or (unique_lane_id == lanes10x[i]) :
				for b in range(len(bcl_directory)):
					bcl_name = bcl_directory[b].rstrip('/').rsplit('/', 1)[1]
					print(bcl_name)
					if any(fname.startswith(fastq_id[i]) for fname in os.listdir(bcl_name+"_fastq")):
						bcl_path = os.path.abspath(bcl_name+"_fastq")
						print("directory: "+bcl_path)
						print(glob.glob(bcl_name+"_fastq/"+fastq_id[i]+"*"))
						if (feature_types[i] != "Multiplexing Capture"  and feature_types[i] != "Gene Expression"):
							w.writerow([fastq_id[i],bcl_path,physical_library_id[i],feature_types[i]])

				
print("## Config File:")
with open('${bamtofastq}_final_${config}', 'r') as f:
	print(f.read())
    
cmd = ["cellranger","multi","--id=${run_id}","--csv=${bamtofastq}_final_${config}","--localcores=16","--localmem=120"]
print(cmd)

p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = p.communicate()
print(stdout.decode())
print(stderr.decode())
if p.returncode != 0:
	sys.exit(p.returncode)

prefix = ""
subprocess.run("mv ${run_id}/outs "+prefix+"${run_id}_outs && rm -rf ${run_id}", shell=True, check=True)
subprocess.run("for i in \$(ls "+prefix+"${run_id}_outs/per_sample_outs); do cp "+prefix+"${run_id}_outs/per_sample_outs/\${i}/web_summary.html "+prefix+"lane10x_${config_id}_\${i}_web_summary.html; done", shell=True, check=True)



"""


}

//* @style @array:{sample,library_types} @multicolumn:{sample,library_types}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 16
    $MEMORY = 120
}
//* platform
//* platform
//* autofill

process cellranger_count {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${run_id}_outs$/) "count/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${run_id}_web_summary.html$/) "cellranger_report/$filename"}
input:
 path reads
 path ref
 path feature_ref
 val bcl_directory
 path config
 val fastq_start_directory
 path fastq_start_reads

output:
 path "${run_id}_outs"  ,emit:g_13_outputDir00 
 path "${run_id}_web_summary.html"  ,emit:g_13_outputHTML11 
 path "${run_id}_raw_feature_bc_matrix.h5"  ,emit:g_13_h5_file22 
 path "${run_id}_filtered_feature_bc_matrix.h5"  ,emit:g_13_h5_file30_g_60 

when:
(params.run_cellranger_count && (params.run_cellranger_count == "yes")) || !params.run_cellranger_count

script:

config_id = config.toString().replace(".csv","").replace("config_","")
run_id = "run_"+config_id
bcl_directory2 = ""
if (bcl_directory){
	bcl_directory2 = bcl_directory.collect{ '"' + it + '"'}
} else if (fastq_start_directory){
	bcl_directory2 = fastq_start_directory.collect{ '"' + it + '"'}
}

"""
#!/usr/bin/env python

import subprocess
import csv  
import os
import sys
import glob

def is_not_blank(s):
	return bool(s and not s.isspace() and s != "null" and not s.startswith('NO_FILE'))

sample = ${sample2}
bcl_directory = ${bcl_directory2}
lane_of_library = ${lane_of_library2}
library_types = ${library_types2}
feature_reference_text = "--feature-ref=${feature_ref}" if is_not_blank("${feature_ref}") else ""
expect_cells_text = "--expect-cells=${cellranger_count_expect_cells}" if is_not_blank("${cellranger_count_expect_cells}") else ""
r1_length_text = "--r1-length=${cellranger_count_r1_length}" if is_not_blank("${cellranger_count_r1_length}") else ""

unique_lane_id = "${config_id}"
print("## All Files:")
for f in glob.glob('**/*', recursive=True):
	print(f)

with open(r'final_${config}', 'a') as f:
	w = csv.writer(f)
	if sample:
		w.writerow(["fastqs","sample","library_type"])
		for i in range(len(sample)):
			if (lane_of_library[i] == "any" or lane_of_library[i] == "") or (unique_lane_id == lane_of_library[i]) :
				for b in range(len(bcl_directory)):
					bcl_name = bcl_directory[b].rstrip('/').rsplit('/', 1)[1]
					## To Support multiple recursive directories
					for path, currentDirectory, files in os.walk(bcl_name+"_fastq"):
						if any(file.startswith(sample[i]+"_S") for file in files):
							cur_path = os.path.abspath(".")
							w.writerow([cur_path+"/"+path,sample[i].strip(),library_types[i]])
					# if any(fname.startswith(sample[i]) for fname in os.listdir(bcl_name+"_fastq")):
					# 	bcl_path = os.path.abspath(bcl_name+"_fastq")
					# 	w.writerow([bcl_path,sample[i],library_types[i]])

print("## Config File:")
with open('final_${config}', 'r') as f:
	print(f.read())

cmd = "cellranger count --localcores=16 --localmem=120 ${cellranger_count_parameters} --id=${run_id} --create-bam=${cellranger_count_create_bam} --chemistry=${cellranger_count_chemistry} --libraries=final_${config} --transcriptome=${ref}"+" "+r1_length_text+" "+feature_reference_text+ " "+ expect_cells_text
print(cmd)

p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = p.communicate()
print(stdout.decode())
print(stderr.decode())
if p.returncode != 0:
	sys.exit(p.returncode)

subprocess.run("mv ${run_id}/outs ${run_id}_outs && rm -rf ${run_id}", shell=True, check=True)
subprocess.run("cp ${run_id}_outs/web_summary.html ${run_id}_web_summary.html", shell=True, check=True)
subprocess.run("cp ${run_id}_outs/*raw*.h5 ${run_id}_raw_feature_bc_matrix.h5", shell=True, check=True)
subprocess.run("cp ${run_id}_outs/*filtered*.h5 ${run_id}_filtered_feature_bc_matrix.h5", shell=True, check=True)



"""


}


process h5_channel_collector {

input:
 path files1
 path files2

output:
 path "h5/*" ,optional:true  ,emit:g_60_h5_file00_g_61 

"""
mkdir h5 && mv *.h5 h5/ 2>/dev/null || true
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

container "quay.io/ummsbiocore/scrna_seurat:2.0"

when:
(params.run_scRNA_Analysis && (params.run_scRNA_Analysis == "yes")) || !params.run_scRNA_Analysis

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

container "quay.io/ummsbiocore/scrna_seurat:2.0"

shell:

'''
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

'''


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

container "quay.io/ummsbiocore/scrna_seurat:2.0"

shell:

varFeatures = params.scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction.varFeatures
selmethod = params.scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction.selmethod
Batch_Effect_Correction = params.scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction.Batch_Effect_Correction
WNN = params.scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction.WNN

//* @style @multicolumn:{varFeatures, selmethod},{Batch_Effect_Correction, WNN}

'''
#!/usr/bin/env Rscript

# libraries
library(Seurat)
library(dplyr)
#install.packages("harmony",repos = "http://cran.us.r-project.org")
library(harmony)

selmethod <- "!{selmethod}"
varFeatures <- "!{varFeatures}"

Data=readRDS("!{seurat_object}")
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
		if (as.logical("!{Batch_Effect_Correction}")){
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
		if (as.logical("!{Batch_Effect_Correction}")){
		Data=RunHarmony(Data,assay.use = DefaultAssay(Data),group.by.vars = "sample",max.iter.harmony = 10000,max.iter.cluster = 10000)
		}

	}
}

if ("!{WNN}"!="") {
original.assay=DefaultAssay(Data)

DefaultAssay(Data)="!{WNN}"

Data=NormalizeData(Data,normalization.method = "CLR",margin=2)

VariableFeatures(Data)=rownames(Data)

Data=ScaleData(Data)

Data=RunPCA(Data,reduction.name = "wpca")

if (Multi_sample==1) {
	Data=RunHarmony(Data,group.by.vars = "sample",assay.use = "!{WNN}",reduction = "wpca",reduction.save = "wharmony")

}

DefaultAssay(Data)=original.assay

}




#if (DefaultAssay(Data)=="SCT"){
#
#if (length(unique(Data$sample))==1) {
#	Data=RunPCA(Data,npcs=100)
#} else {
#	Data <- SplitObject(Data, split.by = "sample")
#	variable.features=SelectIntegrationFeatures(object.list = Data, nfeatures = as.numeric(varFeatures))
#	Data <- merge(Data[[1]],Data[-1])
#	VariableFeatures(Data) <- variable
#	if (all(Data[["percent.mt"]]==0)) {
#		Data <- ScaleData(Data)
#	} else {
#		Data <- ScaleData(Data,vars.to.regress="percent.mt")
#	}
#
#	Data=RunPCA(Data,npcs=100)
#	Data=RunHarmony(Data,assay.use = DefaultAssay(Data),group.by.vars = "sample",max.iter.harmony = 10000,max.iter.cluster = 10000)
#}
#} else {
#	Data <- FindVariableFeatures(Data,selection.method=selmethod,nfeatures=as.numeric(varFeatures))
#	if (all(Data[["percent.mt"]]==0)) {
#		Data <- ScaleData(Data)
#	} else {
#		Data <- ScaleData(Data,vars.to.regress="percent.mt")
#	}
#	Data=RunPCA(Data,npcs=100)
#	if (length(unique(Data$sample))>1) {
#		Data=RunHarmony(Data,assay.use = DefaultAssay(Data),group.by.vars = "sample",max.iter.harmony = 10000,max.iter.cluster = 10000)
#	}
#}
saveRDS(Data,"Reduced_and_Corrected.rds")

'''


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
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /Final_Analysis.rds$/) "scViewer/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tsv$/) "ClusterMarkers/$filename"}
input:
 path seurat_object

output:
 path "final_report.html"  ,emit:g51_19_outputHTML00 
 path "Final_Analysis.rds"  ,emit:g51_19_rdsFile10_g51_25 
 path "*.tsv"  ,emit:g51_19_outFileTSV22 

container "quay.io/ummsbiocore/scrna_seurat:2.0"

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
 path "*.h5ad"  ,emit:g51_22_h5ad_file00 

container "quay.io/ummsbiocore/scrna_seurat:2.0"

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

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 30
}
//* platform
//* platform
//* autofill

process scRNA_Analysis_Module_seurat_to_sce {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /sce_obj.rds$/) "Analysis_Apps/$filename"}
input:
 path seurat_object

output:
 path "sce_obj.rds"  ,emit:g51_25_rdsFile00_g51_27 

label 'scrna_seurat'

script:
"""
#!/usr/bin/env Rscript

library(Seurat)
library(SingleCellExperiment) 

seurat = readRDS("${seurat_object}")

sce = as.SingleCellExperiment(seurat)
saveRDS(sce, file='sce_obj.rds')
"""
}


process scRNA_Analysis_Module_launch_isee {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /launch_iSEE.R$/) "Analysis_Apps/$filename"}
input:
 path sce_object

output:
 path 'launch_iSEE.R'  ,emit:g51_27_outputFileOut00 


shell:

'''
#!/usr/bin/env perl

my $script = <<'EOF';

library(iSEE)

getPath <- function() {
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    needle <- "--file="
    match <- grep(needle, cmdArgs)
    if (length(match) > 0) {
        path <- dirname(normalizePath(sub(needle, "", cmdArgs[match]))[1])
        return(path)
    } else {
        return(normalizePath(getwd()))
    }
}

path = normalizePath(paste0(getPath()))
sce_file = "!{sce_object}"

sce = readRDS(normalizePath(paste0(path, '/', sce_file)))
app = iSEE(sce)
options(shiny.host = "0.0.0.0", shiny.port = 8789)
shiny::runApp(app)

EOF

open OUT, ">launch_iSEE.R";
print OUT $script;
close OUT;

'''
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
 path "Data.loom" ,optional:true  ,emit:g51_30_outputFileOut00 

container "quay.io/ummsbiocore/scrna_seurat:2.0"

shell:
Generate_loom_file = params.scRNA_Analysis_Module_SCEtoLOOM.Generate_loom_file

'''
#!/usr/bin/env Rscript

#library
library(Seurat)

if (as.logical("!{Generate_loom_file}")) {
	
Data=readRDS("!{seurat_object}")



annotation=read.csv("https://huggingface.co/datasets/ctheodoris/Genecorpus-30M/raw/main/example_input_files/gene_info_table.csv",header = T,row.names = 1)
annotation=annotation[annotation$gene_name%in%names(table(annotation$gene_name))[table(annotation$gene_name)==1],]
rownames(annotation)=annotation$gene_name
metadata=Data@meta.data
matrix=Data@assays$RNA@counts

matrix=matrix[rowSums(matrix)>0,]

matrix=matrix[intersect(rownames(matrix),rownames(annotation)),]

annotation=annotation[intersect(rownames(matrix),rownames(annotation)),]

rownames(matrix)=annotation$ensembl_id
matrix=matrix[rownames(matrix)[order(rownames(matrix),decreasing = F)],]

NewData=CreateSeuratObject(matrix,meta.data = metadata)

NewData.loom <- SeuratDisk::as.loom(NewData, filename = "Data.loom", verbose = FALSE)


}
'''


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

container "quay.io/ummsbiocore/scrna_seurat:2.0"

script:
	
"""
build_filtration_report.py --input-dir .

mkdir output
mv by_criteria_summary.tsv output
mv filtration_summary_report.Rmd output
mv overall_filtration_summary.tsv output
"""
}


workflow {


mkfastq_prep()
g_4_bcl00_g_0 = mkfastq_prep.out.g_4_bcl00_g_0
g_4_bcl_directory13_g_13 = mkfastq_prep.out.g_4_bcl_directory13_g_13
(g_4_bcl_directory15_g_55,g_4_bcl_directory16_g_20) = [g_4_bcl_directory13_g_13,g_4_bcl_directory13_g_13]


mkfastq(g_4_bcl00_g_0.flatten())
g_0_reads00_g_13 = mkfastq.out.g_0_reads00_g_13
(g_0_reads00_g_20,g_0_reads00_g_33) = [g_0_reads00_g_13,g_0_reads00_g_13]
g_0_outputDir11 = mkfastq.out.g_0_outputDir11
g_0_outputHTML22 = mkfastq.out.g_0_outputHTML22


cellranger_multi_prep()
g_5_csvFile05_g_20 = cellranger_multi_prep.out.g_5_csvFile05_g_20


cellranger_fastq_prep(g_23_0_g_22,g_24_1_g_22)
g_22_reads00_g_25 = cellranger_fastq_prep.out.g_22_reads00_g_25


cellranger_fastq_collect(g_22_reads00_g_25.collect())
g_25_bcl_directory05_g_13 = cellranger_fastq_collect.out.g_25_bcl_directory05_g_13
(g_25_bcl_directory07_g_55,g_25_bcl_directory07_g_20) = [g_25_bcl_directory05_g_13,g_25_bcl_directory05_g_13]
g_25_reads16_g_13 = cellranger_fastq_collect.out.g_25_reads16_g_13
(g_25_reads10_g_55,g_25_reads18_g_20) = [g_25_reads16_g_13,g_25_reads16_g_13]


cellranger_count_prep()
g_27_csvFile04_g_13 = cellranger_count_prep.out.g_27_csvFile04_g_13


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

Check_and_Build_Module_Add_custom_seq_to_genome_gtf(g50_21_genome00_g50_58,g50_57_gtfFile01_g50_58,g_40_2_g50_58)
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
g50_54_genome10_g_7 = Check_and_Build_Module_check_files.out.g50_54_genome10_g_7
g50_54_genomeSizes22 = Check_and_Build_Module_check_files.out.g50_54_genomeSizes22
g50_54_bed33 = Check_and_Build_Module_check_files.out.g50_54_bed33


cellranger_mkref(g50_54_genome10_g_7,g50_54_gtfFile01_g_7)
g_7_reference00_g_18 = cellranger_mkref.out.g_7_reference00_g_18

g_7_reference00_g_18= g_7_reference00_g_18.ifEmpty(ch_empty_file_1) 


cellranger_ref_checker(g_7_reference00_g_18)
g_18_reference01_g_13 = cellranger_ref_checker.out.g_18_reference01_g_13
(g_18_reference01_g_55,g_18_reference01_g_20) = [g_18_reference01_g_13,g_18_reference01_g_13]

g_0_reads00_g_20= g_0_reads00_g_20.ifEmpty(ch_empty_file_1) 
g_4_bcl_directory16_g_20= g_4_bcl_directory16_g_20.ifEmpty("") 
g_25_bcl_directory07_g_20= g_25_bcl_directory07_g_20.ifEmpty("") 
g_25_reads18_g_20= g_25_reads18_g_20.ifEmpty(ch_empty_file_7) 


if (!((params.run_cellranger_multi && (params.run_cellranger_multi == "yes")) || !params.run_cellranger_multi)){
g_5_csvFile05_g_20.set{g_20_csvFile22}
g_20_outputDir01_g_54 = Channel.empty()
g_20_outputDir00_g_59 = Channel.empty()
g_20_outputHTML11 = Channel.empty()
g_20_runFile30_g_54 = Channel.empty()
} else {

cellranger_multi(g_0_reads00_g_20.collect(),g_18_reference01_g_20,g_11_2_g_20,g_12_3_g_20,g_10_4_g_20,g_5_csvFile05_g_20.flatten(),g_4_bcl_directory16_g_20,g_25_bcl_directory07_g_20,g_25_reads18_g_20.collect())
g_20_outputDir01_g_54 = cellranger_multi.out.g_20_outputDir01_g_54
(g_20_outputDir00_g_59) = [g_20_outputDir01_g_54]
g_20_outputHTML11 = cellranger_multi.out.g_20_outputHTML11
g_20_csvFile22 = cellranger_multi.out.g_20_csvFile22
g_20_runFile30_g_54 = cellranger_multi.out.g_20_runFile30_g_54
g_20_csvFile44 = cellranger_multi.out.g_20_csvFile44
}


Multi_h5_explorer(g_20_outputDir00_g_59)
g_59_h5_file01_g_60 = Multi_h5_explorer.out.g_59_h5_file01_g_60


bamtofastq_10x(g_20_runFile30_g_54,g_20_outputDir01_g_54)
g_54_config_reads06_g_55 = bamtofastq_10x.out.g_54_config_reads06_g_55

g_25_reads10_g_55= g_25_reads10_g_55.ifEmpty(ch_empty_file_1) 
g_4_bcl_directory15_g_55= g_4_bcl_directory15_g_55.ifEmpty("") 
g_25_bcl_directory07_g_55= g_25_bcl_directory07_g_55.ifEmpty("") 


cellranger_multi_after_demultiplexing(g_25_reads10_g_55.collect(),g_18_reference01_g_55,g_11_2_g_55,g_12_3_g_55,g_10_4_g_55,g_4_bcl_directory15_g_55,g_54_config_reads06_g_55.transpose(),g_25_bcl_directory07_g_55)
g_55_outputDir00 = cellranger_multi_after_demultiplexing.out.g_55_outputDir00
g_55_outputHTML11 = cellranger_multi_after_demultiplexing.out.g_55_outputHTML11
g_55_csvFile22 = cellranger_multi_after_demultiplexing.out.g_55_csvFile22

g_0_reads00_g_13= g_0_reads00_g_13.ifEmpty(ch_empty_file_1) 
g_4_bcl_directory13_g_13= g_4_bcl_directory13_g_13.ifEmpty("") 
g_25_bcl_directory05_g_13= g_25_bcl_directory05_g_13.ifEmpty("") 
g_25_reads16_g_13= g_25_reads16_g_13.ifEmpty(ch_empty_file_5) 


cellranger_count(g_0_reads00_g_13.collect(),g_18_reference01_g_13,g_17_2_g_13,g_4_bcl_directory13_g_13,g_27_csvFile04_g_13.flatten(),g_25_bcl_directory05_g_13,g_25_reads16_g_13.collect())
g_13_outputDir00 = cellranger_count.out.g_13_outputDir00
g_13_outputHTML11 = cellranger_count.out.g_13_outputHTML11
g_13_h5_file22 = cellranger_count.out.g_13_h5_file22
g_13_h5_file30_g_60 = cellranger_count.out.g_13_h5_file30_g_60

g_13_h5_file30_g_60= g_13_h5_file30_g_60.ifEmpty(ch_empty_file_1) 
g_59_h5_file01_g_60= g_59_h5_file01_g_60.ifEmpty(ch_empty_file_2) 


h5_channel_collector(g_13_h5_file30_g_60.collect(),g_59_h5_file01_g_60.collect())
g_60_h5_file00_g_61 = h5_channel_collector.out.g_60_h5_file00_g_61


file_to_set_conversion_for_h5(g_60_h5_file00_g_61.flatten())
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
g51_19_rdsFile10_g51_25 = scRNA_Analysis_Module_Clustering_and_Find_Markers.out.g51_19_rdsFile10_g51_25
(g51_19_rdsFile10_g51_22,g51_19_rdsFile10_g51_30) = [g51_19_rdsFile10_g51_25,g51_19_rdsFile10_g51_25]
g51_19_outFileTSV22 = scRNA_Analysis_Module_Clustering_and_Find_Markers.out.g51_19_outFileTSV22


scRNA_Analysis_Module_Create_h5ad(g51_19_rdsFile10_g51_22)
g51_22_h5ad_file00 = scRNA_Analysis_Module_Create_h5ad.out.g51_22_h5ad_file00


scRNA_Analysis_Module_seurat_to_sce(g51_19_rdsFile10_g51_25)
g51_25_rdsFile00_g51_27 = scRNA_Analysis_Module_seurat_to_sce.out.g51_25_rdsFile00_g51_27


scRNA_Analysis_Module_launch_isee(g51_25_rdsFile00_g51_27)
g51_27_outputFileOut00 = scRNA_Analysis_Module_launch_isee.out.g51_27_outputFileOut00


scRNA_Analysis_Module_SCEtoLOOM(g51_19_rdsFile10_g51_30)
g51_30_outputFileOut00 = scRNA_Analysis_Module_SCEtoLOOM.out.g51_30_outputFileOut00


scRNA_Analysis_Module_filter_summary(g51_0_outFileTSV20_g51_34.collect())
g51_34_outputFileHTML00 = scRNA_Analysis_Module_filter_summary.out.g51_34_outputFileHTML00


}

workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
