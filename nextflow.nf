$HOSTNAME = ""
params.outdir = 'results'  

params.publishdict = [:]

def pathChecker(input, path, type){
	cmd = "mkdir -p check && mv ${input} check/. "
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

g_10_4_g_20 = params.cmo_set && file(params.cmo_set, type: 'any').exists() ? file(params.cmo_set, type: 'any') : ch_empty_file_4
g_10_4_g_55 = params.cmo_set && file(params.cmo_set, type: 'any').exists() ? file(params.cmo_set, type: 'any') : ch_empty_file_4
g_11_2_g_20 = params.feature_reference && file(params.feature_reference, type: 'any').exists() ? file(params.feature_reference, type: 'any') : ch_empty_file_2
g_11_2_g_55 = params.feature_reference && file(params.feature_reference, type: 'any').exists() ? file(params.feature_reference, type: 'any') : ch_empty_file_2
g_12_3_g_20 = params.VDJ_reference && file(params.VDJ_reference, type: 'any').exists() ? file(params.VDJ_reference, type: 'any') : ch_empty_file_3
g_12_3_g_55 = params.VDJ_reference && file(params.VDJ_reference, type: 'any').exists() ? file(params.VDJ_reference, type: 'any') : ch_empty_file_3
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
 path "${bcl_files}_fastq"  ,emit:g_0_reads00_g_33 
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
cmo_ids = params.cellranger_multi_prep.cmo_ids
lane_of_sample_10x = params.cellranger_multi_prep.lane_of_sample_10x
description = params.cellranger_multi_prep.description

fastq_id2 = fastq_id.collect{ '"' + it + '"'}
lanes_10x2 = lanes_10x.collect{ '"' + it + '"'}
physical_library_id2 = physical_library_id.collect{ '"' + it + '"'}
feature_types2 = feature_types.collect{ '"' + it + '"'}
sample_id2 = sample_id.collect{ '"' + it + '"'}
cmo_ids2 = cmo_ids.collect{ '"' + it + '"'}
lane_of_sample2 = []
if (lane_of_sample_10x instanceof List) {
	lane_of_sample2 = lane_of_sample_10x.collect{ '"' + it + '"'}	
} 
description2 = description.collect{ '"' + it + '"'}

//* @style @spreadsheet:{fastq_id,lanes_10x,physical_library_id,feature_types},{sample_id,cmo_ids,lane_of_sample_10x,description} 

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
 path "allreads/*/*_{R1,R2}_001.fastq.gz", includeInputs: true   ,emit:g_33_reads00_g_34 
 val "single"  ,emit:g_33_mate10_g_35 

stageInMode 'copy'

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
//* platform
//* platform
//* autofill

process MultiQC {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /multiqc_report.html$/) "multiqc/$filename"}
input:
 path "fastqc/*"

output:
 path "multiqc_report.html" ,optional:true  ,emit:g_52_outputHTML00 

errorStrategy 'ignore'

script:
multiqc_parameters = params.MultiQC.multiqc_parameters
"""
multiqc ${multiqc_parameters} -e general_stats -d -dd 2 .
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

bcl_directory = ${bcl_directory2}
fastq_id = ${fastq_id2}
lanes10x = ${lanes_10x2}
physical_library_id = ${physical_library_id2}	
feature_types = ${feature_types2}
sample_id = ${sample_id2}
cmo_ids = ${cmo_ids2} 
lane_of_samples = ${lane_of_sample2}
description = ${description2}

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

	if sample_id: 
		w.writerow(["[samples]"])
		w.writerow(["sample_id","cmo_ids","description"])
		for i in range(len(sample_id)):
			if (lane_of_samples[i] == "any" or lane_of_samples[i] == "") or (unique_lane_id == lane_of_samples[i]) :
				w.writerow([sample_id[i],cmo_ids[i],description[i]])
				
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

process scRNA_Analysis_Module_Load_Data_h5 {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_filtering_report.html$/) "QC_Reports/$filename"}
input:
 tuple val(name), file(Input)
 path Metadata

output:
 path "${name}.rds"  ,emit:g51_0_rdsFile00_g51_14 
 tuple val(name),file("${name}_filtering_report.html")  ,emit:g51_0_outputFileHTML11 

label 'scrna_seurat'

when:
(params.run_scRNA_Analysis && (params.run_scRNA_Analysis == "yes")) || !params.run_scRNA_Analysis

shell:
minUMI = params.scRNA_Analysis_Module_Load_Data_h5.minUMI
maxUMI = params.scRNA_Analysis_Module_Load_Data_h5.maxUMI
minFeature = params.scRNA_Analysis_Module_Load_Data_h5.minFeature
maxFeature = params.scRNA_Analysis_Module_Load_Data_h5.maxFeature
percent_mt = params.scRNA_Analysis_Module_Load_Data_h5.percent_mt
percent_ribo = params.scRNA_Analysis_Module_Load_Data_h5.percent_ribo
varFeatures = params.scRNA_Analysis_Module_Load_Data_h5.varFeatures
normMethod = params.scRNA_Analysis_Module_Load_Data_h5.normMethod
DoubletRemoval = params.scRNA_Analysis_Module_Load_Data_h5.DoubletRemoval
RemoveMitoGenes = params.scRNA_Analysis_Module_Load_Data_h5.RemoveMitoGenes
RemoveRiboGenes = params.scRNA_Analysis_Module_Load_Data_h5.RemoveRiboGenes

'''
#!/usr/bin/env perl

my $script = <<'EOF';
---
title: "Single Cell RNA-Seq QC and Filtering Report"
author: "Zhaorong Li"
output: 
  html_document:
    toc: true
    toc_float:
      toc_collapsed: true
      toc_depth: 3
    number_sections: true
    fig_caption: yes
    theme: cerulean
    code_folding: hide
editor_options: 
  chunk_output_type: console
params:
  SampleName: ''
  SamplePath: ''
  RawInput: FALSE
  nFeature_RNA_lower_threshold: 0.01
  nFeature_RNA_higher_threshold: 0.99
  nCount_RNA_lower_threshold: 0.01
  nCount_RNA_higher_threshold: 0.99
  ribosomal_contents_threshold: 25
  mitochondrial_contents_threshold: 50
  doubletpercentage: 0.01
  filtered_output: ""
  normMethod: "LogNormalize"
  varFeatures: 3000
  DoubletRemoval: TRUE
  Metadata: ''
  RemoveMitoGenes: FALSE
  RemoveRiboGenes: FALSE

---


```{r setup, include=FALSE}

suppressPackageStartupMessages({
library(dplyr)
library(Seurat)
library(ggplot2)
library(DropletUtils)
library(DoubletFinder)
library(clustree)

})

```

```{r helper_functions, include=FALSE}

get_present_pseudogenes <- function(seurat_object, pseudogenes) {

  # Ensure that the Seurat object has rownames (gene names)

  if (is.null(rownames(seurat_object))) {

    stop("The Seurat object does not have gene names as rownames.")

  }

  

  gene_names <- rownames(seurat_object)

  pseudogenes_in_dataset <- intersect(pseudogenes, gene_names)

  

  return(pseudogenes_in_dataset)

}


```

# Read in samples

If the data is a raw Count Matrix, meaning that the empty droplets are not filtered out by cellranger pipeline or Drop-Seq pipeline, the emptyDrops algorithm from DropUtils will be run in order to remove the empty droplets.

```{r read in samples}
Data=Read10X_h5(params$SamplePath)

MultiModal=F
if (is.list(Data)) {
  MultiModal=T
  Raw=Data
  Data=Data[["Gene Expression"]]
  
}
#If any cells with 0 in the dataset, raw input is true
if (any(colSums(Data)==0)) {
	RawInput=TRUE
} else {
	RawInput=FALSE
}


if (params$RemoveRiboGenes) {
	
	
if (any(grepl("^RP[SL]",rownames(Data)))) {
	rb.genes <- rownames(Data)[grep("^RP[SL]",rownames(Data))]
	Data=Data[!rownames(Data)%in%rb.genes,]
}
else if (any(grepl("^Rp[sl]",rownames(Data)))) {
	rb.genes <- rownames(Data)[grep("^Rp[sl]",rownames(Data))]
	GTgenes=c("Gm42418","AY036118")
	rb.genes <-c(rb.genes,GTgenes)
	Data=Data[!rownames(Data)%in%rb.genes,]
}
}

mt.genes=c()
if (params$RemoveMitoGenes) {


if (any(grepl("^MT-",rownames(Data)))) {
	mt.genes <- rownames(Data)[grep("^MT-",rownames(Data))]
	Data=Data[!rownames(Data)%in%mt.genes,]
}
else if (any(grepl("^mt-",rownames(Data)))) {
	mt.genes <- rownames(Data)[grep("^mt-",rownames(Data))]
	Data=Data[!rownames(Data)%in%mt.genes,]
}
}

if (RawInput) {
  empty=emptyDrops(Data[!rownames(Data)%in%mt.genes,])
  empty=data.frame(empty)
  empty=empty[!is.na(empty$FDR),]
  empty$DropIdentity=ifelse(empty$FDR<0.001,yes="Non Empty",
                            no="Empty")
  ggplot(empty,aes(x=DropIdentity,y=Total))+geom_bar(stat = 'identity')+xlab("Droplet classification")+ylab("Number of UMIs per cell")+ggtitle("Empty Droplet classification")
  Data=Data[,rownames(empty)[empty$FDR<0.05]]

}


```

# QC {.tabset}

In this section the number of genes, number of UMIs and the percentage of mitochondrial contents and ribosomal contents will be visualized.

Mitochondrial contents and ribosomal contents will be calculated based on the mitochondrial and ribosomal genes.

Cells with very high mitochondrial and ribosomal contents will bias the downstream clustering and differential expression analysis.

```{r QC, fig.width=5,fig.height=5}

Data=CreateSeuratObject(Data)
Data$sample=params$SampleName
Data$orig.ident=params$SampleName

if (any(grepl("^MT-",rownames(Data)))|any(grepl("^mt-",rownames(Data)))) {
	if (any(grepl("^MT-",rownames(Data)))) {
	Data[["percent.mt"]]=PercentageFeatureSet(Data,pattern="^MT-")
	rb.genes <- rownames(Data)[grep("^RP[SL]",rownames(Data))]
	Data[["percent.ribo"]] <- PercentageFeatureSet(Data, features = rb.genes)	


}
else if (any(grepl("^mt-",rownames(Data)))) {
	Data[["percent.mt"]]=PercentageFeatureSet(Data,pattern="^mt-")
	
	rb.genes <- rownames(Data)[grep("^Rp[sl]",rownames(Data))]
	present_pseudogenes <- get_present_pseudogenes(Data, c("Gm42418","AY036118") )
    rb.genes <- c( rb.genes, present_pseudogenes )
	
	Data[["percent.ribo"]] <- PercentageFeatureSet(Data, features = rb.genes)	

}} else {
	Data[["percent.mt"]]=0
		Data[["percent.ribo"]]=0

}

Unfiltered = Data

```

## Violin plot of number of gene per cell

``` {r vnFeature, fig.width=5,fig.height=5}
VlnPlot(Unfiltered,features = c("nFeature_RNA"),pt.size = 0)
```

## Violin plot of number of UMIs per cell

``` {r vnCount, fig.width=5,fig.height=5}
VlnPlot(Unfiltered,features = c("nCount_RNA"),pt.size = 0)
```

## Violin plot of mitochondrial contents per cell

``` {r vmt, fig.width=5,fig.height=5}
VlnPlot(Unfiltered,features = c("percent.mt"),pt.size = 0)
```

## Violin plot of ribosomal contents per cell

``` {r vrb, fig.width=5,fig.height=5}
VlnPlot(Unfiltered,features = c("percent.ribo"),pt.size = 0)
```

## Scatter plot of number of genes per cell vs number of UMIs per cell

``` {r sfu, fig.width=5,fig.height=5}
FeatureScatter(Unfiltered,
               feature1 = "nFeature_RNA",
               feature2 = "nCount_RNA")+ NoLegend()
```

## Scatter plot of mitochondrial contents

``` {r smfu, fig.width=5,fig.height=5}
if (isFALSE(all(Unfiltered[["percent.mt"]]==0))){
  FeatureScatter(Unfiltered,
               feature1 = "nFeature_RNA",
               feature2 = "percent.mt")+ NoLegend()

  FeatureScatter(Unfiltered,
               feature1 = "nCount_RNA",
               feature2 = "percent.mt")+ NoLegend()
  }

```

## Scatter plot of mitochondrial contents

``` {r srfu, fig.width=5,fig.height=5}
if (isFALSE(all(Unfiltered[["percent.ribo"]]==0))){
  FeatureScatter(Unfiltered,
               feature1 = "nFeature_RNA",
               feature2 = "percent.ribo")+ NoLegend()

  FeatureScatter(Unfiltered,
               feature1 = "nCount_RNA",
               feature2 = "percent.ribo")+ NoLegend()
}

```

# Doublet Classification {.tabset}

Doublets, or sometimes called multiplets, are the droplets which include two or more cells. Including these droplets in the downstream analysis will bias the results because these droplets include gene expression profiles of more than 1 cell.

DoubletFinder is used to classify the doublet. The violin plot in this section will show the number features and UMIs of the doublets vs that of non-doublets.

## Doublet Simulation and detection 
```{r Doublet Removal, error=FALSE, fig.height=5, fig.width=5, message=FALSE, warning=FALSE,results = FALSE}

DoubletRemovalHandle=as.logical(params$DoubletRemoval)

Data$Doublet.Classification="SingleLet"

if (is.na(DoubletRemovalHandle)) {
	if (MultiModal) {
		if (any(grepl("Multiplexing",names(Raw)))) {
		
			DoubletRemovalHandle=FALSE
			
		} else {
		
			DoubletRemovalHandle=TRUE
			
		}
	}
	else {
	DoubletRemovalHandle=TRUE

	}
} 

if (DoubletRemovalHandle) {
DoubletRemoval <- NormalizeData(Data)
DoubletRemoval <- FindVariableFeatures(DoubletRemoval)
DoubletRemoval <- ScaleData(DoubletRemoval)

DoubletRemoval <- RunPCA(DoubletRemoval,npcs=100)

pc.changes=diff(diff(DoubletRemoval@reductions$pca@stdev))
pc.changes=abs(pc.changes)
pc.changes=which(pc.changes>=mean(pc.changes))

DoubletRemoval=FindNeighbors(DoubletRemoval,dims = 1:(max(pc.changes)+2),reduction = "pca")

DoubletRemoval=FindClusters(DoubletRemoval,resolution = seq(2.0,0.1,-0.1))

names=paste0(DefaultAssay(DoubletRemoval),"_snn_res.")
SC3_Stability=clustree(DoubletRemoval,prefix = names)
SC3_Stability.results=SC3_Stability$data
SC3_Stability.results=SC3_Stability.results[,c(names,"sc3_stability")]
colnames(SC3_Stability.results)[1]="resolution"
SC3_Stability.results.mean=aggregate(sc3_stability~resolution,SC3_Stability.results,mean)
colnames(SC3_Stability.results.mean)[2]="sc3_stability_mean"
Idents(DoubletRemoval)=paste0(DefaultAssay(DoubletRemoval),"_snn_res.",max(as.numeric(as.character(SC3_Stability.results.mean$resolution))[SC3_Stability.results.mean$sc3_stability_mean==max(SC3_Stability.results.mean$sc3_stability_mean)]))

DoubletRemoval$seurat_clusters=Idents(DoubletRemoval)



params_selection <- paramSweep_v3(DoubletRemoval, PCs = 1:(max(pc.changes)+2), sct = F)
params_selection <- summarizeSweep(params_selection, GT = FALSE)
params_selection <- find.pK(params_selection)



annotations <- DoubletRemoval@meta.data$DoubletRemoval$seurat_clusters

homotypic.prop <- modelHomotypic(annotations) 
nExp_poi <- round(params$doubletpercentage*nrow(DoubletRemoval@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


DoubletRemoval <- doubletFinder_v3(DoubletRemoval, PCs = 1:(max(pc.changes)+2), pN = 0.25, pK =as.numeric(as.character(params_selection$pK))[params_selection$BCmetric==max(params_selection$BCmetric)], nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
doublet.classification.name=colnames(DoubletRemoval@meta.data)[ncol(DoubletRemoval@meta.data)]
DoubletRemoval$Doublet.Classification=DoubletRemoval@meta.data[,doublet.classification.name]

Data$Doublet.Classification=DoubletRemoval$Doublet.Classification
}

```

## Doublet Visualization

```{r Doublet Visualization, error=FALSE, fig.height=5, fig.width=5, message=FALSE, warning=FALSE,results = FALSE }

if (DoubletRemovalHandle) {
VlnPlot(DoubletRemoval,c("nCount_RNA","nFeature_RNA"),group.by = "Doublet.Classification",pt.size = 0)+NoLegend()
}

```

# Filtering {.tabset}

Based on the quantile information, cells with too low or too high number of features and UMIs are filtered out, which can be observed on the violin plot.

## Violin Plots of number of genes and UMIs of doublets and singlets

```{r Filtering, fig.width=5,fig.height=5,error=FALSE,warning=FALSE,message=FALSE}

if (DoubletRemovalHandle){
Data$Filtering = ifelse(
  Data$Doublet.Classification!='Doublet'&
  Data$nFeature_RNA>quantile(Data$nFeature_RNA,params$nFeature_RNA_lower_threshold)&
  Data$nFeature_RNA<quantile(Data$nFeature_RNA,params$nFeature_RNA_higher_threshold)&
  Data$nCount_RNA>quantile(Data$nCount_RNA,params$nCount_RNA_lower_threshold)&
  Data$nCount_RNA<quantile(Data$nCount_RNA,params$nCount_RNA_higher_threshold)&
  Data$percent.mt<params$mitochondrial_contents_threshold&
  Data$percent.ribo<params$ribosomal_contents_threshold,
  yes="Keep",
  no="Drop"
)

} else {
Data$Filtering = ifelse(
  Data$nFeature_RNA>quantile(Data$nFeature_RNA,params$nFeature_RNA_lower_threshold)&
  Data$nFeature_RNA<quantile(Data$nFeature_RNA,params$nFeature_RNA_higher_threshold)&
  Data$nCount_RNA>quantile(Data$nCount_RNA,params$nCount_RNA_lower_threshold)&
  Data$nCount_RNA<quantile(Data$nCount_RNA,params$nCount_RNA_higher_threshold)&
  Data$percent.mt<params$mitochondrial_contents_threshold&
  Data$percent.ribo<params$ribosomal_contents_threshold,
  yes="Keep",
  no="Drop")
}

VlnPlot(Data,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo",pt.size = 0),group.by = "Filtering",ncol = 2,pt.size = 0)+NoLegend()

Data=subset(Data,subset=Filtering=="Keep")

```

## Barplot of filtering results

```{r Filtering barplot, fig.width=5,fig.height=5,error=FALSE,warning=FALSE,message=FALSE}

Filtering_statistics=data.frame(Category=c("Unfiltered","Filtered"),CellNumber=c(ncol(Unfiltered),ncol(Data)))
Filtering_statistics$Category=factor(Filtering_statistics$Category,levels = c("Unfiltered","Filtered"))

ggplot(Filtering_statistics,aes(x=Category,y=CellNumber,label=CellNumber))+geom_bar(stat="identity")+geom_text(size = 3, position = position_stack(vjust = 0.75))

```


# Normalization

As the violin plots shown in the QC section, the sequencing depth and coverage of each cell in a single cell RNA-Seq dataset vary significantly.

The normalization step normalize the gene expression profile of each cell, which makes them comparable to each other in the downstream analysis.

The SCTransform is recommended as it enhances the biological signature in the data, however it is quite time-consuming and memory-consuming. 

The LogNormalize is very standard practice time-efficient.

```{r Normalization, error=FALSE, fig.height=5, fig.width=5, message=FALSE, warning=FALSE,results = FALSE}

tryCatch({if (file.exists(params$Metadata)) {
	Metadata=read.table(params$Metadata,sep="\t",check.names=F,header=T,row.names=NULL)
	if ("Sample"%in%colnames(Metadata)) {
	AttributeList=colnames(Metadata)[colnames(Metadata)!="Sample"]
	if (params$SampleName%in%Metadata$Sample) {
		for (i in 1:length(AttributeList)) {
			Data[[AttributeList[i]]]=as.character(Metadata[Metadata$Sample==params$SampleName,AttributeList[i]])
		}
	} else {
		for (i in 1:length(AttributeList)) {
			Data[[AttributeList[i]]]=""
		}
	}

}
}
},
error=function(err) {
	print("No metadata information is added to dataset")
})

if (params$normMethod=="SCT") {
	if (all(Data[["percent.mt"]]==0)) {
		Data <- SCTransform(Data,variable.features.n=params$varFeatures)
	} else {
		Data <- SCTransform(Data,variable.features.n=params$varFeatures,vars.to.regress="percent.mt")

	}
	} else {
	Data <- NormalizeData(object = Data,normalization.method=params$normMethod)
}

if (MultiModal) {
  for (name in names(Raw)[names(Raw)!="Gene Expression"]) {
    Data[[gsub(" ","",name)]]=CreateAssayObject((Raw[[name]][,colnames(Data)]))
  }

}

saveRDS(Data,params$filtered_output)


```

EOF

open OUT, ">!{name}_filtering_rmark.rmd";
print OUT $script;
close OUT;

runCommand("Rscript -e 'rmarkdown::render(\\"!{name}_filtering_rmark.rmd\\",\\"html_document\\", output_file = \\"!{name}_filtering_report.html\\",
params = list(SampleName=\\"!{name}\\",
SamplePath=\\"!{Input}\\",
nFeature_RNA_lower_threshold=as.numeric(!{minFeature}),
nFeature_RNA_higher_threshold=as.numeric(!{maxFeature}),
nCount_RNA_lower_threshold=as.numeric(!{minUMI}),
nCount_RNA_higher_threshold=as.numeric(!{maxUMI}),
ribosomal_contents_threshold=as.numeric(!{percent_ribo}),
mitochondrial_contents_threshold=as.numeric(!{percent_mt}),
filtered_output=\\"!{name}.rds\\",
normMethod=\\"!{normMethod}\\",
varFeatures=as.numeric(!{varFeatures}),
DoubletRemoval=as.character(\\"!{DoubletRemoval}\\"),
Metadata=\\"!{Metadata}\\",
RemoveMitoGenes=as.logical(\\"!{RemoveMitoGenes}\\"),
RemoveRiboGenes=as.logical(\\"!{RemoveRiboGenes}\\")
))'");

sub runCommand {
            my ($com) = @_;
            my $error = system($com);
            if   ($error) { die "Command failed: $error $com\\n"; }
            else          { print "Command successful: $com\\n"; }
          }

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

process scRNA_Analysis_Module_Merge_Seurat_Objects {

input:
 path seurat_obj

output:
 path "merged_filtered_seurat.rds"  ,emit:g51_14_rdsFile00_g51_17 

label 'scrna_seurat'

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

label 'scrna_seurat'


shell:

varFeatures = params.scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction.varFeatures

selmethod = params.scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction.selmethod

Batch_Effect_Correction = params.scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction.Batch_Effect_Correction

WNN = params.scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction.WNN

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
    $CPU  = 1
    $MEMORY = 140
}
//* platform
//* platform
//* autofill

process scRNA_Analysis_Module_Clustering_and_Find_Markers {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /Final_Report.html$/) "Final_Report/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /Final_Analysis.rds$/) "scViewer/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tsv$/) "ClusterMarkers/$filename"}
input:
 path seurat_object

output:
 path "Final_Report.html"  ,emit:g51_19_outputHTML00 
 path "Final_Analysis.rds"  ,emit:g51_19_rdsFile10_g51_22 
 path "*.tsv"  ,emit:g51_19_outFileTSV22 

label 'scrna_seurat'

shell:
minRes = params.scRNA_Analysis_Module_Clustering_and_Find_Markers.minRes
maxRes = params.scRNA_Analysis_Module_Clustering_and_Find_Markers.maxRes
npcs = params.scRNA_Analysis_Module_Clustering_and_Find_Markers.npcs
runCellFindR = params.scRNA_Analysis_Module_Clustering_and_Find_Markers.runCellFindR
findClusterforallResolution = params.scRNA_Analysis_Module_Clustering_and_Find_Markers.findClusterforallResolution

'''
#!/usr/bin/env perl

my $script = <<'EOF';

---
title: "Single Cell RNA-Seq Clustering Report"
author: "Zhaorong Li"
output: 
  html_document:
    toc: true
    toc_float:
      toc_collapsed: true
      toc_depth: 3
    number_sections: true
    fig_caption: yes
    theme: cerulean
    code_folding: hide
editor_options: 
  chunk_output_type: console
params:
  
  SamplePath: ""
  npcs: 25
  minRes: 0.1
  maxRes: 2.0
  filtered_output: ""
  algorithm: 2
---


```{r setup, include=FALSE}

suppressPackageStartupMessages({
library(dplyr)
library(Seurat)
library(ggplot2)
#if(!require(remotes)) install.packages("remotes",repos = "http://cran.us.r-project.org")
#if(!require(data.table)) install.packages("remotes",repos = "http://cran.us.r-project.org")
#if(!require(DT)) install.packages("DT",repos = "http://cran.us.r-project.org")
#if(!require(clustree)) install.packages("clustree",repos = "http://cran.us.r-project.org")

#remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
library(DoubletFinder)
library(data.table)
library(DT)
library(clustree)

is_cluster <- function(tenx, thresh_genes = 10, thresh_val = log(2), pval = 1e-4){
  val = 0 # groups that does not satisfy threshold genes
  counter = 0 # groups that satisfy threshold genes 
  # loop through the identitiy
  matrix_output <- data.frame(row.names = row.names(tenx))
  tenx=NormalizeData(tenx,assay="RNA")
  for (j in sort(unique(Idents(tenx)))){
    if (sum(tenx@active.ident == j) < 5){
      return(FALSE)
    }
    markers <- FindMarkers(tenx, ident.1 = j, min.pct = 0.25,assay="RNA",logfc.threshold=log(2))
    markers <- markers[markers$p_val_adj < pval,]
    #find if the 10th biggest is less than log2, sum 
    print(sort(markers$avg_log2FC, decreasing = TRUE)[thresh_genes])
    # if less than 10 significant genes
    
    if (length((markers$avg_log2FC)) < 10){
      val <- val + 1
    } else if (sort(markers$avg_log2FC, decreasing = TRUE)[thresh_genes] < thresh_val){
      #print(val)
      val <- val + 1
    } else{
      counter = counter + 1
    }
    if (val > 1){
      return(FALSE)
    }
  }
  if (val > 1){
    return(FALSE)
  }
  else{
    return(TRUE)
  }
}

# finds resolution that satisfy
# input: tenx object
# initial resolution of starting clustering
# how much to increment up 
# threshold of genes
# value of the threshold 
find_res <- function(tenx, initial_res = params$minRes, jump = 0.1, thresh_genes = 10, thresh_val = log(2),...) {
  RES_POST <- initial_res # keeping
  RES_IT <- initial_res # iterative
  
  while(TRUE){
    print(paste(RES_IT, sep = ' '))
    tenx <- FindClusters(tenx, resolution = RES_IT,...)
    Idents(tenx)="seurat_clusters"
    # also check if theres only 1 cluster/ then can go up higher es
    # Find number of clusters
    length_group <- length(unique(Idents(tenx)))
    # if only one group then need to look deeper
    if (length_group == 1){
      print(paste(RES_IT, sep = ' '))
      # still not groups at 0.7 res stop and just step as 1
      if (RES_IT == params$maxRes){
        print(paste(RES_IT, sep = ' '))
        break
      }
    } else{
      testing <- is_cluster(tenx)
      if (testing == FALSE){ # if not real group
        print(paste(RES_IT, sep = ' '))
        RES_IT <- RES_IT - jump
        RES_POST <- RES_IT
        print(RES_POST)
        break
      } else{ # valid groups
        RES_POST <- RES_IT
        print(paste(RES_IT, sep = ' '))
      }
    }
    RES_IT <- RES_IT + jump
  }
  # if there is only 1 group, return 0,
  return(RES_POST)
}


})
library(cluster)
Data=readRDS(params$SamplePath)
```

# Sample Statistics {.tabset}

In this section violin plots will show the number of genes, UMIs, mitochondrial percentages and ribosomal percentages of cells after filtering.

If there is more than 1 sample in the dataset, the violin plots will group the cells by samples. This is a very good way to compare quality of cells.

## Number of genes per cell

```{r Number of genes per cell by sample, fig.width=10,fig.height =10}

VlnPlot(Data,"nFeature_RNA",group.by="sample",pt.size = 0)+NoLegend()

```

## Number of UMIs per cell

```{r Number of UMIs per cell by sample, fig.width=10,fig.height =10}

VlnPlot(Data,"nCount_RNA",group.by="sample",pt.size = 0)+NoLegend()

```

## mitochondrial percentage per cell

```{r mitochondrial percentage per cell by sample, fig.width=10,fig.height =10}

VlnPlot(Data,"percent.mt",group.by="sample",pt.size = 0)+NoLegend()

```

## ribosomal percentage per cell

```{r ribosomal percentage per cell by sample, fig.width=10,fig.height =10}

VlnPlot(Data,"percent.ribo",group.by="sample",pt.size = 0)+NoLegend()

```

# Visualize PCA results

The results of principle component analysis give more insight than people usually realize. For example, the genes that contribute the most to the top principle components can help people to do sanity check of the data: ideally these genes will match the gene markers of the cell sub-populations in the data. This means that the cell heterogeneity is being captured.

```{r Dimension Reduction heatmap, fig.width=10,fig.height =10}

DimHeatmap(Data,dims = 1:10,nfeatures = 9,balanced = T,cells = 500,ncol = 3)

```

The genes that contribute the most to the top principle components can be visualized using Dimention Reduction heatmap shown above.


Another good way to visualize the PCA results is elbow plot, which plot the standard deviation of cells on each principle components. 

```{r Elbow Plot, fig.width=10,fig.height =10}

ElbowpotData=data.frame(stdev=Data@reductions[[ifelse("harmony"%in%names(Data@reductions),yes="harmony",no="pca")]]@stdev,PCs=seq(1,length(Data@reductions[[ifelse("harmony"%in%names(Data@reductions),yes="harmony",no="pca")]]@stdev)))

pc.changes=(diff((ElbowpotData$stdev)))*(-1)
pc.changes.raw=pc.changes
pc.changes=which(pc.changes>=mean(pc.changes.raw[pc.changes.raw>0]))
ggplot(ElbowpotData,aes(x=PCs,y=stdev,label=PCs))+geom_point()+theme_bw()+geom_vline(xintercept = max(pc.changes)+2,color="darkred")+geom_vline(xintercept = params$npcs,color="green")


```

# Visualize the UMAP and TSNE {.tabset}

Before doing any clustering, let us first use tSNE and UMAP to see if the batch effects between samples are removed.

If you only have one sample in this analysis, then please just enjoy the beautiful tSNE and UMAP. 

## tSNE {.tabset}

```{r tSNE, fig.width=10,fig.height =10,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}

Data=RunTSNE(Data,dims = 1:ifelse(params$npcs==0,yes=max(pc.changes)+1,no=params$npcs)
,reduction = ifelse("harmony"%in%names(Data@reductions),
                   yes="harmony",
                   no="pca"))
```

### tSNE colored by sample

```{r tSNE colored by sample, fig.width=10,fig.height =10,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}

if (length(unique(Data$sample))==1){
DimPlot(Data,reduction = "tsne")+NoLegend()
} else {
  DimPlot(Data,reduction = "tsne",group.by="sample")+NoLegend()

}
```

### tSNE colored by sample without batch effect correction

```{r tSNE without batch effect correction, fig.width=10,fig.height =10,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}

  Temp=Data
  Temp=RunTSNE(Temp,dims = 1:ifelse(params$npcs==0,yes=max(pc.changes)+1,no=params$npcs),reduction = "pca")
  DimPlot(Temp,reduction = "tsne",group.by="sample")+NoLegend()

```

### tSNE split by sample

```{r tSNE split by sample, fig.width=10,fig.height =10,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}

if (length(unique(Data$sample))==1){
DimPlot(Data,reduction = "tsne")+NoLegend()
} else {
  DimPlot(Data,reduction = "tsne",split.by="sample",ncol=2)+NoLegend()
}
```

### tSNE split by sample without batch effect correction

```{r tSNE split by sample without batch effect correction, fig.width=10,fig.height =10,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}

  DimPlot(Temp,reduction = "tsne",split.by="sample",ncol=2)+NoLegend()
  rm(Temp)

```

## UMAP {.tabset}

```{r umap, fig.width=10,fig.height =10,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}

Data=RunUMAP(Data,dims = 1:ifelse(params$npcs==0,yes=max(pc.changes)+1,no=params$npcs),reduction = ifelse("harmony"%in%names(Data@reductions),
                   yes="harmony",
                   no="pca"))
```

### UMAP colored by sample

```{r UMAP colored by sample, fig.width=10,fig.height =10,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}

if (length(unique(Data$sample))==1){
DimPlot(Data,reduction = "umap")+NoLegend()
} else {
  DimPlot(Data,reduction = "umap",group.by="sample")+NoLegend()

}
```

### UMAP colored by sample without batch effect correction

```{r UMAP without batch effect correction, fig.width=10,fig.height =10,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}
  Temp=Data
  Temp=RunUMAP(Temp,dims = 1:(ifelse(params$npcs==0,yes=max(pc.changes)+1,no=params$npcs)),reduction = "pca")
  DimPlot(Temp,reduction = "umap",group.by="sample")+NoLegend()
```

### UMAP split by sample

```{r UMAP split by sample, fig.width=10,fig.height =10,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}

if (length(unique(Data$sample))==1){
DimPlot(Data,reduction = "umap")+NoLegend()
} else {
  DimPlot(Data,reduction = "umap",split.by="sample",ncol=2)+NoLegend()
}
```

### UMAP split by sample without batch effect correction

```{r UMAP split by sample without batch effect correction, fig.width=10,fig.height =10,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}

  DimPlot(Temp,reduction = "umap",split.by="sample",ncol=2)+NoLegend()
  rm(Temp)
  
```



## Note on tSNE and UMAP.

Both visualizations are widely used. Feel free to choose the one you like.

One thing to take in mind is that it is faster to generate UMAP reduction than to generate tSNE reduction.

# Clustering

## Building the nearest neighborhood graph

In order to cluster the cells, the shared nearest neighborhood graph of cells are constructed using the top principle components (default is 25).


```{r Build snn graph, error=FALSE, fig.height =10, fig.width=10, message=FALSE, warning=FALSE,results = FALSE}
if ("wpca"%in%names(Data@reductions)) {
	ElbowpotData=data.frame(stdev=Data@reductions[[ifelse("wharmony"%in%names(Data@reductions),yes="wharmony",no="wpca")]]@stdev,PCs=seq(1,length(Data@reductions[[ifelse("wharmony"%in%names(Data@reductions),yes="wharmony",no="wpca")]]@stdev)))
	wpc.changes=(diff((ElbowpotData$stdev)))*(-1)
	wpc.changes.raw=wpc.changes
	wpc.changes=which(wpc.changes>=mean(wpc.changes.raw[wpc.changes.raw>0]))
	Data <- FindMultiModalNeighbors(Data,
	reduction.list = list(ifelse("harmony"%in%names(Data@reductions),yes="harmony",no="pca"), 
  ifelse("wharmony"%in%names(Data@reductions),yes="wharmony",no="wpca")
  ), 
  dims.list = list(1:max(pc.changes+1), 1:max(wpc.changes+1)), modality.weight.name = "multi_modal_weight"
)
Data <- RunUMAP(Data, nn.name = "weighted.nn")
Data <- RunTSNE(Data, nn.name = "weighted.nn")

} else {
	Data=FindNeighbors(Data,dims = 1:ifelse(params$npcs==0,yes=max(pc.changes)+1,no=params$npcs),reduction = ifelse("harmony"%in%names(Data@reductions),
                   yes="harmony",
                   no="pca"))

}



```


## Clustering {.tabset}

And then Graph Based Community Detection Algorithm is used to cluster the cells.

In order to select the best clustering resolution, the sc3 stability index is calculate for each resolution. The resolution with the highest mean sc3 stability index (marked by red line in the figure below).

```{r Clustering, fig.width=10,fig.height =10,error=FALSE,warning=FALSE,message=FALSE,results=FALSE,echo=FALSE}

if ("wsnn"%in%names(Data@graphs)) {
	Data <- FindClusters(Data, graph.name = "wsnn",algorithm = params$algorithm,resolution = seq(params$maxRes,params$minRes,-0.1))

} else {
	Data=FindClusters(Data,algorithm = params$algorithm,resolution = seq(params$maxRes,params$minRes,-0.1))

}


if ("wsnn"%in%names(Data@graphs)) {
names=paste0("wsnn","_res.")

} else {
names=paste0(DefaultAssay(Data),"_snn_res.")

}

SC3_Stability=clustree(Data,prefix = names)
SC3_Stability.results=SC3_Stability$data
SC3_Stability.results=SC3_Stability.results[,c(names,"sc3_stability")]
colnames(SC3_Stability.results)[1]="resolution"
SC3_Stability.results.mean=aggregate(sc3_stability~resolution,SC3_Stability.results,mean)
colnames(SC3_Stability.results.mean)[2]="sc3_stability_mean"
if ("wsnn"%in%names(Data@graphs)) {
Idents(Data)=paste0("wsnn","_res.",max(as.numeric(as.character(SC3_Stability.results.mean$resolution))[SC3_Stability.results.mean$sc3_stability_mean==max(SC3_Stability.results.mean$sc3_stability_mean)]))

} else {
Idents(Data)=paste0(DefaultAssay(Data),"_snn_res.",max(as.numeric(as.character(SC3_Stability.results.mean$resolution))[SC3_Stability.results.mean$sc3_stability_mean==max(SC3_Stability.results.mean$sc3_stability_mean)]))
}
Data$seurat_clusters=Idents(Data)

Cluster.distribution=data.frame(table(Data$seurat_clusters,Data$sample))
colnames(Cluster.distribution)=c("Cluster","Sample","CellNumber")



```

### Cluster stability assessment

```{r Cluster stability assessment, fig.width=10,fig.height=10,error=FALSE,warning=FALSE,message=FALSE,results=FALSE,echo=FALSE}
ggplot(SC3_Stability.results,aes(x=resolution,y=sc3_stability))+geom_boxplot()+geom_line(data=SC3_Stability.results.mean,aes(x=resolution,y=sc3_stability_mean,group=1))+geom_vline(xintercept = SC3_Stability.results.mean$resolution[SC3_Stability.results.mean$sc3_stability_mean==max(SC3_Stability.results.mean$sc3_stability_mean)],color='red')
```

### Clustree assessment

This figure shows how the cells are assigned as the resolution changes. The color of the arrow shows the amount of cells going into the cluster in the next level and the direction of the arrow shows the identity of cluster that the cells are going to.

As the resolution increases, the arrows will start to appear "messy". This means that the clustering algorithm is having trouble assigning cells.

```{r Clustree assessment, fig.width=10,fig.height=10,error=FALSE,warning=FALSE,message=FALSE,results=FALSE,echo=FALSE}
SC3_Stability
```

### CellFindR assessment

This is an algorithm developed by Kevin Yu in the Tward lab at UCSF. The algorithm tries to find the optimal clustering results from the single cell RNA-Seq data.

```{r CellFindR, fig.width=10,fig.height=10,error=FALSE,warning=FALSE,message=FALSE,results=FALSE,echo=FALSE}

if (as.logical("!{runCellFindR}")) {
	if ("wsnn"%in%names(Data@graphs)) {
	CellFindR.resolution=find_res(Data,graph.name="wsnn")

	CellFindR=FindClusters(Data,resolution=CellFindR.resolution,graph.name="wsnn")

} else {
	CellFindR.resolution=find_res(Data)

	CellFindR=FindClusters(Data,resolution=CellFindR.resolution)

}

Data$CellFindR.clustering.results=CellFindR$seurat_clusters

}


```








### Sample Distribution over Cluster
```{r Sample distribution, fig.width=10,fig.height=10,error=FALSE,warning=FALSE,message=FALSE,results=FALSE,echo=FALSE}



ggplot(Cluster.distribution,aes(y=Cluster,x=CellNumber,fill=Sample))+geom_bar(stat="identity",position="fill")

```

### Cluster Distribution over Sample
```{r Cluster distribution, fig.width=10,fig.height=10,error=FALSE,warning=FALSE,message=FALSE,results=FALSE,echo=FALSE}

Cluster.distribution=data.frame(table(Data$seurat_clusters,Data$sample))
colnames(Cluster.distribution)=c("Cluster","Sample","CellNumber")

ggplot(Cluster.distribution,aes(y=Sample,x=CellNumber,fill=Cluster))+geom_bar(stat="identity",position="fill")

```

### tSNE Visualization of the Cluster

```{r tSNE Visualization of the Cluster, fig.width=10,fig.height=10,error=FALSE,warning=FALSE,message=FALSE,results=FALSE,echo=FALSE}
DimPlot(Data,reduction = "tsne",label = T)

```

### tSNE distributed by sample

```{r tSNE distributed by sample, fig.width=10,fig.height =10,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}

if (length(unique(Data$sample))==1){
DimPlot(Data,reduction = "tsne")+NoLegend()
} else {
  DimPlot(Data,reduction = "tsne",split.by="sample",ncol=2)+NoLegend()
}
```

### UMAP Visualization of the Cluster

```{r UMAP Visualization of the Cluster, fig.width=10,fig.height =10,error=FALSE,warning=FALSE,message=FALSE,results=FALSE,echo=FALSE}
DimPlot(Data,reduction = "umap",label = T)

```

### UMAP distributed by sample

```{r UMAP distributed by sample, fig.width=10,fig.height =10,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}

if (length(unique(Data$sample))==1){
DimPlot(Data,reduction = "umap")+NoLegend()
} else {
  DimPlot(Data,reduction = "umap",split.by="sample",ncol=2)+NoLegend()
}
```



## Cluster Markers {.tabset}

Use differential expression analysis to find Gene markers for each cluster. These gene markers are very helpful in identifying Cell types.

```{r Find Cluster markers, fig.width=15,fig.height=20,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}

Data=NormalizeData(Data,assay = "RNA")


Data.markers=FindAllMarkers(Data,only.pos = T,assay = "RNA")
write.table(Data.markers,"Cluster.Markers.tsv",quote=F,sep="\t")
if (as.logical("!{findClusterforallResolution}")) {
	for (i in colnames(Data@meta.data)[grepl("snn_res.",colnames(Data@meta.data))]) {
	temp=Data
	Idents(temp)=i
	temp=FindAllMarkers(temp,only.pos = T)
	write.table(temp,paste0(i,".Cluster.Markers.tsv"),quote=F,sep="\t")

}
}



if ("cluster" %in% colnames(Data.markers)){
top10 = Data.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
} else {
	top10=NULL
}
saveRDS(Data,"Final_Analysis.rds")

```

### Top gene markers for clusters

```{r Top gene markers for clusters, fig.width=15,fig.height=20,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}
#if (!is.null(top10)){
DT::datatable(top10,filter = "top",options = list(autoWidth = TRUE))
#}
```

### Heatmap of top gene markers for clusters

```{r heatmap of Top gene markers for clusters, fig.width=15,fig.height=20,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}
if (!is.null(top10)){
Vis=ScaleData(Data,features = top10$gene,assay = "RNA")
DoHeatmap(Vis, features = top10$gene,assay = "RNA") + NoLegend()
}
```


## Cluster Results Quality Control {.tabset}

In this section the number of genes, UMIs, mitochondrial percentages and ribosomal percentages will be plotted by cluster. This step is to check whether the clustering results are significantly biased by the sequencing depth, sequencing coverage and cell viability.

However, researches have shown (insert reference here later) that number of genes, UMIs and mitochondrial contents will vary between cell types and sub-populations. 

### Violin Plots {.tabset}

#### number of genes per cluster

```{r number of genes per cluster v, fig.width=10,fig.height =10}

VlnPlot(Data,"nFeature_RNA",group.by="seurat_clusters",pt.size = 0)+NoLegend()

```

#### number of UMIs per cluster

```{r number of UMIs per cluster v, fig.width=10,fig.height =10}

VlnPlot(Data,"nCount_RNA",group.by="seurat_clusters",pt.size = 0)+NoLegend()

```

#### mitochondrial percentages per cluster

```{r mitochondrial percentages per cluster v, fig.width=10,fig.height =10}

VlnPlot(Data,"percent.mt",group.by="seurat_clusters",pt.size = 0)+NoLegend()

```

#### ribosomal percentages per cluster

```{r ribosomal percentages per cluster v, fig.width=10,fig.height =10}

VlnPlot(Data,"percent.ribo",group.by="seurat_clusters",pt.size = 0)+NoLegend()

```

### Ridge Plots {.tabset}

#### number of genes per cluster

```{r number of genes per cluster r, fig.width=10,fig.height =10}

RidgePlot(Data,"nFeature_RNA",group.by="seurat_clusters")+NoLegend()

```

#### number of UMIs per cluster

```{r number of UMIs per cluster r, fig.width=10,fig.height =10}

RidgePlot(Data,"nCount_RNA",group.by="seurat_clusters")+NoLegend()

```

#### mitochondrial percentages per cluster

```{r mitochondrial percentages per cluster r, fig.width=10,fig.height =10}

RidgePlot(Data,"percent.mt",group.by="seurat_clusters")+NoLegend()

```

#### ribosomal percentages per cluster

```{r ribosomal percentages per cluster r, fig.width=10,fig.height =10}

RidgePlot(Data,"percent.ribo",group.by="seurat_clusters")+NoLegend()

```




EOF

open OUT, "> ClusterandMarker.rmd";
print OUT $script;
close OUT;

runCommand("Rscript -e 'rmarkdown::render(\\"ClusterandMarker.rmd\\",\\"html_document\\", output_file = \\"Final_Report.html\\",
params = list(SamplePath=\\"!{seurat_object}\\",minRes=as.numeric(!{minRes}),
maxRes=as.numeric(!{maxRes}),npcs=as.numeric(!{npcs})))'");

sub runCommand {
            my ($com) = @_;
            my $error = system($com);
            if   ($error) { die "Command failed: $error $com\\n"; }
            else          { print "Command successful: $com\\n"; }
          }

'''



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
 path "*.h5ad"  ,emit:g51_22_h5_file00 

label 'scrna_seurat'

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

label 'scrna_seurat'



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


workflow {


mkfastq_prep()
g_4_bcl00_g_0 = mkfastq_prep.out.g_4_bcl00_g_0
g_4_bcl_directory13_g_13 = mkfastq_prep.out.g_4_bcl_directory13_g_13
(g_4_bcl_directory16_g_20,g_4_bcl_directory15_g_55) = [g_4_bcl_directory13_g_13,g_4_bcl_directory13_g_13]


mkfastq(g_4_bcl00_g_0.flatten())
g_0_reads00_g_33 = mkfastq.out.g_0_reads00_g_33
(g_0_reads00_g_13,g_0_reads00_g_20) = [g_0_reads00_g_33,g_0_reads00_g_33]
g_0_outputDir11 = mkfastq.out.g_0_outputDir11
g_0_outputHTML22 = mkfastq.out.g_0_outputHTML22


cellranger_multi_prep()
g_5_csvFile05_g_20 = cellranger_multi_prep.out.g_5_csvFile05_g_20


cellranger_fastq_prep(g_23_0_g_22,g_24_1_g_22)
g_22_reads00_g_25 = cellranger_fastq_prep.out.g_22_reads00_g_25


cellranger_fastq_collect(g_22_reads00_g_25.collect())
g_25_bcl_directory05_g_13 = cellranger_fastq_collect.out.g_25_bcl_directory05_g_13
(g_25_bcl_directory07_g_20,g_25_bcl_directory07_g_55) = [g_25_bcl_directory05_g_13,g_25_bcl_directory05_g_13]
g_25_reads16_g_13 = cellranger_fastq_collect.out.g_25_reads16_g_13
(g_25_reads18_g_20,g_25_reads10_g_55) = [g_25_reads16_g_13,g_25_reads16_g_13]


cellranger_count_prep()
g_27_csvFile04_g_13 = cellranger_count_prep.out.g_27_csvFile04_g_13


if (!((params.run_FastQC && (params.run_FastQC == "yes")))){
g_0_reads00_g_33.set{g_33_reads00_g_34}
g_33_mate10_g_35 = Channel.empty()
} else {

flatten_cellranger_reads(g_0_reads00_g_33.collect())
g_33_reads00_g_34 = flatten_cellranger_reads.out.g_33_reads00_g_34
g_33_mate10_g_35 = flatten_cellranger_reads.out.g_33_mate10_g_35
}


file_to_set_conversion_for_reads(g_33_reads00_g_34.flatten())
g_34_reads01_g_35 = file_to_set_conversion_for_reads.out.g_34_reads01_g_35


if (!((params.run_FastQC && (params.run_FastQC == "yes")))){
g_34_reads01_g_35.set{g_35_reads11}
g_35_FastQCout04_g_52 = Channel.empty()
} else {

FastQC_after_mkfastq(g_33_mate10_g_35,g_34_reads01_g_35)
g_35_FastQCout04_g_52 = FastQC_after_mkfastq.out.g_35_FastQCout04_g_52
g_35_reads11 = FastQC_after_mkfastq.out.g_35_reads11
}


MultiQC(g_35_FastQCout04_g_52.flatten().toList())
g_52_outputHTML00 = MultiQC.out.g_52_outputHTML00


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
(g_18_reference01_g_20,g_18_reference01_g_55) = [g_18_reference01_g_13,g_18_reference01_g_13]

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



scRNA_Analysis_Module_Load_Data_h5(g_61_h5_file00_g51_0,g_43_1_g51_0)
g51_0_rdsFile00_g51_14 = scRNA_Analysis_Module_Load_Data_h5.out.g51_0_rdsFile00_g51_14
g51_0_outputFileHTML11 = scRNA_Analysis_Module_Load_Data_h5.out.g51_0_outputFileHTML11


scRNA_Analysis_Module_Merge_Seurat_Objects(g51_0_rdsFile00_g51_14.collect())
g51_14_rdsFile00_g51_17 = scRNA_Analysis_Module_Merge_Seurat_Objects.out.g51_14_rdsFile00_g51_17


scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction(g51_14_rdsFile00_g51_17)
g51_17_rdsFile00_g51_19 = scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction.out.g51_17_rdsFile00_g51_19


scRNA_Analysis_Module_Clustering_and_Find_Markers(g51_17_rdsFile00_g51_19)
g51_19_outputHTML00 = scRNA_Analysis_Module_Clustering_and_Find_Markers.out.g51_19_outputHTML00
g51_19_rdsFile10_g51_22 = scRNA_Analysis_Module_Clustering_and_Find_Markers.out.g51_19_rdsFile10_g51_22
(g51_19_rdsFile10_g51_25,g51_19_rdsFile10_g51_30) = [g51_19_rdsFile10_g51_22,g51_19_rdsFile10_g51_22]
g51_19_outFileTSV22 = scRNA_Analysis_Module_Clustering_and_Find_Markers.out.g51_19_outFileTSV22


scRNA_Analysis_Module_Create_h5ad(g51_19_rdsFile10_g51_22)
g51_22_h5_file00 = scRNA_Analysis_Module_Create_h5ad.out.g51_22_h5_file00


scRNA_Analysis_Module_seurat_to_sce(g51_19_rdsFile10_g51_25)
g51_25_rdsFile00_g51_27 = scRNA_Analysis_Module_seurat_to_sce.out.g51_25_rdsFile00_g51_27


scRNA_Analysis_Module_launch_isee(g51_25_rdsFile00_g51_27)
g51_27_outputFileOut00 = scRNA_Analysis_Module_launch_isee.out.g51_27_outputFileOut00


scRNA_Analysis_Module_SCEtoLOOM(g51_19_rdsFile10_g51_30)
g51_30_outputFileOut00 = scRNA_Analysis_Module_SCEtoLOOM.out.g51_30_outputFileOut00


}

workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
