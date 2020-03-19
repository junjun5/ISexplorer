# ISexplorer

## Introduction
mapping-based tool for the identification of IS positions in bacterial genomes

## Installation
ISexplorer is currelntly available only on Linux. Users can download ISexplorer and also need to download required dependencies programs below.
ISexplorer is shell program, that's why you can use it without compiling.

## required dependencies
- GNU Awk (3.1.7)
- Bedtools (2.28.0)
https://bedtools.readthedocs.io/en/latest/
- BWA (0.7.17)
http://bio-bwa.sourceforge.net/
- Samtools (1.2)
http://www.htslib.org/

## Usage
In order to run ISexplorer with default settings, execute the following command:

```
isexplorer --is [Index_target_IS] --reference [Index_target_IS] --read1 [paired_end_read1] --read2 [paired_end_read2] --islength [IS_length_of_target_IS] --insert [insert_size] --readlength [length_of_paired_end_read] --seqcoverage [sequence_coverage_of_paired_end_read]
```

Options:
- ```--is (required)```  
    path to bwa index of query IS
- ```--reference (required)```  
    path to bwa index of reference genome
- ```--read1 (required)```  
    paired-end read1
- ```--read2 (required)```  
    paired-end read2
- ```--islength (required)```  
    length of target IS
- ```--insert (required)```  
    mean insert size of paired-end read
- ```--readlength (required)```  
    length of paired end read
- ```--seqcoverage (required)```  
    sequence coverage of paired end read
- ```--besthit```  
    rate of BestHit Block threshold (default : 0.4)
- ```--multihit```  
    rate of Multihit Block threshold (default : 0.8)
- ```--uniqhit```  
    rage of UniqueHit Block threshold (default : 0.1)
- ```--threads```  
    threads num (default : 1)
- ```--outputdir```  
    path to output directory (default : default)
- ```--no_clean```  
    ISexplorer does not remove intermediate files
- ```--help```  
    print detailed descriptions of command line arguments

  

## Test
We prepared test. Users and developers can download the data, and check whether ISexplorer works properly by the following commands:
