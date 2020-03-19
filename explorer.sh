#!/bin/sh

# sh IShunter.sh
# --is              :Index_IS
# --reference       :Index_Reference
# --Read1           :Read1
# --Read2           :Read2
# --islength        :IS_Length
# --insert          :Insert_Size
# --readlength      :Read_Length
# --seqcoverage     :Seq_Coverage
# -o --outputdir    :Output_Dir
# -t --threads      :Threads
# -h --help         :help
flag_1=false;flag_2=false;flag_3=false;flag_4=false;flag_5=false;flag_6=false;flag_7=false;flag_8=false

usage_exit() {
    echo "USAGE:IShunter.sh [--is Index_IS] [--reference  Index_Reference] [--read1 Read1] [--read2 Read2] [--islength IS_Length] [--insert Insert_Size] [--readlength Read_Length] [--seqcoverage Seq_Coverage] [--outputdir Output_Dir]" \
        1>&2
    echo "Options:"
    echo "      --is: path to query IS's bwa index (required)"
    echo "      --reference: path to reference genome's bwa index (required)"
    echo "      --read1: paired end read 1 (required)"
    echo "      --read2: paired end read 2 (required)"
    echo "      --islength: query IS length (required)"
    echo "      --insert: paired end read's mean insert size (required)"
    echo "      --readlength: paired end read's read length (required)"
    echo "      --seqcoverage: paired end read seqence coverage (required)"
    echo "      --besthit: rate of BestHit Block threshold (default: 0.4)"
    echo "      --multihit: rate of MultiHit Block threshold (default: 0.8)"
    echo "      --uniqhit: rate of UniqueHit BLock threshold (default: 0.1)"
    echo "      --threads: threads num (default: 1)"
    echo "      --outputdir: path to output directory (default: Default)"
    echo "      --no_clean: ISexplorer does not remove intermediate files"
    echo "      --help: print detailed descriptions of command line arguments"
    exit 1
}
Output_Dir="Default"
Threads=1
BestHit=0.4
MultiHit=0.8
UniqHit=0.1
No_Clean=0

while getopts ":-:" OPT
do
    optarg="${!OPTIND}"
    [[ "$OPT" = - ]] && OPT="-$OPTARG"
    case "-$OPT" in
        --is) Index_IS="$optarg";flag_1=true;echo $Index_IS;shift;;
        --reference) Index_Reference="$optarg";flag_2=true;shift;;
        --read1) Read1="$optarg";echo $Read1;flag_3=true;shift;;
        --read2) Read2="$optarg";flag_4=true;shift;;
        --islength) IS_length="$optarg";flag_5=true;shift;;
        --insert) Insert_Size="$optarg";flag_6=true;shift;;
        --readlength) Read_Length="$optarg";flag_7=true;shift;;
        --seqcoverage) Seq_Coverage="$optarg";flag_8=true;shift;;
        --besthit) BestHit="$optarg";shift;;
        --multihit) MultiHit="$optarg";shift;;
        --uniqhit) UniqHit="$optarg";shift;;
        --threads) Threads="$optarg";shift;;
        --outputdir) Output_Dir="$optarg";shift;;
        --no_clean) No_Clean=1;shift;;
        --help) usage_exit ;;
    esac
done
if  ! $flag_1 || ! $flag_2 || ! $flag_3 || ! $flag_4 || ! $flag_5 || ! $flag_6 || ! $flag_7 || ! $flag_8  ; then
    echo "--is,--reference,--read1,--read2,--islength,--insert,--readlength,--seqcoverage must be included"
    usage_exit
    exit 1
fi

echo ${Output_Dir}/tmp

if [ ! -d ${Output_Dir}/tmp ];then
    mkdir -p ${Output_Dir}/tmp
fi

threshold_def=$(((Insert_Size-Read_Length)*Seq_Coverage/Read_Length / 2))
block_gap=$(( Insert_Size - Read_Length * 3))
echo sh IShunter.sh,$Index_IS,$Index_Reference,$Read1,$Read2,$IS_length,$Insert_Size,$Read_Length,$Seq_Coverage,$threshold_def | tee ${Output_Dir}/Logfile.txt
echo -e "[Index_IS]\t"${Index_IS}  | tee -a ${Output_Dir}/Logfile.txt
echo -e "[Index_Reference]\t"${Index_Reference}  | tee -a ${Output_Dir}/Logfile.txt
echo -e "[Read1]\t"${Read1}  | tee -a ${Output_Dir}/Logfile.txt
echo -e "[Read2]\t"${Read2}  | tee -a ${Output_Dir}/Logfile.txt
echo -e "[IS_Length]\t"${IS_length}  | tee -a ${Output_Dir}/Logfile.txt
echo -e "[Insert_Size]\t"${Insert_Size}  | tee -a ${Output_Dir}/Logfile.txt
echo -e "[Read_Length]\t"${Read_Length}  | tee -a ${Output_Dir}/Logfile.txt
echo -e "[Seq_Coverage]\t"${Seq_Coverage}  | tee -a ${Output_Dir}/Logfile.txt
echo -e "[BestHit_Block_rate]\t"${BestHit} | tee -a ${Output_Dir}/Logfile.txt
echo -e "[MultiHit_Block_rate]\t"${MultiHit} | tee -a ${Output_Dir}/Logfile.txt
echo -e "[UniqHit_Block_rate]\t"${UniqHit} | tee -a ${Output_Dir}/Logfile.txt
echo -e "[Threads]\t"${Threads} | tee -a ${Output_Dir}/Logfile.txt

if [ ! $block_gap -gt 0 ]; then
    block_gap=0;
fi

echo -e "[block_Gap_Length]\t"${block_gap}  | tee -a ${Output_Dir}/Logfile.txt
BestHit_threshold=`echo ${threshold_def} | awk -v besthit=${BestHit} '{print $1*besthit}'` ;BestHit_threshold=${BestHit_threshold%.*}
MultiHit_threshold=`echo ${threshold_def} | awk -v multihit=${MultiHit} '{print $1*multihit}'` ;MultiHit_threshold=${MultiHit_threshold%.*}
UniqueHit_threshold=`echo ${threshold_def} | awk -v uniqhit=${UniqHit} '{print $1*uniqhit}'` ;UniqueHit_threshold=${UniqueHit_threshold%.*}

if [ ${UniqueHit_threshold} -eq 0 ]; then
    UniqueHit_threshold=1;
fi

echo -e "[Threshold default]\t" ${threshold_def}  | tee -a ${Output_Dir}/Logfile.txt
echo -e "[BestHit_threshold]\t"${BestHit_threshold}  | tee -a ${Output_Dir}/Logfile.txt
echo -e "[MultiHit_threshold]\t"${MultiHit_threshold}  | tee -a ${Output_Dir}/Logfile.txt
echo -e "[UniqueHit_threshold]\t"${UniqueHit_threshold}  | tee -a ${Output_Dir}/Logfile.txt

# 検出したいISにpair readをマッピング
# Target IS mapped to paired-end reads
bwa mem -t ${Threads} ${Index_IS} ${Read1} ${Read2} 1> ${Output_Dir}/tmp/paired.sam 2>> ${Output_Dir}/Logfile.txt

echo sequence read mapped to target IS ..... done

# マッピング結果からISのright side、left sideをmapped、unmappedで別々に出力
# split mapping result into right side and left side
# left side read is mapped to IS's left side sequences
# right siede read is mapped to IS's right side sequences
samtools view -@ ${Threads} ${Output_Dir}/tmp/paired.sam -f 24 -F 4  > ${Output_Dir}/tmp/left_mapped_toIS.sam      # ISの左端にマップされたものを抽出
samtools view -@ ${Threads} ${Output_Dir}/tmp/paired.sam -f 36 -F 8  > ${Output_Dir}/tmp/left_unmapped_toIS.sam    # ISの左端にマップされたもののpair readを抽出
samtools view -@ ${Threads} ${Output_Dir}/tmp/paired.sam -f 8  -F 20 > ${Output_Dir}/tmp/right_mapped_toIS.sam     # ISの右端にマップされたものを抽出
samtools view -@ ${Threads} ${Output_Dir}/tmp/paired.sam -f 4  -F 40 > ${Output_Dir}/tmp/right_unmapped_toIS.sam   # ISの右端にマップされたもののpair readを抽出

if [ ${No_Clean} -eq 1 ]; then
    rm -f ${Output_Dir}/tmp/paired.sam;
fi

if [ ! -s ${Output_Dir}/tmp/left_mapped_toIS.sam -o ! -s ${Output_Dir}/tmp/right_mapped_toIS.sam ];then
    echo This family-subgroup-IS does not exist in this sequence genome.;
    if [ ${No_Clean} -eq 1 ]; then
        rm -f ${Output_Dir}/tmp/;
    fi
    exit
fi

# mappedされたマッピング結果からSoftclipped(S)、Hardclipped(H)を取り除く
# readが全長マッピングされたリードを取得
# Extract softclipped read and hardclipped read from mapping result

# 両末端がSoft clippedなreadのIDを抽出
cat ${Output_Dir}/tmp/left_mapped_toIS.sam  \
    | awk '{if($6 ~ /S/){num=match($6,"S");if(substr($6,num+1) ~ /S/) {print $1}}}' > ${Output_Dir}/tmp/left_mapped_bothS.id
cat ${Output_Dir}/tmp/right_mapped_toIS.sam  \
    | awk '{if($6 ~ /S/){num=match($6,"S");if(substr($6,num+1) ~ /S/) {print $1}}}' > ${Output_Dir}/tmp/right_mapped_bothS.id

# complete match to target IS
cat ${Output_Dir}/tmp/left_mapped_toIS.sam  \
    | awk '{if($6 !~ /S/){print $1} else{num=match($6,"S");if(substr($6,num+1) !~ /S/) {print $1}}}' > ${Output_Dir}/tmp/left_mapped_match.id
cat ${Output_Dir}/tmp/right_mapped_toIS.sam \
    | awk '{if($6 !~ /S/){print $1} else{num=match($6,"S");if(substr($6,num+1) !~ /S/) {print $1}}}' > ${Output_Dir}/tmp/right_mapped_match.id

# level3 : flanking read (min_length >= 60)
samtools view -@ ${Threads} ${Output_Dir}/tmp/paired.sam -F 4 \
    | awk '{if($6 ~ /S/){num=match($6,"S");if(substr($6,num+1) !~ /S/) {print $0}}}' \
    | awk '{cigar=$6;split(cigar,cigar_list,"S");if(cigar_list[1] !~ /M/ && int(cigar_list[1])>60) {target=substr($10,0,cigar_list[1]);score=substr($11,0,cigar_list[1]);printf "@%s%s:p\n",$1,NR;printf "%s\n+\n%s\n",target,score}}'\
    > ${Output_Dir}/tmp/left_unmapped_p.fq
samtools view -@ ${Threads} ${Output_Dir}/tmp/paired.sam -F 4 \
    | awk '{if($6 ~ /S/){num=match($6,"S");if(substr($6,num+1) !~ /S/) {print $0}}}' \
    | awk '{cigar=$6;split(cigar,cigar_list,"M");len_cigar=length(cigar_list);
    if(cigar_list[len_cigar] ~ /S/) {split(cigar_list[len_cigar],num_S,"S");
    if(int(num_S[1]) > 60){target=substr($10,length($10)-num_S[1]+1,num_S[1]);target_score=substr($11,length($11)-num_S[1]+1,num_S[1]);
        printf "@%s%s:p\n",$1,NR;gsub("A","t",target);gsub("T","a",target);gsub("C","g",target);gsub("G","c",target);
        target=toupper(target);split(target,target_list,//);for(i=length(target_list);i>0;i--){printf "%s",target_list[i]} printf "\n+\n";split(target_score,target_score_list,//);
        for(i=length(target_score_list);i>0;i--){printf "%s",target_score_list[i]}printf "\n"}}}'\
    > ${Output_Dir}/tmp/right_unmapped_p.fq

# unmapped readからpairが全長マッピングされたリード(level12)を取得
# get unmapped read the pair of which is mapped to reference with level1 or level2
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN{file1=Output_Dir"/left_unmapped_toIS.sam";file2=Output_Dir"/left_mapped_match.id";
    while((getline < file1 )>0) {while((getline var < file2 )>0) {if($1 == var) print $0} close(file2)}}' \
    > ${Output_Dir}/tmp/left_unmapped.sam
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN{file1=Output_Dir"/right_unmapped_toIS.sam";file2=Output_Dir"/right_mapped_match.id";while((getline < file1 )>0)
    {while((getline var < file2 )>0) {if($1 == var) print $0} close(file2)}}' \
    > ${Output_Dir}/tmp/right_unmapped.sam

# headerを付与
# get samfile's header
samtools view -H ${Output_Dir}/tmp/paired.sam > ${Output_Dir}/tmp/header_IS.sam
rm -f ${Output_Dir}/tmp/paired.sam
echo Removed sam file.

cat ${Output_Dir}/tmp/header_IS.sam ${Output_Dir}/tmp/left_unmapped.sam  > ${Output_Dir}/tmp/left_unmapped_head.sam
cat ${Output_Dir}/tmp/header_IS.sam ${Output_Dir}/tmp/right_unmapped.sam > ${Output_Dir}/tmp/right_unmapped_head.sam

# flanking readをfastqfileで作成
# From level1 and level2 flanking read, get fastqfile
samtools view -@ ${Threads} -Sb ${Output_Dir}/tmp/left_unmapped_head.sam -o ${Output_Dir}/tmp/left_unmapped.bam
samtools view -@ ${Threads} -Sb ${Output_Dir}/tmp/right_unmapped_head.sam -o ${Output_Dir}/tmp/right_unmapped.bam
bedtools bamtofastq -i ${Output_Dir}/tmp/left_unmapped.bam -fq ${Output_Dir}/tmp/left_unmapped_buf.fq 2>> ${Output_Dir}/Logfile.txt
bedtools bamtofastq -i ${Output_Dir}/tmp/right_unmapped.bam -fq ${Output_Dir}/tmp/right_unmapped_buf.fq 2>> ${Output_Dir}/Logfile.txt
cat ${Output_Dir}/tmp/left_unmapped_buf.fq  ${Output_Dir}/tmp/left_unmapped_p.fq > ${Output_Dir}/tmp/left_unmapped.fq
cat ${Output_Dir}/tmp/right_unmapped_buf.fq ${Output_Dir}/tmp/right_unmapped_p.fq > ${Output_Dir}/tmp/right_unmapped.fq

# flanking readをreference genomeにmapping
# level1 and level2 flanking read map to reference genome
bwa mem -t ${Threads} ${Index_Reference} ${Output_Dir}/tmp/left_unmapped.fq 1> ${Output_Dir}/tmp/left.sam 2>> ${Output_Dir}/Logfile.txt
bwa mem -t ${Threads} ${Index_Reference} ${Output_Dir}/tmp/right_unmapped.fq 1> ${Output_Dir}/tmp/right.sam 2>> ${Output_Dir}/Logfile.txt
# level1 and level2 flanking read map to reference genome with multimapping
bwa mem -a -t ${Threads} ${Index_Reference} ${Output_Dir}/tmp/left_unmapped.fq 1> ${Output_Dir}/tmp/left_multi.sam 2>> ${Output_Dir}/Logfile.txt
bwa mem -a -t ${Threads} ${Index_Reference} ${Output_Dir}/tmp/right_unmapped.fq 1> ${Output_Dir}/tmp/right_multi.sam 2>> ${Output_Dir}/Logfile.txt

echo flanking region read mapped to reference genome ..... done

# extract unique reads mapped to reference genome
cat ${Output_Dir}/tmp/left.sam  | awk '{atmark=substr($1,1,1);if(atmark=="@"){print $0}else{AS=substr($14,6);XS=substr($15,6);if(int(AS)>int(XS)){print $0}}}' > ${Output_Dir}/tmp/buf
samtools view -@ ${Threads} -bS ${Output_Dir}/tmp/buf -o ${Output_Dir}/tmp/left_unique.bam

samtools view -@ ${Threads} -h -F 20 ${Output_Dir}/tmp/left_unique.bam -o ${Output_Dir}/tmp/buf
samtools sort -@ ${Threads} ${Output_Dir}/tmp/buf -T hoge -O bam -o ${Output_Dir}/tmp/left_forward_unique.bam
samtools view -@ ${Threads} -h -f 16 ${Output_Dir}/tmp/left_unique.bam -o ${Output_Dir}/tmp/buf
samtools sort -@ ${Threads} ${Output_Dir}/tmp/buf -T hoge -O bam -o ${Output_Dir}/tmp/left_reverse_unique.bam
bedtools bamtobed -i ${Output_Dir}/tmp/left_forward_unique.bam 1> ${Output_Dir}/tmp/left_forward_unique.bed 2>> ${Output_Dir}/Logfile.txt
bedtools bamtobed -i ${Output_Dir}/tmp/left_reverse_unique.bam 1> ${Output_Dir}/tmp/left_reverse_unique.bed 2>> ${Output_Dir}/Logfile.txt

cat ${Output_Dir}/tmp/right.sam  | awk '{atmark=substr($1,1,1);if(atmark=="@"){print $0}else{AS=substr($14,6);XS=substr($15,6);if(int(AS)>int(XS)){print $0}}}' > ${Output_Dir}/tmp/buf
samtools view -@ ${Threads} -bS ${Output_Dir}/tmp/buf -o ${Output_Dir}/tmp/right_unique.bam

samtools view -@ ${Threads} -h -f 16 ${Output_Dir}/tmp/right_unique.bam -o ${Output_Dir}/tmp/buf
samtools sort -@ ${Threads} ${Output_Dir}/tmp/buf -T hoge -O bam -o ${Output_Dir}/tmp/right_forward_unique.bam
samtools view -@ ${Threads} -h -F 20 ${Output_Dir}/tmp/right_unique.bam -o ${Output_Dir}/tmp/buf
samtools sort -@ ${Threads} ${Output_Dir}/tmp/buf -T hoge -O bam -o ${Output_Dir}/tmp/right_reverse_unique.bam
bedtools bamtobed -i ${Output_Dir}/tmp/right_forward_unique.bam 1> ${Output_Dir}/tmp/right_forward_unique.bed 2>> ${Output_Dir}/Logfile.txt
bedtools bamtobed -i ${Output_Dir}/tmp/right_reverse_unique.bam 1> ${Output_Dir}/tmp/right_reverse_unique.bed 2>> ${Output_Dir}/Logfile.txt

# directionでmap dataを分ける
# divide mapping result with directions
samtools view -@ ${Threads} -bS ${Output_Dir}/tmp/left.sam -F 21 -o ${Output_Dir}/tmp/buf
samtools sort -@ ${Threads} ${Output_Dir}/tmp/buf -T hoge -O bam -o ${Output_Dir}/tmp/left_forward.sorted.bam
samtools view -@ ${Threads} -bS ${Output_Dir}/tmp/left.sam -f 16 -o ${Output_Dir}/tmp/buf
samtools sort -@ ${Threads} ${Output_Dir}/tmp/buf -T hoge -O bam -o ${Output_Dir}/tmp/left_reverse.sorted.bam
samtools view -@ ${Threads} -bS ${Output_Dir}/tmp/right.sam -f 16 -o ${Output_Dir}/tmp/buf
samtools sort -@ ${Threads} ${Output_Dir}/tmp/buf -T hoge -O bam -o ${Output_Dir}/tmp/right_forward.sorted.bam
samtools view -@ ${Threads} -bS ${Output_Dir}/tmp/right.sam -F 20 -o ${Output_Dir}/tmp/buf
samtools sort -@ ${Threads} ${Output_Dir}/tmp/buf -T hoge -O bam -o ${Output_Dir}/tmp/right_reverse.sorted.bam

samtools index ${Output_Dir}/tmp/left_forward.sorted.bam
samtools index ${Output_Dir}/tmp/left_reverse.sorted.bam
samtools index ${Output_Dir}/tmp/right_forward.sorted.bam
samtools index ${Output_Dir}/tmp/right_reverse.sorted.bam

# divide multi mapping result with directions
samtools view -@ ${Threads} -bS ${Output_Dir}/tmp/left_multi.sam -F 20 -o ${Output_Dir}/tmp/buf
samtools sort -@ ${Threads} ${Output_Dir}/tmp/buf -T hoge -O bam -o ${Output_Dir}/tmp/left_forward_multi.sorted.bam
samtools view -@ ${Threads} -bS ${Output_Dir}/tmp/left_multi.sam -f 16 -o ${Output_Dir}/tmp/buf
samtools sort -@ ${Threads} ${Output_Dir}/tmp/buf -T hoge -O bam -o ${Output_Dir}/tmp/left_reverse_multi.sorted.bam
samtools view -@ ${Threads} -bS ${Output_Dir}/tmp/right_multi.sam -f 16 -o ${Output_Dir}/tmp/buf
samtools sort -@ ${Threads} ${Output_Dir}/tmp/buf -T hoge -O bam -o ${Output_Dir}/tmp/right_forward_multi.sorted.bam
samtools view -@ ${Threads} -bS ${Output_Dir}/tmp/right_multi.sam -F 20 -o ${Output_Dir}/tmp/buf
samtools sort -@ ${Threads} ${Output_Dir}/tmp/buf -T hoge -O bam -o ${Output_Dir}/tmp/right_reverse_multi.sorted.bam

samtools index ${Output_Dir}/tmp/left_forward_multi.sorted.bam
samtools index ${Output_Dir}/tmp/left_reverse_multi.sorted.bam
samtools index ${Output_Dir}/tmp/right_forward_multi.sorted.bam
samtools index ${Output_Dir}/tmp/right_reverse_multi.sorted.bam

if [ ! -s ${Output_Dir}/tmp/left_forward.sorted.bam -o ! -s ${Output_Dir}/tmp/right_forward.sorted.bam ];then
    echo This family-subgroup-IS does not exist in this sequence genome.;
    exit
fi

# blockに含まれるreadの本数を計算
# insert size - real_length*2 -Read_Length
# calculate mapped reads numbers mapped to reference genome with best hit in each block
bedtools bamtobed -i ${Output_Dir}/tmp/left_forward.sorted.bam  1> ${Output_Dir}/tmp/left_forward.bed 2>> ${Output_Dir}/Logfile.txt
bedtools bamtobed -i ${Output_Dir}/tmp/left_reverse.sorted.bam  1> ${Output_Dir}/tmp/left_reverse.bed 2>> ${Output_Dir}/Logfile.txt
bedtools bamtobed -i ${Output_Dir}/tmp/right_forward.sorted.bam 1> ${Output_Dir}/tmp/right_forward.bed 2>> ${Output_Dir}/Logfile.txt
bedtools bamtobed -i ${Output_Dir}/tmp/right_reverse.sorted.bam 1> ${Output_Dir}/tmp/right_reverse.bed 2>> ${Output_Dir}/Logfile.txt
# compile mapped reads region into block alollowing 225 gap lengths
bedtools merge -d ${block_gap} -i ${Output_Dir}/tmp/left_forward.bed 1> ${Output_Dir}/tmp/left_forward_block.bed
bedtools merge -d ${block_gap} -i ${Output_Dir}/tmp/left_reverse.bed 1> ${Output_Dir}/tmp/left_reverse_block.bed
bedtools merge -d ${block_gap} -i ${Output_Dir}/tmp/right_forward.bed 1> ${Output_Dir}/tmp/right_forward_block.bed
bedtools merge -d ${block_gap} -i ${Output_Dir}/tmp/right_reverse.bed 1> ${Output_Dir}/tmp/right_reverse_block.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN{file1=Output_Dir"/left_forward_block.bed";file2=Output_Dir"/left_forward.bed";OFS="\t";
    while((getline<file1)>0){count=0; while((getline var <file2)>0) { split(var,array); if(int($2)<=int(array[2]) && int(array[3])<=int($3)){count++;}}close(file2);print $1,$2,$3,count}}'\
    > ${Output_Dir}/tmp/left_forward_block_sumRead.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN{file1=Output_Dir"/left_reverse_block.bed";file2=Output_Dir"/left_reverse.bed";OFS="\t";
    while((getline<file1)>0){count=0; while((getline var <file2)>0) { split(var,array); if(int($2)<=int(array[2]) && int(array[3])<=int($3)){count++;}}close(file2);print $1,$2,$3,count}}'\
    > ${Output_Dir}/tmp/left_reverse_block_sumRead.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN{file1=Output_Dir"/right_forward_block.bed";file2=Output_Dir"/right_forward.bed";OFS="\t";
    while((getline<file1)>0){count=0; while((getline var <file2)>0) { split(var,array); if(int($2)<=int(array[2]) && int(array[3])<=int($3)){count++;}}close(file2);print $1,$2,$3,count}}'\
    > ${Output_Dir}/tmp/right_forward_block_sumRead.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN{file1=Output_Dir"/right_reverse_block.bed";file2=Output_Dir"/right_reverse.bed";OFS="\t";
    while((getline<file1)>0){count=0; while((getline var <file2)>0) { split(var,array); if(int($2)<=int(array[2]) && int(array[3])<=int($3)){count++;}}close(file2);print $1,$2,$3,count}}'\
    > ${Output_Dir}/tmp/right_reverse_block_sumRead.bed
# calculate multi mapped reads numbers mapped to reference genome with multi hit in each block
bedtools bamtobed -i ${Output_Dir}/tmp/left_forward_multi.sorted.bam   1> ${Output_Dir}/tmp/left_forward_multi.bed 2>> ${Output_Dir}/Logfile.txt
bedtools bamtobed -i ${Output_Dir}/tmp/left_reverse_multi.sorted.bam   1> ${Output_Dir}/tmp/left_reverse_multi.bed 2>> ${Output_Dir}/Logfile.txt
bedtools bamtobed -i ${Output_Dir}/tmp/right_forward_multi.sorted.bam  1> ${Output_Dir}/tmp/right_forward_multi.bed 2>> ${Output_Dir}/Logfile.txt
bedtools bamtobed -i ${Output_Dir}/tmp/right_reverse_multi.sorted.bam  1> ${Output_Dir}/tmp/right_reverse_multi.bed 2>> ${Output_Dir}/Logfile.txt
# compile multi mapped reads region into block alollowing ${block_gap} gap lengths
bedtools merge -d ${block_gap} -i ${Output_Dir}/tmp/left_forward_multi.bed 1> ${Output_Dir}/tmp/left_forward_multi_block.bed
bedtools merge -d ${block_gap} -i ${Output_Dir}/tmp/left_reverse_multi.bed 1> ${Output_Dir}/tmp/left_reverse_multi_block.bed
bedtools merge -d ${block_gap} -i ${Output_Dir}/tmp/right_forward_multi.bed 1> ${Output_Dir}/tmp/right_forward_multi_block.bed
bedtools merge -d ${block_gap} -i ${Output_Dir}/tmp/right_reverse_multi.bed 1> ${Output_Dir}/tmp/right_reverse_multi_block.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN{file1=Output_Dir"/left_forward_multi_block.bed";file2=Output_Dir"/left_forward_multi.bed";OFS="\t";
    while((getline<file1)>0){count=0; while((getline var <file2)>0) { split(var,array); if(int($2)<=int(array[2]) && int(array[3])<=int($3)){count++;}}close(file2);print $1,$2,$3,count}}'\
    > ${Output_Dir}/tmp/left_forward_multi_block_sumRead.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN{file1=Output_Dir"/left_reverse_multi_block.bed";file2=Output_Dir"/left_reverse_multi.bed";OFS="\t";
    while((getline<file1)>0){count=0; while((getline var <file2)>0) { split(var,array); if(int($2)<=int(array[2]) && int(array[3])<=int($3)){count++;}}close(file2);print $1,$2,$3,count}}'\
    > ${Output_Dir}/tmp/left_reverse_multi_block_sumRead.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN{file1=Output_Dir"/right_forward_multi_block.bed";file2=Output_Dir"/right_forward_multi.bed";OFS="\t";
    while((getline<file1)>0){count=0; while((getline var <file2)>0) { split(var,array); if(int($2)<=int(array[2]) && int(array[3])<=int($3)){count++;}}close(file2);print $1,$2,$3,count}}'\
    > ${Output_Dir}/tmp/right_forward_multi_block_sumRead.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN{file1=Output_Dir"/right_reverse_multi_block.bed";file2=Output_Dir"/right_reverse_multi.bed";OFS="\t";
    while((getline<file1)>0){count=0; while((getline var <file2)>0) { split(var,array); if(int($2)<=int(array[2]) && int(array[3])<=int($3)){count++;}}close(file2);print $1,$2,$3,count}}'\
    > ${Output_Dir}/tmp/right_reverse_multi_block_sumRead.bed

# あとでcutoffをする閾値を決定
# cutoff blocks with threshold
cat ${Output_Dir}/tmp/left_forward_block_sumRead.bed | awk -v BestHit_threshold=${BestHit_threshold} 'BEGIN{OFS="\t"}{if(BestHit_threshold<=int($4)) {print $0}}' \
    > ${Output_Dir}/tmp/left_forward_block_sumRead_cutoff.bed
cat ${Output_Dir}/tmp/left_reverse_block_sumRead.bed | awk -v BestHit_threshold=${BestHit_threshold} 'BEGIN{OFS="\t"}{if(BestHit_threshold<=int($4)) {print $0}}' \
    > ${Output_Dir}/tmp/left_reverse_block_sumRead_cutoff.bed
cat ${Output_Dir}/tmp/right_forward_block_sumRead.bed | awk -v BestHit_threshold=${BestHit_threshold} 'BEGIN{OFS="\t"}{if(BestHit_threshold<=int($4)) {print $0}}' \
    > ${Output_Dir}/tmp/right_forward_block_sumRead_cutoff.bed
cat ${Output_Dir}/tmp/right_reverse_block_sumRead.bed | awk -v BestHit_threshold=${BestHit_threshold} 'BEGIN{OFS="\t"}{if(BestHit_threshold<=int($4)) {print $0}}' \
    > ${Output_Dir}/tmp/right_reverse_block_sumRead_cutoff.bed

# cutoff multi read blocks with threshold
cat ${Output_Dir}/tmp/left_forward_multi_block_sumRead.bed | awk -v MultiHit_threshold=${MultiHit_threshold} 'BEGIN{OFS="\t"}{if(MultiHit_threshold<=int($4)) {print $0}}' \
    > ${Output_Dir}/tmp/left_forward_multi_block_sumRead_cutoff.bed
cat ${Output_Dir}/tmp/left_reverse_multi_block_sumRead.bed | awk -v MultiHit_threshold=${MultiHit_threshold} 'BEGIN{OFS="\t"}{if(MultiHit_threshold<=int($4)) {print $0}}' \
    > ${Output_Dir}/tmp/left_reverse_multi_block_sumRead_cutoff.bed
cat ${Output_Dir}/tmp/right_forward_multi_block_sumRead.bed | awk -v MultiHit_threshold=${MultiHit_threshold} 'BEGIN{OFS="\t"}{if(MultiHit_threshold<=int($4)) {print $0}}' \
    > ${Output_Dir}/tmp/right_forward_multi_block_sumRead_cutoff.bed
cat ${Output_Dir}/tmp/right_reverse_multi_block_sumRead.bed | awk -v MultiHit_threshold=${MultiHit_threshold} 'BEGIN{OFS="\t"}{if(MultiHit_threshold<=int($4)) {print $0}}' \
    > ${Output_Dir}/tmp/right_reverse_multi_block_sumRead_cutoff.bed

bedtools bamtobed -i ${Output_Dir}/tmp/left_forward_unique.bam 1> ${Output_Dir}/tmp/left_forward_unique.bed 2>> ${Output_Dir}/Logfile.txt
bedtools bamtobed -i ${Output_Dir}/tmp/left_reverse_unique.bam 1> ${Output_Dir}/tmp/left_reverse_unique.bed 2>> ${Output_Dir}/Logfile.txt
bedtools bamtobed -i ${Output_Dir}/tmp/right_forward_unique.bam 1> ${Output_Dir}/tmp/right_forward_unique.bed 2>> ${Output_Dir}/Logfile.txt
bedtools bamtobed -i ${Output_Dir}/tmp/right_reverse_unique.bam 1> ${Output_Dir}/tmp/right_reverse_unique.bed 2>> ${Output_Dir}/Logfile.txt

# exclude block not included unique reads (threshold * 0.1)
awk -v Output_Dir=${Output_Dir}/tmp -v UniqueHit_threshold=${UniqueHit_threshold} 'BEGIN{file1=Output_Dir"/left_forward_block_sumRead_cutoff.bed";file2=Output_Dir"/left_forward_unique.bed";OFS="\t";
    while((getline<file1)>0){count=0;while((getline var < file2)>0){split(var,array);if(int($2)<=int(array[2]) && int(array[3]) <= int($3)){count++}}close(file2);if(count >= UniqueHit_threshold){print $0}}}' \
    > ${Output_Dir}/tmp/left_forward_block_sumRead_cutoff_unique.bed
awk -v Output_Dir=${Output_Dir}/tmp -v UniqueHit_threshold=${UniqueHit_threshold} 'BEGIN{file1=Output_Dir"/left_reverse_block_sumRead_cutoff.bed";file2=Output_Dir"/left_reverse_unique.bed";OFS="\t";
    while((getline<file1)>0){count=0;while((getline var < file2)>0){split(var,array);if(int($2)<=int(array[2]) && int(array[3]) <= int($3)){count++}}close(file2);if(count >= UniqueHit_threshold){print $0}}}' \
    > ${Output_Dir}/tmp/left_reverse_block_sumRead_cutoff_unique.bed
awk -v Output_Dir=${Output_Dir}/tmp -v UniqueHit_threshold=${UniqueHit_threshold} 'BEGIN{file1=Output_Dir"/right_forward_block_sumRead_cutoff.bed";file2=Output_Dir"/right_forward_unique.bed";OFS="\t";
    while((getline<file1)>0){count=0;while((getline var < file2)>0){split(var,array);if(int($2)<=int(array[2]) && int(array[3]) <= int($3)){count++}}close(file2);if(count >= UniqueHit_threshold){print $0}}}' \
    > ${Output_Dir}/tmp/right_forward_block_sumRead_cutoff_unique.bed
awk -v Output_Dir=${Output_Dir}/tmp -v UniqueHit_threshold=${UniqueHit_threshold} 'BEGIN{file1=Output_Dir"/right_reverse_block_sumRead_cutoff.bed";file2=Output_Dir"/right_reverse_unique.bed";OFS="\t";
    while((getline<file1)>0){count=0;while((getline var < file2)>0){split(var,array);if(int($2)<=int(array[2]) && int(array[3]) <= int($3)){count++}}close(file2);if(count >= UniqueHit_threshold){print $0}}}' \
    > ${Output_Dir}/tmp/right_reverse_block_sumRead_cutoff_unique.bed

# exclude multi block not included unique reads (threshold * 0.1)
awk -v Output_Dir=${Output_Dir}/tmp -v UniqueHit_threshold=${UniqueHit_threshold} 'BEGIN{file1=Output_Dir"/left_forward_multi_block_sumRead_cutoff.bed";file2=Output_Dir"/left_forward_unique.bed";OFS="\t";
    while((getline<file1)>0){count=0;while((getline var < file2)>0){split(var,array);if(int($2)<=int(array[2]) && int(array[3]) <= int($3)){count++}}close(file2);if(count >= UniqueHit_threshold){print $0}}}' \
    > ${Output_Dir}/tmp/left_forward_multi_block_sumRead_cutoff_unique.bed
awk -v Output_Dir=${Output_Dir}/tmp -v UniqueHit_threshold=${UniqueHit_threshold} 'BEGIN{file1=Output_Dir"/left_reverse_multi_block_sumRead_cutoff.bed";file2=Output_Dir"/left_reverse_unique.bed";OFS="\t";
    while((getline<file1)>0){count=0;while((getline var < file2)>0){split(var,array);if(int($2)<=int(array[2]) && int(array[3]) <= int($3)){count++}}close(file2);if(count >= UniqueHit_threshold){print $0}}}' \
    > ${Output_Dir}/tmp/left_reverse_multi_block_sumRead_cutoff_unique.bed
awk -v Output_Dir=${Output_Dir}/tmp -v UniqueHit_threshold=${UniqueHit_threshold} 'BEGIN{file1=Output_Dir"/right_forward_multi_block_sumRead_cutoff.bed";file2=Output_Dir"/right_forward_unique.bed";OFS="\t";
    while((getline<file1)>0){count=0;while((getline var < file2)>0){split(var,array);if(int($2)<=int(array[2]) && int(array[3]) <= int($3)){count++}}close(file2);if(count >= UniqueHit_threshold){print $0}}}' \
    > ${Output_Dir}/tmp/right_forward_multi_block_sumRead_cutoff_unique.bed
awk -v Output_Dir=${Output_Dir}/tmp -v UniqueHit_threshold=${UniqueHit_threshold} 'BEGIN{file1=Output_Dir"/right_reverse_multi_block_sumRead_cutoff.bed";file2=Output_Dir"/right_reverse_unique.bed";OFS="\t";
    while((getline<file1)>0){count=0;while((getline var < file2)>0){split(var,array);if(int($2)<=int(array[2]) && int(array[3]) <= int($3)){count++}}close(file2);if(count >= UniqueHit_threshold){print $0}}}' \
    > ${Output_Dir}/tmp/right_reverse_multi_block_sumRead_cutoff_unique.bed
echo Block calculation ..... done


# Detect IS in the first phase
bedtools closest -a ${Output_Dir}/tmp/left_forward_block_sumRead_cutoff_unique.bed -b ${Output_Dir}/tmp/right_forward_block_sumRead_cutoff_unique.bed -D a -iu \
    | awk -v IS_length=${IS_length} 'BEGIN{OFS="\t";ex_position="";ex_length=0;num=0}
    {if($9 <= IS_length*1.5 && $5!="." && $2 <= $6 && $3 <= $7){if($9 !="-1" && $6 != ex_position){ex_block_num=$4;ex_position=$6;num++;array[num]=$0}else{if($4 >= ex_block_num){ex_block_num=$4;ex_position=$6;array[num]=$0}}}}
    END{for(i=1;i<=length(array);i++){print array[i]}}' \
    | awk '{if($9 <= 10 || $9 >= 60) {print $0}}' \
    | awk 'BEGIN{OFS="\t"}{print $0,"F",1}' \
    > ${Output_Dir}/tmp/closest_forward.bed
bedtools closest -a ${Output_Dir}/tmp/right_reverse_block_sumRead_cutoff_unique.bed -b ${Output_Dir}/tmp/left_reverse_block_sumRead_cutoff_unique.bed -D a -iu \
    | awk -v IS_length=${IS_length} 'BEGIN{OFS="\t";ex_position="";ex_length=0;num=0}
    {if($9 <= IS_length*1.5 && $5!="." && $2 <= $6 && $3 <= $7){if($9 !="-1" && $6 != ex_position){ex_block_num=$4;ex_position=$6;num++;array[num]=$0}else{if($4 >= ex_block_num){ex_block_num=$4;ex_position=$6;array[num]=$0}}}}
    END{for(i=1;i<=length(array);i++){print array[i]}}' \
    | awk '{if($9 <= 10 || $9 >= 60) {print $0}}' \
    | awk 'BEGIN{OFS="\t"}{print $0,"R",1}' \
    > ${Output_Dir}/tmp/closest_reverse.bed

cat ${Output_Dir}/tmp/closest_forward.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' \
    > ${Output_Dir}/tmp/closest_forward_left.bed
cat ${Output_Dir}/tmp/closest_forward.bed | awk 'BEGIN{OFS="\t"}{print $5,$6,$7,$8}' \
    > ${Output_Dir}/tmp/closest_forward_right.bed
cat ${Output_Dir}/tmp/closest_reverse.bed | awk 'BEGIN{OFS="\t"}{print $5,$6,$7,$8}' \
    > ${Output_Dir}/tmp/closest_reverse_left.bed
cat ${Output_Dir}/tmp/closest_reverse.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' \
    > ${Output_Dir}/tmp/closest_reverse_right.bed

cat ${Output_Dir}/tmp/closest_forward.bed ${Output_Dir}/tmp/closest_reverse.bed \
    | awk 'BEGIN{OFS="\t"}{if($9==0) {print $1,$2,int(($3+$6)/2),$4,$5,int(($3+$6)/2),$7,$8,$9,$10,$11}else{print $0}}' \
    | sort -nk 3 \
    | awk 'BEGIN{flag=0;num=0}{flag=0;if(num!=0)
    {for(i=1;i<=length(result);i++) {split(result[i],array,"\t");sum_blockread=array[4]+array[8];sum_blockread_new=$4+$8;
        if(array[3]<=$3 && $3<=array[6]){if(sum_blockread_new>=sum_blockread){flag=1;result[i]=$0;break} else{flag=2;break}}} if(flag==0){num++;result[num]=$0}}else{num++;result[num]=$0}}
        END{for(i=1;i<=length(result);i++){print result[i]}}' \
    > ${Output_Dir}/closest_FirstStage.bed
echo First Stage ..... done

# best hitで検出できた箇所のblockを除去
# from best-hit mapped read block, extract IS positions detected first IS detection phase
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN {file1=Output_Dir"/left_forward_block_sumRead_cutoff_unique.bed";file2=Output_Dir"/closest_forward_left.bed";
    while((getline<file1)>0) {condition=0;
        while((getline var < file2)>0) {split(var,var_list,"\t");
            if($2<=var_list[2] && var_list[2]<=$3){condition=1;break}else if($2<=var_list[3] && var_list[3]<=$3){condition=1;break}else if(var_list[2]<=$2 && $3<=var_list[3]){condition=1;break}}
            if(condition==0) {print $0} ;close(file2)}}'\
    > ${Output_Dir}/tmp/left_forward_block_sumRead_cutoff_unique_exclude.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN {file1=Output_Dir"/left_reverse_block_sumRead_cutoff_unique.bed";file2=Output_Dir"/closest_reverse_left.bed";
    while((getline<file1)>0) {condition=0;
        while((getline var < file2)>0) {split(var,var_list,"\t");if($2<=var_list[2] && var_list[2]<=$3){condition=1;break} else if($2<=var_list[3] && var_list[3]<=$3){condition=1;break} else if(var_list[2]<=$2 && $3<=var_list[3]){condition=1;break}}
            if(condition==0) {print $0} ; close(file2)}}'\
    > ${Output_Dir}/tmp/left_reverse_block_sumRead_cutoff_unique_exclude.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN {file1=Output_Dir"/right_forward_block_sumRead_cutoff_unique.bed";file2=Output_Dir"/closest_forward_right.bed";
    while((getline<file1)>0) {condition=0;
        while((getline var < file2)>0) {split(var,var_list,"\t");if($2<=var_list[2] && var_list[2]<=$3){condition=1;break} else if($2<=var_list[3] && var_list[3]<=$3){condition=1;break} else if(var_list[2]<=$2 && $3<=var_list[3]){condition=1;break}}
            if(condition==0) {print $0} ; close(file2)}}'\
    > ${Output_Dir}/tmp/right_forward_block_sumRead_cutoff_unique_exclude.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN {file1=Output_Dir"/right_reverse_block_sumRead_cutoff_unique.bed";file2=Output_Dir"/closest_reverse_right.bed";
    while((getline<file1)>0) {condition=0;
        while((getline var < file2)>0) {split(var,var_list,"\t");if($2<=var_list[2] && var_list[2]<=$3){condition=1;break} else if($2<=var_list[3] && var_list[3]<=$3){condition=1;break} else if(var_list[2]<=$2 && $3<=var_list[3]){condition=1;break}}
        if(condition==0) {print $0} ; close(file2)}}'\
    > ${Output_Dir}/tmp/right_reverse_block_sumRead_cutoff_unique_exclude.bed

# unique hitで検出できた箇所のblockをmultihitから除去
# from multi mapped read block, extract IS positions detected first IS detection phase
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN {file1=Output_Dir"/left_forward_multi_block_sumRead_cutoff.bed";file2=Output_Dir"/closest_forward_left.bed";
    while((getline < file1)>0) {condition=0;
        while((getline var < file2)>0) {split(var,var_list,"\t");if($2<=var_list[2] && var_list[3]<=$3){condition=1;break} else if($2<=var_list[3] && var_list[3]<=$3){condition=1;break} else if(var_list[2]<=$2 && $3<=var_list[3]){condition=1;break}}
             if(condition==0) {print $0} ;close (file2)}}'\
    > ${Output_Dir}/tmp/left_forward_multi_block_sumRead_cutoff_exclude.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN {file1=Output_Dir"/left_reverse_multi_block_sumRead_cutoff.bed";file2=Output_Dir"/closest_reverse_left.bed";
    while((getline < file1)>0) {condition=0;
        while((getline var < file2)>0) {split(var,var_list,"\t");if($2<=var_list[2] && var_list[3]<=$3){condition=1;break} else if($2<=var_list[3] && var_list[3]<=$3){condition=1;break} else if(var_list[2]<=$2 && $3<=var_list[3]){condition=1;break}}
            if(condition==0) {print $0} ;close (file2)}}'\
    > ${Output_Dir}/tmp/left_reverse_multi_block_sumRead_cutoff_exclude.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN {file1=Output_Dir"/right_forward_multi_block_sumRead_cutoff.bed";file2=Output_Dir"/closest_forward_right.bed";
    while((getline < file1)>0) {condition=0;
        while((getline var < file2)>0) {split(var,var_list,"\t");if($2<=var_list[2] && var_list[3]<=$3){condition=1;break} else if($2<=var_list[3] && var_list[3]<=$3){condition=1;break} else if(var_list[2]<=$2 && $3<=var_list[3]){condition=1;break}}
            if(condition==0) {print $0} ;close (file2)}}'\
    > ${Output_Dir}/tmp/right_forward_multi_block_sumRead_cutoff_exclude.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN {file1=Output_Dir"/right_reverse_multi_block_sumRead_cutoff.bed";file2=Output_Dir"/closest_reverse_right.bed";
    while((getline < file1)>0) {condition=0;
        while((getline var < file2)>0) {split(var,var_list,"\t");if($2<=var_list[2] && var_list[3]<=$3){condition=1;break} else if($2<=var_list[3] && var_list[3]<=$3){condition=1;break} else if(var_list[2]<=$2 && $3<=var_list[3]){condition=1;break}}
            if(condition==0) {print $0} ;close (file2)}}'\
    > ${Output_Dir}/tmp/right_reverse_multi_block_sumRead_cutoff_exclude.bed

# from multi mapped read unique block, extract IS positons detected first IS detection phase
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN {file1=Output_Dir"/left_forward_multi_block_sumRead_cutoff_unique.bed";file2=Output_Dir"/closest_forward_left.bed";
    while((getline < file1)>0) {condition=0;
        while((getline var < file2)>0) {split(var,var_list,"\t");if($2<=var_list[2] && var_list[3]<=$3){condition=1;break} else if($2<=var_list[3] && var_list[3]<=$3){condition=1;break} else if(var_list[2]<=$2 && $3<=var_list[3]){condition=1;break}}
        if(condition==0) {print $0} ;close (file2)}}'\
    > ${Output_Dir}/tmp/left_forward_multi_block_sumRead_cutoff_unique_exclude.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN {file1=Output_Dir"/left_reverse_multi_block_sumRead_cutoff_unique.bed";file2=Output_Dir"/closest_reverse_left.bed";
    while((getline < file1)>0) {condition=0;
        while((getline var < file2)>0) {split(var,var_list,"\t");if($2<=var_list[2] && var_list[3]<=$3){condition=1;break} else if($2<=var_list[3] && var_list[3]<=$3){condition=1;break} else if(var_list[2]<=$2 && $3<=var_list[3]){condition=1;break}}
        if(condition==0) {print $0} ;close (file2)}}'\
    > ${Output_Dir}/tmp/left_reverse_multi_block_sumRead_cutoff_unique_exclude.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN {file1=Output_Dir"/right_forward_multi_block_sumRead_cutoff_unique.bed";file2=Output_Dir"/closest_forward_right.bed";
    while((getline < file1)>0) {condition=0;
        while((getline var < file2)>0) {split(var,var_list,"\t");if($2<=var_list[2] && var_list[3]<=$3){condition=1;break} else if($2<=var_list[3] && var_list[3]<=$3){condition=1;break} else if(var_list[2]<=$2 && $3<=var_list[3]){condition=1;break}}
        if(condition==0) {print $0} ;close (file2)}}'\
    > ${Output_Dir}/tmp/right_forward_multi_block_sumRead_cutoff_unique_exclude.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN {file1=Output_Dir"/right_reverse_multi_block_sumRead_cutoff_unique.bed";file2=Output_Dir"/closest_reverse_right.bed";
    while((getline < file1)>0) {condition=0;
        while((getline var < file2)>0) {split(var,var_list,"\t");if($2<=var_list[2] && var_list[3]<=$3){condition=1;break} else if($2<=var_list[3] && var_list[3]<=$3){condition=1;break} else if(var_list[2]<=$2 && $3<=var_list[3]){condition=1;break}}
        if(condition==0) {print $0} ;close (file2)}}'\
    > ${Output_Dir}/tmp/right_reverse_multi_block_sumRead_cutoff_unique_exclude.bed

# Detect IS in the second phase
bedtools closest -a ${Output_Dir}/tmp/left_forward_block_sumRead_cutoff_unique_exclude.bed -b ${Output_Dir}/tmp/right_forward_multi_block_sumRead_cutoff_exclude.bed -D a -iu \
    | awk -v IS_length=${IS_length} 'BEGIN{OFS="\t";ex_position="";ex_length=0;num=0}
    {if($9 <= IS_length*1.5 && $5!="." && $2 <= $6 && $3 <= $7){if($5 !="." && $2 <= $6 && $3 <= $7 && $9 !="-1" && $6 != ex_position){ex_block_num=$4;ex_position=$6;num++;array[num]=$0}else{if($4 >= ex_block_num){ex_block_num=$4;ex_position=$6;array[num]=$0}}}}
    END{for(i=1;i<=length(array);i++){print array[i]}}' \
    | awk '{if($9 <= 10 || $9 >= 60) {print $0}}' \
    | awk 'BEGIN{OFS="\t"}{print $0,"F",2}' \
    > ${Output_Dir}/tmp/closest_forward_multi1.bed
bedtools closest -a ${Output_Dir}/tmp/left_forward_multi_block_sumRead_cutoff_exclude.bed -b ${Output_Dir}/tmp/right_forward_block_sumRead_cutoff_unique_exclude.bed -D a -iu \
    | awk -v IS_length=${IS_length} 'BEGIN{OFS="\t";ex_position="";ex_length=0;num=0}
    {if($9 <= IS_length*1.5 && $5!="." && $2 <= $6 && $3 <= $7){if($5 !="." && $2 <= $6 && $3 <= $7 && $9 !="-1" && $6 != ex_position){ex_block_num=$4;ex_position=$6;num++;array[num]=$0}else{if($4 >= ex_block_num){ex_block_num=$4;ex_position=$6;array[num]=$0}}}}
    END{for(i=1;i<=length(array);i++){print array[i]}}' \
    | awk '{if($9 <= 10 || $9 >= 60) {print $0}}' \
    | awk 'BEGIN{OFS="\t"}{print $0,"F",2}' \
    > ${Output_Dir}/tmp/closest_forward_multi2.bed
bedtools closest -a ${Output_Dir}/tmp/right_reverse_multi_block_sumRead_cutoff_exclude.bed -b ${Output_Dir}/tmp/left_reverse_block_sumRead_cutoff_unique_exclude.bed -D a -iu \
    | awk -v IS_length=${IS_length} 'BEGIN{OFS="\t";ex_position="";ex_length=0;num=0}
    {if($9 <= IS_length*1.5 && $5!="." && $2 <= $6 && $3 <= $7){if($5 !="." && $2 <= $6 && $3 <= $7 && $9 !="-1" && $6 != ex_position){ex_block_num=$4;ex_position=$6;num++;array[num]=$0}else{if($4 >= ex_block_num){ex_block_num=$4;ex_position=$6;array[num]=$0}}}}
    END{for(i=1;i<=length(array);i++){print array[i]}}' \
    | awk '{if($9 <= 10 || $9 >= 60) {print $0}}' \
    | awk 'BEGIN{OFS="\t"}{print $0,"R",2}' \
    > ${Output_Dir}/tmp/closest_reverse_multi1.bed
bedtools closest -a ${Output_Dir}/tmp/right_reverse_block_sumRead_cutoff_unique_exclude.bed -b ${Output_Dir}/tmp/left_reverse_multi_block_sumRead_cutoff_exclude.bed -D a -iu \
    | awk -v IS_length=${IS_length} 'BEGIN{OFS="\t";ex_position="";ex_length=0;num=0}
    {if($9 <= IS_length*1.5 && $5!="." && $2 <= $6 && $3 <= $7){if($5 !="." && $2 <= $6 && $3 <= $7 && $9 !="-1" && $6 != ex_position){ex_block_num=$4;ex_position=$6;num++;array[num]=$0}else{if($4 >= ex_block_num){ex_block_num=$4;ex_position=$6;array[num]=$0}}}}
    END{for(i=1;i<=length(array);i++){print array[i]}}' \
    | awk '{if($9 <= 10 || $9 >= 60) {print $0}}' \
    | awk 'BEGIN{OFS="\t"}{print $0,"R",2}' \
    > ${Output_Dir}/tmp/closest_reverse_multi2.bed

cat ${Output_Dir}/tmp/closest_forward_multi1.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' \
    > ${Output_Dir}/tmp/closest_forward_multi1_left.bed
cat ${Output_Dir}/tmp/closest_forward_multi1.bed | awk 'BEGIN{OFS="\t"}{print $5,$6,$7,$8}' \
    > ${Output_Dir}/tmp/closest_forward_multi1_right.bed
cat ${Output_Dir}/tmp/closest_forward_multi2.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' \
    > ${Output_Dir}/tmp/closest_forward_multi2_left.bed
cat ${Output_Dir}/tmp/closest_forward_multi2.bed | awk 'BEGIN{OFS="\t"}{print $5,$6,$7,$8}' \
    > ${Output_Dir}/tmp/closest_forward_multi2_right.bed
cat ${Output_Dir}/tmp/closest_reverse_multi1.bed | awk 'BEGIN{OFS="\t"}{print $5,$6,$7,$8}' \
    > ${Output_Dir}/tmp/closest_reverse_multi1_left.bed
cat ${Output_Dir}/tmp/closest_reverse_multi1.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' \
    > ${Output_Dir}/tmp/closest_reverse_multi1_right.bed
cat ${Output_Dir}/tmp/closest_reverse_multi2.bed | awk 'BEGIN{OFS="\t"}{print $5,$6,$7,$8}' \
    > ${Output_Dir}/tmp/closest_reverse_multi2_left.bed
cat ${Output_Dir}/tmp/closest_reverse_multi2.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' \
    > ${Output_Dir}/tmp/closest_reverse_multi2_right.bed

cat ${Output_Dir}/tmp/closest_forward_multi1.bed ${Output_Dir}/tmp/closest_forward_multi2.bed ${Output_Dir}/tmp/closest_reverse_multi1.bed ${Output_Dir}/tmp/closest_reverse_multi2.bed \
    | awk 'BEGIN{OFS="\t"}{if($9==0) {print $1,$2,int(($3+$6)/2),$4,$5,int(($3+$6)/2),$7,$8,$9,$10,$11}else{print $0}}' \
    | sort -nk 3 \
    | awk 'BEGIN{flag=0;num=0}
    {flag=0;if(num!=0){
        for(i=1;i<=length(result);i++) {split(result[i],array,"\t");sum_blockread=array[4]+array[8];sum_blockread_new=$4+$8;if(array[3]<=$3 && $3<=array[6]){if(sum_blockread_new>=sum_blockread){flag=1;result[i]=$0;break} else{flag=2;break}}}
        if(flag==0){num++;result[num]=$0}}else{num++;result[num]=$0}}END{for(i=1;i<=length(result);i++){print result[i]}}' \
    > ${Output_Dir}/closest_SecondStage.bed
echo Second Stage ..... done

# multi hitで検出できた箇所のblockを除去
# from multi-hit mapped read block, extract IS positions detected second IS detection phase
awk -v Output_Dir=${Output_Dir}/tmp/ 'BEGIN {file1=Output_Dir"/left_forward_multi_block_sumRead_cutoff_exclude.bed";file2=Output_Dir"/closest_forward_multi1_left.bed";
    while((getline<file1)>0) {condition=0;
        while((getline var < file2)>0) {split(var,var_list,"\t");if($2<=var_list[2] && var_list[2]<=$3){condition=1;break} else if($2<=var_list[3] && var_list[3]<=$3){condition=1;break} else if(var_list[2]<=$2 && $3<=var_list[3]){condition=1;break}}
        if(condition==0) {print $0} ; close(file2)}}'\
    > ${Output_Dir}/tmp/left_forward_multi_block_sumRead_cutoff_exclude2.bed
awk -v Output_Dir=${Output_Dir}/tmp/ 'BEGIN {file1=Output_Dir"/left_forward_multi_block_sumRead_cutoff_exclude2.bed";file2=Output_Dir"/closest_forward_multi2_left.bed";
    while((getline<file1)>0) {condition=0;
        while((getline var < file2)>0) {split(var,var_list,"\t");if($2<=var_list[2] && var_list[2]<=$3){condition=1;break} else if($2<=var_list[3] && var_list[3]<=$3){condition=1;break} else if(var_list[2]<=$2 && $3<=var_list[3]){condition=1;break}}
        if(condition==0) {print $0} ; close(file2)}}'\
    > ${Output_Dir}/tmp/left_forward_multi_block_sumRead_cutoff_exclude3.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN {file1=Output_Dir"/left_reverse_multi_block_sumRead_cutoff_exclude.bed";file2=Output_Dir"/closest_reverse_multi1_left.bed";
    while((getline<file1)>0) {condition=0;
        while((getline var < file2)>0) {split(var,var_list,"\t");if($2<=var_list[2] && var_list[2]<=$3){condition=1;break} else if($2<=var_list[3] && var_list[3]<=$3){condition=1;break} else if(var_list[2]<=$2 && $3<=var_list[3]){condition=1;break}}
        if(condition==0) {print $0} ; close(file2)}}'\
    > ${Output_Dir}/tmp/left_reverse_multi_block_sumRead_cutoff_exclude2.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN {file1=Output_Dir"/left_reverse_multi_block_sumRead_cutoff_exclude2.bed";file2=Output_Dir"/closest_reverse_multi2_left.bed";
    while((getline<file1)>0) {condition=0;
        while((getline var < file2)>0) {split(var,var_list,"\t");if($2<=var_list[2] && var_list[2]<=$3){condition=1;break} else if($2<=var_list[3] && var_list[3]<=$3){condition=1;break} else if(var_list[2]<=$2 && $3<=var_list[3]){condition=1;break}}
        if(condition==0) {print $0} ; close(file2)}}'\
    > ${Output_Dir}/tmp/left_reverse_multi_block_sumRead_cutoff_exclude3.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN {file1=Output_Dir"/right_forward_multi_block_sumRead_cutoff_exclude.bed";file2=Output_Dir"/closest_forward_multi1_right.bed";
    while((getline<file1)>0) {condition=0;
        while((getline var < file2)>0) {split(var,var_list,"\t");if($2<=var_list[2] && var_list[2]<=$3){condition=1;break} else if($2<=var_list[3] && var_list[3]<=$3){condition=1;break} else if(var_list[2]<=$2 && $3<=var_list[3]){condition=1;break}}
        if(condition==0) {print $0} ; close(file2)}}'\
    > ${Output_Dir}/tmp/right_forward_multi_block_sumRead_cutoff_exclude2.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN {file1=Output_Dir"/right_forward_multi_block_sumRead_cutoff_exclude2.bed";file2=Output_Dir"/closest_forward_multi2_right.bed";
    while((getline<file1)>0) {condition=0;
        while((getline var < file2)>0) {split(var,var_list,"\t");if($2<=var_list[2] && var_list[2]<=$3){condition=1;break} else if($2<=var_list[3] && var_list[3]<=$3){condition=1;break} else if(var_list[2]<=$2 && $3<=var_list[3]){condition=1;break}}
        if(condition==0) {print $0} ; close(file2)}}'\
    > ${Output_Dir}/tmp/right_forward_multi_block_sumRead_cutoff_exclude3.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN {file1=Output_Dir"/right_reverse_multi_block_sumRead_cutoff_exclude.bed";file2=Output_Dir"/closest_reverse_multi1_right.bed";
    while((getline<file1)>0) {condition=0;
        while((getline var < file2)>0) {split(var,var_list,"\t");if($2<=var_list[2] && var_list[2]<=$3){condition=1;break} else if($2<=var_list[3] && var_list[3]<=$3){condition=1;break} else if(var_list[2]<=$2 && $3<=var_list[3]){condition=1;break}}
        if(condition==0) {print $0} ; close(file2)}}'\
    > ${Output_Dir}/tmp/right_reverse_multi_block_sumRead_cutoff_exclude2.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN {file1=Output_Dir"/right_reverse_multi_block_sumRead_cutoff_exclude2.bed";file2=Output_Dir"/closest_reverse_multi2_right.bed";
    while((getline<file1)>0) {condition=0;
        while((getline var < file2)>0) {split(var,var_list,"\t");if($2<=var_list[2] && var_list[2]<=$3){condition=1;break} else if($2<=var_list[3] && var_list[3]<=$3){condition=1;break} else if(var_list[2]<=$2 && $3<=var_list[3]){condition=1;break}}
        if(condition==0) {print $0} ; close(file2)}}'\
    > ${Output_Dir}/tmp/right_reverse_multi_block_sumRead_cutoff_exclude3.bed

# from multi-hit mapped read unique block, extract IS positions detected second IS detection phase
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN {file1=Output_Dir"/left_forward_multi_block_sumRead_cutoff_unique_exclude.bed";file2=Output_Dir"/closest_forward_multi1_left.bed";
    while((getline<file1)>0) {condition=0;while((getline var < file2)>0) {split(var,var_list,"\t");
    if($2<=var_list[2] && var_list[2]<=$3){condition=1;break} else if($2<=var_list[3] && var_list[3]<=$3){condition=1;break} else if(var_list[2]<=$2 && $3<=var_list[3]){condition=1;break}} if(condition==0) {print $0};
    close(file2)}}'\
    > ${Output_Dir}/tmp/left_forward_multi_block_sumRead_cutoff_unique_exclude2.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN {file1=Output_Dir"/left_forward_multi_block_sumRead_cutoff_unique_exclude2.bed";file2=Output_Dir"/closest_forward_multi2_left.bed";
    while((getline<file1)>0) {condition=0;while((getline var < file2)>0) {split(var,var_list,"\t");
    if($2<=var_list[2] && var_list[2]<=$3){condition=1;break} else if($2<=var_list[3] && var_list[3]<=$3){condition=1;break} else if(var_list[2]<=$2 && $3<=var_list[3]){condition=1;break}} if(condition==0) {print $0};
    close(file2)}}'\
    > ${Output_Dir}/tmp/left_forward_multi_block_sumRead_cutoff_unique_exclude3.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN{file1=Output_Dir"/left_reverse_multi_block_sumRead_cutoff_unique_exclude.bed";file2=Output_Dir"/closest_reverse_multi1_left.bed";
    while((getline<file1)>0) {condition=0;while((getline var < file2)>0) {split(var,var_list,"\t");
    if($2<=var_list[2] && var_list[2]<=$3){condition=1;break} else if($2<=var_list[3] && var_list[3]<=$3){condition=1;break} else if(var_list[2]<=$2 && $3<=var_list[3]){condition=1;break}} if(condition==0) {print $0} ;
    close(file2)}}'\
    > ${Output_Dir}/tmp/left_reverse_multi_block_sumRead_cutoff_unique_exclude2.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN {file1=Output_Dir"/left_reverse_multi_block_sumRead_cutoff_unique_exclude2.bed";file2=Output_Dir"/closest_reverse_multi2_left.bed";
    while((getline<file1)>0) {condition=0;while((getline var < file2)>0) {split(var,var_list,"\t");
    if($2<=var_list[2] && var_list[2]<=$3){condition=1;break} else if($2<=var_list[3] && var_list[3]<=$3){condition=1;break} else if(var_list[2]<=$2 && $3<=var_list[3]){condition=1;break}} if(condition==0) {print $0} ;
    close(file2)}}'\
    > ${Output_Dir}/tmp/left_reverse_multi_block_sumRead_cutoff_unique_exclude3.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN {file1=Output_Dir"/right_forward_multi_block_sumRead_cutoff_unique_exclude.bed";file2=Output_Dir"/closest_forward_multi1_right.bed";
    while((getline<file1)>0) {condition=0;while((getline var < file2)>0) {split(var,var_list,"\t");
    if($2<=var_list[2] && var_list[2]<=$3){condition=1;break} else if($2<=var_list[3] && var_list[3]<=$3){condition=1;break} else if(var_list[2]<=$2 && $3<=var_list[3]){condition=1;break}} if(condition==0) {print $0} ;
    close(file2)}}'\
    > ${Output_Dir}/tmp/right_forward_multi_block_sumRead_cutoff_unique_exclude2.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN {file1=Output_Dir"/right_forward_multi_block_sumRead_cutoff_unique_exclude2.bed";file2=Output_Dir"/closest_forward_multi2_right.bed";
    while((getline<file1)>0) {condition=0;while((getline var < file2)>0) {split(var,var_list,"\t");
    if($2<=var_list[2] && var_list[2]<=$3){condition=1;break} else if($2<=var_list[3] && var_list[3]<=$3){condition=1;break} else if(var_list[2]<=$2 && $3<=var_list[3]){condition=1;break}} if(condition==0) {print $0} ;
    close(file2)}}'\
    > ${Output_Dir}/tmp/right_forward_multi_block_sumRead_cutoff_unique_exclude3.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN {file1=Output_Dir"/right_reverse_multi_block_sumRead_cutoff_unique_exclude.bed";file2=Output_Dir"/closest_reverse_multi1_right.bed";
    while((getline<file1)>0) {condition=0;while((getline var < file2)>0) {split(var,var_list,"\t");
    if($2<=var_list[2] && var_list[2]<=$3){condition=1;break} else if($2<=var_list[3] && var_list[3]<=$3){condition=1;break} else if(var_list[2]<=$2 && $3<=var_list[3]){condition=1;break}} if(condition==0) {print $0} ;
    close(file2)}}'\
    > ${Output_Dir}/tmp/right_reverse_multi_block_sumRead_cutoff_unique_exclude2.bed
awk -v Output_Dir=${Output_Dir}/tmp 'BEGIN {file1=Output_Dir"/right_reverse_multi_block_sumRead_cutoff_unique_exclude2.bed";file2=Output_Dir"/closest_reverse_multi2_right.bed";
    while((getline<file1)>0) {condition=0;while((getline var < file2)>0) {split(var,var_list,"\t");
    if($2<=var_list[2] && var_list[2]<=$3){condition=1;break} else if($2<=var_list[3] && var_list[3]<=$3){condition=1;break} else if(var_list[2]<=$2 && $3<=var_list[3]){condition=1;break}} if(condition==0) {print $0} ;
    close(file2)}}'\
    > ${Output_Dir}/tmp/right_reverse_multi_block_sumRead_cutoff_unique_exclude3.bed

# Detect IS in the third phase
bedtools closest -a ${Output_Dir}/tmp/left_forward_multi_block_sumRead_cutoff_unique_exclude3.bed -b ${Output_Dir}/tmp/right_forward_multi_block_sumRead_cutoff_exclude3.bed -D a -iu \
    | awk -v IS_length=${IS_length} 'BEGIN{OFS="\t";ex_position="";ex_length=0;num=0}
    {if($9 <= IS_length*1.5 && $5!="." && $2 <= $6 && $3 <= $7){if($5 !="." && $2 <= $6 && $3 <= $7 && $9 !="-1" && $6 != ex_position){ex_block_num=$4;ex_position=$6;num++;array[num]=$0}else{if($4 >= ex_block_num){ex_block_num=$4;ex_position=$6;array[num]=$0}}}}
    END{for(i=1;i<=length(array);i++){print array[i]}}' \
    | awk '{if($9 <= 10 || $9 >= 60) {print $0}}' \
    | awk 'BEGIN{OFS="\t"}{print $0,"F",3}' \
    > ${Output_Dir}/tmp/closest_forward_multi_multi1.bed
bedtools closest -a ${Output_Dir}/tmp/left_forward_multi_block_sumRead_cutoff_exclude3.bed -b ${Output_Dir}/tmp/right_forward_multi_block_sumRead_cutoff_unique_exclude3.bed -D a -iu \
    | awk -v IS_length=${IS_length} 'BEGIN{OFS="\t";ex_position="";ex_length=0;num=0}
    {if($9 <= IS_length*1.5 && $5!="." && $2 <= $6 && $3 <= $7){if($5 !="." && $2 <= $6 && $3 <= $7 && $9 !="-1" && $6 != ex_position){ex_block_num=$4;ex_position=$6;num++;array[num]=$0}else{if($4 >= ex_block_num){ex_block_num=$4;ex_position=$6;array[num]=$0}}}}
    END{for(i=1;i<=length(array);i++){print array[i]}}' \
    | awk '{if($9 <= 10 || $9 >= 60) {print $0}}' \
    | awk 'BEGIN{OFS="\t"}{print $0,"F",3}' \
    > ${Output_Dir}/tmp/closest_forward_multi_multi2.bed
bedtools closest -a ${Output_Dir}/tmp/right_reverse_multi_block_sumRead_cutoff_exclude3.bed -b ${Output_Dir}/tmp/left_reverse_multi_block_sumRead_cutoff_unique_exclude3.bed -D a -iu \
    | awk -v IS_length=${IS_length} 'BEGIN{OFS="\t";ex_position="";ex_length=0;num=0}
    {if($9 <= IS_length*1.5 && $5!="." && $2 <= $6 && $3 <= $7){if($5 !="." && $2 <= $6 && $3 <= $7 && $9 !="-1" && $6 != ex_position){ex_block_num=$4;ex_position=$6;num++;array[num]=$0}else{if($4 >= ex_block_num){ex_block_num=$4;ex_position=$6;array[num]=$0}}}}
    END{for(i=1;i<=length(array);i++){print array[i]}}' \
    | awk '{if($9 <= 10 || $9 >= 60) {print $0}}' \
    | awk 'BEGIN{OFS="\t"}{print $0,"R",3}' \
    > ${Output_Dir}/tmp/closest_reverse_multi_multi1.bed
bedtools closest -a ${Output_Dir}/tmp/right_reverse_multi_block_sumRead_cutoff_unique_exclude3.bed -b ${Output_Dir}/tmp/left_reverse_multi_block_sumRead_cutoff_exclude3.bed -D a -iu \
    | awk -v IS_length=${IS_length} '
    BEGIN{OFS="\t";ex_position="";ex_length=0;num=0}
    {if($9 <= IS_length*1.5 && $5!="." && $2 <= $6 && $3 <= $7){if($5 !="." && $2 <= $6 && $3 <= $7 && $9 !="-1" && $6 != ex_position){ex_block_num=$4;ex_position=$6;num++;array[num]=$0}else{if($4 >= ex_block_num){ex_block_num=$4;ex_position=$6;array[num]=$0}}}}
    END{for(i=1;i<=length(array);i++){print array[i]}}' \
    | awk '{if($9 <= 10 || $9 >= 60) {print $0}}' \
    | awk 'BEGIN{OFS="\t"}{print $0,"R",3}' \
    > ${Output_Dir}/tmp/closest_reverse_multi_multi2.bed

cat ${Output_Dir}/tmp/closest_forward_multi_multi1.bed ${Output_Dir}/tmp/closest_forward_multi_multi2.bed ${Output_Dir}/tmp/closest_reverse_multi_multi1.bed ${Output_Dir}/tmp/closest_reverse_multi_multi2.bed \
    | awk 'BEGIN{OFS="\t"}{if($9==0) {print $1,$2,int(($3+$6)/2),$4,$5,int(($3+$6)/2),$7,$8,$9,$10,$11}else{print $0}}' \
    | sort -nk 3 \
    | awk 'BEGIN{flag=0;num=0}
    {flag=0;if(num!=0)
        {for(i=1;i<=length(result);i++) {split(result[i],array,"\t");sum_blockread=array[4]+array[8];sum_blockread_new=$4+$8;
            if(array[3]<=$3 && $3<=array[6]){
                if(sum_blockread_new>=sum_blockread){flag=1;result[i]=$0;break}else{flag=2;break}}}
    if(flag==0){num++;result[num]=$0}}else{num++;result[num]=$0}}
    END{for(i=1;i<=length(result);i++){print result[i]}}' \
    > ${Output_Dir}/closest_ThirdStage.bed
echo Third Stage ..... done

# Final output
cat ${Output_Dir}/closest_FirstStage.bed ${Output_Dir}/closest_SecondStage.bed ${Output_Dir}/closest_ThirdStage.bed \
    | sort -nk 3 \
    | awk 'BEGIN{flag=0;num=0}
    {flag=0;if(num!=0)
        {for(i=1;i<=length(result);i++) {split(result[i],array,"\t");phase=array[11];phase_new=$11;sum_blockread=array[4]+array[8];sum_blockread_new=$4+$8;
            if(array[3]<=$3 && $3<=array[6]){
                if(phase_new < phase) {flag=1;result[i]=$0;break}
                else if( phase_new == phase) {if(sum_blockread_new>=sum_blockread ){flag=1;result[i]=$0;break}}else{flag=2;break}}}
    if(flag==0){num++;result[num]=$0}}else{num++;result[num]=$0}}
    END{for(i=1;i<=length(result);i++){print result[i]}}' \
    | uniq \
    > ${Output_Dir}/closest_final.bed

if [ ${No_Clean} -eq 1 ]; then
    rm -f ${Output_Dir}/tmp/;
fi
