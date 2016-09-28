{\rtf1\ansi\ansicpg1252\cocoartf1404\cocoasubrtf470
{\fonttbl\f0\fmodern\fcharset0 Courier;\f1\fmodern\fcharset0 Courier-Bold;\f2\fnil\fcharset0 TrebuchetMS;
\f3\fnil\fcharset0 Menlo-Regular;}
{\colortbl;\red255\green255\blue255;\red243\green235\blue0;\red27\green29\blue31;}
\margl1440\margr1440\vieww18180\viewh12860\viewkind0
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0

\f0\fs26 \cf0 \CocoaLigature0 # ahcg_pipeline\
## 
\f1\b Variant calling pipeline for genomic data analysis
\f0\b0 \
\
## Requirements\
\
1. [Python3 - version 3.4.1](https://www.python.org/download/releases/3.4.1/)\
2. [Trimomatic - version 0.36](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip)\
3. [Bowtie2 - version 2.2.9](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/)\
4. [Picard tools - version 2.6.0](https://github.com/broadinstitute/picard/releases/download/2.6.0/picard.jar)\
5. [GATK - version 3.4](https://software.broadinstitute.org/gatk/download/)\
\
## Reference genome\
\
Reference genomes can be downloaded from [Illumina iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.html)\
\
## Test data\
\
Use the following protocol to download and prepare test dataset from NIST sample NA12878\
\
	```\{sh\}\
	wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/	NIST7035_TAAGGCGA_L001_R1_001.fastq.gz\
	wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/	NIST7035_TAAGGCGA_L001_R2_001.fastq.gz\
	gunzip NIST7035_TAAGGCGA_L001_R1_001.fastq.gz\
	gunzip NIST7035_TAAGGCGA_L001_R2_001.fastq.gz\
	head -100000 NIST7035_TAAGGCGA_L001_R1_001.fastq > test_r1.fastq\
	head -100000 NIST7035_TAAGGCGA_L001_R2_001.fastq > test_r2.fastq	```\
\
## Help\
\
To access help use the following command:\
\
```\{sh\}\
python3 ahcg_pipeline.py -h\
```\
\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\qc\partightenfactor0
\cf0 Git\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0
\cf0 \
###### Changing the git URL (from Shashi\'92s to mine)\
My github folder Somineni21 was forked from shashidhar22/ahcg_pipeline. Then the URL (obtained by clicking on the \'93clone or download\'94) was changed from \'93https://github.com/shashidhar22/ahcg_pipeline.git\'94 to \'93https://github.com/Somineni21/ahcg_pipeline.git\'94\
\
## Set up pipeline dependencies\
- Download reference genome and dbSNP \
\
\pard\pardeftab720\partightenfactor0
\cf0 \expnd0\expndtw0\kerning0
\CocoaLigature1 	```\{sh\}\
	wget www.prism.gatech.edu/~sravishankar9/resources.tar.gz\
	tar -zxvf ./resources.tar.gz\
	gunzip ./resources/dbsnp/dbsnp_138.hg19.vcf.gz\
	```\kerning1\expnd0\expndtw0 \CocoaLigature0 \
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\qc\partightenfactor0
\cf0 \
Extra files needed\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0
\cf0 \

\f1\b ** 1
\f0\b0 . 
\f1\b Genome index file - Indexing the reference genome using _samtools faidx_ **\
\
	
\f0\b0 \'91\'92\'92\{sh\}\
	\expnd0\expndtw0\kerning0
\CocoaLigature1 samtools faidx <path to ref.fa>\
	\'91\'92\'92\

\f1\b \kerning1\expnd0\expndtw0 \CocoaLigature0 \
** 2. Genome dictionary - Building genome.dict file using _Picard tools_ **\
\
	
\f0\b0 \'91\'92\'92\{sh\}\
	\expnd0\expndtw0\kerning0
\CocoaLigature1 java -jar picard.jar CreateSequenceDictionary R=ref.fa O=ref.dict\
	\'91\'92\'92\
\kerning1\expnd0\expndtw0 \CocoaLigature0 Note: _-java don\'92t need path as long as you have it somewhere on your local machine. However, picard.jar needs to be given the path_\

\f1\b \
** 3. Bowtie index file - using _bowtie2/bowtie2-build_ **\
\
	
\f0\b0 \'91\'92\'92\{sh\}\
	\expnd0\expndtw0\kerning0
\CocoaLigature1 bowtie2-build <path to ref.fa> <output prefix>\
	\'91\'92\'92\kerning1\expnd0\expndtw0 \CocoaLigature0 \
\
* All the above 3 files were made by using the ref genome file, hg19.fa as input file *\
\
** 1. Genome index file **\
This file is needed for the GATK (it informs GATK, about the position of your read in the reference genome).\
In order to make the Genome index file we need samtools:\
\
** Downloading & Installing samtools **\
\
Samtools were downloaded from the respective website and were installed using the following commands:\
\
	\'91\'92\'92\{sh\}\
	./configure\
	make\
	make install\
	\'91\'92\'92\
Then using the samtools, the ref genome file: hg19.fa was indexed and saved as hg19.fi.fasta by using the command: \
\
	\'91\'92\'92\{sh\}\
	/path/to/samtools-1.3.1/samtools faidx hg19.fa\
	\'91\'92\'92\
\pard\pardeftab720\partightenfactor0
\cf0 ** 2. Genome dictionary **\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0
\cf0 This file is needed for the GATK (this file tells GATK, what the length of the ref genome, start and end positions of all the chromosomes etc).\
\
In order to make the Genome dictionary, we need Picard tools (provided by Shashi)\
\
- Creating genome dictionary file using Picard tools\
\
Genome dictionary of the ref genome, hg19.fa was created using the Picard tools as follows:\
\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardeftab720\ri0\partightenfactor0
\cf0 \cb2 \CocoaLigature1 ahcg_pipeline/resources/genome$ java -jar ../../lib/picard.jar CreateSequenceDictionary R=hg19.fa O=hg19_reference.dict\cb1 \
\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0
\cf0 \CocoaLigature0 3. Bowtie index file\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0
\cf0 \cb2 This file is needed for the samtools. It informs samtools that the following read is in line X of the ref genome index - not sure, check with Shashi?\cb1 \
\
Prebuilt Bowtie index files are available, however, I couldn\'92t use it cause it was not compatible with the hg19 ref genome version.\
To circumvent the issue, I got the appropriate bowtie index files from Wilson (flash drive - copied on to my desktop). Then I moved it from my desktop to the remote server (vagrant) using the following command:\
\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0
\cf3 \CocoaLigature1 	\'91\'92\'92\{sh\}\
	rsync -av -e 'ssh -p 2222' ./Desktop/bowtie_index vagrant@localhost:~/ahcg_pipeline/\cf0 \
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardeftab720\ri0\partightenfactor0
\cf0 	\'91\'92\'92\
## \CocoaLigature0 Running the .py script to call variants\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0
\cf0 Since I now have all the required files in the right place, I then ran Shashi\'92s .py script as follows:\
\
	\'91\'92\'92\{sh\}\
	/home/vagrant/ahcg_pipeline/lib/Python3/bin/python3.4 ahcg_pipeline.py -t /home/vagrant/ahcg_pipeline/lib/	Trimmomatic-0.36/trimmomatic-0.36.jar -b /home/vagrant/ahcg_pipeline/lib/bowtie2-2.2.9/bowtie2 -p /home/vagrant/	ahcg_pipeline/lib/picard.jar -g /home/vagrant/ahcg_pipeline/lib/GenomeAnalysisTK.jar -i /home/vagrant/	ahcg_pipeline/resources/testdata/test_r1.fastq /home/vagrant/ahcg_pipeline/resources/testdata/test_r2.fastq -w /	home/vagrant/ahcg_pipeline/resources/bowtie_index/hg19 -d /home/vagrant/ahcg_pipeline/resources/dbsnp/	dbsnp_138.hg19.vcf -r /home/vagrant/ahcg_pipeline/resources/genome/hg19.fa -a /home/vagrant/ahcg_pipeline/lib/	Trimmomatic-0.36/adapters/NexteraPE-PE.fa -o /home/vagrant/ahcg_pipeline/out/\
	\'91\'92\'92\
\pard\pardeftab720\partightenfactor0
\cf0 \
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0
\cf0 ## To ignore the unwanted files from pushing into the github\
go to vi .gitignore,\
\
then add the name of the file/folder/directory to the list and save it. So next time when you push your stuff to the course (master) github repository, all the stuff under the .gitignore will be excluded\
\
\
\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\qc\partightenfactor0

\f1\b \cf0 09/08/2016
\f0\b0 \
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0
\cf0 ## __NA12878__\
\
Downloaded all the 4 Bam files (from the same subject, NA12878, sequenced at 4 different times) from \'93https://github.com/genome-in-a-bottle/giab_data_indexes/blob/master/NA12878 alignment.index.NA12878_HiSeq_Exome_Garvan_GRCh37_09252015\'94\
\
	\'91\'92\'92\{sh\}\
	wget file.name\
	\'91\'92\'92\
Then all the 4 files were merged to increase the number of reads which thus improves the accuracy of variant calling by using the command\
\
	\'91\'92\'92\{sh\}\
	samtools merge <output.bam> <input.bam1> \'85\'85 <input.bamN>\
	\'91\'92\'92\
Extracting reads mapping to a region of interest (BRCA1, using genomic coordinates as BED file)\
\
	\'91\'92\'92\{sh\}\
	samtools view -L <bed file> -b -o <output bam file> <input bam file>\
	\'91\'92\'92\
Note: _-b just specifies that the output needs to be a bam file_\
\
\pard\pardeftab720\partightenfactor0

\f2 \cf0 \expnd0\expndtw0\kerning0
\CocoaLigature1 Then 
\f0 \kerning1\expnd0\expndtw0 \CocoaLigature0 the output.bam file from above (BRCA1 region extracted from the merged bam file) is then converted to fastq using bedtools\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0
\cf0 \
	\'91\'92\'92\{sh\}\
	bedtools bamtofastq -i <input file.bam> -fq <output file.fastq>\
	\'91\'92\'92\
\
## __Create ssh key for github__\
\
	\'91\'92\'92\{sh\}\
	ssh -add ~/.ssh/id_rsa\
	cd ~/.ssh/\
	vi id_rsa.pub\
	\'91\'92\'92\
copy the ssh key into your github account\
follow instructions at:\
https://help.github.com/articles/adding-a-new-ssh-key-to-your-github-account/\
\
\
09/13/2016\
\
	```\{sh\}\
	/home/vagrant/ahcg_pipeline/lib/Python3/bin/python3.4 ahcg_pipeline.py -t /home/vagrant/ahcg_pipeline/lib/	Trimmomatic-0.36/trimmomatic-0.36.jar -b /home/vagrant/ahcg_pipeline/lib/bowtie2-2.2.9/bowtie2 -p /home/vagrant/	ahcg_pipeline/lib/picard.jar -g /home/vagrant/ahcg_pipeline/lib/GenomeAnalysisTK.jar -i /home/vagrant/	ahcg_pipeline/resources/synthetic_data/reads_B1_27000x150bp_0S_0I_0D_0U_0N_1.fq /home/vagrant/ahcg_pipeline/	resources/synthetic_data/reads_B1_27000x150bp_0S_0I_0D_0U_0N_2.fq -w /home/vagrant/ahcg_pipeline/resources/	genome/chr17_ref_bowtie_indexed -d /home/vagrant/ahcg_pipeline/resources/dbsnp/dbsnp_138.hg19.vcf -r /home/	vagrant/ahcg_pipeline/resources/genome/chr17.fa -a /home/vagrant/ahcg_pipeline/lib/Trimmomatic-0.36/adapters/	NexteraPE-PE.fa -o /home/vagrant/ahcg_pipeline/out/\
	```\
	```\{sh\}\
	/ahcg_pipeline/lib/samtools-1.3.1/samtools faidx /ahcg_pipeline/resources/genome/chr17.fa\
	```\
\
	```\{sh\}\
	java -jar ./lib/picard.jar CreateSequenceDictionary R=./resources/genome/chr17.fa O=chr17.dict\
	```\
\
	```\{sh\}\
	./lib/bowtie2-2.2.9/bowtie2 ./resources/genome/chr17.fa ./resources/genome/chr17_ref_bowtie_indexed
\f3\fs22 \

\f0\fs26 	```\
\

\f1\b 09/14/2016\

\f0\b0 \
## File conversions\
** Convert sam file to bam file **\
\
	\'91\'92\'92\{sh\}\
	samtools view -Sb file.sam > file.bam\
\
qsub -q sh.q script.sh\
09/15/16\
\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0

\fs24 \cf0 wget http://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_nextera/nexterarapidcapture/nexterarapidcapture_exome_targetedregions_v1.2.bed\
\
\

\f1\b 09/28/2016 (Assignment 4)
\f0\b0 \
\
__Verifying the variants called (using our pipeline) with the GIAB_highconf_vcf calls__\
\
Variants that fall into the regions of interest (gene panel) were pulled out from the NA12878_GIAB_highconf.vcf were extracted using the following command\
\
	```\{sh\}\
	/home/vagrant/ahcg_pipeline/lib/bedtools2/bin/bedtools intersect -header -a NA12878_GIAB_highconf.vcf -b ./cancer_geneList_coordinates_noChr.txt > GIAB_variantsCalled\
	\'91\'92\'92\
\
From our pipeline:\
\
	```\{sh\}\
	/home/vagrant/ahcg_pipeline/lib/bedtools2/bin/bedtools intersect -header -a NA12878_exome_variants.vcf -b ./	cancer_geneList_coordinates.txt > NA12878_exome_variants_called\
	\'91\'92\'92\
To get rid off \'93chr\'94 from each line of the VCF file:\
\
	```\{sh\}\
	awk '\{gsub(/^chr/,""); print\}' NA12878_exome_variants_called > no_chr.NA12878_exome.vcf\
	```\
\
To pull out the variants that are common to both the VCF files (ours vs GIAM_highconf): \
\
	```\{sh\}\
	/home/vagrant/ahcg_pipeline/lib/bedtools2/bin/bedtools intersect -header -a GIAB_variantsCalled -b 	no_chr.NA12878_exome.vcf > intersect.vcf\
	```}