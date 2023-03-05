FROM conda/miniconda2@sha256:d30ea858b9d3425127c2e835d39653cdca74368af832d65748cf912f6c1da67b

MAINTAINER loipf

RUN apt-get update && apt install -y wget unzip gzip procps

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge


### bioinformatic tools
### strelka needs python 2.7 - conflict with multiqc (replace with older version)
RUN conda install \
		gatk4==4.2.2.0 \
		lofreq==2.1.5 \
		manta==1.6.0 \
		strelka==2.9.10 \
		samtools==1.12.0 \
		bedtools==2.30.0 \
		vt==2015.11.10
RUN conda install -y numpy scipy
RUN apt install -y r-base 

RUN conda install bwa==0.7.17 
RUN conda install cutadapt==1.18


### download bcftools and set up - with conda package it leads to conflicts
RUN wget https://github.com/samtools/bcftools/releases/download/1.14/bcftools-1.14.tar.bz2
RUN tar -xjf bcftools-1.14.tar.bz2 -C /usr/src/
RUN rm bcftools-1.14.tar.bz2
WORKDIR "/usr/src/bcftools-1.14/"
RUN ./configure --prefix /usr/src/bcftools/
RUN make
RUN make install
ENV PATH="/usr/src/bcftools/bin:${PATH}"
WORKDIR "/"

### download VarDictJava - faster than perl conda version
RUN wget https://github.com/AstraZeneca-NGS/VarDictJava/releases/download/v1.8.2/VarDict-1.8.2.tar
RUN tar -xf VarDict-1.8.2.tar -C /usr/src/
RUN rm VarDict-1.8.2.tar


### download SomaticCombiner
RUN wget https://github.com/mingyi-wang/somatic-combiner/raw/03ebc539d9a037662fde42f05854b179d2249479/example/somaticCombiner.jar
RUN mv somaticCombiner.jar /usr/src/somaticCombiner.jar


### install Conpair (compatible version with gatk4)
RUN wget https://github.com/vladsaveliev/Conpair/archive/0f4d8fe8eb6d86f093ae4c75f0feb406311fb6a9.zip -O conpair.zip
RUN unzip conpair.zip
RUN mv Conpair-* conpair 
RUN mv conpair /usr/src/
RUN pip install /usr/src/conpair


### download mutect2 helper files - even though they are pretty outdated its better than nothing
RUN mkdir /usr/src/mutect2_genome
WORKDIR "/usr/src/mutect2_genome"
RUN wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz -O 1000g_pon.hg38.vcf.gz
#RUN wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi -O /usr/src/mutect2_genome/1000g_pon.hg38.vcf.gz.tbi
RUN wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz -O af-only-gnomad.hg38.vcf.gz
#RUN wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi -O /usr/src/mutect2_genome/af-only-gnomad.hg38.vcf.gz.tbi

RUN wget -N ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt
RUN grep -e '^[^#]' GCF_000001405.39_GRCh38.p13_assembly_report.txt | awk -F "\t" '{ print $1,",",$10 }' | sed -r 's/\s+//g' | awk -F "," '{ print $2, $1}' > UCSC_to_ensembl.chrnames
RUN bcftools annotate --rename-chrs UCSC_to_ensembl.chrnames --threads 10 -Oz -o 1000g_pon_ensembl.hg38.vcf.gz 1000g_pon.hg38.vcf.gz
RUN bcftools annotate --rename-chrs UCSC_to_ensembl.chrnames --threads 10 -Oz -o af-only-gnomad_ensembl.hg38.vcf.gz af-only-gnomad.hg38.vcf.gz
RUN tabix 1000g_pon_ensembl.hg38.vcf.gz
RUN tabix af-only-gnomad_ensembl.hg38.vcf.gz
RUN rm GCF_000001405.39_GRCh38.p13_assembly_report.txt 1000g_pon.hg38.vcf.gz af-only-gnomad.hg38.vcf.gz UCSC_to_ensembl.chrnames

WORKDIR "/"







