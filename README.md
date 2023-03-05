# SOMvc-pipeline_hg19_to_hg38
normal-tumor pair somatic variant calling pipeline - special version for hg19 mapped files to hg38


before running, you have to set up the attached Docker images:
```sh
docker build -t somvc-pipeline_hg19 https://raw.githubusercontent.com/loipf/SOMvc-pipeline_hg19_to_hg38/master/docker/Dockerfile
```

---
### run mapping pipeline

```sh
nextflow run loipf/SOMvc-pipeline_hg19_to_hg38 --project_dir /path/to/folder --sample_match_file normal_tumor_pairs.csv --reference_genome genome.fa.gz -resume --bed_file bed_file_filtered.bed --num_threads 10 -with-docker somvc-pipeline_hg19_to_hg38




