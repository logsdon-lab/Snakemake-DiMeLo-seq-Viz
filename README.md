# Snakemake_DiMeLo-seq_Viz
Workflow to visualize DiMeLo-seq data.
> WIP

```yaml
align:
  # Number of zscores (stdev) to filter abnormal m6A reads.
  filter_lt_m6a_zscore:
  - 2.0
  # Minimum read length
  min_read_length: 0
  samples:
    - name: HG00513
      # Assembly to align to
      asm_fa: /project/logsdon_shared/projects/HGSVC3/new_65_asms_renamed/HG00513-asm-renamed-reort.fa
      treatment:
        name: CENPA
        reads:
        # Reads
        - /project/logsdon_shared/long_read_archive/unsorted/20250729_DMLseq_HG00513_ULK114/D5C/20250729_1253_3B_PAY70585_ab403a2d/pod5/basecalling/20250729_DMLseq_HG00513_ULK114.bam
      control:
        name: "IgG"
        reads:
        - /project/logsdon_shared/long_read_archive/unsorted/20250729_DMLseq_HG00513_ULK114/D5I/20250729_1253_3D_PBC11952_e0f51ee7/pod5/basecalling/20250729_DMLseq_HG00513_ULK114.bam
        - /project/logsdon_shared/projects/HGSVC3/NA20355_Dimelo/tmp_HG00513_IgG_bam/20250730_DMLseq_HG00513_ULK114.bam

  aligner: "minimap2"
  aligner_opts: "-y -a --eqx --cs -x lr:hqae -I8g -s 4000"
  output_dir: "results_HG00513/align"
  logs_dir: "logs_HG00513/align"
  benchmarks_dir: "benchmarks_HG00513/align"
  threads_aln: 24
  mem_aln: 50G
  samtools_view_flag: 2308
```

```bash
snakemake -p --configfile config.yaml --sdm conda --executor lsf --rerun-triggers mtime -j 20 -n
```
