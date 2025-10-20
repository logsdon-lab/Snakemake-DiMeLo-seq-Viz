# Snakemake_DiMeLo-seq_Viz
Workflow to visualize DiMeLo-seq data. For HGSVC centromere paper.

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

To run on UPenn's LPC.
```bash
snakemake -p --configfile config.yaml --sdm conda --executor lsf --rerun-triggers mtime -j 20 -n
```

## `mbamstats`
A binary for a rust script to filter abnormal m6a reads is provided in `workflow/scripts/mbamstats`. This will filter reads by:
1. Collect reads m6a modification stats for each read:
    * Count all possible AT sites.
    * Get all possible m6A sites and filter based on if meet modification threshold of 80%.
    * Calculate the valid proportion of m6A sites by dividing the total number of valid sites over all possible sites.
2. Calculate the mean and stdev valid proportion across all reads.
3. Then, for each read, calculate the valid proporition z-score to be filtered downstream.

To recompile, `rust` is required.
```bash
pushd workflow/scripts/mbamstats
cargo build --release
rm ../bin_mbamstats && cp target/release/mbamstats ../bin_mbamstats
popd
```

## `analysis_workflow_unfilter30.py`
Script contributed by Shenghan Gao that intersects each pileup from the control and treatment with the centromere region, normalizes the treatment signal by the control, and generates `cenplot` images.

> [!NOTE]
> Only functions on UPenn's LPC. Must be rewritten to fit use-case.

Briefly, it:
* Generates 5 kbp non-overlapping windows of the control and treatment pileup.
* Intersects the pileups with a given centromere array bed file.
* Subtracts the control pileup from the treatment pileup.
* Generate a cenplot config toml file.
* Plots the CENP-A signal.
