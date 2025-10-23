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
