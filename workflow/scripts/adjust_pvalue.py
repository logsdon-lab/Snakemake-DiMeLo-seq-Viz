import sys
import polars as pl
from scipy.stats import false_discovery_control


def main():
    infile = sys.argv[1]
    cohen_h = float(sys.argv[2])
    cohen_h_lower_bound = float(sys.argv[3])
    p_value = float(sys.argv[4])

    df_dmr = pl.read_csv(infile, has_header=True, separator="\t")
    p_values_adj = false_discovery_control(df_dmr["map_pvalue"])
    (
        df_dmr.with_columns(map_pvalue_adj=pl.Series(p_values_adj))
        .filter(
            pl.col("map_pvalue_adj").lt(p_value)
            & (
                pl.col("cohen_h").gt(cohen_h)
                & pl.col("cohen_h_low").gt(cohen_h_lower_bound)
            )
            | (
                pl.col("cohen_h").lt(-cohen_h)
                & pl.col("cohen_h_low").lt(-cohen_h_lower_bound)
            )
        )
        .sort("#chrom", "start")
        .select("#chrom", "start", "end")
        .write_csv(file=sys.stdout, include_header=False, separator="\t")
    )


if __name__ == "__main__":
    raise SystemExit(main())
