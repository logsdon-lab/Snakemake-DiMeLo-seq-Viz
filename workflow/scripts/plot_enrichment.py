import os
import json
import toml
import argparse
from typing import TextIO

import polars as pl
import cenplot as cplt


DEF_TOML = """
[settings]
format = ["png", "pdf"]
transparent = false
dim = [ 10, 8]
legend_prop = 0.05
axis_h_pad = 0.01
dpi = 600
legend_pos = "left"
layout = "constrained"

[[tracks]]
position = "relative"
type = "bar"
proportion = 0.05
path = "dimelo"

[tracks.options]
hide_x = true
ymin = 0
legend_title = "Normalized CENP-A\\nDiMeLo-Seq signals\\n(relative to IgG)"
legend_title_only = true
"""

def generate_window_pileup(
    pileup: TextIO,
    window: int
) -> pl.DataFrame:
    df = pl.read_csv(
        pileup,
        has_header=False,
        separator="\t",
        columns=[0, 1, 2, 3],
        new_columns=["chrom", "st", "end", "value"]
    )
    return (
        df.with_columns(
            min_st=pl.col("st").min().over("chrom"),
            max_end=pl.col("st").max().over("chrom")
        ).with_columns(
            pl.col("st") - pl.col("min_st"),
            pl.col("end") - pl.col("min_st"),
        )
        .with_columns(idx=pl.col("st") // window)
        .group_by("chrom", "idx")
        .agg(
            value=pl.col("value").mean(),
            min_st=pl.col("min_st").first(),
            max_end=pl.col("max_end").first()
        )
    )

def main():
    ap = argparse.ArgumentParser()
    # Both must be sorted by first column.
    ap.add_argument("-t", "--treatment_pileup", type=argparse.FileType("rb"), required=True)
    ap.add_argument("-c", "--control_pileup", type=argparse.FileType("rb"), required=True)
    ap.add_argument("-w", "--window", type=int, default=5000)
    ap.add_argument("-j", "--json_tracks", type=str, default=None)
    ap.add_argument("-b", "--base_cenplot_cfg", type=argparse.FileType("rb"), default=None)
    ap.add_argument("-o", "--output_dir", type=str, default="./output")

    args = ap.parse_args()
    output_dir = os.path.dirname(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    df_treatment = generate_window_pileup(args.treatment_pileup, args.window)
    df_control = generate_window_pileup(args.control_pileup, args.window)
    
    df_both = (
        df_treatment
        .join(df_control, on=["chrom", "idx"], how="full")
        .with_columns(
            pl.col("idx").fill_null(pl.col("idx_right")),
            pl.col("chrom").fill_null(pl.col("chrom_right")),
            pl.col("min_st").fill_null(pl.col("min_st_right")),
            pl.col("value").fill_null(pl.lit(0.0))
        )
        .with_columns(
            value_diff=pl.col("value") - pl.col("value_right"),
            st=pl.col("min_st") + (pl.col("idx") * args.window),
            end=pl.col("min_st") + ((pl.col("idx") + 1) * args.window),
        )
        .select("chrom", "st", "end", "value_diff")
        .sort("chrom", "st")
    )

    if args.json_tracks:
        dfs_tracks = json.loads(args.json_tracks)
        assert isinstance(tracks, dict), "Invalid tracks format. Expects {dtype: path}"
        dfs_tracks = {
            dtype: pl.read_csv(path, separator="\t", has_header=False).rename({"column_0": "chrom"})
            for dtype, path in tracks.items()
        }
    else:
        dfs_tracks = {}

    for chrom, df_cenpa in df_both.partition_by("chrom", as_dict=True).items():
        chrom = chrom[0]
        if args.base_cenplot_cfg:
            cfg = toml.load(args.base_cenplot_cfg)
        else:
            cfg = toml.loads(DEF_TOML)
        df_dtypes = {
            **{dtype: df.filter(pl.col("chrom") == chrom) for dtype, df in dfs_tracks.items()},
            "dimelo": df_cenpa
        }
        tracks_chrom = []
        for track in cfg["tracks"]:
            dtype = track.get("path")
            if dtype and isinstance(df_dtypes.get(dtype), pl.DataFrame):
                outfile = os.path.join(output_dir, f"{chrom}_{dtype}.csv")
                df_dtypes[dtype].write_csv(outfile, separator="\t", include_header=False)
                new_track = track.copy()
                new_track["path"] = outfile
                tracks_chrom.append(new_track)
            else:
                tracks_chrom.append(track)

        cfg["tracks"] = tracks_chrom
        output_cfg = os.path.join(output_dir, f"{chrom}.toml")
        with open(output_cfg, "wt") as fh:
            toml.dump(cfg, fh)
        
        with open(output_cfg, "rb") as fh:
            tracks, settings = cplt.read_tracks(fh)

        _ = cplt.plot_tracks(tracks.tracks, settings=settings, outdir=output_dir, chrom=chrom)

if __name__ == "__main__":
    raise SystemExit(main())
