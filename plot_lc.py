"""
Plot light curves from a ztf_bulk_lc.py output directory.

  python plot_lc.py ztf_alerce                 # first 9 objects, one panel each
  python plot_lc.py ztf_dr --n 12 -o dr.png
  python plot_lc.py ztf_alerce --id ZTF22aabxfxz

Reads every data/*.parquet in the directory and works out the schema itself
(DR nested lightcurves, ALeRCE detection/forced_photometry, or Fink).
"""
import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

BANDS = {1: ("g", "tab:green"), 2: ("r", "tab:red"), 3: ("i", "tab:orange")}


def load(outdir):
    """-> long table: object, mjd, band, mag, magerr, kind."""
    files = sorted(glob.glob(os.path.join(outdir, "data", "*.parquet")))
    if not files:
        sys.exit(f"no parquet files under {outdir}/data/")
    frames = []
    for f in files:
        df = pd.read_parquet(f)
        if "lightcurve" in df.columns:                      # --source dr
            lc = df["lightcurve"]
            if not hasattr(lc, "nest"):   # plain pandas reads the struct as dicts
                lc = lc.map(lambda v: {k: list(x) for k, x in v.items()})
                flat = pd.DataFrame({
                    "object": df["target_id"], "band": df["filterid"],
                    "mjd": lc.map(lambda v: v["hmjd"]),
                    "mag": lc.map(lambda v: v["mag"]),
                    "magerr": lc.map(lambda v: v["magerr"]),
                    "catflags": lc.map(lambda v: v["catflags"]),
                }).explode(["mjd", "mag", "magerr", "catflags"])
            else:
                flat = lc.nest.to_flat().reset_index()
                flat = flat.join(df[["target_id", "filterid"]], on="index").rename(
                    columns={"target_id": "object", "filterid": "band", "hmjd": "mjd"})
            flat = flat[flat["catflags"].astype("int64") == 0]   # IRSA's clean-epoch cut
            frames.append(pd.DataFrame({
                "object": flat["object"], "mjd": flat["mjd"].astype(float),
                "band": flat["band"], "mag": flat["mag"].astype(float),
                "magerr": flat["magerr"].astype(float), "kind": "dr"}))
        elif "i:objectId" in df.columns:                    # --source fink
            tag = df.get("d:tag", pd.Series("valid", index=df.index))
            frames.append(pd.DataFrame({
                "object": df["i:objectId"], "mjd": df["mjd"], "band": df["i:fid"],
                "mag": df["i:magpsf"].where(tag != "upperlim", df["i:diffmaglim"]),
                "magerr": df.get("i:sigmapsf"),
                "kind": tag.map(lambda t: "upperlim" if t == "upperlim" else "alert")}))
        else:                                               # --source alerce
            is_forced = df.get("source_table", "") == "ztf.forced_photometry"
            mag = df["mag"] if "mag" in df.columns else df["magpsf"]
            err = df["e_mag"] if "e_mag" in df.columns else df["sigmapsf"]
            if "magpsf" in df.columns and "mag" in df.columns:
                mag = df["mag"].where(is_forced, df["magpsf"])
                err = df["e_mag"].where(is_forced, df["sigmapsf"])
            frames.append(pd.DataFrame({
                "object": df["oid"], "mjd": df["mjd"], "band": df["fid"],
                "mag": mag, "magerr": err,
                "kind": pd.Series("forced", index=df.index).where(is_forced, "alert")}))
    lc = pd.concat(frames, ignore_index=True)
    lc = lc[lc["mag"].notna() & (lc["mag"] > 0) & (lc["mag"] < 30)]
    # display cut: drop the very low-S/N forced-photometry points whose error
    # bars would otherwise set the y-scale (they are kept in the parquet files)
    err = pd.to_numeric(lc["magerr"], errors="coerce")
    lc = lc[err.isna() | ((err > 0) & (err < 1))]
    return lc.sort_values(["object", "mjd"])


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("outdir")
    p.add_argument("--n", type=int, default=9, help="how many objects to plot")
    p.add_argument("--id", default=None, help="plot just this object")
    p.add_argument("-o", "--out", default=None, help="default: <outdir>_lc.png")
    a = p.parse_args()

    lc = load(a.outdir)
    objs = [a.id] if a.id else list(lc.groupby("object").size().sort_values(ascending=False).index[:a.n])
    ncol = 1 if len(objs) == 1 else 3
    nrow = -(-len(objs) // ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 3.1 * nrow), squeeze=False)

    for ax, obj in zip(axes.ravel(), objs):
        d = lc[lc["object"] == obj]
        for fid, (name, col) in BANDS.items():
            b = d[d["band"] == fid]
            if not len(b):
                continue
            for kind, style in (("upperlim", dict(marker="v", ls="none", ms=4, alpha=.35)),
                                ("forced", dict(marker=".", ls="none", ms=5, alpha=.75)),
                                ("alert", dict(marker="o", ls="none", ms=4.5)),
                                ("dr", dict(marker="o", ls="none", ms=3.5))):
                k = b[b["kind"] == kind]
                if not len(k):
                    continue
                lbl = f"{name} {kind}" if kind in ("forced", "upperlim") else name
                if kind == "upperlim":
                    ax.plot(k["mjd"], k["mag"], color=col, label=lbl, **style)
                else:
                    ax.errorbar(k["mjd"], k["mag"], yerr=k["magerr"], color=col,
                                ecolor=col, elinewidth=.7, capsize=0, label=lbl, **style)
        lo, hi = d["mag"].quantile([0.01, 0.99])
        pad = max(0.3, 0.15 * (hi - lo))
        ax.set_ylim(hi + pad, lo - pad)   # inverted
        ax.set_title(f"{obj}  ({len(d)} pts)", fontsize=9)
        ax.set_xlabel("MJD", fontsize=8)
        ax.set_ylabel("mag", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=6, ncol=2, framealpha=.6)
    for ax in axes.ravel()[len(objs):]:
        ax.axis("off")

    out = a.out or a.outdir.rstrip("/") + "_lc.png"
    fig.suptitle(f"{a.outdir}: {lc['object'].nunique()} objects, {len(lc):,} points", fontsize=10)
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    print(f"{out}  ({lc['object'].nunique()} objects, {len(lc):,} points)")


if __name__ == "__main__":
    main()
