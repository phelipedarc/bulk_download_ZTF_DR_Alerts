"""
Bulk download of ZTF light curves for a target list (ZTF ID, RA, Dec).

Sources (--source):

  dr      ZTF Data Release PSF photometry (default DR24) from IRSA's HATS/Parquet
          light-curve catalog on public AWS S3, cross-matched to your targets with
          lsdb. This is IRSA's recommended route for >~1000 light curves; it reads
          Parquet partitions directly instead of one HTTP query per object.
          Output keeps the DR's nested format: one row per (target, DR objectid,
          filter) with the epochs in a nested `lightcurve` column
          (hmjd, mag, magerr, clrcoeff, catflags, ...).

  alerce  Difference-imaging photometry from the ALeRCE TAP service (ADQL,
          many IDs per query): ztf.detection and, optionally,
          ztf.forced_photometry (--alerce-tables). ALeRCE documents TAP as the
          interface for bulk selection, and it is the fastest route here:
          ~25-30 targets/s, so ~57k targets in well under an hour.
          Note this is forced photometry at ZTF *alert* positions, so it only
          covers objects that produced alerts. For arbitrary positions use
          IRSA's batch forced photometry service (registration required,
          1500 positions per submission, 15000 in the system at once).

  fink    Difference-imaging alert photometry from the Fink REST API,
          detections + upper limits. Fink documents this endpoint for simple
          queries, not bulk: it returns HTTP 500 once a request would return
          too much data (~20 objectIds, sometimes fewer), and roughly 15% of
          chunks fail whatever the chunk size. Fink's own bulk route is their
          separate Data Transfer service (login + fink-client). Prefer alerce.

Targets are processed in chunks (--chunk-size); every chunk is written as one
Parquet file and recorded in chunks.csv, so an interrupted run resumes where it
stopped. At the end everything is packed into <outdir>.zip.

Plan / estimate before the full run:
  python ztf_bulk_lc.py targets.csv --source dr --plan        # bytes DR crossmatch must read
  python ztf_bulk_lc.py targets.csv --source dr --sample 500  # real run on a subset + extrapolation

Input: CSV/TSV with a header; columns matched case-insensitively:
  id: ztf_id, ztfid, objectid, oid, id, name | ra: ra, ra_deg | dec: dec, dec_deg, decl
  RA/Dec in decimal degrees (ICRS). Broker sources need alert-style IDs (ZTF18abcdefg).

Dependencies:
  dr     : lsdb>=0.8.1 s3fs pyarrow       (tested with lsdb 0.11.0)
  alerce : pyvo pyarrow pandas
  fink   : requests pyarrow pandas
"""
import argparse
import csv
import json
import logging
import os
import sys
import time
import zipfile
from datetime import datetime, timezone

import pandas as pd

DR_CATALOG = "s3://ipac-irsa-ztf/ztf/enhanced/dr24/lc/hats"
DR_COLUMNS = ["objectid", "objra", "objdec", "filterid", "nepochs", "lightcurve"]
DR_FRAC = 0.935  # measured share of a partition's bytes held by DR_COLUMNS
ALERCE_TAP = "https://tap.alerce.online/tap"
FINK_URL = "https://api.ztf.fink-portal.org"

log = logging.getLogger("ztf_bulk_lc")


# ------------------------------------------------------------------ helpers
def date_to_mjd(s):
    dt = datetime.fromisoformat(s)
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return dt.timestamp() / 86400.0 + 40587.0  # MJD of the Unix epoch


def read_targets(path):
    with open(path, newline="") as f:
        head = f.read(8192)
        f.seek(0)
        try:
            dialect = csv.Sniffer().sniff(head, delimiters=",;\t")
        except csv.Error:
            dialect = csv.excel
        reader = csv.DictReader(f, dialect=dialect)
        cols = {c.strip().lower(): c for c in reader.fieldnames}

        def pick(*names):
            for n in names:
                if n in cols:
                    return cols[n]
            sys.exit(f"Input needs one of the columns {names}; found {reader.fieldnames}")

        idc = pick("ztf_id", "ztfid", "objectid", "oid", "id", "name")
        rac, decc = pick("ra", "ra_deg"), pick("dec", "dec_deg", "decl")
        rows = [(r[idc].strip(), float(r[rac]), float(r[decc]))
                for r in reader if (r[idc] or "").strip()]
    df = pd.DataFrame(rows, columns=["target_id", "ra", "dec"])
    df = df.drop_duplicates("target_id").reset_index(drop=True)
    return df


def sky_sort(df):
    """Order targets along a HEALPix curve so each chunk covers a compact sky area
    (each DR partition is then read by as few chunks as possible)."""
    from hats.pixel_math.spatial_index import compute_spatial_index
    idx = compute_spatial_index(df["ra"].to_numpy(), df["dec"].to_numpy())
    return df.assign(_sidx=idx).sort_values("_sidx").drop(columns="_sidx").reset_index(drop=True)


# values the services put in numeric columns to mean "missing"
JUNK = {"", "nan", "NaN", "null", "NULL", "None", "NA", "\x00", "\ufffd"}


def _is_junk(v):
    return isinstance(v, str) and (v.strip() in JUNK or any(ord(c) < 32 or ord(c) == 0xFFFD for c in v))


def clean_for_parquet(df):
    """Give every column one Arrow-writable type.

    The services return columns with mixed python types: ALeRCE sends nullable
    BIGINTs as int/str, Fink fills non-detection rows with placeholders
    ('nan', '\\x00', '\\ufffd') in columns that are numeric for detections.

    A column is numeric when every value that is not one of those placeholders
    parses as a number -- a rule that does not depend on how many placeholders
    a given chunk happens to hold, so every chunk of a run gets the same
    schema. Integers never pass through float, so 19-digit candids keep every
    digit. Nested values become JSON; anything else stays text.
    """
    for c in df.columns:
        s = df[c]
        if s.dtype != object:
            continue
        vals = s[s.notna()]
        if vals.map(lambda v: isinstance(v, (dict, list))).any():
            df[c] = s.map(lambda v: json.dumps(v) if isinstance(v, (dict, list)) else v).astype("string")
            continue
        real = vals[~vals.map(_is_junk)]
        if not real.size:
            df[c] = s.astype("string")
            continue
        if {type(v) for v in real} <= {bool}:
            df[c] = s.where(~s.map(_is_junk), None).astype("boolean")
            continue
        num = pd.to_numeric(real, errors="coerce")
        if not num.notna().all():                 # some real value is not a number
            df[c] = s.astype("string")
            continue
        num = pd.to_numeric(s, errors="coerce")   # placeholders -> NaN
        ok = num.notna()
        if (num[ok] % 1 == 0).all():              # integral -> exact Int64
            if num[ok].abs().max() <= 2 ** 53:
                df[c] = num.astype("Int64")
            else:                                 # too big for float64: parse the text
                try:
                    df[c] = s.where(~s.map(_is_junk), None).astype("string").astype("Int64")
                except (TypeError, ValueError):
                    df[c] = s.astype("string")
        else:
            df[c] = num
    return df


def load_chunk_log(path):
    if not os.path.exists(path):
        return {}
    with open(path, newline="") as f:
        return {r["chunk"]: r for r in csv.DictReader(f) if r["status"] == "ok"}


def append_chunk_log(path, row):
    new = not os.path.exists(path)
    with open(path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["chunk", "status", "ntargets", "nmatched",
                                          "nrows", "bytes", "seconds", "error"])
        if new:
            w.writeheader()
        w.writerow(row)


# ------------------------------------------------------------------ DR (HATS)
def open_dr(args):
    import lsdb
    so = {"anon": True} if args.dr_catalog.startswith("s3://") else None
    return lsdb.open_catalog(args.dr_catalog, columns=DR_COLUMNS, storage_options=so)


def dr_plan(args, targets, zcat):
    """Which DR partitions the crossmatch touches and how many bytes they hold."""
    import re
    import numpy as np
    from hats.pixel_math.spatial_index import compute_spatial_index
    idx = compute_spatial_index(targets["ra"].to_numpy(), targets["dec"].to_numpy())
    pixels = zcat.get_healpix_pixels()
    by_order = {}
    for o in {p.order for p in pixels}:
        by_order[o] = set(np.unique(idx >> (2 * (29 - o))).tolist())
    touched = [p for p in pixels if p.pixel in by_order[p.order]]

    # one recursive listing of the dataset, then add up the files per partition
    # (partitions are Norder=<o>/Dir=<d>/Npix=<n>/part*.parquet in recent HATS)
    base = zcat.hc_structure.catalog_base_dir
    pat = re.compile(r"Norder=(\d+)/Dir=\d+/Npix=(\d+)")
    sizes = {}
    for path, info in base.fs.find(f"{base.path.rstrip('/')}/dataset", detail=True).items():
        m = pat.search(path)
        if m and path.endswith(".parquet"):
            k = (int(m.group(1)), int(m.group(2)))
            sizes[k] = sizes.get(k, 0) + info["size"]
    tot = sum(sizes.get((p.order, p.pixel), 0) for p in touched)
    whole = sum(sizes.values())

    print(f"\nDR catalog       : {args.dr_catalog}")
    print(f"targets          : {len(targets)}")
    print(f"partitions       : {len(touched)} of {len(pixels)} touched")
    print(f"Parquet to read  : {tot / 1e12:.2f} TB of {whole / 1e12:.2f} TB in the catalog")
    print(f"                   ({DR_FRAC:.0%} of that is actually fetched: the {len(DR_COLUMNS)} "
          f"requested columns, of which the nested `lightcurve` alone is ~86%)")
    eff = tot * DR_FRAC
    print(f"                   => {eff / 1e12:.2f} TB over the wire, + the 10\" margin cache")
    for label, mbps in (("100 Mbit/s", 100), ("500 Mbit/s", 500), ("1 Gbit/s", 1000)):
        print(f"  transfer @ {label:<12}: {eff * 8 / (mbps * 1e6) / 3600:7.1f} h")
    print("Output is far smaller than this (only matched rows are kept); run --sample to measure it.")


def dr_chunk(args, chunk_df, zcat):
    import lsdb
    tc = lsdb.from_dataframe(chunk_df, ra_column="ra", dec_column="dec",
                             catalog_name="targets", margin_threshold=None)
    x = tc.crossmatch(zcat, radius_arcsec=args.radius, n_neighbors=args.n_neighbors,
                      suffixes=("", "_ztf"), suffix_method="overlapping_columns")
    cond = f"lightcurve.hmjd >= {args.mjd_min:.6f}"
    if args.mjd_max is not None:
        cond += f" and lightcurve.hmjd <= {args.mjd_max:.6f}"
    res = x.query(cond).compute()
    res = res.reset_index(drop=True)
    nrows = int(res["lightcurve"].nest.to_flat().shape[0]) if len(res) else 0
    return res, res["target_id"].nunique() if len(res) else 0, nrows


# ------------------------------------------------------------------ ALeRCE (TAP)
def alerce_chunk(args, chunk_df, tap):
    ids = chunk_df["target_id"].tolist()
    id_list = ",".join("'" + i.replace("'", "''") + "'" for i in ids)
    frames = []
    for table in args.alerce_tables:
        tcol = args._alerce_tcol.get(table)
        where = f"oid IN ({id_list})"
        if tcol:
            where += f" AND {tcol} >= {args.mjd_min:.6f}"
            if args.mjd_max is not None:
                where += f" AND {tcol} <= {args.mjd_max:.6f}"
        adql = f"SELECT * FROM {table} WHERE {where}"
        job = tap.run_async(adql, maxrec=args.maxrec) if args.tap_async else tap.search(adql, maxrec=args.maxrec)
        df = job.to_table().to_pandas()
        if len(df) >= args.maxrec:
            raise RuntimeError(f"{table}: hit maxrec={args.maxrec}; lower --chunk-size")
        df.insert(0, "source_table", table)
        frames.append(df)
    res = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    nmatch = res["oid"].nunique() if len(res) else 0
    return res, nmatch, len(res)


def alerce_setup(args):
    import pyvo
    tap = pyvo.dal.TAPService(args.alerce_tap)
    args._alerce_tcol = {}
    for table in args.alerce_tables:  # discover columns once
        cols = tap.search(f"SELECT TOP 1 * FROM {table}").to_table().colnames
        tcol = next((c for c in ("mjd", "MJD") if c in cols), None)
        if tcol is None:
            log.warning("%s: no 'mjd' column (%s...); no time cut applied", table, cols[:12])
        if "oid" not in cols:
            sys.exit(f"{table} has no 'oid' column: {cols}")
        args._alerce_tcol[table] = tcol
        log.info("%s columns: %s", table, ", ".join(cols))
    return tap


# ------------------------------------------------------------------ Fink (API)
def fink_chunk(args, chunk_df, sess):
    payload = {"objectId": ",".join(chunk_df["target_id"]), "withupperlim": True,
               "output-format": "json"}
    if args.fink_columns:
        payload["columns"] = args.fink_columns
    r = sess.post(f"{args.fink_url}/api/v1/objects", json=payload, timeout=args.timeout)
    r.raise_for_status()
    df = pd.DataFrame(r.json())
    if len(df):
        df.insert(0, "mjd", df["i:jd"] - 2400000.5)
        keep = df["mjd"] >= args.mjd_min
        if args.mjd_max is not None:
            keep &= df["mjd"] <= args.mjd_max
        df = df[keep].reset_index(drop=True)
    nmatch = df["i:objectId"].nunique() if len(df) else 0
    return df, nmatch, len(df)


def fink_setup(args):
    import requests
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry
    s = requests.Session()
    s.mount("https://", HTTPAdapter(max_retries=Retry(
        total=args.retries, backoff_factor=1, status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset(["POST"]), respect_retry_after_header=True)))
    return s


# ------------------------------------------------------------------ driver
def run(args):
    targets = read_targets(args.targets)
    n_full = len(targets)
    log.info("%d unique targets", n_full)
    if args.sample:
        targets = targets.sample(n=min(args.sample, n_full), random_state=args.seed)

    if args.source == "dr":
        zcat = open_dr(args)
        if args.plan:
            return dr_plan(args, targets, zcat)
        targets, ctx, fn = sky_sort(targets), zcat, dr_chunk
    elif args.source == "alerce":
        ctx, fn = alerce_setup(args), alerce_chunk
    else:
        ctx, fn = fink_setup(args), fink_chunk

    data_dir = os.path.join(args.outdir, "data")
    os.makedirs(data_dir, exist_ok=True)
    targets.to_csv(os.path.join(args.outdir, "targets_used.csv"), index=False)
    log_path = os.path.join(args.outdir, "chunks.csv")
    done = load_chunk_log(log_path)

    cs = args.chunk_size
    chunks = [targets.iloc[i:i + cs] for i in range(0, len(targets), cs)]
    todo = [(k, c) for k, c in enumerate(chunks) if f"{k:05d}" not in done]
    log.info("%d chunks of <=%d targets; %d done, %d to run", len(chunks), cs,
             len(chunks) - len(todo), len(todo))

    t0 = time.time()
    for n, (k, cdf) in enumerate(todo, 1):
        name = f"{k:05d}"
        path = os.path.join(data_dir, f"{args.source}_chunk{name}.parquet")
        tc = time.time()
        try:
            res, nmatch, nrows = fn(args, cdf.reset_index(drop=True), ctx)
            nbytes = 0
            if len(res):
                clean_for_parquet(res).to_parquet(path + ".part", compression="zstd")
                os.replace(path + ".part", path)
                nbytes = os.path.getsize(path)
            row = dict(chunk=name, status="ok", ntargets=len(cdf), nmatched=nmatch,
                       nrows=nrows, bytes=nbytes, seconds=f"{time.time() - tc:.1f}", error="")
        except Exception as e:  # noqa: BLE001 - log and continue; rerun retries it
            row = dict(chunk=name, status="error", ntargets=len(cdf), nmatched=0, nrows=0,
                       bytes=0, seconds=f"{time.time() - tc:.1f}", error=f"{type(e).__name__}: {e}"[:300])
            log.error("chunk %s failed: %s", name, row["error"])
        append_chunk_log(log_path, row)
        el = time.time() - t0
        log.info("chunk %s: %s targets matched, %s rows, %.1f MB, %ss | %d/%d | ETA %.1f h",
                 name, row["nmatched"], row["nrows"], int(row["bytes"]) / 1e6, row["seconds"],
                 n, len(todo), el / n * (len(todo) - n) / 3600)
        if args.pause:
            time.sleep(args.pause)
    wall = time.time() - t0

    zip_path = args.outdir.rstrip("/") + ".zip"
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_STORED, allowZip64=True) as zf:  # Parquet is already zstd-compressed
        for f in ("chunks.csv", "targets_used.csv"):
            zf.write(os.path.join(args.outdir, f), f)
        for f in sorted(os.listdir(data_dir)):
            if f.endswith(".parquet"):
                zf.write(os.path.join(data_dir, f), f"data/{f}")
    log.info("Wrote %s (%.2f GB)", zip_path, os.path.getsize(zip_path) / 1e9)

    ok = list(load_chunk_log(log_path).values())
    nerr = len(chunks) - len(ok)
    if nerr:
        log.warning("%d chunks failed; rerun the same command to retry them", nerr)
    if args.sample and ok and todo:
        nt = sum(int(r["ntargets"]) for r in ok)
        f = n_full / nt
        b = sum(int(r["bytes"]) for r in ok)
        print(f"\n=== Extrapolated from {nt} sampled targets to {n_full} ===")
        print(f"  targets with data : {sum(int(r['nmatched']) for r in ok) / nt:.1%}")
        print(f"  rows (epochs)     : {sum(int(r['nrows']) for r in ok) * f:,.0f}")
        print(f"  Parquet / zip     : {b * f / 1e9:.2f} GB (zip ~ same; peak disk ~2x)")
        tt = wall * f
        if args.source == "dr":
            print(f"  wall time         : <= {tt / 3600:.1f} h (random sample is spread over the "
                  "sky, so it overestimates: the full sky-sorted run reads each partition ~once; "
                  "see --plan)")
        else:
            print(f"  wall time         : {tt / 3600:.1f} h")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("targets")
    p.add_argument("--source", choices=["dr", "alerce", "fink"], default="dr")
    p.add_argument("-o", "--outdir", default=None, help="default: ztf_<source>")
    p.add_argument("--start", default="2019-01-01", help="UTC date, default 2019-01-01")
    p.add_argument("--end", default=None, help="UTC date, default: no limit")
    p.add_argument("--chunk-size", type=int, default=None,
                   help="targets per chunk (default dr 5000, alerce 1000, fink 10; "
                        "the Fink API 500s on requests that return too much data)")
    p.add_argument("--sample", type=int, default=0, help="run on N random targets and extrapolate")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--pause", type=float, default=0, help="seconds between chunks (be polite to brokers)")
    g = p.add_argument_group("dr")
    g.add_argument("--plan", action="store_true", help="only report partitions/bytes the crossmatch reads")
    g.add_argument("--dr-catalog", default=DR_CATALOG,
                   help="HATS collection (default DR24 on S3; use a local copy or a later DR if available)")
    g.add_argument("--radius", type=float, default=1.5, help="crossmatch radius, arcsec")
    g.add_argument("--n-neighbors", type=int, default=20,
                   help="max DR objectids kept per target (one per filter x field/CCD-quadrant)")
    g = p.add_argument_group("alerce")
    g.add_argument("--alerce-tap", default=ALERCE_TAP)
    g.add_argument("--alerce-tables", default="ztf.detection",
                   help="comma list, e.g. ztf.detection,ztf.forced_photometry")
    g.add_argument("--maxrec", type=int, default=2_000_000)
    g.add_argument("--async", dest="tap_async", action="store_true",
                   help="submit TAP jobs asynchronously (sync is the default and is much faster here)")
    g = p.add_argument_group("fink")
    g.add_argument("--fink-url", default=FINK_URL)
    g.add_argument("--fink-columns", default=None,
                   help="comma list to shrink output (default: all columns)")
    g.add_argument("--timeout", type=float, default=600)
    g.add_argument("--retries", type=int, default=3)
    args = p.parse_args()

    args.outdir = args.outdir or f"ztf_{args.source}"
    args.chunk_size = args.chunk_size or {"dr": 5000, "alerce": 1000, "fink": 10}[args.source]
    args.alerce_tables = [t.strip() for t in args.alerce_tables.split(",") if t.strip()]
    args.mjd_min = date_to_mjd(args.start)
    args.mjd_max = date_to_mjd(args.end) if args.end else None
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S")
    run(args)


if __name__ == "__main__":
    main()
