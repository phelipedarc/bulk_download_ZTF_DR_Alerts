# ZTF bulk light curves

Download ZTF light curves for a list of targets (ID, RA, Dec), in chunks, resumable,
written as Parquet.

Three sources, because they are three different measurements:

| `--source` | what you get | photometry |
|---|---|---|
| `alerce` | ALeRCE TAP: `ztf.detection` and `ztf.forced_photometry` | difference image (host subtracted) |
| `dr` | ZTF DR24 PSF photometry from IRSA's HATS/Parquet catalog on AWS S3 | science image (host included) |
| `fink` | Fink REST API: detections and upper limits | difference image |

For fitting transients, you want difference-image photometry, so `alerce` is usually the
right choice. `dr` is PSF photometry on the science images, which includes the host galaxy
flux; it suits variable stars and AGN, not a supernova on a bright host.

## Install

```bash
pip install -r requirements.txt
```

Only the dependencies for the source you use are needed. `alerce` needs `pyvo`, `fink`
needs `requests`, `dr` needs `lsdb` and `s3fs`.

## Input

CSV or TSV with a header. Column names are matched case-insensitively:

- id: `ztf_id`, `ztfid`, `objectid`, `oid`, `id`, `name`
- RA: `ra`, `ra_deg`
- Dec: `dec`, `dec_deg`, `decl`

RA and Dec in decimal degrees (ICRS). See `targets_example.csv`.

`alerce` and `fink` match on the ID, so they need alert-style IDs (`ZTF18abcdefg`).
`dr` matches on coordinates and works with any ID.

## Usage

```bash
# difference-image photometry, detections plus forced photometry
python ztf_bulk_lc.py targets.csv --source alerce --alerce-tables ztf.detection,ztf.forced_photometry

# DR24 PSF photometry
python ztf_bulk_lc.py targets.csv --source dr --chunk-size 200

# how many bytes the DR crossmatch would read, without running it
python ztf_bulk_lc.py targets.csv --source dr --plan

# try 500 random targets first and extrapolate to the full list
python ztf_bulk_lc.py targets.csv --source alerce --sample 500
```

Useful options:

- `--start` / `--end` UTC dates to cut the time range (default from 2019-01-01)
- `--chunk-size` targets per chunk (default: dr 5000, alerce 1000, fink 10)
- `--pause` seconds between chunks
- `-o` output directory (default `ztf_<source>`)

## Output

```
ztf_alerce/
  targets_used.csv
  chunks.csv                       one row per chunk: status, counts, bytes, seconds
  data/alerce_chunk00000.parquet
  data/alerce_chunk00001.parquet
ztf_alerce.zip
```

Each chunk is written as a separate Parquet file and logged in `chunks.csv`. Rerunning
the same command skips chunks already done and retries the failed ones, so an interrupted
run resumes where it stopped.

Read it back with:

```python
import pandas as pd, glob
df = pd.concat([pd.read_parquet(f) for f in glob.glob("ztf_alerce/data/*.parquet")])
```

For `--source dr` the epochs are in a nested `lightcurve` column, stored as a Parquet
struct of lists (`hmjd`, `mag`, `magerr`, `clrcoeff`, `catflags`). Plain pandas reads each
row as a dict of arrays; `nested_pandas.read_parquet` restores it as a nested column.
Cut on `catflags == 0` for clean epochs.

## Plotting

```bash
python plot_lc.py ztf_alerce --n 9
python plot_lc.py ztf_dr --id ZTF22aabxfxz
```

Writes `<outdir>_lc.png`. Works out the schema from the files, so it takes the output of
any of the three sources. It drops points with a magnitude error above 1 mag so the very
low signal-to-noise forced photometry does not set the y-scale; those points are still in
the Parquet files.

## Scale

Measured on a 2000-target list with median 20 detections per object:

| source | rate | 57,000 targets | output |
|---|---|---|---|
| `alerce` | 25-30 targets/s | 30-40 min | ~350 MB |
| `fink` | ~16 targets/s | ~1 h | ~1 GB |
| `dr` | limited by bytes read | 1-2.5 days | a few GB |

`dr` is a different kind of job. The DR24 catalog is 12,485 Parquet partitions totalling
9.47 TB. A crossmatch reads whole partitions, so 57,000 targets spread over the
extragalactic sky touch about 7,600 of them, roughly 4 TB over the wire to extract a few
GB. Column selection saves little: the nested `lightcurve` column alone is 86% of each
file, and the six columns requested are 93.5% of the bytes. Run `--plan` on your own list
for the exact figure.

Two notes for `dr`:

- Memory. 30 targets in a single chunk reached 48 GB resident. Use a small
  `--chunk-size` (a few hundred), not the 5000 default.
- Do not size the job from `--sample`. Targets are processed in sky order so each
  partition is read about once, but a random sample is scattered and reads a partition per
  target, which overestimates badly. Use `--plan` instead.

## Service limits

- **ALeRCE** documents TAP as the interface for bulk selection and publishes no rate
  limit. Queries of 1000 IDs return in a few seconds and are faster than splitting them.
  Direct database access was retired in favour of TAP.
- **Fink's** REST API is documented for simple queries, not bulk. It returns HTTP 500 once
  a request would return too much data, which in practice is around 20 object IDs and
  sometimes fewer, and about 15% of chunks fail at any chunk size. Coverage was also
  partial: roughly 60% of a test list returned data. Fink's own bulk route is their
  separate Data Transfer service, which needs a login and `fink-client`. Prefer `alerce`.
- **IRSA batch forced photometry (ZFPS)** is not covered by this script. It is the
  official ZTF forced photometry and the only route that works at arbitrary positions,
  including objects that never produced an alert. It needs registration, accepts 1500
  positions per submission with 15,000 positions in the system at once, and is queued:
  a published 12,825-position run took about 2.5 days.

A practical plan for a large list: run `alerce` first, which costs half an hour, check how
many targets came back empty, and submit only those to ZFPS.
