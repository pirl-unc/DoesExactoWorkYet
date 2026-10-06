# Frozen osteosarc benchmark inputs

This panel uses osteosarc **0.15.0** and snapshot **2026-09-28**,
`efb65d4b683bda162879c86c0ea889afcc64704b5228224c061109f306503b91`.
Dataset license: CC0-1.0; upstream code: Apache-2.0.

`inputs.json.gz` is the exact corrected comprehensive catalogue input from the
upstream commit and path recorded in `manifest.json`. Its source receipts,
corrections, catalogue entries, and independent vaccine membership assertions
are preserved. `pipeline.fetch_osteosarc` calls osteosarc's membership join to
select the union of all three sources: **51 variants**, all with resolved alleles
in this snapshot. The earlier overlap-only panel contained 37 mutations.

`epitopes.json` freezes the two pVACtools sequence tables fetched with
`Dataset.download` from the same source inventory. The manifest records both
original table checksums and the parsed JSON checksum. These sequence checks
remain separate from vaccine membership assertions.

The five BAM identities were exported from `Dataset.file` in this snapshot.
Their original URLs, size, modification date, index URLs, sample claims and
provenance are retained. Whole-gene windows, grown to contain overlapping
GENCODE transcripts, are acquired with osteosarc's verified regional extraction
API. Small balanced OpenVax fixtures do not substitute for the gene context
needed by the assembly methods.

`pipeline.prepare_inputs` creates each sample's FASTQs and reference once, pins
every required file's checksum, and rejects incomplete or changed preparations.
CI shares those prepared files across methods and caches them by the frozen data
and preparation code. Original regional BAM acquisition receipts are retained
in each extraction JSON. Methods reject input changes, and the scorer rejects
mixed catalogue/read identities. Older scores remain in history but are not
attached to the expanded panel.

Updating this input set is a reviewed change: export a new corrected upstream
catalogue and source inventory, acquire and freeze the matching epitope tables,
update the manifest checksums/version, regenerate the panel, and run the tests.
Do not silently sync live metadata during a scheduled benchmark or site build.

Osteosarc checks source size and stability during extraction. This adapter also
checks modification dates against the frozen inventory; the upstream gap is
tracked in [osteosarc #112](https://github.com/iskandr/osteosarc/issues/112).

The full expanded Exacto matrix is run after this integration is merged. Local
validation covers the complete variant panel and reference windows, plus a real
MAP2 PacBio extraction, cached reuse, and identical FASTQ rebuilds.
