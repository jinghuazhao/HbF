# HbF

## Commands

The commands are in sequence. 

File | Description
-----|-------------------------
0_github.sh | GitHub batch file
1_make_variant_list.R | Data preparation
2_lookup.sh | Direct pQTL/Gene lookup
3_region.sh | Regional lookup
4_aggregate.sh | Meta-data
5_coloc.sh<sup>1</sup> | Colocalisation analysis

<sup>1</sup> It is to be done formally for regional associations.

## Lookup

The pQTL/Gene lookup also accumulates results, so we do
```bash
grep '#' 2_lookup.sh | sed 's/#//'
```
upon completion.

## Flanking windows

This is set by `M`, and results for M=1 are contained in the regions/ directory. When M=0, we have another a direct pQTL lookup.

Results are in tsv format of named cohorts.
