# HbF

## Commands

The commands are in sequence. 

File | Description
-----|-------------------------
0_github.sh | GitHub batch file
1_make_variant_list.R | Data preparation
2_lookup.sh | Direct pQTL/Gene lookup
3_region.sh | Regional lookup
4_annotate.sh | Meta-data
5_coloc.sh* | Colocalisation analysis

## Lookup

Upon completion we have
```bash
grep '#' 2_lookup.sh | sed 's/#//'
```

## Flanking windows

This is set by `M`.

Results are in tsv format of named cohorts. When M=0, we have another version of direct lookup.

* It is to be done formally for regional associations.
