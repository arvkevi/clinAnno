# clinAnno
`clinAnno` is a collection of Python scripts designed to annotate genetic variants in a .vcf file.  `clinAnno` is currently designed to provide annotations for `PS1` and `PM5` variants according to criteria set forth in the [ACMG recommendations](http://www.nature.com/gim/journal/v17/n5/full/gim201530a.html).

1. A variant is annotated with `PS1` when it has the same amino acid change as an established pathogenic variant, regardless of nucleotide change.   
2. A variant is annotated `PM5` when "a novel missense amino acid change occurs at the same position as another pathogenic missense change (e.g., Trp38Ser and Trp38Leu)".

The [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar/) database is used to identify *previously established, pathogenic* variants. 

## Installation
Install using github:
```sh
git clone https://github.com/arvkevi/clinAnno.git
```

## Requirements
1. An annotated .vcf file that contains amino acid change information according to [HGVS](http://www.hgvs.org/mutnomen/examplesAA.html) nomenclature (p.Trp26Cys).  [Variant Effect Predictor](http://www.ensembl.org/Tools/VEP))

## Usage
### clinVar_parser.py -- parse & save clinVar
`clinVar_parser.py` should be executed once, before any .vcf annotations.
Runtime exceeds 15 minutes.

```sh
clinAnno$ python clinVar_parser.py
```

It will save a snapshot of all "pathogenic" and/or "conflicting" SNVs and Indels from clinVar as a
Python dictionary object (pickle file) in the current working directory.

After the initial execution of clinVar_parser.py it's uneccessary to run it
again, until clinVar releases a new update, or whenever you'd like.

### clinAnno.py -- .vcf annotations
After executing `clinVar_parser.py`, the file `clinVar_obj.p` should be in your repo.

Now, annotate a .vcf:
The example, `clinvar.chr1.anno.vcf`, contains chromosome 1 variants from [clinVar's .vcf](ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/).

```sh
clinAnno$ python clinAnno.py --vcf_in=clinvar.chr1.anno.vcf --vcf_out=clinvar.chr1.anno.clinAnno.vcf
```

Processing for large .vcf files can be sped up with the `--nproc` parameter:
```sh
clinAnno$ python clinAnno.py --vcf_in=clinvar.chr1.anno.vcf --vcf_out=clinvar.chr1.anno.clinAnno.vcf --nproc=4
```

`clinAnno` **has been tested on multi-sample and gzipped .vcf's with success.**

## Results
The output of `clinAnno` is a copy of the original .vcf with the following changes:
If no `PS1` or `PM5` variants are found in clinVar, the record (variant) will be unchanged.
Otherwse, the record (variant) will have either `PS1=varid(s)` or `PM5=varid(s)` prepended to the INFO field, with multiple entries separated by `|`.

The first two records from the `clinAnno` annotated .vcf:
```sh
clinAnno$ grep -v '^##' clinvar.chr1.anno.clinAnno.vcf | head -3
```

| #CHROM | POS    | ID          | REF | ALT | QUAL | FILTER | INFO                                                   |
|--------|--------|-------------|-----|-----|------|--------|--------------------------------------------------------|
| 1      | 949523 | rs786201005 | C   | T   | .    | .      | PS1=183381|RS=786201005;RSPOS=949523;dbSNPBuildID=144; |
| 1      | 949696 | rs672601345 | C   | CG  | .    | .      | RS=672601345;RSPOS=949699;dbSNPBuildID=142;            |

^^^Still working on the markdown syntax for the table^^^
The first variant should show additional information after the `PS1` annotation.

Notice the second variant was not annotated with either `PS1=` or `PM5=`, indicating that the variant did not meet the criteria.  *The INFO field has been truncated*

The integer after `PS1=` is clinVar's unique variation identifier.
This directs you to the variant landing page in clinVar:
<http://www.ncbi.nlm.nih.gov/clinvar/variation/183381/>

## Citation
Richards, Sue, et al.
"Standards and guidelines for the interpretation of sequence variants:
a joint consensus recommendation of the American College of
Medical Genetics and Genomics and the Association for Molecular Pathology."
Genetics in Medicine (2015).
