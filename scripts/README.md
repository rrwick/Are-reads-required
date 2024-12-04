# Scripts

At the core of this paper's methods are these scripts which produce a VCF from an assembly and reference, each in a different manner:
* `vcf_from_assembly_shred.sh`: produces a VCF from an assembly using the 'Shred' method.
* `vcf_from_assembly_ska.sh`: produces a VCF from an assembly using the 'SKA' method.
* `vcf_from_assembly_mummer.sh`: produces a VCF from an assembly using the 'MUMmer' method.

These are miscellaneous scripts created for other parts of the paper's methods:
* `rename_contigs.py`: renames contigs in a FASTA file. Used when producing my reference genome assemblies.
* `msa_with_trycycler.sh`: produces a multiple-sequence alignment using [Trycycler msa](https://github.com/rrwick/Trycycler/wiki/Multiple-sequence-alignment).
* `drop_invariant_sites_and_count_diffs.py`: produces an invariant-free MSA and reports total pairwise differences.
* `missing_once_multiple.py`: reports how many bp in an assembly were missing in the reference, once in the reference and multiple times in the reference. Used when investigating the high post-Medaka error counts.
