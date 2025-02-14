# Prepare ONT reads

Tools used:
* [Dorado](https://github.com/nanoporetech/dorado) v0.7.3
* [Samtools](https://github.com/samtools) v1.21

Basecall (on Spartan):
```bash
cd /data/scratch/projects/punim1894/O2024-060

sbatch --job-name=dorado --time=24:00:00 --ntasks=1 --mem=64000 --cpus-per-task=8 -p gpu-a100 --gres=gpu:1 --wrap "/home/rrwick/programs/dorado-0.7.3-linux-x64/bin/dorado basecaller --kit-name SQK-NBD114-96 sup pod5s > reads.bam"

export PATH=/home/rrwick/programs/dorado-0.7.3-linux-x64/bin:"$PATH"
dorado summary reads.bam > summary.tsv
dorado demux --output-dir reads --no-classify reads.bam

cd reads
rm unclassified.bam *_barcode[01234567]*.bam *_barcode80.bam  # unused barcodes
for b in *.bam; do
    samtools fastq -T '*' "$b" | tr '\t' ' ' | paste - - - - | sort | tr '\t' '\n' > ${b/bam/fastq}
done
for f in SQK*.bam SQK*.fastq; do
    new_f=${f/SQK-NBD114-96_/}
    mv "$f" "$new_f"
done
```

Get a qscore distribution:
```bash
{ cat barcode*.fastq | grep -oP "qs:f:[\d\.]+" | grep -oP "[\d\.]+"; echo "0.0"; } | hist -w 0.25 -x -p '0'
```

```
 64099|                                                                                                 0000000000                                                                                               
 60725|                                                                                              000000000000000                                                                                             
 57352|                                                                                           00000000000000000000                                                                                           
 53978|                                                                                        000000000000000000000000                                                                                          
 50605|                                                                                      000000000000000000000000000                                                                                         
 47231|                                                                                    0000000000000000000000000000000                                                                                       
 43857|                                                                                  0000000000000000000000000000000000                                                                                      
 40484|                                                                               00000000000000000000000000000000000000                                                                                     
 37110|                                                                             000000000000000000000000000000000000000000                                                                                   
 33737|                                                                          0000000000000000000000000000000000000000000000                                                                                  
 30363|                                                                        00000000000000000000000000000000000000000000000000                                                                                
 26989|                                                                     0000000000000000000000000000000000000000000000000000000                                                                              
 23616|                                      000                         0000000000000000000000000000000000000000000000000000000000000                                                                           
 20242|                                    000000                    00000000000000000000000000000000000000000000000000000000000000000000                                                                        
 16869|                                   000000000             0000000000000000000000000000000000000000000000000000000000000000000000000000000                                                                  
 13495|                                 0000000000000      0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000                                                        
 10121|                                0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000                                     
  6748|                             0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000                   
  3374|                           00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000        
     1| 0       0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
       ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
       0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 5 
         .   .   .   .   .   .   .   .   .   . 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 0 
         5   5   5   5   5   5   5   5   5   5   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   
                                                 5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   
```

A qscore of 12.5 seems to be a good breakpoint:
```bash
cd /data/scratch/projects/punim1894/O2024-060/reads
for b in barcode*.fastq; do
    pass=${b/.fastq/_pass.fastq}
    cat "$b" | paste - - - - | awk -F'\t' '{split($1, a, " "); for (i in a) {if (a[i] ~ /^qs:f:/) {split(a[i], b, ":"); if (b[3] >= 12.5) print $0}}}' | tr '\t' '\n' > "$pass"
done
```

Zip the pass reads and clean up the pre-QC FASTQs (can just keep the bams):
```bash
cd /data/scratch/projects/punim1894/O2024-060/reads
rm barcode??.fastq
for f in *.fastq; do
    sbatch --job-name=gzip --time=0:30:00 --ntasks=1 --mem=16000 --cpus-per-task=16 --wrap "pigz -9 -p16 $f"
done
```

Move to Roosta:
```bash
mkdir ~/2024-08_SNP_calling_from_assemblies
cd ~/2024-08_SNP_calling_from_assemblies

for s in RES22-01837 RES22-01838 RES22-01839 RES22-01840 RES22-01841 RES22-01842 RES22-01843 RES22-01844 RES22-01845; do
    mkdir "$s"
    mkdir "$s"/reads
    mkdir "$s"/reads_qc
done

cd ~/2024-08_SNP_calling_from_assemblies/RES22-01837/reads
scp spartan:/data/scratch/projects/punim1894/O2024-060/reads/barcode81_pass.fastq.gz nanopore.fastq.gz

cd ~/2024-08_SNP_calling_from_assemblies/RES22-01838/reads
scp spartan:/data/scratch/projects/punim1894/O2024-060/reads/barcode82_pass.fastq.gz nanopore.fastq.gz

cd ~/2024-08_SNP_calling_from_assemblies/RES22-01839/reads
scp spartan:/data/scratch/projects/punim1894/O2024-060/reads/barcode83_pass.fastq.gz nanopore.fastq.gz

cd ~/2024-08_SNP_calling_from_assemblies/RES22-01840/reads
scp spartan:/data/scratch/projects/punim1894/O2024-060/reads/barcode84_pass.fastq.gz nanopore.fastq.gz

cd ~/2024-08_SNP_calling_from_assemblies/RES22-01841/reads
scp spartan:/data/scratch/projects/punim1894/O2024-060/reads/barcode85_pass.fastq.gz nanopore.fastq.gz

cd ~/2024-08_SNP_calling_from_assemblies/RES22-01842/reads
scp spartan:/data/scratch/projects/punim1894/O2024-060/reads/barcode86_pass.fastq.gz nanopore.fastq.gz

cd ~/2024-08_SNP_calling_from_assemblies/RES22-01843/reads
scp spartan:/data/scratch/projects/punim1894/O2024-060/reads/barcode87_pass.fastq.gz nanopore.fastq.gz

cd ~/2024-08_SNP_calling_from_assemblies/RES22-01844/reads
scp spartan:/data/scratch/projects/punim1894/O2024-060/reads/barcode88_pass.fastq.gz nanopore.fastq.gz

cd ~/2024-08_SNP_calling_from_assemblies/RES22-01845/reads
scp spartan:/data/scratch/projects/punim1894/O2024-060/reads/barcode89_pass.fastq.gz nanopore.fastq.gz
```



# Prepare Illumina reads

```bash
cd ~/2024-08_SNP_calling_from_assemblies/RES22-01837/reads
cp /home/damg/data/I2024-018/DMG2404030_*_R1_*.fastq.gz illumina_1.fastq.gz
cp /home/damg/data/I2024-018/DMG2404030_*_R2_*.fastq.gz illumina_2.fastq.gz

cd ~/2024-08_SNP_calling_from_assemblies/RES22-01838/reads
cp /home/damg/data/I2024-018/DMG2404031_*_R1_*.fastq.gz illumina_1.fastq.gz
cp /home/damg/data/I2024-018/DMG2404031_*_R2_*.fastq.gz illumina_2.fastq.gz

cd ~/2024-08_SNP_calling_from_assemblies/RES22-01839/reads
cp /home/damg/data/I2024-018/DMG2404032_*_R1_*.fastq.gz illumina_1.fastq.gz
cp /home/damg/data/I2024-018/DMG2404032_*_R2_*.fastq.gz illumina_2.fastq.gz

cd ~/2024-08_SNP_calling_from_assemblies/RES22-01840/reads
cp /home/damg/data/I2024-018/DMG2404033_*_R1_*.fastq.gz illumina_1.fastq.gz
cp /home/damg/data/I2024-018/DMG2404033_*_R2_*.fastq.gz illumina_2.fastq.gz

cd ~/2024-08_SNP_calling_from_assemblies/RES22-01841/reads
cp /home/damg/for_transfer/I2024-022/DMG2404034_*_R1_*.fastq.gz illumina_1.fastq.gz
cp /home/damg/for_transfer/I2024-022/DMG2404034_*_R2_*.fastq.gz illumina_2.fastq.gz

cd ~/2024-08_SNP_calling_from_assemblies/RES22-01842/reads
cp /home/damg/data/I2024-018/DMG2404035_*_R1_*.fastq.gz illumina_1.fastq.gz
cp /home/damg/data/I2024-018/DMG2404035_*_R2_*.fastq.gz illumina_2.fastq.gz

cd ~/2024-08_SNP_calling_from_assemblies/RES22-01843/reads
cp /home/damg/data/I2024-018/DMG2404036_*_R1_*.fastq.gz illumina_1.fastq.gz
cp /home/damg/data/I2024-018/DMG2404036_*_R2_*.fastq.gz illumina_2.fastq.gz

cd ~/2024-08_SNP_calling_from_assemblies/RES22-01844/reads
cp /home/damg/for_transfer/I2024-022/DMG2404037_*_R1_*.fastq.gz illumina_1.fastq.gz
cp /home/damg/for_transfer/I2024-022/DMG2404037_*_R2_*.fastq.gz illumina_2.fastq.gz

cd ~/2024-08_SNP_calling_from_assemblies/RES22-01845/reads
cp /home/damg/data/I2024-018/DMG2404038_*_R1_*.fastq.gz illumina_1.fastq.gz
cp /home/damg/data/I2024-018/DMG2404038_*_R2_*.fastq.gz illumina_2.fastq.gz

cd ~/2024-08_SNP_calling_from_assemblies
chmod 664 RES*/reads/*.fastq.gz
```



# Read QC

Tools used:
* [fastp](https://github.com/OpenGene/fastp) v0.23.4
* [Filtlong](https://github.com/rrwick/Filtlong) v0.2.1

```bash
cd ~/2024-08_SNP_calling_from_assemblies
for s in RES*; do
    cd ~/2024-08_SNP_calling_from_assemblies/"$s"/reads_qc
    fastp --in1 ../reads/illumina_1.fastq.gz --in2 ../reads/illumina_2.fastq.gz --out1 illumina_1.fastq.gz --out2 illumina_2.fastq.gz --unpaired1 illumina_u.fastq.gz --unpaired2 illumina_u.fastq.gz
done

cd ~/2024-08_SNP_calling_from_assemblies
for s in RES*; do
    cd ~/2024-08_SNP_calling_from_assemblies/"$s"/reads_qc
    filtlong --min_length 1000 ../reads/nanopore.fastq.gz | pigz -p8 > nanopore.fastq.gz
done
```



# Unicycler assemblies

Tools used:
* [Unicycler](https://github.com/rrwick/Unicycler) v0.5.1
* [SPAdes](https://github.com/ablab/spades) v4.0.0

Unicycler Illumina-only:
```bash
cd ~/2024-08_SNP_calling_from_assemblies
for s in RES*; do
    cd ~/2024-08_SNP_calling_from_assemblies/"$s"
    unicycler -1 reads_qc/illumina_1.fastq.gz -2 reads_qc/illumina_2.fastq.gz -o unicycler_illumina -t 32
done
```

Unicycler hybrid:
```bash
cd ~/2024-08_SNP_calling_from_assemblies
for s in RES*; do
    cd ~/2024-08_SNP_calling_from_assemblies/"$s"
    mkdir unicycler_hybrid
    cp unicycler_illumina/002_depth_filter.gfa unicycler_hybrid
done

cd ~/2024-08_SNP_calling_from_assemblies
for s in RES*; do
    cd ~/2024-08_SNP_calling_from_assemblies/"$s"
    unicycler -1 reads_qc/illumina_1.fastq.gz -2 reads_qc/illumina_2.fastq.gz -l reads_qc/nanopore.fastq.gz -o unicycler_hybrid -t 64
done
```

Note: I'm probably not going to use these Unicycler assemblies in making my reference for each genome, but they can be a useful comparison point, e.g. for small plasmids that might be troublesome in a long-read-first assembly.



# Reference assemblies

Tools used:
* [Canu](https://github.com/marbl/canu) v2.2
* [Flye](https://github.com/mikolmogorov/Flye) v2.9.5
* [miniasm](https://github.com/lh3/miniasm) v0.3
* [Minipolish](https://github.com/rrwick/Minipolish) v0.1.3
* [Racon](https://github.com/lbcb-sci/racon) v1.5.0
* [any2fasta](https://github.com/tseemann/any2fasta) v0.4.2
* [NECAT](https://github.com/xiaochuanle/NECAT) v20200803
* [NextDenovo](https://github.com/Nextomics/NextDenovo) v2.5.2
* [NextPolish](https://github.com/Nextomics/NextPolish) v1.4.1
* [Raven](https://github.com/lbcb-sci/raven) v1.8.3
* [Trycycler](https://github.com/rrwick/Trycycler) v0.5.5
* [Polypolish](https://github.com/rrwick/Polypolish) v0.6.0
* [Pypolca](https://github.com/gbouras13/pypolca) v0.3.1
* [Dnaapler](https://github.com/gbouras13/dnaapler) v0.8.1
* [Sniffles](https://github.com/fritzsedlazeck/Sniffles) v2.4
* [Seqtk](https://github.com/lh3/seqtk) v1.4
* [Samtools](https://github.com/samtools) v1.21
* [Bcftools](https://github.com/samtools/bcftools) v1.21

Generate input assemblies using [Trycycler's extra-thorough method](https://github.com/rrwick/Trycycler/wiki/Generating-assemblies#extra-thorough-assembly):
```bash
genome_size=2880000

cd ~/2024-08_SNP_calling_from_assemblies
for s in RES*; do
    
    cd ~/2024-08_SNP_calling_from_assemblies/"$s"

    trycycler subsample --reads reads_qc/nanopore.fastq.gz --out_dir read_subsets --count 24 --genome_size "$genome_size"
    mkdir assemblies

    for i in 01 07 13 19; do
        canu -p canu -d canu_temp -fast genomeSize="$genome_size" useGrid=false maxThreads="$threads" -nanopore read_subsets/sample_"$i".fastq
        canu_trim.py canu_temp/canu.contigs.fasta > assemblies/assembly_"$i".fasta
        rm -rf canu_temp
    done
    for i in 02 08 14 20; do
        flye --nano-hq read_subsets/sample_"$i".fastq --threads 32 --out-dir flye_temp
        cp flye_temp/assembly.fasta assemblies/assembly_"$i".fasta
        cp flye_temp/assembly_graph.gfa assemblies/assembly_"$i".gfa
        rm -r flye_temp
    done

    for i in 03 09 15 21; do
        miniasm_and_minipolish.sh read_subsets/sample_"$i".fastq 32 > assemblies/assembly_"$i".gfa
        any2fasta assemblies/assembly_"$i".gfa > assemblies/assembly_"$i".fasta
    done

    for i in 04 10 16 22; do
        ~/programs/NECAT/Linux-amd64/bin/necat.pl config config.txt
        realpath read_subsets/sample_"$i".fastq > read_list.txt
        sed -i "s/PROJECT=/PROJECT=necat/" config.txt
        sed -i "s/ONT_READ_LIST=/ONT_READ_LIST=read_list.txt/" config.txt
        sed -i "s/GENOME_SIZE=/GENOME_SIZE="$genome_size"/" config.txt
        sed -i "s/THREADS=4/THREADS=12/" config.txt
        ~/programs/NECAT/Linux-amd64/bin/necat.pl bridge config.txt
        cp necat/6-bridge_contigs/polished_contigs.fasta assemblies/assembly_"$i".fasta
        rm -r necat config.txt read_list.txt
    done

    for i in 05 11 17 23; do
        echo read_subsets/sample_"$i".fastq > input.fofn
        cp ~/programs/NextDenovo/doc/run.cfg nextdenovo_run.cfg
        sed -i "s/genome_size = 1g/genome_size = "$genome_size"/" nextdenovo_run.cfg
        sed -i "s/parallel_jobs = 20/parallel_jobs = 1/" nextdenovo_run.cfg
        sed -i "s/read_type = clr/read_type = ont/" nextdenovo_run.cfg
        sed -i "s/pa_correction = 3/pa_correction = 1/" nextdenovo_run.cfg
        sed -i "s/correction_options = -p 15/correction_options = -p 32/" nextdenovo_run.cfg
        sed -i "s/-t 8/-t 32/" nextdenovo_run.cfg
        nextDenovo nextdenovo_run.cfg
        cp 01_rundir/03.ctg_graph/nd.asm.fasta nextdenovo_temp.fasta
        rm -r 01_rundir nextdenovo_run.cfg input.fofn
        echo read_subsets/sample_"$i".fastq > lgs.fofn
        cat ~/programs/NextPolish/doc/run.cfg | grep -v "sgs" | grep -v "hifi" > nextpolish_run.cfg
        sed -i "s/parallel_jobs = 6/parallel_jobs = 1/" nextpolish_run.cfg
        sed -i "s/multithread_jobs = 5/multithread_jobs = 32/" nextpolish_run.cfg
        sed -i "s|genome = ./raw.genome.fasta|genome = nextdenovo_temp.fasta|" nextpolish_run.cfg
        sed -i "s|-x map-ont|-x map-ont -t 32|" nextpolish_run.cfg
        nextPolish nextpolish_run.cfg
        cp 01_rundir/genome.nextpolish.fasta assemblies/assembly_"$i".fasta
        rm -r 01_rundir pid*.log.info nextpolish_run.cfg lgs.fofn nextdenovo_temp.fasta
    done

    for i in 06 12 18 24; do
        raven --threads 32 --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_"$i".gfa read_subsets/sample_"$i".fastq > assemblies/assembly_"$i".fasta
    done
done

cd ~/2024-08_SNP_calling_from_assemblies
rm -r RES*/read_subsets
```

Trycycler cluster:
```bash
cd ~/2024-08_SNP_calling_from_assemblies
for s in RES*; do
    cd ~/2024-08_SNP_calling_from_assemblies/"$s"
    trycycler cluster --assemblies assemblies/*.fasta --reads reads_qc/nanopore.fastq.gz --out_dir trycycler --threads 64
done
```

And then Trycycler reconcile:
```bash
trycycler reconcile --reads reads_qc/nanopore.fastq.gz --cluster_dir trycycler/cluster_xxx
```

And the remaining Trycycler steps:
```bash
trycycler msa --threads 64 --cluster_dir trycycler/cluster_xxx
trycycler partition --reads reads_qc/nanopore.fastq.gz --cluster_dirs trycycler/cluster_* --threads 96
for c in trycycler/cluster_*; do
    trycycler consensus --cluster_dir "$c"
done
cat trycycler/cluster_*/7_final_consensus.fasta > trycycler.fasta
```

Illumina polishing:
```bash
cd ~/2024-08_SNP_calling_from_assemblies
for s in RES*; do
    cd ~/2024-08_SNP_calling_from_assemblies/"$s"
    bwa index trycycler.fasta
    bwa mem -t 24 -a trycycler.fasta reads_qc/illumina_1.fastq.gz > alignments_1.sam
    bwa mem -t 24 -a trycycler.fasta reads_qc/illumina_2.fastq.gz > alignments_2.sam
    polypolish filter --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam
    polypolish polish trycycler.fasta filtered_1.sam filtered_2.sam > trycycler_polypolish.fasta
    rm *.amb *.ann *.bwt *.pac *.sa *.sam
    pypolca run --careful -a trycycler_polypolish.fasta -1 reads_qc/illumina_1.fastq.gz -2 reads_qc/illumina_2.fastq.gz -t 24 -o pypolca
    cp pypolca/pypolca_corrected.fasta trycycler_polypolish_pypolca.fasta
    rm -r pypolca
done
```

Rotate to starting gene and rename contigs:
```bash
cd ~/2024-08_SNP_calling_from_assemblies
for s in RES*; do
    dnaapler all -i "$s"/trycycler_polypolish_pypolca.fasta -o temp -p "$s"
    seqtk seq temp/"$s"_reoriented.fasta > "$s"/reference.fasta
    scripts/rename_contigs.py "$s"/reference.fasta
    rm -r temp
done
```
Dnaapler didn't always rotate the plasmids (not sure why), so I had to manually rotate a few to match the rest.

Check for SVs with Sniffles:
```bash
cd ~/2024-08_SNP_calling_from_assemblies
for s in RES*; do
    cd ~/2024-08_SNP_calling_from_assemblies/"$s"
    minimap2 -a -x map-ont -t 48 reference.fasta reads_qc/nanopore.fastq.gz | samtools sort > nanopore.bam
    samtools index nanopore.bam
    sniffles -i nanopore.bam -v "$s".vcf
    rm nanopore.bam*
done
```

Check for SNPs with Clair3:
```bash
conda activate clair3

cd ~/2024-08_SNP_calling_from_assemblies
model=/home/wickr/programs/rerio/clair3_models/r1041_e82_400bps_sup_v500
for s in RES*; do
    cd ~/2024-08_SNP_calling_from_assemblies/"$s"
    minimap2 -a -x map-ont -t 32 reference.fasta reads_qc/nanopore.fastq.gz | samtools sort > nanopore.bam
    samtools index nanopore.bam
    samtools faidx reference.fasta
    run_clair3.sh --bam_fn=nanopore.bam --ref_fn=reference.fasta --threads=32 --platform="ont" --model_path="$model" --output=clair3_reference --sample_name="$s" --include_all_ctgs --haploid_precise --no_phasing_for_fa --enable_long_indel
    gunzip -c clair3_reference/merge_output.vcf.gz | bcftools view -f 'PASS,.' | bcftools view -e 'GT="0/1"' > clair3_reference.vcf
    rm nanopore.bam* reference.fasta.fai
done
```

Check for SNPs with FreeBayes:
```bash
conda activate snippy

cd ~/2024-08_SNP_calling_from_assemblies
for s in RES*; do
    cd ~/2024-08_SNP_calling_from_assemblies/"$s"
    bwa index reference.fasta
    bwa mem -t 32 reference.fasta reads_qc/illumina_1.fastq.gz reads_qc/illumina_2.fastq.gz | samtools sort > illumina.bam
    samtools index illumina.bam
    samtools faidx reference.fasta
    freebayes -f reference.fasta --haplotype-length -1 -m 10 -q 10 -p 1 --min-coverage 2 illumina.bam | bcftools view -i 'QUAL>=100' | bcftools reheader -s <(echo "unknown\t$s") > freebayes_reference.vcf
    rm illumina.bam* reference.fasta.*
done
```

Found nothing, so I'm confident these assemblies are rock-solid!
```bash
cd ~/2024-08_SNP_calling_from_assemblies
rm RES*/*.vcf
```

Tidy-up files for making and checking the reference sequence:
```bash
cd ~/2024-08_SNP_calling_from_assemblies
rm RES*/trycycler/*cluster*/4_reads.fastq
for s in RES*; do
    cd ~/2024-08_SNP_calling_from_assemblies/"$s"
    mkdir reference
    mv assemblies trycycler* unicycler* reference
    rm -r clair3_reference clair3_reference.vcf freebayes_reference.vcf
done
```



# Read subsampling

Tools used:
* [Rasusa](https://github.com/mbhall88/rasusa) v2.1.0

For both ONT and Illumina reads, I'd subsampling the QCed reads to various depths for assembly: 20x, 50x and 100x. All genomes are essentially the same size (with 1 bp of 2886598 bp), so I can use the same genome size parameter for each.

```bash
conda activate rasusa

cd ~/2024-08_SNP_calling_from_assemblies
for s in RES*; do
    cd ~/2024-08_SNP_calling_from_assemblies/"$s"
    mkdir reads_subsampled

    rasusa reads --coverage 20 --genome-size 2886598 -o reads_subsampled/nanopore_020x.fastq reads_qc/nanopore.fastq.gz
    rasusa reads --coverage 50 --genome-size 2886598 -o reads_subsampled/nanopore_050x.fastq reads_qc/nanopore.fastq.gz
    rasusa reads --coverage 100 --genome-size 2886598 -o reads_subsampled/nanopore_100x.fastq reads_qc/nanopore.fastq.gz

    rasusa reads --coverage 20 --genome-size 2886598 -o reads_subsampled/illumina_020x_1.fastq -o reads_subsampled/illumina_020x_2.fastq reads_qc/illumina_1.fastq.gz reads_qc/illumina_2.fastq.gz
    rasusa reads --coverage 50 --genome-size 2886598 -o reads_subsampled/illumina_050x_1.fastq -o reads_subsampled/illumina_050x_2.fastq reads_qc/illumina_1.fastq.gz reads_qc/illumina_2.fastq.gz
    rasusa reads --coverage 100 --genome-size 2886598 -o reads_subsampled/illumina_100x_1.fastq -o reads_subsampled/illumina_100x_2.fastq reads_qc/illumina_1.fastq.gz reads_qc/illumina_2.fastq.gz
done

cd ~/2024-08_SNP_calling_from_assemblies
pigz -p16 RES*/reads_subsampled/*.fastq
```



# Unicycler (short-read) assembly

Tools used:
* [Unicycler](https://github.com/rrwick/Unicycler) v0.5.1
* [SPAdes](https://github.com/ablab/spades) v4.0.0

```bash
conda activate unicycler

cd ~/2024-08_SNP_calling_from_assemblies
for s in RES*; do
    cd ~/2024-08_SNP_calling_from_assemblies/"$s"
    for d in 020 050 100; do
        r1=reads_subsampled/illumina_"$d"x_1.fastq.gz
        r2=reads_subsampled/illumina_"$d"x_2.fastq.gz
        unicycler -1 "$r1" -2 "$r2" -o unicycler_short_"$d"x -t 32
        cp unicycler_short_"$d"x/assembly.fasta unicycler_short_"$d"x.fasta
    done
done
```



# Shovill (short-read) assembly

Tools used:
* [Shovill](https://github.com/tseemann/shovill) v1.1.0
* [SPAdes](https://github.com/ablab/spades) v4.0.0

```bash
conda activate shovill

cd ~/2024-08_SNP_calling_from_assemblies
for s in RES*; do
    cd ~/2024-08_SNP_calling_from_assemblies/"$s"
    for d in 020 050 100; do
        r1=reads_subsampled/illumina_"$d"x_1.fastq.gz
        r2=reads_subsampled/illumina_"$d"x_2.fastq.gz
        shovill --outdir shovill_"$d"x --R1 "$r1" --R2 "$r2" --cpus 32
        cp shovill_"$d"x/contigs.fa shovill_"$d"x.fasta
    done
done
```



# SKESA (short-read) assembly

Tools used:
* [SKESA](https://github.com/ncbi/SKESA) v2.5.1

```bash
cd ~/2024-08_SNP_calling_from_assemblies
for s in RES*; do
    cd ~/2024-08_SNP_calling_from_assemblies/"$s"
    for d in 020 050 100; do
        r1=reads_subsampled/illumina_"$d"x_1.fastq.gz
        r2=reads_subsampled/illumina_"$d"x_2.fastq.gz
        skesa --reads "$r1","$r2" --cores 32 --memory 128 > skesa_"$d"x.fasta
    done
done
```



# Flye (long-read) assembly

Tools used:
* [Flye](https://github.com/mikolmogorov/Flye) v2.9.5

```bash
conda activate long-read-assemblers

cd ~/2024-08_SNP_calling_from_assemblies
for s in RES*; do
    cd ~/2024-08_SNP_calling_from_assemblies/"$s"
    for d in 020 050 100; do
        ont=reads_subsampled/nanopore_"$d"x.fastq.gz
        flye --nano-hq "$ont" --threads 32 --out-dir flye_"$d"x
        cp flye_"$d"x/assembly.fasta flye_"$d"x.fasta
    done
done
```



# Canu (long-read) assembly

Tools used:
* [Canu](https://github.com/marbl/canu) v2.2

```bash
conda activate long-read-assemblers

cd ~/2024-08_SNP_calling_from_assemblies
for s in RES*; do
    cd ~/2024-08_SNP_calling_from_assemblies/"$s"
    for d in 020 050 100; do
        ont=reads_subsampled/nanopore_"$d"x.fastq.gz
        canu -p canu -d canu_"$d"x -fast genomeSize=2886598 useGrid=false maxThreads=32 -nanopore "$ont"
        cp canu_"$d"x/canu.contigs.fasta canu_"$d"x.fasta
    done
done
```



# Raven (long-read) assembly

Tools used:
* [Raven](https://github.com/lbcb-sci/raven) v1.8.3

```bash
conda activate long-read-assemblers

cd ~/2024-08_SNP_calling_from_assemblies
for s in RES*; do
    cd ~/2024-08_SNP_calling_from_assemblies/"$s"
    for d in 020 050 100; do
        ont=reads_subsampled/nanopore_"$d"x.fastq.gz
        raven --threads 32 --disable-checkpoints "$ont" > raven_"$d"x.fasta
    done
done
```



# Medaka

Tools used:
* [Medaka](https://github.com/nanoporetech/medaka)  v2.0.0

```bash
conda activate medaka

cd ~/2024-08_SNP_calling_from_assemblies
for s in RES*; do
    cd ~/2024-08_SNP_calling_from_assemblies/"$s"
    for d in 020 050 100; do
        ont=reads_subsampled/nanopore_"$d"x.fastq.gz
        for a in canu flye raven; do
            medaka_consensus -i "$ont" -d "$a"_"$d"x.fasta -o medaka -t 32 --bacteria
            cp medaka/consensus.fasta "$a"_medaka_"$d"x.fasta
            rm -r medaka
            rm "$a"_"$d"x.fasta.fai "$a"_"$d"x.fasta.map-ont.mmi
        done
    done
done
```



# Unicycler (hybrid) assembly

Tools used:
* [Unicycler](https://github.com/rrwick/Unicycler) v0.5.1
* [SPAdes](https://github.com/ablab/spades) v4.0.0

```bash
conda activate unicycler

cd ~/2024-08_SNP_calling_from_assemblies
for s in RES*; do
    cd ~/2024-08_SNP_calling_from_assemblies/"$s"
    for d in 020 050 100; do
        r1=reads_subsampled/illumina_"$d"x_1.fastq.gz
        r2=reads_subsampled/illumina_"$d"x_2.fastq.gz
        ont=reads_subsampled/nanopore_"$d"x.fastq.gz
        unicycler -1 "$r1" -2 "$r2" -l "$ont" -o unicycler_hybrid_"$d"x -t 32
        cp unicycler_hybrid_"$d"x/assembly.fasta unicycler_hybrid_"$d"x.fasta
    done
done
```



# Hybracter (hybrid) assembly

Tools used:
* [Hybracter](https://github.com/gbouras13/hybracter) v0.9.0

```bash
conda activate hybracter

cd ~/2024-08_SNP_calling_from_assemblies
for s in RES*; do
    cd ~/2024-08_SNP_calling_from_assemblies/"$s"
    for d in 020 050 100; do
        r1=reads_subsampled/illumina_"$d"x_1.fastq.gz
        r2=reads_subsampled/illumina_"$d"x_2.fastq.gz
        ont=reads_subsampled/nanopore_"$d"x.fastq.gz
        hybracter hybrid-single -l "$ont" -1 "$r1" -2 "$r2" -s "$s" -c 2800000 -o hybracter_hybrid_"$d"x -t 32
        cp hybracter_hybrid_"$d"x/FINAL_OUTPUT/*complete/"$s"_final.fasta hybracter_hybrid_"$d"x.fasta
    done
done
```



# Prepare reference

I'm using RES22-01844 as the reference since it's the NRS384 wild type genome.

```bash
cd ~/2024-08_SNP_calling_from_assemblies
cp RES22-01844/reference.fasta ref.fasta
samtools faidx ref.fasta
bwa index ref.fasta
```



# Call SNPs from ONT reads

Tools used:
* [Clair3](https://github.com/HKU-BAL/Clair3) v1.0.10
* [Bcftools](https://github.com/samtools/bcftools) v1.21

Notes:
* I based my Clair3 command on what Michael used [here](https://github.com/mbhall88/NanoVarBench/blob/main/workflow/scripts/callers/clair3.sh).
* I only take variants with a `FILTER` value of `PASS` and I exclude any het calls (genotype of `0/1`). That second check shouldn't be necessary with Clair3 in haploid mode, but I noticed a strange behaviour: if no variants are found, it copies an earlier file to final results, which still includes the het calls.

```bash
conda activate clair3

cd ~/2024-08_SNP_calling_from_assemblies
ref=/home/wickr/2024-08_SNP_calling_from_assemblies/ref.fasta
model=/home/wickr/programs/rerio/clair3_models/r1041_e82_400bps_sup_v500
for s in RES*; do
    cd ~/2024-08_SNP_calling_from_assemblies/"$s"
    for d in 020 050 100; do
        ont=reads_subsampled/nanopore_"$d"x.fastq.gz
        minimap2 -a -x map-ont -t 32 "$ref" "$ont" | samtools sort > nanopore.bam
        samtools index nanopore.bam
        run_clair3.sh --bam_fn=nanopore.bam --ref_fn="$ref" --threads=32 --platform="ont" --model_path="$model" --output=clair3_"$d" --sample_name="$s" --include_all_ctgs --haploid_precise --no_phasing_for_fa --enable_long_indel
        gunzip -c clair3_"$d"/merge_output.vcf.gz | bcftools view -f 'PASS,.' | bcftools view -e 'GT="0/1"' > nanopore_"$d"x.vcf
        rm nanopore.bam*
    done
done
```



# Call SNPs from Illumina reads

Tools used:
* [BWA](https://github.com/lh3/bwa) v0.7.18
* [Freebayes](https://github.com/freebayes/freebayes) v1.3.8
* [Samtools](https://github.com/samtools) v1.21
* [Bcftools](https://github.com/samtools/bcftools) v1.21

Notes:
* I based my `freebayes` command using what Michael used [here](https://github.com/mbhall88/NanoVarBench/blob/main/workflow/scripts/callers/freebayes.sh).
* I use a minimum quality score of 100, which is what Snippy uses (via its `--minqual` parameter).
* Freebayes didn't seem to have a way to include a sample name (it just used `unknown`) so I needed to use `bcftools reheader` to add the sample name to the VCF (important for the merge step later).

```bash
conda activate snippy

cd ~/2024-08_SNP_calling_from_assemblies
ref=/home/wickr/2024-08_SNP_calling_from_assemblies/ref.fasta
for s in RES*; do
    cd ~/2024-08_SNP_calling_from_assemblies/"$s"
    for d in 020 050 100; do
        r1=reads_subsampled/illumina_"$d"x_1.fastq.gz
        r2=reads_subsampled/illumina_"$d"x_2.fastq.gz
        bwa mem -t 32 "$ref" "$r1" "$r2" | samtools sort > illumina.bam
        samtools index illumina.bam
        freebayes -f "$ref" --haplotype-length -1 -m 10 -q 10 -p 1 --min-coverage 2 illumina.bam | bcftools view -i 'QUAL>=100' | bcftools reheader -s <(echo "unknown\t$s") > illumina_"$d"x.vcf
        rm illumina.bam*
    done
done
```



# Create VCFs from assemblies

Tools used:
* [all2vcf](https://github.com/rrwick/all2vcf) v0.7.8

When I first used all2vcf to convert from MUMmer's `.snp` file to a `.vcf`, I encountered some errors with indels. I therefore had to fix the tool (see [this pull request](https://github.com/MatteoSchiavinato/all2vcf/pull/7)), and this code used my fixed all2vcf.

Using my reference assemblies:
```bash
conda activate trycycler
export PATH=/home/wickr/programs/all2vcf:"$PATH"

cd ~/2024-08_SNP_calling_from_assemblies
for s in RES*; do
    a="$s"/reference
    scripts/vcf_from_assembly_shred.sh ref.fasta "$a".fasta 32 > "$a"_shred.vcf
    scripts/vcf_from_assembly_ska.sh ref.fasta "$a".fasta > "$a"_ska.vcf
    scripts/vcf_from_assembly_mummer.sh ref.fasta "$a".fasta > "$a"_mummer.vcf
done
```

Using other assemblies:
```bash
conda activate trycycler
export PATH=/home/wickr/programs/all2vcf:"$PATH"

cd ~/2024-08_SNP_calling_from_assemblies
for assembler in unicycler_short shovill skesa flye canu raven flye_medaka canu_medaka raven_medaka unicycler_hybrid hybracter_hybrid; do
    for d in 020 050 100; do
        for s in RES*; do
            a="$s"/"$assembler"_"$d"x
            scripts/vcf_from_assembly_shred.sh ref.fasta "$a".fasta 32 > "$a"_shred.vcf
            scripts/vcf_from_assembly_ska.sh ref.fasta "$a".fasta > "$a"_ska.vcf
            scripts/vcf_from_assembly_mummer.sh ref.fasta "$a".fasta > "$a"_mummer.vcf
        done
    done
done
```



# Create whole-genome MSAs

Tools used:
* [Bcftools](https://github.com/samtools/bcftools) v1.21
* [Trycycler](https://github.com/rrwick/Trycycler) v0.5.5

```bash
cd ~/2024-08_SNP_calling_from_assemblies
mkdir alignments
```

Run bgzip and tabix (for bcftools):
```bash
cd ~/2024-08_SNP_calling_from_assemblies
for v in RES*/*.vcf; do
    bgzip "$v"
    tabix "$v".gz
done
```

Using Illumina reads:
```bash
conda activate trycycler

cd ~/2024-08_SNP_calling_from_assemblies
for d in 020 050 100; do
    rm -f unaligned.fasta
    for s in RES*; do
        echo ">$s" >> unaligned.fasta
        bcftools consensus -f ref.fasta "$s"/illumina_"$d"x.vcf.gz | seqtk seq | grep -vP "^>" | tr -d '\n' >> unaligned.fasta
        printf "\n" >> unaligned.fasta
    done
    scripts/msa_with_trycycler.sh unaligned.fasta 32 > alignments/illumina_"$d".fasta
    rm unaligned.fasta
done
```

Using ONT reads:
```bash
conda activate trycycler

cd ~/2024-08_SNP_calling_from_assemblies
for d in 020 050 100; do
    rm -f unaligned.fasta
    for s in RES*; do
        echo ">$s" >> unaligned.fasta
        bcftools consensus -f ref.fasta "$s"/nanopore_"$d"x.vcf.gz | seqtk seq | grep -vP "^>" | tr -d '\n' >> unaligned.fasta
        printf "\n" >> unaligned.fasta
    done
    scripts/msa_with_trycycler.sh unaligned.fasta 32 > alignments/nanopore_"$d".fasta
    rm unaligned.fasta
done
```

Using my reference assemblies (direct):
```bash
conda activate trycycler

cd ~/2024-08_SNP_calling_from_assemblies
rm -f unaligned.fasta
for s in RES*; do
    echo ">$s" >> unaligned.fasta
    cat "$s"/reference.fasta | seqtk seq | grep -vP "^>" | tr -d '\n' >> unaligned.fasta
    printf "\n" >> unaligned.fasta
done
scripts/msa_with_trycycler.sh unaligned.fasta 32 > alignments/reference_direct.fasta
rm unaligned.fasta
```

Using my reference assemblies (from VCF):
```bash
conda activate trycycler

cd ~/2024-08_SNP_calling_from_assemblies
for m in shred ska mummer; do
    rm -f unaligned.fasta
    for s in RES*; do
        a="$s"/reference
        echo ">$s" >> unaligned.fasta
        bcftools consensus -f ref.fasta "$a"_"$m".vcf.gz | seqtk seq | grep -vP "^>" | tr -d '\n' >> unaligned.fasta
        printf "\n" >> unaligned.fasta
    done
    scripts/msa_with_trycycler.sh unaligned.fasta 32 > alignments/reference_"$m".fasta
    rm unaligned.fasta
done
```

Using other assemblies:
```bash
conda activate trycycler

cd ~/2024-08_SNP_calling_from_assemblies
for assembler in unicycler_short shovill skesa flye canu raven flye_medaka canu_medaka raven_medaka unicycler_hybrid hybracter_hybrid; do
    for d in 020 050 100; do
        for m in shred ska mummer; do
            rm -f unaligned.fasta
            for s in RES*; do
                a="$s"/"$assembler"_"$d"x
                echo ">$s" >> unaligned.fasta
                bcftools consensus -f ref.fasta "$a"_"$m".vcf.gz | seqtk seq | grep -vP "^>" | tr -d '\n' >> unaligned.fasta
                printf "\n" >> unaligned.fasta
            done
            scripts/msa_with_trycycler.sh unaligned.fasta 32 > alignments/"$assembler"_"$d"x_"$m".fasta
            rm unaligned.fasta
        done
    done
done
```

Clean up:
```bash
cd ~/2024-08_SNP_calling_from_assemblies
rm RES*/*.vcf.gz.tbi
gunzip RES*/*.vcf.gz
```



# Establishing truth

```bash
cd ~/2024-08_SNP_calling_from_assemblies/alignments
md5sum *.fasta | sort
```

The following methods all produced the exact same alignments:
* Reference assemblies (direct, shred, MUMmer)
* Illumina 20x
* Illumina 50x
* Illumina 100x
* ONT 100x
* Hybracter 20x (shred, MUMmer)
* Hybracter 50x (shred, MUMmer)
* Hybracter 100x (shred, MUMmer)

This matches my expectations: reference assemblies are solid, Hybracter assemblies are solid, Illumina variant calling is solid and Nanopore variant calling is good (at least at higher depths). So I'm confident in considering these to be the truth.

```bash
cd ~/2024-08_SNP_calling_from_assemblies
for s in RES*; do
    cd ~/2024-08_SNP_calling_from_assemblies/"$s"
    cp reference_shred.vcf truth.vcf
done
```



# Build tree

Tools used:
* [IQ-TREE](https://github.com/iqtree/iqtree2) v2.3.6

I used the alignment based directly on reference assemblies to generate the tree, since this was consistent with my truth variants (see above).

```bash
conda activate iqtree

cd ~/2024-08_SNP_calling_from_assemblies/alignments
../scripts/drop_invariant_sites_and_count_diffs.py reference_direct.fasta > reference_direct_no_invariant.fasta



# Manually replace gaps with a base (A) so IQ-TREE uses them:
cat reference_direct_no_invariant.fasta | sed 's/-/A/g' | sed 's/2A0/2-0/' > reference_direct_no_invariant_bases_for_gaps.fasta

iqtree2 -s reference_direct_no_invariant_bases_for_gaps.fasta --polytomy
```



# Counting errors

Tools used:
* [vcfdist](https://github.com/TimD1/vcfdist) v2.5.3


Run vcfdist on each VCF against the truth:
```bash
cd ~/2024-08_SNP_calling_from_assemblies
for s in RES*; do
    cd ~/2024-08_SNP_calling_from_assemblies/"$s"
    for v in *.vcf; do
        if [[ "$v" == "truth.vcf" ]]; then continue; fi
        vcfdist "$v" truth.vcf ../ref.fasta -p "$v".vcfdist/
    done
done
```

While vcfdist can potentially report different values for `TRUTH_TP` and `QUERY_TP` (see https://github.com/TimD1/vcfdist/issues/26), I checked my results to see if those two values are always the same:
```bash
cd ~/2024-08_SNP_calling_from_assemblies
for t in RES*/*.vcfdist/precision-recall-summary.tsv; do
    echo "$t"
    tail -n+2 "$t" | awk '{if ($4!=$5) print $0}'
done
```
They almost always were, except there were some differences for RES22-01840 Illumina 20x, which for some reason had a variant as `GCT`→`GTT` instead of just `C`→`T`. So I'll just use `TRUTH_TP` and ignore `QUERY_TP`.

I also checked my results for `SV` (structural variant) differences:
```bash
cd ~/2024-08_SNP_calling_from_assemblies
for t in RES*/*.vcfdist/precision-recall-summary.tsv; do
    echo "$t"
    grep -P "^SV\t" "$t" | awk '{if ($4+$5+$6+$7 > 0) print $0}'
done
```
And there were only a few cases (all 20x Illumina assemblies), so I will ignore that metric and just focus on SNP and INDEL stats.

Also, since I already did quality filtering on my VCFs (`QUAL>=100` for freebayes, `PASS` and not a het for Clair3), I'm using the `NONE` threshold results.

I then produced a results table:
```bash
cd ~/2024-08_SNP_calling_from_assemblies
printf "genome\tdepth\tread_method\tassembly_method\tvariant_call_method\tassembly_contigs\tassembly_size\tSNP_TP\tSNP_FN\tSNP_FP\tindel_TP\tindel_FN\tindel_FP\n" > results.tsv

# Read-based variant calling:
for r in illumina nanopore; do
    for t in RES*/"$r"_*.vcfdist/precision-recall-summary.tsv; do
        genome=$(echo "$t" | grep -oP 'RES22-\d+')
        depth=$(echo "$t" | grep -oP '_([0-9]{3})x' | grep -oP '[0-9]{3}' | sed 's/^0//')
        printf "$genome\t$depth\t$r\t\t\t\t\t" >> results.tsv
        grep -P "SNP\tNONE" "$t" | cut -f4,6,7 | tr '\n' '\t' >> results.tsv
        grep -P "INDEL\tNONE" "$t" | cut -f4,6,7 >> results.tsv
    done
done

# Reference assemblies (using full read sets):
for t in RES*/reference_*.vcfdist/precision-recall-summary.tsv; do
    genome=$(echo "$t" | grep -oP 'RES22-\d+')
    method=$(echo "$t" | grep -oP '_mummer|_shred|_ska' | tr -d '_')
    printf "$genome\t\t\treference\t$method\t" >> results.tsv
    seqtk size "$genome"/reference.fasta | tr '\n' '\t' >> results.tsv
    grep -P "SNP\tNONE" "$t" | cut -f4,6,7 | tr '\n' '\t' >> results.tsv
    grep -P "INDEL\tNONE" "$t" | cut -f4,6,7 >> results.tsv
done

# Other assemblies (using subsampled read sets):
for t in RES*/*_*_*.vcfdist/precision-recall-summary.tsv; do
    genome=$(echo "$t" | grep -oP 'RES22-\d+')
    method=$(echo "$t" | grep -oP '_mummer|_shred|_ska' | tr -d '_')
    depth=$(echo "$t" | grep -oP '_([0-9]{3})x' | grep -oP '[0-9]{3}' | sed 's/^0//')
    assembly=$(echo "$t" | grep -oP '/.+_\d+x' | tr -d '/')
    assembler=$(echo "$t" | grep -oP '/.+_\d+x'); assembler="${assembler:1:-5}"
    printf "$genome\t$depth\t\t$assembler\t$method\t" >> results.tsv
    seqtk size "$genome"/"$assembly".fasta | tr '\n' '\t' >> results.tsv
    grep -P "SNP\tNONE" "$t" | cut -f4,6,7 | tr '\n' '\t' >> results.tsv
    grep -P "INDEL\tNONE" "$t" | cut -f4,6,7 >> results.tsv
done
```



# Investigating Medaka errors

Check to see which long-read assemblies are missing plasmid 1. This code prints 'yes' if the plasmid is missing, 'no' if it is not:
```bash
cd ~/2024-08_SNP_calling_from_assemblies
grep -A1 "plasmid_1" ref.fasta > plasmid_1.fasta

for a in canu flye raven; do
    for s in RES22-01837 RES22-01839 RES22-01840 RES22-01841 RES22-01842 RES22-01844 RES22-01845; do
        for d in 020 050 100; do
            printf "$a $s $d\t"
            minimap2 -c -x asm5 plasmid_1.fasta "$s"/"$a"_"$d"x.fasta 2> /dev/null | awk '{if ($10 > 1000 && $10/$11 > 0.99) print $0;}' | { read && echo "no" || echo "yes"; }
        done
    done
done

rm plasmid_1.fasta
```

Check to see how many bases were missed, covered once and covered multiple times:
```bash
cd ~/2024-08_SNP_calling_from_assemblies
for a in canu flye raven; do
    for s in RES22-01837 RES22-01839 RES22-01840 RES22-01841 RES22-01842 RES22-01844 RES22-01845; do
        for d in 020 050 100; do
            printf "$a $s $d\t"
            minimap2 -c -x asm5 ref.fasta "$s"/"$a"_"$d"x.fasta 2> /dev/null | scripts/missing_once_multiple.py
        done
    done
done
```

These commands let me investigate a particular assembly, showing before-Medaka errors, after-Medaka errors, assembly-to-ref alignments and assembly contig headers:
```bash
a=canu; d=100; s=RES22-01839

cd ~/2024-08_SNP_calling_from_assemblies/"$s"/"$a"_"$d"x_shred.vcf.vcfdist; printf "\nBEFORE MEDAKA:\n"; comm -23 <(sort query.tsv | cut -f1-5) <(sort truth.tsv | cut -f1-5) | sort -k1,1 -k2,2n;  cd ~/2024-08_SNP_calling_from_assemblies/"$s"/"$a"_medaka_"$d"x_shred.vcf.vcfdist; printf "\nAFTER MEDAKA:\n"; comm -23 <(sort query.tsv | cut -f1-5) <(sort truth.tsv | cut -f1-5) | sort -k1,1 -k2,2n
printf "\n"; cd ~/2024-08_SNP_calling_from_assemblies; minimap2 -c -x asm5 ref.fasta "$s"/"$a"_"$d"x.fasta 2> /dev/null
printf "\n"; grep ">" "$s"/"$a"_"$d"x.fasta
```



# Tarball data for public repo

For the public data, I used the strain names (e.g. IMAL014) instead of the sample names (e.g. RES22-01837).

```bash
cd ~/2024-08_SNP_calling_from_assemblies
mkdir tarballs
cd tarballs
mkdir assemblies vcfs
for s in NRS384_wildtype NRS384_walKT389A IMAL014 IMAL031 IMAL058 IMAL065 IMAL070; do
    mkdir assemblies/"$s"
    mkdir vcfs/"$s"
done

cp ../RES22-01844/*.fasta assemblies/NRS384_wildtype
cp ../RES22-01845/*.fasta assemblies/NRS384_walKT389A
cp ../RES22-01837/*.fasta assemblies/IMAL014
cp ../RES22-01839/*.fasta assemblies/IMAL031
cp ../RES22-01840/*.fasta assemblies/IMAL058
cp ../RES22-01841/*.fasta assemblies/IMAL065
cp ../RES22-01842/*.fasta assemblies/IMAL070

cp ../RES22-01844/*.vcf vcfs/NRS384_wildtype
cp ../RES22-01845/*.vcf vcfs/NRS384_walKT389A
cp ../RES22-01837/*.vcf vcfs/IMAL014
cp ../RES22-01839/*.vcf vcfs/IMAL031
cp ../RES22-01840/*.vcf vcfs/IMAL058
cp ../RES22-01841/*.vcf vcfs/IMAL065
cp ../RES22-01842/*.vcf vcfs/IMAL070

tar -Jcvf vcfs.tar.xz --owner=0 --group=0 vcfs
tar -Jcvf assemblies.tar.xz --owner=0 --group=0 assemblies
```
