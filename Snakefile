#+ path for files 
SRRIDs = [x.strip() for x in open('data/samples').readlines()]
INDEX  = "ref/GCF_000633955.1_Cs_genomic"                      # reference fasta file need to be indexed with hisat2-index
FA     = "ref/GCF_000633955.1_Cs_genomic.fna"                  # reference fasta file
PFAM   = "/ref/analysis/References/Pfam-A.hmm"                 # pfam hmm for transdecoder
GFF    = "ref/GCF_000633955.1_Cs_genomic.gff"                  # original gff file
UNIPROT = "/ref/analysis/References/uniprot/uniprot-all.fasta" # uniprot for annotation
#-

#+ path for softwares 
HISAT2       = "/programs/hisat2-2.0.4/"
STRINGTIE    = "/programs/stringtie-1.2.4.Linux_x86_64//"
Transdecoder = "/programs/TransDecoder-3.0.0/"
AUGUSTUS     = "/programs/augustus-3.2.2/"
SPECIES      = "arabidopsis"
CUFFLINK     = "/programs/cufflinks-2.2.1.Linux_x86_64/"
SRATK        = "/programs/sratoolkit.2.7.0-ubuntu64/bin/"
#-

rule end:
     input : "finalout/my_csv.csv.addgene.gff3.sort.gff3.merge.all.gff3.sort.gff3","finalout/newgene.annot"

rule NCBIdownload:
     params : "{SRRID}",SRRID=SRRIDs
     output : "{SRRID}.sra" 
     shell  : "python2.7 ./py/ncbi_download.py {params}"

rule fastqdump:
     input  : "{SRRID}.sra"
     output : "data/{SRRID}_1.fastq.gz","data/{SRRID}_2.fastq.gz" # paired_end
     params : cmd=SRATK
     shell  : '''{params.cmd}fastq-dump  --origfmt -I  --split-files --gzip {input}
                 mv {wildcards.SRRID}_?.fastq.gz data/
                 rm {input}'''

           
# Single end 
rule Hisat2:
     input  : 
             #single="data/{SRRID}.fastq.gz", #single end
             fwd="data/{SRRID}_1.fastq.gz",rev="data/{SRRID}_2.fastq.gz" # paired end 
     params : ix=INDEX,cmd=HISAT2
     output : 
             "mapped/{SRRID}.pre.bam"
     threads : 2
     shell  : 
             #"{params.cmd}hisat2 --max-intronlen 30000 -p {threads} -x {params.ix} -U {input.single}  | sambamba view -f bam -o {output} -S /dev/stdin" # single end 
             "{params.cmd}hisat2 --max-intronlen 30000 -p {threads} -x {params.ix} -1 {input.fwd} -2 {input.rev}  | sambamba view -f bam -o {output} -S /dev/stdin" # paired-end
# add "--rna-strandness FR #RF for dUTP protocol" if your library is constructed based on strandness 

rule sambamba_sort : 
     input  : 
            "mapped/{SRRID}.pre.bam"
     output :
            "mapped/{SRRID}.sorted.bam"
            #"mapped/{SRRID}.sorted.bam.bai"
     threads : 10
     shell  : 
            '''sambamba sort -t {threads} -o {output}  {input}
               '''

#Usage: sambamba-merge [options] <output.bam> <input1.bam> <input2.bam> [...]
rule sambamba_merge : 
     input  : 
            bam=expand("mapped/{SRRID}.sorted.bam",SRRID=SRRIDs),
            #bai=expand("mapped/{SRRID}.sorted.bam.bai",SRRID=SRRIDs)
     output : 
            "merged/all.merged.bam"
     threads : 10
     shell  : "sambamba merge -t {threads} {output} {input.bam}"
            
rule stringtie:
     input  :
            "merged/all.merged.bam"
     output : 
            "st_gff_out/all.merged.bam.stringtie.gff"
     threads : 10 
     params : cmd=STRINGTIE
     shell  : 
            "{params.cmd}stringtie -p {threads} -o {output} {input}"

rule stringtie2cuffcompare:
     input : "st_gff_out/all.merged.bam.stringtie.gff"
     output : "st_gff_out/cuffcmp.all.merged.bam.stringtie.gff.tmap"
     params : cmd=CUFFLINK, ref=GFF
     shell  : "(params.cmd)cuffcompare -r {params.ref} {input}"

rule stringtie2gff3:
     input  :
            "st_gff_out/all.merged.bam.stringtie.gff"
     output : 
            "st_gff_out/all.merged.bam.stringtie.gff3"
     params : cmd=Transdecoder
     shell  : 
            "{params.cmd}/util/cufflinks_gtf_to_alignment_gff3.pl {input} > {output}"

rule stringtie2cds:
     input  : "st_gff_out/all.merged.bam.stringtie.gff"
     output : "st_gff_out/all.merged.bam.stringtie.gff.fa"
     params : fa=FA,cmd=Transdecoder
     shell  : 
            "{params.cmd}/util/cufflinks_gtf_genome_to_cdna_fasta.pl {input}  {params.fa} > {output}"


### transdecoder 

rule transdecoder_longOrfs:
     input  : "st_gff_out/all.merged.bam.stringtie.gff.fa"
     output : "all.merged.bam.stringtie.gff.fa.transdecoder_dir/longest_orfs.pep" 
     params : cmd=Transdecoder
     shell  : 
            "{params.cmd}/TransDecoder.LongOrfs -t {input}"

rule hmm:
     input  : "all.merged.bam.stringtie.gff.fa.transdecoder_dir/longest_orfs.pep" 
     output : "all.merged.bam.stringtie.gff.fa.transdecoder_dir/pfam.domtblout"
     threads: 8
     params : pfam=PFAM
     shell  : 
            "hmmscan --cpu {threads} --domtblout {output} {params.pfam} {input}"

rule transdecoder_predict:
     input  : fa="st_gff_out/all.merged.bam.stringtie.gff.fa", pfam="all.merged.bam.stringtie.gff.fa.transdecoder_dir/pfam.domtblout"
     output : "predicted/all.merged.bam.stringtie.gff.fa.transdecoder.gff3" 
     threads: 10
     params : cmd=Transdecoder
     shell  : 
            '''{params.cmd}/TransDecoder.Predict --cpu {threads} -t {input.fa} --retain_pfam_hits {input.pfam}        
               mv all.merged.bam.stringtie.gff.fa.transdecoder.* predicted/'''
        
rule transdecoder_togenome:
     input  : sgff="st_gff_out/all.merged.bam.stringtie.gff", tgff="predicted/all.merged.bam.stringtie.gff.fa.transdecoder.gff3"
     output : "predicted/all.merged.bam.stringtie.gff.fa.transdecoder.gff3.genome.gff"
     shell  : "python2.7 ./py/stringtie.addcds.py {input.sgff} {input.tgff}"

rule transdecoder_noiseremove:
     input  : "predicted/all.merged.bam.stringtie.gff.fa.transdecoder.gff3.genome.gff"
     output : "predicted/selected_mRNA_v4.gff"
     shell  : '''python2.7 ./py/noiseremove.py {input} 
                 mv selected_mRNA_v4.gff predicted/'''
###

### Augustus 

rule augustus:
     input  : "st_gff_out/all.merged.bam.stringtie.gff.fa"
     output : "predicted/all.merged.bam.stringtie.gff.fa.augustus.gff3"
     params : cmd=AUGUSTUS,species=SPECIES
     shell  : "{params.cmd}/bin/augustus --species={params.species} --genemodel=complete --gff3=on --strand=forward {input} > {output}"

rule augustus_togenome:
     input  : sgff="st_gff_out/all.merged.bam.stringtie.gff",tgff="predicted/all.merged.bam.stringtie.gff.fa.augustus.gff3"
     output : "predicted/all.merged.bam.stringtie.gff.fa.augustus.gff3.genome.v1.gff"
     shell  :  '''python2.7 ./py/stringtie.augustus.addcds.py {input.sgff} {input.tgff} {output}'''

###

rule integrate:
     input  : fa="st_gff_out/all.merged.bam.stringtie.gff.fa", ag="predicted/all.merged.bam.stringtie.gff.fa.augustus.gff3.genome.v1.gff", td="predicted/selected_mRNA_v4.gff"
     output : "predicted/my_csv.csv"
     shell  : '''python2.7 ./py/integrate.py {input.fa} {input.ag} {input.td}
                 mv my_csv.csv predicted/'''


rule addgenes:
     input  : "predicted/my_csv.csv"
     output : "predicted/my_csv.csv.addgene.gff3"
     shell  : "python2.7 ./py/addgenefeature.py {input}" 
                
rule gt_sort:
     input  : "predicted/my_csv.csv.addgene.gff3"
     output : "predicted/my_csv.csv.addgene.gff3.sort.gff3"
     shell  : "gt gff3 -sort -retainids {input} > {output}" 


rule cuffcompare:
     input  : "predicted/my_csv.csv.addgene.gff3.sort.gff3"
     output : "predicted/cuffcmp.my_csv.csv.addgene.gff3.sort.gff3.tmap"
     params : rgff=GFF, cmd=CUFFLINK
     shell  : '''{params.cmd}cuffcompare -r {params.rgff} {input}
                 mv cuffcmp.* predicted/'''

rule merge:
     input  : gff="predicted/my_csv.csv.addgene.gff3.sort.gff3", tmap="predicted/cuffcmp.my_csv.csv.addgene.gff3.sort.gff3.tmap", gff_old=GFF, pepfa_new="finalout/my_csv.csv.addgene.gff3.sort.gff3.pep.fa"
     output : gff="predicted/my_csv.csv.addgene.gff3.sort.gff3.merge.all.gff3", fa="finalout/my_csv.csv.addgene.gff3.sort.gff3.pep.fa.new_gene.fa"          
     params : rgff=GFF 
     shell  : '''python2.7 ./py/merge.py {input.tmap} {input.gff} {input.gff_old} {input.pepfa_new}
               cat {params.rgff} predicted/my_csv.csv.addgene.gff3.sort.gff3.merge.gff3 > {output.gff}'''

rule gt_sort2:
     input  : "predicted/my_csv.csv.addgene.gff3.sort.gff3.merge.all.gff3"
     output : "finalout/my_csv.csv.addgene.gff3.sort.gff3.merge.all.gff3.sort.gff3"
     shell  : "gt gff3 -sort -retainids {input} > {output}"


rule gff2cds:
     input  : gff="predicted/my_csv.csv.addgene.gff3.sort.gff3"
     params : fa=FA
     output : "finalout/my_csv.csv.addgene.gff3.sort.gff3.cds.fa"
     shell  : "python2.7 ./py/gff2cds.py {params.fa} {input.gff} {output}"

rule cds2pep:
     input  : "finalout/my_csv.csv.addgene.gff3.sort.gff3.cds.fa"
     output : "finalout/my_csv.csv.addgene.gff3.sort.gff3.pep.fa"
     shell  : "python2.7 ./py/translation.py {input} {output}"

rule newgenes2uniprot:
     input  : "finalout/my_csv.csv.addgene.gff3.sort.gff3.pep.fa.new_gene.fa"
     output : "blastout/newgene.na1.ev1e2.bo"
     params : uniprot=UNIPROT
     threads: 5 
     shell  : "blastp -query {input} -db {params.uniprot} -out {output} -outfmt 7 -num_threads {threads} -num_alignments 1 -evalue 1e-2"  

rule uniprot2annot:
     input  : "blastout/newgene.na1.ev1e2.bo"
     output : "finalout/newgene.annot"
     params : uniprot=UNIPROT
     shell  : "python2.7 ./py/uniprot_annot.py {input} {params.uniprot}; mv blastout/newgene.na1.ev1e2.bo.annot {output} "
     
