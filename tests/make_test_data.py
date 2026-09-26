#!/usr/bin/env python3
"""
Builds a tiny synthetic test set for bottleneck_prep (no real data needed):

  ref.fa          30 kb random "NC_000962.3"
  <sample>.bam    ~100x paired reads with iSNVs at KNOWN frequencies (bwa + samtools)
  <sample>.vcf.gz hand-written VCF in mtb-call2's LoFreq format (what LoFreq would call)
  mask.bed.gz     one masked interval (stands in for PE/PPE)
  samplesheet.tsv, pairs.tsv

Truth (pos: donor AF -> recipient AF):
  P1  D1 -> R1   5000: 0.30 -> 1.00   (transmitted, fixed)
                10000: 0.10 -> 0.00   (lost; recipient has no LoFreq call)
                15000: 0.50 -> 0.40
                20000: 0.25 -> 0.00   (inside mask)
                25000: 0.00 -> 0.20   (recipient-only / de novo)
  P2  D2 -> R2   8000: 0.20 -> 0.20
Requires: bwa, samtools, bgzip on PATH.
"""
import os, random, subprocess, gzip

random.seed(1)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(out, exist_ok=True)
CHROM, L = "NC_000962.3", 30000
ref = "".join(random.choice("ACGT") for _ in range(L))
alt_of = {"A": "G", "C": "T", "G": "A", "T": "C"}

truth = {
    "D1": {5000: 0.30, 10000: 0.10, 15000: 0.50, 20000: 0.25},
    "R1": {5000: 1.00, 15000: 0.40, 25000: 0.20},
    "D2": {8000: 0.20},
    "R2": {8000: 0.20},
}

with open(f"{out}/ref.fa", "w") as f:
    f.write(f">{CHROM}\n" + "\n".join(ref[i:i + 60] for i in range(0, L, 60)) + "\n")
subprocess.run(["samtools", "faidx", f"{out}/ref.fa"], check=True)
subprocess.run(["bwa", "index", f"{out}/ref.fa"], check=True, capture_output=True)

def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))

READ, FRAG, DEPTH = 150, 350, 100
for sample, sites in truth.items():
    n_frag = DEPTH * L // (2 * READ)
    with open(f"{out}/{sample}_1.fq", "w") as r1, open(f"{out}/{sample}_2.fq", "w") as r2:
        for i in range(n_frag):
            s = random.randint(0, L - FRAG)
            frag = list(ref[s:s + FRAG])
            for pos, af in sites.items():          # pos is 1-based
                if s < pos <= s + FRAG and random.random() < af:
                    frag[pos - 1 - s] = alt_of[ref[pos - 1]]
            frag = "".join(frag)
            q = "I" * READ
            r1.write(f"@{sample}_{i}/1\n{frag[:READ]}\n+\n{q}\n")
            r2.write(f"@{sample}_{i}/2\n{revcomp(frag[-READ:])}\n+\n{q}\n")
    bam = f"{out}/{sample}.bam"
    bwa = subprocess.run(["bwa", "mem", "-R", f"@RG\\tID:{sample}\\tSM:{sample}", f"{out}/ref.fa",
                          f"{out}/{sample}_1.fq", f"{out}/{sample}_2.fq"],
                         check=True, capture_output=True)
    subprocess.run(["samtools", "sort", "-o", bam, "-"], input=bwa.stdout, check=True)
    subprocess.run(["samtools", "index", bam], check=True)
    os.remove(f"{out}/{sample}_1.fq"); os.remove(f"{out}/{sample}_2.fq")

    # VCF in the shape mtb-call2's VariantsLoFreq writes (INFO DP/AF/SB/DP4 + GT column)
    vcf = f"{out}/{sample}_lofreq_unfilt.vcf"
    with open(vcf, "w") as f:
        f.write("##fileformat=VCFv4.0\n"
                f"##contig=<ID={CHROM},length={L}>\n"
                '##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">\n'
                '##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">\n'
                '##INFO=<ID=SB,Number=1,Type=Integer,Description="Phred-scaled strand bias">\n'
                '##INFO=<ID=DP4,Number=4,Type=Integer,Description="DP4">\n'
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
        for pos in sorted(sites):
            r = ref[pos - 1]
            f.write(f"{CHROM}\t{pos}\t.\t{r}\t{alt_of[r]}\t500\tPASS\t"
                    f"DP={DEPTH};AF={sites[pos]:.4f};SB=0;DP4=0,0,0,0\tGT\t1\n")
    subprocess.run(["bgzip", "-f", vcf], check=True)

with gzip.open(f"{out}/mask.bed.gz", "wt") as f:
    f.write(f"{CHROM}\t19900\t20100\tfake_PPE\n")

with open(f"{out}/samplesheet.tsv", "w") as f:
    f.write("sample_id\tbam\tlofreq_vcf\n")
    for s in truth:
        f.write(f"{s}\t{out}/{s}.bam\t{out}/{s}_lofreq_unfilt.vcf.gz\n")

with open(f"{out}/pairs.tsv", "w") as f:
    f.write("pair_id\tcohort\tdonor_sample_id\trecipient_sample_id\n"
            "P1\ttest\tD1\tR1\nP2\ttest\tD2\tR2\n")
print("test data written to", out)
