# Vignette Data Generation

The data used in the vignette is available in this directory. The VCFs
are taken from the most recent 1000G whole genome sequencing release
(<https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/>)

The code used to generate the VCFs, annotation files, and summary
statistics for the *NEB* gene are available in the file
1000G\_commands.txt.

We used the following process to generate the summary statistic data:

1.  We first removed all variants with 0 &lt; MAF &lt; 0.01 using allele
    frequencies from the European population of 1000G only.

2.  Using the rare-variant annotation file from 1000G, we assigned to
    each PTV a probability of being causal of 0.80; to each predicted
    deleterious missense variant a probability of being causal of 0.30;
    and to all other missense variants a probability of being causal of
    0.01. We defined a missese variant to be predicted deleterious if
    each of MutationTaster, Polyphen2\_HDIV, Polyphen2\_HVAR, and SIFT
    predicted the variant to be deleterious.

3.  Using *runif()* we then randomly assigned variants to be causal
    based on the probabilities in the previous step.

4.  For each causal variant, we set its effect size to be proportional
    to the -log10 of its MAF.

5.  We set h<sup>2</sup> equal to 2 times the sum of the true causal
    effect squared times MAF(1-MAF) for each causal variant.

6.  To ensure some truly associated variants, we set nh<sup>2</sup>
    = 50. Since n = 503, this meant assigning h<sup>2</sup> = 0.0994. We
    then calculate the proportionality constant for the true effect
    sizes.

7.  We simulated the phenotype from a normal distribution with mean
    equal to the linear combination of the true effect size and
    individual level genotype, and the variance equal to
    1-h<sup>2</sup>.
