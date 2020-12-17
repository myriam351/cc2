CC2 - Écogénomique - Myriam FERBLANTIER - N°22000007
================

  - [Commentaire : Recharge des packages dada2 et
    Rcpp.](#commentaire-recharge-des-packages-dada2-et-rcpp.)
  - [Commentaire : Retraçage du fichier de
    données.](#commentaire-retraçage-du-fichier-de-données.)
  - [Commentaire: Nous avons obtenus les listes correspondantes des
    fichiers fastaq pour les amorces Foward et
    Reverse.](#commentaire-nous-avons-obtenus-les-listes-correspondantes-des-fichiers-fastaq-pour-les-amorces-foward-et-reverse.)
  - [Inspect read quality profiles](#inspect-read-quality-profiles)
  - [Filter and trim](#filter-and-trim)
  - [Learn the Error Rates](#learn-the-error-rates)
  - [Sample Inference](#sample-inference)
  - [Merge paired reads](#merge-paired-reads)
  - [Construct sequence table](#construct-sequence-table)
  - [Revome chimeras](#revome-chimeras)
  - [Track reads through the
    pipeline](#track-reads-through-the-pipeline)
  - [Assign taxonomy](#assign-taxonomy)

``` r
library(Rcpp)
library(dada2)
```

### Commentaire : Recharge des packages dada2 et Rcpp.

``` r
path <- "~/CC2/donnees" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
```

    ##  [1] "filtered"                            "Station5_Fond1_10sept14_R1.fastq"   
    ##  [3] "Station5_Fond1_10sept14_R2.fastq"    "Station5_Fond1_11mars15_R1.fastq"   
    ##  [5] "Station5_Fond1_11mars15_R2.fastq"    "Station5_Fond2_10sept14_R1.fastq"   
    ##  [7] "Station5_Fond2_10sept14_R2.fastq"    "Station5_Fond2_11mars15_R1.fastq"   
    ##  [9] "Station5_Fond2_11mars15_R2.fastq"    "Station5_Fond3_10sept14_R1.fastq"   
    ## [11] "Station5_Fond3_10sept14_R2.fastq"    "Station5_Median1_10sept14_R1.fastq" 
    ## [13] "Station5_Median1_10sept14_R2.fastq"  "Station5_Median2_10sept14_R1.fastq" 
    ## [15] "Station5_Median2_10sept14_R2.fastq"  "Station5_Surface1_10sept14_R1.fastq"
    ## [17] "Station5_Surface1_10sept14_R2.fastq" "Station5_Surface1_11mars15_R1.fastq"
    ## [19] "Station5_Surface1_11mars15_R2.fastq" "Station5_Surface2_10sept14_R1.fastq"
    ## [21] "Station5_Surface2_10sept14_R2.fastq" "Station5_Surface2_11mars15_R1.fastq"
    ## [23] "Station5_Surface2_11mars15_R2.fastq"

### Commentaire : Retraçage du fichier de données.

``` r
# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnRs), "R"), '[', 1)
```

### Commentaire: Nous avons obtenus les listes correspondantes des fichiers fastaq pour les amorces Foward et Reverse.

## Inspect read quality profiles

``` r
plotQualityProfile(fnFs[1:2])
```

![](02_stat-analysiscc2_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
plotQualityProfile(fnRs[1:2])
```

![](02_stat-analysiscc2_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->
\#\#\#Commentaire: Inspection par contrôle qualité des amorces Foward et
Reverse. On peut observer que les amorces Reverse présente un score de
qualité qui diminue à un intervalle entre 200 et 250. En comparaison, on
peut voir que pour les amorces Forward on a un bon score de qualité,
avec une légère diminution vers le 250 ème cycle.

## Filter and trim

``` r
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

### Commentaire: Attributions des noms de fichiers pour les fichiers fastq.gz filtrés.

``` r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,200), trimLeft = c(21,21),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
```

    ##                                    reads.in reads.out
    ## Station5_Fond1_10sept14_R1.fastq     159971    145031
    ## Station5_Fond1_11mars15_R1.fastq     175993    159807
    ## Station5_Fond2_10sept14_R1.fastq     197039    176477
    ## Station5_Fond2_11mars15_R1.fastq      87585     79696
    ## Station5_Fond3_10sept14_R1.fastq     117140    105805
    ## Station5_Median1_10sept14_R1.fastq   116519    106244

### Commentaire : Nous avons utilisé des paramètres de filtrage pour filtrer nos données (soit maxN=0, truncQ=2,trimLeft = c(21,21), rm.phix=TRUE, maxEE=2). Par exemple: Le paramètre maxEE définit le nombre maximum d’“erreurs attendues” autorisées dans une lecture.

### Commentaire : On a sélectionné les amorces tronquées et on les a filtrés en précisent les zones (soit 250 à 200).

## Learn the Error Rates

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 110221135 total bases in 481315 reads from 3 samples will be used for learning the error rates.

``` r
errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 100420969 total bases in 561011 reads from 4 samples will be used for learning the error rates.

### Commentaire : Ici on apprend via la méthode learnErrors les taux d’erreurs à partir des amorces. Pour les amorces Foward on a utilisé 3 échantillons et pour les Reverse on a utilisé 4 échantillons.

``` r
plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    
    ## Warning: Transformation introduced infinite values in continuous y-axis

![](02_stat-analysiscc2_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->
\#\#\# Commentaire : On peut Visualiser des taux d’erreurs estimés, via
la fonction plotErrors().

## Sample Inference

``` r
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 145031 reads in 39978 unique sequences.
    ## Sample 2 - 159807 reads in 37782 unique sequences.
    ## Sample 3 - 176477 reads in 49703 unique sequences.
    ## Sample 4 - 79696 reads in 21434 unique sequences.
    ## Sample 5 - 105805 reads in 31881 unique sequences.
    ## Sample 6 - 106244 reads in 30070 unique sequences.
    ## Sample 7 - 98411 reads in 26954 unique sequences.
    ## Sample 8 - 106995 reads in 28021 unique sequences.
    ## Sample 9 - 70842 reads in 18914 unique sequences.
    ## Sample 10 - 78294 reads in 21347 unique sequences.
    ## Sample 11 - 91238 reads in 25826 unique sequences.

### Commentaire : Visualisation des taux d’erreurs pour les Forward.

``` r
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 145031 reads in 45252 unique sequences.
    ## Sample 2 - 159807 reads in 41345 unique sequences.
    ## Sample 3 - 176477 reads in 55267 unique sequences.
    ## Sample 4 - 79696 reads in 23050 unique sequences.
    ## Sample 5 - 105805 reads in 34435 unique sequences.
    ## Sample 6 - 106244 reads in 31383 unique sequences.
    ## Sample 7 - 98411 reads in 28878 unique sequences.
    ## Sample 8 - 106995 reads in 28735 unique sequences.
    ## Sample 9 - 70842 reads in 21298 unique sequences.
    ## Sample 10 - 78294 reads in 21877 unique sequences.
    ## Sample 11 - 91238 reads in 28105 unique sequences.

### Commentaire : Visualisation des taux d’erreurs pour les Reverse.

``` r
dadaFs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 1016 sequence variants were inferred from 39978 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

### Commentaire : Ici on a effectué une inspection des données, avec 1016 variantes de séquences à partir des séquences uniques de 39978 dans le premier échantillon.

``` r
dadaFs[[3]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 1170 sequence variants were inferred from 49703 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

### Commentaire : On a fait la même chose que précédemment, mais avec l’échantillon 3. On a effectué une inspection des données, avec 1170 variantes de séquences à partir des séquences uniques de 49703 dans le premier échantillon.

## Merge paired reads

``` r
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

    ## 115473 paired-reads (in 4595 unique pairings) successfully merged out of 140326 (in 21160 pairings) input.

    ## 136539 paired-reads (in 3685 unique pairings) successfully merged out of 155618 (in 15520 pairings) input.

    ## 139935 paired-reads (in 6236 unique pairings) successfully merged out of 170590 (in 26605 pairings) input.

    ## 66207 paired-reads (in 2312 unique pairings) successfully merged out of 77396 (in 9469 pairings) input.

    ## 82090 paired-reads (in 3117 unique pairings) successfully merged out of 101553 (in 16083 pairings) input.

    ## 85833 paired-reads (in 3244 unique pairings) successfully merged out of 102767 (in 13959 pairings) input.

    ## 80200 paired-reads (in 2622 unique pairings) successfully merged out of 95148 (in 12056 pairings) input.

    ## 89039 paired-reads (in 3012 unique pairings) successfully merged out of 103729 (in 11950 pairings) input.

    ## 58701 paired-reads (in 1682 unique pairings) successfully merged out of 68495 (in 7924 pairings) input.

    ## 65924 paired-reads (in 1731 unique pairings) successfully merged out of 76263 (in 8159 pairings) input.

    ## 73150 paired-reads (in 2621 unique pairings) successfully merged out of 87831 (in 11884 pairings) input.

### Commentaire : Fusion des lectures d’amorces Forward et Reverse. Un alignement est effectué au préalable. Par exemple, pour le premier, 115473 paires de lectures (dans 4595 paires uniques) ont été fusionnées avec succès à partir de 140326 (dans 21160 paires) entrées.

``` r
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
```

    ##                                                                                                                                                                                                                                                                                                                                                                                sequence
    ## 1     TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTAGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATTAAAACTTTTCAGCTAGAGTATGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCATTACTGACACTGAGGAACGAAAGCATGGGTAGCGAAGAGGATTAGATACCCTCGTAGTCCATGCCGTAAACGATGTGTGTTAGACGTTGGAAATTTATTTTCAGTGTCGCAGGGAAACCGATAAACACACCGCCTGGGGAGTACGACCGCAAGGTT
    ## 2     TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTGGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATCAAAACTTTTCAGCTAGAGTATGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCATTACTGACACTGAGGAACGAAAGCATGGGTAGCGAAGAGGATTAGATACCCTCGTAGTCCATGCCGTAAACGATGTGTGTTAGACGTTGGAAATTTATTTTCAGTGTCGCAGCGAAAGCGATAAACACACCGCCTGGGGAGTACGACCGCAAGGTT
    ## 3     TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTGGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATCAAAACTTTTCAGCTAGAGTTTGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCATTACTGACACTGAGGAACGAAAGCATGGGTAGCGAAGAGGATTAGATACCCTCGTAGTCCATGCCGTAAACGATGTGTGTTAGACGTTGGAAATTTATTTTCAGTGTCGCAGCGAAAGCGATAAACACACCGCCTGGGGAGTACGACCGCAAGGTT
    ## 4     TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTAGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATTAAAACTTTTCAGCTAGAGTATGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCATTACTGACACTGAGGAACGAAAGCATGGGTAGCGAAGAGGATTAGATACCCTCGTAGTCCATGCCGTAAACGATGTGTGTTAGACGTTGGAAATTTATTTTCAGTGTCGCAGCGAAAGCGATAAACACACCGCCTGGGGAGTACGACCGCAAGGTT
    ## 5     TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTGGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATCAAAACTTTTCAGCTAGAGTATGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCATTACTGACACTGAGGAACGAAAGCATGGGTAGCGAAGAGGATTAGATACCCTCGTAGTCCATGCCGTAAACGATGTGTGTTAGACGTTGGAAATTTATTTTCAGTGTCGCAGGGAAACCGATAAACACACCGCCTGGGGAGTACGACCGCAAGGTT
    ## 6 TACGAGGGGTCCTAGCGTTGTCCGGATTTACTGGGCGTAAAGGGTACGTAGGCGTTTTAATAAGTTGTATGTTAAATATCTTAGCTTAACTAAGAAAGTGCATACAAAACTGTTAAGATAGAGTTTGAGAGAGGAACGCAGAATTCATGGTGGAGCGGTGACATGCGTAGATATCATGAGGAAAGTCAAATGCGAAGGCAGCCTTCTGGCTCAAAACTGACGCTGAGGTACGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTATTTGGTGCTGGGGGATTCGACCCTTTCAGTGCCGTAGCTAACGCGATAAATACTCCGCCTGGGGACTACGATCGCAAGATT
    ##   abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 1      5170       1       2     39         0      0      2   TRUE
    ## 2      4127       2       1     39         0      0      2   TRUE
    ## 3      3781       3       1     39         0      0      2   TRUE
    ## 4      2481       1       1     39         0      0      2   TRUE
    ## 5      2176       2       2     39         0      0      2   TRUE
    ## 6      2130       5       9     35         0      0      1   TRUE

## Construct sequence table

``` r
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

    ## [1]    11 17389

### Commentaire :Ici la fonction makeSequenceTable () permet de construire une table de séquences (analogue à une table OTU) à partir de la liste d’échantillons fournie.

``` r
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
```

    ## 
    ##  352  353  362  363  364  365  366  367  368  369  370  371  372  373  374  375 
    ##    1    1    1    1    4  187   23  161  170 4913 3125 2167 2285 2492  101 1651 
    ##  376  377  378  382  386  387  389  392 
    ##   91    4    1    1    2    1    1    5

### Commentaire : Construction d’une table de variantes de séquences d’amplicon de 1 117 389.

## Revome chimeras

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```

    ## Identified 15867 bimeras out of 17389 input sequences.

### Commentaire : On a identifié des séquences chimériques, en utilisant la fonction removeBimeraDenovo(). On a pu identifier 15 867 de chimère sur 17 389 de séquences.

``` r
dim(seqtab.nochim)
```

    ## [1]   11 1522

### Commentaire : la fonction dim() permet de définir la dimension de objet.

``` r
sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.7840742

``` r
1-sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.2159258

### Commentaire : Il y a 2.2 % de séquence chimérique dans notre séquence unique.

## Track reads through the pipeline

``` r
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

    ##                             input filtered denoisedF denoisedR merged nonchim
    ## Station5_Fond1_10sept14_   159971   145031    142239    142879 115473   87651
    ## Station5_Fond1_11mars15_   175993   159807    157275    157884 136539  110945
    ## Station5_Fond2_10sept14_   197039   176477    172714    174073 139935  103398
    ## Station5_Fond2_11mars15_    87585    79696     78289     78639  66207   54402
    ## Station5_Fond3_10sept14_   117140   105805    103125    104001  82090   63805
    ## Station5_Median1_10sept14_ 116519   106244    104130    104717  85833   65324

### Commentaire : Vérification de tout ce qu’on a fait depuis le début. Par exemple: on peut voir qu’on a conservé la majorité de nos lectures brutes.

## Assign taxonomy

``` bash
wget https://zenodo.org/record/3986799/files/silva_nr99_v138_train_set.fa.gz
```

    ## --2020-12-17 13:57:02--  https://zenodo.org/record/3986799/files/silva_nr99_v138_train_set.fa.gz
    ## Resolving zenodo.org (zenodo.org)... 137.138.76.77
    ## Connecting to zenodo.org (zenodo.org)|137.138.76.77|:443... connected.
    ## HTTP request sent, awaiting response... 200 OK
    ## Length: 137973851 (132M) [application/octet-stream]
    ## Saving to: ‘silva_nr99_v138_train_set.fa.gz.3’
    ## 
    ##      0K .......... .......... .......... .......... ..........  0% 6.74M 20s
    ##     50K .......... .......... .......... .......... ..........  0% 13.0M 15s
    ##    100K .......... .......... .......... .......... ..........  0% 4.51M 20s
    ##    150K .......... .......... .......... .......... ..........  0% 13.6M 17s
    ##    200K .......... .......... .......... .......... ..........  0% 8.03M 17s
    ##    250K .......... .......... .......... .......... ..........  0% 87.5M 14s
    ##    300K .......... .......... .......... .......... ..........  0% 98.4M 13s
    ##    350K .......... .......... .......... .......... ..........  0% 18.5M 12s
    ##    400K .......... .......... .......... .......... ..........  0% 37.2M 11s
    ##    450K .......... .......... .......... .......... ..........  0% 39.6M 10s
    ##    500K .......... .......... .......... .......... ..........  0% 73.5M 9s
    ##    550K .......... .......... .......... .......... ..........  0% 75.2M 9s
    ##    600K .......... .......... .......... .......... ..........  0% 84.2M 8s
    ##    650K .......... .......... .......... .......... ..........  0% 28.7M 8s
    ##    700K .......... .......... .......... .......... ..........  0% 93.8M 7s
    ##    750K .......... .......... .......... .......... ..........  0% 11.3M 8s
    ##    800K .......... .......... .......... .......... ..........  0%  102M 7s
    ##    850K .......... .......... .......... .......... ..........  0% 86.1M 7s
    ##    900K .......... .......... .......... .......... ..........  0%  113M 7s
    ##    950K .......... .......... .......... .......... ..........  0% 97.3M 6s
    ##   1000K .......... .......... .......... .......... ..........  0% 96.7M 6s
    ##   1050K .......... .......... .......... .......... ..........  0% 91.7M 6s
    ##   1100K .......... .......... .......... .......... ..........  0% 50.6M 6s
    ##   1150K .......... .......... .......... .......... ..........  0% 36.8M 6s
    ##   1200K .......... .......... .......... .......... ..........  0% 15.5M 6s
    ##   1250K .......... .......... .......... .......... ..........  0% 78.5M 6s
    ##   1300K .......... .......... .......... .......... ..........  1%  108M 6s
    ##   1350K .......... .......... .......... .......... ..........  1%  106M 5s
    ##   1400K .......... .......... .......... .......... ..........  1% 65.9M 5s
    ##   1450K .......... .......... .......... .......... ..........  1% 83.6M 5s
    ##   1500K .......... .......... .......... .......... ..........  1% 86.5M 5s
    ##   1550K .......... .......... .......... .......... ..........  1% 83.2M 5s
    ##   1600K .......... .......... .......... .......... ..........  1% 86.8M 5s
    ##   1650K .......... .......... .......... .......... ..........  1% 78.6M 5s
    ##   1700K .......... .......... .......... .......... ..........  1% 97.4M 5s
    ##   1750K .......... .......... .......... .......... ..........  1% 47.5M 5s
    ##   1800K .......... .......... .......... .......... ..........  1%  110M 4s
    ##   1850K .......... .......... .......... .......... ..........  1% 58.7M 4s
    ##   1900K .......... .......... .......... .......... ..........  1% 66.2M 4s
    ##   1950K .......... .......... .......... .......... ..........  1% 74.9M 4s
    ##   2000K .......... .......... .......... .......... ..........  1% 94.0M 4s
    ##   2050K .......... .......... .......... .......... ..........  1% 47.4M 4s
    ##   2100K .......... .......... .......... .......... ..........  1%  102M 4s
    ##   2150K .......... .......... .......... .......... ..........  1% 39.3M 4s
    ##   2200K .......... .......... .......... .......... ..........  1% 79.6M 4s
    ##   2250K .......... .......... .......... .......... ..........  1% 94.9M 4s
    ##   2300K .......... .......... .......... .......... ..........  1% 81.8M 4s
    ##   2350K .......... .......... .......... .......... ..........  1% 53.9M 4s
    ##   2400K .......... .......... .......... .......... ..........  1% 71.4M 4s
    ##   2450K .......... .......... .......... .......... ..........  1% 99.6M 4s
    ##   2500K .......... .......... .......... .......... ..........  1% 37.7M 4s
    ##   2550K .......... .......... .......... .......... ..........  1% 95.2M 4s
    ##   2600K .......... .......... .......... .......... ..........  1% 34.4M 4s
    ##   2650K .......... .......... .......... .......... ..........  2% 60.9M 4s
    ##   2700K .......... .......... .......... .......... ..........  2%  117M 4s
    ##   2750K .......... .......... .......... .......... ..........  2% 85.0M 4s
    ##   2800K .......... .......... .......... .......... ..........  2% 79.0M 4s
    ##   2850K .......... .......... .......... .......... ..........  2% 95.2M 4s
    ##   2900K .......... .......... .......... .......... ..........  2% 86.3M 3s
    ##   2950K .......... .......... .......... .......... ..........  2% 94.2M 3s
    ##   3000K .......... .......... .......... .......... ..........  2% 92.5M 3s
    ##   3050K .......... .......... .......... .......... ..........  2% 78.7M 3s
    ##   3100K .......... .......... .......... .......... ..........  2% 89.9M 3s
    ##   3150K .......... .......... .......... .......... ..........  2% 65.5M 3s
    ##   3200K .......... .......... .......... .......... ..........  2% 94.9M 3s
    ##   3250K .......... .......... .......... .......... ..........  2% 73.4M 3s
    ##   3300K .......... .......... .......... .......... ..........  2% 94.2M 3s
    ##   3350K .......... .......... .......... .......... ..........  2% 53.3M 3s
    ##   3400K .......... .......... .......... .......... ..........  2% 83.8M 3s
    ##   3450K .......... .......... .......... .......... ..........  2%  106M 3s
    ##   3500K .......... .......... .......... .......... ..........  2% 66.8M 3s
    ##   3550K .......... .......... .......... .......... ..........  2% 73.1M 3s
    ##   3600K .......... .......... .......... .......... ..........  2% 78.2M 3s
    ##   3650K .......... .......... .......... .......... ..........  2% 94.1M 3s
    ##   3700K .......... .......... .......... .......... ..........  2% 57.7M 3s
    ##   3750K .......... .......... .......... .......... ..........  2% 69.6M 3s
    ##   3800K .......... .......... .......... .......... ..........  2% 81.2M 3s
    ##   3850K .......... .......... .......... .......... ..........  2% 66.8M 3s
    ##   3900K .......... .......... .......... .......... ..........  2% 95.1M 3s
    ##   3950K .......... .......... .......... .......... ..........  2% 83.9M 3s
    ##   4000K .......... .......... .......... .......... ..........  3% 74.5M 3s
    ##   4050K .......... .......... .......... .......... ..........  3% 56.2M 3s
    ##   4100K .......... .......... .......... .......... ..........  3% 82.4M 3s
    ##   4150K .......... .......... .......... .......... ..........  3% 84.8M 3s
    ##   4200K .......... .......... .......... .......... ..........  3% 94.5M 3s
    ##   4250K .......... .......... .......... .......... ..........  3% 85.6M 3s
    ##   4300K .......... .......... .......... .......... ..........  3% 23.9M 3s
    ##   4350K .......... .......... .......... .......... ..........  3%  102M 3s
    ##   4400K .......... .......... .......... .......... ..........  3%  117M 3s
    ##   4450K .......... .......... .......... .......... ..........  3%  115M 3s
    ##   4500K .......... .......... .......... .......... ..........  3% 89.5M 3s
    ##   4550K .......... .......... .......... .......... ..........  3%  122M 3s
    ##   4600K .......... .......... .......... .......... ..........  3% 99.8M 3s
    ##   4650K .......... .......... .......... .......... ..........  3%  101M 3s
    ##   4700K .......... .......... .......... .......... ..........  3% 26.2M 3s
    ##   4750K .......... .......... .......... .......... ..........  3% 41.4M 3s
    ##   4800K .......... .......... .......... .......... ..........  3% 98.6M 3s
    ##   4850K .......... .......... .......... .......... ..........  3% 88.7M 3s
    ##   4900K .......... .......... .......... .......... ..........  3%  122M 3s
    ##   4950K .......... .......... .......... .......... ..........  3% 41.1M 3s
    ##   5000K .......... .......... .......... .......... ..........  3%  123M 3s
    ##   5050K .......... .......... .......... .......... ..........  3% 96.4M 3s
    ##   5100K .......... .......... .......... .......... ..........  3%  116M 3s
    ##   5150K .......... .......... .......... .......... ..........  3% 72.7M 3s
    ##   5200K .......... .......... .......... .......... ..........  3% 48.4M 3s
    ##   5250K .......... .......... .......... .......... ..........  3% 38.0M 3s
    ##   5300K .......... .......... .......... .......... ..........  3%  133M 3s
    ##   5350K .......... .......... .......... .......... ..........  4%  122M 3s
    ##   5400K .......... .......... .......... .......... ..........  4%  154M 3s
    ##   5450K .......... .......... .......... .......... ..........  4% 23.5M 3s
    ##   5500K .......... .......... .......... .......... ..........  4% 70.4M 3s
    ##   5550K .......... .......... .......... .......... ..........  4%  115M 3s
    ##   5600K .......... .......... .......... .......... ..........  4% 97.9M 3s
    ##   5650K .......... .......... .......... .......... ..........  4%  116M 3s
    ##   5700K .......... .......... .......... .......... ..........  4%  129M 3s
    ##   5750K .......... .......... .......... .......... ..........  4%  117M 3s
    ##   5800K .......... .......... .......... .......... ..........  4% 81.0M 3s
    ##   5850K .......... .......... .......... .......... ..........  4% 48.1M 3s
    ##   5900K .......... .......... .......... .......... ..........  4% 90.6M 3s
    ##   5950K .......... .......... .......... .......... ..........  4% 82.7M 3s
    ##   6000K .......... .......... .......... .......... ..........  4% 42.0M 3s
    ##   6050K .......... .......... .......... .......... ..........  4% 16.8M 3s
    ##   6100K .......... .......... .......... .......... ..........  4% 68.0M 3s
    ##   6150K .......... .......... .......... .......... ..........  4% 96.8M 3s
    ##   6200K .......... .......... .......... .......... ..........  4% 91.9M 3s
    ##   6250K .......... .......... .......... .......... ..........  4% 82.6M 3s
    ##   6300K .......... .......... .......... .......... ..........  4% 95.0M 3s
    ##   6350K .......... .......... .......... .......... ..........  4% 88.9M 3s
    ##   6400K .......... .......... .......... .......... ..........  4% 63.3M 3s
    ##   6450K .......... .......... .......... .......... ..........  4% 47.4M 3s
    ##   6500K .......... .......... .......... .......... ..........  4% 94.4M 3s
    ##   6550K .......... .......... .......... .......... ..........  4% 91.6M 3s
    ##   6600K .......... .......... .......... .......... ..........  4%  103M 2s
    ##   6650K .......... .......... .......... .......... ..........  4% 93.8M 2s
    ##   6700K .......... .......... .......... .......... ..........  5% 79.9M 2s
    ##   6750K .......... .......... .......... .......... ..........  5%  103M 2s
    ##   6800K .......... .......... .......... .......... ..........  5% 17.0M 3s
    ##   6850K .......... .......... .......... .......... ..........  5% 69.4M 2s
    ##   6900K .......... .......... .......... .......... ..........  5%  109M 2s
    ##   6950K .......... .......... .......... .......... ..........  5% 97.9M 2s
    ##   7000K .......... .......... .......... .......... ..........  5% 92.8M 2s
    ##   7050K .......... .......... .......... .......... ..........  5% 99.6M 2s
    ##   7100K .......... .......... .......... .......... ..........  5%  124M 2s
    ##   7150K .......... .......... .......... .......... ..........  5% 93.3M 2s
    ##   7200K .......... .......... .......... .......... ..........  5% 58.4M 2s
    ##   7250K .......... .......... .......... .......... ..........  5% 30.4M 2s
    ##   7300K .......... .......... .......... .......... ..........  5% 85.4M 2s
    ##   7350K .......... .......... .......... .......... ..........  5% 62.7M 2s
    ##   7400K .......... .......... .......... .......... ..........  5%  109M 2s
    ##   7450K .......... .......... .......... .......... ..........  5% 58.4M 2s
    ##   7500K .......... .......... .......... .......... ..........  5% 79.4M 2s
    ##   7550K .......... .......... .......... .......... ..........  5% 86.0M 2s
    ##   7600K .......... .......... .......... .......... ..........  5% 89.7M 2s
    ##   7650K .......... .......... .......... .......... ..........  5% 75.7M 2s
    ##   7700K .......... .......... .......... .......... ..........  5%  100M 2s
    ##   7750K .......... .......... .......... .......... ..........  5% 75.0M 2s
    ##   7800K .......... .......... .......... .......... ..........  5% 86.0M 2s
    ##   7850K .......... .......... .......... .......... ..........  5% 75.0M 2s
    ##   7900K .......... .......... .......... .......... ..........  5% 12.8M 2s
    ##   7950K .......... .......... .......... .......... ..........  5% 70.0M 2s
    ##   8000K .......... .......... .......... .......... ..........  5%  119M 2s
    ##   8050K .......... .......... .......... .......... ..........  6%  100M 2s
    ##   8100K .......... .......... .......... .......... ..........  6%  111M 2s
    ##   8150K .......... .......... .......... .......... ..........  6% 89.4M 2s
    ##   8200K .......... .......... .......... .......... ..........  6% 82.7M 2s
    ##   8250K .......... .......... .......... .......... ..........  6% 76.8M 2s
    ##   8300K .......... .......... .......... .......... ..........  6% 99.3M 2s
    ##   8350K .......... .......... .......... .......... ..........  6% 75.6M 2s
    ##   8400K .......... .......... .......... .......... ..........  6% 54.6M 2s
    ##   8450K .......... .......... .......... .......... ..........  6% 65.3M 2s
    ##   8500K .......... .......... .......... .......... ..........  6% 66.7M 2s
    ##   8550K .......... .......... .......... .......... ..........  6%  101M 2s
    ##   8600K .......... .......... .......... .......... ..........  6%  126M 2s
    ##   8650K .......... .......... .......... .......... ..........  6% 93.7M 2s
    ##   8700K .......... .......... .......... .......... ..........  6% 99.5M 2s
    ##   8750K .......... .......... .......... .......... ..........  6% 70.7M 2s
    ##   8800K .......... .......... .......... .......... ..........  6% 31.6M 2s
    ##   8850K .......... .......... .......... .......... ..........  6% 73.8M 2s
    ##   8900K .......... .......... .......... .......... ..........  6%  124M 2s
    ##   8950K .......... .......... .......... .......... ..........  6%  100M 2s
    ##   9000K .......... .......... .......... .......... ..........  6% 30.6M 2s
    ##   9050K .......... .......... .......... .......... ..........  6% 66.0M 2s
    ##   9100K .......... .......... .......... .......... ..........  6%  102M 2s
    ##   9150K .......... .......... .......... .......... ..........  6% 96.5M 2s
    ##   9200K .......... .......... .......... .......... ..........  6%  126M 2s
    ##   9250K .......... .......... .......... .......... ..........  6% 22.6M 2s
    ##   9300K .......... .......... .......... .......... ..........  6% 30.4M 2s
    ##   9350K .......... .......... .......... .......... ..........  6% 83.7M 2s
    ##   9400K .......... .......... .......... .......... ..........  7% 89.3M 2s
    ##   9450K .......... .......... .......... .......... ..........  7% 59.1M 2s
    ##   9500K .......... .......... .......... .......... ..........  7% 96.5M 2s
    ##   9550K .......... .......... .......... .......... ..........  7% 89.2M 2s
    ##   9600K .......... .......... .......... .......... ..........  7%  100M 2s
    ##   9650K .......... .......... .......... .......... ..........  7% 19.6M 2s
    ##   9700K .......... .......... .......... .......... ..........  7% 56.1M 2s
    ##   9750K .......... .......... .......... .......... ..........  7% 67.5M 2s
    ##   9800K .......... .......... .......... .......... ..........  7% 70.0M 2s
    ##   9850K .......... .......... .......... .......... ..........  7%  107M 2s
    ##   9900K .......... .......... .......... .......... ..........  7%  121M 2s
    ##   9950K .......... .......... .......... .......... ..........  7%  103M 2s
    ##  10000K .......... .......... .......... .......... ..........  7% 52.1M 2s
    ##  10050K .......... .......... .......... .......... ..........  7% 71.1M 2s
    ##  10100K .......... .......... .......... .......... ..........  7% 38.3M 2s
    ##  10150K .......... .......... .......... .......... ..........  7%  106M 2s
    ##  10200K .......... .......... .......... .......... ..........  7% 3.44M 2s
    ##  10250K .......... .......... .......... .......... ..........  7% 94.2M 2s
    ##  10300K .......... .......... .......... .......... ..........  7% 41.0M 2s
    ##  10350K .......... .......... .......... .......... ..........  7% 75.1M 2s
    ##  10400K .......... .......... .......... .......... ..........  7% 64.2M 2s
    ##  10450K .......... .......... .......... .......... ..........  7% 76.7M 2s
    ##  10500K .......... .......... .......... .......... ..........  7% 90.1M 2s
    ##  10550K .......... .......... .......... .......... ..........  7% 94.5M 2s
    ##  10600K .......... .......... .......... .......... ..........  7%  120M 2s
    ##  10650K .......... .......... .......... .......... ..........  7% 75.5M 2s
    ##  10700K .......... .......... .......... .......... ..........  7% 17.8M 2s
    ##  10750K .......... .......... .......... .......... ..........  8% 90.6M 2s
    ##  10800K .......... .......... .......... .......... ..........  8% 93.8M 2s
    ##  10850K .......... .......... .......... .......... ..........  8%  130M 2s
    ##  10900K .......... .......... .......... .......... ..........  8%  114M 2s
    ##  10950K .......... .......... .......... .......... ..........  8% 95.9M 2s
    ##  11000K .......... .......... .......... .......... ..........  8% 61.7M 2s
    ##  11050K .......... .......... .......... .......... ..........  8% 44.0M 2s
    ##  11100K .......... .......... .......... .......... ..........  8% 62.3M 2s
    ##  11150K .......... .......... .......... .......... ..........  8% 44.6M 2s
    ##  11200K .......... .......... .......... .......... ..........  8% 58.2M 2s
    ##  11250K .......... .......... .......... .......... ..........  8% 99.8M 2s
    ##  11300K .......... .......... .......... .......... ..........  8%  118M 2s
    ##  11350K .......... .......... .......... .......... ..........  8%  117M 2s
    ##  11400K .......... .......... .......... .......... ..........  8% 98.3M 2s
    ##  11450K .......... .......... .......... .......... ..........  8%  101M 2s
    ##  11500K .......... .......... .......... .......... ..........  8%  120M 2s
    ##  11550K .......... .......... .......... .......... ..........  8% 94.1M 2s
    ##  11600K .......... .......... .......... .......... ..........  8% 34.6M 2s
    ##  11650K .......... .......... .......... .......... ..........  8% 69.3M 2s
    ##  11700K .......... .......... .......... .......... ..........  8% 25.4M 2s
    ##  11750K .......... .......... .......... .......... ..........  8% 82.5M 2s
    ##  11800K .......... .......... .......... .......... ..........  8%  106M 2s
    ##  11850K .......... .......... .......... .......... ..........  8% 26.1M 2s
    ##  11900K .......... .......... .......... .......... ..........  8%  123M 2s
    ##  11950K .......... .......... .......... .......... ..........  8%  113M 2s
    ##  12000K .......... .......... .......... .......... ..........  8%  127M 2s
    ##  12050K .......... .......... .......... .......... ..........  8%  134M 2s
    ##  12100K .......... .......... .......... .......... ..........  9% 50.2M 2s
    ##  12150K .......... .......... .......... .......... ..........  9% 94.1M 2s
    ##  12200K .......... .......... .......... .......... ..........  9%  114M 2s
    ##  12250K .......... .......... .......... .......... ..........  9%  116M 2s
    ##  12300K .......... .......... .......... .......... ..........  9%  106M 2s
    ##  12350K .......... .......... .......... .......... ..........  9%  110M 2s
    ##  12400K .......... .......... .......... .......... ..........  9% 26.5M 2s
    ##  12450K .......... .......... .......... .......... ..........  9% 87.6M 2s
    ##  12500K .......... .......... .......... .......... ..........  9% 20.3M 2s
    ##  12550K .......... .......... .......... .......... ..........  9% 56.7M 2s
    ##  12600K .......... .......... .......... .......... ..........  9%  115M 2s
    ##  12650K .......... .......... .......... .......... ..........  9% 28.6M 2s
    ##  12700K .......... .......... .......... .......... ..........  9%  114M 2s
    ##  12750K .......... .......... .......... .......... ..........  9%  114M 2s
    ##  12800K .......... .......... .......... .......... ..........  9% 94.4M 2s
    ##  12850K .......... .......... .......... .......... ..........  9%  105M 2s
    ##  12900K .......... .......... .......... .......... ..........  9%  103M 2s
    ##  12950K .......... .......... .......... .......... ..........  9%  118M 2s
    ##  13000K .......... .......... .......... .......... ..........  9%  101M 2s
    ##  13050K .......... .......... .......... .......... ..........  9%  118M 2s
    ##  13100K .......... .......... .......... .......... ..........  9% 32.5M 2s
    ##  13150K .......... .......... .......... .......... ..........  9% 71.3M 2s
    ##  13200K .......... .......... .......... .......... ..........  9% 91.6M 2s
    ##  13250K .......... .......... .......... .......... ..........  9%  101M 2s
    ##  13300K .......... .......... .......... .......... ..........  9% 79.2M 2s
    ##  13350K .......... .......... .......... .......... ..........  9% 19.0M 2s
    ##  13400K .......... .......... .......... .......... ..........  9%  117M 2s
    ##  13450K .......... .......... .......... .......... .......... 10%  106M 2s
    ##  13500K .......... .......... .......... .......... .......... 10%  134M 2s
    ##  13550K .......... .......... .......... .......... .......... 10% 94.8M 2s
    ##  13600K .......... .......... .......... .......... .......... 10%  123M 2s
    ##  13650K .......... .......... .......... .......... .......... 10%  121M 2s
    ##  13700K .......... .......... .......... .......... .......... 10%  127M 2s
    ##  13750K .......... .......... .......... .......... .......... 10%  112M 2s
    ##  13800K .......... .......... .......... .......... .......... 10% 17.9M 2s
    ##  13850K .......... .......... .......... .......... .......... 10%  105M 2s
    ##  13900K .......... .......... .......... .......... .......... 10% 90.1M 2s
    ##  13950K .......... .......... .......... .......... .......... 10% 72.9M 2s
    ##  14000K .......... .......... .......... .......... .......... 10% 82.4M 2s
    ##  14050K .......... .......... .......... .......... .......... 10%  137M 2s
    ##  14100K .......... .......... .......... .......... .......... 10% 9.94M 2s
    ##  14150K .......... .......... .......... .......... .......... 10% 94.8M 2s
    ##  14200K .......... .......... .......... .......... .......... 10%  124M 2s
    ##  14250K .......... .......... .......... .......... .......... 10%  123M 2s
    ##  14300K .......... .......... .......... .......... .......... 10%  110M 2s
    ##  14350K .......... .......... .......... .......... .......... 10%  140M 2s
    ##  14400K .......... .......... .......... .......... .......... 10%  127M 2s
    ##  14450K .......... .......... .......... .......... .......... 10%  139M 2s
    ##  14500K .......... .......... .......... .......... .......... 10% 50.5M 2s
    ##  14550K .......... .......... .......... .......... .......... 10% 42.2M 2s
    ##  14600K .......... .......... .......... .......... .......... 10% 32.3M 2s
    ##  14650K .......... .......... .......... .......... .......... 10%  105M 2s
    ##  14700K .......... .......... .......... .......... .......... 10% 79.2M 2s
    ##  14750K .......... .......... .......... .......... .......... 10%  132M 2s
    ##  14800K .......... .......... .......... .......... .......... 11% 59.7M 2s
    ##  14850K .......... .......... .......... .......... .......... 11% 52.7M 2s
    ##  14900K .......... .......... .......... .......... .......... 11% 87.8M 2s
    ##  14950K .......... .......... .......... .......... .......... 11%  114M 2s
    ##  15000K .......... .......... .......... .......... .......... 11%  128M 2s
    ##  15050K .......... .......... .......... .......... .......... 11%  123M 2s
    ##  15100K .......... .......... .......... .......... .......... 11%  101M 2s
    ##  15150K .......... .......... .......... .......... .......... 11% 32.3M 2s
    ##  15200K .......... .......... .......... .......... .......... 11% 74.3M 2s
    ##  15250K .......... .......... .......... .......... .......... 11%  122M 2s
    ##  15300K .......... .......... .......... .......... .......... 11% 43.8M 2s
    ##  15350K .......... .......... .......... .......... .......... 11% 65.0M 2s
    ##  15400K .......... .......... .......... .......... .......... 11%  107M 2s
    ##  15450K .......... .......... .......... .......... .......... 11% 59.1M 2s
    ##  15500K .......... .......... .......... .......... .......... 11% 36.1M 2s
    ##  15550K .......... .......... .......... .......... .......... 11% 90.1M 2s
    ##  15600K .......... .......... .......... .......... .......... 11% 83.8M 2s
    ##  15650K .......... .......... .......... .......... .......... 11% 55.8M 2s
    ##  15700K .......... .......... .......... .......... .......... 11% 48.1M 2s
    ##  15750K .......... .......... .......... .......... .......... 11% 80.1M 2s
    ##  15800K .......... .......... .......... .......... .......... 11%  113M 2s
    ##  15850K .......... .......... .......... .......... .......... 11% 73.3M 2s
    ##  15900K .......... .......... .......... .......... .......... 11% 46.5M 2s
    ##  15950K .......... .......... .......... .......... .......... 11%  101M 2s
    ##  16000K .......... .......... .......... .......... .......... 11%  118M 2s
    ##  16050K .......... .......... .......... .......... .......... 11% 54.9M 2s
    ##  16100K .......... .......... .......... .......... .......... 11% 44.8M 2s
    ##  16150K .......... .......... .......... .......... .......... 12% 69.2M 2s
    ##  16200K .......... .......... .......... .......... .......... 12% 99.6M 2s
    ##  16250K .......... .......... .......... .......... .......... 12% 43.0M 2s
    ##  16300K .......... .......... .......... .......... .......... 12% 85.6M 2s
    ##  16350K .......... .......... .......... .......... .......... 12% 98.1M 2s
    ##  16400K .......... .......... .......... .......... .......... 12% 97.4M 2s
    ##  16450K .......... .......... .......... .......... .......... 12% 55.2M 2s
    ##  16500K .......... .......... .......... .......... .......... 12% 58.2M 2s
    ##  16550K .......... .......... .......... .......... .......... 12% 68.3M 2s
    ##  16600K .......... .......... .......... .......... .......... 12% 98.0M 2s
    ##  16650K .......... .......... .......... .......... .......... 12% 83.6M 2s
    ##  16700K .......... .......... .......... .......... .......... 12% 92.3M 2s
    ##  16750K .......... .......... .......... .......... .......... 12%  107M 2s
    ##  16800K .......... .......... .......... .......... .......... 12% 63.3M 2s
    ##  16850K .......... .......... .......... .......... .......... 12% 40.8M 2s
    ##  16900K .......... .......... .......... .......... .......... 12% 27.3M 2s
    ##  16950K .......... .......... .......... .......... .......... 12% 63.0M 2s
    ##  17000K .......... .......... .......... .......... .......... 12% 61.2M 2s
    ##  17050K .......... .......... .......... .......... .......... 12% 34.2M 2s
    ##  17100K .......... .......... .......... .......... .......... 12% 66.2M 2s
    ##  17150K .......... .......... .......... .......... .......... 12% 83.7M 2s
    ##  17200K .......... .......... .......... .......... .......... 12% 99.3M 2s
    ##  17250K .......... .......... .......... .......... .......... 12%  103M 2s
    ##  17300K .......... .......... .......... .......... .......... 12% 78.9M 2s
    ##  17350K .......... .......... .......... .......... .......... 12% 64.6M 2s
    ##  17400K .......... .......... .......... .......... .......... 12% 54.2M 2s
    ##  17450K .......... .......... .......... .......... .......... 12% 56.6M 2s
    ##  17500K .......... .......... .......... .......... .......... 13% 50.4M 2s
    ##  17550K .......... .......... .......... .......... .......... 13% 56.0M 2s
    ##  17600K .......... .......... .......... .......... .......... 13% 85.0M 2s
    ##  17650K .......... .......... .......... .......... .......... 13%  120M 2s
    ##  17700K .......... .......... .......... .......... .......... 13% 42.4M 2s
    ##  17750K .......... .......... .......... .......... .......... 13% 87.1M 2s
    ##  17800K .......... .......... .......... .......... .......... 13% 79.8M 2s
    ##  17850K .......... .......... .......... .......... .......... 13% 42.2M 2s
    ##  17900K .......... .......... .......... .......... .......... 13% 57.8M 2s
    ##  17950K .......... .......... .......... .......... .......... 13% 65.0M 2s
    ##  18000K .......... .......... .......... .......... .......... 13% 93.8M 2s
    ##  18050K .......... .......... .......... .......... .......... 13% 39.7M 2s
    ##  18100K .......... .......... .......... .......... .......... 13% 79.0M 2s
    ##  18150K .......... .......... .......... .......... .......... 13% 76.6M 2s
    ##  18200K .......... .......... .......... .......... .......... 13% 92.1M 2s
    ##  18250K .......... .......... .......... .......... .......... 13% 41.9M 2s
    ##  18300K .......... .......... .......... .......... .......... 13% 59.0M 2s
    ##  18350K .......... .......... .......... .......... .......... 13% 82.4M 2s
    ##  18400K .......... .......... .......... .......... .......... 13% 84.2M 2s
    ##  18450K .......... .......... .......... .......... .......... 13% 40.7M 2s
    ##  18500K .......... .......... .......... .......... .......... 13% 89.3M 2s
    ##  18550K .......... .......... .......... .......... .......... 13% 81.9M 2s
    ##  18600K .......... .......... .......... .......... .......... 13% 66.9M 2s
    ##  18650K .......... .......... .......... .......... .......... 13% 72.6M 2s
    ##  18700K .......... .......... .......... .......... .......... 13% 88.0M 2s
    ##  18750K .......... .......... .......... .......... .......... 13% 57.2M 2s
    ##  18800K .......... .......... .......... .......... .......... 13% 43.1M 2s
    ##  18850K .......... .......... .......... .......... .......... 14% 81.3M 2s
    ##  18900K .......... .......... .......... .......... .......... 14%  117M 2s
    ##  18950K .......... .......... .......... .......... .......... 14% 65.0M 2s
    ##  19000K .......... .......... .......... .......... .......... 14% 45.9M 2s
    ##  19050K .......... .......... .......... .......... .......... 14% 79.1M 2s
    ##  19100K .......... .......... .......... .......... .......... 14%  130M 2s
    ##  19150K .......... .......... .......... .......... .......... 14% 39.8M 2s
    ##  19200K .......... .......... .......... .......... .......... 14%  116M 2s
    ##  19250K .......... .......... .......... .......... .......... 14%  110M 2s
    ##  19300K .......... .......... .......... .......... .......... 14% 83.4M 2s
    ##  19350K .......... .......... .......... .......... .......... 14% 57.3M 2s
    ##  19400K .......... .......... .......... .......... .......... 14% 48.9M 2s
    ##  19450K .......... .......... .......... .......... .......... 14%  116M 2s
    ##  19500K .......... .......... .......... .......... .......... 14% 55.9M 2s
    ##  19550K .......... .......... .......... .......... .......... 14% 95.9M 2s
    ##  19600K .......... .......... .......... .......... .......... 14%  118M 2s
    ##  19650K .......... .......... .......... .......... .......... 14% 69.5M 2s
    ##  19700K .......... .......... .......... .......... .......... 14% 91.9M 2s
    ##  19750K .......... .......... .......... .......... .......... 14% 52.2M 2s
    ##  19800K .......... .......... .......... .......... .......... 14% 74.3M 2s
    ##  19850K .......... .......... .......... .......... .......... 14% 96.5M 2s
    ##  19900K .......... .......... .......... .......... .......... 14% 48.3M 2s
    ##  19950K .......... .......... .......... .......... .......... 14%  106M 2s
    ##  20000K .......... .......... .......... .......... .......... 14% 93.6M 2s
    ##  20050K .......... .......... .......... .......... .......... 14% 75.5M 2s
    ##  20100K .......... .......... .......... .......... .......... 14%  102M 2s
    ##  20150K .......... .......... .......... .......... .......... 14% 56.9M 2s
    ##  20200K .......... .......... .......... .......... .......... 15% 88.4M 2s
    ##  20250K .......... .......... .......... .......... .......... 15%  140M 2s
    ##  20300K .......... .......... .......... .......... .......... 15% 73.6M 2s
    ##  20350K .......... .......... .......... .......... .......... 15% 98.8M 2s
    ##  20400K .......... .......... .......... .......... .......... 15%  104M 2s
    ##  20450K .......... .......... .......... .......... .......... 15% 1021K 2s
    ##  20500K .......... .......... .......... .......... .......... 15% 21.5M 2s
    ##  20550K .......... .......... .......... .......... .......... 15% 49.0M 2s
    ##  20600K .......... .......... .......... .......... .......... 15% 39.1M 2s
    ##  20650K .......... .......... .......... .......... .......... 15% 32.3M 2s
    ##  20700K .......... .......... .......... .......... .......... 15% 85.6M 2s
    ##  20750K .......... .......... .......... .......... .......... 15% 16.9M 2s
    ##  20800K .......... .......... .......... .......... .......... 15%  104M 2s
    ##  20850K .......... .......... .......... .......... .......... 15%  108M 2s
    ##  20900K .......... .......... .......... .......... .......... 15%  108M 2s
    ##  20950K .......... .......... .......... .......... .......... 15% 74.5M 2s
    ##  21000K .......... .......... .......... .......... .......... 15% 75.1M 2s
    ##  21050K .......... .......... .......... .......... .......... 15%  114M 2s
    ##  21100K .......... .......... .......... .......... .......... 15% 93.1M 2s
    ##  21150K .......... .......... .......... .......... .......... 15%  115M 2s
    ##  21200K .......... .......... .......... .......... .......... 15% 79.7M 2s
    ##  21250K .......... .......... .......... .......... .......... 15% 42.7M 2s
    ##  21300K .......... .......... .......... .......... .......... 15% 60.3M 2s
    ##  21350K .......... .......... .......... .......... .......... 15% 34.4M 2s
    ##  21400K .......... .......... .......... .......... .......... 15% 94.8M 2s
    ##  21450K .......... .......... .......... .......... .......... 15% 29.2M 2s
    ##  21500K .......... .......... .......... .......... .......... 15% 80.3M 2s
    ##  21550K .......... .......... .......... .......... .......... 16% 73.3M 2s
    ##  21600K .......... .......... .......... .......... .......... 16% 89.1M 2s
    ##  21650K .......... .......... .......... .......... .......... 16% 92.1M 2s
    ##  21700K .......... .......... .......... .......... .......... 16% 93.0M 2s
    ##  21750K .......... .......... .......... .......... .......... 16% 89.3M 2s
    ##  21800K .......... .......... .......... .......... .......... 16%  100M 2s
    ##  21850K .......... .......... .......... .......... .......... 16% 96.5M 2s
    ##  21900K .......... .......... .......... .......... .......... 16%  105M 2s
    ##  21950K .......... .......... .......... .......... .......... 16% 83.9M 2s
    ##  22000K .......... .......... .......... .......... .......... 16% 45.6M 2s
    ##  22050K .......... .......... .......... .......... .......... 16% 61.6M 2s
    ##  22100K .......... .......... .......... .......... .......... 16% 77.7M 2s
    ##  22150K .......... .......... .......... .......... .......... 16% 80.3M 2s
    ##  22200K .......... .......... .......... .......... .......... 16%  107M 2s
    ##  22250K .......... .......... .......... .......... .......... 16% 20.0M 2s
    ##  22300K .......... .......... .......... .......... .......... 16%  107M 2s
    ##  22350K .......... .......... .......... .......... .......... 16%  110M 2s
    ##  22400K .......... .......... .......... .......... .......... 16% 82.6M 2s
    ##  22450K .......... .......... .......... .......... .......... 16% 94.3M 2s
    ##  22500K .......... .......... .......... .......... .......... 16% 84.8M 2s
    ##  22550K .......... .......... .......... .......... .......... 16% 92.7M 2s
    ##  22600K .......... .......... .......... .......... .......... 16%  116M 2s
    ##  22650K .......... .......... .......... .......... .......... 16%  109M 2s
    ##  22700K .......... .......... .......... .......... .......... 16% 6.83M 2s
    ##  22750K .......... .......... .......... .......... .......... 16% 93.2M 2s
    ##  22800K .......... .......... .......... .......... .......... 16% 35.4M 2s
    ##  22850K .......... .......... .......... .......... .......... 16% 65.8M 2s
    ##  22900K .......... .......... .......... .......... .......... 17% 71.4M 2s
    ##  22950K .......... .......... .......... .......... .......... 17% 86.8M 2s
    ##  23000K .......... .......... .......... .......... .......... 17%  101M 2s
    ##  23050K .......... .......... .......... .......... .......... 17% 10.7M 2s
    ##  23100K .......... .......... .......... .......... .......... 17% 93.2M 2s
    ##  23150K .......... .......... .......... .......... .......... 17% 88.8M 2s
    ##  23200K .......... .......... .......... .......... .......... 17%  104M 2s
    ##  23250K .......... .......... .......... .......... .......... 17%  124M 2s
    ##  23300K .......... .......... .......... .......... .......... 17%  101M 2s
    ##  23350K .......... .......... .......... .......... .......... 17% 77.1M 2s
    ##  23400K .......... .......... .......... .......... .......... 17% 92.3M 2s
    ##  23450K .......... .......... .......... .......... .......... 17%  100M 2s
    ##  23500K .......... .......... .......... .......... .......... 17%  100M 2s
    ##  23550K .......... .......... .......... .......... .......... 17%  130M 2s
    ##  23600K .......... .......... .......... .......... .......... 17% 39.5M 2s
    ##  23650K .......... .......... .......... .......... .......... 17% 80.0M 2s
    ##  23700K .......... .......... .......... .......... .......... 17% 60.0M 2s
    ##  23750K .......... .......... .......... .......... .......... 17% 31.0M 2s
    ##  23800K .......... .......... .......... .......... .......... 17%  103M 2s
    ##  23850K .......... .......... .......... .......... .......... 17%  121M 2s
    ##  23900K .......... .......... .......... .......... .......... 17%  101M 2s
    ##  23950K .......... .......... .......... .......... .......... 17% 81.3M 2s
    ##  24000K .......... .......... .......... .......... .......... 17% 90.2M 2s
    ##  24050K .......... .......... .......... .......... .......... 17% 81.5M 2s
    ##  24100K .......... .......... .......... .......... .......... 17% 95.2M 2s
    ##  24150K .......... .......... .......... .......... .......... 17%  102M 2s
    ##  24200K .......... .......... .......... .......... .......... 17% 24.9M 2s
    ##  24250K .......... .......... .......... .......... .......... 18% 52.5M 2s
    ##  24300K .......... .......... .......... .......... .......... 18%  105M 2s
    ##  24350K .......... .......... .......... .......... .......... 18%  119M 2s
    ##  24400K .......... .......... .......... .......... .......... 18% 54.4M 2s
    ##  24450K .......... .......... .......... .......... .......... 18%  110M 2s
    ##  24500K .......... .......... .......... .......... .......... 18%  114M 2s
    ##  24550K .......... .......... .......... .......... .......... 18%  133M 2s
    ##  24600K .......... .......... .......... .......... .......... 18% 33.3M 2s
    ##  24650K .......... .......... .......... .......... .......... 18%  112M 2s
    ##  24700K .......... .......... .......... .......... .......... 18% 89.5M 2s
    ##  24750K .......... .......... .......... .......... .......... 18% 90.3M 2s
    ##  24800K .......... .......... .......... .......... .......... 18% 64.7M 2s
    ##  24850K .......... .......... .......... .......... .......... 18% 88.6M 2s
    ##  24900K .......... .......... .......... .......... .......... 18% 88.6M 2s
    ##  24950K .......... .......... .......... .......... .......... 18% 22.8M 2s
    ##  25000K .......... .......... .......... .......... .......... 18% 61.7M 2s
    ##  25050K .......... .......... .......... .......... .......... 18% 38.1M 2s
    ##  25100K .......... .......... .......... .......... .......... 18% 63.6M 2s
    ##  25150K .......... .......... .......... .......... .......... 18% 53.6M 2s
    ##  25200K .......... .......... .......... .......... .......... 18% 65.4M 2s
    ##  25250K .......... .......... .......... .......... .......... 18%  110M 2s
    ##  25300K .......... .......... .......... .......... .......... 18%  114M 2s
    ##  25350K .......... .......... .......... .......... .......... 18% 10.8M 2s
    ##  25400K .......... .......... .......... .......... .......... 18% 58.4M 2s
    ##  25450K .......... .......... .......... .......... .......... 18% 38.8M 2s
    ##  25500K .......... .......... .......... .......... .......... 18%  103M 2s
    ##  25550K .......... .......... .......... .......... .......... 18%  115M 2s
    ##  25600K .......... .......... .......... .......... .......... 19%  212K 3s
    ##  25650K .......... .......... .......... .......... .......... 19% 82.6M 3s
    ##  25700K .......... .......... .......... .......... .......... 19% 68.5M 3s
    ##  25750K .......... .......... .......... .......... .......... 19% 86.2M 3s
    ##  25800K .......... .......... .......... .......... .......... 19% 83.5M 3s
    ##  25850K .......... .......... .......... .......... .......... 19% 9.10M 3s
    ##  25900K .......... .......... .......... .......... .......... 19% 57.7M 3s
    ##  25950K .......... .......... .......... .......... .......... 19% 82.4M 3s
    ##  26000K .......... .......... .......... .......... .......... 19% 84.8M 3s
    ##  26050K .......... .......... .......... .......... .......... 19%  106M 3s
    ##  26100K .......... .......... .......... .......... .......... 19% 88.1M 3s
    ##  26150K .......... .......... .......... .......... .......... 19% 48.1M 3s
    ##  26200K .......... .......... .......... .......... .......... 19% 55.0M 3s
    ##  26250K .......... .......... .......... .......... .......... 19% 42.2M 3s
    ##  26300K .......... .......... .......... .......... .......... 19% 53.7M 3s
    ##  26350K .......... .......... .......... .......... .......... 19% 80.1M 3s
    ##  26400K .......... .......... .......... .......... .......... 19% 64.7M 3s
    ##  26450K .......... .......... .......... .......... .......... 19% 89.0M 3s
    ##  26500K .......... .......... .......... .......... .......... 19% 74.8M 3s
    ##  26550K .......... .......... .......... .......... .......... 19% 69.4M 3s
    ##  26600K .......... .......... .......... .......... .......... 19% 80.1M 3s
    ##  26650K .......... .......... .......... .......... .......... 19% 27.2M 3s
    ##  26700K .......... .......... .......... .......... .......... 19% 44.8M 3s
    ##  26750K .......... .......... .......... .......... .......... 19%  101M 3s
    ##  26800K .......... .......... .......... .......... .......... 19% 80.9M 3s
    ##  26850K .......... .......... .......... .......... .......... 19% 86.8M 3s
    ##  26900K .......... .......... .......... .......... .......... 20% 85.3M 3s
    ##  26950K .......... .......... .......... .......... .......... 20% 16.3M 3s
    ##  27000K .......... .......... .......... .......... .......... 20% 70.8M 3s
    ##  27050K .......... .......... .......... .......... .......... 20% 37.3M 3s
    ##  27100K .......... .......... .......... .......... .......... 20% 50.4M 3s
    ##  27150K .......... .......... .......... .......... .......... 20% 78.5M 3s
    ##  27200K .......... .......... .......... .......... .......... 20% 61.9M 3s
    ##  27250K .......... .......... .......... .......... .......... 20% 73.1M 3s
    ##  27300K .......... .......... .......... .......... .......... 20% 96.3M 3s
    ##  27350K .......... .......... .......... .......... .......... 20%  115M 3s
    ##  27400K .......... .......... .......... .......... .......... 20% 31.7M 3s
    ##  27450K .......... .......... .......... .......... .......... 20% 70.0M 3s
    ##  27500K .......... .......... .......... .......... .......... 20% 24.5M 3s
    ##  27550K .......... .......... .......... .......... .......... 20% 42.2M 3s
    ##  27600K .......... .......... .......... .......... .......... 20% 66.5M 3s
    ##  27650K .......... .......... .......... .......... .......... 20% 82.8M 3s
    ##  27700K .......... .......... .......... .......... .......... 20% 87.7M 3s
    ##  27750K .......... .......... .......... .......... .......... 20% 99.0M 3s
    ##  27800K .......... .......... .......... .......... .......... 20% 10.1M 3s
    ##  27850K .......... .......... .......... .......... .......... 20% 82.2M 3s
    ##  27900K .......... .......... .......... .......... .......... 20% 87.6M 3s
    ##  27950K .......... .......... .......... .......... .......... 20% 31.2M 3s
    ##  28000K .......... .......... .......... .......... .......... 20% 60.7M 3s
    ##  28050K .......... .......... .......... .......... .......... 20% 73.5M 3s
    ##  28100K .......... .......... .......... .......... .......... 20% 95.0M 3s
    ##  28150K .......... .......... .......... .......... .......... 20% 81.7M 3s
    ##  28200K .......... .......... .......... .......... .......... 20% 73.7M 3s
    ##  28250K .......... .......... .......... .......... .......... 21% 31.8M 3s
    ##  28300K .......... .......... .......... .......... .......... 21% 97.6M 3s
    ##  28350K .......... .......... .......... .......... .......... 21% 83.3M 3s
    ##  28400K .......... .......... .......... .......... .......... 21%  101M 3s
    ##  28450K .......... .......... .......... .......... .......... 21% 36.8M 3s
    ##  28500K .......... .......... .......... .......... .......... 21% 93.4M 3s
    ##  28550K .......... .......... .......... .......... .......... 21%  110M 3s
    ##  28600K .......... .......... .......... .......... .......... 21% 26.2M 3s
    ##  28650K .......... .......... .......... .......... .......... 21% 95.5M 3s
    ##  28700K .......... .......... .......... .......... .......... 21% 87.5M 3s
    ##  28750K .......... .......... .......... .......... .......... 21%  113M 3s
    ##  28800K .......... .......... .......... .......... .......... 21% 93.6M 3s
    ##  28850K .......... .......... .......... .......... .......... 21% 6.29M 3s
    ##  28900K .......... .......... .......... .......... .......... 21% 21.7M 3s
    ##  28950K .......... .......... .......... .......... .......... 21% 36.4M 3s
    ##  29000K .......... .......... .......... .......... .......... 21% 85.3M 3s
    ##  29050K .......... .......... .......... .......... .......... 21% 97.1M 3s
    ##  29100K .......... .......... .......... .......... .......... 21% 96.3M 3s
    ##  29150K .......... .......... .......... .......... .......... 21% 26.6M 3s
    ##  29200K .......... .......... .......... .......... .......... 21% 66.1M 3s
    ##  29250K .......... .......... .......... .......... .......... 21% 85.5M 3s
    ##  29300K .......... .......... .......... .......... .......... 21% 71.4M 3s
    ##  29350K .......... .......... .......... .......... .......... 21% 99.3M 3s
    ##  29400K .......... .......... .......... .......... .......... 21% 69.8M 3s
    ##  29450K .......... .......... .......... .......... .......... 21% 82.4M 3s
    ##  29500K .......... .......... .......... .......... .......... 21% 68.6M 3s
    ##  29550K .......... .......... .......... .......... .......... 21%  112M 3s
    ##  29600K .......... .......... .......... .......... .......... 22% 91.6M 3s
    ##  29650K .......... .......... .......... .......... .......... 22% 9.43M 3s
    ##  29700K .......... .......... .......... .......... .......... 22% 57.5M 3s
    ##  29750K .......... .......... .......... .......... .......... 22% 20.5M 3s
    ##  29800K .......... .......... .......... .......... .......... 22% 66.2M 3s
    ##  29850K .......... .......... .......... .......... .......... 22% 85.0M 3s
    ##  29900K .......... .......... .......... .......... .......... 22% 64.9M 3s
    ##  29950K .......... .......... .......... .......... .......... 22% 98.1M 3s
    ##  30000K .......... .......... .......... .......... .......... 22% 72.1M 3s
    ##  30050K .......... .......... .......... .......... .......... 22% 66.1M 3s
    ##  30100K .......... .......... .......... .......... .......... 22% 68.9M 3s
    ##  30150K .......... .......... .......... .......... .......... 22%  123M 3s
    ##  30200K .......... .......... .......... .......... .......... 22% 46.9M 3s
    ##  30250K .......... .......... .......... .......... .......... 22% 45.3M 3s
    ##  30300K .......... .......... .......... .......... .......... 22% 61.0M 3s
    ##  30350K .......... .......... .......... .......... .......... 22% 85.0M 3s
    ##  30400K .......... .......... .......... .......... .......... 22% 64.3M 3s
    ##  30450K .......... .......... .......... .......... .......... 22% 30.6M 3s
    ##  30500K .......... .......... .......... .......... .......... 22% 88.8M 3s
    ##  30550K .......... .......... .......... .......... .......... 22% 77.5M 3s
    ##  30600K .......... .......... .......... .......... .......... 22% 78.3M 3s
    ##  30650K .......... .......... .......... .......... .......... 22% 83.0M 3s
    ##  30700K .......... .......... .......... .......... .......... 22% 89.3M 3s
    ##  30750K .......... .......... .......... .......... .......... 22%  108M 3s
    ##  30800K .......... .......... .......... .......... .......... 22% 94.5M 3s
    ##  30850K .......... .......... .......... .......... .......... 22% 64.9M 3s
    ##  30900K .......... .......... .......... .......... .......... 22% 45.2M 3s
    ##  30950K .......... .......... .......... .......... .......... 23% 86.9M 3s
    ##  31000K .......... .......... .......... .......... .......... 23% 89.9M 3s
    ##  31050K .......... .......... .......... .......... .......... 23% 43.3M 3s
    ##  31100K .......... .......... .......... .......... .......... 23% 30.5M 3s
    ##  31150K .......... .......... .......... .......... .......... 23% 93.5M 3s
    ##  31200K .......... .......... .......... .......... .......... 23% 97.9M 3s
    ##  31250K .......... .......... .......... .......... .......... 23% 37.7M 3s
    ##  31300K .......... .......... .......... .......... .......... 23% 19.2M 3s
    ##  31350K .......... .......... .......... .......... .......... 23% 84.8M 3s
    ##  31400K .......... .......... .......... .......... .......... 23% 67.0M 3s
    ##  31450K .......... .......... .......... .......... .......... 23% 93.5M 3s
    ##  31500K .......... .......... .......... .......... .......... 23% 46.8M 3s
    ##  31550K .......... .......... .......... .......... .......... 23% 87.6M 3s
    ##  31600K .......... .......... .......... .......... .......... 23% 68.5M 3s
    ##  31650K .......... .......... .......... .......... .......... 23% 36.3M 3s
    ##  31700K .......... .......... .......... .......... .......... 23% 14.3M 3s
    ##  31750K .......... .......... .......... .......... .......... 23% 86.3M 3s
    ##  31800K .......... .......... .......... .......... .......... 23%  105M 3s
    ##  31850K .......... .......... .......... .......... .......... 23%  110M 3s
    ##  31900K .......... .......... .......... .......... .......... 23%  113M 3s
    ##  31950K .......... .......... .......... .......... .......... 23%  139M 3s
    ##  32000K .......... .......... .......... .......... .......... 23%  114M 3s
    ##  32050K .......... .......... .......... .......... .......... 23%  120M 3s
    ##  32100K .......... .......... .......... .......... .......... 23%  110M 3s
    ##  32150K .......... .......... .......... .......... .......... 23%  102M 3s
    ##  32200K .......... .......... .......... .......... .......... 23% 36.2M 3s
    ##  32250K .......... .......... .......... .......... .......... 23% 82.0M 3s
    ##  32300K .......... .......... .......... .......... .......... 24% 27.3M 3s
    ##  32350K .......... .......... .......... .......... .......... 24% 94.3M 3s
    ##  32400K .......... .......... .......... .......... .......... 24% 81.7M 3s
    ##  32450K .......... .......... .......... .......... .......... 24% 83.9M 3s
    ##  32500K .......... .......... .......... .......... .......... 24% 75.7M 3s
    ##  32550K .......... .......... .......... .......... .......... 24% 11.4M 3s
    ##  32600K .......... .......... .......... .......... .......... 24% 98.7M 3s
    ##  32650K .......... .......... .......... .......... .......... 24%  116M 3s
    ##  32700K .......... .......... .......... .......... .......... 24%  109M 3s
    ##  32750K .......... .......... .......... .......... .......... 24%  138M 3s
    ##  32800K .......... .......... .......... .......... .......... 24%  114M 3s
    ##  32850K .......... .......... .......... .......... .......... 24%  146M 3s
    ##  32900K .......... .......... .......... .......... .......... 24% 45.0M 3s
    ##  32950K .......... .......... .......... .......... .......... 24% 80.0M 3s
    ##  33000K .......... .......... .......... .......... .......... 24% 67.5M 3s
    ##  33050K .......... .......... .......... .......... .......... 24% 76.7M 3s
    ##  33100K .......... .......... .......... .......... .......... 24% 67.8M 3s
    ##  33150K .......... .......... .......... .......... .......... 24% 87.6M 3s
    ##  33200K .......... .......... .......... .......... .......... 24% 34.3M 3s
    ##  33250K .......... .......... .......... .......... .......... 24% 51.0M 3s
    ##  33300K .......... .......... .......... .......... .......... 24% 69.0M 3s
    ##  33350K .......... .......... .......... .......... .......... 24% 94.1M 3s
    ##  33400K .......... .......... .......... .......... .......... 24% 78.5M 3s
    ##  33450K .......... .......... .......... .......... .......... 24%  109M 3s
    ##  33500K .......... .......... .......... .......... .......... 24% 64.5M 3s
    ##  33550K .......... .......... .......... .......... .......... 24% 77.6M 3s
    ##  33600K .......... .......... .......... .......... .......... 24% 73.4M 3s
    ##  33650K .......... .......... .......... .......... .......... 25% 82.1M 3s
    ##  33700K .......... .......... .......... .......... .......... 25% 73.4M 3s
    ##  33750K .......... .......... .......... .......... .......... 25% 20.3M 3s
    ##  33800K .......... .......... .......... .......... .......... 25% 81.5M 3s
    ##  33850K .......... .......... .......... .......... .......... 25% 87.4M 3s
    ##  33900K .......... .......... .......... .......... .......... 25% 48.8M 3s
    ##  33950K .......... .......... .......... .......... .......... 25% 95.2M 3s
    ##  34000K .......... .......... .......... .......... .......... 25% 84.2M 3s
    ##  34050K .......... .......... .......... .......... .......... 25% 82.6M 3s
    ##  34100K .......... .......... .......... .......... .......... 25% 65.1M 3s
    ##  34150K .......... .......... .......... .......... .......... 25% 69.6M 3s
    ##  34200K .......... .......... .......... .......... .......... 25% 76.0M 3s
    ##  34250K .......... .......... .......... .......... .......... 25%  107M 3s
    ##  34300K .......... .......... .......... .......... .......... 25% 81.9M 3s
    ##  34350K .......... .......... .......... .......... .......... 25%  106M 3s
    ##  34400K .......... .......... .......... .......... .......... 25% 73.3M 3s
    ##  34450K .......... .......... .......... .......... .......... 25% 97.8M 3s
    ##  34500K .......... .......... .......... .......... .......... 25% 11.4M 3s
    ##  34550K .......... .......... .......... .......... .......... 25% 80.0M 3s
    ##  34600K .......... .......... .......... .......... .......... 25% 78.9M 3s
    ##  34650K .......... .......... .......... .......... .......... 25%  116M 3s
    ##  34700K .......... .......... .......... .......... .......... 25% 93.5M 3s
    ##  34750K .......... .......... .......... .......... .......... 25% 93.5M 3s
    ##  34800K .......... .......... .......... .......... .......... 25% 95.3M 3s
    ##  34850K .......... .......... .......... .......... .......... 25%  108M 3s
    ##  34900K .......... .......... .......... .......... .......... 25% 79.9M 3s
    ##  34950K .......... .......... .......... .......... .......... 25%  113M 3s
    ##  35000K .......... .......... .......... .......... .......... 26% 69.0M 3s
    ##  35050K .......... .......... .......... .......... .......... 26% 95.3M 3s
    ##  35100K .......... .......... .......... .......... .......... 26% 38.0M 3s
    ##  35150K .......... .......... .......... .......... .......... 26% 95.6M 3s
    ##  35200K .......... .......... .......... .......... .......... 26% 70.7M 3s
    ##  35250K .......... .......... .......... .......... .......... 26% 75.5M 3s
    ##  35300K .......... .......... .......... .......... .......... 26% 81.1M 3s
    ##  35350K .......... .......... .......... .......... .......... 26% 96.3M 3s
    ##  35400K .......... .......... .......... .......... .......... 26% 70.1M 3s
    ##  35450K .......... .......... .......... .......... .......... 26% 32.1M 3s
    ##  35500K .......... .......... .......... .......... .......... 26% 37.0M 3s
    ##  35550K .......... .......... .......... .......... .......... 26%  111M 3s
    ##  35600K .......... .......... .......... .......... .......... 26% 84.8M 3s
    ##  35650K .......... .......... .......... .......... .......... 26% 99.1M 3s
    ##  35700K .......... .......... .......... .......... .......... 26% 88.6M 3s
    ##  35750K .......... .......... .......... .......... .......... 26%  118M 3s
    ##  35800K .......... .......... .......... .......... .......... 26%  255K 3s
    ##  35850K .......... .......... .......... .......... .......... 26% 48.4M 3s
    ##  35900K .......... .......... .......... .......... .......... 26% 34.7M 3s
    ##  35950K .......... .......... .......... .......... .......... 26% 32.9M 3s
    ##  36000K .......... .......... .......... .......... .......... 26% 75.8M 3s
    ##  36050K .......... .......... .......... .......... .......... 26% 26.3M 3s
    ##  36100K .......... .......... .......... .......... .......... 26% 43.3M 3s
    ##  36150K .......... .......... .......... .......... .......... 26% 96.5M 3s
    ##  36200K .......... .......... .......... .......... .......... 26% 89.2M 3s
    ##  36250K .......... .......... .......... .......... .......... 26% 94.7M 3s
    ##  36300K .......... .......... .......... .......... .......... 26% 67.6M 3s
    ##  36350K .......... .......... .......... .......... .......... 27% 76.7M 3s
    ##  36400K .......... .......... .......... .......... .......... 27% 88.5M 3s
    ##  36450K .......... .......... .......... .......... .......... 27% 14.7M 3s
    ##  36500K .......... .......... .......... .......... .......... 27% 71.5M 3s
    ##  36550K .......... .......... .......... .......... .......... 27% 88.5M 3s
    ##  36600K .......... .......... .......... .......... .......... 27% 88.8M 3s
    ##  36650K .......... .......... .......... .......... .......... 27% 53.5M 3s
    ##  36700K .......... .......... .......... .......... .......... 27% 89.2M 3s
    ##  36750K .......... .......... .......... .......... .......... 27%  106M 3s
    ##  36800K .......... .......... .......... .......... .......... 27% 63.8M 3s
    ##  36850K .......... .......... .......... .......... .......... 27% 27.9M 3s
    ##  36900K .......... .......... .......... .......... .......... 27% 81.5M 3s
    ##  36950K .......... .......... .......... .......... .......... 27% 71.1M 3s
    ##  37000K .......... .......... .......... .......... .......... 27% 99.4M 3s
    ##  37050K .......... .......... .......... .......... .......... 27%  103M 3s
    ##  37100K .......... .......... .......... .......... .......... 27% 69.8M 3s
    ##  37150K .......... .......... .......... .......... .......... 27% 24.5M 3s
    ##  37200K .......... .......... .......... .......... .......... 27% 90.3M 3s
    ##  37250K .......... .......... .......... .......... .......... 27% 64.3M 3s
    ##  37300K .......... .......... .......... .......... .......... 27% 94.4M 3s
    ##  37350K .......... .......... .......... .......... .......... 27% 98.1M 3s
    ##  37400K .......... .......... .......... .......... .......... 27% 91.5M 3s
    ##  37450K .......... .......... .......... .......... .......... 27% 25.3M 3s
    ##  37500K .......... .......... .......... .......... .......... 27% 87.3M 3s
    ##  37550K .......... .......... .......... .......... .......... 27% 88.6M 3s
    ##  37600K .......... .......... .......... .......... .......... 27% 91.4M 3s
    ##  37650K .......... .......... .......... .......... .......... 27% 84.5M 3s
    ##  37700K .......... .......... .......... .......... .......... 28% 87.8M 3s
    ##  37750K .......... .......... .......... .......... .......... 28% 92.3M 3s
    ##  37800K .......... .......... .......... .......... .......... 28% 94.9M 3s
    ##  37850K .......... .......... .......... .......... .......... 28% 87.5M 3s
    ##  37900K .......... .......... .......... .......... .......... 28% 33.9M 3s
    ##  37950K .......... .......... .......... .......... .......... 28% 66.6M 3s
    ##  38000K .......... .......... .......... .......... .......... 28% 70.8M 3s
    ##  38050K .......... .......... .......... .......... .......... 28% 81.4M 3s
    ##  38100K .......... .......... .......... .......... .......... 28% 15.4M 3s
    ##  38150K .......... .......... .......... .......... .......... 28% 47.1M 3s
    ##  38200K .......... .......... .......... .......... .......... 28% 87.4M 3s
    ##  38250K .......... .......... .......... .......... .......... 28%  109M 3s
    ##  38300K .......... .......... .......... .......... .......... 28% 98.4M 3s
    ##  38350K .......... .......... .......... .......... .......... 28% 93.8M 3s
    ##  38400K .......... .......... .......... .......... .......... 28% 79.4M 3s
    ##  38450K .......... .......... .......... .......... .......... 28% 82.9M 3s
    ##  38500K .......... .......... .......... .......... .......... 28% 98.8M 3s
    ##  38550K .......... .......... .......... .......... .......... 28% 36.9M 3s
    ##  38600K .......... .......... .......... .......... .......... 28% 72.9M 3s
    ##  38650K .......... .......... .......... .......... .......... 28% 87.5M 3s
    ##  38700K .......... .......... .......... .......... .......... 28% 48.1M 3s
    ##  38750K .......... .......... .......... .......... .......... 28% 30.1M 3s
    ##  38800K .......... .......... .......... .......... .......... 28% 55.1M 3s
    ##  38850K .......... .......... .......... .......... .......... 28% 16.0M 3s
    ##  38900K .......... .......... .......... .......... .......... 28% 77.3M 3s
    ##  38950K .......... .......... .......... .......... .......... 28%  101M 3s
    ##  39000K .......... .......... .......... .......... .......... 28% 75.1M 3s
    ##  39050K .......... .......... .......... .......... .......... 29% 95.5M 3s
    ##  39100K .......... .......... .......... .......... .......... 29% 98.8M 3s
    ##  39150K .......... .......... .......... .......... .......... 29%  113M 3s
    ##  39200K .......... .......... .......... .......... .......... 29%  102M 3s
    ##  39250K .......... .......... .......... .......... .......... 29% 96.0M 3s
    ##  39300K .......... .......... .......... .......... .......... 29% 85.2M 3s
    ##  39350K .......... .......... .......... .......... .......... 29% 82.2M 3s
    ##  39400K .......... .......... .......... .......... .......... 29% 89.4M 3s
    ##  39450K .......... .......... .......... .......... .......... 29% 37.2M 3s
    ##  39500K .......... .......... .......... .......... .......... 29% 97.7M 3s
    ##  39550K .......... .......... .......... .......... .......... 29% 63.9M 3s
    ##  39600K .......... .......... .......... .......... .......... 29%  102M 3s
    ##  39650K .......... .......... .......... .......... .......... 29% 55.0M 3s
    ##  39700K .......... .......... .......... .......... .......... 29% 94.3M 3s
    ##  39750K .......... .......... .......... .......... .......... 29% 52.9M 3s
    ##  39800K .......... .......... .......... .......... .......... 29% 88.7M 3s
    ##  39850K .......... .......... .......... .......... .......... 29% 74.6M 3s
    ##  39900K .......... .......... .......... .......... .......... 29% 96.9M 3s
    ##  39950K .......... .......... .......... .......... .......... 29%  102M 3s
    ##  40000K .......... .......... .......... .......... .......... 29% 20.0M 3s
    ##  40050K .......... .......... .......... .......... .......... 29%  105M 3s
    ##  40100K .......... .......... .......... .......... .......... 29%  114M 3s
    ##  40150K .......... .......... .......... .......... .......... 29%  118M 3s
    ##  40200K .......... .......... .......... .......... .......... 29% 95.7M 3s
    ##  40250K .......... .......... .......... .......... .......... 29% 95.4M 3s
    ##  40300K .......... .......... .......... .......... .......... 29% 97.4M 3s
    ##  40350K .......... .......... .......... .......... .......... 29%  116M 3s
    ##  40400K .......... .......... .......... .......... .......... 30%  112M 3s
    ##  40450K .......... .......... .......... .......... .......... 30%  117M 3s
    ##  40500K .......... .......... .......... .......... .......... 30% 32.6M 3s
    ##  40550K .......... .......... .......... .......... .......... 30% 82.9M 3s
    ##  40600K .......... .......... .......... .......... .......... 30%  103M 3s
    ##  40650K .......... .......... .......... .......... .......... 30% 57.2M 3s
    ##  40700K .......... .......... .......... .......... .......... 30% 97.0M 3s
    ##  40750K .......... .......... .......... .......... .......... 30%  112M 3s
    ##  40800K .......... .......... .......... .......... .......... 30% 25.7M 3s
    ##  40850K .......... .......... .......... .......... .......... 30% 95.2M 3s
    ##  40900K .......... .......... .......... .......... .......... 30%  110M 3s
    ##  40950K .......... .......... .......... .......... .......... 30%  222K 3s
    ##  41000K .......... .......... .......... .......... .......... 30% 78.1M 3s
    ##  41050K .......... .......... .......... .......... .......... 30% 27.6M 3s
    ##  41100K .......... .......... .......... .......... .......... 30% 77.3M 3s
    ##  41150K .......... .......... .......... .......... .......... 30%  103M 3s
    ##  41200K .......... .......... .......... .......... .......... 30% 55.1M 3s
    ##  41250K .......... .......... .......... .......... .......... 30% 96.3M 3s
    ##  41300K .......... .......... .......... .......... .......... 30% 72.8M 3s
    ##  41350K .......... .......... .......... .......... .......... 30% 8.78M 3s
    ##  41400K .......... .......... .......... .......... .......... 30% 69.6M 3s
    ##  41450K .......... .......... .......... .......... .......... 30% 21.5M 3s
    ##  41500K .......... .......... .......... .......... .......... 30% 88.2M 3s
    ##  41550K .......... .......... .......... .......... .......... 30% 40.0M 3s
    ##  41600K .......... .......... .......... .......... .......... 30% 70.3M 3s
    ##  41650K .......... .......... .......... .......... .......... 30% 77.0M 3s
    ##  41700K .......... .......... .......... .......... .......... 30% 97.8M 3s
    ##  41750K .......... .......... .......... .......... .......... 31% 56.5M 3s
    ##  41800K .......... .......... .......... .......... .......... 31% 95.5M 3s
    ##  41850K .......... .......... .......... .......... .......... 31% 41.3M 3s
    ##  41900K .......... .......... .......... .......... .......... 31% 60.7M 3s
    ##  41950K .......... .......... .......... .......... .......... 31% 84.5M 3s
    ##  42000K .......... .......... .......... .......... .......... 31%  102M 3s
    ##  42050K .......... .......... .......... .......... .......... 31% 87.5M 3s
    ##  42100K .......... .......... .......... .......... .......... 31% 83.8M 3s
    ##  42150K .......... .......... .......... .......... .......... 31%  112M 3s
    ##  42200K .......... .......... .......... .......... .......... 31% 27.8M 3s
    ##  42250K .......... .......... .......... .......... .......... 31% 36.8M 3s
    ##  42300K .......... .......... .......... .......... .......... 31% 86.7M 3s
    ##  42350K .......... .......... .......... .......... .......... 31% 82.7M 3s
    ##  42400K .......... .......... .......... .......... .......... 31% 94.9M 3s
    ##  42450K .......... .......... .......... .......... .......... 31% 30.6M 3s
    ##  42500K .......... .......... .......... .......... .......... 31% 43.1M 3s
    ##  42550K .......... .......... .......... .......... .......... 31% 50.0M 3s
    ##  42600K .......... .......... .......... .......... .......... 31% 34.0M 3s
    ##  42650K .......... .......... .......... .......... .......... 31% 81.6M 3s
    ##  42700K .......... .......... .......... .......... .......... 31% 70.9M 3s
    ##  42750K .......... .......... .......... .......... .......... 31% 82.6M 3s
    ##  42800K .......... .......... .......... .......... .......... 31% 71.4M 3s
    ##  42850K .......... .......... .......... .......... .......... 31% 82.4M 3s
    ##  42900K .......... .......... .......... .......... .......... 31% 72.4M 3s
    ##  42950K .......... .......... .......... .......... .......... 31% 81.1M 3s
    ##  43000K .......... .......... .......... .......... .......... 31% 62.7M 3s
    ##  43050K .......... .......... .......... .......... .......... 31% 37.2M 3s
    ##  43100K .......... .......... .......... .......... .......... 32% 65.8M 3s
    ##  43150K .......... .......... .......... .......... .......... 32% 89.0M 3s
    ##  43200K .......... .......... .......... .......... .......... 32% 63.6M 3s
    ##  43250K .......... .......... .......... .......... .......... 32% 56.6M 3s
    ##  43300K .......... .......... .......... .......... .......... 32% 59.3M 3s
    ##  43350K .......... .......... .......... .......... .......... 32%  104M 3s
    ##  43400K .......... .......... .......... .......... .......... 32% 55.8M 3s
    ##  43450K .......... .......... .......... .......... .......... 32% 94.2M 3s
    ##  43500K .......... .......... .......... .......... .......... 32% 46.1M 3s
    ##  43550K .......... .......... .......... .......... .......... 32% 90.2M 3s
    ##  43600K .......... .......... .......... .......... .......... 32% 89.4M 3s
    ##  43650K .......... .......... .......... .......... .......... 32% 64.7M 3s
    ##  43700K .......... .......... .......... .......... .......... 32% 64.0M 3s
    ##  43750K .......... .......... .......... .......... .......... 32% 70.1M 3s
    ##  43800K .......... .......... .......... .......... .......... 32% 82.1M 3s
    ##  43850K .......... .......... .......... .......... .......... 32% 83.4M 3s
    ##  43900K .......... .......... .......... .......... .......... 32% 46.8M 3s
    ##  43950K .......... .......... .......... .......... .......... 32% 74.4M 3s
    ##  44000K .......... .......... .......... .......... .......... 32% 83.5M 3s
    ##  44050K .......... .......... .......... .......... .......... 32% 85.0M 3s
    ##  44100K .......... .......... .......... .......... .......... 32% 98.4M 3s
    ##  44150K .......... .......... .......... .......... .......... 32% 72.9M 3s
    ##  44200K .......... .......... .......... .......... .......... 32% 66.2M 3s
    ##  44250K .......... .......... .......... .......... .......... 32%  110M 3s
    ##  44300K .......... .......... .......... .......... .......... 32% 53.7M 3s
    ##  44350K .......... .......... .......... .......... .......... 32%  111M 3s
    ##  44400K .......... .......... .......... .......... .......... 32% 49.8M 3s
    ##  44450K .......... .......... .......... .......... .......... 33%  105M 3s
    ##  44500K .......... .......... .......... .......... .......... 33%  107M 3s
    ##  44550K .......... .......... .......... .......... .......... 33% 64.3M 3s
    ##  44600K .......... .......... .......... .......... .......... 33% 81.3M 3s
    ##  44650K .......... .......... .......... .......... .......... 33%  102M 3s
    ##  44700K .......... .......... .......... .......... .......... 33%  103M 3s
    ##  44750K .......... .......... .......... .......... .......... 33% 25.7M 3s
    ##  44800K .......... .......... .......... .......... .......... 33% 46.0M 3s
    ##  44850K .......... .......... .......... .......... .......... 33% 62.9M 3s
    ##  44900K .......... .......... .......... .......... .......... 33% 65.7M 3s
    ##  44950K .......... .......... .......... .......... .......... 33% 87.2M 3s
    ##  45000K .......... .......... .......... .......... .......... 33% 88.6M 3s
    ##  45050K .......... .......... .......... .......... .......... 33%  136M 3s
    ##  45100K .......... .......... .......... .......... .......... 33% 82.8M 3s
    ##  45150K .......... .......... .......... .......... .......... 33% 74.8M 3s
    ##  45200K .......... .......... .......... .......... .......... 33%  102M 3s
    ##  45250K .......... .......... .......... .......... .......... 33%  116M 3s
    ##  45300K .......... .......... .......... .......... .......... 33%  101M 3s
    ##  45350K .......... .......... .......... .......... .......... 33% 8.32M 3s
    ##  45400K .......... .......... .......... .......... .......... 33% 28.7M 3s
    ##  45450K .......... .......... .......... .......... .......... 33% 75.6M 3s
    ##  45500K .......... .......... .......... .......... .......... 33% 55.8M 3s
    ##  45550K .......... .......... .......... .......... .......... 33% 99.4M 3s
    ##  45600K .......... .......... .......... .......... .......... 33% 80.5M 3s
    ##  45650K .......... .......... .......... .......... .......... 33%  118M 3s
    ##  45700K .......... .......... .......... .......... .......... 33% 92.9M 3s
    ##  45750K .......... .......... .......... .......... .......... 33% 88.3M 3s
    ##  45800K .......... .......... .......... .......... .......... 34% 85.4M 3s
    ##  45850K .......... .......... .......... .......... .......... 34%  120M 3s
    ##  45900K .......... .......... .......... .......... .......... 34% 96.3M 3s
    ##  45950K .......... .......... .......... .......... .......... 34%  122M 3s
    ##  46000K .......... .......... .......... .......... .......... 34%  107M 3s
    ##  46050K .......... .......... .......... .......... .......... 34%  174K 3s
    ##  46100K .......... .......... .......... .......... .......... 34% 20.5M 3s
    ##  46150K .......... .......... .......... .......... .......... 34% 56.7M 3s
    ##  46200K .......... .......... .......... .......... .......... 34% 80.2M 3s
    ##  46250K .......... .......... .......... .......... .......... 34% 54.9M 3s
    ##  46300K .......... .......... .......... .......... .......... 34% 80.9M 3s
    ##  46350K .......... .......... .......... .......... .......... 34% 36.3M 3s
    ##  46400K .......... .......... .......... .......... .......... 34% 58.9M 3s
    ##  46450K .......... .......... .......... .......... .......... 34% 41.2M 3s
    ##  46500K .......... .......... .......... .......... .......... 34% 61.8M 3s
    ##  46550K .......... .......... .......... .......... .......... 34% 81.0M 3s
    ##  46600K .......... .......... .......... .......... .......... 34% 64.5M 3s
    ##  46650K .......... .......... .......... .......... .......... 34% 77.3M 3s
    ##  46700K .......... .......... .......... .......... .......... 34% 71.5M 3s
    ##  46750K .......... .......... .......... .......... .......... 34% 81.3M 3s
    ##  46800K .......... .......... .......... .......... .......... 34% 24.1M 3s
    ##  46850K .......... .......... .......... .......... .......... 34% 33.2M 3s
    ##  46900K .......... .......... .......... .......... .......... 34% 32.9M 3s
    ##  46950K .......... .......... .......... .......... .......... 34% 46.9M 3s
    ##  47000K .......... .......... .......... .......... .......... 34% 62.6M 3s
    ##  47050K .......... .......... .......... .......... .......... 34% 77.2M 3s
    ##  47100K .......... .......... .......... .......... .......... 34% 61.5M 3s
    ##  47150K .......... .......... .......... .......... .......... 35% 73.0M 3s
    ##  47200K .......... .......... .......... .......... .......... 35% 39.5M 3s
    ##  47250K .......... .......... .......... .......... .......... 35% 48.2M 3s
    ##  47300K .......... .......... .......... .......... .......... 35% 34.5M 3s
    ##  47350K .......... .......... .......... .......... .......... 35% 72.3M 3s
    ##  47400K .......... .......... .......... .......... .......... 35% 33.1M 3s
    ##  47450K .......... .......... .......... .......... .......... 35% 72.0M 3s
    ##  47500K .......... .......... .......... .......... .......... 35% 68.5M 3s
    ##  47550K .......... .......... .......... .......... .......... 35% 64.7M 3s
    ##  47600K .......... .......... .......... .......... .......... 35% 69.6M 3s
    ##  47650K .......... .......... .......... .......... .......... 35% 42.1M 3s
    ##  47700K .......... .......... .......... .......... .......... 35% 60.1M 3s
    ##  47750K .......... .......... .......... .......... .......... 35% 69.3M 3s
    ##  47800K .......... .......... .......... .......... .......... 35% 72.5M 3s
    ##  47850K .......... .......... .......... .......... .......... 35% 30.2M 3s
    ##  47900K .......... .......... .......... .......... .......... 35% 63.7M 3s
    ##  47950K .......... .......... .......... .......... .......... 35% 58.5M 3s
    ##  48000K .......... .......... .......... .......... .......... 35% 61.7M 3s
    ##  48050K .......... .......... .......... .......... .......... 35% 67.2M 3s
    ##  48100K .......... .......... .......... .......... .......... 35% 74.7M 3s
    ##  48150K .......... .......... .......... .......... .......... 35% 17.6M 3s
    ##  48200K .......... .......... .......... .......... .......... 35% 62.5M 3s
    ##  48250K .......... .......... .......... .......... .......... 35% 79.3M 3s
    ##  48300K .......... .......... .......... .......... .......... 35% 43.9M 3s
    ##  48350K .......... .......... .......... .......... .......... 35% 49.4M 3s
    ##  48400K .......... .......... .......... .......... .......... 35% 74.2M 3s
    ##  48450K .......... .......... .......... .......... .......... 35% 55.4M 3s
    ##  48500K .......... .......... .......... .......... .......... 36% 79.7M 3s
    ##  48550K .......... .......... .......... .......... .......... 36% 98.7M 3s
    ##  48600K .......... .......... .......... .......... .......... 36% 9.76M 3s
    ##  48650K .......... .......... .......... .......... .......... 36% 90.4M 3s
    ##  48700K .......... .......... .......... .......... .......... 36%  100M 3s
    ##  48750K .......... .......... .......... .......... .......... 36% 72.8M 3s
    ##  48800K .......... .......... .......... .......... .......... 36% 80.3M 3s
    ##  48850K .......... .......... .......... .......... .......... 36% 79.3M 3s
    ##  48900K .......... .......... .......... .......... .......... 36% 80.6M 3s
    ##  48950K .......... .......... .......... .......... .......... 36% 97.5M 3s
    ##  49000K .......... .......... .......... .......... .......... 36% 73.1M 3s
    ##  49050K .......... .......... .......... .......... .......... 36% 88.1M 3s
    ##  49100K .......... .......... .......... .......... .......... 36% 94.2M 3s
    ##  49150K .......... .......... .......... .......... .......... 36% 34.4M 3s
    ##  49200K .......... .......... .......... .......... .......... 36% 75.9M 3s
    ##  49250K .......... .......... .......... .......... .......... 36% 20.5M 3s
    ##  49300K .......... .......... .......... .......... .......... 36% 93.3M 3s
    ##  49350K .......... .......... .......... .......... .......... 36% 91.9M 3s
    ##  49400K .......... .......... .......... .......... .......... 36% 90.4M 3s
    ##  49450K .......... .......... .......... .......... .......... 36% 92.7M 3s
    ##  49500K .......... .......... .......... .......... .......... 36% 94.7M 3s
    ##  49550K .......... .......... .......... .......... .......... 36%  102M 3s
    ##  49600K .......... .......... .......... .......... .......... 36%  107M 3s
    ##  49650K .......... .......... .......... .......... .......... 36% 55.0M 3s
    ##  49700K .......... .......... .......... .......... .......... 36% 59.6M 3s
    ##  49750K .......... .......... .......... .......... .......... 36% 98.8M 3s
    ##  49800K .......... .......... .......... .......... .......... 36% 85.6M 3s
    ##  49850K .......... .......... .......... .......... .......... 37% 97.6M 3s
    ##  49900K .......... .......... .......... .......... .......... 37% 26.0M 3s
    ##  49950K .......... .......... .......... .......... .......... 37% 57.7M 3s
    ##  50000K .......... .......... .......... .......... .......... 37% 94.8M 3s
    ##  50050K .......... .......... .......... .......... .......... 37% 96.9M 3s
    ##  50100K .......... .......... .......... .......... .......... 37% 90.5M 3s
    ##  50150K .......... .......... .......... .......... .......... 37% 86.7M 3s
    ##  50200K .......... .......... .......... .......... .......... 37%  109M 3s
    ##  50250K .......... .......... .......... .......... .......... 37% 88.2M 3s
    ##  50300K .......... .......... .......... .......... .......... 37% 76.3M 3s
    ##  50350K .......... .......... .......... .......... .......... 37%  111M 3s
    ##  50400K .......... .......... .......... .......... .......... 37% 62.2M 3s
    ##  50450K .......... .......... .......... .......... .......... 37% 99.9M 3s
    ##  50500K .......... .......... .......... .......... .......... 37% 95.7M 3s
    ##  50550K .......... .......... .......... .......... .......... 37% 24.9M 3s
    ##  50600K .......... .......... .......... .......... .......... 37% 15.6M 3s
    ##  50650K .......... .......... .......... .......... .......... 37% 76.1M 3s
    ##  50700K .......... .......... .......... .......... .......... 37% 98.1M 3s
    ##  50750K .......... .......... .......... .......... .......... 37%  114M 3s
    ##  50800K .......... .......... .......... .......... .......... 37%  117M 3s
    ##  50850K .......... .......... .......... .......... .......... 37%  103M 3s
    ##  50900K .......... .......... .......... .......... .......... 37%  113M 3s
    ##  50950K .......... .......... .......... .......... .......... 37%  126M 3s
    ##  51000K .......... .......... .......... .......... .......... 37% 91.1M 3s
    ##  51050K .......... .......... .......... .......... .......... 37%  120M 3s
    ##  51100K .......... .......... .......... .......... .......... 37%  110M 3s
    ##  51150K .......... .......... .......... .......... .......... 37%  108M 3s
    ##  51200K .......... .......... .......... .......... .......... 38% 61.9K 4s
    ##  51250K .......... .......... .......... .......... .......... 38% 12.1M 4s
    ##  51300K .......... .......... .......... .......... .......... 38% 74.1M 4s
    ##  51350K .......... .......... .......... .......... .......... 38% 15.8M 4s
    ##  51400K .......... .......... .......... .......... .......... 38% 47.5M 4s
    ##  51450K .......... .......... .......... .......... .......... 38% 40.4M 4s
    ##  51500K .......... .......... .......... .......... .......... 38% 91.7M 4s
    ##  51550K .......... .......... .......... .......... .......... 38% 10.1M 4s
    ##  51600K .......... .......... .......... .......... .......... 38% 63.5M 4s
    ##  51650K .......... .......... .......... .......... .......... 38% 65.2M 4s
    ##  51700K .......... .......... .......... .......... .......... 38% 59.3M 4s
    ##  51750K .......... .......... .......... .......... .......... 38% 62.2M 4s
    ##  51800K .......... .......... .......... .......... .......... 38% 72.5M 4s
    ##  51850K .......... .......... .......... .......... .......... 38% 96.2M 4s
    ##  51900K .......... .......... .......... .......... .......... 38% 83.6M 4s
    ##  51950K .......... .......... .......... .......... .......... 38% 84.1M 4s
    ##  52000K .......... .......... .......... .......... .......... 38% 43.9M 4s
    ##  52050K .......... .......... .......... .......... .......... 38% 90.9M 4s
    ##  52100K .......... .......... .......... .......... .......... 38% 97.5M 4s
    ##  52150K .......... .......... .......... .......... .......... 38%  105M 4s
    ##  52200K .......... .......... .......... .......... .......... 38% 30.8M 4s
    ##  52250K .......... .......... .......... .......... .......... 38% 82.3M 4s
    ##  52300K .......... .......... .......... .......... .......... 38% 77.4M 4s
    ##  52350K .......... .......... .......... .......... .......... 38% 98.2M 4s
    ##  52400K .......... .......... .......... .......... .......... 38% 92.4M 4s
    ##  52450K .......... .......... .......... .......... .......... 38% 93.3M 4s
    ##  52500K .......... .......... .......... .......... .......... 39% 73.7M 4s
    ##  52550K .......... .......... .......... .......... .......... 39% 64.6M 4s
    ##  52600K .......... .......... .......... .......... .......... 39% 81.7M 4s
    ##  52650K .......... .......... .......... .......... .......... 39%  108M 4s
    ##  52700K .......... .......... .......... .......... .......... 39% 75.3M 4s
    ##  52750K .......... .......... .......... .......... .......... 39% 58.4M 4s
    ##  52800K .......... .......... .......... .......... .......... 39% 20.7M 4s
    ##  52850K .......... .......... .......... .......... .......... 39% 79.1M 4s
    ##  52900K .......... .......... .......... .......... .......... 39% 99.8M 4s
    ##  52950K .......... .......... .......... .......... .......... 39% 98.6M 4s
    ##  53000K .......... .......... .......... .......... .......... 39% 88.5M 4s
    ##  53050K .......... .......... .......... .......... .......... 39% 30.4M 4s
    ##  53100K .......... .......... .......... .......... .......... 39% 70.8M 4s
    ##  53150K .......... .......... .......... .......... .......... 39% 73.1M 4s
    ##  53200K .......... .......... .......... .......... .......... 39% 55.1M 4s
    ##  53250K .......... .......... .......... .......... .......... 39% 71.1M 4s
    ##  53300K .......... .......... .......... .......... .......... 39% 96.5M 4s
    ##  53350K .......... .......... .......... .......... .......... 39%  109M 4s
    ##  53400K .......... .......... .......... .......... .......... 39% 96.5M 4s
    ##  53450K .......... .......... .......... .......... .......... 39% 15.6M 4s
    ##  53500K .......... .......... .......... .......... .......... 39% 71.7M 4s
    ##  53550K .......... .......... .......... .......... .......... 39% 22.1M 4s
    ##  53600K .......... .......... .......... .......... .......... 39% 96.3M 4s
    ##  53650K .......... .......... .......... .......... .......... 39% 49.6M 4s
    ##  53700K .......... .......... .......... .......... .......... 39% 96.8M 4s
    ##  53750K .......... .......... .......... .......... .......... 39% 95.4M 4s
    ##  53800K .......... .......... .......... .......... .......... 39% 94.0M 4s
    ##  53850K .......... .......... .......... .......... .......... 40% 90.3M 4s
    ##  53900K .......... .......... .......... .......... .......... 40% 74.7M 4s
    ##  53950K .......... .......... .......... .......... .......... 40% 27.4M 4s
    ##  54000K .......... .......... .......... .......... .......... 40% 49.0M 4s
    ##  54050K .......... .......... .......... .......... .......... 40% 41.3M 4s
    ##  54100K .......... .......... .......... .......... .......... 40% 85.3M 4s
    ##  54150K .......... .......... .......... .......... .......... 40% 41.6M 4s
    ##  54200K .......... .......... .......... .......... .......... 40% 81.9M 4s
    ##  54250K .......... .......... .......... .......... .......... 40%  102M 4s
    ##  54300K .......... .......... .......... .......... .......... 40% 98.2M 4s
    ##  54350K .......... .......... .......... .......... .......... 40% 21.6M 4s
    ##  54400K .......... .......... .......... .......... .......... 40% 99.7M 4s
    ##  54450K .......... .......... .......... .......... .......... 40%  111M 4s
    ##  54500K .......... .......... .......... .......... .......... 40% 99.2M 4s
    ##  54550K .......... .......... .......... .......... .......... 40%  111M 4s
    ##  54600K .......... .......... .......... .......... .......... 40% 94.2M 4s
    ##  54650K .......... .......... .......... .......... .......... 40%  103M 4s
    ##  54700K .......... .......... .......... .......... .......... 40% 26.0M 4s
    ##  54750K .......... .......... .......... .......... .......... 40% 83.8M 4s
    ##  54800K .......... .......... .......... .......... .......... 40% 37.3M 4s
    ##  54850K .......... .......... .......... .......... .......... 40% 61.1M 4s
    ##  54900K .......... .......... .......... .......... .......... 40% 64.3M 4s
    ##  54950K .......... .......... .......... .......... .......... 40% 46.4M 4s
    ##  55000K .......... .......... .......... .......... .......... 40% 29.6M 4s
    ##  55050K .......... .......... .......... .......... .......... 40% 88.6M 4s
    ##  55100K .......... .......... .......... .......... .......... 40%  108M 4s
    ##  55150K .......... .......... .......... .......... .......... 40% 96.5M 4s
    ##  55200K .......... .......... .......... .......... .......... 41% 96.2M 4s
    ##  55250K .......... .......... .......... .......... .......... 41% 74.7M 4s
    ##  55300K .......... .......... .......... .......... .......... 41% 68.4M 4s
    ##  55350K .......... .......... .......... .......... .......... 41% 85.9M 4s
    ##  55400K .......... .......... .......... .......... .......... 41% 74.8M 4s
    ##  55450K .......... .......... .......... .......... .......... 41%  107M 4s
    ##  55500K .......... .......... .......... .......... .......... 41%  104M 4s
    ##  55550K .......... .......... .......... .......... .......... 41% 91.5M 4s
    ##  55600K .......... .......... .......... .......... .......... 41% 88.9M 4s
    ##  55650K .......... .......... .......... .......... .......... 41% 96.8M 4s
    ##  55700K .......... .......... .......... .......... .......... 41% 94.2M 4s
    ##  55750K .......... .......... .......... .......... .......... 41% 86.3M 4s
    ##  55800K .......... .......... .......... .......... .......... 41%  106M 4s
    ##  55850K .......... .......... .......... .......... .......... 41% 20.8M 4s
    ##  55900K .......... .......... .......... .......... .......... 41% 51.6M 4s
    ##  55950K .......... .......... .......... .......... .......... 41%  111M 4s
    ##  56000K .......... .......... .......... .......... .......... 41%  117M 4s
    ##  56050K .......... .......... .......... .......... .......... 41%  118M 4s
    ##  56100K .......... .......... .......... .......... .......... 41%  111M 4s
    ##  56150K .......... .......... .......... .......... .......... 41% 92.3M 4s
    ##  56200K .......... .......... .......... .......... .......... 41%  113M 4s
    ##  56250K .......... .......... .......... .......... .......... 41%  118M 4s
    ##  56300K .......... .......... .......... .......... .......... 41% 80.3K 5s
    ##  56350K .......... .......... .......... .......... .......... 41% 11.7M 5s
    ##  56400K .......... .......... .......... .......... .......... 41% 88.0M 5s
    ##  56450K .......... .......... .......... .......... .......... 41% 14.7M 5s
    ##  56500K .......... .......... .......... .......... .......... 41% 83.9M 5s
    ##  56550K .......... .......... .......... .......... .......... 42% 50.6M 5s
    ##  56600K .......... .......... .......... .......... .......... 42% 73.5M 5s
    ##  56650K .......... .......... .......... .......... .......... 42% 94.3M 5s
    ##  56700K .......... .......... .......... .......... .......... 42% 37.1M 5s
    ##  56750K .......... .......... .......... .......... .......... 42% 94.2M 5s
    ##  56800K .......... .......... .......... .......... .......... 42% 34.2M 5s
    ##  56850K .......... .......... .......... .......... .......... 42% 47.3M 5s
    ##  56900K .......... .......... .......... .......... .......... 42% 85.9M 5s
    ##  56950K .......... .......... .......... .......... .......... 42% 81.7M 5s
    ##  57000K .......... .......... .......... .......... .......... 42% 92.7M 5s
    ##  57050K .......... .......... .......... .......... .......... 42%  101M 5s
    ##  57100K .......... .......... .......... .......... .......... 42% 87.5M 5s
    ##  57150K .......... .......... .......... .......... .......... 42% 69.9M 5s
    ##  57200K .......... .......... .......... .......... .......... 42% 43.5M 5s
    ##  57250K .......... .......... .......... .......... .......... 42% 47.7M 5s
    ##  57300K .......... .......... .......... .......... .......... 42% 66.5M 5s
    ##  57350K .......... .......... .......... .......... .......... 42% 53.4M 5s
    ##  57400K .......... .......... .......... .......... .......... 42% 31.1M 5s
    ##  57450K .......... .......... .......... .......... .......... 42% 79.7M 5s
    ##  57500K .......... .......... .......... .......... .......... 42% 9.45M 5s
    ##  57550K .......... .......... .......... .......... .......... 42% 90.0M 5s
    ##  57600K .......... .......... .......... .......... .......... 42% 77.8M 5s
    ##  57650K .......... .......... .......... .......... .......... 42% 91.1M 5s
    ##  57700K .......... .......... .......... .......... .......... 42% 89.4M 5s
    ##  57750K .......... .......... .......... .......... .......... 42% 94.4M 5s
    ##  57800K .......... .......... .......... .......... .......... 42%  105M 5s
    ##  57850K .......... .......... .......... .......... .......... 42%  105M 5s
    ##  57900K .......... .......... .......... .......... .......... 43% 15.3M 5s
    ##  57950K .......... .......... .......... .......... .......... 43% 83.0M 5s
    ##  58000K .......... .......... .......... .......... .......... 43% 72.7M 5s
    ##  58050K .......... .......... .......... .......... .......... 43% 88.9M 5s
    ##  58100K .......... .......... .......... .......... .......... 43% 95.2M 5s
    ##  58150K .......... .......... .......... .......... .......... 43%  112M 5s
    ##  58200K .......... .......... .......... .......... .......... 43% 98.2M 5s
    ##  58250K .......... .......... .......... .......... .......... 43%  107M 5s
    ##  58300K .......... .......... .......... .......... .......... 43% 90.5M 4s
    ##  58350K .......... .......... .......... .......... .......... 43% 39.1M 4s
    ##  58400K .......... .......... .......... .......... .......... 43% 78.4M 4s
    ##  58450K .......... .......... .......... .......... .......... 43% 80.5M 4s
    ##  58500K .......... .......... .......... .......... .......... 43% 81.8M 4s
    ##  58550K .......... .......... .......... .......... .......... 43% 85.0M 4s
    ##  58600K .......... .......... .......... .......... .......... 43% 86.9M 4s
    ##  58650K .......... .......... .......... .......... .......... 43% 79.3M 4s
    ##  58700K .......... .......... .......... .......... .......... 43% 93.7M 4s
    ##  58750K .......... .......... .......... .......... .......... 43% 96.5M 4s
    ##  58800K .......... .......... .......... .......... .......... 43%  106M 4s
    ##  58850K .......... .......... .......... .......... .......... 43% 90.4M 4s
    ##  58900K .......... .......... .......... .......... .......... 43% 92.2M 4s
    ##  58950K .......... .......... .......... .......... .......... 43% 86.9M 4s
    ##  59000K .......... .......... .......... .......... .......... 43% 84.1M 4s
    ##  59050K .......... .......... .......... .......... .......... 43% 80.3M 4s
    ##  59100K .......... .......... .......... .......... .......... 43% 85.2M 4s
    ##  59150K .......... .......... .......... .......... .......... 43%  103M 4s
    ##  59200K .......... .......... .......... .......... .......... 43% 29.3M 4s
    ##  59250K .......... .......... .......... .......... .......... 44%  103M 4s
    ##  59300K .......... .......... .......... .......... .......... 44% 84.6M 4s
    ##  59350K .......... .......... .......... .......... .......... 44%  109M 4s
    ##  59400K .......... .......... .......... .......... .......... 44% 82.8M 4s
    ##  59450K .......... .......... .......... .......... .......... 44% 99.4M 4s
    ##  59500K .......... .......... .......... .......... .......... 44% 83.4M 4s
    ##  59550K .......... .......... .......... .......... .......... 44% 84.8M 4s
    ##  59600K .......... .......... .......... .......... .......... 44%  112M 4s
    ##  59650K .......... .......... .......... .......... .......... 44% 86.1M 4s
    ##  59700K .......... .......... .......... .......... .......... 44% 75.3M 4s
    ##  59750K .......... .......... .......... .......... .......... 44%  111M 4s
    ##  59800K .......... .......... .......... .......... .......... 44% 74.6M 4s
    ##  59850K .......... .......... .......... .......... .......... 44% 62.4M 4s
    ##  59900K .......... .......... .......... .......... .......... 44% 11.8M 4s
    ##  59950K .......... .......... .......... .......... .......... 44%  111M 4s
    ##  60000K .......... .......... .......... .......... .......... 44%  112M 4s
    ##  60050K .......... .......... .......... .......... .......... 44% 93.3M 4s
    ##  60100K .......... .......... .......... .......... .......... 44% 96.4M 4s
    ##  60150K .......... .......... .......... .......... .......... 44%  110M 4s
    ##  60200K .......... .......... .......... .......... .......... 44%  105M 4s
    ##  60250K .......... .......... .......... .......... .......... 44%  126M 4s
    ##  60300K .......... .......... .......... .......... .......... 44% 95.1M 4s
    ##  60350K .......... .......... .......... .......... .......... 44%  103M 4s
    ##  60400K .......... .......... .......... .......... .......... 44%  103M 4s
    ##  60450K .......... .......... .......... .......... .......... 44%  104M 4s
    ##  60500K .......... .......... .......... .......... .......... 44% 82.2M 4s
    ##  60550K .......... .......... .......... .......... .......... 44% 63.1M 4s
    ##  60600K .......... .......... .......... .......... .......... 45% 78.8M 4s
    ##  60650K .......... .......... .......... .......... .......... 45% 88.7M 4s
    ##  60700K .......... .......... .......... .......... .......... 45% 83.9M 4s
    ##  60750K .......... .......... .......... .......... .......... 45%  108M 4s
    ##  60800K .......... .......... .......... .......... .......... 45% 22.3M 4s
    ##  60850K .......... .......... .......... .......... .......... 45% 98.6M 4s
    ##  60900K .......... .......... .......... .......... .......... 45% 48.1M 4s
    ##  60950K .......... .......... .......... .......... .......... 45%  105M 4s
    ##  61000K .......... .......... .......... .......... .......... 45% 68.6M 4s
    ##  61050K .......... .......... .......... .......... .......... 45% 78.6M 4s
    ##  61100K .......... .......... .......... .......... .......... 45% 72.1M 4s
    ##  61150K .......... .......... .......... .......... .......... 45% 68.1M 4s
    ##  61200K .......... .......... .......... .......... .......... 45% 94.4M 4s
    ##  61250K .......... .......... .......... .......... .......... 45%  106M 4s
    ##  61300K .......... .......... .......... .......... .......... 45% 91.4M 4s
    ##  61350K .......... .......... .......... .......... .......... 45% 88.9M 4s
    ##  61400K .......... .......... .......... .......... .......... 45% 52.1K 5s
    ##  61450K .......... .......... .......... .......... .......... 45% 11.4M 5s
    ##  61500K .......... .......... .......... .......... .......... 45% 11.6M 5s
    ##  61550K .......... .......... .......... .......... .......... 45% 83.7M 5s
    ##  61600K .......... .......... .......... .......... .......... 45% 14.1M 5s
    ##  61650K .......... .......... .......... .......... .......... 45%  108M 5s
    ##  61700K .......... .......... .......... .......... .......... 45% 35.7M 5s
    ##  61750K .......... .......... .......... .......... .......... 45% 78.8M 5s
    ##  61800K .......... .......... .......... .......... .......... 45% 29.6M 5s
    ##  61850K .......... .......... .......... .......... .......... 45% 63.9M 5s
    ##  61900K .......... .......... .......... .......... .......... 45% 80.1M 5s
    ##  61950K .......... .......... .......... .......... .......... 46% 57.9M 5s
    ##  62000K .......... .......... .......... .......... .......... 46% 36.9M 5s
    ##  62050K .......... .......... .......... .......... .......... 46% 59.4M 5s
    ##  62100K .......... .......... .......... .......... .......... 46% 69.5M 5s
    ##  62150K .......... .......... .......... .......... .......... 46%  100M 5s
    ##  62200K .......... .......... .......... .......... .......... 46% 4.73M 5s
    ##  62250K .......... .......... .......... .......... .......... 46% 42.8M 5s
    ##  62300K .......... .......... .......... .......... .......... 46% 28.2M 5s
    ##  62350K .......... .......... .......... .......... .......... 46% 65.8M 5s
    ##  62400K .......... .......... .......... .......... .......... 46% 82.1M 5s
    ##  62450K .......... .......... .......... .......... .......... 46% 92.3M 5s
    ##  62500K .......... .......... .......... .......... .......... 46% 75.4M 5s
    ##  62550K .......... .......... .......... .......... .......... 46%  112M 5s
    ##  62600K .......... .......... .......... .......... .......... 46% 91.9M 5s
    ##  62650K .......... .......... .......... .......... .......... 46% 66.6M 5s
    ##  62700K .......... .......... .......... .......... .......... 46% 60.5M 5s
    ##  62750K .......... .......... .......... .......... .......... 46% 91.3M 5s
    ##  62800K .......... .......... .......... .......... .......... 46% 24.9M 5s
    ##  62850K .......... .......... .......... .......... .......... 46%  121M 5s
    ##  62900K .......... .......... .......... .......... .......... 46% 41.5M 5s
    ##  62950K .......... .......... .......... .......... .......... 46% 86.8M 5s
    ##  63000K .......... .......... .......... .......... .......... 46% 94.6M 5s
    ##  63050K .......... .......... .......... .......... .......... 46%  108M 5s
    ##  63100K .......... .......... .......... .......... .......... 46% 83.9M 5s
    ##  63150K .......... .......... .......... .......... .......... 46% 87.3M 5s
    ##  63200K .......... .......... .......... .......... .......... 46% 87.4M 5s
    ##  63250K .......... .......... .......... .......... .......... 46% 76.3M 5s
    ##  63300K .......... .......... .......... .......... .......... 47% 66.9M 5s
    ##  63350K .......... .......... .......... .......... .......... 47% 34.9M 5s
    ##  63400K .......... .......... .......... .......... .......... 47% 63.3M 5s
    ##  63450K .......... .......... .......... .......... .......... 47%  118M 5s
    ##  63500K .......... .......... .......... .......... .......... 47% 95.0M 5s
    ##  63550K .......... .......... .......... .......... .......... 47% 62.5M 5s
    ##  63600K .......... .......... .......... .......... .......... 47% 76.2M 5s
    ##  63650K .......... .......... .......... .......... .......... 47%  108M 5s
    ##  63700K .......... .......... .......... .......... .......... 47% 81.3M 5s
    ##  63750K .......... .......... .......... .......... .......... 47%  108M 5s
    ##  63800K .......... .......... .......... .......... .......... 47% 75.3M 5s
    ##  63850K .......... .......... .......... .......... .......... 47% 72.2M 5s
    ##  63900K .......... .......... .......... .......... .......... 47% 66.9M 5s
    ##  63950K .......... .......... .......... .......... .......... 47%  108M 5s
    ##  64000K .......... .......... .......... .......... .......... 47% 15.4M 5s
    ##  64050K .......... .......... .......... .......... .......... 47% 88.1M 5s
    ##  64100K .......... .......... .......... .......... .......... 47% 95.6M 5s
    ##  64150K .......... .......... .......... .......... .......... 47% 99.4M 5s
    ##  64200K .......... .......... .......... .......... .......... 47% 99.2M 5s
    ##  64250K .......... .......... .......... .......... .......... 47%  114M 5s
    ##  64300K .......... .......... .......... .......... .......... 47% 96.9M 5s
    ##  64350K .......... .......... .......... .......... .......... 47% 89.7M 5s
    ##  64400K .......... .......... .......... .......... .......... 47% 89.3M 5s
    ##  64450K .......... .......... .......... .......... .......... 47% 57.2M 5s
    ##  64500K .......... .......... .......... .......... .......... 47% 96.0M 5s
    ##  64550K .......... .......... .......... .......... .......... 47%  115M 5s
    ##  64600K .......... .......... .......... .......... .......... 47% 12.1M 5s
    ##  64650K .......... .......... .......... .......... .......... 48%  113M 5s
    ##  64700K .......... .......... .......... .......... .......... 48% 95.9M 5s
    ##  64750K .......... .......... .......... .......... .......... 48% 99.7M 5s
    ##  64800K .......... .......... .......... .......... .......... 48% 80.8M 5s
    ##  64850K .......... .......... .......... .......... .......... 48%  100M 5s
    ##  64900K .......... .......... .......... .......... .......... 48%  108M 5s
    ##  64950K .......... .......... .......... .......... .......... 48%  127M 5s
    ##  65000K .......... .......... .......... .......... .......... 48% 19.3M 5s
    ##  65050K .......... .......... .......... .......... .......... 48% 70.7M 5s
    ##  65100K .......... .......... .......... .......... .......... 48% 99.2M 5s
    ##  65150K .......... .......... .......... .......... .......... 48% 78.5M 5s
    ##  65200K .......... .......... .......... .......... .......... 48% 81.9M 5s
    ##  65250K .......... .......... .......... .......... .......... 48% 92.3M 5s
    ##  65300K .......... .......... .......... .......... .......... 48% 85.6M 5s
    ##  65350K .......... .......... .......... .......... .......... 48% 94.7M 5s
    ##  65400K .......... .......... .......... .......... .......... 48% 96.3M 5s
    ##  65450K .......... .......... .......... .......... .......... 48% 99.4M 5s
    ##  65500K .......... .......... .......... .......... .......... 48% 52.1M 5s
    ##  65550K .......... .......... .......... .......... .......... 48% 41.4M 5s
    ##  65600K .......... .......... .......... .......... .......... 48% 67.8M 5s
    ##  65650K .......... .......... .......... .......... .......... 48% 91.7M 5s
    ##  65700K .......... .......... .......... .......... .......... 48% 65.0M 5s
    ##  65750K .......... .......... .......... .......... .......... 48% 97.5M 5s
    ##  65800K .......... .......... .......... .......... .......... 48% 80.3M 5s
    ##  65850K .......... .......... .......... .......... .......... 48% 90.7M 5s
    ##  65900K .......... .......... .......... .......... .......... 48% 64.5M 5s
    ##  65950K .......... .......... .......... .......... .......... 48% 26.5M 5s
    ##  66000K .......... .......... .......... .......... .......... 49% 72.0M 5s
    ##  66050K .......... .......... .......... .......... .......... 49% 77.5M 5s
    ##  66100K .......... .......... .......... .......... .......... 49% 83.8M 5s
    ##  66150K .......... .......... .......... .......... .......... 49% 85.0M 5s
    ##  66200K .......... .......... .......... .......... .......... 49% 79.3M 5s
    ##  66250K .......... .......... .......... .......... .......... 49% 87.9M 5s
    ##  66300K .......... .......... .......... .......... .......... 49% 80.9M 5s
    ##  66350K .......... .......... .......... .......... .......... 49%  105M 5s
    ##  66400K .......... .......... .......... .......... .......... 49% 53.9M 5s
    ##  66450K .......... .......... .......... .......... .......... 49% 90.0M 5s
    ##  66500K .......... .......... .......... .......... .......... 49% 84.9M 5s
    ##  66550K .......... .......... .......... .......... .......... 49% 1.36M 5s
    ##  66600K .......... .......... .......... .......... .......... 49% 84.4M 5s
    ##  66650K .......... .......... .......... .......... .......... 49% 29.9M 5s
    ##  66700K .......... .......... .......... .......... .......... 49% 19.0M 5s
    ##  66750K .......... .......... .......... .......... .......... 49% 91.5M 5s
    ##  66800K .......... .......... .......... .......... .......... 49% 86.0M 5s
    ##  66850K .......... .......... .......... .......... .......... 49% 24.6M 5s
    ##  66900K .......... .......... .......... .......... .......... 49% 64.0M 5s
    ##  66950K .......... .......... .......... .......... .......... 49% 95.0M 5s
    ##  67000K .......... .......... .......... .......... .......... 49% 72.9M 5s
    ##  67050K .......... .......... .......... .......... .......... 49% 99.5M 5s
    ##  67100K .......... .......... .......... .......... .......... 49% 90.0M 5s
    ##  67150K .......... .......... .......... .......... .......... 49% 71.9M 5s
    ##  67200K .......... .......... .......... .......... .......... 49% 66.2M 5s
    ##  67250K .......... .......... .......... .......... .......... 49% 71.5M 5s
    ##  67300K .......... .......... .......... .......... .......... 49% 85.5M 5s
    ##  67350K .......... .......... .......... .......... .......... 50%  105M 5s
    ##  67400K .......... .......... .......... .......... .......... 50% 69.7M 5s
    ##  67450K .......... .......... .......... .......... .......... 50% 87.0M 5s
    ##  67500K .......... .......... .......... .......... .......... 50% 62.1M 5s
    ##  67550K .......... .......... .......... .......... .......... 50% 62.4M 5s
    ##  67600K .......... .......... .......... .......... .......... 50% 73.2M 5s
    ##  67650K .......... .......... .......... .......... .......... 50% 86.8M 5s
    ##  67700K .......... .......... .......... .......... .......... 50% 74.3M 5s
    ##  67750K .......... .......... .......... .......... .......... 50% 94.3M 5s
    ##  67800K .......... .......... .......... .......... .......... 50% 76.9M 5s
    ##  67850K .......... .......... .......... .......... .......... 50% 92.1M 5s
    ##  67900K .......... .......... .......... .......... .......... 50% 81.1M 5s
    ##  67950K .......... .......... .......... .......... .......... 50% 27.6M 5s
    ##  68000K .......... .......... .......... .......... .......... 50% 50.7M 4s
    ##  68050K .......... .......... .......... .......... .......... 50% 72.7M 4s
    ##  68100K .......... .......... .......... .......... .......... 50%  106M 4s
    ##  68150K .......... .......... .......... .......... .......... 50% 91.9M 4s
    ##  68200K .......... .......... .......... .......... .......... 50% 92.2M 4s
    ##  68250K .......... .......... .......... .......... .......... 50% 36.2M 4s
    ##  68300K .......... .......... .......... .......... .......... 50% 20.0M 4s
    ##  68350K .......... .......... .......... .......... .......... 50% 21.1M 4s
    ##  68400K .......... .......... .......... .......... .......... 50% 23.9M 4s
    ##  68450K .......... .......... .......... .......... .......... 50% 22.7M 4s
    ##  68500K .......... .......... .......... .......... .......... 50% 23.3M 4s
    ##  68550K .......... .......... .......... .......... .......... 50% 23.2M 4s
    ##  68600K .......... .......... .......... .......... .......... 50% 22.8M 4s
    ##  68650K .......... .......... .......... .......... .......... 50% 21.5M 4s
    ##  68700K .......... .......... .......... .......... .......... 51% 24.8M 4s
    ##  68750K .......... .......... .......... .......... .......... 51% 21.8M 4s
    ##  68800K .......... .......... .......... .......... .......... 51% 27.2M 4s
    ##  68850K .......... .......... .......... .......... .......... 51% 28.7M 4s
    ##  68900K .......... .......... .......... .......... .......... 51% 26.7M 4s
    ##  68950K .......... .......... .......... .......... .......... 51% 23.5M 4s
    ##  69000K .......... .......... .......... .......... .......... 51% 24.9M 4s
    ##  69050K .......... .......... .......... .......... .......... 51% 30.8M 4s
    ##  69100K .......... .......... .......... .......... .......... 51% 85.4M 4s
    ##  69150K .......... .......... .......... .......... .......... 51%  104M 4s
    ##  69200K .......... .......... .......... .......... .......... 51% 93.6M 4s
    ##  69250K .......... .......... .......... .......... .......... 51% 90.2M 4s
    ##  69300K .......... .......... .......... .......... .......... 51% 91.2M 4s
    ##  69350K .......... .......... .......... .......... .......... 51%  104M 4s
    ##  69400K .......... .......... .......... .......... .......... 51% 89.1M 4s
    ##  69450K .......... .......... .......... .......... .......... 51%  101M 4s
    ##  69500K .......... .......... .......... .......... .......... 51% 75.3M 4s
    ##  69550K .......... .......... .......... .......... .......... 51% 88.2M 4s
    ##  69600K .......... .......... .......... .......... .......... 51% 93.3M 4s
    ##  69650K .......... .......... .......... .......... .......... 51% 89.4M 4s
    ##  69700K .......... .......... .......... .......... .......... 51% 82.7M 4s
    ##  69750K .......... .......... .......... .......... .......... 51% 98.6M 4s
    ##  69800K .......... .......... .......... .......... .......... 51% 66.7M 4s
    ##  69850K .......... .......... .......... .......... .......... 51% 57.5M 4s
    ##  69900K .......... .......... .......... .......... .......... 51% 79.8M 4s
    ##  69950K .......... .......... .......... .......... .......... 51% 98.7M 4s
    ##  70000K .......... .......... .......... .......... .......... 51% 89.4M 4s
    ##  70050K .......... .......... .......... .......... .......... 52%  101M 4s
    ##  70100K .......... .......... .......... .......... .......... 52% 63.1M 4s
    ##  70150K .......... .......... .......... .......... .......... 52%  101M 4s
    ##  70200K .......... .......... .......... .......... .......... 52% 96.6M 4s
    ##  70250K .......... .......... .......... .......... .......... 52% 88.4M 4s
    ##  70300K .......... .......... .......... .......... .......... 52%  101M 4s
    ##  70350K .......... .......... .......... .......... .......... 52% 33.3M 4s
    ##  70400K .......... .......... .......... .......... .......... 52% 83.6M 4s
    ##  70450K .......... .......... .......... .......... .......... 52%  103M 4s
    ##  70500K .......... .......... .......... .......... .......... 52% 95.4M 4s
    ##  70550K .......... .......... .......... .......... .......... 52% 82.4M 4s
    ##  70600K .......... .......... .......... .......... .......... 52% 94.2M 4s
    ##  70650K .......... .......... .......... .......... .......... 52% 70.3M 4s
    ##  70700K .......... .......... .......... .......... .......... 52% 94.2M 4s
    ##  70750K .......... .......... .......... .......... .......... 52% 92.8M 4s
    ##  70800K .......... .......... .......... .......... .......... 52% 93.2M 4s
    ##  70850K .......... .......... .......... .......... .......... 52%  115M 4s
    ##  70900K .......... .......... .......... .......... .......... 52% 99.3M 4s
    ##  70950K .......... .......... .......... .......... .......... 52%  117M 4s
    ##  71000K .......... .......... .......... .......... .......... 52% 65.0M 4s
    ##  71050K .......... .......... .......... .......... .......... 52% 79.2M 4s
    ##  71100K .......... .......... .......... .......... .......... 52% 71.3M 4s
    ##  71150K .......... .......... .......... .......... .......... 52%  113M 4s
    ##  71200K .......... .......... .......... .......... .......... 52% 95.9M 4s
    ##  71250K .......... .......... .......... .......... .......... 52%  121M 4s
    ##  71300K .......... .......... .......... .......... .......... 52% 94.9M 4s
    ##  71350K .......... .......... .......... .......... .......... 52%  105M 4s
    ##  71400K .......... .......... .......... .......... .......... 53% 72.5M 4s
    ##  71450K .......... .......... .......... .......... .......... 53%  113M 4s
    ##  71500K .......... .......... .......... .......... .......... 53% 82.4M 4s
    ##  71550K .......... .......... .......... .......... .......... 53%  107M 4s
    ##  71600K .......... .......... .......... .......... .......... 53% 84.8M 4s
    ##  71650K .......... .......... .......... .......... .......... 53%  107M 4s
    ##  71700K .......... .......... .......... .......... .......... 53% 97.2M 4s
    ##  71750K .......... .......... .......... .......... .......... 53% 67.7M 4s
    ##  71800K .......... .......... .......... .......... .......... 53% 98.9M 4s
    ##  71850K .......... .......... .......... .......... .......... 53%  125M 4s
    ##  71900K .......... .......... .......... .......... .......... 53% 83.8M 4s
    ##  71950K .......... .......... .......... .......... .......... 53% 94.4M 4s
    ##  72000K .......... .......... .......... .......... .......... 53%  104M 4s
    ##  72050K .......... .......... .......... .......... .......... 53%  112M 4s
    ##  72100K .......... .......... .......... .......... .......... 53%  103M 4s
    ##  72150K .......... .......... .......... .......... .......... 53%  129M 4s
    ##  72200K .......... .......... .......... .......... .......... 53% 99.2M 4s
    ##  72250K .......... .......... .......... .......... .......... 53%  113M 4s
    ##  72300K .......... .......... .......... .......... .......... 53% 88.9M 4s
    ##  72350K .......... .......... .......... .......... .......... 53%  105M 4s
    ##  72400K .......... .......... .......... .......... .......... 53%  102M 4s
    ##  72450K .......... .......... .......... .......... .......... 53%  107M 4s
    ##  72500K .......... .......... .......... .......... .......... 53% 98.8M 4s
    ##  72550K .......... .......... .......... .......... .......... 53% 73.4M 4s
    ##  72600K .......... .......... .......... .......... .......... 53%  105M 4s
    ##  72650K .......... .......... .......... .......... .......... 53%  122M 4s
    ##  72700K .......... .......... .......... .......... .......... 53%  104M 4s
    ##  72750K .......... .......... .......... .......... .......... 54% 52.5M 4s
    ##  72800K .......... .......... .......... .......... .......... 54%  107M 4s
    ##  72850K .......... .......... .......... .......... .......... 54%  112M 4s
    ##  72900K .......... .......... .......... .......... .......... 54%  110M 4s
    ##  72950K .......... .......... .......... .......... .......... 54%  120M 4s
    ##  73000K .......... .......... .......... .......... .......... 54%  109M 4s
    ##  73050K .......... .......... .......... .......... .......... 54% 93.9M 4s
    ##  73100K .......... .......... .......... .......... .......... 54% 56.8M 4s
    ##  73150K .......... .......... .......... .......... .......... 54%  107M 4s
    ##  73200K .......... .......... .......... .......... .......... 54%  108M 4s
    ##  73250K .......... .......... .......... .......... .......... 54%  124M 4s
    ##  73300K .......... .......... .......... .......... .......... 54%  108M 4s
    ##  73350K .......... .......... .......... .......... .......... 54%  140M 4s
    ##  73400K .......... .......... .......... .......... .......... 54% 91.1M 4s
    ##  73450K .......... .......... .......... .......... .......... 54%  104M 4s
    ##  73500K .......... .......... .......... .......... .......... 54% 13.6M 4s
    ##  73550K .......... .......... .......... .......... .......... 54% 83.8M 4s
    ##  73600K .......... .......... .......... .......... .......... 54% 80.3M 4s
    ##  73650K .......... .......... .......... .......... .......... 54%  121M 4s
    ##  73700K .......... .......... .......... .......... .......... 54% 87.7M 4s
    ##  73750K .......... .......... .......... .......... .......... 54% 84.8M 4s
    ##  73800K .......... .......... .......... .......... .......... 54% 92.5M 4s
    ##  73850K .......... .......... .......... .......... .......... 54% 48.7M 4s
    ##  73900K .......... .......... .......... .......... .......... 54% 96.8M 4s
    ##  73950K .......... .......... .......... .......... .......... 54% 42.3M 4s
    ##  74000K .......... .......... .......... .......... .......... 54%  103M 4s
    ##  74050K .......... .......... .......... .......... .......... 54%  104M 4s
    ##  74100K .......... .......... .......... .......... .......... 55% 85.3M 4s
    ##  74150K .......... .......... .......... .......... .......... 55% 94.4M 4s
    ##  74200K .......... .......... .......... .......... .......... 55% 99.6M 4s
    ##  74250K .......... .......... .......... .......... .......... 55%  122M 4s
    ##  74300K .......... .......... .......... .......... .......... 55%  108M 4s
    ##  74350K .......... .......... .......... .......... .......... 55%  110M 4s
    ##  74400K .......... .......... .......... .......... .......... 55% 40.2M 4s
    ##  74450K .......... .......... .......... .......... .......... 55% 46.8M 4s
    ##  74500K .......... .......... .......... .......... .......... 55% 62.4M 4s
    ##  74550K .......... .......... .......... .......... .......... 55% 93.8M 4s
    ##  74600K .......... .......... .......... .......... .......... 55%  105M 4s
    ##  74650K .......... .......... .......... .......... .......... 55% 6.16M 4s
    ##  74700K .......... .......... .......... .......... .......... 55% 36.6M 4s
    ##  74750K .......... .......... .......... .......... .......... 55% 22.0M 4s
    ##  74800K .......... .......... .......... .......... .......... 55%  100M 4s
    ##  74850K .......... .......... .......... .......... .......... 55% 49.3M 4s
    ##  74900K .......... .......... .......... .......... .......... 55% 61.4M 4s
    ##  74950K .......... .......... .......... .......... .......... 55% 76.9M 4s
    ##  75000K .......... .......... .......... .......... .......... 55%  106M 4s
    ##  75050K .......... .......... .......... .......... .......... 55%  122M 4s
    ##  75100K .......... .......... .......... .......... .......... 55% 23.9M 4s
    ##  75150K .......... .......... .......... .......... .......... 55%  132M 4s
    ##  75200K .......... .......... .......... .......... .......... 55% 44.9M 4s
    ##  75250K .......... .......... .......... .......... .......... 55% 24.6M 4s
    ##  75300K .......... .......... .......... .......... .......... 55% 69.6M 4s
    ##  75350K .......... .......... .......... .......... .......... 55% 58.3M 4s
    ##  75400K .......... .......... .......... .......... .......... 55% 20.4M 4s
    ##  75450K .......... .......... .......... .......... .......... 56%  121M 4s
    ##  75500K .......... .......... .......... .......... .......... 56% 60.1M 4s
    ##  75550K .......... .......... .......... .......... .......... 56%  108M 4s
    ##  75600K .......... .......... .......... .......... .......... 56% 41.2M 4s
    ##  75650K .......... .......... .......... .......... .......... 56%  114M 4s
    ##  75700K .......... .......... .......... .......... .......... 56%  112M 4s
    ##  75750K .......... .......... .......... .......... .......... 56%  139M 4s
    ##  75800K .......... .......... .......... .......... .......... 56%  106M 4s
    ##  75850K .......... .......... .......... .......... .......... 56%  120M 4s
    ##  75900K .......... .......... .......... .......... .......... 56% 21.7M 4s
    ##  75950K .......... .......... .......... .......... .......... 56% 81.8M 4s
    ##  76000K .......... .......... .......... .......... .......... 56% 33.3M 4s
    ##  76050K .......... .......... .......... .......... .......... 56% 23.8M 4s
    ##  76100K .......... .......... .......... .......... .......... 56% 93.6M 4s
    ##  76150K .......... .......... .......... .......... .......... 56% 60.4M 4s
    ##  76200K .......... .......... .......... .......... .......... 56% 70.6M 4s
    ##  76250K .......... .......... .......... .......... .......... 56% 84.5M 4s
    ##  76300K .......... .......... .......... .......... .......... 56% 99.0M 4s
    ##  76350K .......... .......... .......... .......... .......... 56%  128M 4s
    ##  76400K .......... .......... .......... .......... .......... 56% 14.1M 4s
    ##  76450K .......... .......... .......... .......... .......... 56% 88.7M 4s
    ##  76500K .......... .......... .......... .......... .......... 56% 22.4M 4s
    ##  76550K .......... .......... .......... .......... .......... 56%  114M 4s
    ##  76600K .......... .......... .......... .......... .......... 56% 90.9M 4s
    ##  76650K .......... .......... .......... .......... .......... 56% 58.6M 4s
    ##  76700K .......... .......... .......... .......... .......... 56% 98.2M 4s
    ##  76750K .......... .......... .......... .......... .......... 56%  109M 4s
    ##  76800K .......... .......... .......... .......... .......... 57% 99.9M 4s
    ##  76850K .......... .......... .......... .......... .......... 57%  126M 4s
    ##  76900K .......... .......... .......... .......... .......... 57% 82.5M 4s
    ##  76950K .......... .......... .......... .......... .......... 57% 29.2M 4s
    ##  77000K .......... .......... .......... .......... .......... 57%  108M 4s
    ##  77050K .......... .......... .......... .......... .......... 57% 41.1M 4s
    ##  77100K .......... .......... .......... .......... .......... 57%  115M 4s
    ##  77150K .......... .......... .......... .......... .......... 57% 98.3M 4s
    ##  77200K .......... .......... .......... .......... .......... 57% 58.3M 4s
    ##  77250K .......... .......... .......... .......... .......... 57%  111M 4s
    ##  77300K .......... .......... .......... .......... .......... 57%  106M 4s
    ##  77350K .......... .......... .......... .......... .......... 57% 77.7M 4s
    ##  77400K .......... .......... .......... .......... .......... 57% 82.5M 4s
    ##  77450K .......... .......... .......... .......... .......... 57%  126M 4s
    ##  77500K .......... .......... .......... .......... .......... 57%  104M 4s
    ##  77550K .......... .......... .......... .......... .......... 57%  110M 3s
    ##  77600K .......... .......... .......... .......... .......... 57% 94.7M 3s
    ##  77650K .......... .......... .......... .......... .......... 57% 27.1M 3s
    ##  77700K .......... .......... .......... .......... .......... 57% 73.5M 3s
    ##  77750K .......... .......... .......... .......... .......... 57%  116M 3s
    ##  77800K .......... .......... .......... .......... .......... 57% 74.9M 3s
    ##  77850K .......... .......... .......... .......... .......... 57% 61.4M 3s
    ##  77900K .......... .......... .......... .......... .......... 57% 62.8M 3s
    ##  77950K .......... .......... .......... .......... .......... 57% 63.6M 3s
    ##  78000K .......... .......... .......... .......... .......... 57% 81.3M 3s
    ##  78050K .......... .......... .......... .......... .......... 57%  135M 3s
    ##  78100K .......... .......... .......... .......... .......... 58% 33.9M 3s
    ##  78150K .......... .......... .......... .......... .......... 58% 77.1M 3s
    ##  78200K .......... .......... .......... .......... .......... 58% 68.6M 3s
    ##  78250K .......... .......... .......... .......... .......... 58% 79.8M 3s
    ##  78300K .......... .......... .......... .......... .......... 58% 82.7M 3s
    ##  78350K .......... .......... .......... .......... .......... 58% 35.4M 3s
    ##  78400K .......... .......... .......... .......... .......... 58% 74.6M 3s
    ##  78450K .......... .......... .......... .......... .......... 58% 94.6M 3s
    ##  78500K .......... .......... .......... .......... .......... 58%  114M 3s
    ##  78550K .......... .......... .......... .......... .......... 58% 53.8M 3s
    ##  78600K .......... .......... .......... .......... .......... 58% 98.5M 3s
    ##  78650K .......... .......... .......... .......... .......... 58% 71.8M 3s
    ##  78700K .......... .......... .......... .......... .......... 58%  116M 3s
    ##  78750K .......... .......... .......... .......... .......... 58% 66.7M 3s
    ##  78800K .......... .......... .......... .......... .......... 58% 30.9M 3s
    ##  78850K .......... .......... .......... .......... .......... 58%  101M 3s
    ##  78900K .......... .......... .......... .......... .......... 58% 96.3M 3s
    ##  78950K .......... .......... .......... .......... .......... 58% 41.0M 3s
    ##  79000K .......... .......... .......... .......... .......... 58% 89.0M 3s
    ##  79050K .......... .......... .......... .......... .......... 58% 79.7M 3s
    ##  79100K .......... .......... .......... .......... .......... 58% 58.2M 3s
    ##  79150K .......... .......... .......... .......... .......... 58% 98.3M 3s
    ##  79200K .......... .......... .......... .......... .......... 58% 90.8M 3s
    ##  79250K .......... .......... .......... .......... .......... 58% 57.6M 3s
    ##  79300K .......... .......... .......... .......... .......... 58% 61.6M 3s
    ##  79350K .......... .......... .......... .......... .......... 58% 79.2M 3s
    ##  79400K .......... .......... .......... .......... .......... 58% 81.4M 3s
    ##  79450K .......... .......... .......... .......... .......... 59% 90.1M 3s
    ##  79500K .......... .......... .......... .......... .......... 59% 48.1M 3s
    ##  79550K .......... .......... .......... .......... .......... 59% 81.5M 3s
    ##  79600K .......... .......... .......... .......... .......... 59%  105M 3s
    ##  79650K .......... .......... .......... .......... .......... 59% 77.4M 3s
    ##  79700K .......... .......... .......... .......... .......... 59% 69.0M 3s
    ##  79750K .......... .......... .......... .......... .......... 59% 45.6M 3s
    ##  79800K .......... .......... .......... .......... .......... 59% 72.1M 3s
    ##  79850K .......... .......... .......... .......... .......... 59% 90.9M 3s
    ##  79900K .......... .......... .......... .......... .......... 59%  135M 3s
    ##  79950K .......... .......... .......... .......... .......... 59% 60.9M 3s
    ##  80000K .......... .......... .......... .......... .......... 59% 88.9M 3s
    ##  80050K .......... .......... .......... .......... .......... 59%  109M 3s
    ##  80100K .......... .......... .......... .......... .......... 59% 42.4M 3s
    ##  80150K .......... .......... .......... .......... .......... 59% 68.8M 3s
    ##  80200K .......... .......... .......... .......... .......... 59% 51.8M 3s
    ##  80250K .......... .......... .......... .......... .......... 59% 89.1M 3s
    ##  80300K .......... .......... .......... .......... .......... 59%  114M 3s
    ##  80350K .......... .......... .......... .......... .......... 59% 14.5M 3s
    ##  80400K .......... .......... .......... .......... .......... 59% 95.8M 3s
    ##  80450K .......... .......... .......... .......... .......... 59% 22.7M 3s
    ##  80500K .......... .......... .......... .......... .......... 59%  121M 3s
    ##  80550K .......... .......... .......... .......... .......... 59% 47.1M 3s
    ##  80600K .......... .......... .......... .......... .......... 59% 91.3M 3s
    ##  80650K .......... .......... .......... .......... .......... 59%  126M 3s
    ##  80700K .......... .......... .......... .......... .......... 59%  122M 3s
    ##  80750K .......... .......... .......... .......... .......... 59%  103M 3s
    ##  80800K .......... .......... .......... .......... .......... 60%  103M 3s
    ##  80850K .......... .......... .......... .......... .......... 60% 30.8M 3s
    ##  80900K .......... .......... .......... .......... .......... 60% 77.8M 3s
    ##  80950K .......... .......... .......... .......... .......... 60% 96.7M 3s
    ##  81000K .......... .......... .......... .......... .......... 60% 52.7M 3s
    ##  81050K .......... .......... .......... .......... .......... 60% 40.2M 3s
    ##  81100K .......... .......... .......... .......... .......... 60% 62.0M 3s
    ##  81150K .......... .......... .......... .......... .......... 60% 92.6M 3s
    ##  81200K .......... .......... .......... .......... .......... 60% 68.1M 3s
    ##  81250K .......... .......... .......... .......... .......... 60% 23.7M 3s
    ##  81300K .......... .......... .......... .......... .......... 60%  109M 3s
    ##  81350K .......... .......... .......... .......... .......... 60% 32.2M 3s
    ##  81400K .......... .......... .......... .......... .......... 60% 89.8M 3s
    ##  81450K .......... .......... .......... .......... .......... 60%  111M 3s
    ##  81500K .......... .......... .......... .......... .......... 60%  105M 3s
    ##  81550K .......... .......... .......... .......... .......... 60% 94.0M 3s
    ##  81600K .......... .......... .......... .......... .......... 60% 83.7M 3s
    ##  81650K .......... .......... .......... .......... .......... 60% 93.1M 3s
    ##  81700K .......... .......... .......... .......... .......... 60% 94.7M 3s
    ##  81750K .......... .......... .......... .......... .......... 60%  122M 3s
    ##  81800K .......... .......... .......... .......... .......... 60% 94.9M 3s
    ##  81850K .......... .......... .......... .......... .......... 60% 94.5M 3s
    ##  81900K .......... .......... .......... .......... .......... 60% 84.6M 3s
    ##  81950K .......... .......... .......... .......... .......... 60% 28.4M 3s
    ##  82000K .......... .......... .......... .......... .......... 60% 34.0M 3s
    ##  82050K .......... .......... .......... .......... .......... 60%  104M 3s
    ##  82100K .......... .......... .......... .......... .......... 60% 77.6M 3s
    ##  82150K .......... .......... .......... .......... .......... 61% 87.3M 3s
    ##  82200K .......... .......... .......... .......... .......... 61% 96.4M 3s
    ##  82250K .......... .......... .......... .......... .......... 61%  115M 3s
    ##  82300K .......... .......... .......... .......... .......... 61% 58.7M 3s
    ##  82350K .......... .......... .......... .......... .......... 61% 77.8M 3s
    ##  82400K .......... .......... .......... .......... .......... 61% 82.7M 3s
    ##  82450K .......... .......... .......... .......... .......... 61% 33.3M 3s
    ##  82500K .......... .......... .......... .......... .......... 61% 59.5M 3s
    ##  82550K .......... .......... .......... .......... .......... 61% 93.1M 3s
    ##  82600K .......... .......... .......... .......... .......... 61% 98.4M 3s
    ##  82650K .......... .......... .......... .......... .......... 61% 38.3M 3s
    ##  82700K .......... .......... .......... .......... .......... 61% 89.5M 3s
    ##  82750K .......... .......... .......... .......... .......... 61% 47.7M 3s
    ##  82800K .......... .......... .......... .......... .......... 61% 81.6M 3s
    ##  82850K .......... .......... .......... .......... .......... 61%  125M 3s
    ##  82900K .......... .......... .......... .......... .......... 61% 87.2M 3s
    ##  82950K .......... .......... .......... .......... .......... 61% 56.7M 3s
    ##  83000K .......... .......... .......... .......... .......... 61% 78.1M 3s
    ##  83050K .......... .......... .......... .......... .......... 61% 79.2M 3s
    ##  83100K .......... .......... .......... .......... .......... 61% 82.7M 3s
    ##  83150K .......... .......... .......... .......... .......... 61% 54.0M 3s
    ##  83200K .......... .......... .......... .......... .......... 61% 64.7M 3s
    ##  83250K .......... .......... .......... .......... .......... 61% 66.3M 3s
    ##  83300K .......... .......... .......... .......... .......... 61% 93.2M 3s
    ##  83350K .......... .......... .......... .......... .......... 61% 50.4M 3s
    ##  83400K .......... .......... .......... .......... .......... 61% 64.2M 3s
    ##  83450K .......... .......... .......... .......... .......... 61% 93.7M 3s
    ##  83500K .......... .......... .......... .......... .......... 62% 77.2M 3s
    ##  83550K .......... .......... .......... .......... .......... 62% 75.2M 3s
    ##  83600K .......... .......... .......... .......... .......... 62% 61.6M 3s
    ##  83650K .......... .......... .......... .......... .......... 62% 65.2M 3s
    ##  83700K .......... .......... .......... .......... .......... 62% 76.7M 3s
    ##  83750K .......... .......... .......... .......... .......... 62% 98.2M 3s
    ##  83800K .......... .......... .......... .......... .......... 62% 58.0M 3s
    ##  83850K .......... .......... .......... .......... .......... 62% 83.0M 3s
    ##  83900K .......... .......... .......... .......... .......... 62% 80.1M 3s
    ##  83950K .......... .......... .......... .......... .......... 62% 87.1M 3s
    ##  84000K .......... .......... .......... .......... .......... 62% 68.7M 3s
    ##  84050K .......... .......... .......... .......... .......... 62% 58.3M 3s
    ##  84100K .......... .......... .......... .......... .......... 62% 78.2M 3s
    ##  84150K .......... .......... .......... .......... .......... 62% 90.7M 3s
    ##  84200K .......... .......... .......... .......... .......... 62%  106M 3s
    ##  84250K .......... .......... .......... .......... .......... 62% 50.5M 3s
    ##  84300K .......... .......... .......... .......... .......... 62% 59.4M 3s
    ##  84350K .......... .......... .......... .......... .......... 62%  105M 3s
    ##  84400K .......... .......... .......... .......... .......... 62% 67.4M 3s
    ##  84450K .......... .......... .......... .......... .......... 62%  115M 3s
    ##  84500K .......... .......... .......... .......... .......... 62% 24.3M 3s
    ##  84550K .......... .......... .......... .......... .......... 62% 78.7M 3s
    ##  84600K .......... .......... .......... .......... .......... 62%  127M 3s
    ##  84650K .......... .......... .......... .......... .......... 62%  134M 3s
    ##  84700K .......... .......... .......... .......... .......... 62%  118M 3s
    ##  84750K .......... .......... .......... .......... .......... 62%  139M 3s
    ##  84800K .......... .......... .......... .......... .......... 62% 69.8M 3s
    ##  84850K .......... .......... .......... .......... .......... 63% 86.6M 3s
    ##  84900K .......... .......... .......... .......... .......... 63%  110M 3s
    ##  84950K .......... .......... .......... .......... .......... 63%  130M 3s
    ##  85000K .......... .......... .......... .......... .......... 63% 97.8M 3s
    ##  85050K .......... .......... .......... .......... .......... 63%  100M 3s
    ##  85100K .......... .......... .......... .......... .......... 63% 43.8M 3s
    ##  85150K .......... .......... .......... .......... .......... 63% 82.1M 3s
    ##  85200K .......... .......... .......... .......... .......... 63% 40.9M 3s
    ##  85250K .......... .......... .......... .......... .......... 63% 26.2M 3s
    ##  85300K .......... .......... .......... .......... .......... 63% 86.1M 3s
    ##  85350K .......... .......... .......... .......... .......... 63% 60.8M 3s
    ##  85400K .......... .......... .......... .......... .......... 63%  105M 3s
    ##  85450K .......... .......... .......... .......... .......... 63% 75.3M 3s
    ##  85500K .......... .......... .......... .......... .......... 63% 91.0M 3s
    ##  85550K .......... .......... .......... .......... .......... 63%  135M 3s
    ##  85600K .......... .......... .......... .......... .......... 63% 65.2M 3s
    ##  85650K .......... .......... .......... .......... .......... 63%  100M 3s
    ##  85700K .......... .......... .......... .......... .......... 63% 32.0M 3s
    ##  85750K .......... .......... .......... .......... .......... 63% 87.8M 3s
    ##  85800K .......... .......... .......... .......... .......... 63%  102M 3s
    ##  85850K .......... .......... .......... .......... .......... 63% 71.9M 3s
    ##  85900K .......... .......... .......... .......... .......... 63% 58.9M 3s
    ##  85950K .......... .......... .......... .......... .......... 63%  122M 3s
    ##  86000K .......... .......... .......... .......... .......... 63%  135M 3s
    ##  86050K .......... .......... .......... .......... .......... 63% 36.0M 3s
    ##  86100K .......... .......... .......... .......... .......... 63% 26.3M 3s
    ##  86150K .......... .......... .......... .......... .......... 63%  123M 3s
    ##  86200K .......... .......... .......... .......... .......... 64% 40.8M 3s
    ##  86250K .......... .......... .......... .......... .......... 64%  119M 3s
    ##  86300K .......... .......... .......... .......... .......... 64% 72.1M 3s
    ##  86350K .......... .......... .......... .......... .......... 64% 66.2M 3s
    ##  86400K .......... .......... .......... .......... .......... 64% 83.7M 3s
    ##  86450K .......... .......... .......... .......... .......... 64%  139M 3s
    ##  86500K .......... .......... .......... .......... .......... 64% 76.1M 3s
    ##  86550K .......... .......... .......... .......... .......... 64%  132M 3s
    ##  86600K .......... .......... .......... .......... .......... 64% 86.1M 3s
    ##  86650K .......... .......... .......... .......... .......... 64% 49.8M 3s
    ##  86700K .......... .......... .......... .......... .......... 64% 51.8M 3s
    ##  86750K .......... .......... .......... .......... .......... 64% 92.0M 3s
    ##  86800K .......... .......... .......... .......... .......... 64% 54.3M 3s
    ##  86850K .......... .......... .......... .......... .......... 64% 59.3M 3s
    ##  86900K .......... .......... .......... .......... .......... 64%  103M 3s
    ##  86950K .......... .......... .......... .......... .......... 64% 68.6M 3s
    ##  87000K .......... .......... .......... .......... .......... 64% 88.4M 3s
    ##  87050K .......... .......... .......... .......... .......... 64%  103M 3s
    ##  87100K .......... .......... .......... .......... .......... 64% 73.8M 3s
    ##  87150K .......... .......... .......... .......... .......... 64% 58.3M 3s
    ##  87200K .......... .......... .......... .......... .......... 64%  104M 3s
    ##  87250K .......... .......... .......... .......... .......... 64% 64.3M 3s
    ##  87300K .......... .......... .......... .......... .......... 64% 69.7M 3s
    ##  87350K .......... .......... .......... .......... .......... 64% 92.5M 3s
    ##  87400K .......... .......... .......... .......... .......... 64% 40.4M 3s
    ##  87450K .......... .......... .......... .......... .......... 64% 74.8M 3s
    ##  87500K .......... .......... .......... .......... .......... 64% 99.0M 3s
    ##  87550K .......... .......... .......... .......... .......... 65%  115M 3s
    ##  87600K .......... .......... .......... .......... .......... 65% 74.5M 3s
    ##  87650K .......... .......... .......... .......... .......... 65%  104M 3s
    ##  87700K .......... .......... .......... .......... .......... 65% 89.4M 3s
    ##  87750K .......... .......... .......... .......... .......... 65%  141M 3s
    ##  87800K .......... .......... .......... .......... .......... 65% 81.6M 3s
    ##  87850K .......... .......... .......... .......... .......... 65% 37.9M 3s
    ##  87900K .......... .......... .......... .......... .......... 65% 95.6M 3s
    ##  87950K .......... .......... .......... .......... .......... 65% 89.8M 3s
    ##  88000K .......... .......... .......... .......... .......... 65%  113M 3s
    ##  88050K .......... .......... .......... .......... .......... 65% 41.2M 3s
    ##  88100K .......... .......... .......... .......... .......... 65% 79.7M 3s
    ##  88150K .......... .......... .......... .......... .......... 65%  101M 3s
    ##  88200K .......... .......... .......... .......... .......... 65%  123M 3s
    ##  88250K .......... .......... .......... .......... .......... 65%  117M 3s
    ##  88300K .......... .......... .......... .......... .......... 65% 48.3M 3s
    ##  88350K .......... .......... .......... .......... .......... 65%  115M 3s
    ##  88400K .......... .......... .......... .......... .......... 65%  122M 3s
    ##  88450K .......... .......... .......... .......... .......... 65% 87.4M 3s
    ##  88500K .......... .......... .......... .......... .......... 65%  108M 3s
    ##  88550K .......... .......... .......... .......... .......... 65% 68.1M 3s
    ##  88600K .......... .......... .......... .......... .......... 65% 86.8M 3s
    ##  88650K .......... .......... .......... .......... .......... 65%  125M 3s
    ##  88700K .......... .......... .......... .......... .......... 65% 78.1M 3s
    ##  88750K .......... .......... .......... .......... .......... 65% 94.2M 3s
    ##  88800K .......... .......... .......... .......... .......... 65% 51.0M 3s
    ##  88850K .......... .......... .......... .......... .......... 65%  125M 3s
    ##  88900K .......... .......... .......... .......... .......... 66% 85.5M 3s
    ##  88950K .......... .......... .......... .......... .......... 66% 95.5M 3s
    ##  89000K .......... .......... .......... .......... .......... 66% 99.9M 3s
    ##  89050K .......... .......... .......... .......... .......... 66% 54.9M 3s
    ##  89100K .......... .......... .......... .......... .......... 66% 92.5M 3s
    ##  89150K .......... .......... .......... .......... .......... 66%  115M 3s
    ##  89200K .......... .......... .......... .......... .......... 66% 75.9M 3s
    ##  89250K .......... .......... .......... .......... .......... 66%  111M 3s
    ##  89300K .......... .......... .......... .......... .......... 66%  141M 2s
    ##  89350K .......... .......... .......... .......... .......... 66% 53.7M 2s
    ##  89400K .......... .......... .......... .......... .......... 66%  106M 2s
    ##  89450K .......... .......... .......... .......... .......... 66% 90.2M 2s
    ##  89500K .......... .......... .......... .......... .......... 66% 85.6M 2s
    ##  89550K .......... .......... .......... .......... .......... 66%  115M 2s
    ##  89600K .......... .......... .......... .......... .......... 66%  112M 2s
    ##  89650K .......... .......... .......... .......... .......... 66% 98.9M 2s
    ##  89700K .......... .......... .......... .......... .......... 66%  113M 2s
    ##  89750K .......... .......... .......... .......... .......... 66% 82.8M 2s
    ##  89800K .......... .......... .......... .......... .......... 66% 45.4M 2s
    ##  89850K .......... .......... .......... .......... .......... 66% 88.4M 2s
    ##  89900K .......... .......... .......... .......... .......... 66%  126M 2s
    ##  89950K .......... .......... .......... .......... .......... 66%  100M 2s
    ##  90000K .......... .......... .......... .......... .......... 66% 49.0M 2s
    ##  90050K .......... .......... .......... .......... .......... 66%  104M 2s
    ##  90100K .......... .......... .......... .......... .......... 66%  107M 2s
    ##  90150K .......... .......... .......... .......... .......... 66%  111M 2s
    ##  90200K .......... .......... .......... .......... .......... 66% 92.4M 2s
    ##  90250K .......... .......... .......... .......... .......... 67%  144M 2s
    ##  90300K .......... .......... .......... .......... .......... 67% 71.3M 2s
    ##  90350K .......... .......... .......... .......... .......... 67% 90.6M 2s
    ##  90400K .......... .......... .......... .......... .......... 67%  107M 2s
    ##  90450K .......... .......... .......... .......... .......... 67%  120M 2s
    ##  90500K .......... .......... .......... .......... .......... 67%  138M 2s
    ##  90550K .......... .......... .......... .......... .......... 67% 68.7M 2s
    ##  90600K .......... .......... .......... .......... .......... 67%  103M 2s
    ##  90650K .......... .......... .......... .......... .......... 67% 96.2M 2s
    ##  90700K .......... .......... .......... .......... .......... 67%  110M 2s
    ##  90750K .......... .......... .......... .......... .......... 67%  148M 2s
    ##  90800K .......... .......... .......... .......... .......... 67% 22.8M 2s
    ##  90850K .......... .......... .......... .......... .......... 67% 38.1M 2s
    ##  90900K .......... .......... .......... .......... .......... 67% 55.8M 2s
    ##  90950K .......... .......... .......... .......... .......... 67% 71.7M 2s
    ##  91000K .......... .......... .......... .......... .......... 67%  142M 2s
    ##  91050K .......... .......... .......... .......... .......... 67% 64.0M 2s
    ##  91100K .......... .......... .......... .......... .......... 67% 61.5M 2s
    ##  91150K .......... .......... .......... .......... .......... 67%  134M 2s
    ##  91200K .......... .......... .......... .......... .......... 67%  140M 2s
    ##  91250K .......... .......... .......... .......... .......... 67% 42.2M 2s
    ##  91300K .......... .......... .......... .......... .......... 67%  110M 2s
    ##  91350K .......... .......... .......... .......... .......... 67%  160M 2s
    ##  91400K .......... .......... .......... .......... .......... 67%  144M 2s
    ##  91450K .......... .......... .......... .......... .......... 67%  164M 2s
    ##  91500K .......... .......... .......... .......... .......... 67% 40.5M 2s
    ##  91550K .......... .......... .......... .......... .......... 67% 31.9M 2s
    ##  91600K .......... .......... .......... .......... .......... 68% 87.9M 2s
    ##  91650K .......... .......... .......... .......... .......... 68% 71.7M 2s
    ##  91700K .......... .......... .......... .......... .......... 68% 91.7M 2s
    ##  91750K .......... .......... .......... .......... .......... 68% 42.6M 2s
    ##  91800K .......... .......... .......... .......... .......... 68% 22.0M 2s
    ##  91850K .......... .......... .......... .......... .......... 68%  101M 2s
    ##  91900K .......... .......... .......... .......... .......... 68%  138M 2s
    ##  91950K .......... .......... .......... .......... .......... 68% 18.6M 2s
    ##  92000K .......... .......... .......... .......... .......... 68%  136M 2s
    ##  92050K .......... .......... .......... .......... .......... 68%  154M 2s
    ##  92100K .......... .......... .......... .......... .......... 68%  146M 2s
    ##  92150K .......... .......... .......... .......... .......... 68% 3.84M 2s
    ##  92200K .......... .......... .......... .......... .......... 68% 67.7M 2s
    ##  92250K .......... .......... .......... .......... .......... 68% 48.7M 2s
    ##  92300K .......... .......... .......... .......... .......... 68%  125M 2s
    ##  92350K .......... .......... .......... .......... .......... 68% 55.0M 2s
    ##  92400K .......... .......... .......... .......... .......... 68%  126M 2s
    ##  92450K .......... .......... .......... .......... .......... 68% 59.6M 2s
    ##  92500K .......... .......... .......... .......... .......... 68%  150M 2s
    ##  92550K .......... .......... .......... .......... .......... 68% 21.9M 2s
    ##  92600K .......... .......... .......... .......... .......... 68% 43.0M 2s
    ##  92650K .......... .......... .......... .......... .......... 68% 34.0M 2s
    ##  92700K .......... .......... .......... .......... .......... 68% 33.1M 2s
    ##  92750K .......... .......... .......... .......... .......... 68% 41.4M 2s
    ##  92800K .......... .......... .......... .......... .......... 68% 17.8M 2s
    ##  92850K .......... .......... .......... .......... .......... 68% 75.0M 2s
    ##  92900K .......... .......... .......... .......... .......... 68%  107M 2s
    ##  92950K .......... .......... .......... .......... .......... 69%  140M 2s
    ##  93000K .......... .......... .......... .......... .......... 69%  139M 2s
    ##  93050K .......... .......... .......... .......... .......... 69%  169M 2s
    ##  93100K .......... .......... .......... .......... .......... 69%  118M 2s
    ##  93150K .......... .......... .......... .......... .......... 69%  130M 2s
    ##  93200K .......... .......... .......... .......... .......... 69% 10.3M 2s
    ##  93250K .......... .......... .......... .......... .......... 69%  110M 2s
    ##  93300K .......... .......... .......... .......... .......... 69%  132M 2s
    ##  93350K .......... .......... .......... .......... .......... 69%  148M 2s
    ##  93400K .......... .......... .......... .......... .......... 69%  108M 2s
    ##  93450K .......... .......... .......... .......... .......... 69%  111M 2s
    ##  93500K .......... .......... .......... .......... .......... 69% 85.5M 2s
    ##  93550K .......... .......... .......... .......... .......... 69%  123M 2s
    ##  93600K .......... .......... .......... .......... .......... 69%  115M 2s
    ##  93650K .......... .......... .......... .......... .......... 69% 91.3M 2s
    ##  93700K .......... .......... .......... .......... .......... 69%  119M 2s
    ##  93750K .......... .......... .......... .......... .......... 69%  118M 2s
    ##  93800K .......... .......... .......... .......... .......... 69%  132M 2s
    ##  93850K .......... .......... .......... .......... .......... 69% 94.4M 2s
    ##  93900K .......... .......... .......... .......... .......... 69%  110M 2s
    ##  93950K .......... .......... .......... .......... .......... 69% 47.8M 2s
    ##  94000K .......... .......... .......... .......... .......... 69% 59.1M 2s
    ##  94050K .......... .......... .......... .......... .......... 69% 86.8M 2s
    ##  94100K .......... .......... .......... .......... .......... 69%  110M 2s
    ##  94150K .......... .......... .......... .......... .......... 69% 48.1M 2s
    ##  94200K .......... .......... .......... .......... .......... 69% 88.9M 2s
    ##  94250K .......... .......... .......... .......... .......... 69%  163M 2s
    ##  94300K .......... .......... .......... .......... .......... 70% 19.7M 2s
    ##  94350K .......... .......... .......... .......... .......... 70% 92.4M 2s
    ##  94400K .......... .......... .......... .......... .......... 70% 10.9M 2s
    ##  94450K .......... .......... .......... .......... .......... 70%  127M 2s
    ##  94500K .......... .......... .......... .......... .......... 70% 26.6M 2s
    ##  94550K .......... .......... .......... .......... .......... 70%  139M 2s
    ##  94600K .......... .......... .......... .......... .......... 70% 48.8M 2s
    ##  94650K .......... .......... .......... .......... .......... 70% 75.6M 2s
    ##  94700K .......... .......... .......... .......... .......... 70% 81.1M 2s
    ##  94750K .......... .......... .......... .......... .......... 70%  101M 2s
    ##  94800K .......... .......... .......... .......... .......... 70%  108M 2s
    ##  94850K .......... .......... .......... .......... .......... 70%  139M 2s
    ##  94900K .......... .......... .......... .......... .......... 70%  117M 2s
    ##  94950K .......... .......... .......... .......... .......... 70%  130M 2s
    ##  95000K .......... .......... .......... .......... .......... 70% 88.3M 2s
    ##  95050K .......... .......... .......... .......... .......... 70% 42.9M 2s
    ##  95100K .......... .......... .......... .......... .......... 70% 89.7M 2s
    ##  95150K .......... .......... .......... .......... .......... 70% 23.5M 2s
    ##  95200K .......... .......... .......... .......... .......... 70% 82.6M 2s
    ##  95250K .......... .......... .......... .......... .......... 70% 80.5M 2s
    ##  95300K .......... .......... .......... .......... .......... 70% 91.0M 2s
    ##  95350K .......... .......... .......... .......... .......... 70% 87.9M 2s
    ##  95400K .......... .......... .......... .......... .......... 70%  129M 2s
    ##  95450K .......... .......... .......... .......... .......... 70% 86.1M 2s
    ##  95500K .......... .......... .......... .......... .......... 70% 46.4M 2s
    ##  95550K .......... .......... .......... .......... .......... 70%  108M 2s
    ##  95600K .......... .......... .......... .......... .......... 70% 91.1M 2s
    ##  95650K .......... .......... .......... .......... .......... 71% 70.6M 2s
    ##  95700K .......... .......... .......... .......... .......... 71% 81.1M 2s
    ##  95750K .......... .......... .......... .......... .......... 71%  136M 2s
    ##  95800K .......... .......... .......... .......... .......... 71% 32.8M 2s
    ##  95850K .......... .......... .......... .......... .......... 71% 96.8M 2s
    ##  95900K .......... .......... .......... .......... .......... 71% 99.8M 2s
    ##  95950K .......... .......... .......... .......... .......... 71%  106M 2s
    ##  96000K .......... .......... .......... .......... .......... 71%  133M 2s
    ##  96050K .......... .......... .......... .......... .......... 71% 5.46M 2s
    ##  96100K .......... .......... .......... .......... .......... 71% 75.7M 2s
    ##  96150K .......... .......... .......... .......... .......... 71% 32.4M 2s
    ##  96200K .......... .......... .......... .......... .......... 71% 20.6M 2s
    ##  96250K .......... .......... .......... .......... .......... 71% 55.5M 2s
    ##  96300K .......... .......... .......... .......... .......... 71% 87.0M 2s
    ##  96350K .......... .......... .......... .......... .......... 71%  121M 2s
    ##  96400K .......... .......... .......... .......... .......... 71% 71.2M 2s
    ##  96450K .......... .......... .......... .......... .......... 71% 47.2M 2s
    ##  96500K .......... .......... .......... .......... .......... 71% 70.0M 2s
    ##  96550K .......... .......... .......... .......... .......... 71% 86.9M 2s
    ##  96600K .......... .......... .......... .......... .......... 71% 56.5M 2s
    ##  96650K .......... .......... .......... .......... .......... 71% 79.0M 2s
    ##  96700K .......... .......... .......... .......... .......... 71% 44.4M 2s
    ##  96750K .......... .......... .......... .......... .......... 71% 59.2M 2s
    ##  96800K .......... .......... .......... .......... .......... 71% 45.9M 2s
    ##  96850K .......... .......... .......... .......... .......... 71% 39.1M 2s
    ##  96900K .......... .......... .......... .......... .......... 71% 28.7M 2s
    ##  96950K .......... .......... .......... .......... .......... 71%  121M 2s
    ##  97000K .......... .......... .......... .......... .......... 72% 75.8M 2s
    ##  97050K .......... .......... .......... .......... .......... 72% 77.2M 2s
    ##  97100K .......... .......... .......... .......... .......... 72% 96.3M 2s
    ##  97150K .......... .......... .......... .......... .......... 72%  119M 2s
    ##  97200K .......... .......... .......... .......... .......... 72% 96.7M 2s
    ##  97250K .......... .......... .......... .......... .......... 72%  128M 2s
    ##  97300K .......... .......... .......... .......... .......... 72% 6.54M 2s
    ##  97350K .......... .......... .......... .......... .......... 72% 71.4M 2s
    ##  97400K .......... .......... .......... .......... .......... 72% 20.5M 2s
    ##  97450K .......... .......... .......... .......... .......... 72% 28.1M 2s
    ##  97500K .......... .......... .......... .......... .......... 72% 32.9M 2s
    ##  97550K .......... .......... .......... .......... .......... 72% 58.6M 2s
    ##  97600K .......... .......... .......... .......... .......... 72% 41.2M 2s
    ##  97650K .......... .......... .......... .......... .......... 72%  100M 2s
    ##  97700K .......... .......... .......... .......... .......... 72% 47.0M 2s
    ##  97750K .......... .......... .......... .......... .......... 72% 29.9M 2s
    ##  97800K .......... .......... .......... .......... .......... 72%  100M 2s
    ##  97850K .......... .......... .......... .......... .......... 72% 50.5M 2s
    ##  97900K .......... .......... .......... .......... .......... 72% 93.1M 2s
    ##  97950K .......... .......... .......... .......... .......... 72% 12.3M 2s
    ##  98000K .......... .......... .......... .......... .......... 72% 66.7M 2s
    ##  98050K .......... .......... .......... .......... .......... 72% 39.5M 2s
    ##  98100K .......... .......... .......... .......... .......... 72% 45.1M 2s
    ##  98150K .......... .......... .......... .......... .......... 72%  109M 2s
    ##  98200K .......... .......... .......... .......... .......... 72% 54.8M 2s
    ##  98250K .......... .......... .......... .......... .......... 72% 73.9M 2s
    ##  98300K .......... .......... .......... .......... .......... 72% 73.8M 2s
    ##  98350K .......... .......... .......... .......... .......... 73% 93.9M 2s
    ##  98400K .......... .......... .......... .......... .......... 73% 52.7M 2s
    ##  98450K .......... .......... .......... .......... .......... 73% 49.7M 2s
    ##  98500K .......... .......... .......... .......... .......... 73% 55.8M 2s
    ##  98550K .......... .......... .......... .......... .......... 73% 47.4M 2s
    ##  98600K .......... .......... .......... .......... .......... 73% 46.0M 2s
    ##  98650K .......... .......... .......... .......... .......... 73% 42.3M 2s
    ##  98700K .......... .......... .......... .......... .......... 73% 63.4M 2s
    ##  98750K .......... .......... .......... .......... .......... 73% 87.6M 2s
    ##  98800K .......... .......... .......... .......... .......... 73% 33.7M 2s
    ##  98850K .......... .......... .......... .......... .......... 73% 82.3M 2s
    ##  98900K .......... .......... .......... .......... .......... 73%  110M 2s
    ##  98950K .......... .......... .......... .......... .......... 73%  123M 2s
    ##  99000K .......... .......... .......... .......... .......... 73%  111M 2s
    ##  99050K .......... .......... .......... .......... .......... 73% 99.9M 2s
    ##  99100K .......... .......... .......... .......... .......... 73% 65.3M 2s
    ##  99150K .......... .......... .......... .......... .......... 73% 85.7M 2s
    ##  99200K .......... .......... .......... .......... .......... 73% 69.5M 2s
    ##  99250K .......... .......... .......... .......... .......... 73%  124M 2s
    ##  99300K .......... .......... .......... .......... .......... 73% 58.6M 2s
    ##  99350K .......... .......... .......... .......... .......... 73%  118M 2s
    ##  99400K .......... .......... .......... .......... .......... 73% 74.3M 2s
    ##  99450K .......... .......... .......... .......... .......... 73% 68.7M 2s
    ##  99500K .......... .......... .......... .......... .......... 73% 76.0M 2s
    ##  99550K .......... .......... .......... .......... .......... 73% 33.1M 2s
    ##  99600K .......... .......... .......... .......... .......... 73%  101M 2s
    ##  99650K .......... .......... .......... .......... .......... 73%  119M 2s
    ##  99700K .......... .......... .......... .......... .......... 74% 58.0M 2s
    ##  99750K .......... .......... .......... .......... .......... 74% 45.4M 2s
    ##  99800K .......... .......... .......... .......... .......... 74% 73.9M 2s
    ##  99850K .......... .......... .......... .......... .......... 74% 49.1M 2s
    ##  99900K .......... .......... .......... .......... .......... 74% 97.1M 2s
    ##  99950K .......... .......... .......... .......... .......... 74%  121M 2s
    ## 100000K .......... .......... .......... .......... .......... 74%  114M 2s
    ## 100050K .......... .......... .......... .......... .......... 74%  135M 2s
    ## 100100K .......... .......... .......... .......... .......... 74%  121M 2s
    ## 100150K .......... .......... .......... .......... .......... 74%  145M 2s
    ## 100200K .......... .......... .......... .......... .......... 74% 4.60M 2s
    ## 100250K .......... .......... .......... .......... .......... 74% 52.7M 2s
    ## 100300K .......... .......... .......... .......... .......... 74% 49.5M 2s
    ## 100350K .......... .......... .......... .......... .......... 74% 77.0M 2s
    ## 100400K .......... .......... .......... .......... .......... 74% 45.3M 2s
    ## 100450K .......... .......... .......... .......... .......... 74% 49.1M 2s
    ## 100500K .......... .......... .......... .......... .......... 74% 48.0M 2s
    ## 100550K .......... .......... .......... .......... .......... 74% 52.3M 2s
    ## 100600K .......... .......... .......... .......... .......... 74%  116M 2s
    ## 100650K .......... .......... .......... .......... .......... 74%  132M 2s
    ## 100700K .......... .......... .......... .......... .......... 74%  118M 2s
    ## 100750K .......... .......... .......... .......... .......... 74%  147M 2s
    ## 100800K .......... .......... .......... .......... .......... 74%  111M 2s
    ## 100850K .......... .......... .......... .......... .......... 74% 21.8M 2s
    ## 100900K .......... .......... .......... .......... .......... 74%  111M 2s
    ## 100950K .......... .......... .......... .......... .......... 74%  123M 2s
    ## 101000K .......... .......... .......... .......... .......... 74% 60.1M 2s
    ## 101050K .......... .......... .......... .......... .......... 75%  111M 2s
    ## 101100K .......... .......... .......... .......... .......... 75% 63.7M 2s
    ## 101150K .......... .......... .......... .......... .......... 75% 98.1M 2s
    ## 101200K .......... .......... .......... .......... .......... 75% 54.9M 2s
    ## 101250K .......... .......... .......... .......... .......... 75% 85.4M 2s
    ## 101300K .......... .......... .......... .......... .......... 75%  113M 2s
    ## 101350K .......... .......... .......... .......... .......... 75% 27.3M 2s
    ## 101400K .......... .......... .......... .......... .......... 75% 25.0M 2s
    ## 101450K .......... .......... .......... .......... .......... 75%  114M 2s
    ## 101500K .......... .......... .......... .......... .......... 75% 21.3M 2s
    ## 101550K .......... .......... .......... .......... .......... 75%  131M 2s
    ## 101600K .......... .......... .......... .......... .......... 75% 78.9M 2s
    ## 101650K .......... .......... .......... .......... .......... 75%  142M 2s
    ## 101700K .......... .......... .......... .......... .......... 75%  111M 2s
    ## 101750K .......... .......... .......... .......... .......... 75% 90.3M 2s
    ## 101800K .......... .......... .......... .......... .......... 75%  113M 2s
    ## 101850K .......... .......... .......... .......... .......... 75% 86.4M 2s
    ## 101900K .......... .......... .......... .......... .......... 75% 74.7M 2s
    ## 101950K .......... .......... .......... .......... .......... 75%  136M 2s
    ## 102000K .......... .......... .......... .......... .......... 75% 30.7M 2s
    ## 102050K .......... .......... .......... .......... .......... 75% 53.3M 2s
    ## 102100K .......... .......... .......... .......... .......... 75% 44.6M 2s
    ## 102150K .......... .......... .......... .......... .......... 75%  147M 2s
    ## 102200K .......... .......... .......... .......... .......... 75% 44.0M 2s
    ## 102250K .......... .......... .......... .......... .......... 75%  136M 2s
    ## 102300K .......... .......... .......... .......... .......... 75% 40.8M 2s
    ## 102350K .......... .......... .......... .......... .......... 75%  106M 2s
    ## 102400K .......... .......... .......... .......... .......... 76% 66.4M 2s
    ## 102450K .......... .......... .......... .......... .......... 76%  142M 2s
    ## 102500K .......... .......... .......... .......... .......... 76%  114M 2s
    ## 102550K .......... .......... .......... .......... .......... 76%  119M 2s
    ## 102600K .......... .......... .......... .......... .......... 76% 69.2M 2s
    ## 102650K .......... .......... .......... .......... .......... 76% 74.9M 2s
    ## 102700K .......... .......... .......... .......... .......... 76% 49.0M 2s
    ## 102750K .......... .......... .......... .......... .......... 76% 11.8M 2s
    ## 102800K .......... .......... .......... .......... .......... 76% 30.8M 2s
    ## 102850K .......... .......... .......... .......... .......... 76% 17.8M 2s
    ## 102900K .......... .......... .......... .......... .......... 76% 27.2M 2s
    ## 102950K .......... .......... .......... .......... .......... 76% 20.7M 2s
    ## 103000K .......... .......... .......... .......... .......... 76% 19.8M 2s
    ## 103050K .......... .......... .......... .......... .......... 76% 21.3M 2s
    ## 103100K .......... .......... .......... .......... .......... 76% 32.7M 2s
    ## 103150K .......... .......... .......... .......... .......... 76% 26.9M 2s
    ## 103200K .......... .......... .......... .......... .......... 76% 21.7M 2s
    ## 103250K .......... .......... .......... .......... .......... 76% 68.5M 2s
    ## 103300K .......... .......... .......... .......... .......... 76% 70.0M 2s
    ## 103350K .......... .......... .......... .......... .......... 76% 87.0M 2s
    ## 103400K .......... .......... .......... .......... .......... 76% 30.7M 2s
    ## 103450K .......... .......... .......... .......... .......... 76% 75.4M 2s
    ## 103500K .......... .......... .......... .......... .......... 76% 66.1M 2s
    ## 103550K .......... .......... .......... .......... .......... 76% 66.3M 2s
    ## 103600K .......... .......... .......... .......... .......... 76% 60.5M 2s
    ## 103650K .......... .......... .......... .......... .......... 76% 81.0M 2s
    ## 103700K .......... .......... .......... .......... .......... 77% 79.5M 2s
    ## 103750K .......... .......... .......... .......... .......... 77% 87.0M 2s
    ## 103800K .......... .......... .......... .......... .......... 77% 72.3M 2s
    ## 103850K .......... .......... .......... .......... .......... 77% 61.5M 2s
    ## 103900K .......... .......... .......... .......... .......... 77% 58.7M 2s
    ## 103950K .......... .......... .......... .......... .......... 77% 78.5M 2s
    ## 104000K .......... .......... .......... .......... .......... 77% 73.3M 2s
    ## 104050K .......... .......... .......... .......... .......... 77%  106M 2s
    ## 104100K .......... .......... .......... .......... .......... 77% 73.8M 2s
    ## 104150K .......... .......... .......... .......... .......... 77% 94.9M 2s
    ## 104200K .......... .......... .......... .......... .......... 77% 86.2M 2s
    ## 104250K .......... .......... .......... .......... .......... 77%  104M 2s
    ## 104300K .......... .......... .......... .......... .......... 77% 80.1M 2s
    ## 104350K .......... .......... .......... .......... .......... 77% 88.4M 2s
    ## 104400K .......... .......... .......... .......... .......... 77% 73.6M 2s
    ## 104450K .......... .......... .......... .......... .......... 77% 69.7M 2s
    ## 104500K .......... .......... .......... .......... .......... 77% 90.2M 2s
    ## 104550K .......... .......... .......... .......... .......... 77%  134M 1s
    ## 104600K .......... .......... .......... .......... .......... 77% 75.3M 1s
    ## 104650K .......... .......... .......... .......... .......... 77%  121M 1s
    ## 104700K .......... .......... .......... .......... .......... 77% 80.6M 1s
    ## 104750K .......... .......... .......... .......... .......... 77%  103M 1s
    ## 104800K .......... .......... .......... .......... .......... 77% 88.0M 1s
    ## 104850K .......... .......... .......... .......... .......... 77%  102M 1s
    ## 104900K .......... .......... .......... .......... .......... 77% 86.0M 1s
    ## 104950K .......... .......... .......... .......... .......... 77% 82.7M 1s
    ## 105000K .......... .......... .......... .......... .......... 77% 91.2M 1s
    ## 105050K .......... .......... .......... .......... .......... 78% 69.3M 1s
    ## 105100K .......... .......... .......... .......... .......... 78% 77.2M 1s
    ## 105150K .......... .......... .......... .......... .......... 78%  121M 1s
    ## 105200K .......... .......... .......... .......... .......... 78% 76.1M 1s
    ## 105250K .......... .......... .......... .......... .......... 78%  121M 1s
    ## 105300K .......... .......... .......... .......... .......... 78% 66.5M 1s
    ## 105350K .......... .......... .......... .......... .......... 78% 93.2M 1s
    ## 105400K .......... .......... .......... .......... .......... 78% 98.8M 1s
    ## 105450K .......... .......... .......... .......... .......... 78%  129M 1s
    ## 105500K .......... .......... .......... .......... .......... 78% 74.1M 1s
    ## 105550K .......... .......... .......... .......... .......... 78%  113M 1s
    ## 105600K .......... .......... .......... .......... .......... 78%  109M 1s
    ## 105650K .......... .......... .......... .......... .......... 78%  119M 1s
    ## 105700K .......... .......... .......... .......... .......... 78% 90.5M 1s
    ## 105750K .......... .......... .......... .......... .......... 78%  135M 1s
    ## 105800K .......... .......... .......... .......... .......... 78% 70.1M 1s
    ## 105850K .......... .......... .......... .......... .......... 78% 86.9M 1s
    ## 105900K .......... .......... .......... .......... .......... 78%  108M 1s
    ## 105950K .......... .......... .......... .......... .......... 78%  113M 1s
    ## 106000K .......... .......... .......... .......... .......... 78%  112M 1s
    ## 106050K .......... .......... .......... .......... .......... 78%  134M 1s
    ## 106100K .......... .......... .......... .......... .......... 78%  115M 1s
    ## 106150K .......... .......... .......... .......... .......... 78%  119M 1s
    ## 106200K .......... .......... .......... .......... .......... 78%  104M 1s
    ## 106250K .......... .......... .......... .......... .......... 78% 19.8M 1s
    ## 106300K .......... .......... .......... .......... .......... 78% 88.2M 1s
    ## 106350K .......... .......... .......... .......... .......... 78% 50.9M 1s
    ## 106400K .......... .......... .......... .......... .......... 79% 34.5M 1s
    ## 106450K .......... .......... .......... .......... .......... 79% 78.5M 1s
    ## 106500K .......... .......... .......... .......... .......... 79% 37.2M 1s
    ## 106550K .......... .......... .......... .......... .......... 79% 69.8M 1s
    ## 106600K .......... .......... .......... .......... .......... 79% 23.4M 1s
    ## 106650K .......... .......... .......... .......... .......... 79%  109M 1s
    ## 106700K .......... .......... .......... .......... .......... 79% 99.2M 1s
    ## 106750K .......... .......... .......... .......... .......... 79% 51.5M 1s
    ## 106800K .......... .......... .......... .......... .......... 79% 84.5M 1s
    ## 106850K .......... .......... .......... .......... .......... 79%  101M 1s
    ## 106900K .......... .......... .......... .......... .......... 79% 87.7M 1s
    ## 106950K .......... .......... .......... .......... .......... 79%  121M 1s
    ## 107000K .......... .......... .......... .......... .......... 79% 95.3M 1s
    ## 107050K .......... .......... .......... .......... .......... 79% 14.8M 1s
    ## 107100K .......... .......... .......... .......... .......... 79% 4.21M 1s
    ## 107150K .......... .......... .......... .......... .......... 79% 88.1M 1s
    ## 107200K .......... .......... .......... .......... .......... 79% 53.6M 1s
    ## 107250K .......... .......... .......... .......... .......... 79%  129M 1s
    ## 107300K .......... .......... .......... .......... .......... 79% 5.57M 1s
    ## 107350K .......... .......... .......... .......... .......... 79% 62.8M 1s
    ## 107400K .......... .......... .......... .......... .......... 79% 31.0M 1s
    ## 107450K .......... .......... .......... .......... .......... 79% 89.8M 1s
    ## 107500K .......... .......... .......... .......... .......... 79% 72.3M 1s
    ## 107550K .......... .......... .......... .......... .......... 79%  126M 1s
    ## 107600K .......... .......... .......... .......... .......... 79%  101M 1s
    ## 107650K .......... .......... .......... .......... .......... 79% 49.6M 1s
    ## 107700K .......... .......... .......... .......... .......... 79% 81.6M 1s
    ## 107750K .......... .......... .......... .......... .......... 80% 94.0M 1s
    ## 107800K .......... .......... .......... .......... .......... 80% 71.8M 1s
    ## 107850K .......... .......... .......... .......... .......... 80% 50.7M 1s
    ## 107900K .......... .......... .......... .......... .......... 80% 51.8M 1s
    ## 107950K .......... .......... .......... .......... .......... 80% 19.1M 1s
    ## 108000K .......... .......... .......... .......... .......... 80%  106M 1s
    ## 108050K .......... .......... .......... .......... .......... 80%  113M 1s
    ## 108100K .......... .......... .......... .......... .......... 80% 38.1M 1s
    ## 108150K .......... .......... .......... .......... .......... 80% 43.4M 1s
    ## 108200K .......... .......... .......... .......... .......... 80% 96.0M 1s
    ## 108250K .......... .......... .......... .......... .......... 80%  111M 1s
    ## 108300K .......... .......... .......... .......... .......... 80%  105M 1s
    ## 108350K .......... .......... .......... .......... .......... 80%  129M 1s
    ## 108400K .......... .......... .......... .......... .......... 80%  110M 1s
    ## 108450K .......... .......... .......... .......... .......... 80% 97.4M 1s
    ## 108500K .......... .......... .......... .......... .......... 80%  105M 1s
    ## 108550K .......... .......... .......... .......... .......... 80%  102M 1s
    ## 108600K .......... .......... .......... .......... .......... 80%  109M 1s
    ## 108650K .......... .......... .......... .......... .......... 80%  112M 1s
    ## 108700K .......... .......... .......... .......... .......... 80% 92.6M 1s
    ## 108750K .......... .......... .......... .......... .......... 80% 86.7M 1s
    ## 108800K .......... .......... .......... .......... .......... 80% 40.2M 1s
    ## 108850K .......... .......... .......... .......... .......... 80% 66.2M 1s
    ## 108900K .......... .......... .......... .......... .......... 80% 49.2M 1s
    ## 108950K .......... .......... .......... .......... .......... 80%  116M 1s
    ## 109000K .......... .......... .......... .......... .......... 80% 38.1M 1s
    ## 109050K .......... .......... .......... .......... .......... 80% 97.3M 1s
    ## 109100K .......... .......... .......... .......... .......... 81%  108M 1s
    ## 109150K .......... .......... .......... .......... .......... 81%  117M 1s
    ## 109200K .......... .......... .......... .......... .......... 81%  108M 1s
    ## 109250K .......... .......... .......... .......... .......... 81%  106M 1s
    ## 109300K .......... .......... .......... .......... .......... 81% 92.0M 1s
    ## 109350K .......... .......... .......... .......... .......... 81% 28.8M 1s
    ## 109400K .......... .......... .......... .......... .......... 81%  110M 1s
    ## 109450K .......... .......... .......... .......... .......... 81%  124M 1s
    ## 109500K .......... .......... .......... .......... .......... 81%  112M 1s
    ## 109550K .......... .......... .......... .......... .......... 81%  119M 1s
    ## 109600K .......... .......... .......... .......... .......... 81% 8.41M 1s
    ## 109650K .......... .......... .......... .......... .......... 81% 19.2M 1s
    ## 109700K .......... .......... .......... .......... .......... 81% 40.0M 1s
    ## 109750K .......... .......... .......... .......... .......... 81% 41.1M 1s
    ## 109800K .......... .......... .......... .......... .......... 81% 46.2M 1s
    ## 109850K .......... .......... .......... .......... .......... 81% 93.9M 1s
    ## 109900K .......... .......... .......... .......... .......... 81% 67.0M 1s
    ## 109950K .......... .......... .......... .......... .......... 81%  116M 1s
    ## 110000K .......... .......... .......... .......... .......... 81% 71.4M 1s
    ## 110050K .......... .......... .......... .......... .......... 81% 86.9M 1s
    ## 110100K .......... .......... .......... .......... .......... 81% 83.6M 1s
    ## 110150K .......... .......... .......... .......... .......... 81% 72.4M 1s
    ## 110200K .......... .......... .......... .......... .......... 81% 44.8M 1s
    ## 110250K .......... .......... .......... .......... .......... 81% 37.8M 1s
    ## 110300K .......... .......... .......... .......... .......... 81% 4.19M 1s
    ## 110350K .......... .......... .......... .......... .......... 81%  155M 1s
    ## 110400K .......... .......... .......... .......... .......... 81%  156M 1s
    ## 110450K .......... .......... .......... .......... .......... 82%  173M 1s
    ## 110500K .......... .......... .......... .......... .......... 82%  165M 1s
    ## 110550K .......... .......... .......... .......... .......... 82%  179M 1s
    ## 110600K .......... .......... .......... .......... .......... 82%  168M 1s
    ## 110650K .......... .......... .......... .......... .......... 82%  178M 1s
    ## 110700K .......... .......... .......... .......... .......... 82%  159M 1s
    ## 110750K .......... .......... .......... .......... .......... 82%  169M 1s
    ## 110800K .......... .......... .......... .......... .......... 82%  153M 1s
    ## 110850K .......... .......... .......... .......... .......... 82% 9.17M 1s
    ## 110900K .......... .......... .......... .......... .......... 82% 42.0M 1s
    ## 110950K .......... .......... .......... .......... .......... 82% 60.0M 1s
    ## 111000K .......... .......... .......... .......... .......... 82% 67.1M 1s
    ## 111050K .......... .......... .......... .......... .......... 82% 66.3M 1s
    ## 111100K .......... .......... .......... .......... .......... 82% 43.9M 1s
    ## 111150K .......... .......... .......... .......... .......... 82%  115M 1s
    ## 111200K .......... .......... .......... .......... .......... 82%  136M 1s
    ## 111250K .......... .......... .......... .......... .......... 82%  156M 1s
    ## 111300K .......... .......... .......... .......... .......... 82%  137M 1s
    ## 111350K .......... .......... .......... .......... .......... 82% 16.2M 1s
    ## 111400K .......... .......... .......... .......... .......... 82% 19.5M 1s
    ## 111450K .......... .......... .......... .......... .......... 82%  132M 1s
    ## 111500K .......... .......... .......... .......... .......... 82% 52.5M 1s
    ## 111550K .......... .......... .......... .......... .......... 82% 65.0M 1s
    ## 111600K .......... .......... .......... .......... .......... 82% 65.4M 1s
    ## 111650K .......... .......... .......... .......... .......... 82% 93.2M 1s
    ## 111700K .......... .......... .......... .......... .......... 82% 89.0M 1s
    ## 111750K .......... .......... .......... .......... .......... 82%  140M 1s
    ## 111800K .......... .......... .......... .......... .......... 83%  126M 1s
    ## 111850K .......... .......... .......... .......... .......... 83%  152M 1s
    ## 111900K .......... .......... .......... .......... .......... 83% 33.8M 1s
    ## 111950K .......... .......... .......... .......... .......... 83% 56.7M 1s
    ## 112000K .......... .......... .......... .......... .......... 83% 20.5M 1s
    ## 112050K .......... .......... .......... .......... .......... 83% 61.9M 1s
    ## 112100K .......... .......... .......... .......... .......... 83% 35.8M 1s
    ## 112150K .......... .......... .......... .......... .......... 83% 60.3M 1s
    ## 112200K .......... .......... .......... .......... .......... 83% 48.9M 1s
    ## 112250K .......... .......... .......... .......... .......... 83%  119M 1s
    ## 112300K .......... .......... .......... .......... .......... 83%  120M 1s
    ## 112350K .......... .......... .......... .......... .......... 83% 40.0M 1s
    ## 112400K .......... .......... .......... .......... .......... 83% 98.7M 1s
    ## 112450K .......... .......... .......... .......... .......... 83%  122M 1s
    ## 112500K .......... .......... .......... .......... .......... 83%  106M 1s
    ## 112550K .......... .......... .......... .......... .......... 83%  133M 1s
    ## 112600K .......... .......... .......... .......... .......... 83%  140M 1s
    ## 112650K .......... .......... .......... .......... .......... 83%  165M 1s
    ## 112700K .......... .......... .......... .......... .......... 83% 95.1M 1s
    ## 112750K .......... .......... .......... .......... .......... 83%  112M 1s
    ## 112800K .......... .......... .......... .......... .......... 83% 92.5M 1s
    ## 112850K .......... .......... .......... .......... .......... 83%  136M 1s
    ## 112900K .......... .......... .......... .......... .......... 83%  110M 1s
    ## 112950K .......... .......... .......... .......... .......... 83%  152M 1s
    ## 113000K .......... .......... .......... .......... .......... 83% 39.3M 1s
    ## 113050K .......... .......... .......... .......... .......... 83% 44.9M 1s
    ## 113100K .......... .......... .......... .......... .......... 83% 76.2M 1s
    ## 113150K .......... .......... .......... .......... .......... 84% 68.2M 1s
    ## 113200K .......... .......... .......... .......... .......... 84% 87.4M 1s
    ## 113250K .......... .......... .......... .......... .......... 84%  110M 1s
    ## 113300K .......... .......... .......... .......... .......... 84%  103M 1s
    ## 113350K .......... .......... .......... .......... .......... 84%  119M 1s
    ## 113400K .......... .......... .......... .......... .......... 84% 59.8M 1s
    ## 113450K .......... .......... .......... .......... .......... 84% 92.8M 1s
    ## 113500K .......... .......... .......... .......... .......... 84% 66.6M 1s
    ## 113550K .......... .......... .......... .......... .......... 84% 31.8M 1s
    ## 113600K .......... .......... .......... .......... .......... 84% 68.2M 1s
    ## 113650K .......... .......... .......... .......... .......... 84%  111M 1s
    ## 113700K .......... .......... .......... .......... .......... 84% 85.4M 1s
    ## 113750K .......... .......... .......... .......... .......... 84% 56.8M 1s
    ## 113800K .......... .......... .......... .......... .......... 84% 72.8M 1s
    ## 113850K .......... .......... .......... .......... .......... 84%  140M 1s
    ## 113900K .......... .......... .......... .......... .......... 84% 44.4M 1s
    ## 113950K .......... .......... .......... .......... .......... 84% 92.1M 1s
    ## 114000K .......... .......... .......... .......... .......... 84%  101M 1s
    ## 114050K .......... .......... .......... .......... .......... 84% 53.9M 1s
    ## 114100K .......... .......... .......... .......... .......... 84% 59.2M 1s
    ## 114150K .......... .......... .......... .......... .......... 84% 93.4M 1s
    ## 114200K .......... .......... .......... .......... .......... 84% 96.3M 1s
    ## 114250K .......... .......... .......... .......... .......... 84% 92.2M 1s
    ## 114300K .......... .......... .......... .......... .......... 84% 76.2M 1s
    ## 114350K .......... .......... .......... .......... .......... 84% 99.7M 1s
    ## 114400K .......... .......... .......... .......... .......... 84% 33.5M 1s
    ## 114450K .......... .......... .......... .......... .......... 84%  119M 1s
    ## 114500K .......... .......... .......... .......... .......... 85%  103M 1s
    ## 114550K .......... .......... .......... .......... .......... 85%  165M 1s
    ## 114600K .......... .......... .......... .......... .......... 85% 43.5M 1s
    ## 114650K .......... .......... .......... .......... .......... 85%  105M 1s
    ## 114700K .......... .......... .......... .......... .......... 85% 44.1M 1s
    ## 114750K .......... .......... .......... .......... .......... 85%  115M 1s
    ## 114800K .......... .......... .......... .......... .......... 85%  106M 1s
    ## 114850K .......... .......... .......... .......... .......... 85%  145M 1s
    ## 114900K .......... .......... .......... .......... .......... 85% 30.9M 1s
    ## 114950K .......... .......... .......... .......... .......... 85% 87.1M 1s
    ## 115000K .......... .......... .......... .......... .......... 85%  116M 1s
    ## 115050K .......... .......... .......... .......... .......... 85%  164M 1s
    ## 115100K .......... .......... .......... .......... .......... 85% 79.0M 1s
    ## 115150K .......... .......... .......... .......... .......... 85%  115M 1s
    ## 115200K .......... .......... .......... .......... .......... 85% 51.3M 1s
    ## 115250K .......... .......... .......... .......... .......... 85% 91.3M 1s
    ## 115300K .......... .......... .......... .......... .......... 85% 85.4M 1s
    ## 115350K .......... .......... .......... .......... .......... 85%  168M 1s
    ## 115400K .......... .......... .......... .......... .......... 85% 40.5M 1s
    ## 115450K .......... .......... .......... .......... .......... 85% 97.0M 1s
    ## 115500K .......... .......... .......... .......... .......... 85%  138M 1s
    ## 115550K .......... .......... .......... .......... .......... 85%  135M 1s
    ## 115600K .......... .......... .......... .......... .......... 85% 65.7M 1s
    ## 115650K .......... .......... .......... .......... .......... 85% 65.3M 1s
    ## 115700K .......... .......... .......... .......... .......... 85%  121M 1s
    ## 115750K .......... .......... .......... .......... .......... 85%  157M 1s
    ## 115800K .......... .......... .......... .......... .......... 85% 9.23M 1s
    ## 115850K .......... .......... .......... .......... .......... 86%  128M 1s
    ## 115900K .......... .......... .......... .......... .......... 86%  139M 1s
    ## 115950K .......... .......... .......... .......... .......... 86%  175M 1s
    ## 116000K .......... .......... .......... .......... .......... 86%  138M 1s
    ## 116050K .......... .......... .......... .......... .......... 86%  136M 1s
    ## 116100K .......... .......... .......... .......... .......... 86%  106M 1s
    ## 116150K .......... .......... .......... .......... .......... 86%  118M 1s
    ## 116200K .......... .......... .......... .......... .......... 86%  124M 1s
    ## 116250K .......... .......... .......... .......... .......... 86%  156M 1s
    ## 116300K .......... .......... .......... .......... .......... 86% 10.5M 1s
    ## 116350K .......... .......... .......... .......... .......... 86% 72.4M 1s
    ## 116400K .......... .......... .......... .......... .......... 86% 62.6M 1s
    ## 116450K .......... .......... .......... .......... .......... 86% 57.0M 1s
    ## 116500K .......... .......... .......... .......... .......... 86% 63.9M 1s
    ## 116550K .......... .......... .......... .......... .......... 86%  105M 1s
    ## 116600K .......... .......... .......... .......... .......... 86% 87.6M 1s
    ## 116650K .......... .......... .......... .......... .......... 86% 67.3M 1s
    ## 116700K .......... .......... .......... .......... .......... 86% 83.4M 1s
    ## 116750K .......... .......... .......... .......... .......... 86% 86.3M 1s
    ## 116800K .......... .......... .......... .......... .......... 86% 56.7M 1s
    ## 116850K .......... .......... .......... .......... .......... 86%  109M 1s
    ## 116900K .......... .......... .......... .......... .......... 86% 86.3M 1s
    ## 116950K .......... .......... .......... .......... .......... 86% 82.4M 1s
    ## 117000K .......... .......... .......... .......... .......... 86% 60.0M 1s
    ## 117050K .......... .......... .......... .......... .......... 86% 89.8M 1s
    ## 117100K .......... .......... .......... .......... .......... 86% 67.9M 1s
    ## 117150K .......... .......... .......... .......... .......... 86% 91.1M 1s
    ## 117200K .......... .......... .......... .......... .......... 87% 67.7M 1s
    ## 117250K .......... .......... .......... .......... .......... 87% 93.0M 1s
    ## 117300K .......... .......... .......... .......... .......... 87% 71.5M 1s
    ## 117350K .......... .......... .......... .......... .......... 87% 79.6M 1s
    ## 117400K .......... .......... .......... .......... .......... 87% 88.3M 1s
    ## 117450K .......... .......... .......... .......... .......... 87% 87.9M 1s
    ## 117500K .......... .......... .......... .......... .......... 87% 99.8M 1s
    ## 117550K .......... .......... .......... .......... .......... 87% 89.8M 1s
    ## 117600K .......... .......... .......... .......... .......... 87% 91.1M 1s
    ## 117650K .......... .......... .......... .......... .......... 87%  102M 1s
    ## 117700K .......... .......... .......... .......... .......... 87% 72.6M 1s
    ## 117750K .......... .......... .......... .......... .......... 87% 71.4M 1s
    ## 117800K .......... .......... .......... .......... .......... 87% 68.4M 1s
    ## 117850K .......... .......... .......... .......... .......... 87% 39.0M 1s
    ## 117900K .......... .......... .......... .......... .......... 87% 69.5M 1s
    ## 117950K .......... .......... .......... .......... .......... 87% 74.7M 1s
    ## 118000K .......... .......... .......... .......... .......... 87% 66.1M 1s
    ## 118050K .......... .......... .......... .......... .......... 87% 76.9M 1s
    ## 118100K .......... .......... .......... .......... .......... 87% 74.1M 1s
    ## 118150K .......... .......... .......... .......... .......... 87% 94.4M 1s
    ## 118200K .......... .......... .......... .......... .......... 87% 84.2M 1s
    ## 118250K .......... .......... .......... .......... .......... 87% 13.2M 1s
    ## 118300K .......... .......... .......... .......... .......... 87% 95.7M 1s
    ## 118350K .......... .......... .......... .......... .......... 87%  116M 1s
    ## 118400K .......... .......... .......... .......... .......... 87% 90.2M 1s
    ## 118450K .......... .......... .......... .......... .......... 87%  103M 1s
    ## 118500K .......... .......... .......... .......... .......... 87% 95.8M 1s
    ## 118550K .......... .......... .......... .......... .......... 88%  102M 1s
    ## 118600K .......... .......... .......... .......... .......... 88%  109M 1s
    ## 118650K .......... .......... .......... .......... .......... 88% 91.0M 1s
    ## 118700K .......... .......... .......... .......... .......... 88%  103M 1s
    ## 118750K .......... .......... .......... .......... .......... 88%  121M 1s
    ## 118800K .......... .......... .......... .......... .......... 88%  114M 1s
    ## 118850K .......... .......... .......... .......... .......... 88%  101M 1s
    ## 118900K .......... .......... .......... .......... .......... 88% 14.2M 1s
    ## 118950K .......... .......... .......... .......... .......... 88% 97.3M 1s
    ## 119000K .......... .......... .......... .......... .......... 88% 99.0M 1s
    ## 119050K .......... .......... .......... .......... .......... 88%  110M 1s
    ## 119100K .......... .......... .......... .......... .......... 88%  114M 1s
    ## 119150K .......... .......... .......... .......... .......... 88% 99.5M 1s
    ## 119200K .......... .......... .......... .......... .......... 88% 99.1M 1s
    ## 119250K .......... .......... .......... .......... .......... 88%  119M 1s
    ## 119300K .......... .......... .......... .......... .......... 88% 98.8M 1s
    ## 119350K .......... .......... .......... .......... .......... 88%  108M 1s
    ## 119400K .......... .......... .......... .......... .......... 88%  126M 1s
    ## 119450K .......... .......... .......... .......... .......... 88% 98.3M 1s
    ## 119500K .......... .......... .......... .......... .......... 88%  120M 1s
    ## 119550K .......... .......... .......... .......... .......... 88% 35.3M 1s
    ## 119600K .......... .......... .......... .......... .......... 88% 83.4M 1s
    ## 119650K .......... .......... .......... .......... .......... 88% 81.6M 1s
    ## 119700K .......... .......... .......... .......... .......... 88% 47.8M 1s
    ## 119750K .......... .......... .......... .......... .......... 88% 35.7M 1s
    ## 119800K .......... .......... .......... .......... .......... 88%  125M 1s
    ## 119850K .......... .......... .......... .......... .......... 88%  121M 1s
    ## 119900K .......... .......... .......... .......... .......... 89% 11.5M 1s
    ## 119950K .......... .......... .......... .......... .......... 89% 91.5M 1s
    ## 120000K .......... .......... .......... .......... .......... 89%  114M 1s
    ## 120050K .......... .......... .......... .......... .......... 89%  108M 1s
    ## 120100K .......... .......... .......... .......... .......... 89%  110M 1s
    ## 120150K .......... .......... .......... .......... .......... 89%  129M 1s
    ## 120200K .......... .......... .......... .......... .......... 89%  106M 1s
    ## 120250K .......... .......... .......... .......... .......... 89%  125M 1s
    ## 120300K .......... .......... .......... .......... .......... 89%  121M 1s
    ## 120350K .......... .......... .......... .......... .......... 89% 7.34M 1s
    ## 120400K .......... .......... .......... .......... .......... 89% 38.7M 1s
    ## 120450K .......... .......... .......... .......... .......... 89% 93.5M 1s
    ## 120500K .......... .......... .......... .......... .......... 89% 19.0M 1s
    ## 120550K .......... .......... .......... .......... .......... 89%  117M 1s
    ## 120600K .......... .......... .......... .......... .......... 89% 9.79M 1s
    ## 120650K .......... .......... .......... .......... .......... 89%  104M 1s
    ## 120700K .......... .......... .......... .......... .......... 89%  100M 1s
    ## 120750K .......... .......... .......... .......... .......... 89%  130M 1s
    ## 120800K .......... .......... .......... .......... .......... 89%  105M 1s
    ## 120850K .......... .......... .......... .......... .......... 89% 15.8M 1s
    ## 120900K .......... .......... .......... .......... .......... 89% 93.5M 1s
    ## 120950K .......... .......... .......... .......... .......... 89%  120M 1s
    ## 121000K .......... .......... .......... .......... .......... 89%  101M 1s
    ## 121050K .......... .......... .......... .......... .......... 89%  110M 1s
    ## 121100K .......... .......... .......... .......... .......... 89% 81.8M 1s
    ## 121150K .......... .......... .......... .......... .......... 89% 40.6M 1s
    ## 121200K .......... .......... .......... .......... .......... 89% 75.7M 1s
    ## 121250K .......... .......... .......... .......... .......... 90%  123M 1s
    ## 121300K .......... .......... .......... .......... .......... 90%  108M 1s
    ## 121350K .......... .......... .......... .......... .......... 90% 96.2M 1s
    ## 121400K .......... .......... .......... .......... .......... 90%  106M 1s
    ## 121450K .......... .......... .......... .......... .......... 90%  112M 1s
    ## 121500K .......... .......... .......... .......... .......... 90% 79.3M 1s
    ## 121550K .......... .......... .......... .......... .......... 90% 35.3M 1s
    ## 121600K .......... .......... .......... .......... .......... 90% 75.6M 1s
    ## 121650K .......... .......... .......... .......... .......... 90%  110M 1s
    ## 121700K .......... .......... .......... .......... .......... 90% 66.9M 1s
    ## 121750K .......... .......... .......... .......... .......... 90% 92.3M 1s
    ## 121800K .......... .......... .......... .......... .......... 90% 25.1M 1s
    ## 121850K .......... .......... .......... .......... .......... 90% 87.4M 1s
    ## 121900K .......... .......... .......... .......... .......... 90% 47.2M 1s
    ## 121950K .......... .......... .......... .......... .......... 90% 75.2M 1s
    ## 122000K .......... .......... .......... .......... .......... 90% 53.7M 1s
    ## 122050K .......... .......... .......... .......... .......... 90% 74.2M 1s
    ## 122100K .......... .......... .......... .......... .......... 90% 41.7M 1s
    ## 122150K .......... .......... .......... .......... .......... 90% 44.0M 1s
    ## 122200K .......... .......... .......... .......... .......... 90% 74.6M 1s
    ## 122250K .......... .......... .......... .......... .......... 90% 90.6M 1s
    ## 122300K .......... .......... .......... .......... .......... 90% 34.4M 1s
    ## 122350K .......... .......... .......... .......... .......... 90% 97.4M 1s
    ## 122400K .......... .......... .......... .......... .......... 90% 85.7M 1s
    ## 122450K .......... .......... .......... .......... .......... 90% 69.7M 1s
    ## 122500K .......... .......... .......... .......... .......... 90% 30.0M 1s
    ## 122550K .......... .......... .......... .......... .......... 90% 81.6M 1s
    ## 122600K .......... .......... .......... .......... .......... 91% 77.2M 1s
    ## 122650K .......... .......... .......... .......... .......... 91% 65.1M 1s
    ## 122700K .......... .......... .......... .......... .......... 91%  104M 1s
    ## 122750K .......... .......... .......... .......... .......... 91% 84.1M 1s
    ## 122800K .......... .......... .......... .......... .......... 91% 18.6M 1s
    ## 122850K .......... .......... .......... .......... .......... 91% 79.8M 1s
    ## 122900K .......... .......... .......... .......... .......... 91% 56.3M 1s
    ## 122950K .......... .......... .......... .......... .......... 91%  108M 1s
    ## 123000K .......... .......... .......... .......... .......... 91% 45.8M 1s
    ## 123050K .......... .......... .......... .......... .......... 91% 72.3M 1s
    ## 123100K .......... .......... .......... .......... .......... 91% 64.0M 1s
    ## 123150K .......... .......... .......... .......... .......... 91% 43.1M 1s
    ## 123200K .......... .......... .......... .......... .......... 91%  117M 1s
    ## 123250K .......... .......... .......... .......... .......... 91%  146M 1s
    ## 123300K .......... .......... .......... .......... .......... 91% 74.4M 1s
    ## 123350K .......... .......... .......... .......... .......... 91% 91.2M 1s
    ## 123400K .......... .......... .......... .......... .......... 91%  108M 1s
    ## 123450K .......... .......... .......... .......... .......... 91%  120M 1s
    ## 123500K .......... .......... .......... .......... .......... 91%  101M 1s
    ## 123550K .......... .......... .......... .......... .......... 91%  145M 0s
    ## 123600K .......... .......... .......... .......... .......... 91%  101M 0s
    ## 123650K .......... .......... .......... .......... .......... 91%  108M 0s
    ## 123700K .......... .......... .......... .......... .......... 91%  107M 0s
    ## 123750K .......... .......... .......... .......... .......... 91%  116M 0s
    ## 123800K .......... .......... .......... .......... .......... 91% 28.2M 0s
    ## 123850K .......... .......... .......... .......... .......... 91% 48.1M 0s
    ## 123900K .......... .......... .......... .......... .......... 91%  102M 0s
    ## 123950K .......... .......... .......... .......... .......... 92% 90.0M 0s
    ## 124000K .......... .......... .......... .......... .......... 92% 80.8M 0s
    ## 124050K .......... .......... .......... .......... .......... 92% 96.1M 0s
    ## 124100K .......... .......... .......... .......... .......... 92%  108M 0s
    ## 124150K .......... .......... .......... .......... .......... 92% 33.2M 0s
    ## 124200K .......... .......... .......... .......... .......... 92% 79.2M 0s
    ## 124250K .......... .......... .......... .......... .......... 92%  106M 0s
    ## 124300K .......... .......... .......... .......... .......... 92% 89.8M 0s
    ## 124350K .......... .......... .......... .......... .......... 92%  116M 0s
    ## 124400K .......... .......... .......... .......... .......... 92% 46.7M 0s
    ## 124450K .......... .......... .......... .......... .......... 92% 43.2M 0s
    ## 124500K .......... .......... .......... .......... .......... 92% 50.0M 0s
    ## 124550K .......... .......... .......... .......... .......... 92% 99.1M 0s
    ## 124600K .......... .......... .......... .......... .......... 92% 71.1M 0s
    ## 124650K .......... .......... .......... .......... .......... 92% 50.0M 0s
    ## 124700K .......... .......... .......... .......... .......... 92% 79.3M 0s
    ## 124750K .......... .......... .......... .......... .......... 92%  113M 0s
    ## 124800K .......... .......... .......... .......... .......... 92% 94.1M 0s
    ## 124850K .......... .......... .......... .......... .......... 92%  105M 0s
    ## 124900K .......... .......... .......... .......... .......... 92% 76.2M 0s
    ## 124950K .......... .......... .......... .......... .......... 92% 81.3M 0s
    ## 125000K .......... .......... .......... .......... .......... 92% 98.8M 0s
    ## 125050K .......... .......... .......... .......... .......... 92% 38.7M 0s
    ## 125100K .......... .......... .......... .......... .......... 92%  107M 0s
    ## 125150K .......... .......... .......... .......... .......... 92%  118M 0s
    ## 125200K .......... .......... .......... .......... .......... 92% 26.2M 0s
    ## 125250K .......... .......... .......... .......... .......... 92% 86.2M 0s
    ## 125300K .......... .......... .......... .......... .......... 93% 78.3M 0s
    ## 125350K .......... .......... .......... .......... .......... 93%  140M 0s
    ## 125400K .......... .......... .......... .......... .......... 93%  118M 0s
    ## 125450K .......... .......... .......... .......... .......... 93% 98.1M 0s
    ## 125500K .......... .......... .......... .......... .......... 93%  106M 0s
    ## 125550K .......... .......... .......... .......... .......... 93% 61.8M 0s
    ## 125600K .......... .......... .......... .......... .......... 93% 99.6M 0s
    ## 125650K .......... .......... .......... .......... .......... 93% 10.3M 0s
    ## 125700K .......... .......... .......... .......... .......... 93%  118M 0s
    ## 125750K .......... .......... .......... .......... .......... 93%  148M 0s
    ## 125800K .......... .......... .......... .......... .......... 93% 96.6M 0s
    ## 125850K .......... .......... .......... .......... .......... 93%  111M 0s
    ## 125900K .......... .......... .......... .......... .......... 93%  149M 0s
    ## 125950K .......... .......... .......... .......... .......... 93%  182M 0s
    ## 126000K .......... .......... .......... .......... .......... 93% 47.3M 0s
    ## 126050K .......... .......... .......... .......... .......... 93%  161M 0s
    ## 126100K .......... .......... .......... .......... .......... 93%  147M 0s
    ## 126150K .......... .......... .......... .......... .......... 93%  151M 0s
    ## 126200K .......... .......... .......... .......... .......... 93%  115M 0s
    ## 126250K .......... .......... .......... .......... .......... 93%  187M 0s
    ## 126300K .......... .......... .......... .......... .......... 93% 39.0M 0s
    ## 126350K .......... .......... .......... .......... .......... 93%  139M 0s
    ## 126400K .......... .......... .......... .......... .......... 93% 64.0M 0s
    ## 126450K .......... .......... .......... .......... .......... 93%  129M 0s
    ## 126500K .......... .......... .......... .......... .......... 93%  148M 0s
    ## 126550K .......... .......... .......... .......... .......... 93%  179M 0s
    ## 126600K .......... .......... .......... .......... .......... 93% 29.9M 0s
    ## 126650K .......... .......... .......... .......... .......... 94%  160M 0s
    ## 126700K .......... .......... .......... .......... .......... 94%  145M 0s
    ## 126750K .......... .......... .......... .......... .......... 94% 13.4M 0s
    ## 126800K .......... .......... .......... .......... .......... 94% 84.7M 0s
    ## 126850K .......... .......... .......... .......... .......... 94% 36.4M 0s
    ## 126900K .......... .......... .......... .......... .......... 94% 49.2M 0s
    ## 126950K .......... .......... .......... .......... .......... 94% 79.4M 0s
    ## 127000K .......... .......... .......... .......... .......... 94% 49.0M 0s
    ## 127050K .......... .......... .......... .......... .......... 94%  148M 0s
    ## 127100K .......... .......... .......... .......... .......... 94% 18.7M 0s
    ## 127150K .......... .......... .......... .......... .......... 94%  196M 0s
    ## 127200K .......... .......... .......... .......... .......... 94%  152M 0s
    ## 127250K .......... .......... .......... .......... .......... 94%  174M 0s
    ## 127300K .......... .......... .......... .......... .......... 94%  156M 0s
    ## 127350K .......... .......... .......... .......... .......... 94%  167M 0s
    ## 127400K .......... .......... .......... .......... .......... 94% 7.59M 0s
    ## 127450K .......... .......... .......... .......... .......... 94% 35.7M 0s
    ## 127500K .......... .......... .......... .......... .......... 94% 29.6M 0s
    ## 127550K .......... .......... .......... .......... .......... 94% 32.1M 0s
    ## 127600K .......... .......... .......... .......... .......... 94% 28.3M 0s
    ## 127650K .......... .......... .......... .......... .......... 94% 32.3M 0s
    ## 127700K .......... .......... .......... .......... .......... 94% 33.9M 0s
    ## 127750K .......... .......... .......... .......... .......... 94% 22.9M 0s
    ## 127800K .......... .......... .......... .......... .......... 94% 33.1M 0s
    ## 127850K .......... .......... .......... .......... .......... 94% 34.4M 0s
    ## 127900K .......... .......... .......... .......... .......... 94% 31.5M 0s
    ## 127950K .......... .......... .......... .......... .......... 94% 29.6M 0s
    ## 128000K .......... .......... .......... .......... .......... 95% 30.5M 0s
    ## 128050K .......... .......... .......... .......... .......... 95% 34.0M 0s
    ## 128100K .......... .......... .......... .......... .......... 95% 31.1M 0s
    ## 128150K .......... .......... .......... .......... .......... 95% 33.9M 0s
    ## 128200K .......... .......... .......... .......... .......... 95% 33.1M 0s
    ## 128250K .......... .......... .......... .......... .......... 95% 22.4M 0s
    ## 128300K .......... .......... .......... .......... .......... 95% 36.2M 0s
    ## 128350K .......... .......... .......... .......... .......... 95% 39.2M 0s
    ## 128400K .......... .......... .......... .......... .......... 95% 34.7M 0s
    ## 128450K .......... .......... .......... .......... .......... 95% 36.3M 0s
    ## 128500K .......... .......... .......... .......... .......... 95% 36.8M 0s
    ## 128550K .......... .......... .......... .......... .......... 95% 39.4M 0s
    ## 128600K .......... .......... .......... .......... .......... 95% 23.9M 0s
    ## 128650K .......... .......... .......... .......... .......... 95% 44.7M 0s
    ## 128700K .......... .......... .......... .......... .......... 95% 88.5M 0s
    ## 128750K .......... .......... .......... .......... .......... 95% 92.8M 0s
    ## 128800K .......... .......... .......... .......... .......... 95% 96.9M 0s
    ## 128850K .......... .......... .......... .......... .......... 95% 87.5M 0s
    ## 128900K .......... .......... .......... .......... .......... 95% 47.5M 0s
    ## 128950K .......... .......... .......... .......... .......... 95%  123M 0s
    ## 129000K .......... .......... .......... .......... .......... 95%  101M 0s
    ## 129050K .......... .......... .......... .......... .......... 95% 99.8M 0s
    ## 129100K .......... .......... .......... .......... .......... 95% 97.9M 0s
    ## 129150K .......... .......... .......... .......... .......... 95% 76.6M 0s
    ## 129200K .......... .......... .......... .......... .......... 95%  111M 0s
    ## 129250K .......... .......... .......... .......... .......... 95%  142M 0s
    ## 129300K .......... .......... .......... .......... .......... 95%  108M 0s
    ## 129350K .......... .......... .......... .......... .......... 96% 85.3M 0s
    ## 129400K .......... .......... .......... .......... .......... 96% 97.4M 0s
    ## 129450K .......... .......... .......... .......... .......... 96%  134M 0s
    ## 129500K .......... .......... .......... .......... .......... 96% 95.6M 0s
    ## 129550K .......... .......... .......... .......... .......... 96%  135M 0s
    ## 129600K .......... .......... .......... .......... .......... 96% 83.1M 0s
    ## 129650K .......... .......... .......... .......... .......... 96%  106M 0s
    ## 129700K .......... .......... .......... .......... .......... 96%  113M 0s
    ## 129750K .......... .......... .......... .......... .......... 96%  118M 0s
    ## 129800K .......... .......... .......... .......... .......... 96%  112M 0s
    ## 129850K .......... .......... .......... .......... .......... 96%  109M 0s
    ## 129900K .......... .......... .......... .......... .......... 96% 83.4M 0s
    ## 129950K .......... .......... .......... .......... .......... 96% 80.3M 0s
    ## 130000K .......... .......... .......... .......... .......... 96% 97.8M 0s
    ## 130050K .......... .......... .......... .......... .......... 96%  107M 0s
    ## 130100K .......... .......... .......... .......... .......... 96% 88.2M 0s
    ## 130150K .......... .......... .......... .......... .......... 96%  117M 0s
    ## 130200K .......... .......... .......... .......... .......... 96% 79.3M 0s
    ## 130250K .......... .......... .......... .......... .......... 96%  106M 0s
    ## 130300K .......... .......... .......... .......... .......... 96% 75.0M 0s
    ## 130350K .......... .......... .......... .......... .......... 96%  107M 0s
    ## 130400K .......... .......... .......... .......... .......... 96%  106M 0s
    ## 130450K .......... .......... .......... .......... .......... 96%  152M 0s
    ## 130500K .......... .......... .......... .......... .......... 96%  122M 0s
    ## 130550K .......... .......... .......... .......... .......... 96% 94.1M 0s
    ## 130600K .......... .......... .......... .......... .......... 96% 75.8M 0s
    ## 130650K .......... .......... .......... .......... .......... 97% 85.5M 0s
    ## 130700K .......... .......... .......... .......... .......... 97% 98.8M 0s
    ## 130750K .......... .......... .......... .......... .......... 97%  158M 0s
    ## 130800K .......... .......... .......... .......... .......... 97%  114M 0s
    ## 130850K .......... .......... .......... .......... .......... 97% 92.0M 0s
    ## 130900K .......... .......... .......... .......... .......... 97%  103M 0s
    ## 130950K .......... .......... .......... .......... .......... 97%  104M 0s
    ## 131000K .......... .......... .......... .......... .......... 97% 93.5M 0s
    ## 131050K .......... .......... .......... .......... .......... 97% 91.0M 0s
    ## 131100K .......... .......... .......... .......... .......... 97%  117M 0s
    ## 131150K .......... .......... .......... .......... .......... 97%  125M 0s
    ## 131200K .......... .......... .......... .......... .......... 97% 91.3M 0s
    ## 131250K .......... .......... .......... .......... .......... 97%  118M 0s
    ## 131300K .......... .......... .......... .......... .......... 97%  106M 0s
    ## 131350K .......... .......... .......... .......... .......... 97%  104M 0s
    ## 131400K .......... .......... .......... .......... .......... 97% 94.1M 0s
    ## 131450K .......... .......... .......... .......... .......... 97%  125M 0s
    ## 131500K .......... .......... .......... .......... .......... 97%  103M 0s
    ## 131550K .......... .......... .......... .......... .......... 97%  132M 0s
    ## 131600K .......... .......... .......... .......... .......... 97%  107M 0s
    ## 131650K .......... .......... .......... .......... .......... 97%  158M 0s
    ## 131700K .......... .......... .......... .......... .......... 97% 90.0M 0s
    ## 131750K .......... .......... .......... .......... .......... 97%  109M 0s
    ## 131800K .......... .......... .......... .......... .......... 97%  107M 0s
    ## 131850K .......... .......... .......... .......... .......... 97%  135M 0s
    ## 131900K .......... .......... .......... .......... .......... 97%  108M 0s
    ## 131950K .......... .......... .......... .......... .......... 97%  105M 0s
    ## 132000K .......... .......... .......... .......... .......... 98%  115M 0s
    ## 132050K .......... .......... .......... .......... .......... 98%  109M 0s
    ## 132100K .......... .......... .......... .......... .......... 98% 92.0M 0s
    ## 132150K .......... .......... .......... .......... .......... 98%  144M 0s
    ## 132200K .......... .......... .......... .......... .......... 98%  130M 0s
    ## 132250K .......... .......... .......... .......... .......... 98%  146M 0s
    ## 132300K .......... .......... .......... .......... .......... 98%  106M 0s
    ## 132350K .......... .......... .......... .......... .......... 98%  117M 0s
    ## 132400K .......... .......... .......... .......... .......... 98%  106M 0s
    ## 132450K .......... .......... .......... .......... .......... 98% 93.5M 0s
    ## 132500K .......... .......... .......... .......... .......... 98% 97.4M 0s
    ## 132550K .......... .......... .......... .......... .......... 98%  157M 0s
    ## 132600K .......... .......... .......... .......... .......... 98%  111M 0s
    ## 132650K .......... .......... .......... .......... .......... 98%  114M 0s
    ## 132700K .......... .......... .......... .......... .......... 98%  134M 0s
    ## 132750K .......... .......... .......... .......... .......... 98%  118M 0s
    ## 132800K .......... .......... .......... .......... .......... 98%  121M 0s
    ## 132850K .......... .......... .......... .......... .......... 98%  112M 0s
    ## 132900K .......... .......... .......... .......... .......... 98% 94.1M 0s
    ## 132950K .......... .......... .......... .......... .......... 98%  105M 0s
    ## 133000K .......... .......... .......... .......... .......... 98%  133M 0s
    ## 133050K .......... .......... .......... .......... .......... 98%  147M 0s
    ## 133100K .......... .......... .......... .......... .......... 98%  109M 0s
    ## 133150K .......... .......... .......... .......... .......... 98%  105M 0s
    ## 133200K .......... .......... .......... .......... .......... 98% 99.1M 0s
    ## 133250K .......... .......... .......... .......... .......... 98%  125M 0s
    ## 133300K .......... .......... .......... .......... .......... 98%  101M 0s
    ## 133350K .......... .......... .......... .......... .......... 99%  143M 0s
    ## 133400K .......... .......... .......... .......... .......... 99%  126M 0s
    ## 133450K .......... .......... .......... .......... .......... 99% 78.7M 0s
    ## 133500K .......... .......... .......... .......... .......... 99% 95.9M 0s
    ## 133550K .......... .......... .......... .......... .......... 99% 94.3M 0s
    ## 133600K .......... .......... .......... .......... .......... 99%  106M 0s
    ## 133650K .......... .......... .......... .......... .......... 99% 55.9M 0s
    ## 133700K .......... .......... .......... .......... .......... 99% 72.5M 0s
    ## 133750K .......... .......... .......... .......... .......... 99%  153M 0s
    ## 133800K .......... .......... .......... .......... .......... 99%  107M 0s
    ## 133850K .......... .......... .......... .......... .......... 99%  120M 0s
    ## 133900K .......... .......... .......... .......... .......... 99%  134M 0s
    ## 133950K .......... .......... .......... .......... .......... 99% 43.3M 0s
    ## 134000K .......... .......... .......... .......... .......... 99% 68.9M 0s
    ## 134050K .......... .......... .......... .......... .......... 99%  105M 0s
    ## 134100K .......... .......... .......... .......... .......... 99%  110M 0s
    ## 134150K .......... .......... .......... .......... .......... 99%  149M 0s
    ## 134200K .......... .......... .......... .......... .......... 99% 13.9M 0s
    ## 134250K .......... .......... .......... .......... .......... 99%  141M 0s
    ## 134300K .......... .......... .......... .......... .......... 99%  155M 0s
    ## 134350K .......... .......... .......... .......... .......... 99%  113M 0s
    ## 134400K .......... .......... .......... .......... .......... 99%  144M 0s
    ## 134450K .......... .......... .......... .......... .......... 99%  194M 0s
    ## 134500K .......... .......... .......... .......... .......... 99%  181M 0s
    ## 134550K .......... .......... .......... .......... .......... 99%  175M 0s
    ## 134600K .......... .......... .......... .......... .......... 99%  193M 0s
    ## 134650K .......... .......... .......... .......... .......... 99%  192M 0s
    ## 134700K .......... .......... .......... ..........           100%  172M=5.7s
    ## 
    ## 2020-12-17 13:57:08 (23.1 MB/s) - ‘silva_nr99_v138_train_set.fa.gz.3’ saved [137973851/137973851]

### Commentaire : On a importé les données Silva.

``` r
taxa <- assignTaxonomy(seqtab.nochim, "~/CC2/silva_nr99_v138_train_set.fa.gz", multithread = TRUE)
```

### Commentaire : La fonction assignTaxonomy() met en œuvre l’algorithme du RDP Naive Bayesian Classifier.

``` r
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
```

    ##      Kingdom    Phylum             Class                 Order            
    ## [1,] "Bacteria" "Proteobacteria"   "Alphaproteobacteria" "SAR11 clade"    
    ## [2,] "Bacteria" "Cyanobacteria"    "Cyanobacteriia"      "Synechococcales"
    ## [3,] "Bacteria" "Proteobacteria"   "Alphaproteobacteria" "SAR11 clade"    
    ## [4,] "Bacteria" "Proteobacteria"   "Alphaproteobacteria" "SAR11 clade"    
    ## [5,] "Bacteria" "Proteobacteria"   "Alphaproteobacteria" "SAR11 clade"    
    ## [6,] "Bacteria" "Actinobacteriota" "Acidimicrobiia"      "Actinomarinales"
    ##      Family             Genus                    
    ## [1,] "Clade I"          "Clade Ia"               
    ## [2,] "Cyanobiaceae"     "Synechococcus CC9902"   
    ## [3,] "Clade I"          "Clade Ia"               
    ## [4,] "Clade I"          "Clade Ia"               
    ## [5,] "Clade II"         NA                       
    ## [6,] "Actinomarinaceae" "Candidatus Actinomarina"

### Commentaire : On a effectué l’assignement taxonomique en utilisant les données SILVA. Ensuite, on a examiné les affectations taxonomiques. Par exemple, on peut voir les différentes phyla, les classes, les ordres, les familles, les genres et les espèces.

``` r
save.image(file = "02_stat-analysiscc2-with-DADA2_FinalENV")
```

### Commentaire : La fonction save.image () permet que les objets du fichier 02\_stat-analysiscc2 soient sauvegardés. On va ensuite charger ce fichier de donnée dans un autre fichier 03\_stat-analysiscc2 en utilisant la fonction load() par la suite.
