# Copy number variation in the Chimpanzee X chromosome

We want to describe the copy number variation of the genes in the Chimpanzee (Pan Troglodytes) X chromosome (Pan Tro 3).

## Phase 1 - Is there any variation?

By constructing dot-plots of the PanTro3 reference X chromosome in windows of ??Mbp, we surveilled the chromosome for genes containing internal duplications.

Fig1, arbitrary dotplot
![what](https://github.com/cmkobel/ChimpX/blob/master/01dotplot/1plots/window_1.png "example dotplot")


Next thing was to extract all the genes with internal duplications and concatenate them into an *Artificial Chromosome*.

We sourced the web for public Chimpanzee (Pan Troglodytes) genomes (all subspecies). We collected 38 such genomes. (See `subjects/subjects.tsv` for furter information).

Then we mapped the reads from 38 subjects (Chimpanzees) onto the artificial chromosome. By reading off the coverage in each position along each concatenated gene in the artificial chromosome, we were able to describe the relative copy number variation in each subject. We decided that the *DMD* gene is a single copy gene. This means that by dividing the coverage for each gene by the coverage for *DMD* in each subject, we can approximate an absolute copy number for each subject. This data is visualized in the following figure:

Fig2, data phase 1
![2](https://github.com/cmkobel/ChimpX/blob/master/visualization/chimpx_region_points.png "variation in selected genes")


The y-axis represents the Copy Number. The x-axis represents each gene sorted by descending copy number variance.  Each dot is a subject. Colored by sex.

After this initial assessment, we decided to go ahead and look for variation in all ~1000 genes in the X chromosome. One caveat of using the aforementioned dotplot-method, is that it is only possibly to recognize duplications of content within the dotplot window. A quick investigation in the structure of the artificial chromosome highligths many repetitions of sequences which are distantly positioned and unlinked to the gene. This quick investigation is visualized in the following figure:

Fig3, artificial chromosome internal structure
![2](https://github.com/cmkobel/ChimpX/blob/master/visualization/ac3_dotplot_no69_debug.pdf "artificial chromosome internal structure")

## Phase 2 - Including all genes in the X chromosome
Let's see what the standard deviation of coverage is after filtering along the chromosome:


![what](https://github.com/cmkobel/ChimpX/blob/master/visualization/sd_across_30.png?raw=true "standard deviaton of coverage across the Chimpanzee X-chromosome")
*Fig 4: standard deviaton of coverage across the Chimpanzee X-chromosome*


## Phase 3 - Relating to meiotic drive?
