# Comparing DESeq2 and edgeR for Gene Expression Analysis: My Experience and Insights

  Recently, I had the opportunity to experiment with two of the most popular R libraries for gene expression analysis: DESeq2 and edgeR. Both tools are widely recognized for their ability to identify differentially expressed genes (DEGs) from RNA sequencing data. Despite their common goal, their underlying methods and workflows differ, and each may be better suited to specific analysis scenarios.
  
## Similarities Between DESeq2 and edgeR
  At their core, both DESeq2 and edgeR take a raw count matrix as input and output key statistics such as log2 fold changes, p-values, and adjusted p-values in the form of False Discovery Rate (FDR). They model RNA-seq data using the negative binomial distribution, which effectively captures both biological variability and technical noise inherent in sequencing experiments.
  However, while they share this statistical foundation, the tools diverge notably in how they approach normalization and dispersion estimation, two critical steps in differential expression analysis.

## Differences in Normalization and Dispersion Estimation

  DESeq2 employs the median-of-ratios normalization method. This approach accounts for differences in sequencing depth and adjusts for RNA composition bias across samples. It is especially robust against outliers and compositional skew.   A caveat, however, is that DESeq2 assumes that most genes are not differentially expressed, which can limit its effectiveness if this assumption is violated.
  On the other hand, edgeR uses the Trimmed Mean of M-values (TMM) normalization method. TMM tends to perform better on datasets where many genes are differentially expressed or where expression changes are asymmetric. This method can be more sensitive to extreme counts if the trimming is not carefully applied, but it offers a powerful way to handle strong global shifts in gene expression.
  Regarding dispersion estimation, DESeq2 uses a shrinkage approach, which borrows information across all genes to estimate dispersions more reliably, particularly beneficial for experiments with small sample sizes. It applies an empirical Bayes shrinkage not only to dispersion estimates but also to fold changes, improving stability and interpretability.
  EdgeR also utilizes empirical Bayes techniques but offers more flexibility. It provides options to estimate dispersions as tagwise, common, or trended, making it well-suited for larger sample sizes or more complex experimental designs.

## Statistical Testing Approaches
  When it comes to hypothesis testing, both tools support likelihood ratio tests (LRT) and quasi-likelihood F-tests. DESeq2 typically relies on Wald tests for straightforward pairwise comparisons, providing a simpler framework for many typical analyses. In contrast, edgeR favors quasi-likelihood F-tests (QLF), which tend to offer improved control over false positive rates, particularly valuable in more complex or noisy datasets.
  
## Usability and Learning Curve
  From a usability standpoint, DESeq2 generally offers a more beginner-friendly experience with high-level wrapper functions that abstract many complexities. EdgeR, conversely, presents a steeper learning curve, requiring more manual configuration and understanding of its detailed options, but rewarding advanced users with greater flexibility and control.
  Both packages are well-documented and supported by active communities, so users at any level can find guidance and examples to help them.
  
## Choosing Between DESeq2 and edgeR: When to Use Which
  Based on the literature, certain scenarios lend themselves better to one tool over the other:
For small sample sizes (e.g., three replicates per group), DESeq2’s shrinkage methods provide greater stability and reliable results.
For large and complex experimental designs such as time-course studies or multi-factor experiments, edgeR’s flexibility in model specification makes it the preferred choice.
When dealing with datasets featuring many differentially expressed genes or global expression shifts, edgeR’s TMM normalization offers advantages.
For simpler pairwise comparisons or straightforward differential expression tasks, DESeq2 is easier to implement and interpret.
When you need complex contrast matrices or fine-grained control over the model, edgeR shines.

## Example Code

library(DESeq2)
Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)  # change 'condition' to your factor

Pre-filtering: remove genes with low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

Run DESeq pipeline: normalization, dispersion estimation, model fitting
dds <- DESeq(dds)

Results extraction (default: condition B vs condition A)
res <- results(dds)

View top differentially expressed genes
head(res[order(res$padj), ])

Plot MA-plot
plotMA(res, main="DESeq2", ylim=c(-5,5))

library(edgeR)
Create DGEList object
dge <- DGEList(counts=counts, group=coldata$condition)

Filter lowly expressed genes
keep <- filterByExpr(dge)
dge <- dge[keep,, keep.lib.sizes=FALSE]

Normalize with TMM
dge <- calcNormFactors(dge)

Estimate dispersion
dge <- estimateDisp(dge)

Fit negative binomial GLM
fit <- glmQLFit(dge, design=model.matrix(~ coldata$condition))

Perform quasi-likelihood F-test for differential expression
qlf <- glmQLFTest(fit, coef=2)  # coef=2 compares second level vs first level

Extract top DE genes
topTags(qlf)

Plot MD plot (similar to MA-plot)
plotMD(qlf)

## Normalization Results

### Before
!before

### EdgeR
!edger

### DESeq2
!before


## Summary and Final Thoughts

  In summary, DESeq2 and edgeR each have strengths tailored to different experimental setups. DESeq2 is excellent for smaller, simpler analyses and provides an accessible entry point for newcomers. EdgeR is better suited for larger, more complicated datasets and users who want granular control.
In practice, many researchers run both tools and compare results to gain confidence in their findings. Ultimately, the choice should be guided by your experimental design, sample size, and specific biological questions.
If you are just getting started with RNA-seq differential expression,it is recommended beginning with DESeq2. However, if your study involves multi-factor designs or you want to maximize flexibility, then edgeR is a powerful tool to master. Both tools are battle-tested and supported by extensive documentation and user communities, so you can’t go wrong with either. The best approach is to understand their differences and choose the one that aligns with your project's needs.
