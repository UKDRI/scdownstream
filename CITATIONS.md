# UK DRI scdownstream: Citations

If you use this pipeline, please cite the tools it ran for your analysis. The tools are grouped below
by whether they are part of the [curated tool set](docs/usage.md#supported-tool-choices) or are
implemented but still pending curation.

## Frameworks

- [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/)

  > Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.

- [nf-core](https://pubmed.ncbi.nlm.nih.gov/32055031/)

  This pipeline is no longer part of nf-core, but it was created with the nf-core template and
  continues to use nf-core framework code and modules. See
  [`ACKNOWLEDGEMENTS.md`](ACKNOWLEDGEMENTS.md).

  > Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x. PubMed PMID: 32055031.

## Base pipelines

Upstream nf-core/scdownstream, and therefore this fork, drew on the learnings and implementations of:

- [panpipes](https://doi.org/10.1101/2023.03.11.532085)

  > Curion F, Rich-Griffin C, Agarwal D, et al. Panpipes: a pipeline for multiomic single-cell and spatial transcriptomic data analysis. Published online December 18, 2023:2023.03.11.532085. doi:10.1101/2023.03.11.532085

- [scFlow](https://doi.org/10.1101/2021.08.16.456499)

  > Khozoie C, Fancy N, Marjaneh MM, Murphy AE, Matthews PM, Skene N. scFlow: A Scalable and Reproducible Analysis Pipeline for Single-Cell RNA Sequencing Data. Published online August 19, 2021:2021.08.16.456499. doi:10.1101/2021.08.16.456499

- [SIMBA](https://github.com/Mye-InfoBank/SIMBA)

  > Trummer, N. et al. 2024. SIMBA (Single-cell Integration Methods pipeline for Big Atlases)

- [YASCP](https://github.com/wtsi-hgi/yascp)

  > Ozols, M. et al. 2023. YASCP (Yet Another Single Cell Pipeline)

## Curated tool set

### Core data handling and analysis

- [Scanpy](https://pubmed.ncbi.nlm.nih.gov/29409532/)

  > Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. Genome Biol. 2018 Feb 6;19(1):15. doi: 10.1186/s13059-017-1382-4. PubMed PMID: 29409532; PubMed Central PMCID: PMC5802054.

- [AnnData](https://doi.org/10.1101/2021.12.16.473007)

  > Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata: Annotated data. Published online December 19, 2021:2021.12.16.473007. doi:10.1101/2021.12.16.473007

- [anndataR](https://anndatar.data-intuitive.com/)

  > Cannoodt R, Zappia L, Morgan M, Deconinck L (2025). anndataR: AnnData interoperability in R. R package version 0.99.0

- [Seurat](https://pubmed.ncbi.nlm.nih.gov/29608179/) — used for RDS input conversion and the
  SingleCellExperiment output

  > Butler A, Hoffman P, Smibert P, Papalexi E, Satija R. Integrating single-cell transcriptomic data across different conditions, technologies, and species. Nat Biotechnol. 2018 Apr;36(5):411-420. doi: 10.1038/nbt.4096. Epub 2018 Mar 12. PubMed PMID: 29608179; PubMed Central PMCID: PMC5965097.

### Doublet detection

- [Scrublet](https://pubmed.ncbi.nlm.nih.gov/30954476/)

  > Wolock SL, Lopez R, Klein AM. Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data. Cell Syst. 2019 Apr 24;8(4):281-291.e9. doi: 10.1016/j.cels.2018.11.005. PubMed PMID: 30954476; PubMed Central PMCID: PMC6625319.

### Ambient RNA correction

- [decontX](https://pubmed.ncbi.nlm.nih.gov/32138773/)

  > Yang S, Corbett SE, Koga Y, Wang Z, Johnson WE, Yajima M, Campbell JD. Decontamination of ambient RNA in single-cell RNA-seq with DecontX. Genome Biol. 2020 Mar 5;21(1):57. doi: 10.1186/s13059-020-1950-6. PubMed PMID: 32138773; PubMed Central PMCID: PMC7059894.

- [SoupX](https://pubmed.ncbi.nlm.nih.gov/33367645/)

  > Young MD, Behjati S. SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data. Gigascience. 2020 Dec 1;9(12):giaa151. doi: 10.1093/gigascience/giaa151. PubMed PMID: 33367645; PubMed Central PMCID: PMC7763177.

- [CellBender](https://pubmed.ncbi.nlm.nih.gov/37550580/)

  > Fleming SJ, Chaffin MD, Arduini A, Akkad AD, Banks E, Marioni JC, Philippakis AA, Ellinor PT, Babadi M. Unsupervised removal of systematic background noise from droplet-based single-cell experiments using CellBender. Nat Methods. 2023 Sep;20(9):1323-1335. doi: 10.1038/s41592-023-01943-7. PubMed PMID: 37550580.

- [scAR](https://doi.org/10.1101/2022.01.14.476312)

  > Sheng C, Lopes R, Li G, et al. Probabilistic modeling of ambient noise in single-cell omics data. Published online January 17, 2022:2022.01.14.476312. doi:10.1101/2022.01.14.476312

### Integration

- [scVI](https://pubmed.ncbi.nlm.nih.gov/30504886/)

  > Lopez R, Regier J, Cole MB, Jordan MI, Yosef N. Deep generative modeling for single-cell transcriptomics. Nat Methods. 2018 Dec;15(12):1053-1058. doi: 10.1038/s41592-018-0229-2. PubMed PMID: 30504886; PubMed Central PMCID: PMC6289068.

- [scvi-tools](https://pubmed.ncbi.nlm.nih.gov/35132262/)

  > Gayoso A, Lopez R, Xing G, et al. A Python library for probabilistic analysis of single-cell omics data. Nat Biotechnol. 2022 Feb;40(2):163-166. doi: 10.1038/s41587-021-01206-w. PubMed PMID: 35132262.

### Cell type annotation

- [CellTypist](https://pubmed.ncbi.nlm.nih.gov/35549406/)

  > Domínguez Conde C, Xu C, Jarvis LB, et al. Cross-tissue immune cell analysis reveals tissue-specific features in humans. Science. 2022 May 13;376(6594):eabl5197. doi: 10.1126/science.abl5197. PubMed PMID: 35549406; PubMed Central PMCID: PMC9835110.

- [SingleR](https://pubmed.ncbi.nlm.nih.gov/30643263/)

  > Aran D, Looney AP, Liu L, Wu E, Fong V, Hsu A, et al. Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. Nat Immunol. 2019;20(2):163-172. doi: 10.1038/s41590-018-0276-y. Epub 2018 Dec 17. PubMed PMID: 30531964; PubMed Central PMCID: PMC6350770.

- [celldex](https://pubmed.ncbi.nlm.nih.gov/30643263/)

  > Aran D, Looney AP, Liu L, Wu E, Fong V, Hsu A, Chak S, Naikawadi RP, Wolters PJ, Abate AR, Butte AJ, Bhattacharya M (2019). "Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage." Nat. Immunol., 20, 163-172. doi:10.1038/s41590-018-0276-y.

### Clustering and dimensionality reduction

- [Leiden](https://pubmed.ncbi.nlm.nih.gov/30914743/)

  > Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden: guaranteeing well-connected communities. Sci Rep. 2019 Mar 26;9(1):5233. doi: 10.1038/s41598-019-41695-z. PubMed PMID: 30914743; PubMed Central PMCID: PMC6435756.

- [UMAP](https://arxiv.org/abs/1802.03426)

  > McInnes L, Healy J, Melville J. UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. arXiv:1802.03426.

### Cell–cell communication

- [LIANA+](https://pubmed.ncbi.nlm.nih.gov/39223377/)

  > Dimitrov D, Schäfer PSL, Farr E, Rodriguez-Mier P, Lobentanzer S, Badia-i-Mompel P, Dugourd A, Tanevski J, Ramirez Flores RO, Saez-Rodriguez J. LIANA+ provides an all-in-one framework for cell-cell communication inference. Nat Cell Biol. 2024 Sep;26(9):1613-1622. doi: 10.1038/s41556-024-01469-w. PubMed PMID: 39223377.

- [HCOP](https://pubmed.ncbi.nlm.nih.gov/25361968/) — ortholog mapping for non-human data

  > Eyre TA, Wright MW, Lush MJ, Bruford EA. HCOP: a searchable database of human orthology predictions. Brief Bioinform. 2007 Mar;8(1):2-5. doi: 10.1093/bib/bbl030.

### Pseudobulk differential expression

- [decoupler](https://pubmed.ncbi.nlm.nih.gov/36699385/)

  > Badia-i-Mompel P, Vélez Santiago J, Braunger J, Geiss C, Dimitrov D, Müller-Dott S, Taus P, Dugourd A, Holland CH, Ramirez Flores RO, Saez-Rodriguez J. decoupleR: ensemble of computational methods to infer biological activities from omics data. Bioinform Adv. 2022 Sep 30;2(1):vbac016. doi: 10.1093/bioadv/vbac016. PubMed PMID: 36699385.

- [PyDESeq2](https://pubmed.ncbi.nlm.nih.gov/37669147/)

  > Muzellec B, Teleńczuk M, Cabeli V, Andreux M. PyDESeq2: a python package for bulk RNA-seq differential expression analysis. Bioinformatics. 2023 Sep 2;39(9):btad547. doi: 10.1093/bioinformatics/btad547. PubMed PMID: 37669147.

- [DESeq2](https://pubmed.ncbi.nlm.nih.gov/25516281/) — the method PyDESeq2 implements

  > Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 2014;15(12):550. doi: 10.1186/s13059-014-0550-8. PubMed PMID: 25516281; PubMed Central PMCID: PMC4302049.

### Reporting

- [MultiQC](https://pubmed.ncbi.nlm.nih.gov/27312411/)

  > Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8. doi: 10.1093/bioinformatics/btw354. Epub 2016 Jun 16. PubMed PMID: 27312411; PubMed Central PMCID: PMC5039924.

- [Quarto](https://quarto.org/)

  > Allaire JJ, Teague C, Scheidegger C, Xie Y, Dervieux C. Quarto. doi:10.5281/zenodo.5960048

## Tools present in the codebase, pending curation

These are implemented but not yet part of the curated tool set — see
[Supported tool choices](docs/usage.md#supported-tool-choices). They are listed here so the citations
are ready as each is validated and enabled.

- [scANVI](https://pubmed.ncbi.nlm.nih.gov/34310650/)

  > Xu C, Lopez R, Mehlman E, Regier J, Jordan MI, Yosef N. Probabilistic harmonization and annotation of single-cell transcriptomics data with deep generative models. Mol Syst Biol. 2021 Jan;17(1):e9620. doi: 10.15252/msb.20209620. PubMed PMID: 33502086.

- [Harmony](https://pubmed.ncbi.nlm.nih.gov/31740819/)

  > Korsunsky I, Millard N, Fan J, Slowikowski K, Zhang F, Wei K, Baglaenko Y, Brenner M, Loh PR, Raychaudhuri S. Fast, sensitive and accurate integration of single-cell data with Harmony. Nat Methods. 2019 Dec;16(12):1289-1296. doi: 10.1038/s41592-019-0619-0. PubMed PMID: 31740819.

- [BBKNN](https://pubmed.ncbi.nlm.nih.gov/31400197/)

  > Polański K, Young MD, Miao Z, Meyer KB, Teichmann SA, Park JE. BBKNN: fast batch alignment of single cell transcriptomes. Bioinformatics. 2020 Feb 1;36(3):964-965. doi: 10.1093/bioinformatics/btz625. PubMed PMID: 31400197.

- [ComBat](https://pubmed.ncbi.nlm.nih.gov/16632515/)

  > Johnson WE, Li C, Rabinovic A. Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostatistics. 2007 Jan;8(1):118-27. doi: 10.1093/biostatistics/kxj037. PubMed PMID: 16632515.

- [SCimilarity](https://doi.org/10.1101/2023.07.18.549537)

  > Heimberg G, Kuo T, DePianto D, et al. Scalable querying of human cell atlases via a foundational model reveals commonalities across fibrosis-associated macrophages. Published online July 19, 2023:2023.07.18.549537. doi:10.1101/2023.07.18.549537

- [SOLO](https://pubmed.ncbi.nlm.nih.gov/32592658/)

  > Bernstein NJ, Fong NL, Lam I, Roy MA, Hendrickson DG, Kelley DR. Solo: Doublet Identification in Single-Cell RNA-Seq via Semi-Supervised Deep Learning. Cell Syst. 2020 Jul 22;11(1):95-101.e5. doi: 10.1016/j.cels.2020.05.010. PubMed PMID: 32592658.

- [DoubletDetection](https://doi.org/10.5281/zenodo.2678041)

  > Gayoso A, Shor J, Carr AJ, Sharma R, Pe'er D. DoubletDetection. doi:10.5281/zenodo.2678041

- [scds](https://pubmed.ncbi.nlm.nih.gov/31501871/)

  > Bais AS, Kostka D. scds: computational annotation of doublets in single-cell RNA sequencing data. Bioinformatics. 2020 Feb 15;36(4):1150-1158. doi: 10.1093/bioinformatics/btz698. PubMed PMID: 31501871.

## Software packaging and containerisation

- [Anaconda](https://anaconda.com)

  > Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web.

- [Bioconda](https://pubmed.ncbi.nlm.nih.gov/29967506/)

  > Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018 Jul;15(7):475-476. doi: 10.1038/s41592-018-0046-7. PubMed PMID: 29967506.

- [BioContainers](https://pubmed.ncbi.nlm.nih.gov/28379341/)

  > da Veiga Leprevost F, Grüning B, Aflitos SA, Röst HL, Uszkoreit J, Barsnes H, Vaudel M, Moreno P, Gatto L, Weber J, Bai M, Jimenez RC, Sachsenberg T, Pfeuffer J, Alvarez RV, Griss J, Nesvizhskii AI, Perez-Riverol Y. BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics. 2017 Aug 15;33(16):2580-2582. doi: 10.1093/bioinformatics/btx192. PubMed PMID: 28379341; PubMed Central PMCID: PMC5870671.

- [Docker](https://dl.acm.org/doi/10.5555/2600239.2600241)

  > Merkel, D. (2014). Docker: lightweight linux containers for consistent development and deployment. Linux Journal, 2014(239), 2. doi: 10.5555/2600239.2600241.

- [Singularity](https://pubmed.ncbi.nlm.nih.gov/28494014/)

  > Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. eCollection 2017. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.

- [Apptainer](https://apptainer.org/)

  > Apptainer: Application containers for Linux. Apptainer Community.
