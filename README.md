# Similarity Projection and Related Programs
## Description
This repository contains stable release(s) of the software developed in the course of my PhD candidature, and subsequent revisions.

The programs included in this repository implement algorithms for alignment-free sequence comparison which have been demonstrated to achieve accuracy superior to that of [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi), with execution speed comparable to that of [USEARCH](http://drive5.com).

Work to date primarily focusses on biological sequence data, including protein and RNA, but the approach should be equally applicable to any sequentially structured symbolic data.

The ultimate outcome of the research is a streamlined procedure which projects biological sequences into a metric space using a technique similar to Vector Quantisation to generate a sparse high-dimensional binary feature vector which represents each sequence. A standard metric in the binary embedding space is then used to estimate the evolutionary distance between the original sequences. 

## Publications

* The Similarity Projection concept, algorithm, and preliminary results were published at the 2017 IEEE eScience Conference:<br />
    Buckingham, L., Chappell, T., Hogan, J. M., & Geva, S. (2017, 24-27 Oct. 2017). Similarity Projection: A Geometric Measure for Comparison of Biological Sequences. Paper presented at the 2017 IEEE 13th International Conference on e-Science (e-Science).
* The Compressed Kmer Vocabulary metric embedding method was introduced at the 2018 Australasian Document Computing Symposium<br /> 
     Buckingham, L., Geva, S., & Hogan, J. M. (2018). Protein database search using compressed k-mer vocabularies. Paper presented at the Proceedings of the 23rd Australasian Document Computing Symposium, Dunedin, New Zealand. https://doi.org/10.1145/3291992.3291997

An extended empirical study is under preparation at present; results and associated software will be posted here in due course. This work was severely impacted by workplace disruption associated with COVID-19.

## Work in progress...

I'm still working on this document!
Content will be uploaded in the near future.
