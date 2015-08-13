
# Hadeda is a tool to identify clusters of similar samples, also known as 'clones'.


## Sound bite version
Clonal strains are identified using a density based clustering algorithm that groups strains that are within 10 SNPs of each other. 

Ester, Martin; Kriegel, Hans-Peter; Sander, Jörg; Xu, Xiaowei (1996). "A density-based algorithm for discovering clusters in large spatial databases with noise". In Simoudis, Evangelos; Han, Jiawei; Fayyad, Usama M. Proceedings of the Second International Conference on Knowledge Discovery and Data Mining (KDD-96). AAAI Press. pp. 226–231. ISBN 1-57735-004-9. CiteSeerX: 10.1.1.71.1980.


## Short version
Clonal strains are identified using a density based clustering algorithm. The algorithms starts by identifying seed clusters. A seed cluster starts by randomly selecting any unclustered strain and counting the number of strains that are with 10 SNPs distance. If this is larger than 3, those starting strains plus its neighbors are considered a core cluster. This core cluster is iteratively expanded using the same procedure. Lastly, when the core no longer grows, all strains that are within SNPs of any single strain within the cluster are added to the cluster. The clonal cluster thus consists of a densely connected core and a looser connected periphery. 


## Longer version
Hadeda uses a density based clustering algorithm that is similar to a hybrid of 'full' and single linkage clustering.  The samples in the cluster core need to be linked tightly together, the periphery can be single linked to the core.

The algorithms has two parameters.
1) E: epsilon: The range to query when expanding a cluster. This is the maximum number of SNPs between two samples for them to be considered linked
2) M: minimum points: The minimum number of samples that need to be within range of a sample to consider that sample a core sample of a cluster.  This value is set to 3 by default.

Linkage is determined using a threshold based approach on the number of SNPs. Samples are linked when they have a pairwise SNP-distance that is smaller than E.

The algorithm processes each sample that is not yet part of a cluster in turn. It looks for samples that are within E SNPs from that sample. The samples that have at least M neighbors are added to the core, the remainder are added to the periphery. Each sample that is not yet examined in the core is then further expanded using the same strategy.

The algorithm makes a final pass and groups all pairs of samples that are within E as clonal pairs. 
 
When E is sufficiently small, the clusters of samples could be considered clones. For higher values, the clusters correspond to spoligotypes, lineages or species

--- 

TL;DR: Hadeda clusters samples into groups that are roughly E SNPs apart, with E being the threshold you see in the figure or file. For X=10, we consider samples within a cluster to be clonal.




