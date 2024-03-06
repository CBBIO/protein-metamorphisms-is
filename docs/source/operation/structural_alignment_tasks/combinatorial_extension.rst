Combinatorial Extension Alignment Task
======================================
Combinatorial Extension (CE) alignment task, a highly efficient algorithm for protein structure alignment. CE is renowned for its ability to rapidly identify alignments between protein structures by extending an initial seed alignment in a combinatorial fashion. This method is particularly effective for comparing protein domains, identifying fold similarities, and analyzing protein structure-function relationships.

The CE algorithm operates under the premise that a good structural alignment can be constructed by piecing together short segments of alignments, known as aligned fragment pairs (AFPs), that are identified across the entire length of the protein structures. This combinatorial approach to extending alignments allows CE to uncover optimal alignments that highlight both the similarities and differences in the 3D structures of proteins, even in the presence of significant structural variations.

More information about CE-Align can be found at the official website: `CE-Align <https://pymolwiki.org/index.php/Cealign/>`_.

.. autofunction:: protein_metamorphisms_is.operations.structural_alignment_tasks.combinatorial_extension.align_task
