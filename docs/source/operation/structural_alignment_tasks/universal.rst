Universal Structural Alignment Task
===================================


US-align (Universal Structural alignment) is a unified protocol to compare 3D structures of different macromolecules (proteins, RNAs and DNAs) in different forms (monomers, oligomers and heterocomplexes) for both pairwise and multiple structure alignments. The core alogrithm of US-align is extended from TM-align and generates optimal structural alignments by maximizing TM-score of compared strucures through heuristic dynamic programming iterations. Large-scale benchmark tests showed that US-align can generate more accurate structural alignments with significantly reduced CPU time, compared to the state-of-the-art methods developed for specific structural alignment tasks. TM-score has values in (0,1] with 1 indicating an identical structure match, where a TM-score ≥0.5 (or 0.45) means the structures share the same global topology for proteins (or RNAs).

More information about US-align can be found at the official website: `US-align <https://zhanggroup.org/US-align/>`_.

.. autofunction:: protein_metamorphisms_is.operations.structural_alignment_tasks.universal.align_task
