FATCAT Alignment Task
=====================

FATCAT (Flexible structure AlignmenT by Chaining Aligned fragment pairs allowing Twists)

Protein structures are flexible and undergo structural rearrangements as part of their function. FATCAT (Flexible structure AlignmenT by Chaining Aligned fragment pairs allowing Twists) is an approach for flexible protein structure comparison. It simultaneously addresses the two major goals of flexible structure alignment; optimizing the alignment and minimizing the number of rigid-body movements (twists) around pivot points (hinges) introduced in the reference structure. In FATCAT, the structure alignment is formulated as the AFPs (aligned fragment pairs) chaining process allowing at most t twists, and the flexible structure alignment is transformed into a rigid structure alignment when t is forced to be 0. Dynamic programming is used to find the optimal chaining.

More information about FATCAT can be found at the official website: `FATCAT <https://fatcat.godziklab.org/>`_.

.. autofunction:: protein_metamorphisms_is.operations.structural_alignment_tasks.fatcat.align_task

