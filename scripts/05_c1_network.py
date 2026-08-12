"""Prepare C1 interaction-network node and edge tables for Cytoscape.

Data required
-------------
- sample-specific MAG relative abundance
- sample-specific raw abundance / coverage
- sample-to-habitat mapping
- curated C1 metabolic capacities
- MAG taxonomy

Biological edges require metabolic compatibility and occurrence of both MAGs
within at least one shared habitat. The interaction score weights habitat-level
abundance and sample overlap; it is not a direct estimate of carbon flux.
"""

import numpy as np


def interaction_score(Rs, Rt, Cs, Ct, shared_fraction=None):
    """Calculate the habitat-specific interaction score.

    Rs, Rt: habitat-summed source/target relative abundance
    Cs, Ct: corresponding raw abundance sums
    shared_fraction: fraction of habitat replicates containing both MAGs.
                     If None, both occur in the habitat but not the same sample.
    """
    relative_term = np.log(1 + 1e6 * np.sqrt(Rs * Rt))
    raw_term = np.log(1 + np.sqrt(Cs * Ct))
    overlap_bonus = 0.65 if shared_fraction is None else 1 + shared_fraction
    return relative_term * raw_term * overlap_bonus


def compatible_edges(source_traits, target_traits):
    """Return allowed directed C1 handoffs for a source-target pair."""
    edges = []

    if source_traits.get("methane_oxidizer", False):
        if target_traits.get("methanol_consumer", False):
            edges.append("methane_to_methanol")
        if target_traits.get("formate_consumer", False):
            edges.append("methane_to_formate")

    source_can_release_downstream_c1 = (
        source_traits.get("methanol_consumer", False)
        or source_traits.get("formate_consumer", False)
        or source_traits.get("formaldehyde_processing", False)
    )
    if source_can_release_downstream_c1 and target_traits.get("formate_consumer", False):
        edges.append("methanol_to_formate")

    return edges

# For each source-target-interaction combination observed in more than one
# habitat, retain the habitat producing the highest interaction score.
# Export node and edge tables for Cytoscape. For visualization, use a
# representative subset of high-scoring edges across interaction classes and
# source/target MAGs to reduce visual complexity.
