from Bio.SeqFeature import SeqFeature
import numpy as np
from numpy.typing import NDArray
from utils import argmax


class FeatureAggregator:
    def __init__(self, mode: str) -> None:
        # Variables for ID creation
        self.gene_indices: dict[str, int] = dict()
        self.gene_counter: int = 0

        self.gene_subfeature_ids: dict[int, int] = dict()

        # Variables for ID creation
        self.feature_indices: dict[str, int] = dict()
        self.feature_counter: int = 0

        self.feature_ids: dict[str, int] = (
            dict()
        )  # Used for checking if the feature has already been seen
        self.feature_counters: dict[str, int] = (
            dict()
        )  # Used for assigning incremental IDs

        # Variables for target selection
        self.has_cds: bool = False
        self.cds_lengths: NDArray
        self.exon_lengths: NDArray

        # Define target selecting function by mode
        match mode:
            case "all":
                self.select_targets = self._select_all_targets
            case "cds_longest":
                self.select_targets = self._select_target_with_longest_aggr_cds_or_exon
            case _:
                raise Exception(f"Invalid target selection mode {mode}")

        return None

    def _select_all_targets(self, gene: SeqFeature) -> list[SeqFeature]:
        """
        Return all the sub-features of a gene.
        """
        return gene.sub_features

    def _update_cds_and_exon_info(self, gene: SeqFeature) -> None:
        """
        Update the information needed for selecting features based on aggregate CDS and exon lengths.
        """
        # Reset selection variables
        self.cds_lengths = np.zeros(len(gene.sub_features), dtype=np.int32)
        self.exon_lengths = np.zeros(len(gene.sub_features), dtype=np.int32)
        self.has_cds: bool = False

        for i, child in enumerate(gene.sub_features):
            for grand_child in child.sub_features:
                if grand_child.type == "exon":
                    self.exon_lengths[i] += len(grand_child)
                if grand_child.type == "CDS":
                    self.has_cds = True
                    self.cds_lengths[i] += len(grand_child)

        return None

    def _select_target_with_longest_aggr_cds_or_exon(
        self, gene: SeqFeature
    ) -> list[SeqFeature]:
        """
        Return the gene sub-feature with the longest aggregate CDS if any sub-feature contains a CDS.
        Elsewise, the sub-feature with the longest aggregate exon length.
        """
        self._update_cds_and_exon_info(gene)

        if self.has_cds:
            return [gene.sub_features[argmax(self.cds_lengths)]]
        else:
            return [gene.sub_features[argmax(self.exon_lengths)]]

    def create_aggregated(self, parent: SeqFeature, feature: SeqFeature) -> SeqFeature:
        new_feature = SeqFeature(
            location=feature.location,
            id=f"{feature.type}-aggregate-{parent.id}",
            type=feature.type + "-aggregate",
            qualifiers={"Parent": [parent.id]},
        )
        # The sub_feature attribute must be created manually because it is deprecated in BioPython
        new_feature.sub_features = []

        new_subfeatures = []

        # Copy the entire hierarchy, but this won't work for aggregating features beyond this point
        for child in feature.sub_features:
            new_subfeatures.append(self.create_aggregated(new_feature, child))

        feature.sub_features += new_subfeatures

        return new_feature

    def aggregate_sub_features(self, gene: SeqFeature) -> None:
        if len(gene.sub_features) == 0:
            return None

        if gene.id not in self.gene_indices:
            self.gene_counter += 1
            self.gene_indices[gene.id] = self.gene_counter

        target_transcripts: list[SeqFeature] = self.select_targets(gene)

        for transcript in target_transcripts:
            new_features = []

            for child in transcript.sub_features:
                new_features.append(
                    self.create_aggregated(transcript, child)
                )

            transcript.sub_features += new_features

        return None
