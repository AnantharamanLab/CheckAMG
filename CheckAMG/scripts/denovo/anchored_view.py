import pst
import torch


class AnchoredView:
    def __init__(self, anchor_idx: int, num_proteins: int, context_size: int):
        self.anchor_idx = anchor_idx

        # TODO: this might pose a problem depending on the fragmentation pattern (if any...)
        if anchor_idx >= num_proteins:
            raise ValueError(
                f"The anchor_idx {anchor_idx} must be less than the number of proteins {num_proteins}"
            )

        self.num_proteins = num_proteins
        self.context_size = min(context_size, num_proteins)

        self._construct_windows()

    def _construct_windows(self):
        # these are all possible self.context_size+1 length contiguous runs of proteins that contain the anchor
        # the anchor should be in each possible position in windows of that size

        self._windows: list[range] = []
        for protein_idx in range(self.num_proteins):
            lookahead = (
                min(self.num_proteins - 1, protein_idx + self.context_size) + 1
            )
            context = range(protein_idx, lookahead)

            if self.anchor_idx in context:
                self._windows.append(context)

    def all_anchored_views(self) -> list[torch.Tensor]:
        base_edge_index = pst.GenomeGraph.create_fully_connected_graph(
            num_nodes=self.context_size + 1
        )

        edge_indices: list[torch.Tensor] = []
        for window in self._windows:
            edge_index = base_edge_index + window.start

            # cap max ptn to be a real ptn from this genome
            valid_mask = (edge_index < self.num_proteins).all(0)
            edge_index = edge_index[:, valid_mask]

            # ignore single ptn views
            if (
                edge_index.shape[-1] > 1
                and len(edge_index.unique()) == self.num_nodes_per_view
            ):
                edge_indices.append(edge_index)

        return edge_indices

    @property
    def num_nodes_per_view(self) -> int:
        # all subgraphs should have self.context_size + 1 nodes
        # self.context_size for the neighbors + 1 for the anchor
        return self.context_size + 1

    @property
    def total_views(self) -> int:
        # context_size = 4
        # num_proteins = 10

        # anchor = 0
        # views: 0-4

        # anchor = 5
        # views 0-5, 1-6, 2-7, 3-8, 4-9, 5-10
        num_views = self.context_size + 1  # max number of views to start with

        lookbehind = self.anchor_idx - self.context_size
        if lookbehind < 0:
            # remove the number of leftward views since they aren't real
            num_views += lookbehind

        lookahead = self.anchor_idx + self.context_size - self.num_proteins + 1
        if lookahead > 0:
            # remove rightward views
            num_views -= lookahead

        return num_views

    @classmethod
    def from_genome_graph(
        cls,
        graph: pst.GenomeGraph,
        anchor_idx: int,
        context_size: int,
    ):
        return cls(
            anchor_idx=anchor_idx,
            num_proteins=int(graph.num_proteins),
            context_size=context_size,
        )
