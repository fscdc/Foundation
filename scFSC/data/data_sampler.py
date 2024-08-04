from typing import Iterable, List, Sequence
import numpy as np
import torch
from torch.utils.data import Sampler, SubsetRandomSampler, BatchSampler

class SequentialSubsetSampler(Sampler):
    """Samples elements sequentially from a given list of indices, without replacement.

    Arguments:
        indices (Sequence[int]): A sequence of indices to sample from.
    """

    def __init__(self, indices: Sequence[int]):
        self.indices = indices

    def __iter__(self) -> Iterable[int]:
        return iter(self.indices)

    def __len__(self) -> int:
        return len(self.indices)


class BatchSubsetsSampler(Sampler[List[int]]):
    """Samples batches of indices from multiple subsets. Each subset is sampled 
    either randomly or sequentially, and each batch contains indices from only one subset.

    Arguments:
        subsets (List[Sequence[int]]): A list of subsets of indices.
        batch_size (int): Size of each batch.
        intra_subset_shuffle (bool): Whether to shuffle indices within each subset.
        inter_subset_shuffle (bool): Whether to shuffle the order of subsets.
        drop_last (bool): Whether to drop the last batch if it has fewer than `batch_size` elements.
    """

    def __init__(
        self,
        subsets: List[Sequence[int]],
        batch_size: int,
        intra_subset_shuffle: bool = True,
        inter_subset_shuffle: bool = True,
        drop_last: bool = False,
    ):
        self.subsets = subsets
        self.batch_size = batch_size
        self.intra_subset_shuffle = intra_subset_shuffle
        self.inter_subset_shuffle = inter_subset_shuffle
        self.drop_last = drop_last

        # Create appropriate samplers for each subset
        if intra_subset_shuffle:
            self.subset_samplers = [SubsetRandomSampler(subset) for subset in subsets]
        else:
            self.subset_samplers = [
                SequentialSubsetSampler(subset) for subset in subsets
            ]

        # Create batch samplers for each subset sampler
        self.batch_samplers = [
            BatchSampler(sampler, batch_size, drop_last)
            for sampler in self.subset_samplers
        ]

        if inter_subset_shuffle:
            # Mapping from batch index to batch sampler
            _id_to_batch_sampler = []
            for i, batch_sampler in enumerate(self.batch_samplers):
                _id_to_batch_sampler.extend([i] * len(batch_sampler))
            self._id_to_batch_sampler = np.array(_id_to_batch_sampler)

            assert len(self._id_to_batch_sampler) == len(self)

            self.batch_sampler_iterators = [
                iter(batch_sampler) for batch_sampler in self.batch_samplers
            ]

    def __iter__(self) -> Iterable[List[int]]:
        if self.inter_subset_shuffle:
            # Randomly sample from batch samplers
            random_idx = torch.randperm(len(self._id_to_batch_sampler))
            batch_sampler_ids = self._id_to_batch_sampler[random_idx]
            for batch_sampler_id in batch_sampler_ids:
                batch_sampler_iter = self.batch_sampler_iterators[batch_sampler_id]
                yield next(batch_sampler_iter)
        else:
            for batch_sampler in self.batch_samplers:
                yield from batch_sampler

    def __len__(self) -> int:
        return sum(len(batch_sampler) for batch_sampler in self.batch_samplers)
