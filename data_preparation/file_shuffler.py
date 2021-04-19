import os
import shutil

import numpy as np

from dto import FileShufflingConfig


class FileShuffler:
    def __init__(self, config: FileShufflingConfig):
        self._configuration = config

    def run(self):
        directory = self._configuration.directory
        n_permutations = self._configuration.n_permutations

        names = os.listdir(directory)
        if self._configuration.save_permutations:
            save_file = os.path.join(directory, "permutation_list.txt")

        for n in range(n_permutations):
            dest_dir = os.path.join(directory, f'{n}_permutation')
            os.mkdir(dest_dir)
            perm = np.random.permutation(names)

            for src, dest in zip(names, perm):
                shutil.copy(os.path.join(directory, src), os.path.join(dest_dir, dest))

            if self._configuration.save_permutations:
                with open(save_file, "a") as file:
                    for name in perm:
                        file.write(name)
                        file.write("\t")
                    file.write("\n")

