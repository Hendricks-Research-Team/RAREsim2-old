from raresim.common.sparse import *
from raresim.common.bins import loadBins
from raresim.common.legend import *
from raresim.engine.pruners import *
from typing import Union

class DefaultRunner:
    def __init__(self, runConfig: RunConfig):
        self.args = runConfig.args
        self.runConfig = runConfig
        self.matrix_reader = SparseMatrixReader()
        self.matrix_writer = SparseMatrixWriter()

    def run(self):
        """
        Execute the main process for running the simulation.

        This method performs the following steps:
        1. Loads the necessary data, including the sparse matrix and legend.
        2. Validates that the legend and haplotype file lengths match.
        3. Verifies that input and output legend files are provided.
        4. Retrieves and applies the appropriate transformer based on the run configuration.
        5. Optionally writes the new variant legend if zeroed rows are removed.
        6. Writes the new haplotype file.

        Raises:
            IllegalArgumentException: If input or output legend files are not provided.
        """
        # Start with loading all the necessary data
        matrix: SparseMatrix = self.matrix_reader.loadSparseMatrix(self.args.sparse_matrix)
        legend: Legend = LegendReaderWriter.load_legend(self.args.input_legend)
        
        bins = None
        if self.runConfig.run_type != "probabilistic":
            bins = self.get_bins()

        # Validate inputs
        if legend.row_count() != matrix.num_rows():
            print(f"Legend and Hap file lengths do not match. \n"
                                           f"Legend: {legend.row_count()}, Haps: {matrix.num_rows()}")

        if self.args.input_legend is None:
            raise IllegalArgumentException("Input legend file not provided")

        transformer = self.get_transformer(bins, legend, matrix)
        transformer.transform()

        if self.runConfig.remove_zeroed_rows:
            print()
            print('Writing new variant legend')
            if self.args.output_legend is None:
                raise IllegalArgumentException("Output legend file not provided when remove_zeroed_rows is True")
            LegendReaderWriter.write_legend(legend, self.args.output_legend)

        print()
        print('Writing new haplotype file')
        self.matrix_writer.writeToHapsFile(matrix, self.args.output_hap)

    def get_bins(self) -> list:
        """
        Retrieve the appropriate bin definitions based on the current run configuration.

        Bins are loaded based on the run mode. In 'func_split' mode, separate bins are loaded for functional and synonymous variants.
        In 'syn_only' and 'fun_only' modes, bins are loaded for only synonymous and functional variants, respectively. In all other modes,
        bins are loaded for all variants.

        :return: A dictionary mapping bin IDs to lists of allele counts, or a dictionary of such dictionaries if run mode is 'func_split'
        :rtype: dict
        """
        mode = self.runConfig.run_type
        if mode == "func_split":
            bins = {'fun': loadBins(self.args.exp_fun_bins), 'syn': loadBins(self.args.exp_syn_bins)}
        elif mode == "syn_only":
            bins = loadBins(self.args.syn_bins_only)
        elif mode == "fun_only":
            bins = loadBins(self.args.fun_bins_only)
        else:
            bins = loadBins(self.args.exp_bins)
        return bins

    def get_transformer(self, bins: Union[dict, list], legend, matrix) -> Pruner:
        """
        Retrieve the appropriate transformer based on the current run configuration.

        The transformer is chosen based on the run mode. In 'func_split' mode, a FunctionalSplitPruner is used. In 'syn_only' and 'fun_only' modes,
        a FunctionalSplitPruner is used with the appropriate bins. In all other modes, a StandardPruner is used.

        Parameters
        ----------
        bins : dict
            A dictionary mapping bin IDs to lists of allele counts, or a dictionary of such dictionaries if run mode is 'func_split'
        legend : Legend
            The Legend object describing the input data
        matrix : SparseMatrix
            The SparseMatrix object containing the input data

        Returns
        -------
        Pruner
            The chosen transformer
        """
        mode = self.runConfig.run_type
        print(f"Running with run mode: {mode}")
        if mode == "standard":
            return StandardPruner(self.runConfig, bins, legend, matrix)
        if mode == "func_split":
            return FunctionalSplitPruner(self.runConfig, bins, legend, matrix)
        if mode == "fun_only" or mode == "syn_only":
            return StandardPruner(self.runConfig, bins, legend, matrix)
        if mode == "syn_only":
            return FunctionalSplitPruner(self.runConfig, bins, legend, matrix)
        if mode == "probabilistic":
            return ProbabilisticPruner(self.runConfig, legend, matrix)
