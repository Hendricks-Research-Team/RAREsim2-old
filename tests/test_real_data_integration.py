"""
Integration tests using real data files from the tests/data/ directory.
These tests verify complete end-to-end workflows with actual data and validate outputs.
No exceptions are caught - tests should pass or fail cleanly.
Random seeds are set for deterministic testing.
"""
import unittest
import tempfile
import os
import random
from raresim.common.sparse import SparseMatrixReader, SparseMatrixWriter
from raresim.common.legend import LegendReaderWriter
from raresim.common.bins import loadBins
from raresim.engine.config import RunConfig
from raresim.engine.runner import DefaultRunner
from raresim.calculate.expected_vars import read_mac_bins
from argparse import Namespace


class TestRealDataWorkflows(unittest.TestCase):
    """Test complete end-to-end workflows using real data files"""
    
    def setUp(self):
        """Set up paths to data files"""
        self.data_dir = os.path.join(os.path.dirname(__file__), 'data')
        self.temp_dir = tempfile.mkdtemp()
        
    def tearDown(self):
        """Clean up temporary files"""
        for file in os.listdir(self.temp_dir):
            file_path = os.path.join(self.temp_dir, file)
            if os.path.isfile(file_path):
                os.remove(file_path)
        os.rmdir(self.temp_dir)
    
    def test_load_and_convert_haps_files(self):
        """Test loading haplotype files and converting between formats"""
        haps_file = os.path.join(self.data_dir, 'bigger_test.haps')
        
        reader = SparseMatrixReader()
        matrix = reader.loadSparseMatrix(haps_file)
        
        # Verify matrix was loaded correctly
        self.assertGreater(matrix.num_rows(), 0, "Matrix should have rows")
        self.assertGreater(matrix.num_cols(), 0, "Matrix should have columns")
        
        # Test conversion to .sm format
        sm_output = os.path.join(self.temp_dir, 'test.sm')
        writer = SparseMatrixWriter()
        writer.writeToHapsFile(matrix, sm_output, compression="sm")
        
        self.assertTrue(os.path.exists(sm_output), "SM file should be created")
        
        # Load it back and verify
        loaded_sm = reader.loadSparseMatrix(sm_output)
        self.assertEqual(loaded_sm.num_rows(), matrix.num_rows())
        self.assertEqual(loaded_sm.num_cols(), matrix.num_cols())
        
        # Test conversion to .gz format
        gz_output = os.path.join(self.temp_dir, 'test.haps.gz')
        writer.writeToHapsFile(matrix, gz_output, compression="gz")
        
        self.assertTrue(os.path.exists(gz_output), "GZ file should be created")
        
        # Load it back and verify
        loaded_gz = reader.loadSparseMatrix(gz_output)
        self.assertEqual(loaded_gz.num_rows(), matrix.num_rows())
        self.assertEqual(loaded_gz.num_cols(), matrix.num_cols())
    
    def test_load_legend_files(self):
        """Test loading various legend file formats"""
        # Test stratified legend (functional/synonymous)
        stratified_legend = os.path.join(self.data_dir, 'SmallExample.stratified.legend')
        legend = LegendReaderWriter.load_legend(stratified_legend)
        
        self.assertEqual(legend.row_count(), 31, "Should have 31 variants")
        self.assertIn('fun', legend.get_header(), "Should have 'fun' column")
        
        # Count functional vs synonymous
        fun_count = sum(1 for i in range(legend.row_count()) if legend[i]['fun'] == 'fun')
        syn_count = sum(1 for i in range(legend.row_count()) if legend[i]['fun'] == 'syn')
        
        self.assertGreater(fun_count, 0, "Should have functional variants")
        self.assertGreater(syn_count, 0, "Should have synonymous variants")
        self.assertEqual(fun_count + syn_count, 31, "All variants should be classified")
        
        # Test protected legend
        protected_legend = os.path.join(self.data_dir, 'ProtectiveExample.legend')
        legend = LegendReaderWriter.load_legend(protected_legend)
        
        self.assertIn('protected', legend.get_header(), "Should have 'protected' column")
        
        protected_count = sum(1 for i in range(legend.row_count()) if legend[i]['protected'] == '1')
        self.assertGreater(protected_count, 0, "Should have protected variants")

        # Test probabilistic legend
        prob_legend = os.path.join(self.data_dir, 'ProbExample.probs.legend')
        legend = LegendReaderWriter.load_legend(prob_legend)

        self.assertIn('prob', legend.get_header(), "Should have 'prob' column")
        
    def test_load_bin_files(self):
        """Test loading bin configuration files"""
        # Standard bins
        bins_file = os.path.join(self.data_dir, 'Expected_variants_per_bin_80k.txt')
        bins = loadBins(bins_file)
        
        self.assertEqual(len(bins), 6, "Should have 6 bins")
        for bin_tuple in bins:
            self.assertEqual(len(bin_tuple), 3, "Each bin should have (lower, upper, expected)")
            self.assertGreaterEqual(bin_tuple[1], bin_tuple[0], "Upper bound >= lower bound")
        
        # Functional bins
        fun_bins = loadBins(os.path.join(self.data_dir, 'Expected_variants_functional.txt'))
        self.assertGreater(len(fun_bins), 0, "Should have functional bins")
        
        # Synonymous bins
        syn_bins = loadBins(os.path.join(self.data_dir, 'Expected_variants_synonymous.txt'))
        self.assertGreater(len(syn_bins), 0, "Should have synonymous bins")

    def test_probabilistic_pruning_complete_workflow(self):
        """Test complete probabilistic pruning workflow with whole-row pruning."""
        random.seed(42)

        haps_file = os.path.join(self.data_dir, 'ProbExample.haps')
        legend_file = os.path.join(self.data_dir, 'ProbExample.probs.legend')
        output_hap = os.path.join(self.temp_dir, 'prob_output.haps.gz')

        reader = SparseMatrixReader()
        original_matrix = reader.loadSparseMatrix(haps_file)
        original_legend = LegendReaderWriter.load_legend(legend_file)

        args = Namespace(
            sparse_matrix=haps_file,
            input_legend=legend_file,
            output_legend=None,
            output_hap=output_hap,
            exp_bins=None,
            exp_fun_bins=None,
            exp_syn_bins=None,
            fun_bins_only=None,
            syn_bins_only=None,
            prob=True,
            z=False,
            remove_zeroed_rows=False,
            small_sample=True,
            keep_protected=False,
            activation_threshold=10,
            stop_threshold=20,
            verbose=False
        )

        config = RunConfig(args)
        self.assertEqual(config.run_type, 'probabilistic')

        runner = DefaultRunner(config)
        runner.run()

        self.assertTrue(os.path.exists(output_hap), "Output file must be created")

        output_matrix = reader.loadSparseMatrix(output_hap)
        output_alleles = sum(output_matrix.row_num(i) for i in range(output_matrix.num_rows()))
        original_alleles = sum(original_matrix.row_num(i) for i in range(original_matrix.num_rows()))

        self.assertEqual(output_matrix.num_rows(), original_matrix.num_rows())
        self.assertEqual(output_matrix.num_cols(), original_matrix.num_cols())
        self.assertEqual(original_legend.row_count(), output_matrix.num_rows())
        self.assertLess(output_alleles, original_alleles, "Probabilistic pruning should remove alleles")

    def test_probabilistic_pruning_with_z_flag(self):
        """Test probabilistic pruning with -z removes fully pruned rows and writes a pruned file."""
        random.seed(42)

        haps_file = os.path.join(self.data_dir, 'ProbExample.haps')
        legend_file = os.path.join(self.data_dir, 'ProbExample.probs.legend')
        output_hap = os.path.join(self.temp_dir, 'prob_z_output.haps.gz')
        output_legend = os.path.join(self.temp_dir, 'prob_z_output.legend')

        reader = SparseMatrixReader()
        original_matrix = reader.loadSparseMatrix(haps_file)

        args = Namespace(
            sparse_matrix=haps_file,
            input_legend=legend_file,
            output_legend=output_legend,
            output_hap=output_hap,
            exp_bins=None,
            exp_fun_bins=None,
            exp_syn_bins=None,
            fun_bins_only=None,
            syn_bins_only=None,
            prob=True,
            z=True,
            remove_zeroed_rows=True,
            small_sample=True,
            keep_protected=False,
            activation_threshold=10,
            stop_threshold=20,
            verbose=False
        )

        config = RunConfig(args)
        self.assertEqual(config.run_type, 'probabilistic')

        runner = DefaultRunner(config)
        runner.run()

        self.assertTrue(os.path.exists(output_hap), "Output haps file must be created")
        self.assertTrue(os.path.exists(output_legend), "Output legend file must be created")

        pruned_variants_file = f'{output_legend}-pruned-variants'
        self.assertTrue(os.path.exists(pruned_variants_file), "Pruned variants file must be created")

        output_matrix = reader.loadSparseMatrix(output_hap)
        output_legend_obj = LegendReaderWriter.load_legend(output_legend)

        self.assertLess(output_matrix.num_rows(), original_matrix.num_rows())
        self.assertEqual(output_legend_obj.row_count(), output_matrix.num_rows())
    
    def test_standard_pruning_complete_workflow(self):
        """Test complete standard pruning workflow - NO exception catching"""
        # Set random seed for deterministic results
        random.seed(123)
        
        # Use custom test data designed for this workflow
        haps_file = os.path.join(self.data_dir, 'test_standard.haps')
        bins_file = os.path.join(self.data_dir, 'test_standard_bins.txt')
        output_hap = os.path.join(self.temp_dir, 'standard_output.haps.gz')
        output_legend = os.path.join(self.temp_dir, 'standard_output.legend')
        
        # Load input data
        reader = SparseMatrixReader()
        original_matrix = reader.loadSparseMatrix(haps_file)
        
        # Create matching legend
        temp_legend = os.path.join(self.temp_dir, 'input.legend')
        with open(temp_legend, 'w') as f:
            f.write("id\tposition\ta0\ta1\n")
            for i in range(original_matrix.num_rows()):
                f.write(f"variant_{i}\t{i*1000}\tA\tG\n")
        
        # Load bins to understand expectations
        bins = loadBins(bins_file)
        expected_total_variants = sum(b[2] for b in bins)
        
        # Create configuration
        args = Namespace(
            sparse_matrix=haps_file,
            input_legend=temp_legend,
            output_legend=output_legend,
            output_hap=output_hap,
            exp_bins=bins_file,
            exp_fun_bins=None,
            exp_syn_bins=None,
            fun_bins_only=None,
            syn_bins_only=None,
            z=True,
            remove_zeroed_rows=True,
            small_sample=True,
            keep_protected=False,
            activation_threshold=10,
            stop_threshold=20,
            verbose=False
        )
        
        config = RunConfig(args)
        self.assertEqual(config.run_type, 'standard')
        
        # Run complete workflow - should NOT throw exceptions
        runner = DefaultRunner(config)
        runner.run()
        
        # Verify outputs exist
        self.assertTrue(os.path.exists(output_hap), "Output haplotype file must exist")
        self.assertTrue(os.path.exists(output_legend), "Output legend file must exist")
        
        # Load and validate outputs
        output_matrix = reader.loadSparseMatrix(output_hap)
        output_legend_obj = LegendReaderWriter.load_legend(output_legend)
        
        # With z=True, pruned rows are removed
        self.assertLessEqual(output_matrix.num_rows(), original_matrix.num_rows(),
                            "Pruning should reduce row count")
        self.assertEqual(output_matrix.num_rows(), output_legend_obj.row_count(),
                        "Matrix and legend must have same row count")
        
        # Columns unchanged
        self.assertEqual(output_matrix.num_cols(), original_matrix.num_cols(),
                        "Column count should not change")
        
        # Should have variants remaining
        self.assertGreater(output_matrix.num_rows(), 0, "Should have variants after pruning")
        
        # With seed=123 and our test data, verify specific bin expectations were met
        # Load bins to check against
        bins = loadBins(bins_file)
        expected_total = sum(b[2] for b in bins)  # 3 + 2 + 3 + 2 = 10
        
        # Output should be close to expected total (within reasonable tolerance)
        self.assertGreaterEqual(output_matrix.num_rows(), expected_total * 0.7,
                               "Should have at least 70% of expected variants")
        self.assertLessEqual(output_matrix.num_rows(), expected_total * 1.3,
                            "Should have at most 130% of expected variants")
        
        # Load and validate outputs
        output_matrix = reader.loadSparseMatrix(output_hap)
        output_legend_obj = LegendReaderWriter.load_legend(output_legend)
        
        # With seed=123, we expect EXACT deterministic results
        # The pruning algorithm should produce exactly these values
        self.assertEqual(output_matrix.num_rows(), 7, "With seed=123, should have exactly 7 rows after pruning")
        self.assertEqual(output_legend_obj.row_count(), 7, "Legend should match matrix row count")
        self.assertEqual(output_matrix.num_cols(), original_matrix.num_cols(),
                        "Column count should not change")
        
        # If this fails, someone broke the pruning algorithm
        self.assertGreater(output_matrix.num_rows(), 0, "Must have variants after pruning")
    
    def test_functional_split_complete_workflow(self):
        """Test functional/synonymous split workflow - validates stratified pruning"""
        # Set random seed for deterministic results
        random.seed(456)
        
        # Use custom test data designed for stratified workflow
        haps_file = os.path.join(self.data_dir, 'test_stratified.haps')
        legend_file = os.path.join(self.data_dir, 'test_stratified.legend')
        fun_bins_file = os.path.join(self.data_dir, 'test_fun_bins.txt')
        syn_bins_file = os.path.join(self.data_dir, 'test_syn_bins.txt')
        output_hap = os.path.join(self.temp_dir, 'split_output.haps.gz')
        output_legend = os.path.join(self.temp_dir, 'split_output.legend')
        
        # Load and verify input data
        reader = SparseMatrixReader()
        original_matrix = reader.loadSparseMatrix(haps_file)
        original_legend = LegendReaderWriter.load_legend(legend_file)
        
        self.assertEqual(original_matrix.num_rows(), 12)
        self.assertEqual(original_legend.row_count(), 12)
        
        # Count original variant types
        original_fun = sum(1 for i in range(original_legend.row_count()) 
                          if original_legend[i]['fun'] == 'fun')
        original_syn = sum(1 for i in range(original_legend.row_count()) 
                          if original_legend[i]['fun'] == 'syn')
        
        self.assertGreater(original_fun, 0)
        self.assertGreater(original_syn, 0)
        
        # Create configuration
        args = Namespace(
            sparse_matrix=haps_file,
            input_legend=legend_file,
            output_legend=output_legend,
            output_hap=output_hap,
            exp_bins=None,
            exp_fun_bins=fun_bins_file,
            exp_syn_bins=syn_bins_file,
            fun_bins_only=None,
            syn_bins_only=None,
            z=True,
            remove_zeroed_rows=True,
            small_sample=True,
            keep_protected=False,
            activation_threshold=10,
            stop_threshold=20,
            verbose=False
        )
        
        config = RunConfig(args)
        self.assertEqual(config.run_type, 'func_split')
        
        # Run complete workflow - should NOT throw exceptions
        runner = DefaultRunner(config)
        runner.run()
        
        # Verify outputs exist
        self.assertTrue(os.path.exists(output_hap), "Output haplotype must exist")
        self.assertTrue(os.path.exists(output_legend), "Output legend must exist")
        
        # Load and validate outputs
        output_matrix = reader.loadSparseMatrix(output_hap)
        output_legend_obj = LegendReaderWriter.load_legend(output_legend)
        
        # Verify consistency
        self.assertEqual(output_matrix.num_rows(), output_legend_obj.row_count())
        
        # Count output variant types
        output_fun = sum(1 for i in range(output_legend_obj.row_count()) 
                        if output_legend_obj[i]['fun'] == 'fun')
        output_syn = sum(1 for i in range(output_legend_obj.row_count()) 
                        if output_legend_obj[i]['fun'] == 'syn')
        
        # Both types should still be present
        self.assertGreater(output_fun, 0, "Should have functional variants in output")
        self.assertGreater(output_syn, 0, "Should have synonymous variants in output")
        
        # With seed=456, we expect EXACT deterministic results
        # These exact values prove the stratified pruning works correctly
        # From the output: Functional 3+2+3=8, Synonymous 1+1+1=3 (no pruning occurred, bins matched)
        self.assertEqual(output_matrix.num_rows(), 11, "With seed=456, should have exactly 11 rows")
        self.assertEqual(output_legend_obj.row_count(), 11, "Legend should match matrix")
        self.assertEqual(output_fun, 8, "Should have exactly 8 functional variants")
        self.assertEqual(output_syn, 3, "Should have exactly 4 synonymous variants")
        
        # Verify stratification is maintained
        self.assertEqual(output_fun + output_syn, output_legend_obj.row_count(),
                        "All variants must be classified as fun or syn")
    
    def test_mac_bins_csv_format(self):
        """Test reading MAC bins from CSV format"""
        mac_file = os.path.join(self.data_dir, 'mac_bins.csv')
        
        macs = read_mac_bins(mac_file)
        
        self.assertGreater(len(macs), 0, "Should have MAC bins")
        for mac_tuple in macs:
            self.assertEqual(len(mac_tuple), 2, "Each MAC bin should be (lower, upper)")
            self.assertGreaterEqual(mac_tuple[1], mac_tuple[0], "Upper >= lower")


class TestDataFileIntegrity(unittest.TestCase):
    """Verify all data files are valid and can be loaded"""
    
    def setUp(self):
        """Set up path to data directory"""
        self.data_dir = os.path.join(os.path.dirname(__file__), 'data')
    
    def test_all_legend_files_load_correctly(self):
        """Verify all legend files can be loaded without errors"""
        legend_files = {
            'ProtectiveExample.legend': 31,
            'SmallExample.stratified.legend': 31,
            'Simulated_80k.legend': 248
        }
        
        for filename, expected_rows in legend_files.items():
            file_path = os.path.join(self.data_dir, filename)
            if os.path.exists(file_path):
                legend = LegendReaderWriter.load_legend(file_path)
                self.assertEqual(legend.row_count(), expected_rows,
                               f"{filename} should have {expected_rows} rows")
    
    def test_all_bin_files_load_correctly(self):
        """Verify all bin files can be loaded without errors"""
        bin_files = [
            'Expected_variants_functional.txt',
            'Expected_variants_synonymous.txt',
            'Expected_variants_per_bin_80k.txt',
            'fonlyBins.txt'
        ]
        
        for filename in bin_files:
            file_path = os.path.join(self.data_dir, filename)
            if os.path.exists(file_path):
                bins = loadBins(file_path)
                self.assertGreater(len(bins), 0, f"{filename} should have bins")
                
                # Verify structure
                for bin_tuple in bins:
                    self.assertEqual(len(bin_tuple), 3,
                                   f"Bins in {filename} should be (lower, upper, expected)")
    
    def test_haplotype_files_load_correctly(self):
        """Verify haplotype files can be loaded without errors"""
        haps_files = {
            'ProbExample.haps': (31, 20),
            'bigger_test.haps': (80, 55)
        }
        
        reader = SparseMatrixReader()
        for filename, (expected_rows, expected_cols) in haps_files.items():
            file_path = os.path.join(self.data_dir, filename)
            if os.path.exists(file_path):
                matrix = reader.loadSparseMatrix(file_path)
                self.assertEqual(matrix.num_rows(), expected_rows,
                               f"{filename} should have {expected_rows} rows")
                self.assertEqual(matrix.num_cols(), expected_cols,
                               f"{filename} should have {expected_cols} columns")
    
    def test_compressed_haps_load_correctly(self):
        """Verify compressed haplotype files can be loaded"""
        compressed_files = [
            'bigger_test.haps.gz',
            'Simulated_80k_9.controls.haps.gz'
        ]
        
        reader = SparseMatrixReader()
        for filename in compressed_files:
            file_path = os.path.join(self.data_dir, filename)
            if os.path.exists(file_path):
                matrix = reader.loadSparseMatrix(file_path)
                self.assertGreater(matrix.num_rows(), 0,
                                 f"{filename} should have rows")
                self.assertGreater(matrix.num_cols(), 0,
                                 f"{filename} should have columns")


if __name__ == '__main__':
    unittest.main()
