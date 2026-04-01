import unittest
import random
import tempfile
import io
from argparse import Namespace
import os
from contextlib import redirect_stdout

from raresim.common.legend import Legend
from raresim.common.sparse import SparseMatrix
from raresim.engine.config import RunConfig
from raresim.engine.pruners import StandardPruner, FunctionalSplitPruner, ProbabilisticPruner


class TestStandardPruner(unittest.TestCase):
    def setUp(self):
        """Set up test data"""
        random.seed(123)
        header = ["id", "position", "a0", "a1", "protected", "fun"]
        legend = Legend(header)
        matrix = SparseMatrix()
        
        # Select protected variants: 5 from each bin
        # Bin 0 (AC=1): rows 5-39
        # Bin 1 (AC=2-5): rows 40-59
        # Bin 2 (AC=6-10): rows 60-79
        # Bin 3 (AC>10): rows 80-99
        self.protected_indices = [
            10, 15, 20, 25, 30,  # 5 from bin 0 (AC=1)
            42, 45, 48, 52, 55,  # 5 from bin 1 (AC=2-5)
            62, 65, 68, 72, 75,  # 5 from bin 2 (AC=6-10)
            82, 85, 88, 92, 95   # 5 from bin 3 (AC>10)
        ]
        
        for i in range(100):
            protected = "1" if i in self.protected_indices else "0"
            fun = "fun" if i % 2 == 0 else "syn"  # Alternate between fun and syn
            legend.add_row([f"allele{i}", f"{i}", "C", "G", protected, fun])
            if i < 5:
                matrix.add_row([])
            if 5 <= i < 40:
                columns = random.sample(range(30), 1)
                matrix.add_row(sorted(columns))
            if 40 <= i < 60:
                num_columns = random.randrange(2, 5)
                columns = random.sample(range(30), num_columns)
                matrix.add_row((sorted(columns)))
            if 60 <= i < 80:
                num_columns = random.randrange(6, 10)
                columns = random.sample(range(30), num_columns)
                matrix.add_row((sorted(columns)))
            if 80 <= i:
                num_columns = random.randrange(11, 30)
                columns = random.sample(range(30), num_columns)
                matrix.add_row((sorted(columns)))

        self.legend = legend
        self.matrix = matrix
        self.bins = [
            [1,1,25],
            [2,5,10],
            [6,10,3]
        ]
        self.temp_dir = tempfile.mkdtemp()
        self.legend_file = os.path.join(self.temp_dir, "legend")
    def tearDown(self):
        """Clean up temporary files"""
        pruned_file = f'{self.legend_file}-pruned-variants'
        if os.path.exists(pruned_file):
            os.remove(pruned_file)
        if os.path.exists(self.legend_file):
            os.remove(self.legend_file)
        os.rmdir(self.temp_dir)
        
    def test_assign_bins(self):
        args = Namespace(z=False, exp_bins='some_bins', small_sample=True, activation_threshold=5, stop_threshold=2, keep_protected=False, output_legend=self.legend_file, input_legend=self.legend_file)
        run_config = RunConfig(args)
        assert(run_config.run_type == "standard")
        
        pruner = StandardPruner(run_config, self.bins, self.legend, self.matrix)
        assigned_bins = pruner.assign_bins()
        
        assert(len(assigned_bins) == 4)
        assert(len(assigned_bins[0]) == 35)
        assert(len(assigned_bins[1]) == 20)
        assert(len(assigned_bins[2]) == 20)
        assert(len(assigned_bins[3]) == 20)
        
    def test_full_transform_no_z(self):
        args = Namespace(z=False, exp_bins='some_bins', small_sample=True, activation_threshold=5, stop_threshold=2, keep_protected=False, output_legend=self.legend_file, input_legend=self.legend_file)
        runConfig = RunConfig(args)
        assert(runConfig.run_type == "standard")
        
        pruner = StandardPruner(runConfig, self.bins, self.legend, self.matrix)
        pruner.transform()
        
        after_assigned_bins = pruner.assign_bins()
        
        assert(abs(len(after_assigned_bins[0]) - self.bins[0][2]) <= 5)
        assert(abs(len(after_assigned_bins[1]) - self.bins[1][2]) <= 5)
        assert(abs(len(after_assigned_bins[2]) - self.bins[2][2]) <= 5)
        assert(len(after_assigned_bins[3]) == 20)
        
        # z flag is false. Assert that the number of rows in the matrix and legend is 100 still.
        assert(self.matrix.num_rows() == 100)
        assert(self.legend.row_count() == 100)
        
        # z flag is false. Count number of rows with 0 columns and assert that when added with the length of all
        # sub-lists in after_assigned_bins is equal to 100.
        num_rows_with_0_columns = 0
        for i in range(self.matrix.num_rows()):
            row_count = self.matrix.row_num(i)
            if row_count == 0:
                num_rows_with_0_columns += 1
        assert(num_rows_with_0_columns + sum(len(after_assigned_bins[entry]) for entry in after_assigned_bins) == 100)
        
        # Verify pruned variants file was created and contains expected data
        pruned_file = f'{self.legend_file}-pruned-variants'
        assert(os.path.exists(pruned_file))
        
        with open(pruned_file, 'r') as f:
            lines = f.readlines()
            # First line is header
            assert(lines[0].strip() == "\t".join(self.legend.get_header()))
            # Number of pruned variants should equal num_rows_with_0_columns
            assert(len(lines) - 1 == num_rows_with_0_columns)
    
    def test_full_transform_with_z(self):
        args = Namespace(z=True, exp_bins='some_bins', small_sample=True, activation_threshold=5, stop_threshold=2, keep_protected=False, output_legend=self.legend_file, input_legend=self.legend_file)
        runConfig = RunConfig(args)
        assert(runConfig.run_type == "standard")
        
        pruner = StandardPruner(runConfig, self.bins, self.legend, self.matrix)
        pruner.transform()
        
        after_assigned_bins = pruner.assign_bins()
        
        assert(abs(len(after_assigned_bins[0]) - self.bins[0][2]) <= 5)
        assert(abs(len(after_assigned_bins[1]) - self.bins[1][2]) <= 5)
        assert(abs(len(after_assigned_bins[2]) - self.bins[2][2]) <= 5)
        assert(len(after_assigned_bins[3]) == 20)
        
        # z flag is true. Assert that the number of rows in the matrix and legend has been reduced.
        # The total should equal the sum of all kept variants in bins.
        total_kept_variants = sum(len(after_assigned_bins[entry]) for entry in after_assigned_bins)
        assert(self.matrix.num_rows() == total_kept_variants)
        assert(self.legend.row_count() == total_kept_variants)
        
        # Verify pruned variants file was created and contains expected data
        pruned_file = f'{self.legend_file}-pruned-variants'
        assert(os.path.exists(pruned_file))
        
        with open(pruned_file, 'r') as f:
            lines = f.readlines()
            # First line is header
            assert(lines[0].strip() == "\t".join(self.legend.get_header()))
            # Number of pruned variants should equal 100 - total_kept_variants
            assert(len(lines) - 1 == 100 - total_kept_variants)
    
    def test_full_transform_with_keep_protected(self):
        args = Namespace(z=False, exp_bins='some_bins', small_sample=True, activation_threshold=5, stop_threshold=2, keep_protected=True, output_legend=self.legend_file, input_legend=self.legend_file, verbose=True)
        runConfig = RunConfig(args)
        assert(runConfig.run_type == "standard")
        
        pruner = StandardPruner(runConfig, self.bins, self.legend, self.matrix)
        pruner.transform()
        
        after_assigned_bins = pruner.assign_bins()
        
        # Verify all protected variants are still in the bins (not pruned)
        all_kept_rows = []
        for bin_id in after_assigned_bins:
            all_kept_rows.extend(after_assigned_bins[bin_id])
        
        for protected_idx in self.protected_indices:
            if self.matrix.row_num(protected_idx) > 0:  # Only check if variant has AC > 0
                assert(protected_idx in all_kept_rows), f"Protected variant {protected_idx} was pruned!"
        
        # Verify pruned variants file doesn't contain protected variants
        pruned_file = f'{self.legend_file}-pruned-variants'
        assert(os.path.exists(pruned_file))
        
        with open(pruned_file, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:  # Skip header
                variant_id = line.split('\t')[0]
                # Extract index from variant_id (e.g., "allele10" -> 10)
                idx = int(variant_id.replace('allele', ''))
                assert(idx not in self.protected_indices), f"Protected variant {idx} was written to pruned file!"
    
    def test_borrowing_from_extra_rows(self):
        """Test that smaller bins can borrow from extra_rows (variants pruned from larger bins)"""
        # Processing order: Bin 2 → Bin 1 → Bin 0 (reverse order)
        # Bin 2 (processed first) will be pruned down, creating extra_rows with AC=6-10
        # Bin 1 (processed second) needs more than it has, so it borrows from extra_rows
        # The borrowed variants (originally AC=6-10) will be pruned down to AC=2-5
        bins_with_borrowing = [
            [1, 1, 25],   # Bin 0: AC=1, need 25 (have 35, processed last)
            [2, 5, 30],   # Bin 1: AC=2-5, need 30 (have 20, need to borrow ~10 from extra_rows)
            [6, 10, 3]    # Bin 2: AC=6-10, need 3 (have 20, will prune ~17 into extra_rows, processed first)
        ]
        
        args = Namespace(z=False, exp_bins='some_bins', small_sample=True, activation_threshold=5, stop_threshold=2, keep_protected=False, output_legend=self.legend_file, input_legend=self.legend_file, verbose=False)
        runConfig = RunConfig(args)
        assert(runConfig.run_type == "standard")
        
        pruner = StandardPruner(runConfig, bins_with_borrowing, self.legend, self.matrix)
        pruner.transform()
        
        after_assigned_bins = pruner.assign_bins()
        
        # Bin 1 should have close to 30 variants (borrowed from extra_rows created by bin 2)
        assert(abs(len(after_assigned_bins[1]) - 30) <= 5), f"Bin 1 should have ~30 variants but has {len(after_assigned_bins[1])}"
        
        # Verify that all variants in bin 1 have AC in the range [2, 5]
        # (borrowed variants from bin 2 were pruned down to fit this bin)
        for row_id in after_assigned_bins[1]:
            row_ac = self.matrix.row_num(row_id)
            assert(2 <= row_ac <= 5), f"Variant {row_id} in bin 1 has AC={row_ac}, expected AC in [2, 5]"
    
    def test_borrowing_from_reserve_pool(self):
        """Test that bins can borrow from reserve_pool when extra_rows is insufficient"""
        # Processing order: Bin 2 → Bin 1 → Bin 0 (reverse order)
        # Bin 2 (processed first) doesn't prune much, so extra_rows will be small
        # Bin 1 (processed second) needs many more variants than extra_rows can provide
        # So it must borrow from reserve_pool (bin 3 with AC > 10)
        bins_with_reserve_borrowing = [
            [1, 1, 25],   # Bin 0: AC=1, need 25 (have 35, processed last)
            [2, 5, 30],   # Bin 1: AC=2-5, need 30 (have 20, need to borrow ~10)
            [6, 10, 18]   # Bin 2: AC=6-10, need 18 (have 20, will prune only ~2 into extra_rows, processed first)
        ]
        # Bin 2 prunes ~2 into extra_rows, but bin 1 needs ~10 more
        # So bin 1 must borrow ~8 from reserve_pool (bin 3 with AC > 10, which has 20 variants)
        
        args = Namespace(z=False, exp_bins='some_bins', small_sample=True, activation_threshold=5, stop_threshold=2, keep_protected=False, output_legend=self.legend_file, input_legend=self.legend_file, verbose=False)
        runConfig = RunConfig(args)
        assert(runConfig.run_type == "standard")
        
        pruner = StandardPruner(runConfig, bins_with_reserve_borrowing, self.legend, self.matrix)
        pruner.transform()
        
        after_assigned_bins = pruner.assign_bins()
        
        # Bin 1 should have close to 30 variants (borrowed from extra_rows AND reserve_pool)
        assert(abs(len(after_assigned_bins[1]) - 30) <= 5), f"Bin 1 should have ~30 variants but has {len(after_assigned_bins[1])}"
        
        # Verify that all variants in bin 1 have AC in the range [2, 5]
        # (borrowed variants from reserve_pool were pruned down to fit this bin)
        for row_id in after_assigned_bins[1]:
            row_ac = self.matrix.row_num(row_id)
            assert(2 <= row_ac <= 5), f"Variant {row_id} in bin 1 has AC={row_ac}, expected AC in [2, 5]"
        
        # Reserve pool (bin 3) should have fewer variants than the original 20
        # because some were borrowed for bin 1
        assert(len(after_assigned_bins[3]) < 20), f"Reserve pool should have been reduced from borrowing, but has {len(after_assigned_bins[3])}"
    
    def test_fun_only_mode(self):
        """Test fun_only mode: prune functional variants, keep all synonymous variants"""
        # In fun_only mode, only functional variants (fun='fun') are binned and pruned
        # All synonymous variants (fun='syn') pass through untouched
        # With alternating fun/syn, we have ~50 of each type
        # Bin 0 (AC=1): ~17-18 functional variants (half of 35)
        # Bin 1 (AC=2-5): ~10 functional variants (half of 20)
        # Bin 2 (AC=6-10): ~10 functional variants (half of 20)
        # Bin 3 (AC>10): ~10 functional variants (half of 20)
        
        bins_fun_only = [
            [1, 1, 12],   # Bin 0: AC=1, need 12 functional variants
            [2, 5, 8],    # Bin 1: AC=2-5, need 8 functional variants
            [6, 10, 2]    # Bin 2: AC=6-10, need 2 functional variants
        ]
        
        args = Namespace(z=False, exp_bins=None, fun_bins_only='some_bins', syn_bins_only=None, exp_fun_bins=None, exp_syn_bins=None, small_sample=True, activation_threshold=5, stop_threshold=2, keep_protected=False, output_legend=self.legend_file, input_legend=self.legend_file, verbose=False)
        runConfig = RunConfig(args)
        assert(runConfig.run_type == "fun_only")
        
        pruner = StandardPruner(runConfig, bins_fun_only, self.legend, self.matrix)
        pruner.transform()
        
        after_assigned_bins = pruner.assign_bins()
        
        # Verify bin assignments match targets
        assert(abs(len(after_assigned_bins[0]) - bins_fun_only[0][2]) <= 5)
        assert(abs(len(after_assigned_bins[1]) - bins_fun_only[1][2]) <= 5)
        assert(abs(len(after_assigned_bins[2]) - bins_fun_only[2][2]) <= 5)
        
        # Count functional and synonymous variants in the final output
        rows_to_keep = pruner.get_all_kept_rows(after_assigned_bins)
        
        fun_count = 0
        syn_count = 0
        for row_id in rows_to_keep:
            if self.legend[row_id]['fun'] == 'fun':
                fun_count += 1
            else:
                syn_count += 1
        
        # All synonymous variants with AC > 0 should be kept (exactly 50)
        assert(syn_count == 50), f"Expected exactly 50 synonymous variants to be kept, but got {syn_count}"
    
    def test_syn_only_mode(self):
        """Test syn_only mode: prune synonymous variants, keep all functional variants"""
        # In syn_only mode, only synonymous variants (fun='syn') are binned and pruned
        # All functional variants (fun='fun') pass through untouched
        # With alternating fun/syn, we have ~50 of each type
        
        bins_syn_only = [
            [1, 1, 12],   # Bin 0: AC=1, need 12 synonymous variants
            [2, 5, 8],    # Bin 1: AC=2-5, need 8 synonymous variants
            [6, 10, 2]    # Bin 2: AC=6-10, need 2 synonymous variants
        ]
        
        args = Namespace(z=False, exp_bins=None, fun_bins_only=None, syn_bins_only='some_bins', exp_fun_bins=None, exp_syn_bins=None, small_sample=True, activation_threshold=5, stop_threshold=2, keep_protected=False, output_legend=self.legend_file, input_legend=self.legend_file, verbose=False)
        runConfig = RunConfig(args)
        assert(runConfig.run_type == "syn_only")
        
        pruner = StandardPruner(runConfig, bins_syn_only, self.legend, self.matrix)
        pruner.transform()
        
        after_assigned_bins = pruner.assign_bins()
        
        # Verify bin assignments match targets
        assert(abs(len(after_assigned_bins[0]) - bins_syn_only[0][2]) <= 5)
        assert(abs(len(after_assigned_bins[1]) - bins_syn_only[1][2]) <= 5)
        assert(abs(len(after_assigned_bins[2]) - bins_syn_only[2][2]) <= 5)
        
        # Count functional and synonymous variants in the final output
        rows_to_keep = pruner.get_all_kept_rows(after_assigned_bins)
        
        fun_count = 0
        syn_count = 0
        for row_id in rows_to_keep:
            if self.legend[row_id]['fun'] == 'fun':
                fun_count += 1
            else:
                syn_count += 1
        
        # All functional variants with AC > 0 should be kept (exactly 50)
        assert(fun_count == 50), f"Expected exactly 50 functional variants to be kept, but got {fun_count}"
    
                
class TestFunctionalSplitPruner(unittest.TestCase):
    def setUp(self):
        """Set up test data"""
        random.seed(789)  # Different seed from StandardPruner tests and integration tests
        header = ["id", "position", "a0", "a1", "protected", "fun"]
        legend = Legend(header)
        matrix = SparseMatrix()
        
        # Select protected variants: 5 from each bin for both fun and syn
        self.protected_indices = [
            10, 15, 20, 25, 30,  # 5 from bin 0 (AC=1)
            42, 45, 48, 52, 55,  # 5 from bin 1 (AC=2-5)
            62, 65, 68, 72, 75,  # 5 from bin 2 (AC=6-10)
            82, 85, 88, 92, 95   # 5 from bin 3 (AC>10)
        ]
        
        for i in range(100):
            protected = "1" if i in self.protected_indices else "0"
            fun = "fun" if i % 2 == 0 else "syn"  # Alternate between fun and syn
            legend.add_row([f"allele{i}", f"{i}", "C", "G", protected, fun])
            if i < 5:
                matrix.add_row([])
            if 5 <= i < 40:
                columns = random.sample(range(30), 1)
                matrix.add_row(sorted(columns))
            if 40 <= i < 60:
                num_columns = random.randrange(2, 5)
                columns = random.sample(range(30), num_columns)
                matrix.add_row((sorted(columns)))
            if 60 <= i < 80:
                num_columns = random.randrange(6, 10)
                columns = random.sample(range(30), num_columns)
                matrix.add_row((sorted(columns)))
            if 80 <= i:
                num_columns = random.randrange(11, 30)
                columns = random.sample(range(30), num_columns)
                matrix.add_row((sorted(columns)))

        self.legend = legend
        self.matrix = matrix
        # Separate bins for functional and synonymous variants
        self.bins = {
            'fun': [
                [1, 1, 12],   # Bin 0: AC=1, need 12 functional variants
                [2, 5, 8],    # Bin 1: AC=2-5, need 8 functional variants
                [6, 10, 2]    # Bin 2: AC=6-10, need 2 functional variants
            ],
            'syn': [
                [1, 1, 15],   # Bin 0: AC=1, need 15 synonymous variants
                [2, 5, 6],    # Bin 1: AC=2-5, need 6 synonymous variants
                [6, 10, 3]    # Bin 2: AC=6-10, need 3 synonymous variants
            ]
        }
        self.temp_dir = tempfile.mkdtemp()
        self.legend_file = os.path.join(self.temp_dir, "legend")
    
    def tearDown(self):
        """Clean up test files"""
        # Remove all files in temp directory
        for file in os.listdir(self.temp_dir):
            os.remove(os.path.join(self.temp_dir, file))
        os.rmdir(self.temp_dir)
    
    def test_assign_bins(self):
        """Test that variants are correctly assigned to functional and synonymous bins"""
        args = Namespace(z=False, exp_bins=None, exp_fun_bins='some_bins', exp_syn_bins='some_bins', small_sample=True, activation_threshold=5, stop_threshold=2, keep_protected=False, output_legend=self.legend_file, input_legend=self.legend_file)
        run_config = RunConfig(args)
        assert(run_config.run_type == "func_split")
        
        pruner = FunctionalSplitPruner(run_config, self.bins, self.legend, self.matrix)
        assigned_bins = pruner.assign_bins()
        
        # Verify structure
        assert('fun' in assigned_bins)
        assert('syn' in assigned_bins)
        
        # Count functional and synonymous variants
        fun_total = sum([len(assigned_bins['fun'][i]) for i in assigned_bins['fun']])
        syn_total = sum([len(assigned_bins['syn'][i]) for i in assigned_bins['syn']])
        
        # Should have ~50 of each (excluding AC=0 variants)
        assert(fun_total >= 45 and fun_total <= 50), f"Expected ~50 functional variants, got {fun_total}"
        assert(syn_total >= 45 and syn_total <= 50), f"Expected ~50 synonymous variants, got {syn_total}"
    
    def test_full_transform_no_z(self):
        """Test full transform without z-flag (keep zeroed rows)"""
        args = Namespace(z=False, exp_bins=None, exp_fun_bins='some_bins', exp_syn_bins='some_bins', small_sample=True, activation_threshold=5, stop_threshold=2, keep_protected=False, output_legend=self.legend_file, input_legend=self.legend_file)
        runConfig = RunConfig(args)
        assert(runConfig.run_type == "func_split")
        
        pruner = FunctionalSplitPruner(runConfig, self.bins, self.legend, self.matrix)
        pruner.transform()
        
        after_assigned_bins = pruner.assign_bins()
        
        # Verify functional bins match targets
        assert(abs(len(after_assigned_bins['fun'][0]) - self.bins['fun'][0][2]) <= 5)
        assert(abs(len(after_assigned_bins['fun'][1]) - self.bins['fun'][1][2]) <= 5)
        assert(abs(len(after_assigned_bins['fun'][2]) - self.bins['fun'][2][2]) <= 5)
        
        # Verify synonymous bins match targets
        assert(abs(len(after_assigned_bins['syn'][0]) - self.bins['syn'][0][2]) <= 5)
        assert(abs(len(after_assigned_bins['syn'][1]) - self.bins['syn'][1][2]) <= 5)
        assert(abs(len(after_assigned_bins['syn'][2]) - self.bins['syn'][2][2]) <= 5)
        
        # z flag is false. Assert that the number of rows in the matrix and legend is 100 still.
        assert(self.legend.row_count() == 100)
        assert(self.matrix.num_rows() == 100)
    
    def test_full_transform_with_z(self):
        """Test full transform with z-flag (remove zeroed rows)"""
        args = Namespace(z=True, exp_bins=None, exp_fun_bins='some_bins', exp_syn_bins='some_bins', small_sample=True, activation_threshold=5, stop_threshold=2, keep_protected=False, output_legend=self.legend_file, input_legend=self.legend_file)
        runConfig = RunConfig(args)
        assert(runConfig.run_type == "func_split")
        
        pruner = FunctionalSplitPruner(runConfig, self.bins, self.legend, self.matrix)
        pruner.transform()
        
        after_assigned_bins = pruner.assign_bins()
        
        # Verify functional bins match targets
        assert(abs(len(after_assigned_bins['fun'][0]) - self.bins['fun'][0][2]) <= 5)
        assert(abs(len(after_assigned_bins['fun'][1]) - self.bins['fun'][1][2]) <= 5)
        assert(abs(len(after_assigned_bins['fun'][2]) - self.bins['fun'][2][2]) <= 5)
        
        # Verify synonymous bins match targets
        assert(abs(len(after_assigned_bins['syn'][0]) - self.bins['syn'][0][2]) <= 5)
        assert(abs(len(after_assigned_bins['syn'][1]) - self.bins['syn'][1][2]) <= 5)
        assert(abs(len(after_assigned_bins['syn'][2]) - self.bins['syn'][2][2]) <= 5)
        
        # z flag is true. Assert that the number of rows has been reduced.
        total_kept_variants = sum([len(after_assigned_bins['fun'][i]) for i in after_assigned_bins['fun']])
        total_kept_variants += sum([len(after_assigned_bins['syn'][i]) for i in after_assigned_bins['syn']])
        
        assert(self.legend.row_count() == total_kept_variants)
        assert(self.matrix.num_rows() == total_kept_variants)
    
    def test_full_transform_with_keep_protected(self):
        """Test that protected variants are not pruned in func_split mode"""
        args = Namespace(z=False, exp_bins=None, exp_fun_bins='some_bins', exp_syn_bins='some_bins', small_sample=True, activation_threshold=5, stop_threshold=2, keep_protected=True, output_legend=self.legend_file, input_legend=self.legend_file, verbose=True)
        runConfig = RunConfig(args)
        assert(runConfig.run_type == "func_split")
        
        pruner = FunctionalSplitPruner(runConfig, self.bins, self.legend, self.matrix)
        pruner.transform()
        
        after_assigned_bins = pruner.assign_bins()
        
        # Verify all protected variants are still in the bins (not pruned)
        all_kept_rows = []
        for bin_id in after_assigned_bins['fun']:
            all_kept_rows.extend(after_assigned_bins['fun'][bin_id])
        for bin_id in after_assigned_bins['syn']:
            all_kept_rows.extend(after_assigned_bins['syn'][bin_id])
        
        for protected_idx in self.protected_indices:
            if self.matrix.row_num(protected_idx) > 0:  # Only check if variant has AC > 0
                assert(protected_idx in all_kept_rows), f"Protected variant {protected_idx} was pruned!"
        
        # Verify pruned variants file doesn't contain protected variants
        pruned_file = f'{self.legend_file}-pruned-variants'
        assert(os.path.exists(pruned_file))
        
        with open(pruned_file, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:  # Skip header
                variant_id = line.split('\t')[0]
                idx = int(variant_id.replace('allele', ''))
                assert(idx not in self.protected_indices), f"Protected variant {idx} was written to pruned file!"
    
    def test_borrowing_no_crossover(self):
        """Test that functional and synonymous variants cannot borrow from each other's extra_rows"""
        # Create scenario where both fun and syn need to borrow
        # But we'll verify they only borrow from their own extra_rows pools
        bins_with_borrowing = {
            'fun': [
                [1, 1, 12],   # Bin 0: AC=1, need 12 (have ~17)
                [2, 5, 15],   # Bin 1: AC=2-5, need 15 (have ~10, need to borrow ~5)
                [6, 10, 2]    # Bin 2: AC=6-10, need 2 (have ~10, will prune ~8 into extra_rows)
            ],
            'syn': [
                [1, 1, 12],   # Bin 0: AC=1, need 12 (have ~17)
                [2, 5, 15],   # Bin 1: AC=2-5, need 15 (have ~10, need to borrow ~5)
                [6, 10, 2]    # Bin 2: AC=6-10, need 2 (have ~10, will prune ~8 into extra_rows)
            ]
        }
        
        args = Namespace(z=False, exp_bins=None, exp_fun_bins='some_bins', exp_syn_bins='some_bins', small_sample=True, activation_threshold=5, stop_threshold=2, keep_protected=False, output_legend=self.legend_file, input_legend=self.legend_file, verbose=False)
        runConfig = RunConfig(args)
        assert(runConfig.run_type == "func_split")
        
        pruner = FunctionalSplitPruner(runConfig, bins_with_borrowing, self.legend, self.matrix)
        pruner.transform()
        
        after_assigned_bins = pruner.assign_bins()
        
        # Verify both fun and syn bin 1 have close to 15 variants (borrowed from their own extra_rows)
        assert(abs(len(after_assigned_bins['fun'][1]) - 15) <= 5), f"Fun bin 1 should have ~15 variants but has {len(after_assigned_bins['fun'][1])}"
        assert(abs(len(after_assigned_bins['syn'][1]) - 15) <= 5), f"Syn bin 1 should have ~15 variants but has {len(after_assigned_bins['syn'][1])}"
        
        # Verify all variants in fun bins are functional
        for bin_id in after_assigned_bins['fun']:
            for row_id in after_assigned_bins['fun'][bin_id]:
                assert(self.legend[row_id]['fun'] == 'fun'), f"Found non-functional variant {row_id} in functional bin {bin_id}"
        
        # Verify all variants in syn bins are synonymous
        for bin_id in after_assigned_bins['syn']:
            for row_id in after_assigned_bins['syn'][bin_id]:
                assert(self.legend[row_id]['fun'] == 'syn'), f"Found non-synonymous variant {row_id} in synonymous bin {bin_id}"


class TestProbabilisticPruner(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.legend_file = os.path.join(self.temp_dir, "legend")
        self.legend = Legend(["id", "position", "a0", "a1", "prob"])
        self.legend.add_row(["var0", "1", "A", "G", "0.0"])
        self.legend.add_row(["var1", "2", "A", "G", "1.0"])
        self.legend.add_row(["var2", "3", "A", "G", "."])

        self.matrix = SparseMatrix(cols=4)
        self.matrix.add_row([0, 1, 2])
        self.matrix.add_row([1, 3])
        self.matrix.add_row([0, 2])

    def tearDown(self):
        pruned_file = f'{self.legend_file}-pruned-variants'
        if os.path.exists(pruned_file):
            os.remove(pruned_file)
        os.rmdir(self.temp_dir)

    def test_probabilistic_pruning_removes_whole_rows(self):
        args = Namespace(
            z=False,
            exp_bins=None,
            exp_fun_bins=None,
            exp_syn_bins=None,
            fun_bins_only=None,
            syn_bins_only=None,
            prob=True,
            small_sample=True,
            activation_threshold=10,
            stop_threshold=20,
            keep_protected=False,
            output_legend=self.legend_file,
            input_legend=self.legend_file,
            verbose=False
        )
        run_config = RunConfig(args)
        pruner = ProbabilisticPruner(run_config, self.legend, self.matrix)

        pruner.transform()

        self.assertEqual(self.matrix.row_num(0), 0, "prob=0.0 should always prune the whole variant")
        self.assertEqual(self.matrix.get_row(0), [0, 0, 0, 0])
        self.assertEqual(self.matrix.get_row(1), [0, 1, 0, 1], "prob=1.0 should always keep the variant")
        self.assertEqual(self.matrix.get_row(2), [1, 0, 1, 0], "'.' should leave the variant unchanged")

        pruned_file = f'{self.legend_file}-pruned-variants'
        self.assertTrue(os.path.exists(pruned_file))
        with open(pruned_file, 'r') as f:
            lines = f.readlines()
        self.assertEqual(len(lines), 2)
        self.assertIn("var0", lines[1])

    def test_probabilistic_pruning_prints_allele_frequency_distributions(self):
        args = Namespace(
            z=False,
            exp_bins=None,
            exp_fun_bins=None,
            exp_syn_bins=None,
            fun_bins_only=None,
            syn_bins_only=None,
            prob=True,
            small_sample=True,
            activation_threshold=10,
            stop_threshold=20,
            keep_protected=False,
            output_legend=self.legend_file,
            input_legend=self.legend_file,
            verbose=False
        )
        run_config = RunConfig(args)
        pruner = ProbabilisticPruner(run_config, self.legend, self.matrix)

        stdout = io.StringIO()
        with redirect_stdout(stdout):
            pruner.transform()

        output = stdout.getvalue()
        self.assertIn("Input allele frequency distribution:", output)
        self.assertIn("New allele frequency distribution:", output)
        self.assertIn("Bin", output)
        self.assertIn("Expected", output)
        self.assertIn("Actual", output)
