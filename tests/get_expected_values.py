"""
Helper script to determine exact expected values for deterministic tests.
Run this to get the expected outputs with specific random seeds.
"""
import random
import tempfile
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from raresim.common.sparse import SparseMatrixReader
from raresim.common.legend import LegendReaderWriter
from raresim.engine.config import RunConfig
from raresim.engine.runner import DefaultRunner
from argparse import Namespace

data_dir = os.path.join(os.path.dirname(__file__), 'data')
temp_dir = tempfile.mkdtemp()

reader = SparseMatrixReader()

print("=" * 60)
print("STANDARD PRUNING (seed=123)")
print("=" * 60)
random.seed(123)
haps_file = os.path.join(data_dir, 'test_standard.haps')
legend_file = os.path.join(temp_dir, 'temp.legend')
output_hap = os.path.join(temp_dir, 'standard_output.haps.gz')
output_legend = os.path.join(temp_dir, 'standard_output.legend')

original_matrix = reader.loadSparseMatrix(haps_file)
with open(legend_file, 'w') as f:
    f.write("id\tposition\ta0\ta1\n")
    for i in range(original_matrix.num_rows()):
        f.write(f"variant_{i}\t{i*1000}\tA\tG\n")

args = Namespace(
    sparse_matrix=haps_file,
    input_legend=legend_file,
    output_legend=output_legend,
    output_hap=output_hap,
    exp_bins=os.path.join(data_dir, 'test_standard_bins.txt'),
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
runner = DefaultRunner(config)
runner.run()

output_matrix = reader.loadSparseMatrix(output_hap)
output_legend_obj = LegendReaderWriter.load_legend(output_legend)
print(f"Output rows: {output_matrix.num_rows()}")
print(f"Output legend rows: {output_legend_obj.row_count()}")

print("\n" + "=" * 60)
print("FUNCTIONAL SPLIT (seed=456)")
print("=" * 60)
random.seed(456)
haps_file = os.path.join(data_dir, 'test_stratified.haps')
legend_file = os.path.join(data_dir, 'test_stratified.legend')
output_hap = os.path.join(temp_dir, 'split_output.haps.gz')
output_legend = os.path.join(temp_dir, 'split_output.legend')

original_matrix = reader.loadSparseMatrix(haps_file)
original_legend = LegendReaderWriter.load_legend(legend_file)

args = Namespace(
    sparse_matrix=haps_file,
    input_legend=legend_file,
    output_legend=output_legend,
    output_hap=output_hap,
    exp_bins=None,
    exp_fun_bins=os.path.join(data_dir, 'test_fun_bins.txt'),
    exp_syn_bins=os.path.join(data_dir, 'test_syn_bins.txt'),
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
runner = DefaultRunner(config)
runner.run()

output_matrix = reader.loadSparseMatrix(output_hap)
output_legend_obj = LegendReaderWriter.load_legend(output_legend)
output_fun = sum(1 for i in range(output_legend_obj.row_count()) 
                if output_legend_obj[i]['fun'] == 'fun')
output_syn = sum(1 for i in range(output_legend_obj.row_count()) 
                if output_legend_obj[i]['fun'] == 'syn')

print(f"Output rows: {output_matrix.num_rows()}")
print(f"Output functional: {output_fun}")
print(f"Output synonymous: {output_syn}")

# Cleanup
import shutil
shutil.rmtree(temp_dir)
