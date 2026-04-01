from raresim.common.sparse import SparseMatrixReader, SparseMatrixWriter
from raresim.engine.runner import DefaultRunner
from raresim.engine.config import RunConfig
from raresim.calculate.expected_vars import calc
import argparse
import random
import os
import gzip


def parseCommand():
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(dest='command', required=True, help='Available commands: sim, convert, extract')

    sim_parser = subparsers.add_parser('sim')
    convert_parser = subparsers.add_parser('convert')
    extract_parser = subparsers.add_parser('extract')
    calc_parser = subparsers.add_parser('calc')

    calc_parser.add_argument('--mac',
                        dest='mac',
                        required=True,
                        help='MAC bin bounds (lower and upper allele counts) for the simulated sample size')
    calc_parser.add_argument('-o',
                        dest='output',
                        required=True,
                        help='Output file name')
    calc_parser.add_argument('-N',
                        dest='n',
                        required=True,
                        help='Number of individuals in the simulated sample')
    calc_parser.add_argument('--pop',
                        dest='pop',
                        help='Population (AFR, EAS, NFE, or SAS) to use default values for if not providing alpha, beta, omega, phi, and b values or target data')
    calc_parser.add_argument('--alpha',
                        dest='alpha',
                        help='Shape parameter to estimate the expected AFS distribution (must be > 0)')
    calc_parser.add_argument('--beta',
                        dest='beta',
                        help='Shape parameter to estimate the expected AFS distribution')
    calc_parser.add_argument('--omega',
                        dest='omega',
                        help='Scaling parameter to estimate the expected number of variants per (Kb) for sample size N (range of 0-1)')
    calc_parser.add_argument('--phi',
                        dest='phi',
                        help='Shape parameter to estimate the expected number of variants per (Kb) for sample size N (must be > 0)')
    calc_parser.add_argument('-b',
                        dest='b',
                        help='Scale parameter to estimate the expected AFS distribution')
    calc_parser.add_argument('--nvar_target_data',
                        dest='nvar_target_data',
                        help='Target downsampling data with the number of variants per Kb to estimate the expected number of variants per Kb for sample size N')
    calc_parser.add_argument('--afs_target_data',
                        dest='afs_target_data',
                        help='Target AFS data with the proportion of variants per MAC bin to estimate the expected AFS distribution')
    calc_parser.add_argument('--reg_size',
                        dest='reg_size',
                        help='Size of simulated genetic region in kilobases (Kb)')
    calc_parser.add_argument('-w',
                        dest='w',
                        default='1.0',
                        help='Weight to multiply the expected number of variants by in non-stratified simulations (default value of 1)')
    calc_parser.add_argument('--w_fun',
                        dest='w_fun',
                        default='1.0',
                        help='Weight to multiply the expected number of functional variants by in stratified simulations (default value of 1)')
    calc_parser.add_argument('--w_syn',
                        dest='w_syn',
                        default='1.0',
                        help='Weight to multiply the expected number of synonymous variants by in stratified simulations (default value of 1)')

    extract_parser.add_argument('-i',
                                dest='input_file',
                                required=True,
                                help='Input haplotype file (gzipped)')
    extract_parser.add_argument('-o',
                                dest='output_file',
                                required=True,
                                help='Output haplotype file name')
    extract_parser.add_argument('-s', '--seed',
                                dest='seed',
                                type=int,
                                help='Optional seed for reproducibility')
    extract_parser.add_argument('-n',
                                dest='num',
                                type=int,
                                required=True,
                                help='Number of haplotypes to extract')

    sim_parser.add_argument('-m',
                        dest='sparse_matrix',
                        required=True,
                        help='Input haplotype file (can be a .haps, .sm, or .gz file)')

    sim_parser.add_argument('-b',
                        dest='exp_bins',
                        help='Expected number of functional and synonymous variants per MAC bin')

    sim_parser.add_argument('--functional_bins',
                        dest='exp_fun_bins',
                        help='Expected number of variants per MAC bin for functional variants (must be used with --synonymous_bins) ')

    sim_parser.add_argument('--synonymous_bins',
                        dest='exp_syn_bins',
                        help='Expected number of variants per MAC bin for synonymous variants (must be used with --functional_bins) ')

    sim_parser.add_argument('-l',
                        dest='input_legend',
                        required=True,
                        help='Input legend file')

    sim_parser.add_argument('-L',
                        dest='output_legend',
                        help='Output legend file (only required when using -z)')

    sim_parser.add_argument('-H',
                        dest='output_hap',
                        required=True,
                        help='Output compressed haplotype file')

    sim_parser.add_argument('--f_only',
                        dest='fun_bins_only',
                        help='Expected number of variants per MAC bin for only functional variants')

    sim_parser.add_argument('--s_only',
                        dest='syn_bins_only',
                        help='Expected number of variants per MAC bin for only synonymous variants')

    sim_parser.add_argument('-z',
                        action='store_true',
                        help='Monomorphic and pruned variants (rows of zeros) are removed from the output haplotype file')

    sim_parser.add_argument('-prob',
                        action='store_true',
                        help='Variants are pruned row by row using the keep probability in the legend file')

    sim_parser.add_argument('--small_sample',
                        action='store_true',
                        help='Overrides error to allow for simulation of small sample sizes (<10,000 haplotypes)')

    sim_parser.add_argument('--keep_protected',
                        action='store_true',
                        help='Variants designated with a 1 in the protected column of the legend file will not be pruned')

    sim_parser.add_argument('--stop_threshold',
                        dest='stop_threshold',
                        default='20',
                        help='Percentage threshold for stopping the pruning process (0-100). Prevents the number of variants from falling below the specified percentage of the expected count for any given MAC bin during pruning (default value of 20)')

    sim_parser.add_argument('--activation_threshold',
                        dest='activation_threshold',
                        default='10',
                        help='Percentage threshold for activating the pruning process (0-100). Requires that the actual number of variants for a MAC bin must be more than the given percentage different from the expected number to activate pruning on the bin')

    sim_parser.add_argument('--verbose',
                        action='store_true',
                        help='when using --keep_protected and this flag, the program will additionally print the before and after Allele Frequency Distributions with the protected variants pulled out')
    convert_parser.add_argument('-i',
                                dest='input_file',
                                required=True,
                                help='Input haplotype file')

    convert_parser.add_argument('-o',
                                dest='output_file',
                                required=True,
                                help='Output haplotype file')

    args = parser.parse_args()

    return args

def extract(args):
    random.seed(args.seed)
    with gzip.open(args.input_file, 'rt') as f:
        line = f.readline()
        columns = line.split()
    size = len(columns)
    columnsToExtract = random.sample(range(0, size), args.num)
    otherColumns = [i for i in range(size) if i not in columnsToExtract]
    columnsToExtract.sort()
    base, ext = os.path.splitext(args.output_file)
    output_file_name = base
    with gzip.open(f'{output_file_name}-sample.gz', 'wb') as s:
        with gzip.open(f'{output_file_name}-remainder.gz', 'wb') as r:
            with gzip.open(args.input_file, 'rt') as input_haps:
                for l in input_haps.readlines():
                    cols = l.split()
                    sampleLine = [cols[i] for i in columnsToExtract]
                    remainderLine = [cols[i] for i in otherColumns]
                    s.write((" ".join(sampleLine) + "\n").encode())
                    r.write((" ".join(remainderLine) + "\n").encode())


def main():
    command = parseCommand()
    if command.command == 'sim':
        runConfig = RunConfig(command)
        runner: DefaultRunner = DefaultRunner(runConfig)
        runner.run()

    elif command.command == 'convert':
        args = command
        reader = SparseMatrixReader()
        writer = SparseMatrixWriter()
        matrix = reader.loadSparseMatrix(args.input_file)
        writer.writeToHapsFile(matrix, args.output_file, "sm")

    elif command.command == 'extract':
        args = command
        extract(args)

    elif command.command == 'calc':
        args = command
        calc(args)


if __name__ == '__main__':
    main()
