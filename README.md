[![PyPI version](https://badge.fury.io/py/raresim.svg)](https://badge.fury.io/py/raresim)
[![Python Version](https://img.shields.io/pypi/pyversions/raresim.svg)](https://pypi.org/project/raresim/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# RAREsim2
Python interface for flexible simulation of rare-variant genetic data using real haplotypes


## Installation

### From PyPI
```bash
pip install raresim
```

### From TestPyPI (for testing pre-releases)
```bash
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ raresim
```

### From Source
```bash
git clone https://github.com/RMBarnard/raresim.git
cd raresim
pip install -e .  # Install in development mode
```

## Main Functions

### CALC
Calculate the expected number of variants per MAC bin using default population parameters, user-provided parameters, or target data.
```
usage: __main__.py calc [-h] --mac MAC -o OUTPUT -N N [--pop POP]
                        [--alpha ALPHA] [--beta BETA] [--omega OMEGA]
                        [--phi PHI] [-b B]
                        [--nvar_target_data NVAR_TARGET_DATA]
                        [--afs_target_data AFS_TARGET_DATA]
                        [--reg_size REG_SIZE] [-w W] [--w_fun W_FUN]
                        [--w_syn W_SYN]

options:
  -h, --help            show this help message and exit
  --mac MAC             MAC bin bounds (lower and upper allele counts) for the simulated sample size
  -o OUTPUT             Output file name
  -N N                  Number of individuals in the simulated sample
  --pop POP             Population (AFR, EAS, NFE, or SAS) to use default values for if not providing
                        alpha, beta, omega, phi, and b values or target data
  --alpha ALPHA         Shape parameter to estimate the expected AFS distribution (must be > 0)
  --beta BETA           Shape parameter to estimate the expected AFS distribution
  --omega OMEGA         Scaling parameter to estimate the expected number of variants per (Kb) for
                        sample size N (range of 0-1)
  --phi PHI             Shape parameter to estimate the expected number of variants per (Kb) for
                        sample size N (must be > 0)
  -b B                  Scale parameter to estimate the expected AFS distribution
  --nvar_target_data NVAR_TARGET_DATA
                        Target downsampling data with the number of variants per Kb to estimate the
                        expected number of variants per Kb for sample size N
  --afs_target_data AFS_TARGET_DATA
                        Target AFS data with the proportion of variants per MAC bin to estimate the
                        expected AFS distribution
  --reg_size REG_SIZE   Size of simulated genetic region in kilobases (Kb)
  -w W                  Weight to multiply the expected number of variants by in non-stratified
                        simulations (default value of 1)
  --w_fun W_FUN         Weight to multiply the expected number of functional variants by in
                        stratified simulations (default value of 1)
  --w_syn W_SYN         Weight to multiply the expected number of synonymous variants by in
                        stratified simulations (default value of 1)
```

#### Default Population Parameters
The expected number of functional and synonymous variants can be estimated using default parameters for the following populations: African (AFR), East Asian (EAS), Non-Finnish European (NFE), and South Asian (SAS).

```text
$ python3 -m raresim calc \
    --mac example/mac_bins.txt \
    -o example/mac_bin_estimates_default.txt \
    -N 10000 \
    --pop NFE \
    --reg_size 19.029

Calculated 842.5888117489534 total variants (accounting for region size)
```

#### Target Data
The user can also use their own target data - this is necessary to calculate the expected number of functional and/or synonymous variants for stratified simulations. Note, the simulation parameters are output if the user wants to use them instead of target data for future simulations.

```text
$ python3 -m raresim calc \
    --mac example/mac_bins.txt \
    -o example/mac_bin_estimates_target.txt \
    -N 10000 \
    --nvar_target_data example/nvar_target.txt \
    --afs_target_data example/afs_target.txt \
    --reg_size 19.029

Calculating synonymous values
Calculated the following params from AFS target data. alpha: 1.9397807693228122, beta: 0.34101610369526514, b: 0.8464846288340953
Calculated the following params from nvar target data. omega: 0.6295595643083463, phi: 0.04392478579419536
Calculated 275.6537313477067 total variants (accounting for region size)

Calculating functional values
Calculated the following params from AFS target data. alpha: 2.1388159441481442, beta: 0.4285647164342115, b: 1.134635990601139
Calculated the following params from nvar target data. omega: 0.6413547202832528, phi: 0.08338724275310817
Calculated 583.3570639000195 total variants (accounting for region size)
```

Note: Two MAC bin estimate files will be output (one for functional variants and another for synonymous variants) if the 
input AFS file is stratified by functional status. If it's not stratified, then just one file will be output.

#### User-Provided Parameters
If parameters are known from previous simulations, the user can provide those instead of having to provide and fit target data.

```text
$ python3 -m raresim calc \
    --mac example/mac_bins.txt \
    -o example/mac_bin_estimates_params.txt \
    -N 10000 \
    --alpha 1.947 \
    --beta 0.118 \
    -b 0.6676 \
    --omega 0.6539 \
    --phi 0.1073 \
    --reg_size 19.029

Calculated 842.5888117489534 total variants (accounting for region size)
```

### SIM
Simulate new allele frequencies by pruning (i.e., removing) certain variants from an input haplotype file given the expected number of variants for the simulated sample size. A list of pruned variants (.legend-pruned-variants) is also output along with the new haplotype file.
```
usage: __main__.py sim [-h] -m SPARSE_MATRIX [-b EXP_BINS]
                       [--functional_bins EXP_FUN_BINS]
                       [--synonymous_bins EXP_SYN_BINS] -l INPUT_LEGEND
                       [-L OUTPUT_LEGEND] -H OUTPUT_HAP
                       [--f_only FUN_BINS_ONLY] [--s_only SYN_BINS_ONLY] [-z]
                       [-prob] [--small_sample] [--keep_protected]
                       [--stop_threshold STOP_THRESHOLD]
                       [--activation_threshold ACTIVATION_THRESHOLD]
                       [--verbose]

options:
  -h, --help            show this help message and exit
  -m SPARSE_MATRIX      Input haplotype file (can be a .haps, .sm, or .gz file)
  -b EXP_BINS           Expected number of functional and synonymous variants per MAC bin
  --functional_bins EXP_FUN_BINS
                        Expected number of variants per MAC bin for functional variants (must be used
                        with --synonymous_bins) 
  --synonymous_bins EXP_SYN_BINS
                        Expected number of variants per MAC bin for synonymous variants (must be used
                        with --functional_bins) 
  -l INPUT_LEGEND       Input legend file
  -L OUTPUT_LEGEND      Output legend file (only required when using -z)
  -H OUTPUT_HAP         Output compressed haplotype file
  --f_only FUN_BINS_ONLY
                        Expected number of variants per MAC bin for only functional variants
  --s_only SYN_BINS_ONLY
                        Expected number of variants per MAC bin for only synonymous variants
  -z                    Monomorphic and pruned variants (rows of zeros) are removed from the output
                        haplotype file
  -prob                 Variants are pruned row by row using the keep probability in the
                        legend file
  --small_sample        Overrides error to allow for simulation of small sample sizes (<10,000
                        haplotypes)
  --keep_protected      Variants designated with a 1 in the protected column of the legend file will
                        not be pruned
  --stop_threshold STOP_THRESHOLD
                        Percentage threshold for stopping the pruning process (0-100). Prevents the
                        number of variants from falling below the specified percentage of the expected
                        count for any given MAC bin during pruning (default value of 20)
  --activation_threshold ACTIVATION_THRESHOLD
                        Percentage threshold for activating the pruning process (0-100). Requires that
                        the actual number of variants for a MAC bin must be more than the given
                        percentage different from the expected number to activate pruning on the bin
                        (default value of 10)
  --verbose             when using --keep_protected and this flag, the program will additionally print
                        the before and after Allele Frequency Distributions with the protected variants
                        pulled out
```

```text
$ python3 -m raresim sim \
    -m example/example.haps.gz \
    -b example/mac_bin_estimates_default.txt \
    -l example/example.legend \
    -L example/output.legend \
    -H example/output.haps.gz \
    -z

Running with run mode: standard
Input allele frequency distribution:
Bin         Expected          Actual
[1,1]       452.7055560068    1002
[2,2]       130.4830742030    484
[3,5]       120.6258509819    768
[6,10]      52.2181585555     663
[11,20]     29.5461366439     681
[21,100]    26.2774091990     856
[101,200]   3.6164427260      79
[201,∞]     N/A               65

New allele frequency distribution:
Bin         Expected          Actual
[1,1]       452.7055560068    472
[2,2]       130.4830742030    119
[3,5]       120.6258509819    110
[6,10]      52.2181585555     48
[11,20]     29.5461366439     28
[21,100]    26.2774091990     47
[101,200]   3.6164427260      3
[201,∞]     N/A               65

Writing new variant legend

Writing new haplotype file
[====================] 100%
```
Note: An updated legend file is only output when using the z flag (when pruned variants are removed from the haplotype file). If not using the z flag,
then the order and amount of rows (i.e., variants) in the haplotype file will remain unchanged and match the input legend file. Also, if the input
haplotype file contains monomorphic variants (i.e., rows of zeros) when using the z flag, then the .legend-pruned-variants file will contain both
monomorphic and actual pruned variants.

#### Stratified (Functional/Synonymous) Pruning
To perform stratified simulations where functional and synonymous variants are pruned separately:
1. add a column to the legend file (`-l`) named "fun", where functional variants have the value "fun" and synonymous variants have the value "syn"
2. provide separate MAC bin files with the expected number of variants per bin for functional (`--functional_bins`) and synonymous (`--synonymous_bins`) variants
   
```text
$ python3 -m raresim sim \
    -m example/example.haps.gz \
    --functional_bins example/mac_bin_estimates_target_fun.txt \
    --synonymous_bins example/mac_bin_estimates_target_syn.txt \
    -l example/example.legend \
    -L example/output_stratified.legend \
    -H example/output_stratified.haps.gz \
    -z

Running with run mode: func_split
Input allele frequency distribution:
Functional
Bin         Expected          Actual
[1,1]       308.6658613719    706
[2,2]       99.2199432898     332
[3,5]       92.6656147375     541
[6,10]      38.2293812491     463
[11,20]     19.9237792915     489
[21,100]    15.1688219483     607
[101,200]   1.6493333218      52
[201,∞]     N/A               46

Synonymous
Bin         Expected          Actual
[1,1]       132.0653670095    296
[2,2]       44.8145869897     152
[3,5]       45.0536145138     227
[6,10]      20.7498071235     200
[11,20]     12.1186468959     192
[21,100]    11.0509676181     249
[101,200]   1.5493808935      27
[201,∞]     N/A               19

New allele frequency distribution:
Functional
Bin         Expected          Actual
[1,1]       308.6658613719    290
[2,2]       99.2199432898     99
[3,5]       92.6656147375     88
[6,10]      38.2293812491     47
[11,20]     19.9237792915     18
[21,100]    15.1688219483     22
[101,200]   1.6493333218      1
[201,∞]     N/A               46

Synonymous
Bin         Expected          Actual
[1,1]       132.0653670095    134
[2,2]       44.8145869897     42
[3,5]       45.0536145138     51
[6,10]      20.7498071235     22
[11,20]     12.1186468959     11
[21,100]    11.0509676181     11
[101,200]   1.5493808935      2
[201,∞]     N/A               19

Writing new variant legend

Writing new haplotype file
[====================] 100%
```

#### Only Functional/Synonymous Variants
To prune only functional or only synonymous variants:
1. add a column to the legend file (`-l`) named "fun", where functional variants have the value "fun" and synonymous variants have the value "syn"
2. provide a MAC bin file with the expected number of variants per bin for only functional (`--f_only`) or only synonymous (`--s_only`) variants

```text
$ python3 -m raresim sim \
    -m example/example.haps.gz \
    --f_only example/mac_bin_estimates_target_fun.txt \
    -l example/example.legend \
    -L example/output_fun_only.legend \
    -H example/output_fun_only.haps.gz \
    -z

Running with run mode: fun_only
Input allele frequency distribution:
Bin         Expected          Actual
[1,1]       308.6658613719    706
[2,2]       99.2199432898     332
[3,5]       92.6656147375     541
[6,10]      38.2293812491     463
[11,20]     19.9237792915     489
[21,100]    15.1688219483     607
[101,200]   1.6493333218      52
[201,∞]     N/A               46

New allele frequency distribution:
Bin         Expected          Actual
[1,1]       308.6658613719    312
[2,2]       99.2199432898     92
[3,5]       92.6656147375     102
[6,10]      38.2293812491     38
[11,20]     19.9237792915     17
[21,100]    15.1688219483     15
[101,200]   1.6493333218      2
[201,∞]     N/A               46

Writing new variant legend

Writing new haplotype file
[====================] 100%
```

#### Given Probabilities
To prune variants using known or given probabilities, add a column to the legend file (`-l`) named `prob`. A single random draw is generated for each variant row, and the row is kept with the probability given in the legend. When using the `-z` flag, fully pruned and monomorphic variants are removed from the output haplotype file, and a pruned-variants file is created.

```text
$ python3 -m raresim sim \
    -m example/example.haps.gz \
    -l example/example.legend \
    -L example/output_probs.legend \
    -H example/output_probs.haps.gz \
    -prob \
    -z

Running with run mode: probabilistic

Writing new variant legend

Writing new haplotype file
[====================] 100%
```


#### Protected Status
To exclude protected variants from the pruning process, add a column to the legend file (`-l`) named "protected". Any row with a 0 in this column will be eligible for pruning while any row with a 1 will still be counted but will not be eligible for pruning.

```text
$ python3 -m raresim sim \
    -m example/example.haps.gz \
    -b example/mac_bin_estimates_default.txt \
    -l example/example.protected.legend \
    -L example/output_protected.legend \
    -H example/output_protected.haps.gz \
    --keep_protected \
    -z

Running with run mode: standard
Input allele frequency distribution:
Bin         Expected          Actual
[1,1]       452.7055560068    1002
[2,2]       130.4830742030    484
[3,5]       120.6258509819    768
[6,10]      52.2181585555     663
[11,20]     29.5461366439     681
[21,100]    26.2774091990     856
[101,200]   3.6164427260      79
[201,∞]     N/A               65

New allele frequency distribution:
Bin         Expected          Actual
[1,1]       452.7055560068    462
[2,2]       130.4830742030    131
[3,5]       120.6258509819    123
[6,10]      52.2181585555     52
[11,20]     29.5461366439     32
[21,100]    26.2774091990     25
[101,200]   3.6164427260      3
[201,∞]     N/A               65

Writing new variant legend

Writing new haplotype file
[====================] 100%
```


### EXTRACT
Randomly extract a subset of haplotypes (.haps-sample.gz) and output the remaining haplotypes separately (.haps-remainder.gz).
```
options:
  -h, --help            show this help message and exit
  -i INPUT_FILE         Input haplotype file (gzipped)
  -o OUTPUT_FILE        Output haplotype file name
  -s SEED, --seed SEED  Optional seed for reproducibility
  -n NUM                Number of haplotypes to extract
```

```bash
$ python3 -m raresim extract \
    -i example/example.haps.gz \
    -o example/example_subset.haps.gz \
    -n 20000 \
    --seed 3
```


## Complete Workflow Demonstration

For a complete end-to-end workflow demonstrating how to use RAREsim2, see the [RAREsim2_demo](https://github.com/JessMurphy/RAREsim2_demo) repository. This repository demonstrates how to:

- Prepare the required input files
- Perform initial simulations with an over-abundance of rare variants using Hapgen2
- Create datasets for multiple case-control simulation scenarios using RAREsim2
- Perform power analyses for rare variant association methods (Burden, SKAT, SKAT-O)


## Additional Resources

- **Contributing**: See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on contributing to the project
- **GitHub Repository**: [https://github.com/RMBarnard/raresim](https://github.com/RMBarnard/raresim)
- **Issues**: Report bugs or request features at [https://github.com/RMBarnard/raresim/issues](https://github.com/RMBarnard/raresim/issues)
