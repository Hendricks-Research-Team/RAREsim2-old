import random
import copy
from typing import Dict, List
from raresim.common.sparse import SparseMatrix
from raresim.common.legend import Legend


def copy_bin_assignments(bin_assignments):
    return copy.deepcopy(bin_assignments)


def _format_expected(expected) -> str:
    if expected == "N/A":
        return expected
    return f"{float(expected):.3f}"


def _bin_label_from_ranges(lower, upper) -> str:
    return f"[{lower},{upper}]"


def _bin_label_from_standard_bins(bins: list, bin_id: int, total_bins: int) -> str:
    if bin_id < len(bins):
        return _bin_label_from_ranges(bins[bin_id][0], bins[bin_id][1])
    return f"[{str(bins[bin_id - 1][1] + 1)},\u221E]"


def print_bin(bins: list, bin_assignments: dict):
    """
    Backward-compatible bin printer used by existing tests and callers.
    """
    print(f"{'Bin':<12}{'Expected':<12}{'Actual':<10}")
    for bin_id in range(len(bin_assignments)):
        expected = bins[bin_id][2] if bin_id < len(bins) else "N/A"
        print(
            f"{_bin_label_from_standard_bins(bins, bin_id, len(bin_assignments)):<12}"
            f"{_format_expected(expected):<12}"
            f"{len(bin_assignments[bin_id]):<10}"
        )


def print_bin_comparison(bins: list, input_assignments: dict, output_assignments: dict):
    print(f"{'Bin':<12}{'Expected':<12}{'Input':<10}{'Output':<10}")
    for bin_id in range(len(input_assignments)):
        expected = bins[bin_id][2] if bin_id < len(bins) else "N/A"
        print(
            f"{_bin_label_from_standard_bins(bins, bin_id, len(input_assignments)):<12}"
            f"{_format_expected(expected):<12}"
            f"{len(input_assignments[bin_id]):<10}"
            f"{len(output_assignments[bin_id]):<10}"
        )


def summarize_observed_afd(matrix: SparseMatrix) -> dict:
    allele_count_bins = {}
    for row_id in range(matrix.num_rows()):
        allele_count = matrix.row_num(row_id)
        if allele_count == 0:
            continue
        allele_count_bins.setdefault(allele_count, 0)
        allele_count_bins[allele_count] += 1
    return allele_count_bins


def print_observed_afd_comparison(input_summary: dict, output_summary: dict) -> None:
    """
    Print observed allele frequency distributions as a single comparison table.
    """
    print(f"{'Bin':<12}{'Expected':<12}{'Input':<10}{'Output':<10}")
    for allele_count in sorted(set(input_summary) | set(output_summary)):
        print(
            f"[{allele_count},{allele_count}]".ljust(12)
            + f"{'N/A':<12}"
            + f"{input_summary.get(allele_count, 0):<10}"
            + f"{output_summary.get(allele_count, 0):<10}"
        )


def build_probabilistic_bins(matrix: SparseMatrix, legend: Legend) -> list:
    """
    Infer MAC bins from the input matrix by grouping rows with the same keep
    probability.

    Each unique probability is assigned a single inferred MAC range based on
    the non-monomorphic rows carrying that probability.
    """
    grouped = {}
    for row_id in range(matrix.num_rows()):
        row_mac = matrix.row_num(row_id)
        if row_mac <= 0:
            continue

        probability = legend[row_id].get('prob', '.')
        if probability == '.':
            continue

        if probability not in grouped:
            grouped[probability] = {
                'probability': probability,
                'rows': [],
                'lower': row_mac,
                'upper': row_mac,
                'expected': 0.0,
            }

        grouped[probability]['rows'].append(row_id)
        grouped[probability]['lower'] = min(grouped[probability]['lower'], row_mac)
        grouped[probability]['upper'] = max(grouped[probability]['upper'], row_mac)

    inferred_bins = []
    for probability, values in grouped.items():
        values['expected'] = float(probability) * len(values['rows'])
        inferred_bins.append(values)

    return sorted(inferred_bins, key=lambda item: (item['lower'], item['upper'], float(item['probability'])))


def print_probabilistic_bin_summary(probability_bins: list, matrix: SparseMatrix) -> None:
    """
    Print probabilistic pruning bins inferred from the input legend.
    """
    print(f"{'Prob':<10}{'Bin':<12}{'Expected':<12}{'Input':<10}{'Output':<10}")
    for probability_bin in probability_bins:
        output = 0
        for row_id in probability_bin['rows']:
            if matrix.row_num(row_id) > 0:
                output += 1

        bin_label = f"[{probability_bin['lower']},{probability_bin['upper']}]"
        print(
            f"{probability_bin['probability']:<10}"
            f"{bin_label:<12}"
            f"{_format_expected(probability_bin['expected']):<12}"
            f"{len(probability_bin['rows']):<10}"
            f"{output:<10}"
        )


def prune_bins(extra_rows: list, bin_assignments: dict, legend: Legend, matrix: SparseMatrix, bins: list,
               activation_threshold: int, stopping_threshold: int) -> None:
    """
    Prune the rows assigned to each bin so that the bins are more closely aligned with the desired allele counts.

    The following steps are performed:

    1. For each bin, if the actual number of variants in the bin is greater than the desired number of variants by more than 3,
       randomly remove variants from the bin until the difference is less than or equal to 3.

    2. For each bin, if the actual number of variants in the bin is less than the desired number of variants by more than 3,
       randomly add variants from the extra_rows list to the bin until the difference is less than or equal to 3.

    3. Once all bins have been pruned, the bin_assignments dictionary is updated to reflect the new row assignments.

    This method is called by the transform method of the DefaultTransformer.

    Parameters
    ----------
    extra_rows : list
        A list of row ids that are not currently assigned to a bin.
    bin_assignments : dict
        A dictionary of lists where each list contains the row numbers of variants assigned to the corresponding bin.
    legend : Legend
        The Legend object describing the input data.
    matrix : SparseMatrix
        The SparseMatrix object containing the input data.
    bins : list
        A list of bin definitions, each as a tuple of (lower, upper, target)
    activation_threshold : int
        The percentage of the bin size that is required to activate the pruning process.
    stopping_threshold : int
        The percentage of the bin size that is required to stop the pruning process.

    This method is called by the transform method of the DefaultTransformer.

    Returns
    -------
    None
    """
    reserve_pool = bin_assignments[len(bin_assignments) - 1]
    # Loop through the bins from largest to smallest
    for bin_id in range(len(bin_assignments) - 1)[::-1]:
        # If there are any rows with too many 1s for the largest bin, they get put into an extra bin,
        # and we do not prune these alleles away
        if bin_id == len(bins):
            continue

        # How many variants to we need and how many do we have
        need = bins[bin_id][2]
        have = len(bin_assignments[bin_id])

        activation = min([(float(activation_threshold) / 100) * need, 10])
        stop = (float(stopping_threshold) / 100) * need

        # If we have more rows in the bin than what we need, remove some rows
        if have - need > activation:
            # Calculate the probability to remove any given row in the bin
            prob_remove = 1 - float(need) / float(have) if have != 0 else 0
            row_ids_to_rem = []
            for i in range(have):
                # If we hit our stopping threshold to stop the pruning process, then stop the pruning process
                if have - len(row_ids_to_rem) <= need - stop:
                    break
                flip = random.uniform(0, 1)
                # If the row is 'chosen' for removal, remove it and add the row to the list of rows that may be used
                # to make up for not having enough rows in later bins
                if flip < prob_remove:
                    row_id = bin_assignments[bin_id][i]
                    if 'protected' in legend[row_id] and legend[row_id]['protected'] == '1':
                        print(
                            f"WARNING: Attempting to prune a row that is protected. This should not happen and is a bug. RowId = {row_id}, binId={bin_id}")
                    # Add the ith row in the bin to the list of row ids to remove and the list of available
                    # rows to pull from if needed later on
                    row_ids_to_rem.append(row_id)
                    extra_rows.append(row_id)
            for row_id in row_ids_to_rem:
                bin_assignments[bin_id].remove(row_id)

        # If we don't have enough rows in the current bin, pull from the list of excess rows
        elif have < round(need - activation):
            pullFromCommons = False
            if len(extra_rows) < round(abs(need - have)):
                if len(reserve_pool) < round(abs(need - have)):
                    raise Exception(f'ERROR: Current bin has {have} variants, but the model needs {need} variants '
                                    f'and only {len(extra_rows)} excess rows are available to use and prune down '
                                    f'from larger bins and only {len(reserve_pool)} common variants are available '
                                    f'to prune down.')
                else:
                    print(
                        f"WARNING: No more extra pruned rows to pull from. Common variants may be pruned down for bin {bin_id + 1} if necessary")
                    pullFromCommons = True

            # Calculate the probability to use any given row from the available list
            prob_add_from_extras = float(need - have) / float(len(extra_rows)) if len(extra_rows) > 0 else 1
            prob_add_from_reserve = float(need - have - len(extra_rows)) / float(len(reserve_pool))
            row_ids_to_add_from_extras = []
            for i in range(len(extra_rows)):
                flip = random.uniform(0, 1)
                # If the current row is 'chosen' for use, we will prune it down to the desired number of variants
                # and add it back in to the current bin and remove it from the list of available rows to pull from
                if flip < prob_add_from_extras:
                    row_id = extra_rows[i]
                    row_ids_to_add_from_extras.append(row_id)
                    num_to_keep = random.randint(bins[bin_id][0], bins[bin_id][1])
                    num_to_rem = matrix.row_num(row_id) - num_to_keep
                    matrix.prune_row(row_id, num_to_rem)
                    if matrix.row_num(row_id) != num_to_keep:
                        print(
                            f"WARNING: Requested to prune row {row_id} down to {num_to_keep} variants, but after the request the row still has {matrix.row_num} variants. This should not happen and the issue likely needs to be raised with a developer.")
                    bin_assignments[bin_id].append(row_id)

            for row_id in row_ids_to_add_from_extras:
                extra_rows.remove(row_id)

            if pullFromCommons:
                row_ids_to_add_from_reserve = []
                for i in range(len(reserve_pool)):
                    flip = random.uniform(0, 1)
                    if flip < prob_add_from_reserve:
                        row_id = reserve_pool[i]
                        row_ids_to_add_from_reserve.append(row_id)
                        num_to_keep = random.randint(bins[bin_id][0], bins[bin_id][1])
                        num_to_rem = matrix.row_num(row_id) - num_to_keep
                        matrix.prune_row(row_id, num_to_rem)
                        if matrix.row_num(row_id) != num_to_keep:
                            print(
                                f"WARNING: Requested to prune row {row_id} down to {num_to_keep} variants, but after the request the row still has {matrix.row_num} variants. This should not happen and the issue likely needs to be raised with a developer.")
                        bin_assignments[bin_id].append(row_id)
                for row_id in row_ids_to_add_from_reserve:
                    reserve_pool.remove(row_id)


def adjust_for_protected_variants(bins, bin_assignments, legend) -> Dict[int, List]:
    ret = {bin_id: [] for bin_id in range(len(bin_assignments))}
    for bin_id in range(len(bins)):
        rows_in_bin = bin_assignments[bin_id].copy()
        for row_id in rows_in_bin:
            if legend[row_id]['protected'] == '1':
                bin_assignments[bin_id].remove(row_id)
                bins[bin_id][2] -= 1
                ret[bin_id].append(row_id)
    return ret


def add_protected_rows_back(bins, bin_assignments, protected_var_counts_per_bin) -> None:
    for bin_id in protected_var_counts_per_bin:
        if bin_id == len(bins): continue
        bins[bin_id][2] += len(protected_var_counts_per_bin[bin_id])
        bin_assignments[bin_id] += protected_var_counts_per_bin[bin_id]
        bin_assignments[bin_id] = sorted(bin_assignments[bin_id])
