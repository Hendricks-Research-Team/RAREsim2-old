from raresim.common.exceptions import IllegalArgumentException

class RunConfig:
    def __init__(self, args):
        self.args = args
        self.run_type: str = self.__determine_run_type()
        self.remove_zeroed_rows: bool = args.z
        self.small_sample: bool = args.small_sample
        self.is_probabilistic: bool = getattr(args, 'prob', False)
        self.activation_threshold: int = args.activation_threshold
        self.stop_threshold: int = args.stop_threshold

    def __determine_run_type(self) -> str:
        if self.args.exp_bins is None and not getattr(self.args, 'prob', False):
            if self.args.exp_fun_bins is not None \
                    and self.args.exp_syn_bins is not None:
                return "func_split"
            elif self.args.fun_bins_only is not None:
                return "fun_only"
            elif self.args.syn_bins_only is not None:
                return "syn_only"
            else:
                raise IllegalArgumentException('If variants are split by functional/synonymous ' +
                                               'files must be provided for --functional_bins ' +
                                               'and --synonymous_bins')
        elif getattr(self.args, 'prob', False):
            return "probabilistic"
        else:
            return "standard"
