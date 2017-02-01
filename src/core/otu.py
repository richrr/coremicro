from probability import row_randomize_probability
from functools import total_ordering


def filter_present(vals, min_val):
    """Returns a list of the values in vals that are above the
    specified minimum abundance"""
    return [v for v in vals if v > min_val]


def filter_absent(vals, min_val):
    """Returns a list of values in vals that are at or below the specified
    minimum abundance"""
    return [v for v in vals if v <= min_val]


@total_ordering             # Only have to implement __lt__ and __eq__ for cmp
class Otu:
    def __init__(self, values, name, min_abundance, i_indexes):
        """Construct from a list of ordered values, otu name, minimum abundance
        to be present, and list of indexes of the interest group"""
        # Store the inputs
        self.i_indexes = i_indexes
        self.min_abundance = min_abundance
        self.name = name
        self.values = values
        # Form interest and out groups
        self.i_values = [v for i, v in enumerate(values) if i in i_indexes]
        self.o_values = [v for i, v in enumerate(values) if i not in i_indexes]
        # Whole table stats
        self.total = len(values)
        self.present = len(filter_present(self.values, self.min_abundance))
        self.absent = len(filter_absent(self.values, self.min_abundance))
        self.interest = len(self.i_values)
        # Interest group stats
        self.i_present = len(filter_present(self.i_values, self.min_abundance))
        self.i_absent = len(filter_absent(self.i_values, self.min_abundance))
        self.i_presence_frac = self.i_present / float(self.interest)
        # Out group stats
        self.out = len(self.o_values)
        self.o_present = len(filter_present(self.o_values, self.min_abundance))
        self.o_absent = len(filter_absent(self.o_values, self.min_abundance))
        self.o_presence_frac = self.o_present / float(self.interest)
        # pval
        self.pval = row_randomize_probability(self)

    def __lt__(self, other):
        """Comparison for sorting. Greatest i_presence_frac is smallest, with
        ties broken first by least p-value and finally by first
        (alphabetically) OTU name"""
        return (self.i_presence_frac > other.i_presence_frac or
                self.corrected_pval < other.corrected_pval or
                self.otu < other.otu)

    def __eq__(self, other):
        """Equality for sorting. Equal if equal i_presence_frac, corrected_pval,
        and OTU name (should never happen in use)"""
        return (self.i_presence_frac == other.i_presence_frac and
                self.corrected_pval == other.corrected_pval and
                self.otu == other.otu)

    def __str__(self):
        """Print summary of OTU as TSV line"""
        return '\t'.join(map(str, [self.name, self.pval, self.corrected_pval,
                                   self.i_presence_frac,
                                   self.o_presence_frac]))
