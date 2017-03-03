import numpy
from functools import total_ordering


def present(values, min_abundance):
    """Replaces all non-present values with None"""
    return [v if v > min_abundance else None for v in values]


def absent(values, min_abundance):
    """Replaces all non-absent values with None"""
    return [v if v <= min_abundance else None for v in values]


def interest(values, i_indexes):
    """Replaces all non-interest values with None"""
    return [v if i in i_indexes else None for i, v in enumerate(values)]


def out(values, i_indexes):
    """Replaces all non-out values with None"""
    return [v if i not in i_indexes else None for i, v in enumerate(values)]


def mean(values):
    """Calculate the mean of the given list of values. None values are treated
    as not existing"""
    return numpy.mean([v for v in values if v is not None])


def count(values):
    """Number of non-None values in values"""
    return len([v for v in values if v is not None])


def standard_error(values):
    """Standard error of the given list of values, treating None values as not
    existing"""
    return (numpy.std([v for v in values if v is not None], ddof=1) /
            numpy.sqrt(len(values)))


@total_ordering             # Only have to implement __lt__ and __eq__ for cmp
class Otu:
    def __init__(self, values, name, min_abundance, i_indexes):
        """Construct from a list of ordered values, otu name, minimum abundance
        to be present, and list of indexes of the interest group"""
        # Store the inputs
        self.name = name
        self.i_indexes = i_indexes
        self.min_abundance = min_abundance
        self.values = values

        interest_values = interest(values, i_indexes)
        out_values = out(values, i_indexes)
        self.interest_present = count(present(interest_values, min_abundance))
        self.interest_absent = count(absent(interest_values, min_abundance))
        self.interest_mean = mean(interest_values)
        self.interest_error = standard_error(interest_values)
        self.out_present = count(present(out_values, min_abundance))
        self.out_absent = count(absent(out_values, min_abundance))
        self.out_mean = mean(out_values)
        self.out_error = standard_error(out_values)

        self.interest = self.interest_present + self.interest_absent
        self.out = self.out_present + self.out_absent
        self.present = self.interest_present + self.out_present
        self.absent = self.interest_absent + self.out_absent
        self.total = self.interest + self.out

        self.interest_frac = self.interest_present / float(self.interest)
        self.out_frac = self.out_present / float(self.out)

    def __lt__(self, other):
        """Comparison for sorting. Greatest interest presence frac is smallest,
        with ties broken first by least p-value and finally by first
        (alphabetically) OTU name. Corrected_pval must have been calculated
        first"""
        return (self.interest_frac > other.interest_frac or
                self.corrected_pval < other.corrected_pval or
                self.name < other.name)

    def __eq__(self, other):
        """Equality for sorting. Equal if equal interest presence frac,
        corrected_pval, and OTU name (should never happen in use)
        Corrected_pval must have been calculated first"""
        return (self.interest_frac == other.interest_frac and
                self.corrected_pval == other.corrected_pval and
                self.name == other.name)
