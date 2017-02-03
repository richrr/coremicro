from functools import total_ordering


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

    def child(self, values):
        "Returns a copy of this with it's values replaced with the given"
        return Otu(values, self.name, self.min_abundance, self.i_indexes)

    def present(self):
        "Replaces all non-present values with None"
        return self.child([v if v > self.min_abundance else None
                           for v in self.values])

    def absent(self):
        "Replaces all non-absent values with None"
        return self.child([v if v <= self.min_abundance else None
                           for v in self.values])

    def interest(self):
        "Replaces all non-interest values with None"
        return self.child([v if i in self.i_indexes else None
                           for i, v in enumerate(self.values)])

    def out(self):
        "Replaces all non-out values with None"
        return self.child([v if i not in self.i_indexes else None
                           for i, v in enumerate(self.values)])

    def count(self):
        return len([v for v in self.values if v is not None])

    def frac(self):
        "Returns the fraction of values that are present"
        return (self.present().count()) / float(self.count())

    def __lt__(self, other):
        """Comparison for sorting. Greatest interest presence frac is smallest,
        with ties broken first by least p-value and finally by first
        (alphabetically) OTU name. Corrected_pval must have been calculated
        first"""
        return ((self.interest().frac() >
                 other.interest().frac()) or
                self.corrected_pval < other.corrected_pval or
                self.name < other.name)

    def __eq__(self, other):
        """Equality for sorting. Equal if equal interest presence frac,
        corrected_pval, and OTU name (should never happen in use)
        Corrected_pval must have been calculated first"""
        return ((self.interest().presence_frac() ==
                 other.interest().presence_frac()) and
                self.corrected_pval == other.corrected_pval and
                self.name == other.name)

    def __str__(self):
        """Print summary of OTU as TSV line. pval and corrected_pval must have
        been added first"""
        return '\t'.join(map(str, [self.name, self.pval, self.corrected_pval,
                                   self.interest().frac(),
                                   self.out().frac()]))
