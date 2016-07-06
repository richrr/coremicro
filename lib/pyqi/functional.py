#!/usr/bin/env python

from pyqi.core.interface import Interface, InterfaceOption
from testing import FilterSamples, FilterObservations
from pyqi.core.command import Command
from types import NoneType
from pyqi.core.factory import general_factory
from biom.table import Table as BIOMTable
from biom.parse import parse_biom_table

class Functional(Interface):
    def _validate_usage_examples(self, foo):
        pass

    def _the_in_validator(self, in_):
        if isinstance(in_, NoneType):
            raise ValueError("wtf")

    def _input_handler(self, in_, *args, **kwargs):
        cmd_in = self._get_inputs()[0]
        if args:
            raise ValueError("wtf")
        if not isinstance(in_, cmd_in.Type):
            raise ValueError("wtf2")

        # assume kwargs values are all in memory already
        kwargs[cmd_in.Name] = in_
        return kwargs

    def _the_out_validator(self, result):
        pass

    def _output_handler(self, result):
        cmd_out = self._get_outputs()[0]
        return result[cmd_out.Name]

class FunctionalOption(InterfaceOption):
    def __init__(self, **kwargs):
        super(FunctionalOption, self).__init__(**kwargs)
        self.Type = self.Parameter.DataType

filter_samples_inputs = [FunctionalOption(Parameter=FilterSamples.CommandIns['table'])]
filter_samples_outputs = [FunctionalOption(Parameter=FilterSamples.CommandOuts['filtered'])]


filter_obs_inputs = [FunctionalOption(Parameter=FilterSamples.CommandIns['table'])]
filter_obs_outputs = [FunctionalOption(Parameter=FilterSamples.CommandOuts['filtered'])]

f = general_factory(FilterSamples, None, filter_samples_inputs, filter_samples_outputs, None, Functional)()
g = general_factory(FilterObservations, None, filter_obs_inputs, filter_obs_outputs, None, Functional)()

t1 = parse_biom_table(open('/Users/mcdonadt/ResearchWork/software/biom-format/examples/rich_dense_otu_table.biom'))
t2 = parse_biom_table(open('/Users/mcdonadt/ResearchWork/software/biom-format/examples/min_dense_otu_table.biom'))
t3 = parse_biom_table(open('/Users/mcdonadt/ResearchWork/software/biom-format/examples/rich_sparse_otu_table.biom'))

res = g(f("asdasdasd", seqs_per_sample=3), seqs_per_observation=1)
res2 = map(lambda x: f(x, seqs_per_sample=3), [t1, t2, t3])
res3 = reduce(lambda x,y: x.merge(f(y, seqs_per_sample=3)), [t1,t2,t3])

print res
for i in res2:
    print i
print "***"
print res3
