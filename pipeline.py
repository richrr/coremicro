from mapreduce.base_handler import PipelineBase
from mapreduce import mapreduce_pipeline
from mapreduce.input_readers import InputReader
from mapreduce.errors import BadReaderParamsError
from mapreduce.output_writers import OutputWriter

from process_data import shuffle_dicts, convert_shuffled_dict_to_str
from process_data import finish_analysis


class RunRandomPipeline(PipelineBase):

    def run(self, factor, group, out_group, data, mapping, numb_to_run):
        yield mapreduce_pipeline.MapreducePipeline(
            'run_random',
            'process_data.map_random_data',
            'reducer_spec',
            'RandomizedDataInputReader',
            output_writer_spec=None,
            mapper_params={
                'factor': factor,
                'group': group,
                'out_group': out_group,
                'data': data,
                'mapping': mapping,
            },
            reducer_params={
            },
            shards=numb_to_run,
        )

# Mapper: return (frac_thresh, otus) for each thresh for each run
# Reducer: compile list of all otus at each thresh


class RandomizedDataInputReader(InputReader):
    '''
    InputReader for randomized data
    '''

    def __init__(self, factor, group, out_group, data, mapping, done=False):
        super(RandomizedDataInputReader, self).__init__()
        self.factor = factor
        self.group = group
        self.out_group = out_group
        self.data = data
        self.mapping = mapping
        self.done = done

    def next(self):
        '''
        Returns the next input as a dictionary
        '''
        if self.done:
            raise StopIteration()
        self.done = True
        return {
            'factor': self.factor,
            'group': self.group,
            'out_group': self.out_group,
            'data': self.data,
            'mapping': self.mapping,
        }

    @classmethod
    def from_json(cls, input_shard_state):
        '''
        Creates an instance of the InputReader for for the given
        input_shard_state
        '''
        return cls(input_shard_state.get('factor'),
                   input_shard_state.get('group'),
                   input_shard_state.get('out_group'),
                   input_shard_state.get('data'),
                   input_shard_state.get('mapping'),
                   input_shard_state.get('done'))

    def to_json(self):
        '''
        Returns an input_shard_state for the remaining inputs
        '''
        return {
            'factor': self.factor,
            'group': self.group,
            'out_group': self.out_group,
            'data': self.data,
            'mapping': self.mapping,
            'done': self.done
        }

    @classmethod
    def validate(cls, mapper_spec):
        '''
        Validates that the given mapper_spec is ok
        '''
        factor = mapper_spec.params.get('factor')
        group = mapper_spec.params.get('group')
        out_group = mapper_spec.params.get('out_group')
        data = mapper_spec.params.get('data')
        mapping = mapper_spec.params.get('mapping')
        if not factor:
            raise BadReaderParamsError('Factor required')
        if not group:
            raise BadReaderParamsError('Group required')
        if not out_group:
            raise BadReaderParamsError('Out group required')
        if not data:
            raise BadReaderParamsError('Data required')
        if not mapping:
            raise BadReaderParamsError('Mapping required')

    @classmethod
    def split_inputs(cls, mapper_spec):
        '''
        Returns a list of input readers, one per shard; each given a randomized
        mapping
        '''
        factor = mapper_spec.params.get('factor')
        group = mapper_spec.params.get('group')
        out_group = mapper_spec.params.get('out_group')
        data = mapper_spec.params.get('data')
        mapping = mapper_spec.params.get('mapping')
        return [cls(factor, group, out_group, data,
                    convert_shuffled_dict_to_str(shuffle_dicts(mapping)))
                for i in range(mapper_spec.shard_count)]


class RandomDataOutputWriter(OutputWriter):
    def __init__(self, core={}, out={}):
        self.core = core
        self.out = out

    @classmethod
    def validate(cls, mapper_spec):
        '''
        Validates mapper specification
        '''

    @classmethod
    def init_job(cls, mapreduce_state):
        '''
        Initialize the job-level state
        '''
        mapreduce_state.writer_state['core'] = dict()
        mapreduce_state.writer_state['out'] = dict()

    @classmethod
    def finalize_job(cls, mapreduce_state):
        '''
        Finalize the job-level state
        '''
        finish_analysis(key,    # Need to add in key to rest of code!
                        mapreduce_state.writer_state['core'],
                        mapreduce_state.writer_state['out'])

    @classmethod
    def from_json(cls, state):
        '''
        Creates an instance of the OutputWriter from the given json state
        '''
        return cls(state.get('core'),
                   state.get('out'))

    def to_json(self):
        '''
        Returns the writer state to serialize to json
        '''
        return {
            'core': self.core,
            'out': self.out,
        }

    @classmethod
    def create(cls, mr_spec, shard_number, shard_attempt):
        '''
        create new writer for a shard
        '''
        return cls()

    def write(self, data, writer_state):
        '''
        Write data
        '''
        key = data[0]
        otus = data[1]
        group, thresh = key.split(' ')
        if thresh in self['group']:
            self['group'].append(otus)
        else:
            self['group'] = [otus]

    def finalize(self, ctx, shard_state):
        '''
        Finalize writer shard-level state
        '''
        core = shard_state.writer_state['core']
        out = shard_state.writer_state['out']
        for frac in self.core:
            if frac in core:
                core[frac].append(self.core[frac])
            else:
                core[frac] = self.core[frac]
        for frac in self.out:
            if frac in out:
                out[frac].append(self.out[frac])
            else:
                out[frac] = self.out[frac]

    @classmethod
    def get_filenames(cls, mapreduce_state):
        '''
        Returns list of filenames this writer wrote to
        '''
