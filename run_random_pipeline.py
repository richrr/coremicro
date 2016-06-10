from mapreduce.base_handler import PipelineBase
from mapreduce import mapreduce_pipeline
from mapreduce.input_readers import InputReader
from mapreduce.errors import BadReaderParamsError
from mapreduce.output_writers import OutputWriter, _get_params


from process_results import ProcessResultsPipeline, perform_sign_calc
from storage import Result_RandomDict, Result_TrueDict, OriginalBiom


class RunRandomPipeline(PipelineBase):
    def run(self, key, num_to_run):
        output = yield mapreduce_pipeline.MapreducePipeline(
            'run_random',
            mapper_spec='process_data.map_random_data',
            reducer_spec='process_results.reduce_random_data',
            input_reader_spec='run_random_pipeline.RandomizedDataInputReader',
            output_writer_spec='run_random_pipeline.RandomizedDataOutputWriter',
            mapper_params={
                'key': key,
            },
            reducer_params={
                'key': key,
            },
            shards=num_to_run,
        )

        yield ProcessResultsPipeline(output)


class RandomizedDataInputReader(InputReader):
    '''
    InputReader for randomized data
    '''

    def __init__(self, key, done=False):
        super(RandomizedDataInputReader, self).__init__()
        self.key = key
        self.done = done

    def next(self):
        '''
        Returns the next input as a dictionary
        '''
        if self.done:
            raise StopIteration()
        self.done = True
        return {
            'key': self.key
        }

    @classmethod
    def from_json(cls, input_shard_state):
        '''
        Creates an instance of the InputReader for for the given
        input_shard_state
        '''
        return cls(input_shard_state.get('key'),
                   input_shard_state.get('done'))

    def to_json(self):
        '''
        Returns an input_shard_state for the remaining inputs
        '''
        return {
            'key': self.key,
            'done': self.done
        }

    @classmethod
    def validate(cls, mapper_spec):
        '''
        Validates that the given mapper_spec is ok
        '''
        key = mapper_spec.params.get('key')
        if not key:
            raise BadReaderParamsError('Key required')

    @classmethod
    def split_input(cls, mapper_spec):
        '''
        Returns a list of input readers, one per shard; each given a randomized
        mapping
        '''
        key = mapper_spec.params.get('key')
        return [cls(key)
                for i in range(mapper_spec.shard_count)]


class RandomizedDataOutputWriter(OutputWriter):
    def __init__(self, key):
        super(RandomizedDataOutputWriter, self).__init__()
        self.key = key

    @classmethod
    def validate(cls, mapper_spec):
        '''
        Validates mapper specification
        '''

    @classmethod
    def from_json(cls, state):
        '''
        Creates an instance of the OutputWriter from the given json state
        '''
        return cls(state.get('key'))

    def to_json(self):
        '''
        Returns the writer state to serialize to json
        '''
        return {
            'key': self.key
        }

    @classmethod
    def create(cls, mr_spec, shard_number, shard_attempt, _writer_state=None):
        '''
        create new writer for a shard
        '''
        return cls(_get_params(mr_spec.mapper).get('key'))

    def write(self, data):
        '''
        Write data
        '''
        print 'Writing data'
        category = data[0]
        compiled = data[1]
        (user_args, to_email, p_val_adj, DELIM, NTIMES,
         otu_table_biom, g_info_list, factor, group, out_group,
         OUTPFILE, categ_samples_dict) = OriginalBiom.get_params(self.key)
        if category == 'core':
            core_true = Result_TrueDict.get_entry(self.key).true_results
            results = perform_sign_calc(self.key, compiled, p_val_adj,
                                        DELIM, core_true, NTIMES)
            Result_RandomDict.add_entry(self.key, results)
        elif category == 'out':
            out_true = Result_TrueDict.get_entry(self.key,
                                                 out_group=True).true_results
            results = perform_sign_calc(self.key, compiled, p_val_adj,
                                        DELIM, out_true, NTIMES)
            Result_RandomDict.add_entry(self.key, results, out_group=True)

    def finalize(self, ctx, shard_state):
        '''
        Finalize writer shard-level state
        '''
        pass

    @classmethod
    def get_filenames(cls, mapreduce_state):
        '''
        Returns list of filenames this writer wrote to
        '''
        return _get_params(mapreduce_state.mapreduce_spec.mapper).get('key')
