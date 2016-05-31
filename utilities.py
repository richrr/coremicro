from storage import OriginalBiom


# get unique elements from last column (otus)
def compile_results(otus, DELIM):
    taxon_list = list()         # this may be a json or list
    result_ = otus              # json.loads(otus)
    for l in result_:
        l = l.strip()
        contents = l.split(DELIM)
        if '#' in l or not l:
            continue
        taxon_list.append(contents[1])
    return list(set(taxon_list))
