'''
This is a modest collection of methods that I regularly reuse.
'''

import sys
#import networkx as nx


######################################################################
# READ FILE, ARG: FILENAME
# RETURN LIST
#####################################################################
def read_file(filename):
    f = open(filename)
    lines = f.readlines()
    f.close()
    print 'Read %d lines from file %s ' %(len(lines), filename)
    return lines


#####################################################################
# WRITE TXT TO FILE, ARG: TXT , FILENAME 
####################################################################
def writeTXT_to_file(txt, filename):
    outfile = open(filename,'w')
    outfile.write("%s\n" % txt)
    outfile.close()
    print 'TXT written to file %s' % (filename)
    return

#####################################################################
# WRITE LIST TO FILE, ARG: LIST , FILENAME 
####################################################################
def writeLIST_to_file(list, filename):
    outfile = open(filename,'w')
    for l in list:
        outfile.write("%s\n" % l)
    outfile.close()
    print 'LIST written to file %s' % (filename)
    return


#####################################################################
# WRITE DICTIONARY TO FILE, ARG: DICTIONARY , FILENAME
######################################################################
def writeDICT_to_file(DICT, filename):
    outfile = open(filename,'w')
    for k, v in DICT.items():
        outfile.write("%s\t%s\n" % (k , v))
    outfile.close()
    print 'DICT written to file %s' % (filename)
    return

######################################################################
# CALC HISTOGRAM OF KEYS WITH # OF VALUES, 
# ARG: DICT , OPTIONAL ARG: strings for key and value
######################################################################
def calc_hist(DICT, k='KEYS' , v='VALUES'):
    hist = {}  # counts the number of kegg ids with particular number of values
    for key in DICT:
        l = len(DICT[key])
        if l not in hist:
            hist[l] = 0
        hist[l]+=1
    for key in sorted(hist):
        print '%d %s have %d %s entries' % (hist[key], k, key, v)
    print '\n'
    return

######################################################################
# CONVERT LIST TO DICTIONARY, 
# ARG: LIST
# OPTIONAL ARG: DELIMITER FOR SPLITTING INTO KEY & VALUE, JOINER for multiple mappings , string for filename & 0-indexed columns to use as key and value
######################################################################
def list_to_dict(List, delim='\t', joiner=',', string="current", key_column=0, val_column=1):
    local_dict = dict()
    multiple_mappings = 'no'
    for l in List:
        if '#' in l or not l: # checks for commented lines and empty lines
            continue
        key = value = ''
        array = l.strip().split(delim)
        #print array
        key = array[key_column]
        value = array[val_column]
        if key in local_dict:
            #print '%s key already exists in %s file. Adding %s to its value in dict ' %(key, string, value)
            multiple_mappings = ''
            current_values = list()
            current_values.append(local_dict[key])
            current_values.append(value)
            current_values = condense_list(current_values , joiner)
            local_dict[key] = joiner.join(current_values)
        else:
            local_dict[key] = value
    print "There is %s multiple mappings, multiple values to the key are separated by %s" % (multiple_mappings, joiner)
    return local_dict


######################################################################
# MAP IDENTIFIERS: TWO FILES into LIST, 
# ARG: TWO FILES (all lines from the first file (data file) are printed, the second file (annotation file) is used for making dict)
# OPTIONAL ARG: DELIMITER FOR SPLITTING, JOINER string , key and value columns, insert newid or replace oldid with newid, 
#               0-indexed column to where value is to be inserted, keep lines with missing value from file 2
'''
Map the identifiers in various columns of a tab-delimited file to a different
namespace, using a supplied mapper. Columns are 0-indexed.
comments begin with '#'

This script takes a single tab-delimited file with one or more columns. This is
the file to be mapped:
--- This is matrix.txt ---
  id1_n1   val1,1   val1,2
  id2_n1   val2,1   val2,2
  id3_n1   val3,1   val3,2

The script also takes a namespace mapper, which is a multi-column,
tab-delimited file. Each row corresponds to a single entity, and each column
contains IDs for that entity from a different namespace. For example:
--- This is mapper.txt ---
  #namespace1  namespace2  namespace3
  id1_n1       id1_n2      id1_n3
  id2_n1       id2_n2      id2_n3
  id3_n1       id3_n2      id3_n3

Lastly, the script requires two integers that specify the column to be used as
the "from" namespace (key) and the column to be used as the "to" namespace (value).

Example: map-identifiers (matrix.txt, mapper.txt, '\t', '\t', 0, 2)

  # comments begin with '#'
  id1_n3   val1,1   val1,2
  id2_n3   val2,1   val2,2
  id3_n3   val3,1   val3,2 
'''


######################################################################
def map_identifiers(file1, file2, delim='\t', joiner='\t', key_column=0, val_column=1, insert=False, insert_column=1, skipmissing=False):
    f1 = read_file(file1) # treated as a list

    f2 = read_file(file2) # used to make a dictionary
    d2 = list_to_dict(f2, delim, joiner, "current", key_column, val_column)

    newlist = list()

    for l in f1: 
        if '#' in l:
            newlist.append(l.strip())
            continue
        array = l.strip().split(delim)
        # insert newid
        if array[0] in d2 and insert:
            array.insert(insert_column , d2[array[0]])
        # replace id with newid 
        elif array[0] in d2 and not insert:
            array[0] = d2[array[0]]
        # skip line if missing value
        elif array[0] not in d2 and skipmissing:
            continue
        # the following if are executed only if skipmissing is False
        # insert missing value
        elif array[0] not in d2 and insert:
            array.insert(insert_column , 'NA')
        # replace id with missing value
        elif array[0] not in d2 and not insert:
            array[0] = 'NA'
        
        newlist.append(joiner.join(array))

    return newlist


######################################################################
# splits the elements of the list with requested delimiter, 
# returns list with unique elements
######################################################################
def condense_list(List , delim=','):
    condensed_list = list()
    for x in List:
        condensed_list.extend(x.split(delim))
    return sorted(list(set(condensed_list)))



def readDict(f, fromCol=1, toCol=2, sep='\t'):
    '''
    Read the dict from the given tab-delimited file. The dict
    maps items in the fromCol to items in the toCol (1-based column index).
    '''
    itemMap = {}
    for line in open(f, 'r').readlines():
        if line=='':
            continue
        if line[0]=='#':
            continue
        items = line.rstrip().split(sep)
        if len(items)<max(fromCol, toCol):
            continue
        key = items[fromCol-1]
        val = items[toCol-1]
        if key=='':
            continue
        if val=='':
            continue
        itemMap[key] = val
    return itemMap

def readColumnsSep(f, sep='\t', *cols):
    '''
    Read multiple columns and return the items from those columns
    in each line as a tuple.
    
    foo.txt:
        a b c
        d e f
        g h i
        
    Calling "readColumnsSep('foo.txt', '\t', 0, 2)" will return:
        [(a, c), (d, f), (g, i)]
    
    '''
    if len(cols)==0:
        return []
    rows = []
    for line in open(f, 'r').readlines():
        if line=='':
            continue
        if line[0]=='#':
            continue
        items = line.rstrip().split(sep)
        if len(items)<=max(cols):
            continue
        rows.append(tuple([items[c] for c in cols]))
    return rows



def readItemList(f, col=0, sep='\t'):
    '''
    Read the given column of the tab-delimited file f
    and return it as a list. Col is the 0-based column index.
    '''
    itemlist = []
    for line in open(f, 'r').readlines():
        if line=='':
            continue
        if line[0]=='#':
            continue
        items = line.rstrip().split(sep)
        if len(items)<=col:
            continue
        itemlist.append(items[col])
    return itemlist

def readItemSet(f, col=0, sep='\t'):
    '''
    Read the given column of the tab-delimited file f
    and return it as a set. Col is the 0-based column index.
    
    A wrapper to readItemList, returning a set instead of a list.
    '''
    return set(readItemList(f, col=col, sep=sep))



def jaccardIndex(A, B, step=1):
    '''
    Compute the Jaccard index between A and B at intervals of the
    given step size.
    
    Returns a list of tuples (x,y) where y is the JI at cutoff x.
    '''
    seenOnce = set()
    seenTwice = set()
    minRange = min(len(A), len(B))

    if step<1:
        step = 1
    
    #JIvalues stores a list of pairs [x,y], where 'y' is the JI at cutoff 'x'
    JIvalues = []
    for k in range(minRange):
        currA = A[k]
        currB = B[k]
        if currA in seenOnce:
            seenOnce.remove(currA)
            seenTwice.add(currA)
        elif currA not in seenTwice: #so currA hasn't been seen at all
            seenOnce.add(currA)
        if currB in seenOnce:
            seenOnce.remove(currB)
            seenTwice.add(currB)
        elif currB not in seenTwice: #so currB hasn't been seen at all
            seenOnce.add(currB)

        #get the various JI values
        if ((k+1)%step==0) or (k+1==minRange):
            currJI = 1.0*len(seenTwice)/(len(seenOnce)+len(seenTwice))
            JIvalues.append( (k+1, currJI) )
    return JIvalues


def generalKT(A, B, step=1):
    '''
    Compute the generalization fo Kentall's tau for top-k lists by Fagin, Kumar,
    and Sivakumar (http://epubs.siam.org/doi/abs/10.1137/S0895480102412856).
    Briefly, this method compute the number of mismatched pair ordering between
    the two lists.
    
    This method is explained nicely in Mccown and Nelson
    (http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.65.3888).
    '''
    mapA = {}
    mapB = {}
    for i,a in enumerate(A):
        mapA[a] = i
    for i,b in enumerate(B):
        mapB[b] = i
    
    minRange = min(len(A), len(B))

    if step<1:
        step = 1
    
    #vals stores a list of pairs [x,y], where 'y' is the generalized Kendall tau at cutoff 'x'
    vals = []
    mismatched = set()
    for k in range(minRange):
        a = A[k]
        b = B[k]
        
        rnkaInB = mapB.get(a, minRange+1)
        rnkbInA = mapA.get(b, minRange+1)
        
        for rnkTempbInB, tempb in enumerate(B[:k+1]):
            rnkTempbInA = mapA.get(tempb, minRange+1)
            if (k<rnkTempbInA and rnkaInB>rnkTempbInB) or (k>rnkTempbInA and rnkaInB<rnkTempbInB):
                mismatched.add( tuple(sorted((a, tempb))) )
        
        # only index to k here, since we compare a to b in the previous for loop
        # by indexing to k+1        
        for rnkTempaInA, tempa in enumerate(A[:k]):
            rnkTempaInB = mapB.get(tempa, minRange+1)
            if (k<rnkTempaInB and rnkbInA>rnkTempaInA) or (k>rnkTempaInB and rnkbInA<rnkTempaInA):
                mismatched.add( tuple(sorted((tempa,b))) )
        
        #update vals
        if ((k+1)%step==0) or (k+1==minRange):
            numpairs = (k+1)*(k+1)
            currVal = 1.0-1.0*len(mismatched)/numpairs
            vals.append( (k+1, currVal) )
    return vals


def computePR(pos, neg, values):
    '''
    pos and neg are sets of positive and negative items, respectively, and
    values is a list of pairs (item, val) where items are ranked in increasing
    order of their val. Precision and recall are computed for every unique val
    in the list of (item,val) pairs, and PR values are returned as a list of
    pairs (precision, recall). 
    '''
    pr = []
    tp = 0.0
    fp = 0.0
    lastVal = None
    for item, val in sorted(values, key=lambda x: x[1]):
        if item in pos:
            tp += 1.0
        elif item in neg:
            fp += 1.0
        else:
            continue
        
        if tp==0 and fp==0:
            continue
        
        precision = tp/(tp+fp)
        recall = tp/len(pos)
        if lastVal==None:
            pr.append( (precision, recall) )
        # if the val is the same as the lastVal, update the last (precision, recall) entry
        elif lastVal==val:
            pr[-1] = ( (precision, recall) )
        # only add a new entry if it is different from the previous
        elif pr[-1]!=(precision, recall):
            pr.append( (precision, recall) )
        lastVal = val
    return pr
