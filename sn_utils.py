#!/usr/bin/env python3

import numpy as np

from Bio import SeqIO


def read_fastafile(fasta_file):
    with open(fasta_file,'r') as fa_fh:
        #fa_dict = SeqIO.to_dict(SeqIO.parse(fa_fh, "fasta"))
        fa_dict = seqrcds_to_ordereddict(SeqIO.parse(fa_fh, "fasta")) 
    return fa_dict

def read_fasta(fasta_text):
    fa_dict = SeqIO.to_dict(SeqIO.parse(fasta_text, "fasta"))
    return fa_dict

def read_tomlfile(toml_file):
    from collections import OrderedDict
    import toml
    with open(toml_file) as tomlfh:
        toml_d = toml.loads(tomlfh.read(), OrderedDict) #取り出す時の順番が毎回同じ方が助かる場合があるから。tomlファイルに記載された順番どおりではない。
    return toml_d

def random_sample_fasta(fa_d, n):
    import random
    fa_sampled_l = random.sample(list(fa_d.items()), n)
    fa_sampled_d = dict()
    for seq in fa_sampled_l: #データ形式をSeqIOで読み込んだものと同じに戻す
        fa_sampled_d.update({seq[0]: seq[1]})
    return fa_sampled_d

def read_fasta_subtype(fasta_text):
    from collections import defaultdict
    fasta_d = read_fasta(fasta_text)
    subtypes_d = defaultdict(dict)
    for k, v in fasta_d.items():
        subtype = k.split(":")[-1].split("-")[0] #ドメイン等の領域が書いてある場合もあるから.
        subtypes_d[subtype].update({k:v})
    return subtypes_d

def seqrcds_to_dict(sequences, key_function=None):  #from BioPython (SeqIO)
    """Turns a sequence iterator or list into a dictionary. 
    
        - sequences  - An iterator that returns SeqRecord objects, 
          or simply a list of SeqRecord objects. 
        - key_function - Optional callback function which when given a 
          SeqRecord should return a unique key for the dictionary. 
    
    e.g. key_function = lambda rec : rec.name 
    or,  key_function = lambda rec : rec.description.split()[0] 
    
    If key_function is omitted then record.id is used, on the assumption 
    that the records objects returned are SeqRecords with a unique id. 
    
    If there are duplicate keys, an error is raised. 
    """
    if key_function is None: 
          key_function = lambda rec: rec.id 
    
    d = dict() 
    for record in sequences: 
        key = key_function(record) 
        if key in d: 
            raise ValueError("Duplicate key '%s'" % key) 
        d[key] = record 
    return d 

def seqrcds_to_ordereddict(sequences, key_function=None):  #modified to return as OrderedDict
    from collections import OrderedDict
    if key_function is None: 
          key_function = lambda rec: rec.id 
    
    d = OrderedDict()
    for record in sequences: 
        key = key_function(record) 
        if key in d: 
            raise ValueError("Duplicate key '%s'" % key) 
        d[key] = record 
    return d 

def draw_annopos(ax, anno_dict, rows=3, readingframe=False, fs=9):
    """
    anno_dict = {name:[start,end]}
    """
    from matplotlib.patches import Rectangle
    y1 ,height, pad = 0, 1, 0.2
    ax.set_ylim([-pad,rows*(height+pad)])
    anno_elements = []
    for name, x in anno_dict.items():
        anno_elements.append({'name': name,
                     'x1': x[0], 'x2': x[1], 'width': x[1] - x[0]})
    anno_elements.sort(key = lambda x:x['x1'])
    for ai, anno in enumerate(anno_elements):
        if readingframe:
            anno['y1'] = y1 + (height + pad) * (2 - (anno['x1'])%3)
        else:
            anno['y1'] = y1 + (height + pad) * (ai%rows)
        anno['y2'] = anno['y1'] + height

    for anno in anno_elements:
        r = Rectangle((anno['x1'], anno['y1']),
                      anno['width'],
                      height,
                      facecolor=[0.8] * 3,
                      edgecolor='k',
                      label=anno['name'])
        
        xt = anno['x1'] + 0.5 * anno['width']
        yt = anno['y1'] + 0.2 * height + height*(anno['width']<500)

        ax.add_patch(r)
        ax.text(xt, yt,
                anno['name'],
                color='k', 
                fontsize=fs,
                ha='center')

#original -https://github.com/neherlab/HIVEVO_figures/blob/master/src/util.py
def draw_genome(ax, annotations,rows=3, readingframe=True,fs=9):
    from matplotlib.patches import Rectangle
    y1 ,height, pad = 0, 1, 0.2
    ax.set_ylim([-pad,rows*(height+pad)])
    anno_elements = []
    for name, feature in annotations.iteritems():
        x = [feature.location.nofuzzy_start, feature.location.nofuzzy_end]
        anno_elements.append({'name': name,
                     'x1': x[0], 'x2': x[1], 'width': x[1] - x[0]})
    anno_elements.sort(key = lambda x:x['x1'])
    for ai, anno in enumerate(anno_elements):
        if readingframe:
            anno['y1'] = y1 + (height + pad) * (2 - (anno['x1'])%3)
        else:
            anno['y1'] = y1 + (height + pad) * (ai%rows)
        anno['y2'] = anno['y1'] + height

    for anno in anno_elements:
        r = Rectangle((anno['x1'], anno['y1']),
                      anno['width'],
                      height,
                      facecolor=[0.8] * 3,
                      edgecolor='k',
                      label=anno['name'])
        
        xt = anno['x1'] + 0.5 * anno['width']
        yt = anno['y1'] + 0.2 * height + height*(anno['width']<500)

        ax.add_patch(r)
        ax.text(xt, yt,
                anno['name'],
                color='k', 
                fontsize=fs,
                ha='center')

#移動平均
def running_average_masked(obs, ws, min_valid_fraction=0.95):
    '''
    calculates a running average via convolution, fixing the edges
    obs     --  observations (a masked array)
    ws      --  window size (number of points to average)

    https://github.com/neherlab/HIVEVO_figures/blob/master/src/evolutionary_rates.py
    '''
    #tmp_vals = np.convolve(np.ones(ws, dtype=float), obs*(1-obs.mask), mode='same')
    tmp_vals = np.convolve(np.ones(ws, dtype=float), obs.filled(0), mode='same')

     # if the array is not masked, edges needs to be explictly fixed due to smaller counts
    if len(obs.mask.shape) == 0:
        tmp_valid = ws*np.ones_like(tmp_vals)
        # fix the edges. using mode='same' assumes zeros outside the range
        if ws%2==0:
            tmp_vals[:ws//2]*=float(ws)/np.arange(ws//2,ws)
            if ws//2>1:
                tmp_vals[-ws//2+1:]*=float(ws)/np.arange(ws-1,ws//2,-1.0)
        else:
            tmp_vals[:ws//2]*=float(ws)/np.arange(ws//2+1,ws)
            tmp_vals[-ws//2:]*=float(ws)/np.arange(ws,ws//2,-1.0)

    # if the array is masked, then we get the normalizer from counting the unmasked values
    else:
        tmp_valid = np.convolve(np.ones(ws, dtype=float), (1-obs.mask), mode='same')

    run_avg = np.ma.array(tmp_vals / tmp_valid)
    run_avg.mask = tmp_valid < ws * min_valid_fraction

    return run_avg
