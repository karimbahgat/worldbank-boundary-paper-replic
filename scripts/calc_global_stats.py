# imports
import boundarytools
import numpy as np

import os
import sys
import json
from urllib.request import urlopen
import csv
import traceback
import multiprocessing as mp
import datetime
from time import time
import itertools

# params

OUTPUT_DIR = 'global_stats'
REPLACE = False
IGNORE_SOURCES = []
MAXPROCS = 4


def loop_country_levels():
    url = 'https://raw.githubusercontent.com/wmgeolab/geoContrast/stable/releaseData/geoContrast-meta.csv'
    raw = urlopen(url).read().decode('utf8')
    reader = csv.DictReader(raw.split('\n'))
    def key(row):
        return row['boundaryISO'], int(row['boundaryType'][-1])
    iso_levels = sorted(set([key(row) for row in reader]))
    for iso,level in iso_levels:
        yield iso,level

def get_country_level_stats(iso, level):
    #results = {}

    # open areas
    try:
        path = 'global_relations/{}-ADM{}-areas.json'.format(iso, level)
        with open(path) as r:
            areas = json.loads(r.read())
    except:
        return None

    # open relations
    try:
        path = 'global_relations/{}-ADM{}-relations.json'.format(iso, level)
        with open(path) as r:
            relations = json.loads(r.read())
    except:
        return None

    # collect pairwise source similarities
    source_similarities = {}
    for src,relations2 in relations.items():
        As = areas[src]
        source_similarities_row = {}
        for src2,featurepairs in relations2.items():
            print('')
            print(src,src2)
            Bs = areas[src2]

            # for each feat add various area measurements
            mainArea = sum(As)
            comparisonArea = sum(Bs)
            isecs = [AB for i1,i2,(Adiff,Bdiff,AB) in featurepairs]
            isecArea = sum(isecs)
            print(f'A {mainArea} B {comparisonArea} isec {isecArea}')

            # decompose into Adiff, Bdiff, and union
            mainDiffArea = (1 - (isecArea / mainArea)) * mainArea
            comparisonDiffArea = (1 - (isecArea / comparisonArea)) * comparisonArea
            unionArea = mainDiffArea + comparisonDiffArea + isecArea
            print(f'Adiff {mainDiffArea} Bdiff {comparisonDiffArea} union {unionArea}')

            # calc perc of matching
            matches = match_features(As, Bs, featurepairs)
            match_areas = [stats['within'] * As[i1]
                           for i1,i2,stats in matches if stats]
            matchArea = sum(match_areas)
            percArea = matchArea / unionArea * 100
            source_similarities_row[src2] = percArea
            print(f'match {matchArea} perc {percArea}')

        source_similarities[src] = source_similarities_row

    # return dict of source pairs, each with a single probability/similarity metric
    return source_similarities

def match_features(As, Bs, relations):
    # this function returns a list of matches for a particular source combination
    # each list entry is the match of each feature with the most similar feature in the other source

    # create relations lookup dict
    key = lambda x: x[0] # group by row index
    relations_lookup = dict([(k,list(group)) for k,group in itertools.groupby(sorted(relations, key=key), key=key)])

    # calc similarities
    similarities = []
    for i,A in enumerate(As):
        related = relations_lookup.get(i, []) # if no related then prob = 0
        # calc similarity of feature A with all features in the other source
        A = As[i] # could also be Adiff+AB

        simils = []
        for _,i2,(Adiff,Bdiff,AB) in related:
            AorB = sum([Adiff,Bdiff,AB])
            B = Bs[i2] # could also be Bdiff+AB
            if A == 0 or B == 0:
                # rare freak case of 0-area geoms
                continue
            equality = AB / AorB
            within = AB / A
            stats = {'equality':equality, 'within':within}
            simils.append((i, i2, stats))

        if simils:
            best = sorted(simils, key=lambda x: x[-1]['equality'])[-1]
            similarities.append(best)
        else:
            similarities.append((i,None,None))

    # do a second pass so that a match is only kept if it's the best in both sources
    similarities2 = []
    for i,i2,stats in similarities:
        similarities2.append((i,i2,stats))
        if i2 is None:
            continue
        othermatches = [(_i,_i2,_stats) for _i,_i2,_stats in similarities if _i2 == i2]
        assert len(othermatches) >= 1
        bestmatch = sorted(othermatches, key=lambda x: x[-1]['equality'])[-1]
        if i == bestmatch[0]:
            # this is the best match, keep it
            pass
        else:
            # this is not the best match, set match to null
            similarities2[i] = (i, None, None)
    
    return similarities2


def process(iso, level):
    # calc stats
    stats = get_country_level_stats(iso, level)
    if not stats:
        # relations havent been calculated, skip
        return

    # dump stats as json
    filename = '{}-ADM{}-stats.json'.format(iso, level)
    path = os.path.join(OUTPUT_DIR, filename)
    with open(path, 'w', encoding='utf8') as w:
        json.dump(stats, w, indent=4)

    # dump errors as json if any
    if 0: #errors:
        filename = '{}-ADM{}-errors.json'.format(iso, level)
        path = os.path.join(OUTPUT_DIR, filename)
        with open(path, 'w', encoding='utf8') as w:
            w.write(json.dumps(errors, indent=4))


def process_logger(logfile, func, **kwargs):
    logger = open(logfile, 'w', encoding='utf8')
    sys.stdout = logger
    sys.stderr = logger
    print('PID:',os.getpid())
    print('time',datetime.datetime.now().isoformat())
    print('working path',os.path.abspath(''))
    # run it
    func(**kwargs)
    # finish
    print('finished!')
    print('time',datetime.datetime.now().isoformat())


if __name__ == '__main__':
    maxprocs = MAXPROCS
    procs = []

    # loop country levels
    for iso,level in loop_country_levels():
        print('')
        print(iso,level)

        # skip if stats already exists
        filename = '{}-ADM{}-stats.json'.format(iso, level)
        path = os.path.join(OUTPUT_DIR, filename)
        if not REPLACE and os.path.lexists(path):
            print('already exists, skipping')
            continue

        # local
        #process(iso, level)
        #continue

        # multiprocessing
        logfile = '{}-ADM{}-log.txt'.format(iso, level)
        logpath = os.path.join(OUTPUT_DIR, logfile)
        p = mp.Process(target=process_logger,
                       args=[logpath,process],
                       kwargs={'iso':iso, 'level':level})
        p.start()
        procs.append((p,time()))

        # Wait in line
        while len(procs) >= maxprocs:
            for p,t in procs:
                if not p.is_alive():
                    procs.remove((p,t))
                elif time()-t > 900:
                    p.terminate()
                    procs.remove((p,t))

    # waiting for last ones
    for p,t in procs:
        p.join()
