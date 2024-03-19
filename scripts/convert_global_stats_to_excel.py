# imports
import boundarytools
import os
import json
import numpy as np
import math
import csv
from urllib.request import urlopen
import pythongis as pg

# params
BRANCH = 'gadm4'

# load country boundaries
#url = 'https://www.geoboundaries.org/data/geoBoundariesCGAZ-3_0_0/ADM0/simplifyRatio_10/geoBoundariesCGAZ_ADM0.geojson'
#geoj = boundarytools.utils.load_geojson_url(url)
with open('data/gb-countries-simple.json') as r:
    geoj = json.loads(r.read())
    
# collect stats
def load_meta():
    url = f'https://raw.githubusercontent.com/wmgeolab/geoContrast/{BRANCH}/releaseData/geoContrast-meta.csv'
    raw = urlopen(url).read().decode('utf8')
    print(len(raw), raw[:100])
    reader = csv.DictReader(raw.split('\n'))
    return list(reader)
      
META = load_meta()
    
def get_country_level_stats(iso, level):
    # open source pair stats
    try:
        path = 'global_stats/{}-ADM{}-stats.json'.format(iso, level)
        with open(path) as r:
            stats = json.loads(r.read())
        return stats
    except:
        return None

def calc_prob(stats):
    # get all sources
    sources = set()
    for src1,stats2 in stats.items():
        sources.add(src1)
        for src2 in stats2.keys():
            sources.add(src2)
    # create source matrix
    rowtots = []
    for src1 in sources:
        row = []
        for src2 in sources:
            if src1==src2: continue
            prob = stats[src1][src2]
            row.append(prob)
        row = [v for v in row if v != 0 and not math.isnan(v)] # exclude 0, ie weird error
        rowtot = np.mean(row) if row else 0
        rowtots.append(rowtot)
    rowtots = [v for v in rowtots if v != 0] # exclude 0, ie weird error
    tottot = float(np.mean(rowtots)) if rowtots else None
    if tottot and tottot > 100:
        tottot = 100.0
    return tottot
      
# collect for each level
for f in geoj['features']:
    props = f['properties']
    iso = props['shapeISO']
    #print(iso)
    for level in range(0, 4+1):
        #print('ADM',level)
        stats = get_country_level_stats(iso, level)
        if stats:
            #print(stats)
            val = calc_prob(stats)
            #print(val)
            props['prob{}'.format(level)] = None if val is None else val
        else:
            props['prob{}'.format(level)] = None

# output to excel
path = 'global_stats.csv'
with open(path, 'w', newline='') as csvfile:
    fieldnames = list(geoj['features'][0]['properties'].keys())
    fieldnames = [f for f in fieldnames if not f.startswith('prob')]
    fieldnames = sorted(fieldnames)
    fieldnames += [f'prob{lvl}' for lvl in range(4+1)]
    print(fieldnames)
    writer = csv.DictWriter(csvfile, fieldnames)
    writer.writeheader()
    sortby = lambda f: f['properties']['shapeISO']
    for f in sorted(geoj['features'], key=sortby):
        props = f['properties']
        writer.writerow(props)
