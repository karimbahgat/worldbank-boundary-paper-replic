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
            
# some figure configs
breaks = [0,50,80,90,95,100]
classes = 5
colors = [(146,3,85), (230,230,230), (45,107,26)] # ['red','white','green']

def save_map(geoj, color_by, title, output, reverse_colors=False):
    d = pg.VectorData()
    d.fields = list(geoj['features'][0]['properties'].keys())
    for f in geoj['features']:
        if f['properties']['shapeName'] == 'Antarctica':
            continue
        d.add_feature(f['properties'], f['geometry'])
    print(d)
    m = pg.renderer.Map(4000,2000,background='lightblue')
    m.zoom_bbox(-180,-70,180,90)
    if reverse_colors:
        _colors = list(reversed(colors))
    else:
        _colors = colors
    m.add_layer(d, fillcolor='gray', outlinewidth='0.3px',
                transparency=0.3, legend=False)
    m.add_layer(d.select(lambda f: f[color_by] != None),
                fillcolor={'key':color_by,
                              'breaks':breaks,
                              'classes':classes,
                              'colors':_colors},
                outlinewidth='2.5px',
                legendoptions={'title':title, 'valueformat':lambda v: format(v,'.0f'), 'direction':'n'})
    m.add_legend({'fillcolor':None,'outlinecolor':None})#, xy=('50%w','99%h'), anchor='s')
    m.save(output)
            
# visualize country year stdev
for level in range(0, 4+1):
    print(level)
    #fig = boundarytools.utils.show_dataset(geoj, color_by='yr_std{}'.format(level))
    #fig.savefig('figures/yr_std{}.png'.format(level))
    save_map(geoj, 'prob{}'.format(level),
             'ADM{} Percent\nSource Agreement'.format(level),
             'figures/agree{}.png'.format(level),
             )#reverse_colors=True)



