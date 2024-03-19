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
SOURCES = ['geoBoundaries (Open)', 'GADM v3.6', 'OSM-Boundaries', 'SALB', 'geoBoundaries (Humanitarian)', 'Natural Earth v4.1', 'IPUMS']

# load country boundaries
#url = 'https://www.geoboundaries.org/data/geoBoundariesCGAZ-3_0_0/ADM0/simplifyRatio_10/geoBoundariesCGAZ_ADM0.geojson'
#geoj = boundarytools.utils.load_geojson_url(url)
with open('data/gb-countries-simple.json') as r:
    geoj = json.loads(r.read())
    
# collect stats
def load_meta():
    url = 'https://raw.githubusercontent.com/wmgeolab/geoContrast/main/releaseData/geoContrast-meta.csv'
    raw = urlopen(url).read().decode('utf8')
    print(len(raw), raw[:100])
    reader = csv.DictReader(raw.split('\n'))
    return list(reader)
      
META = load_meta()
    
def get_country_level_stats(iso, level):
    # open source pair stats
    countrylevelrows = [r
                        for r in META
                        if r['boundaryISO']==iso
                        and r['boundaryType']=='ADM{}'.format(level)
                       ]
    if SOURCES:
        countrylevelrows = [r for r in countrylevelrows
                            if r['boundarySource-1'] in SOURCES]
    lineres = [r['boundaryCount'] for r in countrylevelrows]
    lineres = [float(v) for v in lineres if v != 'Unknown']
    lineres = [v for v in lineres if not math.isnan(v)]
    stats = {}
    stats['admincounts_list'] = lineres
    stats['admincounts_mean'] = np.mean(lineres) if lineres else None
    stats['admincounts_std'] = np.std(lineres) if lineres else None
    stats['admincounts_stdperc'] = (stats['admincounts_std']/stats['admincounts_mean'])*100 if lineres else None
    return stats
      
# collect for each level
for f in geoj['features']:
    props = f['properties']
    iso = props['shapeISO']
    #print(iso)
    for level in range(0, 4+1):
        #print('ADM',level)
        stats = get_country_level_stats(iso, level)
        for k,v in stats.items():
            k = k + str(level)
            props[k] = v

        # special stats comparing 
            
# some figure configs
import matplotlib.pyplot as plt
plt.rcParams['axes.grid'] = False
plt.rcParams['lines.linewidth'] = 0.5

breaks = 'natural'
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
                legendoptions={'title':title}) #, 'valueformat':lambda v: format(v,'.0f')}) #, 'direction':'e'})
    m.add_legend({'fillcolor':None,'outlinecolor':None})#, xy=('50%w','99%h'), anchor='s')
    m.save(output)
            
# visualize country year stdev
for level in range(0, 4+1):
    print(level)
    save_map(geoj, 'admincounts_mean{}'.format(level),
             'ADM{} Average\n# of Divisions'.format(level),
             'figures/bcount_mean{}.png'.format(level),
             reverse_colors=True)
    save_map(geoj, 'admincounts_std{}'.format(level),
             'ADM{} Variation in\n# of Divisions'.format(level),
             'figures/bcount_std{}.png'.format(level),
             reverse_colors=True)
    save_map(geoj, 'admincounts_stdperc{}'.format(level),
             'ADM{} Variation in\n# of Divisions (%)'.format(level),
             'figures/bcount_stdperc{}.png'.format(level),
             reverse_colors=True)
    


    
