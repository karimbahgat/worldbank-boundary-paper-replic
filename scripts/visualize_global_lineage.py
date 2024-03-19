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
SOURCES = ['geoBoundaries (Open)', 'OCHA', 'SALB'] #'GADM', 'OpenStreetMap', 'SALB', 'OCHA', 'Natural_Earth']

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
    
def get_source_class(source, sources):
    sources = [s.lower() for s in sources if s]
    # gadm = contains "gadm"
    test = lambda s: 'gadm' in s
    if any([src for src in sources if test(src)]):
        return 'GADM'
    # gaul = contains "gaul"
    test = lambda s: 'gaul' in s
    if any([src for src in sources if test(src)]):
        return 'GAUL'
    # salb = contains "salb"
    test = lambda s: 'salb' in s
    if any([src for src in sources if test(src)]):
        return 'SALB'
    # osm = contains "osm"
    test = lambda s: 'osm' in s or 'openstreetmap' in s.replace(' ','')
    if any([src for src in sources if test(src)]):
        return 'OSM'
    # ocha = contains "ocha"
    # note, ignore listing ocha as a source for the ocha dataset/collection
    test = lambda s: 'ocha' in s
    if source != 'OCHA' and any([src for src in sources if test(src)]):
        return 'OCHA'
    # gov sources = contains "statist*","census","agency","bureau","ministr*","department","government*","national"
    # note, salb dataset is automatically authoritative
    test = lambda s: 'statist' in s \
                        or 'census' in s \
                        or 'agency' in s \
                        or 'bureau' in s \
                        or 'ministr' in s \
                        or 'department' in s \
                        or 'government' in s \
                        or 'national' in s
    if source == 'SALB' or any([src for src in sources if test(src)]):
        return 'Authoritative'
    # other
    return 'Other'


def get_country_source_stats(iso, level, source):
    # open source pair stats
    countrylevelrows = [r
                        for r in META
                        if r['boundaryISO']==iso
                        and r['boundaryType']=='ADM{}'.format(level)
                        and r['boundaryCollection'] == source
                       ]
    if len(countrylevelrows) == 0:
        return {'source_type': None}

    elif len(countrylevelrows) == 1:
        row = countrylevelrows[0]

        # get source classification based on all secondary or higher source fields
        # (source-1 lists the source dataset itself)
        stats = {}
        sources = [val for field,val in row.items() 
                    if field.startswith('boundarySource-') and not field.endswith('-1')]
        print(sources)
        stats['source_type'] = get_source_class(source, sources)

        #print(iso,source,levels,stats)
        return stats

    else:
        raise Exception('More than one row entry')
      
# collect for each source
level = 2
for f in geoj['features']:
    props = f['properties']
    iso = props['shapeISO']
    #print(iso)
    for source in SOURCES:
        #print('ADM',level)
        stats = get_country_source_stats(iso, level, source)
        #print(stats)
        for k,v in stats.items():
            k = k + '_' + str(source).replace(' ', '_')
            props[k] = v
        #print(props)
    print(iso,props)
            
# some figure configs
#import matplotlib.pyplot as plt
#plt.rcParams['axes.grid'] = False
#plt.rcParams['lines.linewidth'] = 0.5

#breaks = [0,1,2,3,4,5]
#classes = 5
#colors = [(146,3,85), (230,230,230), (45,107,26)] # ['red','white','green']

def save_map(geoj, color_by, title, output):
    d = pg.VectorData()
    d.fields = list(geoj['features'][0]['properties'].keys())
    for f in geoj['features']:
        if f['properties']['shapeName'] == 'Antarctica':
            continue
        d.add_feature(f['properties'], f['geometry'])
    print(d)
    m = pg.renderer.Map(4000,2000,background='lightblue')
    m.zoom_bbox(-180,-70,180,90)
    m.add_layer(d, fillcolor='gray', outlinewidth='0.3px',
                transparency=0.3, legend=False)
    _colors = {'Authoritative':'blue',
                'GADM':'purple',
                'GAUL':'red',
                'OSM':'green',
                'SALB':'pink',
                'OCHA':'orange',
                'Other':'brown'}
    m.add_layer(d.select(lambda f: f[color_by] != None),
                fillcolor={'key':color_by,
                              'breaks':'unique',
                              #'classes':classes,
                              'colors':_colors,
                              },
                outlinewidth='2.5px',
                legendoptions={'title':title}) #, 'valueformat':lambda v: format(v,'.0f')}) #, 'direction':'e'})
    m.add_legend({'fillcolor':None,'outlinecolor':None})#, xy=('50%w','99%h'), anchor='s')
    m.save(output)
            
# visualize
for source in SOURCES:
    source = source.replace(' ', '_')
    print(source)
    save_map(geoj, 'source_type_{}'.format(source),
             'Sourced from'.format(source),
             'figures/source_type_{}.png'.format(source),
             )
    


    
