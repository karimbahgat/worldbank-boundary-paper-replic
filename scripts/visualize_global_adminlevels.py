# imports
import boundarytools
import os
import json
import numpy as np
import math
import csv
from urllib.request import urlopen
import itertools
import pythongis as pg

# params
BRANCH = 'gadm4'
SOURCES = ['geoBoundaries (Open)', 'GADM', 'OpenStreetMap', 'SALB', 'OCHA', 'Natural_Earth']

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

def calc_source_level_matches(level, countryrows):
    # admin matching algorithm
    # for iso,level
    #   for each reference source
    #     for each comparison source
    #       find the comparison level that's closest to the reference
    #       eg gadm1=13, vs salb0=1, salb1=10 (MATCH), salb2=25
    #       if the closest comparison has same level as reference, set dummy to 1
    #   calc percentage of all source combinations considered to be same level (share of dummies)
    key = lambda r: r['boundaryCollection']
    sourcegroups = {src:list(group)
                    for src,group in itertools.groupby(sorted(countryrows, key=key), key=key)
                    }

    match_matrix = {}
    for src1,group1 in sourcegroups.items():
        match_row = {}
        match_matrix[src1] = match_row

        row_at_level = None
        for r in group1:
            if r['boundaryType']=='ADM{}'.format(level):
                row_at_level = r

        if not row_at_level:
            continue

        for src2,group2 in sourcegroups.items():
            if src1 == src2: continue
            levelcounts = [(r['boundaryType'],int(r['boundaryCount']))
                            for r in group2]
            if not f'ADM{level}' in [r['boundaryType'] for r in group2]:
                # it must be possible to match with comparison source 
                # ie, both must have data for the same admin level
                continue
            key = lambda r: abs( int(row_at_level['boundaryCount']) - r[1] )
            nearest_match = sorted(levelcounts, key=key)[0]
            match_at_same_level = f'ADM{level}' == nearest_match[0]
            match_row[src2] = match_at_same_level
    
    return match_matrix

def get_country_level_stats(iso, level):
    # open source pair stats
    countryrows = [r
                    for r in META
                    if r['boundaryISO']==iso
                    ]
    if SOURCES:
        countryrows = [r for r in countryrows
                            if r['boundaryCollection'] in SOURCES]
    stats = {}
    match_matrix = calc_source_level_matches(level, countryrows)
    print(match_matrix)
    vals = [v for row in match_matrix.values() for v in row.values()]
    stats['percent_admincount_same_level'] = np.mean(vals) * 100 if vals else None
    print(iso, level, stats)

    #if iso=='DZA' and level == 2:
    #    x=1 # set breakpoint here for debugging

    #counts = [r['boundaryCount'] for r in countryrows]
    #counts = [float(v) for v in counts if v != 'Unknown']
    #counts = [v for v in counts if not math.isnan(v)]
    #stats = {}
    #stats['admincounts_list'] = counts
    #stats['admincounts_mean'] = np.mean(counts) if counts else None
    #stats['admincounts_std'] = np.std(counts) if counts else None
    #stats['admincounts_stdperc'] = (stats['admincounts_std']/stats['admincounts_mean'])*100 if counts else None

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
breaks = [0,25,50,75,95,100] #'natural'
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
            
# visualize country admin match percentage
for level in range(0, 4+1):
    print(level)
    save_map(geoj, 'percent_admincount_same_level{}'.format(level),
             'ADM{} Same Level\nAcross Sources (%)'.format(level),
             'figures/bcount_match{}.png'.format(level),
             reverse_colors=False)
    


    
