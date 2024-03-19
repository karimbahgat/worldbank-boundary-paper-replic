# imports
import boundarytools
import os
import json
import numpy as np
import math
import csv
import itertools
from urllib.request import urlopen

# params (note that gbHum = UN OCHA)
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

# collect misc stats for each country and level
countrystats = {}
#isos = [f['properties']['shapeISO'] for f in geoj['features']]
key = lambda r: r['boundaryISO']
for iso,group in itertools.groupby(sorted(META, key=key), key=key):
    #print(iso)
    group = list(group)
    group = [r for r in group if r['boundaryCollection'] in SOURCES]
    countrystats[iso] = levelstats = {}
    for level in range(0, 3+1):
        levelstats[level] = indics = {}
        # lineres
        sourcestats = dict([(r['boundaryCollection'],r['statsLineResolution']) for r in group 
                            if r['boundaryType']=='ADM{}'.format(level)])
        indics['lineres'] = sourcestats
        # year
        sourcestats = dict([(r['boundaryCollection'],r['boundaryYearRepresented']) for r in group 
                            if r['boundaryType']=='ADM{}'.format(level)])
        indics['year'] = sourcestats
    #print(json.dumps(levelstats, indent=4))

# create region lookup
regions = 'World,Northern America,Latin America and the Caribbean,Europe,Africa,Asia,Oceania'.split(',')
isoregions = {}
for iso,group in itertools.groupby(sorted(META, key=key), key=key):
    reg = list(group)[0]['Continent']
    isoregions[iso] = reg

# aggregate across countries, so we get level-indic-source
levelstats = {}
for level in range(0, 3+1):
    levelstats[level] = indics = {}
    for regi,reg in enumerate(regions):
        # lineres
        indics[f'lineres_{regi}'] = sourcestats = {}
        for src in SOURCES:
            vals = [_levelstats[level]['lineres'].get(src, None) for _iso,_levelstats in countrystats.items()
                    if regi == 0 or isoregions[_iso]==reg]
            vals = [float(v) for v in vals if v not in (None,'nan')]
            sourcestats[src] = '{:,.0f}'.format(np.mean(vals)) if vals else '-'
        # year
        indics[f'year_{regi}'] = sourcestats = {}
        for src in SOURCES:
            vals = [_levelstats[level]['year'].get(src, None) for _iso,_levelstats in countrystats.items()
                    if regi == 0 or isoregions[_iso]==reg]
            vals = [float(v) for v in vals if v not in (None,'Unknown')]
            sourcestats[src] = '{:.1f}'.format(np.mean(vals)) if vals else '-'

# print latex table
def print_table(indic):
    rows = []
    for level,indics in levelstats.items():
        rows.append('\\hline')
        rows.append(f'\\textbf{{ADM{level}}} \\\\')
        for src in SOURCES:
            indicvals = []
            indicvals += [indics[indic+f'_{regi}'][src] for regi,reg in enumerate(regions)]
            valstrings = ' & '.join(indicvals)
            src = src.replace('_',' ').replace(' (Open)','').replace('SALB','UN SALB').replace('geoBoundaries (Humanitarian)','UN OCHA').replace('OpenStreetMap','OSM')
            row = f'\hspace{{0.1cm}} {src} & {valstrings} \\\\'
            rows.append(row)
    print('\n'.join(rows))

print('YEAR')
print_table('year')

print('\nLINERES')
print_table('lineres')
