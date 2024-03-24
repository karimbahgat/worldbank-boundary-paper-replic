# imports
import os
import json
import numpy as np
import math
import csv
import itertools
from urllib.request import urlopen

# set all paths relative to this script
curdir = os.path.dirname(__file__)
os.chdir(curdir)

# params
BRANCH = 'worldbank-replication' # which branch of the boundary data repo (geocontrast)
SOURCES = ['geoBoundaries (Open)', 'GADM v4.0.4', 'OSM-Boundaries', 'UN SALB', 'OCHA', 'Natural Earth v5.0.1']
OUTPUT_DIR = os.path.abspath('../output/tables')

# load country boundaries
#url = 'https://www.geoboundaries.org/data/geoBoundariesCGAZ-3_0_0/ADM0/simplifyRatio_10/geoBoundariesCGAZ_ADM0.geojson'
#geoj = boundarytools.utils.load_geojson_url(url)
with open('../data/gb-countries-simple.json') as r:
    geoj = json.loads(r.read())

# load source metadata
def load_meta():
    url = f'https://raw.githubusercontent.com/wmgeolab/geoContrast/{BRANCH}/releaseData/geoContrast-meta.csv'
    raw = urlopen(url).read().decode('utf8')
    print(len(raw), raw[:100])
    reader = csv.DictReader(raw.split('\n'))
    return list(reader)

META = load_meta()

########################

def collect_countrystats():
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
    return countrystats

def define_isoregions():
    regions = 'World,Northern America,Latin America and the Caribbean,Europe,Africa,Asia,Oceania'.split(',')
    isoregions = {}
    key = lambda r: r['boundaryISO']
    for iso,group in itertools.groupby(sorted(META, key=key), key=key):
        reg = list(group)[0]['Continent']
        isoregions[iso] = reg
    return regions,isoregions

def aggregate_levelstats(countrystats, isoregions, regions):
    levelstats = {}
    for level in range(0, 3+1):
        levelstats[level] = indics = {}
        for regi,reg in enumerate(regions):
            # lineres
            indics[f'lineres_{regi}'] = sourcestats = {}
            for src in SOURCES:
                vals = [_levelstats[level]['lineres'].get(src, None) 
                        for _iso,_levelstats in countrystats.items()
                        if regi == 0 or isoregions[_iso]==reg]
                vals = [float(v) for v in vals if v not in (None,'nan')]
                sourcestats[src] = '{:,.0f}'.format(np.mean(vals)) if vals else '-'
            # year
            indics[f'year_{regi}'] = sourcestats = {}
            for src in SOURCES:
                vals = [_levelstats[level]['year'].get(src, None) 
                        for _iso,_levelstats in countrystats.items()
                        if regi == 0 or isoregions[_iso]==reg]
                vals = [float(v) for v in vals if v not in (None,'Unknown')]
                sourcestats[src] = '{:.1f}'.format(np.mean(vals)) if vals else '-'
    
    return levelstats

def latex_table(levelstats, indic):
    rows = []
    for level,indics in levelstats.items():
        rows.append('\\hline')
        rows.append(f'\\textbf{{ADM{level}}} \\\\')
        for src in SOURCES:
            indicvals = []
            indicvals += [indics[indic+f'_{regi}'][src] for regi,reg in enumerate(regions)]
            valstrings = ' & '.join(indicvals)
            row = f'\hspace{{0.1cm}} {src} & {valstrings} \\\\'
            rows.append(row)
    return '\n'.join(rows)


if __name__ == '__main__':

    # collect misc stats for each country and level
    countrystats = collect_countrystats()

    # create region lookup
    regions,isoregions = define_isoregions()

    # aggregate across countries, so we get level-indic-source
    levelstats = aggregate_levelstats(countrystats, isoregions, regions)

    # generate year table
    with open(f'{OUTPUT_DIR}/year_stats_latex_table.txt', 'w') as f:
        f.write( latex_table(levelstats, 'year') )

    # generate lineres table
    with open(f'{OUTPUT_DIR}/lineres_stats_latex_table.txt', 'w') as f:
        f.write( latex_table(levelstats, 'lineres') )
