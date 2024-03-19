# imports
import boundarytools as bt
import os
import json
import numpy as np
import math
import csv
from urllib.request import urlopen
import pythongis as pg

# params
branch = 'gadm4'
iso = 'TCD'
lvl = 2

# load country boundaries
#url = 'https://www.geoboundaries.org/data/geoBoundariesCGAZ-3_0_0/ADM0/simplifyRatio_10/geoBoundariesCGAZ_ADM0.geojson'
#geoj = boundarytools.utils.load_geojson_url(url)
with open('data/gb-countries-simple.json') as r:
    geoj = json.loads(r.read())
countries = pg.VectorData()
countries.fields = list(geoj['features'][0]['properties'].keys())
print(countries.fields)
for f in geoj['features']:
    countries.add_feature(f['properties'], f['geometry'])

# define map
def makemap(geoj, country, source, level):
    print(country, source, level)
    print('count', len(geoj['features']))

    #crs = '+proj=aea +lat_1=27 +lat_2=45 +lat_0=35 +lon_0=105 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +no_defs'
    
    m = pg.renderer.Map(1200,1000,background='white') #,crs=crs)
    m.add_layer(countries, fillcolor='lightgray', outlinewidth='0.5px')
    
    d = pg.VectorData()
    for f in geoj['features']:
        d.add_feature([], f['geometry'])
    #color = pg.renderer.rgb('blue')[:3] + (200,)
    color = (30, 144, 255, 200)
    m.add_layer(d, fillcolor=color)
    
    m.zoom_bbox(*d.bbox, geographic=True)
    m.save('figures/temporal-{}-ADM{}-{}.png'.format(country, level, source))

# load nigeria sources
##iso,lvl = 'CIV',1
##sources = bt.utils.find_geocontrast_sources(iso, lvl)
##for src,url in sources.items():
##    if src == 'OSM-Boundaries': continue
##    geoj = bt.utils.load_topojson_url(url)
##    makemap(geoj, iso, src, lvl)

# load germany sources
##iso,lvl = 'DEU', 1
##sources = bt.utils.find_geocontrast_sources(iso, lvl)
##for src,url in sources.items():
##    if src == 'OSM-Boundaries': continue
##    geoj = bt.utils.load_topojson_url(url)
##    makemap(geoj, iso, src, lvl)

# load country sources
sources = bt.utils.find_geocontrast_sources(iso, lvl, branch=branch)
for src,url in sources.items():
    #if src == 'OSM-Boundaries': continue
    url = url.replace('/stable/', f'/{branch}/')
    geoj = bt.utils.load_topojson_url(url)
    makemap(geoj, iso, src, lvl)





    
