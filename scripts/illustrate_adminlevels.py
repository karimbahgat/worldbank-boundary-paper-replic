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
iso = 'MWI'

# identify possible

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

# load main country bounds
for f in countries:
    if f['shapeISO'] == iso:
        break
polys = f.geometry['coordinates']
largest = sorted(polys, key=lambda poly: len(poly[0]))[-1]
xs,ys = zip(*largest[0])
bbox = min(xs),min(ys),max(xs),max(ys)

# define map
def makemap(geoj, source, level):
    print(source, level)
    print('count', len(geoj['features']))

    crs = '+proj=aea +lat_1=27 +lat_2=45 +lat_0=35 +lon_0=105 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +no_defs'
    
    m = pg.renderer.Map(1000,1000,background='white') #,crs=crs)
    m.add_layer(countries, fillcolor='lightgray', outlinewidth='0.5px')
    
    d = pg.VectorData()
    for f in geoj['features']:
        d.add_feature([], f['geometry'])
    color = (30, 144, 255, 200)
    m.add_layer(d, fillcolor=color) #(45,107,26))
    
    m.zoom_bbox(*bbox, geographic=True) #*d.bbox, geographic=True)
    m.zoom_out(1.1)
    m.save('figures/adminlevels-{}-ADM{}.png'.format(source, level))

# load china sources
for lvl in range(1, 3+1):
    print(lvl)
    sources = bt.utils.find_geocontrast_sources(iso, lvl, branch=branch)
    for src,url in sources.items():
        if 'Authoritative' in src: 
            continue
        url = url.replace('/stable/', f'/{branch}/')
        geoj = bt.utils.load_topojson_url(url)
        makemap(geoj, src, lvl)





    
