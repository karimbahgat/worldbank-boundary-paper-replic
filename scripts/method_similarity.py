# imports
import boundarytools as bt
import os
import json
import numpy as np
import math
import csv
from urllib.request import urlopen
import pythongis as pg

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

def geoj2data(geoj):
    d = pg.VectorData()
    d.fields = list(geoj['features'][0]['properties'].keys())
    for f in geoj['features']:
        d.add_feature(f['properties'], f['geometry'])
    return d

def similarity(f1, f2):
    print(f1,f2)
    shp1 = f1.get_shapely()
    shp2 = f2.get_shapely()
    isec = shp1.intersection(shp2)
    union = shp1.union(shp2)
    simil = (isec.area/union.area) * 100
    print('simil',simil)
    return isec,union,simil

# define map
def makemap(feat1, feat2, data2, mapname):
    print(mapname)

    #crs = '+proj=aea +lat_1=27 +lat_2=45 +lat_0=35 +lon_0=105 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +no_defs'
    
    m = pg.renderer.Map(1100,1000,background='white') #,crs=crs)

    # countries? 
    #m.add_layer(countries, fillcolor='lightgray', outlinewidth='0.5px')

    # background gray
    m.add_layer(data2, fillcolor='lightgray', outlinecolor=None, 
                legend=False, #legendoptions={'title':'source B comparisons'}
                )

    # calc similarity
    isec,union,simil = similarity(feat1, feat2)
    print('isec area', isec.area)
    print('union area', union.area)
    print('% similarity', simil)

    # feat 1 nonmatch
    color = pg.renderer.rgb('blue')
    m.add_layer(feat1.dataset(), fillcolor=color[:3]+(100,), outlinecolor=None,
                legend=False, #legendoptions={'title':'source A unit'},
                )

    # feat 2 nonmatch
    color = pg.renderer.rgb('red')
    m.add_layer(feat2.dataset(), fillcolor=color[:3]+(100,), outlinecolor=None,
                legend=False, #legendoptions={'title':'source B match = {}'.format(feat2['NAME_1'])},
                )

    # isec
    _d = pg.VectorData()
    _d.add_feature([], isec.__geo_interface__)
    #fisec = pg.vector.data.Feature(_d, [], isec.__geo_interface__)
    #color = (100,100,100)
    color = (30, 144, 255, 200)
    m.add_layer(_d, fillcolor=color, outlinecolor=color, outlinewidth='4px',
                legendoptions={'title':'intersection = {} km2'.format(round(isec.area,2))},
                )

    # source B red outlines
    color = pg.renderer.rgb('red')
    m.add_layer(data2, fillcolor=None, outlinecolor=color, outlinewidth='1px',
                legendoptions={'title':'source B comparisons'}
                )

    # feat 2
    m.add_layer(feat2.dataset(), fillcolor=None, outlinecolor='red', outlinewidth='5px',
                legendoptions={'title':'source B match = {}'.format(feat2['NAME_1'])},
                )

    # feat 1
    m.add_layer(feat1.dataset(), fillcolor=None, outlinecolor='black', outlinewidth='5px',
                legendoptions={'title':'source A unit = {}'.format(feat1['shapeName'])},
                )
    
    m.zoom_auto()
    m.zoom_out(1.2)

    title = '{}% Match'.format(round(simil, 1))
    titleoptions = {'fillcolor':None, 'outlinecolor':None, 'xy':('1%w','1%h'), 'anchor':'nw'}
    m.title = title
    m.titleoptions = titleoptions
    legend_layers = [m.layers[i] for i in [-1,-2,-4,-3]]
    m.add_legend({'layers':legend_layers, 'fillcolor':None, 'outlinecolor':None, 'direction':'s'}, #, 'title':title, 'titleoptions':titleoptions},
                 xy=('1%w','7%h'), anchor='nw')
    m.save('figures/similarity-{}.png'.format(mapname))

# load sources
iso,lvl = 'OMN', 1
sources = bt.utils.find_geocontrast_sources(iso, lvl)

geoj = bt.utils.load_topojson_url(sources['geoBoundaries (Open)'])
d = geoj2data(geoj)

geoj2 = bt.utils.load_topojson_url(sources['GADM v3.6'])
d2 = geoj2data(geoj2)

for i,f1 in enumerate(d):
    i += 1 # 1-based
    name1 = f1['shapeName']
    # get best match
    sortby = lambda feat2: similarity(f1,feat2)[-1]
    f2 = sorted(d2, key=sortby)[-1]
    if similarity(f1, f2)[-1] == 0:
        continue
    name2 = f2['NAME_1']
    print(i, name1, name2)
    # make map
    mapname = '{}-ADM{}-{}'.format(iso, lvl, i)
    makemap(f1, f2, d2, mapname)





    
