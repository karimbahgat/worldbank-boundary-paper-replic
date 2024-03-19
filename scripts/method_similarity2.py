# imports
import boundarytools as bt
import os
import json
import numpy as np
import math
import csv
from urllib.request import urlopen
import pythongis as pg

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

# load sources
iso,lvl = 'OMN', 1
sources = bt.utils.find_geocontrast_sources(iso, lvl)

geoj = bt.utils.load_topojson_url(sources['geoBoundaries (Open)'])
d = geoj2data(geoj)
for f in d:
    if f['shapeName'] == 'Al Wusta':
        print('found', f['shapeName'])
        f1 = f

geoj2 = bt.utils.load_topojson_url(sources['GADM v4.0.4'])
d2 = geoj2data(geoj2)
print(d2.fields)
for f in d2:
    if f['NAME_1'] == 'Al Wusta':
        print('found', f['NAME_1'])
        f2 = f

# render source 1 unit
m = pg.renderer.Map(2000,2000,background='white') #,crs=crs)
m.add_layer(f1.dataset(), fillcolor='lightgray')
m.zoom_auto()
m.zoom_out(1.1)
m.save('figures/similarity-simple-source1.png')

# render source 2 unit
m = pg.renderer.Map(2000,2000,background='white') #,crs=crs)
m.add_layer(f2.dataset(), fillcolor='lightgray')
m.zoom_auto()
m.zoom_out(1.1)
m.save('figures/similarity-simple-source2.png')

# render both intersection
isec,union,simil = similarity(f1, f2)
print('isec area', isec.area)
print('union area', union.area)
print('% similarity', simil)
_isec = pg.VectorData()
_isec.add_feature([], isec.__geo_interface__)
m = pg.renderer.Map(2000,2000,background='white') #,crs=crs)
m.add_layer(f1.dataset(), fillcolor='lightgray')
m.add_layer(f2.dataset(), fillcolor='lightgray')
m.add_layer(_isec, fillcolor='gray')
m.zoom_auto()
m.zoom_out(1.1)
m.save('figures/similarity-simple-both.png')




    
