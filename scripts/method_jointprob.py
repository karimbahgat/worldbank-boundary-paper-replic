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
    #print(f1,f2)
    shp1 = f1.get_shapely()
    shp2 = f2.get_shapely()
    if shp1.intersects(shp2):
        try:
            isec = shp1.intersection(shp2)
            union = shp1.union(shp2)
            simil = (isec.area/union.area) * 100
            print('simil',simil)
            return isec,union,simil
        except:
            pass

    # no overlaps
    return None,None,0






    

# load sources
iso,lvl = 'OMN', 1
sources = bt.utils.find_geocontrast_sources(iso, lvl)

geoj = bt.utils.load_topojson_url(sources['geoBoundaries (Open)'])
d = geoj2data(geoj)

geoj2 = bt.utils.load_topojson_url(sources['GADM v3.6'])
d2 = geoj2data(geoj2)







# stage 1: prob of point in each unit
d.compute('area', lambda f: f.get_shapely().area)
d.compute('areatot', lambda f: f['area'], by=lambda f: 1, stat='sum')
d.compute('areaprob', lambda f: f['area']/f['areatot'])

def makemap(data, mapname):
    print(mapname)

    #crs = '+proj=aea +lat_1=27 +lat_2=45 +lat_0=35 +lon_0=105 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +no_defs'
    
    m = pg.renderer.Map(1200,1000,background='white') #,crs=crs)

    # countries? 
    #m.add_layer(countries, fillcolor='lightgray', outlinewidth='0.5px')

    # data
    color = pg.renderer.rgb('blue')
    m.add_layer(data, fillcolor={'key':'areaprob', 'breaks':'proportional', 'colors':[color[:3]+(100,),color]},
                outlinecolor='black', outlinewidth='3px',
                text=lambda f: '{}%'.format(round(f['areaprob']*100,1)),
                textoptions={'textsize':10},
                legend=False, #legendoptions={'title':'source B comparisons'}
                )
    
    m.zoom_auto()
    m.zoom_out(1.1)

    #title = 'Match Similarity = {}%'.format(round(simil, 1))
    #titleoptions = {'fillcolor':None, 'outlinecolor':None, 'xy':('1%w','1%h'), 'anchor':'nw'}
    #m.title = title
    #m.titleoptions = titleoptions
    #m.add_legend({'fillcolor':None, 'outlinecolor':None, 'direction':'s'}, #, 'title':title, 'titleoptions':titleoptions},
    #             xy=('1%w','7%h'), anchor='nw')
    m.save('figures/jointprob-{}.png'.format(mapname))

mapname = '{}-ADM{}-stage1'.format(iso,lvl)
makemap(d, mapname)






# stage 2: prob of point in each unit of A matching with same unit in B

def findmatch(feat1, data2):
    candidates = []

    # calc all simils
    for feat2 in data2:
        isec,union,simil = similarity(feat1, feat2)
        if simil:
            candidates.append((feat2,isec,union,simil))

    # get closest match
    key = lambda e: e[-1]
    matchfeat,isec,union,simil = sorted(candidates, key=key, reverse=True)[0]

    # calc prob same
    probsame = isec.area / feat1.get_shapely().area

    return matchfeat,isec,probsame

def makemap(data, matches, shared, mapname):
    print(mapname)

    #crs = '+proj=aea +lat_1=27 +lat_2=45 +lat_0=35 +lon_0=105 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +no_defs'
    
    m = pg.renderer.Map(1200,1000,background='white') #,crs=crs)

    # countries? 
    #m.add_layer(countries, fillcolor='lightgray', outlinewidth='0.5px')

    # data background
    color = pg.renderer.rgb('blue')
    m.add_layer(data, fillcolor=color[:3]+(100,),
                outlinecolor=None,
                text=lambda f: '{}%'.format(round(f['probsame']*100,1)),
                textoptions={'textsize':10},
                legend=False, #legendoptions={'title':'Different in source B', 'titleoptions':{'textsize':10}},
                )

    # shared
    sharedarea = sum([pg.vector.geography.Geography(f.geometry).area for f in shared])
    title = 'AB overlap = {:,.0f} km2'.format(sharedarea)
    color = pg.renderer.rgb('blue')
    m.add_layer(shared, fillcolor=color[:3]+(200,), outlinecolor=None,
                legendoptions={'title':title} #, 'titleoptions':{'textsize':10}},
                )

    # data outlines
    totarea = sum([pg.vector.geography.Geography(f.geometry).area for f in data])
    title = 'Source A = {:,.0f} km2'.format(totarea)
    color = pg.renderer.rgb('black')
    m.add_layer(data, fillcolor=None,
                outlinecolor=color, outlinewidth='3px',
                text=lambda f: '{}%'.format(round(f['probsame']*100,1)),
                textoptions={'textsize':10},
                legendoptions={'title':title} #, 'titleoptions':{'textsize':10}}
                )

    # matches
    #color = pg.renderer.rgb('red')
    #m.add_layer(matches, fillcolor=None, outlinecolor=color, outlinewidth='3px')
    
    m.zoom_auto()
    m.zoom_out(1.1)

    simil = sharedarea / totarea * 100
    title = 'P(A=B) = {}%'.format(round(simil, 1))
    titleoptions = {'fillcolor':None, 'outlinecolor':None, 'xy':('3%w','3%h'), 'anchor':'nw'}
    m.title = title
    m.titleoptions = titleoptions
    legend_layers = [m.layers[-1],m.layers[-2]]
    m.add_legend({'layers':legend_layers, 'fillcolor':None, 'outlinecolor':None, 'direction':'s'}, #, 'title':title, 'titleoptions':titleoptions},
                 xy=('5%w','10%h'), anchor='nw')
    m.save('figures/jointprob-{}.png'.format(mapname))

# create layer of nearest matches
d.compute('probsame', lambda f: None)
matches = pg.VectorData(fields=d2.fields)
shared = pg.VectorData(fields=d2.fields)
for f in d:
    matchfeat,isec,probsame = findmatch(f, d2)
    f['probsame'] = probsame
    matches.add_feature(matchfeat.row, matchfeat.geometry)
    shared.add_feature(matchfeat.row, isec.__geo_interface__)

mapname = '{}-ADM{}-stage2'.format(iso,lvl)
makemap(d, matches, shared, mapname)







# stage 3: final probability

def makemap(data, mapname):
    print(mapname)

    #crs = '+proj=aea +lat_1=27 +lat_2=45 +lat_0=35 +lon_0=105 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +no_defs'
    
    m = pg.renderer.Map(1200,1000,background='white') #,crs=crs)

    # countries? 
    #m.add_layer(countries, fillcolor='lightgray', outlinewidth='0.5px')

    # data
    color = pg.renderer.rgb('blue')
    m.add_layer(data, fillcolor={'key':'finalprob', 'breaks':'proportional', 'colors':[color[:3]+(100,),color]},
                outlinecolor='black', outlinewidth='3px',
                text=lambda f: '{}%'.format(round(f['finalprob']*100,1)),
                textoptions={'textsize':10},
                legend=False, #legendoptions={'title':'source B comparisons'}
                )
    
    m.zoom_auto()
    m.zoom_out(1.1)

    #title = 'Match Similarity = {}%'.format(round(simil, 1))
    #titleoptions = {'fillcolor':None, 'outlinecolor':None, 'xy':('1%w','1%h'), 'anchor':'nw'}
    #m.title = title
    #m.titleoptions = titleoptions
    #m.add_legend({'fillcolor':None, 'outlinecolor':None, 'direction':'s'}, #, 'title':title, 'titleoptions':titleoptions},
    #             xy=('1%w','7%h'), anchor='nw')
    m.save('figures/jointprob-{}.png'.format(mapname))

d.compute('finalprob', lambda f: f['areaprob'] * f['probsame'])
mapname = '{}-ADM{}-stage3'.format(iso,lvl)
makemap(d, mapname)




    
