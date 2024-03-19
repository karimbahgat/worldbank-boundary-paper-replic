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
        f = d.add_feature(f['properties'], f['geometry'])
        f.geometry = f.get_shapely().simplify(0.01).buffer(0).__geo_interface__
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
            #print('simil',simil)
            return isec,union,simil
        except:
            pass

    # no overlaps
    return None,None,0

def makemap(data, matches, shared, mapname):
    print(mapname)

    #crs = '+proj=aea +lat_1=27 +lat_2=45 +lat_0=35 +lon_0=105 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +no_defs'
    
    m = pg.renderer.Map(1100,1000,background='white') #,crs=crs)

    # countries? 
    #m.add_layer(countries, fillcolor='lightgray', outlinewidth='0.5px')

    # data background
    color = pg.renderer.rgb('blue')
    m.add_layer(data, fillcolor=color[:3]+(100,),
                outlinecolor=None,
                #text=lambda f: '{}%'.format(round(f['probsame']*100,1)),
                #textoptions={'textsize':10},
                #legend=False, #legendoptions={'title':'Different in source B', 'titleoptions':{'textsize':10}},
                )

    # shared
    sharedarea = sum([pg.vector.geography.Geography(f.geometry).area for f in shared])
    title = 'AB overlap = {:,.0f} km2'.format(sharedarea)
    color = pg.renderer.rgb('blue')
    m.add_layer(shared, fillcolor=color[:3]+(200,), outlinecolor=None,
                #legendoptions={'title':title} #, 'titleoptions':{'textsize':10}},
                )

    # data outlines
    totarea = sum([pg.vector.geography.Geography(f.geometry).area for f in data])
    title = 'Source A = {:,.0f} km2'.format(totarea)
    color = pg.renderer.rgb('black')
    m.add_layer(data, fillcolor=None,
                outlinecolor=color, outlinewidth='3px',
                #text=lambda f: '{}%'.format(round(f['probsame']*100,1)),
                #textoptions={'textsize':10},
                #legendoptions={'title':title} #, 'titleoptions':{'textsize':10}}
                )

    # matches
    #color = pg.renderer.rgb('red')
    #m.add_layer(matches, fillcolor=None, outlinecolor=color, outlinewidth='3px')
    
    m.zoom_auto()
    m.zoom_out(1.1)

    simil = sharedarea / totarea * 100
    print('source agreement', round(simil,1))
    #title = 'P(A=B) = {}%'.format(round(simil, 1))
    #titleoptions = {'fillcolor':None, 'outlinecolor':None, 'xy':('3%w','3%h'), 'anchor':'nw'}
    #m.title = title
    #m.titleoptions = titleoptions
    #legend_layers = [m.layers[-1],m.layers[-2]]
    #m.add_legend({'layers':legend_layers, 'fillcolor':None, 'outlinecolor':None, 'direction':'s'}, #, 'title':title, 'titleoptions':titleoptions},
    #             xy=('5%w','10%h'), anchor='nw')
    m.save('figures/sourcematrix-{}.png'.format(mapname))






    

# load sources
iso,lvl = 'OMN', 1
sources = bt.utils.find_geocontrast_sources(iso, lvl)
keepsources = ['Natural Earth v4.1', 'GADM v3.6', 'OSM-Boundaries', 'geoBoundaries (Open)']
sources = [(src,sources[src]) for src in keepsources]

for i,(src,url) in enumerate(sources):
    i += 1 # 1 based index
    geoj = bt.utils.load_topojson_url(url)
    d = geoj2data(geoj)

    for i2,(src2,url2) in enumerate(sources):
        i2 += 1 # 1 based index
        if src == src2:
            continue

        print(i,i2)
        print(src,src2)
        geoj2 = bt.utils.load_topojson_url(url2)
        d2 = geoj2data(geoj2)




        # match

        def findmatch(feat1, data2):
            #print('finding match')
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

        # create layer of nearest matches
        d.compute('probsame', lambda f: None)
        matches = pg.VectorData(fields=d2.fields)
        shared = pg.VectorData(fields=d2.fields)
        for f in d:
            matchfeat,isec,probsame = findmatch(f, d2)
            f['probsame'] = probsame
            _=matches.add_feature(matchfeat.row, matchfeat.geometry)
            _=shared.add_feature(matchfeat.row, isec.__geo_interface__)




        # make map
        mapname = '{}-ADM{}-{}-{}'.format(iso,lvl,i,i2)
        makemap(d, matches, shared, mapname)

# get row/col totals
import numpy as np
mat = np.array([[100, 94.5, 77.0, 81.1],
                [95.2, 100, 78.9, 81.2],
                [77.7, 79.1, 100, 95.0],
                [73.4, 76.3, 86.9, 100]])
rowtots = [mat[i,:][np.arange(mat.shape[0])!=i].mean() for i in range(mat.shape[0])]
print('rowtots',rowtots)
coltots = [mat[:,i][np.arange(mat.shape[1])!=i].mean() for i in range(mat.shape[1])]
print('coltots',coltots)
print('final',np.array(rowtots).mean())
