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

'''
def get_landsat_raster(bbox):
    from landsatxplore.api import API
    user = 'kbahgat'
    pw = input('Please input Earth Explorer password:\n>>> ')
    landsat_api = API(user, pw)

    print(bbox)
    dataset = 'landsat_8_c1'
    scenes = landsat_api.search(
        dataset=dataset,
        bbox=bbox,
        max_cloud_cover=10,
    )
    scenes = list(scenes)

    from landsatxplore.earthexplorer import EarthExplorer, EE_DOWNLOAD_URL, DATA_PRODUCTS
    ee = EarthExplorer(user, pw)
    for scene in scenes[:1]:
        print(scene)
        #path = ee.download(scene['entity_id'], 'temp')
        url = EE_DOWNLOAD_URL.format(
            data_product_id=DATA_PRODUCTS[dataset], entity_id=scene['entity_id']
        )
        path = ee._download(url, 'temp', timeout=300, chunk_size=1024*10)
        r = pg.RasterData(path)
        yield r
'''

# preload global satims
if 0:
    import urllib
    print('dl west')
    urllib.request.urlretrieve('https://eoimages.gsfc.nasa.gov/images/imagerecords/57000/57752/land_shallow_topo_west.tif',
                       'temp/land_shallow_topo_west.tif'
                      )
    print('dl east')
    urllib.request.urlretrieve('https://eoimages.gsfc.nasa.gov/images/imagerecords/57000/57752/land_shallow_topo_east.tif',
                       'temp/land_shallow_topo_east.tif'
                      )
    
# load global sat ims
print('open global sats')
import PIL.Image
PIL.Image.MAX_IMAGE_PIXELS = 466560000*2

scale = 180.0/21600
#affine = [scale,0,-180,
#         0,-scale,90]
#sat_west = pg.RasterData('temp/land_shallow_topo_west.tif', affine=affine)
#sat_west.set_geotransform(affine=affine)
#affine = [scale,0,0,
#         0,-scale,90]
#sat_east = pg.RasterData('temp/land_shallow_topo_east.tif', affine=affine)
#sat_east.set_geotransform(affine=affine)

def get_sat(bbox, padding=0):
    if padding:
        xmin,ymax,xmax,ymin = bbox
        w,h = abs(xmax-xmin),abs(ymax-ymin)
        xpad,ypad = w*padding, h*padding
        bbox = [xmin-xpad,ymax+ypad,xmax+xpad,ymin-ypad]
    if bbox[0] < 0:
        _bbox = list(bbox)
        _bbox[2] = min(0, _bbox[2])
        print(_bbox)
        yield sat_west.manage.crop(_bbox)
    if bbox[2] >= 0:
        _bbox = list(bbox)
        _bbox[0] = max(0, _bbox[0])
        print(_bbox)
        yield sat_east.manage.crop(_bbox)

# define sat tile loading
# some parts from https://jimmyutterstrom.com/blog/2019/06/05/map-tiles-to-geotiff/
import io
from math import log, tan, radians, cos, pi, floor, ceil
from PIL import Image
from urllib.request import urlopen

def xyz2merc(x, y, zoom):
    n = (2 ** zoom) or 1 # only 1 tile at zoom level 0
    # get global mercator bounds
    # https://www.maptiler.com/google-maps-coordinates-tile-bounds-projection/#4/87.28/18.02
    #xmin,ymin,xmax,ymax = [-20037508.342789244, -20037508.342789244, 20037508.342789244, 20037508.342789244]
    xmin,ymin,xmax,ymax = [-20026376.39, -20048966.10, 20026376.39, 20048966.10]
    width,height = xmax-xmin, ymax-ymin
    # convert xy to 0-1 range
    x /= n
    y /= n
    # convert to mercator
    x = xmin + width * x
    y = ymin + height * (1-y)
    return x,y

def sec(x):
    return(1/cos(x))

def latlon_to_xyz(lat, lon, z):
    tile_count = pow(2, z)
    x = (lon + 180) / 360
    y = (1 - log(tan(radians(lat)) + sec(radians(lat))) / pi) / 2
    return(tile_count*x, tile_count*y)

def bbox_to_xyz(lon_min, lon_max, lat_min, lat_max, z):
    x_min, y_max = latlon_to_xyz(lat_min, lon_min, z)
    x_max, y_min = latlon_to_xyz(lat_max, lon_max, z)
    return(floor(x_min), floor(x_max),
           floor(y_min), floor(y_max))

def fetch_tile(x, y, z, tile_server):
    url = tile_server.replace(
        "{x}", str(x)).replace(
        "{y}", str(y)).replace(
        "{z}", str(z))
    print(url)
    fobj = io.BytesIO(urlopen(url).read())
    img = Image.open(fobj)
    return img

def get_sat_tiles(mapp, padding=0):
    '''Retrieve and merge xyz tiles into a single RasterData
    Only works for wgs84 maps so far.
    '''
    bbox = mapp.bbox
    xs,ys = [bbox[0],bbox[2]], [bbox[1],bbox[3]]
    bbox = min(xs),min(ys),max(xs),max(ys)
    #tile_server = 'https://api.maptiler.com/tiles/satellite/{z}/{x}/{y}.jpg?key=' + 'aknzJQRnZg32XVVPrcYH'
    #tile_server = 'https://services.arcgisonline.com/arcgis/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}'
    tile_server = 'http://mt0.google.com/vt/lyrs=s&hl=en&x={x}&y={y}&z={z}'

    # convert to latlon if needed
    # ... 

    # determine zoom level from map extent and size
    # how many 256 tiles fit in map
    from math import floor, ceil, log10
    w,h = mapp.width,mapp.height
    tilexfit,tileyfit = w/256, h/256
    print(tilexfit,tileyfit)
    # how many map extents fit in world
    dx,dy = bbox[2]-bbox[0], bbox[3]-bbox[1]
    extxfit,extyfit = 360/dx, 180/dy
    print(extxfit,extyfit)
    # how many tiles total needed to cover entire world
    tilextot,tileytot = tilexfit*extxfit, tileyfit*extyfit
    print(tilextot,tileytot)
    # n = 2**z  ->  log n = z * log 2  ->  log n / log 2 = z
    n = tilextot
    z = floor( log10(n) / log10(2) )
    print(z)
    
    # padding
    if padding:
        xmin,ymin,xmax,ymax = bbox
        w,h = abs(xmax-xmin),abs(ymax-ymin)
        xpad,ypad = w*padding, h*padding
        bbox = [xmin-xpad,ymin-ypad,xmax+xpad,ymax+ypad]

    # loop range of xy coords within bbox
    print('bbox',bbox)
    #xmax,ymax,xmin,ymin = bbox_to_xyz(*bbox, z) # xy coords
    xmin,ymin = latlon_to_xyz(bbox[3], bbox[0], z) # note that input is ymax,xmin
    xmin,ymin = map(floor, (xmin,ymin))
    xmax,ymax = xmin+floor(tilexfit), ymin+floor(tileyfit)
    print('tile ranges',xmin,ymin,xmax,ymax)
    for x in range(xmin, xmax + 1):
        for y in range(ymin, ymax + 1):
            print(f"{x},{y},{z}")
            img = fetch_tile(x, y, z, tile_server)
            crs = '+init=EPSG:3857'
            _x1,_y1 = xyz2merc(x, y+1, z)
            _x2,_y2 = xyz2merc(x+1, y, z)
            #georef = {'bbox':(_x1,_y1,_x2,_y2)}
            georef = {'xscale':(_x2-_x1)/(img.size[0]-1), # -1 only to cover gaps
                      'yscale':-(_y2-_y1)/(img.size[1]-1), # -1 only to cover gaps
                      'xoffset':_x1,
                      'yoffset':_y2,
                      } 
            print(georef)
            r = pg.RasterData(image=img, crs=crs)
            r.set_geotransform(**georef)
            print(r)
            yield r
    

# define map
def makemap(geoj, geoj2, zoomfunc, zoomfunc2, name):
    print(name)

    #crs = '+proj=aea +lat_1=27 +lat_2=45 +lat_0=35 +lon_0=105 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +no_defs'
    
    m = pg.renderer.Map(1200,1000,background='lightgray') #,crs=crs)
    textoptions = {'textsize': 15}

    # countries? 
    #m.add_layer(countries, fillcolor='lightgray', outlinewidth='0.5px')

    # source 1
    d = pg.VectorData()
    d.fields = list(geoj['features'][0]['properties'].keys())
    for f in geoj['features']:
        d.add_feature(f['properties'], f['geometry'])
    d = d.select(zoomfunc)
    m.add_layer(d, fillcolor=None, outlinecolor='blue', outlinewidth='10px', legend=False)
    m.add_layer(d.convert.to_points('vertex'), fillcolor='blue', fillsize='14px', outlinecolor=None, legendoptions={'title':'geoBoundaries'})

    # source2
    d = pg.VectorData()
    d.fields = list(geoj2['features'][0]['properties'].keys())
    for f in geoj2['features']:
        d.add_feature(f['properties'], f['geometry'])
    d = d.select(zoomfunc2)
    m.add_layer(d, fillcolor=None, outlinecolor='red', outlinewidth='10px', legend=False)
    m.add_layer(d.convert.to_points('vertex'), fillcolor='red', fillsize='14px', outlinecolor=None, legendoptions={'title':'Natural Earth'})
    
    # zoom
    m.zoom_auto()
    #m.zoom_out(1.1)
    #m.zoom_in(3)

    # old
    #m.offset('-35%w', '30%h')
    #m.zoom_in(3)
    #m.offset('-33%w', 0)
    #m.zoom_in(3.5)

    # new
    #m.offset('8%w', '-2%h')
    m.offset('-8%w', '-4%h')
    m.zoom_in(3)
    m.zoom_in(2)

    # add satellite base layer
    bbox = m.bbox
    print(bbox)
    for sat in get_sat_tiles(m):
    #for sat in get_sat(bbox, padding=0.1):
    #for sat in [pg.RasterData('temp/localsat.tif')]: 
        #for b in sat.bands:
        #    b.compute('min(val+10, 255)')
        m.add_layer(sat)
        m.move_layer(-1, 0)

    # save zoomed in map
    m.add_legend({'direction':'s', 'outlinecolor':None}) #, xy=('2%w','93%h'), anchor='sw')
    opts = dict(length=0.09,
                labeloptions={'textcolor':(255,255,255)},
                symboloptions={'fillcolor':(255,255,255)})
    m.add_scalebar(opts) #, xy=('4%w','98%h'), anchor='sw') 
    m.save('figures/resolution-{}-zoomed.png'.format(name))

    # create rectangle layer for current map view
    x1,y1,x2,y2 = bbox
    ring = [(x1,y1),(x2,y1),(x2,y2),(x1,y2),(x1,y1)]
    poly = {'type':'Polygon', 'coordinates':[ring]}
    d = pg.VectorData()
    d.add_feature([], poly)
    m.add_layer(d, fillcolor=None, outlinecolor="black", outlinewidth='8px', legend=False)

    # zoom out
    m.zoom_auto()
    m.zoom_out(1.1)

    # add satellite base layer for new zoom
    bbox = m.bbox
    print(bbox)
    import tiler
    for sat,footprint in tiler.get_sat_tiles(m):
    #for sat in get_sat(bbox, padding=0.1):
    #for sat in [pg.RasterData('temp/localsat.tif')]: 
        #for b in sat.bands:
        #    b.compute('min(val+10, 255)')
        m.add_layer(sat)
        m.move_layer(-1, 0)

    # save zoomed out map with zoom rectangle
    m.render()
    m.save('figures/resolution-{}.png'.format(name))

# load data sources
iso,lvl = 'ETH', 1
sources = bt.utils.find_geocontrast_sources(iso, lvl)
geoj = bt.utils.load_topojson_url(sources['geoBoundaries (Open)'])
geoj2 = bt.utils.load_topojson_url(sources['Natural Earth v5.0.1'])
zoomfunc = lambda f: f['shapeName']=='Tigray'
zoomfunc2 = lambda f: f['name']=='Tigray'
name = '{}-{}-Tigray-gB vs NE'.format(iso, lvl)
makemap(geoj, geoj2, zoomfunc, zoomfunc2, name)





    
