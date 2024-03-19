
import pythongis as pg

# define sat tile loading
# some parts from https://jimmyutterstrom.com/blog/2019/06/05/map-tiles-to-geotiff/
import io
import os
from math import log, tan, radians, cos, pi, floor, ceil, degrees, atan, sinh
from PIL import Image
from urllib.request import urlopen
import requests
import hashlib

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
    if x < 0: 
        x = 0
    elif x > 1:
        x = 1
    if lat <= -90:
        y = 1
    elif lat >= 90:
        y = 0
    else:
        y = (1 - log(tan(radians(lat)) + sec(radians(lat))) / pi) / 2
    return((tile_count-0)*x, (tile_count-0)*y)

# def latlon_to_xyz(lat, lon, z):
#     # according to chatgpt
#     lat_rad = radians(lat)
#     n = 2.0 ** z
#     xtile = int((lon + 180.0) / 360.0 * n)
#     ytile = int((1.0 - log(tan(lat_rad) + (1 / cos(lat_rad))) / pi) / 2.0 * n)
#     return (xtile, ytile)

def bbox_to_xyz(lon_min, lat_min, lon_max, lat_max, z):
    x_min, y_max = latlon_to_xyz(lat_min, lon_min, z)
    x_max, y_min = latlon_to_xyz(lat_max, lon_max, z)
    return(floor(x_min), floor(y_min),
           floor(x_max), floor(y_max))

def fetch_tile(x, y, z, tile_server):
    url = tile_server.replace(
        "{x}", str(x)).replace(
        "{y}", str(y)).replace(
        "{z}", str(z))
    print('fetching from', url)
    #fobj = io.BytesIO(urlopen(url).read())
    headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.36'}
    fobj = io.BytesIO(requests.get(url, headers=headers).content)
    img = Image.open(fobj)
    return img

def fetch_tilecache(x, y, z, tile_server):
    folder = 'temp'
    if not os.path.exists(folder):
        os.mkdir(folder)
    urlhash = hashlib.md5(tile_server.encode('utf8')).hexdigest()
    tilename = f'tile_cache_{urlhash}_x{x}_y{y}_z{z}.png'
    path = f'{folder}/{tilename}'
    print(path)
    # fetch from cache or url
    if os.path.exists(path):
        # fetch cached tile image
        img = Image.open(path)
    else:
        # fetch img from url
        img = fetch_tile(x, y, z, tile_server)
        # store in cache
        img.save(path)
    # return
    return img

# from https://github.com/jimutt/tiles-to-tiff/blob/master/tiles_to_tiff/tile_convert.py

def mercatorToLat(mercatorY):
    return(degrees(atan(sinh(mercatorY))))


def y_to_lat_edges(y, z):
    tile_count = pow(2, z)
    unit = 1 / tile_count
    relative_y1 = y * unit
    relative_y2 = relative_y1 + unit
    lat1 = mercatorToLat(pi * (1 - 2 * relative_y1))
    lat2 = mercatorToLat(pi * (1 - 2 * relative_y2))
    return(lat1, lat2)


def x_to_lon_edges(x, z):
    tile_count = pow(2, z)
    unit = 360 / tile_count
    lon1 = -180 + x * unit
    lon2 = lon1 + unit
    return(lon1, lon2)


def tile_edges(x, y, z):
    lat2, lat1 = y_to_lat_edges(y, z)
    lon1, lon2 = x_to_lon_edges(x, z)
    return[lon1, lat1, lon2, lat2]

# custom

def as_raster(img, _x1, _y1, _x2, _y2):
    crs = None #'+init=EPSG:3857'
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
    return r

def as_vector(_x1, _y1, _x2, _y2):
    crs = None #'+init=EPSG:3857'
    d = pg.VectorData(crs=crs)
    poly = [[(_x1,_y1),(_x2,_y1),(_x2,_y2),(_x1,_y2)]]
    geom = {'type':'Polygon', 'coordinates':poly}
    d.add_feature([], geom)
    return d

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
    #if mapp.crs != '+init=EPSG:4326':
    #    _transform = pg.renderer.get_crs_transformer(mapp.crs, '+init=EPSG:4326')
    #    bbox = pg.renderer.reproject_bbox(bbox, _transform)

    # determine zoom level from map extent and size
    # how many 256 tiles fit in map
    from math import floor, ceil, log10, log2
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
    # n = number of tiles needed to cover world in each direction at zoom level z
    # n = 2**z  ->  log n = z * log 2  ->  log n / log 2 = z
    n = tilextot
    z = floor( log2(n) )
    #z = floor( log10(n) / log10(2) )
    #z = int(floor(log(n) / log(2)))
    print(z)
    
    # padding
    if padding:
        xmin,ymin,xmax,ymax = bbox
        w,h = abs(xmax-xmin),abs(ymax-ymin)
        xpad,ypad = w*padding, h*padding
        bbox = [xmin-xpad,ymin-ypad,xmax+xpad,ymax+ypad]

    # loop range of xy coords within bbox
    print('bbox',bbox)
    xmin,ymin,xmax,ymax = bbox_to_xyz(*bbox, z) # xy coords
    #xmin,ymin = latlon_to_xyz(bbox[3], bbox[0], z) # note that input is ymax,xmin
    #xmax,ymax = floor(xmin+tilexfit), floor(ymin+tileyfit)
    #xmin,ymin = map(floor, (xmin,ymin))
    print('tile ranges',xmin,ymin,xmax,ymax)
    for x in range(xmin, xmax + 1):
        for y in range(ymin, ymax + 1):
            print(f"{x},{y},{z}")
            # get merc coorrds
            #_x1,_y1 = xyz2merc(x, y+1, z)
            #_x2,_y2 = xyz2merc(x+1, y, z)
            _x1,_y1,_x2,_y2 = tile_edges(x, y, z)
            # fetch tile
            img = fetch_tilecache(x, y, z, tile_server)
            # create raster
            r = as_raster(img, _x1, _y1, _x2, _y2)
            # create vector
            v = as_vector(_x1, _y1, _x2, _y2)
            yield r,v