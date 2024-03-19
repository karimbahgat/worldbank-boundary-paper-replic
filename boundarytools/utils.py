from zipfile import ZipFile
import io

import numpy as np
# NOTE: must use shape rather than proxy asShape which will be deprecated
# and lead to repetative warnings
# https://github.com/shapely/shapely/issues/1207
from shapely.geometry import shape

import sys
sys.path.append(r'C:\Users\kimok\Documents\GitHub\geoContrast')
import topojson_simple

def get_shapely(feat):
    if 'shapely' in feat:
        return feat['shapely']
    else:
        return shape(feat['geometry'])

def get_bbox(feat):
    #if 'shapely' in feat:
    #    return feat['shapely'].bounds
    if 'bbox' in feat:
        return feat['bbox']
    coords = (p 
            for ring in iter_rings(feat['geometry'])
            for p in ring)
    xs,ys = zip(*coords)
    return min(xs),min(ys),max(xs),max(ys)

def bbox_union(*bboxes):
    xmins,ymins,xmaxs,ymaxs = zip(*bboxes)
    xmin,ymin,xmax,ymax = min(xmins),min(ymins),max(xmaxs),max(ymaxs)
    return xmin,ymin,xmax,ymax

def iter_rings(geoj):
    '''Iterates through all rings of a polygon/multipolygon
    '''
    if geoj['type'] == 'Polygon':
        for ring in geoj['coordinates']:
            yield ring
    elif geoj['type'] == 'MultiPolygon':
        for poly in geoj['coordinates']:
            for ring in poly:
                yield ring

def topo2geoj(data):
    '''Data is the toplevel topojson dict containing: type, objects, arcs'''
    #from topojson.utils import geometry
    #from topojson.ops import dequantize
    # lyr = list(data['objects'].keys())[0]
    # tfeatures = data['objects'][lyr]['geometries']
    # data['arcs'] = [np.array(arc) for arc in data['arcs']] # topojson.utils.geometry() assumes np arrays
    # transform = data.get('transform')
    # if transform:
    #     raise Exception('Quantized coordinates not currently supported')
    #     # transform arg doesnt do anything for polygons, so need to do this manually
    #     # NOT WORKING YET
    #     print(str(data['arcs'])[:1000])
    #     data['arcs'] = [dequantize(arc.T, **transform).T for arc in data['arcs']]
    #     print(str(data['arcs'])[:1000])
    # geoj = {'type': "FeatureCollection", 'features': []}
    # for tfeat in tfeatures:
    #     #print(tfeat['type'], tfeat['properties']) #, tfeat['arcs']) 
    #     feat = {'type': "Feature"}
    #     feat['properties'] = tfeat['properties'].copy()
    #     feat['geometry'] = geometry(tfeat, data['arcs'], transform)
    #     coords = feat['geometry']['coordinates']
    #     geoj['features'].append(feat)
    
    #from topojson.utils import serialize_as_geojson
    #lyr = list(data['objects'].keys())[0]
    #geoj = serialize_as_geojson(data, objectname=lyr)
    geoj = topojson_simple.decode.geojson(data)

    # ignore null-geom features
    geoj['features'] = [feat for feat in geoj['features']
                        if feat['geometry'] != None]

    return geoj

def iter_geocontrast_metatable(branch='stable'):
    from urllib.request import urlopen
    import csv
    url = f'https://raw.githubusercontent.com/wmgeolab/geoContrast/{branch}/releaseData/geoContrast-meta.csv'
    raw = urlopen(url).read().decode('utf8').split('\n')
    reader = csv.DictReader(raw, delimiter=',')
    for row in reader:
        yield row

def find_geocontrast_sources(iso, level, branch='stable'):
    sources = {}
    for row in iter_geocontrast_metatable(branch=branch):
        if row['boundaryISO'] == iso and row['boundaryType'] == 'ADM{}'.format(level):
            source = row['boundarySource-1']
            apiURL = row['apiURL']
            sources[source] = apiURL
    
    return sources

def load_topojson_url(url, load_shapely=False):
    from urllib.request import urlopen
    import json
    if url.endswith('.zip'):
        fobj = io.BytesIO(urlopen(url).read())
        with ZipFile(fobj) as archive:
            filename = url.split('/')[-1].replace('.zip','')
            topoj = json.loads(archive.open(filename).read())
    else:
        topoj = json.loads(urlopen(url).read())
    coll = topo2geoj(topoj)
    if load_shapely:
        for feat in coll['features']:
            feat['shapely'] = shape(feat['geometry'])
            feat['bbox'] = get_bbox(feat)
    return coll

def load_geojson_url(url, load_shapely=False):
    from urllib.request import urlopen
    import json
    geoj = json.loads(urlopen(url).read())
    if load_shapely:
        for feat in geoj['features']:
            feat['shapely'] = shape(feat['geometry'])
            feat['bbox'] = get_bbox(feat)
    return geoj

def morphology(arr, kernel, op, dtype):
    count = 0
    output = np.zeros(arr.shape, dtype=dtype)
    # should be much faster... 
    kernel_half = (kernel.shape[0] - 1) // 2
    for ky in range(kernel.shape[0]):
        ky_offset = kernel_half - ky
        y1 = max(ky_offset, 0)
        y2 = min(arr.shape[0]+ky_offset, arr.shape[0])
        #oy1 = max(ky, 0)
        #oy2 = min(ky+(y2-y1), output.shape[0])
        for kx in range(kernel.shape[1]):
            kx_offset = kernel_half - kx
            x1 = max(kx_offset, 0)
            x2 = min(arr.shape[1]+kx_offset, arr.shape[1])

            kval = kernel[ky,kx]
            kval_extract = kval * arr[y1:y2, x1:x2]

            slices = slice(ky,kval_extract.shape[0]), slice(kx,kval_extract.shape[1])
            #print(slices,output[slices].shape)
            output[slices] = op([output[slices], kval_extract[slices]])
    count = output.sum()
    return count,output

def burn(val, geometry, drawer, transform):
    from PIL import ImagePath
    geotype = geometry["type"]

    # set the coordspace to vectordata bbox
    a,b,c,d,e,f = transform

    # if fill val is None, then draw binary outline
    fill = val
    outline = 255 if val is None else None
    holefill = 0 if val is not None else None
    holeoutline = 255 if val is None else None
    #print ["burnmain",fill,outline,"burnhole",holefill,holeoutline]

    # make all multis so can treat all same
    coords = geometry["coordinates"]
    if not "Multi" in geotype:
        coords = [coords]

    # polygon, basic black fill, no outline
    if "Polygon" in geotype:
        for poly in coords:
            # exterior
            exterior = [tuple(p) for p in poly[0]]
            path = ImagePath.Path(exterior)
            #print list(path)[:10]
            path.transform((a,b,c,d,e,f))
            #print list(path)[:10]
            drawer.polygon(path, fill=fill, outline=outline)
            # holes
            if len(poly) > 1:
                for hole in poly[1:]:
                    hole = [tuple(p) for p in hole]
                    path = ImagePath.Path(hole)
                    path.transform((a,b,c,d,e,f))
                    drawer.polygon(path, fill=holefill, outline=holeoutline)
                    
    # line, 1 pixel line thickness
    elif "LineString" in geotype:
        for line in coords:
            line = [tuple(p) for p in line]
            path = ImagePath.Path(line)
            path.transform((a,b,c,d,e,f))
            drawer.line(path, fill=outline)
        
    # point, 1 pixel square size
    elif "Point" in geotype:
        for p in coords:
            p = tuple(p)
            path = ImagePath.Path([p])
            path.transform((a,b,c,d,e,f))
            drawer.point(path, fill=outline)

def show_surface(surf, minval=None, maxval=None, flipy=True, cmap=None):
    import matplotlib.pyplot as plt
    # setup plot
    plt.clf()
    ax = plt.gca()
    ax.set_aspect('equal', 'datalim')
    # add surface
    opts = {}
    if not cmap:
        cmap = 'PiYG'
    if cmap:
        opts['cmap'] = cmap
    if flipy:
        opts['origin'] = 'lower'
    plt.imshow(surf, **opts)
    if minval is not None or maxval is not None:
        plt.clim(minval, maxval)
    # show
    plt.colorbar()
    return plt.gcf()

def show_dataset(data, color_by=None, cmap=None, minval=None, maxval=None, flipy=True):
    import matplotlib.pyplot as plt
    # setup plot
    plt.clf()
    ax = plt.gca()
    ax.set_aspect('equal', 'datalim')
    # calc the colors
    import classypie as cp
    key = lambda f: f['properties'][color_by]
    if not cmap:
        cmap = 'PiYG'
    from matplotlib import cm
    cmap_obj = cm.get_cmap(cmap)
    color_by_colors = [tuple(col) for col in cmap_obj(range(cmap_obj.N+1))]
    breaks = cp.breaks([0,1], 'equal', classes=100,
                            minval=minval, maxval=maxval)
    classif = cp.Classifier(data['features'], breaks=breaks, #classes=100,
                            key=key, classvalues=color_by_colors,
                            minval=minval, maxval=maxval)
    # plot the data
    for feat,col in classif:
        for ring in iter_rings(feat['geometry']):
            ring = list(ring)
            if ring[0]!=ring[-1]: ring.append(ring[-1])
            x,y = zip(*ring)
            plt.fill(x, y, color=col)
            plt.plot(x, y, color='black', marker='')
    # show
    plt.colorbar(cm.ScalarMappable(cmap=cmap))
    return plt.gcf()

def show_datasets(data1, data2, surf=None, bbox=None, minval=None, maxval=None, flipy=True, cmap=None):
    import matplotlib.pyplot as plt
    # setup plot
    plt.clf()
    ax = plt.gca()
    ax.set_aspect('equal', 'datalim')
    # set the data bbox
    bboxes = [get_bbox(feat1) for feat1 in data1['features']] + [get_bbox(feat2) for feat2 in data2['features']]
    bbox = bbox_union(*bboxes)
    xmin,ymin,xmax,ymax = bbox
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax) # maybe depends on flipy? 
    # add surface
    if surf is not None:
        opts = {}
        if cmap:
            opts['cmap'] = cmap
        if flipy:
            opts['origin'] = 'lower'

        if bbox:
            xmin,ymin,xmax,ymax = bbox
            opts['extent'] = [xmin,xmax,ymax,ymin] # maybe depends on flipy? 
        else:
            raise Exception('bbox arg is required in order to know where the surf represents')            

        plt.imshow(surf, **opts)

        if minval is not None or maxval is not None:
            plt.clim(minval, maxval)
    # data1
    for feat in data1['features']:
        for ring in iter_rings(feat['geometry']):
            ring = list(ring)
            if ring[0]!=ring[-1]: ring.append(ring[-1])
            x,y = zip(*ring)
            plt.plot(x, y, color='tab:blue', marker='')
    # data2
    for feat in data2['features']:
        for ring in iter_rings(feat['geometry']):
            ring = list(ring)
            if ring[0]!=ring[-1]: ring.append(ring[-1])
            x,y = zip(*ring)
            plt.plot(x, y, color='tab:red', marker='')
    # show
    if surf is not None:
        plt.colorbar()
    if flipy:
        pass #ax.invert_yaxis() # maybe not necessary? 
    return plt.gcf()
        
def show_boundaries(boundaries, surf=None, bbox=None, flipy=True):
    import matplotlib.pyplot as plt
    # setup plot
    plt.clf()
    ax = plt.gca()
    ax.set_aspect('equal', 'datalim')
    if not bbox:
        bboxes = [b.uncertainty_bbox() for b in boundaries]
        bbox = bbox_union(*bboxes)
    xmin,ymin,xmax,ymax = bbox
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymax, ymin)
    # add in surface
    if surf is not None:
        extent = [xmin,xmax,ymax,ymin]
        plt.imshow(surf, extent=extent)
    # and boundaries
    for i,bnd in enumerate(boundaries):
        #print(i, bnd)
        # main shape
        color = None
        #for ring in iter_rings(bnd.geom):
        #    ring = list(ring)
        #    if ring[0]!=ring[-1]: ring.append(ring[-1])
        #    x,y = zip(*ring)
        #    p = plt.plot(x, y, color=color, linestyle='-.', linewidth=2, label=label) #, marker='o')
        #    color = p[0].get_color()
        # inner
        inner = shape(bnd.geom).buffer(-bnd.precision_range_max).__geo_interface__
        for ring in iter_rings(inner):
            ix,iy = zip(*ring)
            label = 'Boundary {}'.format(i+1) if color is None else None # only the first ring gets a label
            p = plt.plot(ix,iy, color=color, linestyle='--', label=label) #, marker='o')
            color = p[0].get_color()
        # outer
        outer = shape(bnd.geom).buffer(bnd.precision_range_max).__geo_interface__
        for ring in iter_rings(outer):
            ox,oy = zip(*ring)
            plt.plot(ox,oy, color=color, linestyle='--', label="_nolegend") #, marker='o')
    # show
    if surf is not None:
        plt.colorbar()
    if flipy:
        ax.invert_yaxis()
    plt.legend()
    return plt.gcf()
