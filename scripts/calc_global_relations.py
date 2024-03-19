# imports
import boundarytools
import numpy as np
from pyproj import Geod
geod = Geod(ellps="WGS84")
from shapely.geometry import shape
from shapely.validation import make_valid

import os
import sys
import json
from urllib.request import urlopen
import csv
import traceback
import multiprocessing as mp
import datetime
from time import time
import tempfile

# params

BRANCH = 'gadm4'
OUTPUT_DIR = 'global_relations'
SOURCES = ['geoBoundaries (Open)', 'GADM v4.0.4', 'OSM-Boundaries', 'UN SALB', 'OCHA', 'Natural Earth v5.0.1']
IGNORE_SOURCES = []
MAXPROCS = 2
COUNTRIES = [] #'CAN','CHL']


def source_in_sources(source, sources):
    # to allow for some variation, match if source matches the beginning of any in sources
    matches = (source.startswith(comp) for comp in sources)
    return any(matches)


def loop_country_levels():
    url = f'https://raw.githubusercontent.com/wmgeolab/geoContrast/{BRANCH}/releaseData/geoContrast-meta.csv'
    raw = urlopen(url).read().decode('utf8')
    reader = csv.DictReader(raw.split('\n'))
    def key(row):
        return row['boundaryISO'], int(row['boundaryType'][-1])
    iso_levels = sorted(set([key(row) for row in reader]))
    for iso,level in iso_levels:
        yield iso,level


def get_country_level_areas(country, level):
    # first get all possible sources
    sourcedict = boundarytools.utils.find_geocontrast_sources(country, level, branch=BRANCH)
    sourcedict = dict([(s,u.replace('/stable/',f'/{BRANCH}/')) for s,u in sourcedict.items()]) # override with chosen branch
    print('available sources:', sourcedict.keys())

    # define how to load and prep each source 
    def load_source(url):
        print('loading', url)
        coll = boundarytools.utils.load_topojson_url(url, load_shapely=True)
        
        # delete geojson repr to reduce memory (only need shapely)
        for feat in coll['features']:
            del feat['geometry']

        # simplify and validate
        print('simplifying and validating')
        for feat in coll['features']:
            # simplify
            simplified = feat['shapely'].simplify(0.01, True) # ca 1km

            # try simple fix if invalid
            if not simplified.is_valid:
                simplified = simplified.buffer(0)

            # try again if the fix didn't work (eg buffer 0 can lead to empty geom)
            if not simplified.is_valid or simplified.is_empty:
                # original shape was probably invalid
                # force validate
                shp = feat['shapely']
                shp = make_valid(shp)

                # fix by splitting original into linestring and polygonizing the lines
                #simplified = feat['shapely']
                #if simplified.geom_type == 'MultiPolygon':
                #    polys = simplified.geoms
                #else:
                #    polys = [simplified]
                #lines = []
                #for poly in polys:
                #    for ring in [poly.exterior] + list(poly.interiors):
                #        line = ring
                #        lines.append(line)
                        
                # then extract only the valid polygons from these lines
                #polys = polygonize(lines)
                #simplified = MultiPolygon(polys)
                #simplified = simplified

                # simplify the validated geom
                simplified = shp.simplify(0.01, True) # ca 1km

            # overwrite with simplified geom
            assert simplified.is_valid and not simplified.is_empty and simplified.geom_type != 'GeometryCollection'
            feat['shapely'] = simplified

        # validating
        #print('validating')
        #for feat in coll['features']:
        #    if not feat['shapely'].is_valid:
        #        feat['shapely'] = feat['shapely'].buffer(0)
        #        assert not feat['shapely'].is_empty
        return coll

    # download and prep all sources in local temp files
    print('downloading and prepping all source data')
    localfiles = {}
    for source,url in sourcedict.items():
        if SOURCES and not source_in_sources(source, SOURCES): continue
        elif source in IGNORE_SOURCES: continue
        print('-->', source)
        coll = load_source(url)
        for f in coll['features']:
            f['geometry'] = f['shapely'].__geo_interface__
            del f['shapely']
        tmp = tempfile.NamedTemporaryFile(mode='w+')
        tmp.write(json.dumps(coll))
        localfiles[url] = tmp
        
    def load_local_source(url):
        tmp = localfiles[url]
        tmp.seek(0)
        coll = json.loads(tmp.read())
        for f in coll['features']:
            f['shapely'] = shape(f['geometry'])
        return coll

    # calc relations bw all pairs of sources
    source_areas = {}
    source_results_matrix = {}
    source_errors_matrix = {}
    #for source1,coll1 in sourcedata.items():
    for source1,url1 in sourcedict.items():
        if SOURCES and not source_in_sources(source1, SOURCES): continue
        elif source1 in IGNORE_SOURCES: continue
        print('-->', source1)
        coll1 = load_local_source(url1)

        print('calc areas', url1)
        areas = []
        for feat in coll1['features']:
            area,perim = geojson_area_perimeter(feat['shapely'].__geo_interface__)
            areas.append(round(area, 1))
        source_areas[source1] = areas
        
        source_results_row = {}
        source_errors_row = {}
        
        for source2,url2 in sourcedict.items():
            if source2 == source1: continue

            if SOURCES and not source_in_sources(source2, SOURCES): continue
            elif source2 in IGNORE_SOURCES: continue
            print('-->', source1, 'vs', source2)
            coll2 = load_local_source(url2)

            # calc feature pair areas
            print('calculating pairwise feature relations')
            feature_pair_areas,feature_pair_errors = get_feature_pair_areas(coll1, coll2)
            source_results_row[source2] = feature_pair_areas
            if feature_pair_errors:
                source_errors_row[source2] = feature_pair_errors
                
            # calc simil stats
            #def similarity(feat1, feat2):
            #    isec_area = feat1['shapely'].intersection(feat2['shapely']).area
            #    union_area = feat1['shapely'].union(feat2['shapely']).area
            #    equality = isec_area / union_area
            #    within = isec_area / feat1['shapely'].area
            #    contains = isec_area / feat2['shapely'].area
            #    stats = {'equality':equality, 'within':within, 'contains':contains}
            #    print(stats)
            #    return stats
            #feat1,feat2 = coll1['features'][0], coll2['features'][0]
            #source_results_row[source2] = similarity(feat1, feat2)
            # first from persepctive of coll1
            #feat1,feat2 = coll1['features'][0], coll2['features'][0]
            #feat1['similarity'] = similarity(feat1, feat2)
            # then from perspective of coll2
            #feat2['similarity'] = similarity(feat2, feat1)
            
        source_results_matrix[source1] = source_results_row
        if source_errors_row:
            source_errors_matrix[source1] = source_errors_row

    # close and clear all local temp files
    for fobj in localfiles.values():
        fobj.close()

    return source_areas, source_results_matrix, source_errors_matrix


def get_feature_pair_areas(coll1, coll2):
    # calculate the relations
    # storing the full matrices is a lot of redundant information
    # and fairly quickly becomes impractical for eg >10,000 features, ie 10,000^2
    # instead store 3-tuples: row,col,[Adiffarea,Bdiffarea,ABarea]
    def bbox_overlap(bbox1, bbox2):
        xmin1,ymin1,xmax1,ymax1 = bbox1
        xmin2,ymin2,xmax2,ymax2 = bbox2
        overlap = (xmin1 <= xmax2 and xmax1 >= xmin2 and ymin1 <= ymax2 and ymax1 >= ymin2)
        return overlap
    relations = []
    errors = []
    for i1,feat1 in enumerate(coll1['features']):
        #print('feat1',i1)
        for i2,feat2 in enumerate(coll2['features']):
            try:
                if bbox_overlap(feat1['bbox'], feat2['bbox']) and feat1['shapely'].intersects(feat2['shapely']):
                    # AB intersection
                    geom = feat1['shapely'].intersection(feat2['shapely'])
                    AB,_ = geojson_area_perimeter(geom.__geo_interface__)
                    if AB == 0:
                        # intersection was probably just tangential, dont add
                        continue
                    # A diff
                    geom = feat1['shapely'].difference(feat2['shapely'])
                    Adiff,_ = geojson_area_perimeter(geom.__geo_interface__)
                    # B diff
                    geom = feat2['shapely'].difference(feat1['shapely'])
                    Bdiff,_ = geojson_area_perimeter(geom.__geo_interface__)
                    # add entry
                    entry = [i1,i2,[round(Adiff,1),round(Bdiff,1),round(AB,1)]]
                    relations.append(entry)
                else:
                    # no intersection, so Adiff and Bdiff are exclusive
                    # no need to store anything about the relations
                    pass
            except:
                err = traceback.format_exc()
                entry = [i1,i2,err]
                errors.append(entry)
    return relations, errors


def geojson_area_perimeter(geoj):
    # area may be negative if incorrect orientation
    # but the abs(area) will be correct as long as ext and holes
    # have opposite orientation
    polys = []
    if geoj['type'] == 'MultiPolygon':
        polys = geoj['coordinates']
    elif geoj['type'] == 'Polygon':
        polys = [geoj['coordinates']]
    elif geoj['type'] == 'GeometryCollection':
        polys = []
        for geom in geoj['geometries']:
            if geom['type'] == 'Polygon':
                polys.append(geom['coordinates'])
            elif geom['type'] == 'MultiPolygon':
                polys.extend(geom['coordinates'])
    if not polys:
        # intersection was prob only along a line or a point
        return 0.0, 0.0
        
    area = 0
    perim = 0
    for poly in polys:
        for i,ring in enumerate(poly):
            coords = np.array(ring)
            lons,lats = coords[:,0],coords[:,1]
            _area,_perim = geod.polygon_area_perimeter(lons, lats)
            if i == 0:
                # exterior
                area += abs(_area)
            else:
                # hole
                area -= abs(_area)
            perim += _perim
    return area, perim


def process(iso, level):
    # calc stats
    areas,relations,errors = get_country_level_areas(iso, level)

    # dump source feature areas as json
    filename = '{}-ADM{}-areas.json'.format(iso, level)
    path = os.path.join(OUTPUT_DIR, filename)
    with open(path, 'w', encoding='utf8') as w:
        json.dump(areas, w) #, indent=4)

    # dump feature pair relations as json
    filename = '{}-ADM{}-relations.json'.format(iso, level)
    path = os.path.join(OUTPUT_DIR, filename)
    with open(path, 'w', encoding='utf8') as w:
        json.dump(relations, w) #, indent=4)

    # dump errors as json if any
    if errors:
        filename = '{}-ADM{}-errors.json'.format(iso, level)
        path = os.path.join(OUTPUT_DIR, filename)
        with open(path, 'w', encoding='utf8') as w:
            w.write(json.dumps(errors, indent=4))


def process_logger(logfile, func, **kwargs):
    logger = open(logfile, 'w', encoding='utf8', buffering=1)
    sys.stdout = logger
    sys.stderr = logger
    print('PID:',os.getpid())
    print('time',datetime.datetime.now().isoformat())
    print('working path',os.path.abspath(''))
    # run it
    try:
        func(**kwargs)
    except:
        traceback.print_exc()
    # finish
    print('finished!')
    print('time',datetime.datetime.now().isoformat())


if __name__ == '__main__':
    maxprocs = MAXPROCS
    procs = []

    # loop country levels
    for iso,level in loop_country_levels():
        if COUNTRIES and iso not in COUNTRIES:
            continue
        
        print('')
        print(iso,level)

        # skip if stats already exists
        filename = '{}-ADM{}-relations.json'.format(iso, level)
        path = os.path.join(OUTPUT_DIR, filename)
        if os.path.lexists(path):
            print('already exists, skipping')
            continue

        # local
        #process(iso, level)
        #continue

        # multiprocessing
        logfile = '{}-ADM{}-log.txt'.format(iso, level)
        logpath = os.path.join(OUTPUT_DIR, logfile)
        p = mp.Process(target=process_logger,
                       args=[logpath,process],
                       kwargs={'iso':iso, 'level':level})
        p.start()
        procs.append((p,time()))

        # Wait in line
        while len(procs) >= maxprocs:
            for p,t in procs:
                if not p.is_alive():
                    procs.remove((p,t))
                elif time()-t > 60*60*2: # max 2 hr
                    p.terminate()
                    procs.remove((p,t))

    # waiting for last ones
    for p,t in procs:
        p.join()
