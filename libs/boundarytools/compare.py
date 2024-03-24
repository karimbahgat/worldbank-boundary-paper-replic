
# NOTE: it's possible the raster surfaces will be too computational/impractical
# so start with pure vector comparisons rather than raster ones
# ... 

from .uncertainty import _line_dists, _line_resolution_min

from .utils import iter_rings, get_shapely, get_bbox, bbox_union

import shapely
from shapely.geometry import asShape, LineString, MultiLineString

import numpy as np

def boundary_distances(geom1, geom2, interval_dist=None, signed_distances=False):
    # signed_distances is whether to consider the sign of the distances in terms of distance inside (negative)
    # ...and outside (positive) geom2. More computationally costly. 

    # determine min line resolution
    if not interval_dist:
        dists = _line_dists([geom1])
        res1 = _line_resolution_min(dists)
        dists = _line_dists([geom2])
        res2 = _line_resolution_min(dists)
        interval_dist = min(res1, res2)
        print('interval dist set to min res', interval_dist)

    # create shapely geoms
    if signed_distances:
        inside = asShape(geom2)
    rings = list(iter_rings(geom2))
    shp2 = MultiLineString(rings)

    # walk and sample points along boundary of geom1
    dists = []
    for ring in iter_rings(geom1):
        shp1 = LineString(ring)
        perim = shp1.length
        pos = 0
        while pos < perim:
            sample = shp1.interpolate(pos)
            dist = sample.distance(shp2)
            if signed_distances:
                _,nearest = shapely.ops.nearest_points(sample, inside)
                if nearest.intersects(inside):
                    dist *= -1
            dists.append(dist)
            pos += interval_dist

    return dists

def joint_probability_surface(boundaries, resolution=None, bbox=None):
    if not bbox:
        bboxes = [b.uncertainty_bbox() for b in boundaries]
        bbox = bbox_union(*bboxes)
    if not resolution:
        dx = bbox[2]-bbox[0]
        dy = bbox[3]-bbox[1]
        import math
        diag = math.hypot(dx, dy)
        resolution = diag / 300.0

    joint = boundaries[0].uncertainty_surface(resolution, bbox)
    for b in boundaries[1:]:
        surf = b.uncertainty_surface(resolution, bbox)
        joint *= surf
    return joint

def disjoint_probability_surface(boundaries, resolution=None, bbox=None):
    if not bbox:
        bboxes = [b.uncertainty_bbox() for b in boundaries]
        bbox = bbox_union(*bboxes)
    if not resolution:
        dx = bbox[2]-bbox[0]
        dy = bbox[3]-bbox[1]
        import math
        diag = math.hypot(dx, dy)
        resolution = diag / 300.0

    joint = 1 - boundaries[0].uncertainty_surface(resolution, bbox) # not
    for b in boundaries[1:]:
        surf = 1 - b.uncertainty_surface(resolution, bbox) # not
        joint *= surf
    return joint

def symmetric_difference_probability_surface(boundaries1, boundaries2, resolution=None, bbox=None):
    # resolution
    if not resolution:
        ranges = [bnd.precision_range_max for bnd in boundaries1] + [bnd.precision_range_max for bnd in boundaries2]
        resolution = min(ranges) / 2.0
        print('resolution', resolution)

    # get bboxes
    if not bbox:
        bboxes = [b.uncertainty_bbox() for b in boundaries1] + [b.uncertainty_bbox() for b in boundaries2]
        bbox = bbox_union(*bboxes)

    # probability of being inside any of the boundaries
    cumul = None
    for i,bnd1 in enumerate(boundaries1):
        print(i+1,bnd1)
        matches = []
        for i2,bnd2 in enumerate(boundaries2):
            #print('--> bnd2',i2+1,bnd2)
            isec = bnd1.bbox_intersection(bnd2)
            if isec:
                inside1 = bnd1.uncertainty_surface(resolution=resolution, bbox=bbox)
                inside2 = bnd2.uncertainty_surface(resolution=resolution, bbox=bbox)
                joint = inside1 * inside2
                union = np.maximum(inside1, inside2)

                equals = joint.sum() / union.sum()
                matches.append((equals,joint,inside1,inside2))

        # only compare to the most equal in source2
        equals,joint,inside1,inside2 = sorted(matches, key=lambda x: -x[0])[0]

        # calc symmetric difference
        symdiff = (inside1-joint) + (inside2-joint)
        #boundarytools.utils.show_boundaries([bnd1,bnd2], uniq, bbox=bbox)

        if cumul is None:
            cumul = symdiff
        else:
            cumul = np.maximum(cumul, symdiff)

    return cumul
    
def difference_probability_surface(boundaries1, boundaries2, resolution=None, bbox=None):
    # resolution
    if not resolution:
        ranges = [bnd.precision_range_max for bnd in boundaries1] + [bnd.precision_range_max for bnd in boundaries2]
        resolution = min(ranges) / 2.0
        print('resolution', resolution)

    # get bboxes
    if not bbox:
        bboxes = [b.uncertainty_bbox() for b in boundaries1] + [b.uncertainty_bbox() for b in boundaries2]
        bbox = bbox_union(*bboxes)

    # probability of being inside any of the boundaries
    cumul = None
    for i,bnd1 in enumerate(boundaries1):
        print(i+1,bnd1)
        matches = []
        for i2,bnd2 in enumerate(boundaries2):
            #print('--> bnd2',i2+1,bnd2)
            isec = bnd1.bbox_intersection(bnd2)
            if isec:
                inside1 = bnd1.uncertainty_surface(resolution=resolution, bbox=bbox)
                inside2 = bnd2.uncertainty_surface(resolution=resolution, bbox=bbox)
                joint = inside1 * inside2
                union = np.maximum(inside1, inside2)

                equals = joint.sum() / union.sum()
                matches.append((equals,joint,inside1,inside2))

        # only compare to the most equal in source2
        equals,joint,inside1,inside2 = sorted(matches, key=lambda x: -x[0])[0]

        # calc what's unique to source1
        uniq = inside1-joint 
        #boundarytools.utils.show_boundaries([bnd1,bnd2], uniq, bbox=bbox)

        if cumul is None:
            cumul = uniq
        else:
            cumul = np.maximum(cumul, uniq)

    return cumul

def similarity_surface(boundaries1, boundaries2, metric='equality', resolution=None, bbox=None):
    # resolution
    if not resolution:
        ranges = [bnd.precision_range_max for bnd in boundaries1] + [bnd.precision_range_max for bnd in boundaries2]
        resolution = min(ranges) / 2.0
        print('resolution', resolution)

    # get bboxes
    if not bbox:
        bboxes = [b.uncertainty_bbox() for b in boundaries1] + [b.uncertainty_bbox() for b in boundaries2]
        bbox = bbox_union(*bboxes)

    # similarity of overlapping boundaries for each pixel
    cumul = None
    maxprob = None
    for i,bnd1 in enumerate(boundaries1):
        print(i+1,bnd1)
        for i2,bnd2 in enumerate(boundaries2):
            #print('--> bnd2',i2+1,bnd2)
            isec = bnd1.bbox_intersection(bnd2)
            if isec:
                stats = bnd1.similarity(bnd2, resolution=resolution, bbox=bbox)
                val = stats[metric]
                
                #from boundarytools.utils import show_boundaries
                #show_boundaries([bnd1,bnd2], bbox=bbox)

                joint = joint_probability_surface([bnd1,bnd2], resolution, bbox)
                #from boundarytools.utils import show_surface
                #show_surface(joint, 0, 1)

                # set to statval where joint prob is higher than any previous prob
                # ie each pixel gets set to the pair similarity with highest prob
                if cumul is not None:
                    cumul = np.where(joint>maxprob, val, cumul)
                    maxprob = np.maximum(maxprob, joint)
                else:
                    cumul = np.where(joint>0, val, float('nan'))
                    maxprob = joint
                
                #show_surface(cumul, 0, 1)

    return cumul

def probability_feature_same_as_features(feat, features):
    geom1 = get_shapely(feat)
    bbox1 = get_bbox(feat)
    xmin1,ymin1,xmax1,ymax1 = bbox1

    cumprob = 0
    for feat2 in features:
        bbox2 = get_bbox(feat2)
        xmin2,ymin2,xmax2,ymax2 = bbox2
        boxoverlap = (xmin1 <= xmax2 and xmax1 >= xmin2 and ymin1 <= ymax2 and ymax1 >= ymin2)
        if boxoverlap:
            geom2 = get_shapely(feat2)
            try:
                isec = geom1.intersection(geom2)
            except:
                # unable to perform intersection, skip
                continue
            if not isec.is_empty:
                try:
                    union = geom1.union(geom2)
                except:
                    # unable to perform union, skip
                    continue
                isec_area = isec.area
                probability_intersects_geom2 = isec_area / geom1.area
                probability_is_same = isec_area / union.area
                probability_both = probability_intersects_geom2 * probability_is_same
                cumprob += probability_both

    return cumprob

def probabilities_aggregate(features, aggregate_to_feat, probability_field):
    geom1 = get_shapely(aggregate_to_feat)
    totarea = geom1.area
    
    cumprob = 0
    for feat in features:
        geom2 = get_shapely(feat)
        probability_intersects_geom2 = geom2.area / totarea
        probability_is_same = feat['properties'][probability_field]
        probability_both = probability_intersects_geom2 * probability_is_same
        cumprob += probability_both

    return cumprob

