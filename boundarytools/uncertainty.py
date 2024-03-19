import numpy as np
from datetime import datetime, timedelta
import math
import copy

from .utils import iter_rings, burn, morphology, bbox_union


def _line_dists(shapes):
    dists = []
    for geoj in shapes:
        for ring in iter_rings(geoj):
            ring = list(ring)
            if ring[0]!=ring[-1]: ring.append(ring[-1])
            for i in range(len(ring)-1):
                v1,v2 = ring[i:i+2]
                dx,dy = v1[0]-v2[0], v1[1]-v2[1]
                d = math.hypot(dx, dy)
                dists.append(d)
    return dists

def _line_resolution_min(dists):
    return sorted(dists) [ len(dists) // 100 ]

def _line_resolution_med(dists):
    return np.median(dists)


# classes

class Boundary(object):

    def __init__(self, geom, precision, precision_range_max=None):
        '''
        - geom is a geojson dict
        - precision is an equation string
            that determines the precision membership function of x as distance from geometry boundary (in meters).
        - precision_range_max is how far from the boundary the precision function reaches
        '''
        self.geom = geom
        self.precision = precision
        self.precision_range_max = precision_range_max or self.line_resolution_med()

        self._uncertainty_surface = {}

    def bbox(self, expand=None):
        xs,ys = [],[]
        for ring in iter_rings(self.geom):
            _xs,_ys = zip(*ring)
            xs.extend(_xs)
            ys.extend(_ys)
        bbox = min(xs), min(ys), max(xs), max(ys)

        if expand:
            return bbox[0]-expand, bbox[1]-expand, bbox[2]+expand, bbox[3]+expand, 
        else:
            return bbox

    def line_resolution_min(self):
        dists = _line_dists([self.geom])
        return _line_resolution_min(dists)

    def line_resolution_med(self):
        dists = _line_dists([self.geom])
        return _line_resolution_med(dists)

    def distance_surface(self, resolution, bbox, maxdist=None, boundary=True):
        # PIL VERSION
        # NOTE: adding 1 pixel padding on all sides, because the distance transform doesnt calc around edges
        # get coord to pixel drawing transform
        xscale,xskew,xoff = resolution,0,bbox[0]
        yskew,yscale,yoff = 0,resolution,bbox[1]
        affine = np.array([[xscale,xskew,xoff],
                           [yskew,yscale,yoff],
                           [0,0,1]])
        inv_affine = np.linalg.inv(affine).flatten()[:6]

        # create image and drawer
        from PIL import Image, ImageDraw
        width = (bbox[2]-bbox[0])
        height = (bbox[3]-bbox[1])
        imwidth = int(width / resolution)
        imheight = int(height / resolution)
        ###print('distance surface',resolution,bbox,imwidth,imheight)

        # exit early
        if (imwidth,imheight) == (1,1) or (imwidth,imheight) == (0,0): 
            return np.array([[1]])
        
        # draw the footprint outline
        outline = Image.new('L', (imwidth, imheight), 0)
        drawer = ImageDraw.Draw(outline)
        burn(None, self.geom, drawer, inv_affine)
        #outline.show()

        # calc pixel distances by iteratively dilating and incrementing a distance counter
        from PIL import ImageMorph, ImageMath, ImageOps
        outline = ImageOps.expand(outline, border=1, fill=0) # add padding for distance transform
        dilation_straight = ImageMorph.MorphOp(patterns = ["4:(... .0. .1.)->1"])
        initval = maxdist+1 if maxdist else 0 
        distances = Image.new('F', (imwidth+2, imheight+2), initval) # add padding for distance transform
        distances.paste(0, mask=outline)
        count = True
        cumul = outline
        dist = resolution
        while count:
            # dilate straight vertical and horizontal
            count1,grown1 = dilation_straight.apply(cumul)
            edge = ImageMath.eval('convert(grown ^ cumul, "1")', {'grown':grown1, 'cumul':cumul})
            distances.paste(dist, mask=edge)

            dist += resolution

            cumul = grown1 
            count = count1
            if maxdist and dist > maxdist:
                break

        # distances inside polygons are made negative
        if not boundary and 'Polygon' in self.geom['type']:
            # draw the footprint fill
            inside = Image.new('L', (imwidth, imheight), 0)
            drawer = ImageDraw.Draw(inside)
            burn(255, self.geom, drawer, inv_affine)
            inside = ImageOps.expand(inside, border=1, fill=0) # add padding for distance transform
            #inside.show()
            neg = ImageMath.eval('-d', {'d':distances})
            distances.paste(neg, mask=inside)

        # voila, return
        distances = np.array(distances)
        distances = distances[1:-1,1:-1] # strip away outer padding
        #Image.fromarray((distances/distances.max()*255).astype('uint8')).resize((1000,1000)).show()
        return distances

    def distance_surface2(self, resolution, bbox, maxdist=None, boundary=True):
        # NUMPY VERSION
        # get coord to pixel drawing transform
        xscale,xskew,xoff = resolution,0,bbox[0]-resolution # 1 pixel padding on each side
        yskew,yscale,yoff = 0,resolution,bbox[1]-resolution # 1 pixel padding on each side
        affine = np.array([[xscale,xskew,xoff],
                           [yskew,yscale,yoff],
                           [0,0,1]])
        inv_affine = np.linalg.inv(affine).flatten()[:6]

        # create image and drawer
        from PIL import Image, ImageDraw, ImagePath
        width = (bbox[2]-bbox[0]) + resolution * 2 # 1 pixel padding on each side
        height = (bbox[3]-bbox[1]) + resolution * 2 # 1 pixel padding on each side
        imwidth = int(width / resolution)
        imheight = int(height / resolution)
        ###print('distance surface',resolution,bbox,imwidth,imheight)

        # exit early
        if (imwidth,imheight) == (2,2) or (imwidth,imheight) == (3,3): # padding of 2 pixels, so 2,2 means size 0, 3,3 means size 1
            # sort of hacky, maybe not best option...
            return np.array([[1]])
        
        # draw the footprint outline
        outline = Image.new('L', (imwidth, imheight), 0)
        drawer = ImageDraw.Draw(outline)
        burn(None, self.geom, drawer, inv_affine)
        #outline.show()

        # calc pixel distances by iteratively dilating and incrementing a distance counter
        initval = maxdist+1 if maxdist else 0 
        distances = Image.new('F', (imwidth, imheight), initval)
        distances.paste(0, mask=outline)
        distances = np.array(distances)
        count = True
        cumul = np.array(outline, dtype=bool)
        dist = resolution
        dilate_kernel = np.array([[0,1,0],[1,1,1],[0,1,0]], dtype=bool)
        prevcount = None
        while count != prevcount:
            # dilate straight vertical and horizontal
            prevcount = count
            count1,grown1 = morphology(cumul, dilate_kernel, lambda res: np.maximum(*res), dtype=bool) 
            edge = grown1 ^ cumul
            distances[edge] = dist
            
            dist += resolution
            cumul = grown1 
            count = count1 
            if maxdist and dist > maxdist:
                break

        # distances inside polygons are made negative
        if not boundary and 'Polygon' in self.geom['type']:
            # draw the footprint fill
            inside = Image.new('L', (imwidth, imheight), 0)
            drawer = ImageDraw.Draw(inside)
            burn(255, self.geom, drawer, inv_affine)
            #inside.show()
            neg = -distances
            inside = np.array(inside, dtype=bool)
            distances[inside] = neg[inside]

        # voila, return
        distances = distances[1:-1,1:-1] # strip away outer padding
        #Image.fromarray((distances/distances.max()*255).astype('uint8')).resize((1000,1000)).show()
        return distances

    def precision_kernel(self, resolution, maxdist=None):
        '''Generates distretized 2D precision kernel at a given resolution,
        representing the membership of binary yes/no overlap at arbitrary point.
        NOTE: BEWARE OF ALIASING EFFECT IF THE BBOX DONT FALL EXACTLY ON THE GEOM...
        TODO: NOT SURE IF SHOULD ALLOW maxdist, SINCE FOOTPRINT ALREADY HAS precision_range_max...
        '''

        # build distance matrix (for one quadrant of a kernel)
        maxdist = maxdist or self.precision_range_max
        kernel_center = {'type':'Point', 'coordinates':[0,0]}
        kernel_footprint = Boundary(kernel_center, self.precision, self.precision_range_max)
        kernel_bbox = [0,0,maxdist,maxdist]
        #cornersurf = kernel_footprint.distance_surface(resolution, kernel_bbox, boundary=True)
        cornersurf = kernel_footprint.precision_surface(resolution, kernel_bbox, maxdist)
        cornerkernel = np.array(cornersurf)

        # expand kernel to all four quadrants
        if cornerkernel.shape[0] > 1:
            bottom_right = cornerkernel
            #print('quad',bottom_right)
            bottom_left = np.fliplr(bottom_right[:,1:]) # minus the kernel center
            bottom = np.concatenate([bottom_left, bottom_right], axis=1)
            #print('bottom',bottom)
            top = np.flipud(bottom[1:,:]) # minus the kernel center
            kernel = np.concatenate([top, bottom], axis=0)
            #print('full kernel',kernel.shape,kernel)
        else:
            kernel = np.array([[1.0]])
        ###print('kernel footprint',resolution,kernel_bbox,kernel.shape)

        # normalize so whole distribution sums to 1? 
        #summ = kernel.sum()
        #kernel = kernel / summ
        #print('kernel',kernel)
        
        return kernel

    def precision_surface(self, resolution, bbox=None, maxdist=None):
        '''Generates distretized 2D accuracy matrix at a given resolution,
        representing the membership of binary yes/no overlap at each point.
        NOTE: BEWARE OF ALIASING EFFECT IF THE BBOX DONT FALL EXACTLY ON THE GEOM...
        '''
        # create distance grid
        maxdist = maxdist or self.precision_range_max
        bbox = bbox or self.bbox(maxdist)
        distgrid = self.distance_surface(resolution, bbox, boundary=True, maxdist=maxdist)

        # eval accuracy equation of the distance grid
        context = {'sqrt':np.sqrt, 'pi':np.pi, 'exp': np.exp}
        context['x'] = distgrid
        outgrid = eval(self.precision, {}, context)

        # clamp to 0-1
        outgrid = np.minimum(outgrid, 1)
        outgrid = np.maximum(outgrid, 0)
            
        return outgrid

    def uncertainty_surface(self, resolution, bbox=None, maxdist=None):
        maxdist = maxdist or self.precision_range_max
        bbox = bbox or self.bbox(maxdist)

        # get from cache? 
        k = maxdist,bbox
        if k in self._uncertainty_surface:
            return copy.deepcopy(self._uncertainty_surface[k])

        # get coord to pixel drawing transform
        xscale,xskew,xoff = resolution,0,bbox[0]-resolution # 1 pixel padding on each side
        yskew,yscale,yoff = 0,resolution,bbox[1]-resolution # 1 pixel padding on each side
        affine = np.array([[xscale,xskew,xoff],
                           [yskew,yscale,yoff],
                           [0,0,1]])
        inv_affine = np.linalg.inv(affine).flatten()[:6]

        # create image and drawer
        from PIL import Image, ImageDraw, ImagePath
        width = (bbox[2]-bbox[0]) + resolution * 2 # 1 pixel padding on each side
        height = (bbox[3]-bbox[1]) + resolution * 2 # 1 pixel padding on each side
        imwidth = int(width / resolution)
        imheight = int(height / resolution)
        ###print('distance surface',resolution,bbox,imwidth,imheight)

        # exit early
        if (imwidth,imheight) == (2,2) or (imwidth,imheight) == (3,3): # padding of 2 pixels, so 2,2 means size 0, 3,3 means size 1
            # sort of hacky, maybe not best option...
            return np.array([[1]])
        
        # draw the footprint insides
        inside = Image.new('L', (imwidth, imheight), 0)
        drawer = ImageDraw.Draw(inside)
        burn(255, self.geom, drawer, inv_affine)
        #outline.show()

        inside_surf = np.array(inside)
        precision_kernel = self.precision_kernel(resolution, maxdist=maxdist)

        # normalize so whole kernel distribution sums to 1? 
        summ = precision_kernel.sum()
        precision_kernel = precision_kernel / summ
        precision_kernel = precision_kernel / 255.0 # not sure why this is needed or if it truly works...

        kernel_offset = (precision_kernel.shape[0] - 1) // 2

        # convolve kernel over inside surface
        output = np.zeros(inside_surf.shape, dtype=np.float32)

        def convolute():
            # slow proper version: iterate each pixel, and paste sum the kernel
            for y in range(inside_surf.shape[0]):
                for x in range(inside_surf.shape[1]):
                    prec = inside_surf[y,x]
                    if not prec:
                        continue
                    # weight kernel by precision membership
                    weighted_kernel = precision_kernel * prec
                    # sum add weighted kernel to output uncertainty
                    b1 = output
                    b2 = weighted_kernel
                    #print('b1',b1.shape)
                    #print('b2',b2.shape)

                    pos_v, pos_h = y-kernel_offset,x-kernel_offset
                    v_range1 = slice(max(0, pos_v), max(min(pos_v + b2.shape[0], b1.shape[0]), 0))
                    h_range1 = slice(max(0, pos_h), max(min(pos_h + b2.shape[1], b1.shape[1]), 0))

                    v_range2 = slice(max(0, -pos_v), min(-pos_v + b1.shape[0], b2.shape[0]))
                    h_range2 = slice(max(0, -pos_h), min(-pos_h + b1.shape[1], b2.shape[1]))

                    b1[v_range1, h_range1] += b2[v_range2, h_range2]

                    #output[y-kernel_offset-1:y+kernel_offset, x-kernel_offset-1:x+kernel_offset] += weighted_kernel #.flatten()

        def convolute2(output=output):
            # faster version: iterate only the kernel values, apply to entire image with offsets
            # not working correctly yet, cant get the offsets and slices right
            kernel_half = (precision_kernel.shape[0] - 1) // 2
            for ky in range(precision_kernel.shape[0]):
                ky_offset = kernel_half - ky
                y1 = max(ky_offset, 0)
                y2 = min(inside_surf.shape[0]+ky_offset, inside_surf.shape[0])
                #oy1 = max(ky, 0)
                #oy2 = min(ky+(y2-y1), output.shape[0])
                for kx in range(precision_kernel.shape[1]):
                    kx_offset = kernel_half - kx
                    x1 = max(kx_offset, 0)
                    x2 = min(inside_surf.shape[1]+kx_offset, inside_surf.shape[1])
                    #ox1 = max(kx, 0)
                    #ox2 = min(kx+(x2-x1), output.shape[1])
                    #print(ky,kx, (y1,y2), (x1,x2))
                    kval = precision_kernel[ky,kx]
                    kval_arr_weighted = kval * inside_surf[y1:y2, x1:x2]
                    #print(output.shape)
                    #print(kval_arr_weighted.shape)
                    # TODO: slices here are being REAL FUNKY, hacky solution, NOT CORRECT
                    #slices = slice(oy1,oy2),slice(ox1,ox2)
                    slices = slice(ky,kval_arr_weighted.shape[0]), slice(kx,kval_arr_weighted.shape[1])
                    #print(slices,output[slices].shape)
                    output[slices] += kval_arr_weighted[slices]
                    #print(output)
            
        convolute()

        # cache results
        self._uncertainty_surface[k] = copy.deepcopy(output)

        return output

    def show(self, surf=None, bbox=None, flipy=True):
        import matplotlib.pyplot as plt
        from shapely.geometry import asShape
        # setup plot
        plt.clf()
        ax = plt.gca()
        ax.set_aspect('equal', 'datalim')
        if not bbox:
            bbox = self.bbox(self.precision_range_max)
        xmin,ymin,xmax,ymax = bbox
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymax, ymin)
        # add in surface
        if surf is not None:
            extent = [xmin,xmax,ymax,ymin]
            plt.imshow(surf, extent=extent)
        # main shape
        for ring in iter_rings(self.geom):
            ring = list(ring)
            if ring[0]!=ring[-1]: ring.append(ring[-1])
            x,y = zip(*ring)
            plt.plot(x, y, color='tab:blue') #, marker='o')
        # outer range
        buf = self.precision_range_max
        outer = asShape(self.geom).buffer(buf).__geo_interface__
        for ring in iter_rings(outer):
            ring = list(ring)
            if ring[0]!=ring[-1]: ring.append(ring[-1])
            x,y = zip(*ring)
            plt.plot(x, y, color='black')
        # inner range
        buf = -buf
        inner = asShape(self.geom).buffer(buf).__geo_interface__
        for ring in iter_rings(inner):
            ring = list(ring)
            if ring[0]!=ring[-1]: ring.append(ring[-1])
            x,y = zip(*ring)
            plt.plot(x, y, color='black')
        # show
        if surf is not None:
            plt.colorbar()
        if flipy:
            ax.invert_yaxis()
        plt.show()

    #################

    def uncertainty_bbox(self):
        if not hasattr(self, '_uncertainty_bbox'):
            bbox = self.bbox(self.precision_range_max)
            self._uncertainty_bbox = list(bbox)
        return self._uncertainty_bbox

    def bbox_intersection(self, other):
        bbox1 = self.uncertainty_bbox()
        bbox2 = other.uncertainty_bbox()

        if bbox2[0] <= bbox1[2] and bbox2[2] >= bbox1[0] and bbox2[1] <= bbox1[3] and bbox2[3] >= bbox1[1]:
            xmin =  max(bbox1[0],bbox2[0])
            ymin =  max(bbox1[1],bbox2[1])
            xmax =  min(bbox1[2],bbox2[2])
            ymax =  min(bbox1[3],bbox2[3])
            isec_bbox = xmin,ymin,xmax,ymax
            return isec_bbox

    def bbox_union(self, other):
        bbox1 = self.uncertainty_bbox()
        bbox2 = other.uncertainty_bbox()
        # get combined bbox
        xmin =  min(bbox1[0],bbox2[0])
        ymin =  min(bbox1[1],bbox2[1])
        xmax =  max(bbox1[2],bbox2[2])
        ymax =  max(bbox1[3],bbox2[3])
        return xmin,ymin,xmax,ymax

    def overlap_surface(self, other, resolution=None, bbox=None):
        if not bbox:
            bbox = self.bbox_union(other)
        if not resolution:
            dx = bbox[2]-bbox[0]
            dy = bbox[3]-bbox[1]
            import math
            diag = math.hypot(dx, dy)
            resolution = diag / 300.0

        ###print('surf1')
        surf1 = self.uncertainty_surface(resolution, bbox)
        ###print('surf2')
        surf2 = other.uncertainty_surface(resolution, bbox)
        #print('surf1')#,surf1.shape,surf1)
        #print('surf2')#,surf2.shape,surf2)

        import numpy as np
        ###print('fuzzy_min',surf1.shape,surf2.shape)
        fuzzy_min = np.minimum(surf1, surf2)
        #print('fuzzy_min',fuzzy_min.shape)#,fuzzy_min)
        return surf1,surf2,fuzzy_min

    def difference_surface(self, other, resolution=None, bbox=None):
        if not bbox:
            bbox = self.bbox_union(other)
        if not resolution:
            dx = bbox[2]-bbox[0]
            dy = bbox[3]-bbox[1]
            import math
            diag = math.hypot(dx, dy)
            resolution = diag / 300.0

        print('surf1')
        surf1 = self.uncertainty_surface(resolution, bbox)
        print('surf2')
        surf2 = other.uncertainty_surface(resolution, bbox)
        #print('surf1')#,surf1.shape,surf1)
        #print('surf2')#,surf2.shape,surf2)

        import numpy as np
        print('fuzzy_min',surf1.shape,surf2.shape)
        fuzzy_diff = abs(surf1 - surf2)
        #print('fuzzy_min',fuzzy_min.shape)#,fuzzy_min)
        return surf1,surf2,fuzzy_diff

    def similarity(self, other, resolution=None, bbox=None):
        bbox1 = self.uncertainty_bbox()
        bbox2 = other.uncertainty_bbox()
        xmin1,ymin1,xmax1,ymax1 = bbox1
        xmin2,ymin2,xmax2,ymax2 = bbox2
        boxoverlap = (xmin1 <= xmax2 and xmax1 >= xmin2 and ymin1 <= ymax2 and ymax1 >= ymin2)
        if not boxoverlap:
            return {'equality':0, 'within':0, 'contains':0}
        
        if not bbox:
            bbox = self.bbox_union(other)
        if not resolution:
            dx = bbox[2]-bbox[0]
            dy = bbox[3]-bbox[1]
            import math
            diag = math.hypot(dx, dy)
            resolution = diag / 300.0

        surf1,surf2,fuzzy_min = self.overlap_surface(other, resolution, bbox)
        #from PIL import Image
        #Image.fromarray((surf1/surf1.max()*255).astype('uint8')).resize((1000,1000)).show()
        #Image.fromarray((surf2/surf2.max()*255).astype('uint8')).resize((1000,1000)).show()
        #Image.fromarray((fuzzy_min/fuzzy_min.max()*255).astype('uint8')).resize((1000,1000)).show()

        import numpy as np
        similarity = {}

        fuzzy_max = np.maximum(surf1, surf2)
        fuzzy_min_sum = fuzzy_min.sum()
        similarity['equality'] = float( fuzzy_min_sum / fuzzy_max.sum() )
        similarity['within'] = float( fuzzy_min_sum / surf1.sum() )
        similarity['contains'] = float( fuzzy_min_sum / surf2.sum() )
        
        return similarity

class NormalBoundary(Boundary):
    def __init__(self, geom, mean=None, stdev=None):
        '''
        - geom is a geojson dict
        '''
        self.geom = geom

        if stdev is None:
            dists = _line_dists([self.geom])
            lineres = _line_resolution_med(dists) # median seems to be the best match for the normal distribution
            stdev = lineres
        if mean is None:
            mean = 0

        self.precision = '1 / (sqrt(2*pi)*{sig}) * exp((-((x-{mu})/{sig})**2)/2.0) / 100.0'.format(mu=mean, sig=stdev)
        self.precision_range_max = stdev * 3 # if the median line resolution = std of the normal distribution, then 3x std = full range

        self._uncertainty_surface = {}


# functions

def probability_inside(boundaries, resolution=None, bbox=None):
    # resolution
    if not resolution:
        resolution = min(bnd.precision_range_max for bnd in boundaries) / 4.0
        print('resolution', resolution)

    # get bboxes
    if not bbox:
        bboxes = [b.uncertainty_bbox() for b in boundaries]
        bbox = bbox_union(*bboxes)

    # probability of being inside any of the boundaries
    cumul = None
    for bnd in boundaries:
        inside = bnd.uncertainty_surface(resolution, bbox=bbox)
        #bnd.show(inside, bbox=bbox)
        # get the probability inside our boundary AND not inside any other boundary
        certinside = inside
        for bnd2 in boundaries:
            if bnd2 is bnd:
                continue
            insideother = bnd2.uncertainty_surface(resolution, bbox=bbox)
            notinsideother = 1 - insideother
            certinside *= notinsideother
        #bnd.show(certinside, bbox=bbox)
        # add to cumul
        if cumul is None:
            cumul = certinside
        else:
            cumul += certinside
        #bnd.show(cumul, bbox=bbox)

    #boundarytools.utils.show_boundaries(boundaries, cumul, bbox=bbox)
    return cumul

def crisp_footprints(boundaries, resolution=None, bbox=None):
    # resolution
    if not resolution:
        resolution = min(bnd.precision_range_max for bnd in boundaries) / 4.0
        print('resolution', resolution)

    # get bboxes
    if not bbox:
        bboxes = [b.uncertainty_bbox() for b in boundaries]
        bbox = bbox_union(*bboxes)

    # pixels that are crisply within any of the boundaries
    cumul = None
    for bnd in boundaries:
        inside = bnd.uncertainty_surface(resolution, bbox=bbox)
        #bnd.show(certinside, bbox=bbox)
        # add to cumul
        if cumul is None:
            cumul = inside
        else:
            cumul = np.maximum(cumul, inside)
        #bnd.show(cumul, bbox=bbox)

    # crispify (0.5 is the breaking point between inside and outside)
    cumul = np.where(cumul>0.5, 1, 0)

    #boundarytools.utils.show_boundaries(boundaries, cumul, bbox=bbox)
    return cumul
