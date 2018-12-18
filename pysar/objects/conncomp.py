############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2018, Zhang Yunjun                          #
# Author:  Zhang Yunjun                                    #
############################################################
# Recommend import:
#   from pysar.objects import connectComponent


import os
import time
import itertools
import numpy as np
from matplotlib import pyplot as plt
from scipy.sparse import csgraph as csg
from scipy.spatial import cKDTree
try:
    from skimage import measure, segmentation as seg, morphology as morph
except ImportError:
    raise ImportError('Could not import skimage!')
from .ramp import deramp


class connectComponent:
    """ Object for bridging connected components.
    
    Example:
        conncomp, atr = readfile.read_attribute('filt_fine.unw.conncomp')
        cc = connectComponent(conncomp=conncomp, metadata=atr)
        cc.label()
        cc.find_mst_bridge()
        
        unw = readfile.read('filt_fine.unw')[0]
        bdg_unw = cc.unwrap_conn_comp(unw, atr, ramp_type='linear')
        writefile.write(bdg_unw, 'bdg_filt_fine.unw', atr)
    """
    def __init__(self, conncomp, metadata):
        """Parameters: conncomp : 2D np.ndarray in np.bool_ format
                       metadata : dict, attributes
        """
        if type(conncomp).__module__ != np.__name__:
            raise ValueError('Input conncomp is not np.ndarray: {}'.format(type(conncomp).__module__))
        self.conncomp = conncomp
        self.metadata = metadata
        self.refY = int(self.metadata['REF_Y'])
        self.refX = int(self.metadata['REF_X'])
        self.length, self.width = self.conncomp.shape

    def label(self, min_area=1e4, erosion_size=5):
        label_image = measure.label(self.conncomp, connectivity=1)
        # take regions with large enough areas
        min_area = min(min_area, label_image.size * 3e-3)
        flag_slabel = np.bincount(label_image.flatten()) < min_area
        flag_slabel[0] = False
        label_small = np.where(flag_slabel)[0]
        for i in label_small:
            label_image[label_image == i] = 0
        # re-label
        self.labelImg, self.numLabel = measure.label(label_image, connectivity=1, return_num=True)
        # reference label, where reference pixel is located
        self.labelRef = self.labelImg[self.refY, self.refX]
        # find label boundaries to facilitate bridge finding
        self.find_boundary(erosion_size=erosion_size)
        return

    def find_boundary(self, erosion_size=5):
        self.labelErosion = morph.erosion(self.labelImg, morph.disk(erosion_size)).astype(np.uint8)
        self.labelBound = seg.find_boundaries(self.labelErosion, mode='thick').astype(np.uint8)
        self.labelBound *= self.labelErosion
        return

    def get_all_connections(self):
        regions = measure.regionprops(self.labelBound)
        if len(regions) < self.numLabel:
            msg = 'Some regions are too small --> lost during erosion.'
            msg += '\n1) decrease erosion_size value, or'
            msg += '\n2) increase min_area value.'
            raise ValueError(msg)

        trees = []
        for i in range(self.numLabel):
            trees.append(cKDTree(regions[i].coords))

        self.connDict = dict()
        self.distMat = np.zeros((self.numLabel, self.numLabel), dtype=np.float32)
        for i, j in itertools.combinations(range(self.numLabel), 2):
            # find shortest bridge
            dist, idx = trees[i].query(regions[j].coords)
            idx_min = np.argmin(dist)
            yxj = regions[j].coords[idx_min,:]
            yxi = regions[i].coords[idx[idx_min],:]
            dist_min = dist[idx_min]
            # save
            n0, n1 = str(i+1), str(j+1)
            conn = dict()
            conn[n0] = yxi
            conn[n1] = yxj
            conn['distance'] = dist_min
            self.connDict['{}{}'.format(n0, n1)] = conn
            self.distMat[i,j] = self.distMat[j,i] = dist_min
        return

    def find_mst_bridge(self):
        if not hasattr(self, 'distMat'):
            self.get_all_connections()

        # MST bridges with breadth_first_order
        distMatMst = csg.minimum_spanning_tree(self.distMat)
        succs, preds = csg.breadth_first_order(distMatMst, i_start=self.labelRef-1, directed=False)

        # save to self.bridges
        self.bridges = []
        for i in range(1, succs.size):
            n0 = preds[succs[i]] + 1
            n1 = succs[i] + 1
            # read conn
            nn = sorted([str(n0), str(n1)])
            conn = self.connDict['{}{}'.format(nn[0], nn[1])]
            y0, x0 = conn[str(n0)]
            y1, x1 = conn[str(n1)]
            # save bdg
            bridge = dict()
            bridge['x0'] = x0
            bridge['y0'] = y0
            bridge['x1'] = x1
            bridge['y1'] = y1
            bridge['label0'] = n0
            bridge['label1'] = n1
            self.bridges.append(bridge)
        self.num_bridge = len(self.bridges)
        return

    def unwrap_conn_comp(self, unw, atr, radius=50, ramp_type=None, print_msg=False):
        start_time = time.time()
        radius = int(min(radius, min(self.conncomp.shape)*0.05))

        ref_y, ref_x = int(atr['REF_Y']), int(atr['REF_X'])
        unw -= unw[ref_y, ref_x]
        unw = np.array(unw, dtype=np.float32)

        if ramp_type is not None:
            if print_msg:
                print('estimate a {} ramp'.format(ramp_type))
            ramp_mask = (self.labelImg == self.labelImg[ref_y, ref_x])
            unw, ramp = deramp(unw, ramp_mask, ramp_type, metadata=atr)

        for bridge in self.bridges:
            # get mask of AOI
            x0, y0 = bridge['x0'], bridge['y0']
            x1, y1 = bridge['x1'], bridge['y1']
            x00 = max(0, x0 - radius); x01 = min(self.width,  x0 + radius)
            y00 = max(0, y0 - radius); y01 = min(self.length, y0 + radius)
            x10 = max(0, x1 - radius); x11 = min(self.width,  x1 + radius)
            y10 = max(0, y1 - radius); y11 = min(self.length, y1 + radius)
            aoi_mask0 = np.zeros(self.labelImg.shape, dtype=np.bool_)
            aoi_mask1 = np.zeros(self.labelImg.shape, dtype=np.bool_)
            aoi_mask0[y00:y01, x00:x01] = True
            aoi_mask1[y10:y11, x10:x11] = True

            label_mask0 = self.labelImg == bridge['label0']
            label_mask1 = self.labelImg == bridge['label1']

            # get phase difference
            value0 = np.nanmedian(unw[aoi_mask0 * label_mask0])
            value1 = np.nanmedian(unw[aoi_mask1 * label_mask1])
            diff_value = value1 - value0

            # estimate integer number of phase jump
            num_jump = (np.abs(diff_value) + np.pi) // (2.*np.pi)
            if diff_value > 0:
                num_jump *= -1

            if print_msg:
                print('phase diff {}-{}: {:04.1f} rad --> num of jump: {}'.format(bridge['label1'],
                                                                                  bridge['label0'],
                                                                                  diff_value,
                                                                                  num_jump))

            # add phase jump
            unw[label_mask1] += 2.* np.pi * num_jump
        
        # add ramp back
        if ramp_type is not None:
            unw += ramp
        if print_msg:
            print('time used: {:.2f} secs.'.format(time.time()-start_time))
        return unw

    def plot_bridge(self):
        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=[8, 3], sharey=True)
        im = axs[0].imshow(self.conncomp)
        im = axs[1].imshow(self.labelImg)
        for bridge in self.bridges:
            axs[1].plot([bridge['x0'], bridge['x1']], [bridge['y0'], bridge['y1']], '-', lw=2)
        axs[1].plot(self.refX, self.refY, 'ks')
        plt.colorbar(im)
        plt.show()
        return


