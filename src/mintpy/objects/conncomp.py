"""Class / utilities for connected components."""
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2018                               #
############################################################
# Recommend import:
#   from mintpy.objects.conncomp import connectComponent


import itertools
import time

import numpy as np
from scipy.sparse import csgraph as csg
from scipy.spatial import cKDTree
from skimage import measure, morphology as morph, segmentation as seg

from mintpy.objects.ramp import deramp


######################################## utilities functions ##############################
def label_conn_comp(mask, min_area=2.5e3, erosion_size=5, print_msg=False):
    """Label / clean up the conn comp (mask)

    Parameters: mask         - 2D np.ndarray of bool/int
                min_area     - float, minimum region/area size
                erosion_size - int (odd number), size of erosion structure
                               set to 0 to turn it off.
    Returns:    label_img    - 2d np.ndarray of int, labeled array where all
                               connected regions are assigned the same value
                num_label    - int, number of labeled regions
    """

    # label
    label_img, num_label = measure.label(mask, connectivity=1, return_num=True)

    ## remove small regions
    min_area = min(min_area, label_img.size * 3e-3)
    if print_msg:
        print(f'remove regions with area < {int(min_area)}')
    mask = morph.remove_small_objects(label_img, min_size=min_area, connectivity=1)
    label_img[mask == 0] = 0
    # update label
    label_img, num_label = measure.label(label_img, connectivity=1, return_num=True) # re-label

    ## remove regions that would disappear after erosion
    # to ensure the consistency between label_img and label_bound
    if erosion_size > 0:
        erosion_structure = np.ones((erosion_size, erosion_size))
        label_erosion_img = morph.erosion(label_img, erosion_structure).astype(np.uint8)

        erosion_regions = measure.regionprops(label_erosion_img)
        if len(erosion_regions) < num_label:
            if print_msg:
                print('regions lost during morphological erosion operation:')

            label_erosion = [reg.label for reg in erosion_regions]
            for orig_reg in measure.regionprops(label_img):
                if orig_reg.label not in label_erosion:
                    label_img[label_img == orig_reg.label] = 0
                    if print_msg:
                        print('label: {}, area: {}, bbox: {}'.format(orig_reg.label,
                                                                     orig_reg.area,
                                                                     orig_reg.bbox))

            # update label
            label_img, num_label = measure.label(label_img, connectivity=1, return_num=True)

    return label_img, num_label


def label_boundary(label_img, num_label, erosion_size=5, print_msg=False):
    """Label the boundary of the labeled array

    Parameters: label_img    - 2d np.ndarray of int, labeled array where all connected regions are assigned the same value
                num_label    - int, number of labeled regions
    Returns:    label_img    - 2d np.ndarray of int, labeled array where all connected regions are assigned the same value
                num_label    - int, number of labeled regions
                label_bound  - 2d np.ndarrary of bool, where True represent a boundary pixel.
    """

    if erosion_size > 0:
        # remove regions that would disappear after erosion
        # to ensure the consistency between label_img and label_bound
        erosion_structure = np.ones((erosion_size, erosion_size))
        label_erosion_img = morph.erosion(label_img, erosion_structure).astype(np.uint8)

        erosion_regions = measure.regionprops(label_erosion_img)
        if len(erosion_regions) < num_label:
            if print_msg:
                print('regions lost during morphological erosion operation:')

            label_erosion = [reg.label for reg in erosion_regions]
            for orig_reg in measure.regionprops(label_img):
                if orig_reg.label not in label_erosion:
                    label_img[label_img == orig_reg.label] = 0
                    if print_msg:
                        print('label: {}, area: {}, bbox: {}'.format(orig_reg.label,
                                                                     orig_reg.area,
                                                                     orig_reg.bbox))

        # update label
        label_img, num_label = measure.label(label_img, connectivity=1, return_num=True) # re-label

    # get label boundaries to facilitate bridge finding
    label_bound = seg.find_boundaries(label_erosion_img, mode='thick').astype(np.uint8)
    label_bound *= label_erosion_img

    return label_img, num_label, label_bound



######################################## beginning of connectComponent class ##############################
class connectComponent:
    """ Object for bridging connected components.

    Example:
        unw_file = 'filt_fine.unw'
        # prepare connectComponent object
        atr = readfile.read_attribute(unw_file)
        conncomp = readfile.read(unw_file+'.conncomp')[0]
        cc = connectComponent(conncomp=conncomp, metadata=atr)
        cc.label()
        cc.find_mst_bridge()

        # run bridging
        unw = readfile.read(unw_file)[0]
        bdg_unw = cc.unwrap_conn_comp(unw, ramp_type='linear')

        # write output file
        writefile.write(bdg_unw, 'bdg_'+unw_file, atr)
    """

    def __init__(self, conncomp, metadata):
        """Parameters: conncomp : 2D np.ndarray in np.bool_ format
                       metadata : dict, attributes
        """
        if type(conncomp).__module__ != np.__name__:
            raise ValueError(f'Input conncomp is not np.ndarray: {type(conncomp).__module__}')
        self.conncomp = conncomp
        self.metadata = metadata
        if 'REF_Y' in metadata.keys():
            self.refY = int(self.metadata['REF_Y'])
            self.refX = int(self.metadata['REF_X'])
        else:
            self.refY = None
            self.refX = None
        self.length, self.width = self.conncomp.shape


    def label(self, min_area=2.5e3, erosion_size=5, print_msg=False):
        """ Label the connected components
        Returns: self.labelImg   - 2D np.ndarray in int64 to mask areas to be corrected
                 self.labelBound - 2D np.ndarray in uint8 for label boundaries to find bridges
        """
        # label conn comp
        (self.labelImg,
         self.numLabel) = label_conn_comp(self.conncomp,
                                          min_area=min_area,
                                          print_msg=print_msg)

        # label conn comp boundaries
        (self.labelImg,
         self.numLabel,
         self.labelBound) = label_boundary(self.labelImg,
                                           self.numLabel,
                                           erosion_size=erosion_size,
                                           print_msg=print_msg)

        # reference label (ref_y/x or the largest one)
        if self.refY is not None:
            self.labelRef = self.labelImg[self.refY, self.refX]
            if self.labelRef == 0:
                raise ValueError('input reference point is NOT included in the connectComponent.')
        else:
            regions = measure.regionprops(self.labelImg)
            idx = np.argmax([region.area for region in regions])
            self.labelRef = regions[idx].label
        return


    def get_all_bridge(self):
        """ Search all possible connections among labeled regions
        Returns:    connDict : dict of connection, i.e.:
                        {'1_2': {'1': array([1232,  345]),
                                 '2': array([868, 239]),
                                 'distance': 379.1200337623956},
                         '1_3': {'1': array([1232,  345]),
                                 '3': array([1089,  191]),
                                 'distance': 210.1547049199708},
                         '1_4': {'1': array([1204, 1143]),
                                 '4': array([1217, 1157]),
                                 'distance': 19.1049731745428},
                         '1_5': {'1': array([1263,  557]),
                                 '5': array([1270,  565]),
                                 'distance': 10.63014581273465},
                         '2_3': {'2': array([868, 239]),
                                 '3': array([891, 249]),
                                 'distance': 25.079872407968907},
                         '2_4': {'2': array([868, 239]),
                                 '4': array([1273, 1103]),
                                 'distance': 954.2122405419037},
                         '2_5': {'2': array([868, 239]),
                                 '5': array([1269,  566]),
                                 'distance': 517.4263232577175},
                         '3_4': {'3': array([996, 275]),
                                 '4': array([1319, 1085]),
                                 'distance': 872.0258023705492},
                         '3_5': {'3': array([1015,  264]),
                                 '5': array([1289,  545]),
                                 'distance': 392.4754769409167},
                         '4_5': {'4': array([1319, 1085]),
                                 '5': array([1305,  670]),
                                 'distance': 415.2360774306587}
                        }
                    distMat : 2D np.array in size of (nLabel, nLabel), i.e.:
                        array([[  0.      , 379.12003 , 210.15471 ,  19.104973,  10.630146],
                               [379.12003 ,   0.      ,  25.079872, 954.2122  , 517.42633 ],
                               [210.15471 ,  25.079872,   0.      , 872.0258  , 392.47546 ],
                               [ 19.104973, 954.2122  , 872.0258  ,   0.      , 415.23608 ],
                               [ 10.630146, 517.42633 , 392.47546 , 415.23608 ,   0.      ]],
                              dtype=float32)
        """
        regions = measure.regionprops(self.labelBound)

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
            self.connDict[f'{n0}_{n1}'] = conn
            self.distMat[i,j] = self.distMat[j,i] = dist_min
        return self.connDict, self.distMat


    def find_mst_bridge(self):
        """Search for bridges to connect all labeled areas using the minimum spanning tree algorithm
        Returns:    bridges : list of dict, i.e.:
                        [{'label0': 1, 'label1': 3, 'x0': 345, 'x1': 191, 'y0': 1232, 'y1': 1089},
                         {'label0': 1, 'label1': 4, 'x0': 1143, 'x1': 1157, 'y0': 1204, 'y1': 1217},
                         {'label0': 1, 'label1': 5, 'x0': 557, 'x1': 565, 'y0': 1263, 'y1': 1270},
                         {'label0': 3, 'label1': 2, 'x0': 249, 'x1': 239, 'y0': 891, 'y1': 868}]
        """
        if not hasattr(self, 'distMat'):
            self.get_all_bridge()

        # MST bridges with breadth_first_order
        distMatMst = csg.minimum_spanning_tree(self.distMat)
        succs, preds = csg.breadth_first_order(distMatMst, i_start=self.labelRef-1, directed=False)

        # save to self.bridges
        self.bridges = []
        for i in range(1, succs.size):
            n0 = preds[succs[i]] + 1
            n1 = succs[i] + 1
            # read conn
            if n0 > n1:
                nn = [str(n1), str(n0)]
            else:
                nn = [str(n0), str(n1)]

            conn = self.connDict[f'{nn[0]}_{nn[1]}']
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
            bridge['distance'] = ((x1 - x0)**2 + (y1 - y0)**2)**0.5
            self.bridges.append(bridge)
        self.num_bridge = len(self.bridges)
        return self.bridges


    def get_bridge_endpoint_aoi_mask(self, bridge, radius=50):
        # get AOI mask
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
        return aoi_mask0, aoi_mask1


    def unwrap_conn_comp(self, unw, radius=50, ramp_type=None, print_msg=False):
        start_time = time.time()
        radius = int(min(radius, min(self.conncomp.shape)*0.05))

        unw = np.array(unw, dtype=np.float32)
        if self.refY is not None:
            unw[unw != 0.] -= unw[self.refY, self.refX]

        if ramp_type is not None:
            if print_msg:
                print(f'estimate a {ramp_type} ramp')
            ramp_mask = (self.labelImg == self.labelRef)
            unw, ramp = deramp(unw, ramp_mask, ramp_type, metadata=self.metadata)

        for bridge in self.bridges:
            # prepare masks
            aoi_mask0, aoi_mask1 = self.get_bridge_endpoint_aoi_mask(bridge, radius=radius)
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

            # add phase jump
            unw[label_mask1] += 2.* np.pi * num_jump

            if print_msg:
                print(('phase diff {}_{}: {:04.1f} rad --> '
                       'num of jump: {}').format(bridge['label1'],
                                                 bridge['label0'],
                                                 diff_value,
                                                 num_jump))

        # add ramp back
        if ramp_type is not None:
            unw += ramp
        if print_msg:
            print(f'time used: {time.time()-start_time:.2f} secs.')
        return unw


    def plot_bridge(self, ax, cmap='jet', radius=50):
        # label background
        ax.imshow(self.labelImg, cmap=cmap, interpolation='nearest')
        # bridges
        for bridge in self.bridges:
            ax.plot([bridge['x0'], bridge['x1']],
                    [bridge['y0'], bridge['y1']], 'w-', lw=1)
            # endpoint window
            if radius > 0:
                aoi_mask0, aoi_mask1 = self.get_bridge_endpoint_aoi_mask(bridge, radius=radius)
                label_mask0 = self.labelImg == bridge['label0']
                label_mask1 = self.labelImg == bridge['label1']
                mask0 = np.ma.masked_where(~(aoi_mask0*label_mask0), np.zeros(self.labelImg.shape))
                mask1 = np.ma.masked_where(~(aoi_mask1*label_mask1), np.zeros(self.labelImg.shape))
                ax.imshow(mask0, cmap='gray', alpha=0.3, vmin=0, vmax=1)
                ax.imshow(mask1, cmap='gray', alpha=0.3, vmin=0, vmax=1)
        # reference pixel
        ax.plot(self.refX, self.refY, 'ks', ms=2)
        return ax

######################################## end of connectComponent class ####################################
