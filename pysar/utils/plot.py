###############################################################################
# Program is part of PySAR v2.0 
# Copyright(c) 2018, Zhang Yunjun 
# Author:  Zhang Yunjun
###############################################################################
# Recommend usage:
#     import pysar.utils.plot as pp

import datetime
import numpy as np
import numpy.matlib as matlib
import scipy.ndimage as ndimage

import matplotlib as mpl
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LightSource
from matplotlib.offsetbox import AnchoredText
from matplotlib.patheffects import withStroke
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap, cm, pyproj

import pysar.utils.datetime as ptime
import pysar.utils.network as pnet
import pysar.utils.utils as ut

minFigSizeSingle = 6.0       # default min size in inch, for single plot
maxFigSizeSingle = 10.0      # default min size in inch, for single plot
defaultFigSizeMulti = [15.0, 8.0]    # default size in inch, for multiple subplots

mplColors = ['#1f77b4',\
             '#ff7f0e',\
             '#2ca02c',\
             '#d62728',\
             '#9467bd',\
             '#8c564b',\
             '#e377c2',\
             '#7f7f7f',\
             '#bcbd22',\
             '#17becf']


############################################ Class Begein ###############################################
class Basemap2(Basemap):
    # add drawscale method to Basemap class. 
    # Basemap.drawmapscale() do not support 'cyl' projection.
    def drawscale(self, lat_c, lon_c, distance, ax=None, font_size=12, yoffset=None, color='k'): 
        """draw a simple map scale from x1,y to x2,y in map projection 
        coordinates, label it with actual distance
        Inputs:
            lat_c/lon_c : float, longitude and latitude of scale bar center, in degree
            distance    : float, distance of scale bar, in m
            yoffset     : float, optional, scale bar length at two ends, in degree
        Example:
            m.drawscale(33.06, 131.18, 2000)
        ref_link: http://matplotlib.1069221.n5.nabble.com/basemap-scalebar-td14133.html
        """
        gc = pyproj.Geod(a=self.rmajor,b=self.rminor) 
        if distance > 1000.0: distance = np.rint(distance/1000.0)*1000.0
        lon_c2, lat_c2, az21 = gc.fwd(lon_c, lat_c, 90, distance)
        length = np.abs(lon_c - lon_c2)
        lon0 = lon_c - length/2.0
        lon1 = lon_c + length/2.0
        if not yoffset:
            yoffset = 0.1*length

        self.plot([lon0,lon1],[lat_c,lat_c],color=color)
        self.plot([lon0,lon0],[lat_c,lat_c+yoffset],color=color)
        self.plot([lon1,lon1],[lat_c,lat_c+yoffset],color=color)
        if not ax:  ax = plt.gca()
        if distance < 1000.0:
            ax.text(lon0+0.5*length, lat_c+yoffset*3, '%d m'%(distance),\
                    verticalalignment='top', horizontalalignment='center',fontsize=font_size, color=color) 
        else:
            ax.text(lon0+0.5*length, lat_c+yoffset*3, '%d km'%(distance/1000.0),\
                    verticalalignment='top', horizontalalignment='center',fontsize=font_size, color=color) 
    
    def auto_lalo_sequence(self, geo_box, lalo_step=None, max_tick_num=4, step_candidate=[1,2,3,4,5]):
        '''Auto calculate lat/lon label sequence based on input geo_box
        Inputs:
            geo_box        : 4-tuple of float, defining UL_lon, UL_lat, LR_lon, LR_lat coordinate
            max_tick_num   : int, rough major tick number along the longer axis
            step_candidate : list of int, candidate list for the significant number of step
        Outputs:
            lats/lons : np.array of float, sequence of lat/lon auto calculated from input geo_box
            lalo_step : float, lat/lon label step
        Example:
            geo_box = (128.0, 37.0, 138.0, 30.0)
            lats, lons, step = m.auto_lalo_sequence(geo_box)
        '''
        max_lalo_dist = max([geo_box[1]-geo_box[3], geo_box[2]-geo_box[0]])

        if not lalo_step:
            # Initial tick step
            lalo_step = ut.round_to_1(max_lalo_dist/max_tick_num)

            # Final tick step - choose from candidate list
            digit = np.int(np.floor(np.log10(lalo_step)))
            lalo_step_candidate = [i*10**digit for i in step_candidate]
            distance = [(i - max_lalo_dist/max_tick_num)**2 for i in lalo_step_candidate]
            lalo_step = lalo_step_candidate[distance.index(min(distance))]
        print('label step - '+str(lalo_step)+' degree')

        # Auto tick sequence
        digit = np.int(np.floor(np.log10(lalo_step)))
        lat_major = np.ceil(geo_box[3]/10**(digit+1))*10**(digit+1)
        lats = np.unique(np.hstack((np.arange(lat_major, lat_major-10.*max_lalo_dist, -lalo_step),\
                                    np.arange(lat_major, lat_major+10.*max_lalo_dist, lalo_step))))
        lats = np.sort(lats[np.where(np.logical_and(lats>=geo_box[3], lats<=geo_box[1]))])

        lon_major = np.ceil(geo_box[0]/10**(digit+1))*10**(digit+1)
        lons = np.unique(np.hstack((np.arange(lon_major, lon_major-10.*max_lalo_dist, -lalo_step),\
                                    np.arange(lon_major, lon_major+10.*max_lalo_dist, lalo_step))))
        lons = np.sort(lons[np.where(np.logical_and(lons>=geo_box[0], lons<=geo_box[2]))])
 
        return lats, lons, lalo_step


    def draw_lalo_label(self, geo_box, ax=None, lalo_step=None, labels=[1,0,0,1], font_size=12, color='k'):
        '''Auto draw lat/lon label/tick based on coverage from geo_box
        Inputs:
            geo_box : 4-tuple of float, defining UL_lon, UL_lat, LR_lon, LR_lat coordinate
            labels  : list of 4 int, positions where the labels are drawn as in [left, right, top, bottom]
                      default: [1,0,0,1]
            ax      : axes object the labels are drawn
            draw    : bool, do not draw if False
        Outputs:
            
        Example:
            geo_box = (128.0, 37.0, 138.0, 30.0)
            m.draw_lalo_label(geo_box)
        '''
        lats, lons, step = self.auto_lalo_sequence(geo_box, lalo_step=lalo_step)

        digit = np.int(np.floor(np.log10(step)))
        fmt = '%.'+'%d'%(abs(min(digit, 0)))+'f'
        # Change the 2 lines below for customized label
        #lats = np.linspace(31.55, 31.60, 2)
        #lons = np.linspace(130.60, 130.70, 3)

        # Plot x/y tick without label
        if not ax:
            ax = plt.gca()
        ax.tick_params(which='both', direction='in', labelsize=font_size, bottom=True, top=True, left=True, right=True)

        ax.set_xticks(lons)
        ax.set_yticks(lats)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        #ax.xaxis.tick_top()
        
        # Plot x/y label
        labels_lat = np.multiply(labels, [1,1,0,0])
        labels_lon = np.multiply(labels, [0,0,1,1])
        self.drawparallels(lats, fmt=fmt, labels=labels_lat, linewidth=0.05, fontsize=font_size, color=color, textcolor=color)
        self.drawmeridians(lons, fmt=fmt, labels=labels_lon, linewidth=0.05, fontsize=font_size, color=color, textcolor=color)



############################################ Plot Utilities #############################################
def add_inner_title(ax, title, loc, size=None, **kwargs):
    if size is None:
        size = dict(size=plt.rcParams['legend.fontsize'])
    at = AnchoredText(title, loc=loc, prop=size,
                      pad=0., borderpad=0.5,
                      frameon=False, **kwargs)
    ax.add_artist(at)
    at.txt._text.set_path_effects([withStroke(foreground="w", linewidth=3)])
    return at


def auto_flip_direction(atr_dict):
    '''Check flip left-right and up-down based on attribute dict, for radar-coded file only'''
    # default value
    flip_lr = False
    flip_ud = False

    # check auto flip only for file in radar coord
    try:
        atr_dict['X_FIRST']
        geocoord = True
    except:
        geocoord = False
        try:
            print(atr_dict['ORBIT_DIRECTION']+' orbit')
            if atr_dict['ORBIT_DIRECTION'][0].upper() == 'A':
                flip_ud = True
            else:
                flip_lr = True
        except: pass
    return flip_lr, flip_ud


def auto_row_col_num(subplot_num, data_shape, fig_size, fig_num=1):
    '''Get optimal row and column number given figure size number of subplots
    
    Inputs:
        subplot_num : int, total number of subplots
        data_shape  : list of 2 float, data size in pixel in row and column direction of each plot
        fig_size    : list of 2 float, figure window size in inches
        fig_num     : int, number of figure windows, optional, default = 1.

    Outputs:
        row_num : number of subplots in row    direction per figure
        col_num : number of subplots in column direction per figure
    '''
    subplot_num_per_fig = int(np.ceil(float(subplot_num)/float(fig_num)))

    data_shape_ratio = float(data_shape[0])/float(data_shape[1])
    num_ratio = fig_size[1]/fig_size[0]/data_shape_ratio
    row_num = np.sqrt(subplot_num_per_fig*num_ratio)
    col_num = np.sqrt(subplot_num_per_fig/num_ratio)
    while np.rint(row_num)*np.rint(col_num) < subplot_num_per_fig:
        if row_num%1 > col_num%1:
            row_num += 0.5
        else:
            col_num += 0.5
    row_num = int(np.rint(row_num))
    col_num = int(np.rint(col_num))
    
    return row_num, col_num


def check_colormap_input(atr_dict, colormap=None, datasetName=None):
    if not colormap:
        if (atr_dict['FILE_TYPE'] in ['coherence','temporal_coherence','.cor','.mli','.slc','.amp','.ramp']\
            or (atr_dict['FILE_TYPE'] == 'ifgramStack' and datasetName.split('-')[0] in ['coherence', 'connectComponent'])):
              colormap = 'gray'
        else: colormap = 'jet'
    print('colormap: '+colormap)

    # Modified hsv colormap by H. Fattahi
    if colormap == 'hsv':
        cdict1 = {'red':   ((0.0, 0.0, 0.0),
                            (0.5, 0.0, 0.0),
                            (0.6, 1.0, 1.0),
                            (0.8, 1.0, 1.0),
                            (1.0, 0.5, 0.5)),
                  'green': ((0.0, 0.0, 0.0),
                            (0.2, 0.0, 0.0),
                            (0.4, 1.0, 1.0),
                            (0.6, 1.0, 1.0),
                            (0.8, 0.0, 0.0),
                            (1.0, 0.0, 0.0)),
                   'blue':  ((0.0, 0.5, .5),
                            (0.2, 1.0, 1.0),
                            (0.4, 1.0, 1.0),
                            (0.5, 0.0, 0.0),
                            (1.0, 0.0, 0.0),)
                 }
        colormap = LinearSegmentedColormap('BlueRed1', cdict1)
    else:
        colormap = plt.get_cmap(colormap)
    return colormap


def auto_adjust_xaxis_date(ax, datevector, fontSize=12, every_year=1):
    '''Adjust X axis
    Input:
        ax : matplotlib figure axes object
        datevector : list of float, date in years
                     i.e. [2007.013698630137, 2007.521917808219, 2007.6463470319634]
    Output:
        ax  - matplotlib figure axes object
        dss - datetime.date object, xmin
        dee - datetime.date object, xmax
    '''

    # Min/Max
    ts=datevector[0] -0.2;  ys=int(ts);  ms=int((ts-ys)*12.0)
    te=datevector[-1]+0.3;  ye=int(te);  me=int((te-ye)*12.0)
    if ms>12:   ys = ys+1;   ms=1
    if me>12:   ye = ye+1;   me=1
    if ms<1:    ys = ys-1;   ms=12
    if me<1:    ye = ye-1;   me=12
    dss=datetime.date(ys,ms,1)
    dee=datetime.date(ye,me,1)
    ax.set_xlim(dss,dee)

    # Label/Tick format
    ax.fmt_xdata = mdates.DateFormatter('%Y-%m-%d %H:%M:%S')
    ax.xaxis.set_major_locator(mdates.YearLocator(every_year))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax.xaxis.set_minor_locator(mdates.MonthLocator())

    # Label font size
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontSize)
    #fig2.autofmt_xdate()     #adjust x overlap by rorating, may enble again
    return ax, dss, dee


def auto_adjust_yaxis(ax, dataList, fontSize=12, ymin=None, ymax=None):
    '''Adjust Y axis
    Input:
        ax       : matplot figure axes object
        dataList : list of float, value in y axis
        fontSize : float, font size
        ymin     : float, lower y axis limit
        ymax     : float, upper y axis limit
    Output:
        ax
    '''
    # Min/Max
    dataRange = max(dataList) - min(dataList)
    if ymin is None:  ymin = min(dataList) - 0.1*dataRange
    if ymax is None:  ymax = max(dataList) + 0.1*dataRange
    ax.set_ylim([ymin, ymax])
    ## Tick/Label setting
    #xticklabels = plt.getp(ax, 'xticklabels')
    #yticklabels = plt.getp(ax, 'yticklabels')
    #plt.setp(yticklabels, 'color', 'k', fontsize=fontSize)
    #plt.setp(xticklabels, 'color', 'k', fontsize=fontSize)

    return ax



####################################### Plot ################################################
def plot_coherence_history(ax, date12_list, coherence_list, plot_dict={}):
    '''Plot min/max Coherence of all interferograms for each date'''
    # Figure Setting
    if not 'fontsize'    in plot_dict.keys():   plot_dict['fontsize']    = 12
    if not 'linewidth'   in plot_dict.keys():   plot_dict['linewidth']   = 2
    if not 'markercolor' in plot_dict.keys():   plot_dict['markercolor'] = 'orange'
    if not 'markersize'  in plot_dict.keys():   plot_dict['markersize']  = 16
    if not 'disp_title'  in plot_dict.keys():   plot_dict['disp_title']  = True
    if not 'every_year'  in plot_dict.keys():   plot_dict['every_year']  = 1

    # Get date list
    m_dates = [date12.split('-')[0] for date12 in date12_list]
    s_dates = [date12.split('-')[1] for date12 in date12_list]
    date8_list = sorted(ptime.yyyymmdd(list(set(m_dates + s_dates))))

    dates, datevector = ptime.date_list2vector(date8_list)
    bar_width = ut.most_common(np.diff(dates).tolist())*3/4
    x_list = [i-bar_width/2 for i in dates]

    coh_mat = pnet.coherence_matrix(date12_list, coherence_list)

    ax.bar(x_list, np.nanmax(coh_mat, axis=0), bar_width.days, label='Max Coherence')
    ax.bar(x_list, np.nanmin(coh_mat, axis=0), bar_width.days, label='Min Coherence')

    if plot_dict['disp_title']:
        ax.set_title('Coherence History of All Related Interferograms')

    ax = auto_adjust_xaxis_date(ax, datevector, plot_dict['fontsize'], every_year=plot_dict['every_year'])[0]
    ax.set_ylim([0.0,1.0])

    ax.set_xlabel('Time [years]',fontsize=plot_dict['fontsize'])
    ax.set_ylabel('Coherence',fontsize=plot_dict['fontsize'])
    ax.legend(loc='lower right')

    return ax


def plot_network(ax, date12_list, date_list, pbase_list, plot_dict={}, date12_list_drop=[], printMsg=True):
    '''Plot Temporal-Perp baseline Network
    Inputs
        ax : matplotlib axes object
        date12_list : list of string for date12 in YYMMDD-YYMMDD format
        date_list   : list of string, for date in YYYYMMDD/YYMMDD format
        pbase_list  : list of float, perp baseline, len=number of acquisition
        plot_dict   : dictionary with the following items:
                      fontsize
                      linewidth
                      markercolor
                      markersize

                      coherence_list : list of float, coherence value of each interferogram, len = number of ifgrams
                      disp_min/max :  float, min/max range of the color display based on coherence_list
                      colormap : string, colormap name
                      coh_thres : float, coherence of where to cut the colormap for display
                      disp_title : bool, show figure title or not, default: True
                      disp_drop: bool, show dropped interferograms or not, default: True
    Output
        ax : matplotlib axes object
    '''
    
    # Figure Setting
    if not 'fontsize'    in plot_dict.keys():   plot_dict['fontsize']    = 12
    if not 'linewidth'   in plot_dict.keys():   plot_dict['linewidth']   = 2
    if not 'markercolor' in plot_dict.keys():   plot_dict['markercolor'] = 'orange'
    if not 'markersize'  in plot_dict.keys():   plot_dict['markersize']  = 16
    # For colorful display of coherence
    if not 'coherence_list' in plot_dict.keys():  plot_dict['coherence_list'] = None
    if not 'cbar_label'     in plot_dict.keys():  plot_dict['cbar_label']     = 'Coherence'
    if not 'disp_min'       in plot_dict.keys():  plot_dict['disp_min']       = 0.2
    if not 'disp_max'       in plot_dict.keys():  plot_dict['disp_max']       = 1.0
    if not 'colormap'       in plot_dict.keys():  plot_dict['colormap']       = 'RdBu'
    if not 'disp_title'     in plot_dict.keys():  plot_dict['disp_title']     = True
    if not 'coh_thres'      in plot_dict.keys():  plot_dict['coh_thres']      = None
    if not 'disp_drop'      in plot_dict.keys():  plot_dict['disp_drop']      = True
    if not 'every_year'     in plot_dict.keys():  plot_dict['every_year']     = 1
    coh_list = plot_dict['coherence_list']
    disp_min = plot_dict['disp_min']
    disp_max = plot_dict['disp_max']
    coh_thres = plot_dict['coh_thres']
    transparency = 0.7

    # Date Convert
    date8_list = ptime.yyyymmdd(sorted(date_list))
    date6_list = ptime.yymmdd(date8_list)
    dates, datevector = ptime.date_list2vector(date8_list)
    tbase_list = ptime.date_list2tbase(date8_list)[0]

    ## maxBperp and maxBtemp
    ifgram_num = len(date12_list)
    pbase12 = np.zeros(ifgram_num)
    tbase12 = np.zeros(ifgram_num)
    for i in range(ifgram_num):
        m_date, s_date = date12_list[i].split('-')
        m_idx = date6_list.index(m_date)
        s_idx = date6_list.index(s_date)
        pbase12[i] = pbase_list[s_idx] - pbase_list[m_idx]
        tbase12[i] = tbase_list[s_idx] - tbase_list[m_idx]
    print(('max perpendicular baseline: %.2f m' % (np.max(np.abs(pbase12)))))
    print(('max temporal      baseline: %d days' % (np.max(tbase12))))

    ## Keep/Drop - date12
    date12_list_keep = sorted(list(set(date12_list) - set(date12_list_drop)))
    idx_date12_keep = [date12_list.index(i) for i in date12_list_keep]
    idx_date12_drop = [date12_list.index(i) for i in date12_list_drop]
    if not date12_list_drop:
        plot_dict['disp_drop'] = False

    ## Keep/Drop - date
    m_dates = [i.split('-')[0] for i in date12_list_keep]
    s_dates = [i.split('-')[1] for i in date12_list_keep]
    date8_list_keep = ptime.yyyymmdd(sorted(list(set(m_dates + s_dates))))
    date8_list_drop = sorted(list(set(date8_list) - set(date8_list_keep)))
    idx_date_keep = [date8_list.index(i) for i in date8_list_keep]
    idx_date_drop = [date8_list.index(i) for i in date8_list_drop]


    # Ploting
    #ax=fig.add_subplot(111)
    ## Colorbar when conherence is colored
    if coh_list is not None:
        data_min = min(coh_list)
        data_max = max(coh_list)
        # Normalize
        normalization = False
        if normalization:
            coh_list = [(coh-data_min) / (data_min-data_min) for coh in coh_list]
            disp_min = data_min
            disp_max = data_max

        if printMsg:
            print('showing coherence')
            print(('colormap: '+plot_dict['colormap']))
            print(('display range: '+str([disp_min, disp_max])))
            print(('data    range: '+str([data_min, data_max])))

        splitColormap = True
        if splitColormap:
            # Use lower/upper part of colormap to emphasis dropped interferograms
            if not coh_thres:
                # Find proper cut percentage so that all keep pairs are blue and drop pairs are red
                coh_list_keep = [coh_list[i] for i in idx_date12_keep]
                coh_list_drop = [coh_list[i] for i in idx_date12_drop]
                if coh_list_drop:
                    coh_thres = max(coh_list_drop)
                else:
                    coh_thres = min(coh_list_keep)
            if coh_thres < disp_min:
                disp_min = 0.0
                if printMsg:
                    print(('data range exceed orginal display range, set new display range to: [0.0, %f]' % (disp_max)))
            c1_num = np.ceil(200.0 * (coh_thres - disp_min) / (disp_max - disp_min)).astype('int')
            coh_thres = c1_num / 200.0 * (disp_max-disp_min) + disp_min
            cmap = plt.get_cmap(plot_dict['colormap'])
            colors1 = cmap(np.linspace(0.0, 0.3, c1_num))
            colors2 = cmap(np.linspace(0.6, 1.0, 200 - c1_num))
            cmap = LinearSegmentedColormap.from_list('truncate_RdBu', np.vstack((colors1, colors2)))
            if printMsg:
                print(('color jump at '+str(coh_thres)))
        else:
            cmap = plt.get_cmap(plot_dict['colormap'])

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", "3%", pad="3%")
        norm = mpl.colors.Normalize(vmin=disp_min, vmax=disp_max)
        cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
        cbar.set_label(plot_dict['cbar_label'], fontsize=plot_dict['fontsize'])

        #plot low coherent ifgram first and high coherence ifgram later
        coh_list_keep = [coh_list[date12_list.index(i)] for i in date12_list_keep]
        date12_list_keep = [x for _,x in sorted(zip(coh_list_keep, date12_list_keep))]

    ## Dot - SAR Acquisition
    if idx_date_keep:
        x_list = [dates[i] for i in idx_date_keep]
        y_list = [pbase_list[i] for i in idx_date_keep]
        ax.plot(x_list, y_list, 'ko', alpha=0.7, ms=plot_dict['markersize'], mfc=plot_dict['markercolor'])
    if idx_date_drop:
        x_list = [dates[i] for i in idx_date_drop]
        y_list = [pbase_list[i] for i in idx_date_drop]
        ax.plot(x_list, y_list, 'ko', alpha=0.7, ms=plot_dict['markersize'], mfc='gray')

    ## Line - Pair/Interferogram        
    # interferograms dropped
    if plot_dict['disp_drop']:
        for date12 in date12_list_drop:
            date1, date2 = date12.split('-')
            idx1 = date6_list.index(date1)
            idx2 = date6_list.index(date2)
            x = np.array([dates[idx1], dates[idx2]])
            y = np.array([pbase_list[idx1], pbase_list[idx2]])
            if coh_list:
                coh = coh_list[date12_list.index(date12)]
                coh_idx = (coh - disp_min) / (disp_max - disp_min)
                ax.plot(x, y, '--', lw=plot_dict['linewidth'], alpha=transparency, c=cmap(coh_idx)) 
            else:
                ax.plot(x, y, '--', lw=plot_dict['linewidth'], alpha=transparency, c='k')

    # interferograms kept
    for date12 in date12_list_keep:
        date1, date2 = date12.split('-')
        idx1 = date6_list.index(date1)
        idx2 = date6_list.index(date2)
        x = np.array([dates[idx1], dates[idx2]])
        y = np.array([pbase_list[idx1], pbase_list[idx2]])
        if coh_list is not None:
            coh = coh_list[date12_list.index(date12)]
            coh_idx = (coh - disp_min) / (disp_max - disp_min)
            ax.plot(x, y, '-', lw=plot_dict['linewidth'], alpha=transparency, c=cmap(coh_idx)) 
        else:
            ax.plot(x, y, '-', lw=plot_dict['linewidth'], alpha=transparency, c='k')

    if plot_dict['disp_title']:
        ax.set_title('Interferogram Network', fontsize=plot_dict['fontsize'])

    # axis format
    ax = auto_adjust_xaxis_date(ax, datevector, plot_dict['fontsize'], every_year=plot_dict['every_year'])[0]
    ax = auto_adjust_yaxis(ax, pbase_list, plot_dict['fontsize'])
    ax.set_xlabel('Time [years]',fontsize=plot_dict['fontsize'])
    ax.set_ylabel('Perp Baseline [m]',fontsize=plot_dict['fontsize'])

    # Legend
    if plot_dict['disp_drop']:
        solid_line = mlines.Line2D([],[],color='k',ls='solid', label='Interferograms')
        dash_line  = mlines.Line2D([],[],color='k',ls='dashed', label='Interferograms dropped')
        ax.legend(handles=[solid_line,dash_line])

    return ax


def plot_perp_baseline_hist(ax, date8_list, pbase_list, plot_dict={}, date8_list_drop=[]):
    ''' Plot Perpendicular Spatial Baseline History
    Inputs
        ax : matplotlib axes object
        date8_list : list of string, date in YYYYMMDD format
        pbase_list : list of float, perp baseline 
        plot_dict : dictionary with the following items:
                    fontsize
                    linewidth
                    markercolor
                    markersize
                    disp_title : bool, show figure title or not, default: True
                    every_year : int, number of years for the major tick on xaxis
        date8_list_drop : list of string, date dropped in YYYYMMDD format
                          e.g. ['20080711', '20081011']
    Output:
        ax : matplotlib axes object
    '''
    # Figure Setting
    if not 'fontsize'    in plot_dict.keys():   plot_dict['fontsize']    = 12
    if not 'linewidth'   in plot_dict.keys():   plot_dict['linewidth']   = 2
    if not 'markercolor' in plot_dict.keys():   plot_dict['markercolor'] = 'orange'
    if not 'markersize'  in plot_dict.keys():   plot_dict['markersize']  = 16
    if not 'disp_title'  in plot_dict.keys():   plot_dict['disp_title']  = True
    if not 'every_year'  in plot_dict.keys():   plot_dict['every_year']  = 1
    transparency = 0.7

    # Date Convert
    dates, datevector = ptime.date_list2vector(date8_list)

    # Get index of date used and dropped
    #date8_list_drop = ['20080711', '20081011']  # for debug
    idx_keep = list(range(len(date8_list)))
    idx_drop = []
    for i in date8_list_drop:
        idx = date8_list.index(i)
        idx_keep.remove(idx)
        idx_drop.append(idx)

    # Plot
    #ax=fig.add_subplot(111)

    # Plot date used
    if idx_keep:
        x_list = [dates[i] for i in idx_keep]
        y_list = [pbase_list[i] for i in idx_keep]
        ax.plot(x_list, y_list, '-ko', alpha=transparency, lw=plot_dict['linewidth'], \
                ms=plot_dict['markersize'], mfc=plot_dict['markercolor'])
    
    # Plot date dropped
    if idx_drop:
        x_list = [dates[i] for i in idx_drop]
        y_list = [pbase_list[i] for i in idx_drop]
        ax.plot(x_list, y_list, 'ko', alpha=transparency, ms=plot_dict['markersize'], mfc='gray')

    if plot_dict['disp_title']:
        ax.set_title('Perpendicular Baseline History',fontsize=plot_dict['fontsize'])

    # axis format
    ax = auto_adjust_xaxis_date(ax, datevector, plot_dict['fontsize'], every_year=plot_dict['every_year'])[0]
    ax = auto_adjust_yaxis(ax, pbase_list, plot_dict['fontsize'])
    ax.set_xlabel('Time [years]',fontsize=plot_dict['fontsize'])
    ax.set_ylabel('Perpendicular Baseline [m]',fontsize=plot_dict['fontsize'])

    return ax


def plot_coherence_matrix(ax, date12_list, coherence_list, date12_list_drop=[], plot_dict={}):
    '''Plot Coherence Matrix of input network
    
    if date12_list_drop is not empty, plot KEPT pairs in the upper triangle and
                                           ALL  pairs in the lower triangle.
    '''
    # Figure Setting
    if not 'fontsize'    in plot_dict.keys():   plot_dict['fontsize']    = 12
    if not 'linewidth'   in plot_dict.keys():   plot_dict['linewidth']   = 2
    if not 'markercolor' in plot_dict.keys():   plot_dict['markercolor'] = 'orange'
    if not 'markersize'  in plot_dict.keys():   plot_dict['markersize']  = 16
    if not 'disp_title'  in plot_dict.keys():   plot_dict['disp_title']  = True
    if not 'cbar_label'  in plot_dict.keys():   plot_dict['cbar_label']  = 'Coherence'

    coh_mat = pnet.coherence_matrix(date12_list, coherence_list)

    if date12_list_drop:
        # Date Convert
        m_dates = [i.split('-')[0] for i in date12_list]
        s_dates = [i.split('-')[1] for i in date12_list]
        date6_list = ptime.yymmdd(sorted(list(set(m_dates + s_dates))))
        # Set dropped pairs' value to nan, in upper triangle only.
        for date12 in date12_list_drop:
            idx1, idx2 = [date6_list.index(i) for i in date12.split('-')]
            coh_mat[idx1, idx2] = np.nan

    #Show diagonal value as black, to be distinguished from un-selected interferograms
    diag_mat = np.diag(np.ones(coh_mat.shape[0]))
    diag_mat[diag_mat == 0.] = np.nan
    im = ax.imshow(diag_mat, cmap='gray_r', vmin=0.0, vmax=1.0, interpolation='nearest')
    #Show coherence matrix
    im = ax.imshow(coh_mat, cmap='jet', vmin=0.0, vmax=1.0, interpolation='nearest')

    date_num = coh_mat.shape[0]
    if date_num < 30:
        tick_list = list(range(0,date_num,5))
    else:
        tick_list = list(range(0,date_num,10))
    ax.get_xaxis().set_ticks(tick_list)
    ax.get_yaxis().set_ticks(tick_list)
    ax.set_xlabel('Image Number', fontsize=plot_dict['fontsize'])
    ax.set_ylabel('Image Number', fontsize=plot_dict['fontsize'])

    if plot_dict['disp_title']:
        ax.set_title('Coherence Matrix')

    # Colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "3%", pad="3%")
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label(plot_dict['cbar_label'], fontsize=plot_dict['fontsize'])

    # Legend
    if date12_list_drop:
        ax.plot([],[],label='Upper: used ifgrams')
        ax.plot([],[],label='Lower: all ifgrams')
        ax.legend(handlelength=0)

    return ax


def plot_dem_lalo(bmap, dem, box, inps_dict):
    '''Plot DEM in geo-coordinate
    Inputs:
        bmap  : basemap object
        dem   : dem data, 2D np.int16 matrix
        box   : geo bounding box, 4-tuple as (urcrnrlon,urcrnrlat,llcrnrlon,llcrnrlat)
        inps_dict : dict with the following 5 items:
                    'disp_dem_shade'    : bool,  True/False
                    'disp_dem_contour'  : bool,  True/False
                    'dem_contour_step'  : float, 200.0
                    'dem_contour_smooth': float, 3.0
    
    Examples:
        dem_disp_dict = {'dem': 'gsi10m_30m.dem', 'disp_dem_shade': True, 'disp_dem_contour': True,\
                         'dem_contour_step': 200.0, 'dem_contour_smooth': 3.0}
        bmap = plot_dem_lalo(bmap,dem,geo_box,dem_inps_dict)
    '''

    if inps_dict['disp_dem_shade']:
        print('show shaded relief DEM')
        ls = LightSource(azdeg=315, altdeg=45)
        dem_shade = ls.shade(dem, vert_exag=1.0, cmap=plt.cm.gray, vmin=-5000, vmax=np.nanmax(dem)+2000)

        ### mask out water
        #dem_shade = ls.shade(dem, vert_exag=1.0, cmap=plt.cm.gray, vmin=-5000, vmax=np.nanmax(dem)+500)
        #mask_file = '/Users/jeromezhang/Documents/insarlab/Kyushu/Velocity/mask_land.h5'
        #mask_mat = readfile.read(mask_file, datasetName='mask', box=inps_dict['dem_pix_box'])[0]
        #dem_shade = mask.mask_matrix(dem_shade, mask_mat)

        bmap.imshow(dem_shade, origin='upper', interpolation='spline16')

    if inps_dict['disp_dem_contour']:
        print('show contour in step - '+str(inps_dict['dem_contour_step'])+' m'+\
              ' with smoothing factor - '+str(inps_dict['dem_contour_smooth']))
        c_x = np.linspace(box[0], box[2], num=dem.shape[1], endpoint='FALSE').reshape(1,dem.shape[1])
        c_y = np.linspace(box[1], box[3], num=dem.shape[0], endpoint='FALSE').reshape(dem.shape[0],1)
        c_xx = matlib.repmat(c_x, dem.shape[0], 1)
        c_yy = matlib.repmat(c_y, 1, dem.shape[1])

        dem_contour = ndimage.gaussian_filter(dem, sigma=inps_dict['dem_contour_smooth'], order=0)
        contour_sequence = np.arange(-6000, 9000, inps_dict['dem_contour_step'])

        bmap.contour(c_xx, c_yy, dem_contour, contour_sequence, origin='upper',colors='black',alpha=0.5, latlon='FALSE')

    return bmap


def plot_dem_yx(ax, dem, inps_dict=dict()):
    '''Plot DEM in radar coordinate
    Inputs:
        ax         : matplotlib axes object
        dem        : dem data, 2D np.int16 matrix
        inps_dict : dict with the following 5 items:
                    'disp_dem_shade'    : bool,  True/False
                    'disp_dem_contour'  : bool,  True/False
                    'dem_contour_step'  : float, 200.0
                    'dem_contour_smooth': float, 3.0

    Examples:
        dem_disp_dict = {'dem': 'gsi10m_30m.dem', 'disp_dem_shade': True, 'disp_dem_contour': True,\
                         'dem_contour_step': 200.0, 'dem_contour_smooth': 3.0}
        ax = plot_dem_yx(ax,dem,dem_disp_dict)
    '''
    if not inps_dict:
        inps_dict['disp_dem_shade']     = True
        inps_dict['disp_dem_contour']   = True
        inps_dict['dem_contour_smooth'] = 3.0
        inps_dict['dem_contour_step']   = 200.0

    if inps_dict['disp_dem_shade']:
        print('show shaded relief DEM')
        ls = LightSource(azdeg=315, altdeg=45)
        dem_shade = ls.shade(dem, vert_exag=1.0, cmap=plt.cm.gray, vmin=-5000, vmax=np.nanmax(dem)+2000)
        ax.imshow(dem_shade, interpolation='spline16')

    if inps_dict['disp_dem_contour']:
        print('show contour in step: '+str(inps_dict['dem_contour_step'])+' m'+\
              ' with smoothing factor - '+str(inps_dict['dem_contour_smooth']))
        dem_contour = ndimage.gaussian_filter(dem, sigma=inps_dict['dem_contour_smooth'], order=0)
        contour_sequence = np.arange(-6000, 9000, inps_dict['dem_contour_step'])
        ax.contour(dem_contour, contour_sequence, origin='lower',colors='black',alpha=0.5)

    return ax




