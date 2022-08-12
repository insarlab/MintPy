############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2016                               #
############################################################


import matplotlib.pyplot as plt
from mintpy.utils import readfile, ptime, utils as ut, plot as pp


#############################  Main Function  ################################
def spatial_average(inps):
    print('\n*************** Spatial Average ******************')
    mean_list, date_list = ut.spatial_average(
        inps.file,
        datasetName=inps.datasetName,
        maskFile=inps.mask_file,
        saveList=True,
    )

    atr = readfile.read_attribute(inps.file)
    if inps.disp_fig and atr['FILE_TYPE'] == 'timeseries':
        dates = ptime.date_list2vector(date_list)[0]
        # plot
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(dates, mean_list, '-o')#, lw=2, ms=16, alpha=0.7) #, mfc='crimson')
        ax.set_title('Spatial Average', fontsize=12)
        ax = pp.auto_adjust_xaxis_date(ax, dates)[0]
        ax.set_xlabel('Time [years]', fontsize=12)
        ax.set_ylabel('Mean', fontsize=12)
        plt.show()
