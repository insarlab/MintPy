############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Dec 2020                           #
############################################################
# Recommend import:
#    from mintpy.objects.progress import progressBar
# OR from mintpy.utils import ptime


import sys
import time
import numpy as np


###########################Simple progress bar######################
class progressBar:
    """Creates a text-based progress bar. Call the object with
    the simple print command to see the progress bar, which looks
    something like this:
    [=======> 22%       ]
    You may specify the progress bar's min and max values on init.

    Code modified from PyAPS release 1.0 (http://earthdef.caltech.edu/projects/pyaps/wiki/Main)
    Code originally from http://code.activestate.com/recipes/168639/

    example:
        from mintpy.objects import timeseries
        from mintpy.utils import ptime
        date_list = timeseries('timeseries.h5').get_date_list()
        num_date = len(date_list)
        prog_bar = ptime.progressBar(maxValue=num_date)
        for i in range(num_date):
            do_work()
            prog_bar.update(i+1, suffix=date_list[i])
        prog_bar.close()
    """

    def __init__(self, maxValue=100, prefix='', minValue=0, totalWidth=70, print_msg=True):
        self.prog_bar = "[]"  # This holds the progress bar string
        self.min = minValue
        self.max = maxValue
        self.span = maxValue - minValue
        self.suffix = ''
        self.prefix = prefix
        self.print_msg = print_msg

        ## calculate total width based on console width
        #rows, columns = os.popen('stty size', 'r').read().split()
        #self.width = round(int(columns) * 0.7 / 10) * 10
        self.width = totalWidth
        self.reset()

    def reset(self):
        self.start_time = time.time()
        self.amount = 0        # When amount == max, we are 100% done
        self.update_amount(0)  # Build progress bar string

    def update_amount(self, newAmount=0, suffix=''):
        """ Update the progress bar with the new amount (with min and max
        values set at initialization; if it is over or under, it takes the
        min or max value as a default. """
        newAmount = max(newAmount, self.min)
        newAmount = min(newAmount, self.max)
        self.amount = newAmount

        # Add whitespace to suffix if turned ON
        if suffix:
            suffix = ' '+suffix

        # Figure out the new percent done (round to an integer)
        diffFromMin = np.float(self.amount - self.min)
        percentDone = (diffFromMin / np.float(self.span)) * 100.0
        percentDone = np.int(np.round(percentDone))

        # Figure out how many hash bars the percentage should be
        allFull = self.width - 2 - 18
        numHashes = (percentDone / 100.0) * allFull
        numHashes = np.int(np.round(numHashes))

        # Build a progress bar with an arrow of equal signs; special cases for empty and full
        if numHashes == 0:
            self.prog_bar = '{}[>{}]'.format(self.prefix, ' '*(allFull-1))

        elif numHashes == allFull:
            self.prog_bar = '{}[{}]'.format(self.prefix, '='*allFull)
            self.prog_bar += suffix

        else:
            self.prog_bar = '[{}>{}]'.format('='*(numHashes-1), ' '*(allFull-numHashes))

            # figure out where to put the percentage (roughly centered)
            percentPlace = int(len(self.prog_bar)/2 - len(str(percentDone)))
            percentString = ' {}% '.format(percentDone)

            # slice the percentage into the bar
            self.prog_bar = ''.join([self.prog_bar[0:percentPlace],
                                     percentString,
                                     self.prog_bar[percentPlace+len(percentString):]])

            # prefix and suffix
            self.prog_bar = self.prefix + self.prog_bar + suffix

            # time info - elapsed time and estimated remaining time
            if percentDone > 0:
                elapse_time = time.time() - self.start_time
                remain_time = int(elapse_time * (100./percentDone-1))
                self.prog_bar += '{:5d}s / {:5d}s'.format(int(elapse_time), int(remain_time))

    def update(self, value, every=1, suffix=''):
        """ Updates the amount, and writes to stdout. Prints a
         carriage return first, so it will overwrite the current
          line in stdout."""
        if value % every == 0 or value >= self.max:
            self.update_amount(newAmount=value, suffix=suffix)
            if self.print_msg:
                sys.stdout.write('\r' + self.prog_bar)
                sys.stdout.flush()

    def close(self):
        """Prints a blank space at the end to ensure proper printing
        of future statements."""
        if self.print_msg:
            print(' ')
################################End of progress bar class####################################


