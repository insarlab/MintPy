import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("TkAgg")

data = h5py.File('/Users/David/code/PySAR/timeseries.h5', 'r')


def is_error(path1, path2):
    keys = ['bperp', 'date', 'timeseries']

    data1 = h5py.File(path1, 'r')
    data2 = h5py.File(path2, 'r')

    for key in keys:
        if key not in data1.keys():
            raise Exception(path1 + "does not contain the following key: " + key)
        if key not in data2.keys():
            raise Exception(path2 + "does not contain the following key: " + key)
    good = True
    if not bperp_is_good(data1['bperp'], data2['bperp']):
        print("BAD BPERP VALS")
        good = False
    if not dates_are_good(data1['date'], data2['date']):
        print(data1['date'][:])
        print(data2['date'][:])
        print( data2['date'][:] == data2['date'][:])

        print("BAD DATE VALS")
        good = False
    if not timeseries_is_good(data1['timeseries'], data2['timeseries']):
        print("BAD TIMESERIES")
        good = False
    if good:
        print("All is good in the world")
    return good

is_error_func = np.vectorize(lambda x: abs(x) > 1e-5)

def bperp_is_good(bperp1, bperp2):
    if bperp1.shape != bperp2.shape:
        raise Exception("Shapes do not match for bperp: " + str(bperp1.shape) + " != " + str(bperp2.shape))
    diff = bperp1[:] - bperp2[:]
    return not is_error_func(diff).any()

def dates_are_good(date1, date2):
    if date1.shape != date2.shape:
        raise Exception("Date dimensions don't match: " + str(date1.shape) + " != " + str(date2.shape))

    return np.array_equal(date1[:], date2[:])

def timeseries_is_good(t1, t2):
    if t1.shape != t2.shape:
        raise Exception("Date dimensions don't match: " + str(t1.shape) + " != " + str(t2.shape))

    diff = t1[:] - t2[:]
    errors = is_error_func(diff)
    one_d_errors = errors.ravel()

    n, bins, patches = plt.hist(one_d_errors, 20, facecolor='blue', alpha=0.5)
    plt.show()
    print(errors.shape)
    print(np.sum(errors))
    return not errors.any()

if __name__ == "__main__":
    import sys
    print(is_error(sys.argv[1], sys.argv[2]))
