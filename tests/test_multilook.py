import numpy as np
import pytest

from mintpy.multilook import get_ampcor_box, multilook_file
from mintpy.utils import readfile, writefile


def test_get_ampcor_box():
    # the grid of pycuampcor.dense_offsets: 64 x 64 windows every 32 x 16 pixels, searched +/- 20 x 10
    # pixels, in a file of 1000 x 500 pixels (width x length)
    box = get_ampcor_box(1000, 500, lks_y=16, lks_x=32, search_win=[20, 10], xcorr_win=[64, 64])
    n_x = (1000 - 2*20 - 64) // 32          # 28 windows across
    n_y = (500 - 2*10 - 64) // 16           # 26 windows down
    # the blocks centered on the windows: the first window starts at the search range
    x0 = 20 + 64 // 2 - 32 // 2
    y0 = 10 + 64 // 2 - 16 // 2
    assert box == (x0, y0, x0 + n_x * 32, y0 + n_y * 16)

    # with a margin, and windows as large as the step: the blocks are the windows
    box = get_ampcor_box(1000, 500, lks_y=64, lks_x=64, search_win=[20, 20], xcorr_win=[64, 64], margin=5)
    assert box[:2] == (25, 25)

    # no window fits
    with pytest.raises(ValueError):
        get_ampcor_box(100, 100, lks_y=16, lks_x=16, search_win=[40, 40], xcorr_win=[64, 64])


# even and odd steps and windows: the nearest pixels are the window centers for any parity
@pytest.mark.parametrize('lks_y, lks_x, xcorr', [(16, 32, [64, 64]), (15, 31, [64, 64]), (15, 32, [65, 63])])
def test_multilook_ampcor_nearest(tmp_path, lks_y, lks_x, xcorr):
    # a full resolution geometry with the column and row numbers as values
    length, width = 300, 400
    rows, cols = np.mgrid[0:length, 0:width].astype(np.float32)
    infile = str(tmp_path / 'geometryRadar.h5')
    meta = {'FILE_TYPE': 'geometry', 'LENGTH': length, 'WIDTH': width}
    writefile.write({'latitude': rows, 'longitude': cols}, out_file=infile, metadata=meta)

    # multilooked to the grid of the dense offsets, with the nearest pixels: the window centers
    search = [20, 10]
    outfile = str(tmp_path / 'geometryRadar_mli.h5')
    multilook_file(infile, lks_y, lks_x, outfile=outfile, method='nearest',
                   search_win=search, xcorr_win=xcorr)
    lat = readfile.read(outfile, datasetName='latitude')[0]
    lon = readfile.read(outfile, datasetName='longitude')[0]
    n_y = (length - 2*search[1] - xcorr[1]) // lks_y
    n_x = (width - 2*search[0] - xcorr[0]) // lks_x
    assert lat.shape == lon.shape == (n_y, n_x)
    # the center (pixel) of the window (i, j) of pycuampcor: start + window // 2 + step * index
    assert np.array_equal(lat[:, 0], search[1] + xcorr[1] // 2 + lks_y * np.arange(n_y))
    assert np.array_equal(lon[0, :], search[0] + xcorr[0] // 2 + lks_x * np.arange(n_x))
