from osgeo import gdal
from pyproj import Proj, Transformer
import argparse
import os
def convert2degree(infile, x, y):
    '''
    convert points in the projection coordinates of the file to lon/lat in degree/
    inputs: infile -- gtiff file
            x1,y1 --scale value, for example x=605750, y=3886250, x, y can be 1/2D numpy array
    pay attention, Transform.from_proj(p1,p2, always_xy=True) make the x,y <-> lon, lat 
    '''
    ds = gdal.Open(infile)

    srs = ds.GetSpatialRef()

    if (not srs.IsProjected()) and (srs.GetAttrValue('unit') == 'degree'):
        return x, y 
    
    p_in = Proj(ds.GetProjection())
    p_out = Proj('epsg:4326')
    transformer = Transformer.from_proj(p_in, p_out, always_xy=True)   
    return transformer.transform(x, y)

def snwe_in_degree(infile, snwe):

    lon, lat = convert2degree(infile, [snwe[2], snwe[3]], [snwe[0], snwe[1]])
    
    snwe_degree=(lat[0], lat[1], lon[0], lon[1] )

    return snwe_degree


def main():
    parser = argparse.ArgumentParser(
    prog=os.path.basename(__file__),
    description=__doc__,
    )

    parser.add_argument(
    "--infile", required=True, 
    help="Geotiff file"   
    )

    parser.add_argument(
        "--x", required=True, 
        help="x-coordinate"
    )

    parser.add_argument(
        "--y", required=True, 
        help="y-coordinate"
    )

    args = parser.parse_args()

    file = args.infile

    x = args.x

    y = args.y

    lon, lat = convert2degree(file, x, y)

    print("lon {}, lat {} in degree:".format(lon, lat))

if __name__ == "__main__":

    main()

