from osgeo import gdal
from pyproj import Proj, Transformer
import argparse
import os
def convert2degree(infile, points):
    '''
    convert points in the projection coordinates of the file to lon/lat in degree/
    inputs: infile -- gtiff file
            points -- 'x1,y1', 'x2,y2',... for example points '3886250, 605750'
    pay attention, lat, lon = transformer.transform(x,y)
    '''
    ds = gdal.Open(infile)

    srs = ds.GetSpatialRef()

    if (not srs.IsProjected()) and (srs.GetAttrValue('unit') == 'degree'):
        return points
    
    p_in = Proj(ds.GetProjection())
    p_out = Proj('epsg:4326')
    transformer = Transformer.from_proj(p_in, p_out)
    points_out = []
    for point in points:   
        lat, lon = transformer.transform(point[0], point[1])  
        points_out.append((lon,lat)) 
    del ds
    return points_out

def snwe_in_degree(infile, snwe):
    points = [(snwe[2], snwe[0]), (snwe[3], snwe[1])]
    points_out = convert2degree(infile, points)   
    snwe_degree=(points_out[0][1], points_out[1][1], points_out[0][0], points_out[1][0] )
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
        "--points", required=True, nargs ="*", 
        help="'x1,y1','x2,y2',..."
    )

    args = parser.parse_args()

    file = args.infile

    points = args.points

    points_list =[]

    for point in points:

        plist = point.split(",")
        ptuple = tuple(map(float, plist))
        points_list.append(ptuple)



    points_out = convert2degree(file, points_list)

    print("points in degree: {}".format(points_out))

if __name__ == "__main__":

    main()

