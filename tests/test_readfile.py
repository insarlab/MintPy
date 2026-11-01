import numpy as np
import pytest
from osgeo import gdal, osr

from mintpy.utils import readfile


@pytest.fixture()
def geotiff(tmp_path):
    """Create a small georeferenced raster without external test data."""
    gdal.UseExceptions()
    fname = tmp_path / 'hyp3_unw_phase.tif'
    ds = gdal.GetDriverByName('GTiff').Create(str(fname), 4, 3, 1, gdal.GDT_Float32)
    ds.SetGeoTransform((-120.0, 0.001, 0.0, 35.0, 0.0, -0.001))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    ds.SetProjection(srs.ExportToWkt())
    ds.GetRasterBand(1).WriteArray(np.arange(12, dtype=np.float32).reshape(3, 4))
    ds = None
    return fname


@pytest.mark.parametrize('sidecar_ext', ['.xml', '.aux.xml'])
@pytest.mark.parametrize('with_rsc', [False, True])
def test_geotiff_with_qgis_sidecar(geotiff, sidecar_ext, with_rsc):
    """QGIS PAM metadata must not be parsed as ISCE metadata (#1410)."""
    if with_rsc:
        geotiff.with_name(geotiff.name + '.rsc').write_text(
            'WIDTH 4\nFILE_LENGTH 3\nPROCESSOR hyp3\nDATA_TYPE float32\n'
            'X_FIRST -120.0\nY_FIRST 35.0\nX_STEP 0.001\nY_STEP -0.001\n',
            encoding='utf-8',
        )
    expected = readfile.read_attribute(str(geotiff))
    geotiff.with_name(geotiff.name + sidecar_ext).write_text(
        '<PAMDataset><PAMRasterBand band="1"><Metadata>'
        '<MDI key="STATISTICS_MEAN">5.5</MDI>'
        '</Metadata></PAMRasterBand></PAMDataset>',
        encoding='utf-8',
    )

    assert readfile.read_attribute(str(geotiff)) == expected
    assert int(expected['WIDTH']) == 4
    assert int(expected['LENGTH']) == 3
    assert expected['PROCESSOR'] == ('hyp3' if with_rsc else 'gdal')
    data, _ = readfile.read(str(geotiff))
    np.testing.assert_array_equal(data, np.arange(12, dtype=np.float32).reshape(3, 4))


@pytest.mark.parametrize('xml_type', ['imageFile', 'PAMDataset'])
def test_isce_binary_with_xml(tmp_path, xml_type):
    """ISCE binary rasters must still read their XML metadata."""
    fname = tmp_path / 'hgt.rdr'
    fname.write_bytes(bytes(4 * 3 * 4))
    if xml_type == 'imageFile':
        suffix = '.xml'
        xml = (
            '<imageFile>'
            '<property name="width"><value>4</value></property>'
            '<property name="length"><value>3</value></property>'
            '<property name="data_type"><value>FLOAT</value></property>'
            '<property name="bands"><value>1</value></property>'
            '</imageFile>'
        )
    else:
        suffix = '.aux.xml'
        xml = (
            '<PAMDataset><Metadata domain="ENVI">'
            '<MDI key="samples">4</MDI><MDI key="lines">3</MDI>'
            '<MDI key="data_type">4</MDI><MDI key="bands">1</MDI>'
            '</Metadata></PAMDataset>'
        )
    fname.with_name(fname.name + suffix).write_text(xml, encoding='utf-8')

    metadata = readfile.read_attribute(str(fname))
    assert metadata['PROCESSOR'] == 'isce'
    assert metadata['WIDTH'] == '4'
    assert metadata['LENGTH'] == '3'
    assert metadata['DATA_TYPE'] == 'float32'
