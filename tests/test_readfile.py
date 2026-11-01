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
@pytest.mark.parametrize('metafile_ext', [None, '.xml'])
def test_isce_binary_with_xml(tmp_path, xml_type, metafile_ext):
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

    if metafile_ext:
        # prep_isce explicitly requests XML to avoid stale dimensions in an RSC.
        fname.with_name(fname.name + '.rsc').write_text(
            'WIDTH 1\nFILE_LENGTH 1\nPROCESSOR isce\nDATA_TYPE float32\n',
            encoding='utf-8',
        )

    metadata = readfile.read_attribute(str(fname), metafile_ext=metafile_ext)
    assert metadata['PROCESSOR'] == 'isce'
    assert metadata['WIDTH'] == '4'
    assert metadata['LENGTH'] == '3'
    assert metadata['DATA_TYPE'] == 'float32'


@pytest.mark.parametrize('with_rsc', [False, True])
@pytest.mark.parametrize('sidecar_ext, metafile_ext', [
    ('.xml', '.xml'),
    ('.aux.xml', '.xml'),
    ('.aux.xml', '.aux.xml'),
])
def test_geotiff_with_explicit_qgis_xml(geotiff, with_rsc, sidecar_ext, metafile_ext):
    """An explicitly selected QGIS XML must fall back to raster metadata."""
    expected = readfile.read_attribute(str(geotiff))
    if with_rsc:
        geotiff.with_name(geotiff.name + '.rsc').write_text(
            'WIDTH 4\nFILE_LENGTH 3\nPROCESSOR hyp3\nDATA_TYPE float32\n',
            encoding='utf-8',
        )
    geotiff.with_name(geotiff.name + sidecar_ext).write_text(
        '<PAMDataset><PAMRasterBand band="1"><Metadata>'
        '<MDI key="STATISTICS_MEAN">5.5</MDI>'
        '</Metadata></PAMRasterBand></PAMDataset>',
        encoding='utf-8',
    )

    assert readfile.read_attribute(str(geotiff), metafile_ext=metafile_ext) == expected


@pytest.mark.parametrize('with_rsc', [False, True])
def test_geotiff_with_gdal_statistics(geotiff, with_rsc):
    """Reproduce the sidecar using GDAL itself rather than hand-written XML."""
    if with_rsc:
        geotiff.with_name(geotiff.name + '.rsc').write_text(
            'WIDTH 4\nFILE_LENGTH 3\nPROCESSOR hyp3\nDATA_TYPE float32\n',
            encoding='utf-8',
        )
    expected = readfile.read_attribute(str(geotiff))
    ds = gdal.Open(str(geotiff), gdal.GA_ReadOnly)
    ds.GetRasterBand(1).ComputeStatistics(False)
    ds = None
    assert geotiff.with_name(geotiff.name + '.aux.xml').is_file()

    assert readfile.read_attribute(str(geotiff)) == expected
