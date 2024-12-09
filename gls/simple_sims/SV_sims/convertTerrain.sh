gdal_translate -of XYZ -a_srs EPSG:4326 ./terrain/exportImage.tiff ./terrain/dem.xyz
xyz2kdt -v ./terrain/dem < ./terrain/dem.xyz
rm ./terrain/dem.xyz
