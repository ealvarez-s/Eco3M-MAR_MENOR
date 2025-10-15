import copernicusmarine

copernicusmarine.subset(
  dataset_id="med-ogs-co2-rean-d",
  variables=["fpco2", "spco2"],
  minimum_longitude=-1,
  maximum_longitude=2,
  minimum_latitude=36,
  maximum_latitude=39,
  start_datetime="2010-01-01T00:00:00",
  end_datetime="2010-01-31T23:59:59",
  minimum_depth=0,
  maximum_depth=200,
  output_filename = "med-ogs-co2-rean-d.nc",
  output_directory = "reanalysis"
)
