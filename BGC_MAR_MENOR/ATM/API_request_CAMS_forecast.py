### conda activate to_download_ERA5
import cdsapi

dataset = "cams-global-greenhouse-gas-forecasts"
request = {
    "variable": ["co2_column_mean_molar_fraction"],
    "date": ["2024-03-01/2026-01-16"],
    "leadtime_hour": ["0"],
    "data_format": "netcdf_zip",
    "area": [39, -2, 35, 2]
}

client = cdsapi.Client()
client.retrieve(dataset, request).download()
