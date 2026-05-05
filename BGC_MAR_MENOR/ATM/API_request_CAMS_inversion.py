### conda activate to_download_ERA5
import cdsapi

dataset = "satellite-carbon-dioxide"
request = {
    "processing_level": ["level_3"],
    "variable": "xco2",
    "sensor_and_algorithm": ["merged_obs4mips"],
    "version": ["4_6"]
}

client = cdsapi.Client()
client.retrieve(dataset, request).download()
