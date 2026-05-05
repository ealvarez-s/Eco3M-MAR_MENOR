### conda activate to_download_ERA5
import cdsapi

dataset = "cams-global-ghg-reanalysis-egg4"
request = {
    "date": ["2011-01-01/2020-12-31"],
    "step": ["12"],
    "data_format": "netcdf_zip",
    "variable": ["co2_column_mean_molar_fraction"],
    "area": [39, -2, 35, 2]
}

client = cdsapi.Client()
client.retrieve(dataset, request).download()


dataset = "cams-global-ghg-reanalysis-egg4-monthly"
request = {
    "year": [
        "2003", "2004", "2005",
        "2006", "2007", "2008",
        "2009", "2010", "2011",
        "2012", "2013", "2014",
        "2015", "2016", "2017",
        "2018", "2019", "2020"
    ],
    "month": [
        "01", "02", "03",
        "04", "05", "06",
        "07", "08", "09",
        "10", "11", "12"
    ],
    "product_type": ["monthly_mean"],
    "data_format": "netcdf_zip",
    "variable": ["co2_column_mean_molar_fraction"],
    "area": [39, -2, 35, 2]
}

client = cdsapi.Client()
client.retrieve(dataset, request).download()
