import math
import json
import csv
import pandas as pd

# Define constants with appropriate units
temp_SLR = 40 + 273.15  # Surface layer air temperature (K)
GHI_slr = 1000  # Global Horizontal Irradiance (W/m²)
pressure_SLR = 101325  # Atmospheric pressure for SLR (Pa)
speed_SLR = 0.61  # Wind speed for SLR (m/s)
k_air = 0.027  # Thermal conductivity of air (W/m·K)
dynvis_air = 1.908e-5  # Dynamic viscosity of air (Pa·s)
avogadros_number = 6.023e23  # Avogadro's number (molecules/mol)
stefan_boltzmann = 5.67e-8  # Stefan-Boltzmann constant (W/m²K⁴)
specific_air = 287.058  # Specific gas constant for dry air (J/kg·K)
absorptivity = 0.8  # Absorptivity of conductor (dimensionless)
emissivity = 0.8  # Emissivity of conductor (dimensionless)
temp_c = 75 + 273.15  # Conductor temperature (K)
resistance_per_meter = 9.391e-5

theta_slr = math.pi/2

def calculate_standard_rating(diameter, resist, voltage):
    """Calculate the standard thermal rating of a transmission line conductor."""

    # Maximum convective heat transfer
    q_c = 82.10
    # Radiative heat loss (W/m)
    q_r = 39.11
    # Solar heat absorption (W/m)
    q_s = 22.45

    # Maximum conductor current rating (A)
    max_rating = math.sqrt((q_c + q_r - q_s) / resistance_per_meter ) 

    return ((max_rating / 1000) * voltage * (3 ** (1/2))) / 100 ### per unit

def haversine(coord1, coord2):
    # Coordinates in [lon, lat]
    lon1, lat1 = coord1
    lon2, lat2 = coord2
    # Convert degrees to radians
    lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])
    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
    c = 2 * math.asin(math.sqrt(a))
    r = 6371.0  # Earth radius in kilometers
    return c * r

# Load transmission line data from JSON file
with open("C:\\Users\\gogre\\Desktop\\191\\Final project\\CATS_lines.json", "r") as f:
    data = json.load(f)

# Initialize dictionary to store transmission line properties
data_dict = {}

# Extract relevant properties from JSON data
for idx, feature in enumerate(data["features"]):
    props = feature["properties"]
    data_dict[idx] = {}
    index_dict = data_dict[idx]
    index_dict["Voltage"] = props["kV"]
    coordinates = feature["geometry"]["coordinates"]

    # Calculate total line distance
    distance = 0.0
    for i in range(1, len(coordinates)):
        distance += haversine(coordinates[i-1], coordinates[i])

    index_dict["Resistance"] = (distance * 1000)  # in kilometers

    index_dict["From"] = props["f_bus"]
    index_dict["To"] = props["t_bus"]
    index_dict["Rating_i"] = props["rate_a"]
    index_dict["Rating_f"] = 0


    ###setting diameters
    if index_dict["Voltage"] == 66.0:
        index_dict["Diameter"] = 1.108 / 39.37  
    elif index_dict["Voltage"] == 115.0:
        index_dict["Diameter"] = 1.108 / 39.37  
    elif index_dict["Voltage"] == 230.0:
        index_dict["Diameter"] = 1.108 / 39.37  
    elif index_dict["Voltage"] == 500.0:
        index_dict["Diameter"] = 1.108 / 39.37  
    else:
        index_dict["Diameter"] = None


    ##rating for transformers
    if index_dict["Diameter"] is None:
        index_dict["Rating_f"] +=  9999
    else:
        resistance = index_dict["Resistance"]
        diameter = index_dict["Diameter"]
        voltage = index_dict["Voltage"]
        rating = calculate_standard_rating(diameter, resistance, voltage)
        index_dict["Rating_f"] += rating
    

    ### adding the data to the dictionary 

    index_dict["U_wind"] = []
    index_dict["V_wind"] = []
    index_dict["Wind"] = []
    index_dict["Pressure"] = []
    index_dict["GHI"] = []



    

# Define file path for saving CSV data
csv_file_path = "test.csv"

# Save the transmission line data to CSV
with open(csv_file_path, mode="w", newline="") as file:
    writer = csv.writer(file)

    # Write CSV header
    writer.writerow(["From", "To", "Rating"])

    # Write data rows
    for idx in data_dict:
        index_dict = data_dict[idx]
        writer.writerow([
            index_dict["From"],
            index_dict["To"],
            index_dict["Rating_f"]
        ])

print(f"CSV file saved successfully at {csv_file_path}")