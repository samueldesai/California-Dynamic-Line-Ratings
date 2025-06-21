import math
import json
import csv
import pandas as pd
import numpy as np

# --- Constants used for thermal rating calculations ---

avogadros_number = 6.023e23  # Avogadro’s number (unused here)
stefan_boltzmann = 5.67e-8  # Stefan–Boltzmann constant (W/m²K⁴)
specific_air = 287.058  # Specific gas constant for dry air (J/kg·K)
absorptivity = 0.8  # Conductor absorptivity (dimensionless)
emissivity = 0.8  # Conductor emissivity (dimensionless)
temp_surface = 75 + 273.15  # Assumed conductor temperature (K)
resistance_per_meter = 9.391e-5
temp_SLR = 40 + 273.15  # Surface layer air temperature (K)
GHI_slr = 1000  # Global Horizontal Irradiance (W/m²)
pressure_SLR = 101325  # Atmospheric pressure for SLR (Pa)
Tc= 100
resistance_per_meter = 9.391e-5

# --- Function to calculate a static (standard) rating ---
def calculate_standard_rating(diameter, resist):
    """Estimate thermal rating using reference (static) conditions."""
    K = 1.194 - math.cos(0) + 0.194 * math.cos(0) + 0.368 * math.sin(0)
    density_air = pressure_SLR / (specific_air * temp_SLR)  # Air density
    reynolds = diameter * density_air * speed_SLR / dynvis_air  # Reynolds number

    # Three empirical convective cooling approximations
    qc_zero = 3.645 * (density_air**0.5) * (diameter**0.75) * ((temp_c - temp_SLR)**1.25)
    qc_low = K * (1.01 + 1.35 * reynolds**0.52) * k_air * (temp_c - temp_SLR)
    qc_high = K * (0.754 * reynolds**0.6) * k_air * (temp_c - temp_SLR)
    q_c = max(qc_zero, qc_low, qc_high)  # Use highest convective cooling

    # Radiative and solar heat exchange
    q_r = math.pi * diameter * stefan_boltzmann * emissivity * ((temp_c**4) - (temp_SLR**4))
    q_s = absorptivity * diameter * GHI_slr  # Solar heat gain

    # Thermal rating in Amps (with unit scaling)
    max_rating = math.sqrt((q_c + q_r - q_s) / resist) / 100
    return max_rating

# --- Function to calculate a dynamic line rating ---
import math

def calculate_dynamic_rating(
    diameter, 
    voltage, 
    dlr_pressure, 
    dlr_temp, 
    dlr_speed, 
    dlr_theta, 
    dlr_ghi,
):
    specific_air = 287.0  # J/(kg·K)
    K = 1.194 - math.cos(dlr_theta) + 0.194 * math.cos(2 * dlr_theta) + 0.368 * math.sin(2 * dlr_theta)
    tfilm = (Tc + dlr_temp) / 2
    delta_temp = Tc - dlr_temp

    # Air properties
    density_air = dlr_pressure / (specific_air * (tfilm + 273.15))
    k_air = 2.424e-2 + 7.477e-5 * tfilm - 4.407e-9 * (tfilm**2)
    dynvisair = (1.458e-6 * (tfilm + 273.15)**1.5) / (tfilm + 383.4)
    reynolds = diameter * density_air * dlr_speed / dynvisair

    # Convective cooling (IEEE 738 compliant selection)
    q1 = 3.645 * density_air**0.5 * diameter**0.75 * delta_temp**1.25
    q2 = K * (1.01 + 1.35 * reynolds**0.52) * k_air * delta_temp
    q3 = K * (0.754 * reynolds**0.6) * k_air * delta_temp

    q_c = max(q1,q2,q3)
    

    # Radiative heat loss (fixed ambient term)
    term_cond = (Tc + 273.15) / 100.0
    term_amb = (dlr_temp + 273.15) / 100.0
    q_r = 17.8 * diameter * emissivity * (term_cond**4 - term_amb**4)

    # Solar heat gain (requires solar geometry adjustment)
    q_s = absorptivity * diameter * dlr_ghi  # Simplified; use effective irradiance
    
    # Current calculation (amperes)
    I = math.sqrt((q_c + q_r - q_s) / resistance_per_meter) / 1000

    return I * voltage * (3** (1/2)) / 100


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


# --- Load JSON line data ---
with open("C:\\Users\\gogre\\Desktop\\191\\Final project\\CATS_lines.json", "r") as f:
    data = json.load(f)

# --- Initialize dictionary to store line metadata ---
data_dict = {}
for idx, feature in enumerate(data["features"]):
    props = feature["properties"]
    data_dict[idx] = {
        "Resistance": props["br_r"],
        "Voltage": props["kV"],
        "From": props["f_bus"],
        "To": props["t_bus"],
        "Rating_i": props["rate_a"],  # Initial rating
        "Rating_f": []  # Dynamic ratings over time
    }
    coordinates = feature["geometry"]["coordinates"]

    # Calculate total line distance
    distance = 0.0
    for i in range(1, len(coordinates)):
        distance += haversine(coordinates[i-1], coordinates[i])
    data_dict[idx]["Resistance"] = resistance_per_meter * (distance * 1000)  # in kilometers


    # Set conductor diameter based on voltage level (converted from inches to meters)
    voltage = props["kV"]
    if voltage == 66.0:
        data_dict[idx]["Diameter"] = 1.108 / 39.37
    elif voltage in [115.0, 230.0, 500.0]:
        data_dict[idx]["Diameter"] = 1.108 / 39.37
    else:
        data_dict[idx]["Diameter"] = None  # Handle unknown voltages

# --- Load weather and angle data from CSVs ---
csv1_data = pd.read_csv("bus_weather_data_parallel.csv")
csv2_data = pd.read_csv("bus_ghi_data.csv")
angles_data = pd.read_csv("line_angles.csv")

# --- Parse weather CSVs ---
csv1_bus_data = {}
for _, row1 in csv1_data.iterrows():
    bus = int(row1[0])
    csv1_bus_data[bus] = {
        "U_wind": row1[3:291].astype(float).tolist(),
        "V_wind": row1[291:579].astype(float).tolist(),
        "Temperature": row1[579:867].astype(float).tolist(),
        "Pressure": row1[867:1155].astype(float).tolist(),
        "GHI": []
    }

for _, row2 in csv2_data.iterrows():
    bus = int(row2[0])
    csv1_bus_data[bus]["GHI"] = row2[3:291].astype(float).tolist()

# --- Helper: return values or default fallback ---
def get_values(data, key, default):
    return np.array(data.get(key, [default] * len(csv1_bus_data[1]["U_wind"])))

# --- Merge and process weather + geometry for each line ---
for idx in range(len(data_dict)):
    loc_dict = data_dict[idx]
    from_bus = loc_dict["From"]
    to_bus = loc_dict["To"]
    
    # Get weather data from each bus end
    from_dict = csv1_bus_data.get(from_bus, {})
    to_dict = csv1_bus_data.get(to_bus, {})

    # Retrieve and average values
    from_u, to_u = get_values(from_dict, "U_wind", 0.61), get_values(to_dict, "U_wind", 0.61)
    from_v, to_v = get_values(from_dict, "V_wind", 0), get_values(to_dict, "V_wind", 0)
    from_temp, to_temp = get_values(from_dict, "Temperature", temp_SLR), get_values(to_dict, "Temperature", temp_SLR)
    from_pressure, to_pressure = get_values(from_dict, "Pressure", pressure_SLR), get_values(to_dict, "Pressure", pressure_SLR)
    from_ghi, to_ghi = get_values(from_dict, "GHI", GHI_slr), get_values(to_dict, "GHI", GHI_slr)


    # PROPER WIND CALCULATIONS
    # 1. Compute vector averages first
    U_avg = (from_u + to_u) / 2
    V_avg = (from_v + to_v) / 2
    wind_speed_avg = np.sqrt(U_avg**2 + V_avg**2)
    wind_dir_avg = np.arctan2(V_avg, U_avg)  # True vector mean direction

    # Store averaged values
    loc_dict["average_wind"] = wind_speed_avg.tolist()
    loc_dict["average_temp"] = ((from_temp + to_temp) / 2).tolist()
    loc_dict["average_pressure"] = ((from_pressure + to_pressure) / 2).tolist()
    loc_dict["average_ghi"] = ((from_ghi + to_ghi) / 2).tolist()

    # LINE ANGLE PROCESSING
    line_angle = angles_data.iloc[idx, 1] if idx < len(angles_data) else None
    loc_dict["Angle (Radians)"] = line_angle  # Single value per line

    if line_angle is not None:
        # Calculate acute angle between wind and line (0-90°)
        angle_diff = (wind_dir_avg - line_angle + np.pi) % (2 * np.pi) - np.pi
        acute_angle = np.minimum(np.abs(angle_diff), np.pi - np.abs(angle_diff))
        loc_dict["Angle_Diff"] = acute_angle.tolist()  # Use acute angle for IEEE 738 K-factor
    else:
        loc_dict["Angle_Diff"] = [np.pi/2] * len(loc_dict["average_wind"]) 

# --- Compute dynamic ratings for each hour ---



for line in data_dict:
    line_data = data_dict[line]
    line = data_dict[0]
    for t in range(len(line_data["average_wind"])):
        if line_data["Diameter"] is None:
            line_data["Rating_f"].append(9999)  # Default for unknown diameter
        else:
            rating = calculate_dynamic_rating(
                line_data["Diameter"],
                line_data["Voltage"],
                line_data["average_pressure"][t],
                (line_data["average_temp"][t]) - 273.15,
                line_data["average_wind"][t],
                line_data["Angle_Diff"][t],
                line_data["average_ghi"][t])
            
            line_data["Rating_f"].append(rating)

# --- Write results to CSV file ---
csv_file_path = "TEST_DLR_Ratings.csv"
months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", 
          "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
hour_labels = [f"{month} {hour}" for month in months for hour in range(24)]

with open(csv_file_path, mode="w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["From", "To"] + hour_labels)
    for idx in data_dict:
        d = data_dict[idx]
        writer.writerow([d["From"], d["To"]] + d.get("Rating_f", [None]*288))

print(f"CSV file saved successfully at {csv_file_path}")
