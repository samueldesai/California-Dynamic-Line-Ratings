import json
import math

# Function to calculate angle using atan2 to handle all directions and vertical lines
def calculateangle(x1, y1, x2, y2):
    dx = x2 - x1
    dy = y2 - y1
    return (math.atan2(dy, dx)) 

# Load your GeoJSON data
with open("C:\\Users\\gogre\\Desktop\\191\\Final project\\CATS_lines.json", "r") as f:
    data = json.load(f)

# Initialize dictionary to store scalar coordinate values
coordinates_dict = {
    "long1": [],
    "lat1": [],
    "long2": [],
    "lat2": []
}

# Extract start and end points from the properties dictionary
for feature in data["features"]:
    props = feature["properties"]  # Accessing the 'properties' dictionary
    
    # Extract individual longitude and latitude values
    long1 = props["Lon1"]
    lat1 = props["Lat1"]
    long2 = props["Lon2"]
    lat2 = props["Lat2"]

    # Append individual float values to respective lists
    coordinates_dict["long1"].append(long1)
    coordinates_dict["lat1"].append(lat1)
    coordinates_dict["long2"].append(long2)
    coordinates_dict["lat2"].append(lat2)

# Compute angles for each line
angles_list = []

for i in range(len(data["features"])):
    x1 = coordinates_dict["long1"][i]
    y1 = coordinates_dict["lat1"][i]
    x2 = coordinates_dict["long2"][i]
    y2 = coordinates_dict["lat2"][i]
    angle_rad = calculateangle(x1, y1, x2, y2)
    angles_list.append(angle_rad)

import csv

# Define the output file name
csv_filename = "angles_output.csv"

# Write angles to a CSV file with index numbers
with open(csv_filename, mode="w", newline="") as file:
    writer = csv.writer(file)
    
    # Write header
    writer.writerow(["Index", "Angle (Radians)"])
    
    # Write each index and angle
    for index, angle in enumerate(angles_list):
        writer.writerow([index, angle])

print(f"Angles successfully exported to {csv_filename}")