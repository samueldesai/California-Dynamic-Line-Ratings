import pandas as pd
import numpy as np
import requests
import json
import os
from datetime import datetime
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import time

# Test mode flag - set to True to only process the first bus
TEST_MODE = False

def create_coordinate_key(longitude, latitude, precision=2):
    return f"{round(float(longitude), precision)},{round(float(latitude), precision)}"

def analyze_unique_coordinates(bus_df):
    groups = defaultdict(list)
    for idx, row in bus_df.iterrows():
        key = create_coordinate_key(row['Lon'], row['Lat'])
        groups[key].append({
            'index': idx,
            'bus_number': row['Bus number'],
            'longitude': row['Lon'],
            'latitude': row['Lat']
        })
    print(f"Unique coordinate groups: {len(groups)} (from {len(bus_df)} buses)")
    return groups

def create_output_dataframe(bus_df):
    months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    ghi_cols = [f"{m}_GHI_H{h:02d}" for m in months for h in range(24)]
    new_df = pd.DataFrame(np.nan, index=bus_df.index, columns=ghi_cols)
    return pd.concat([bus_df, new_df], axis=1), ghi_cols

def populate_solar_data(output_df, ghi_cols, bus_idx, solar_data):
    if solar_data is None:
        return
    for _, row in solar_data.iterrows():
        m = int(row['month']) - 1
        h = int(row['hour'])
        ghi_col = ghi_cols[m*24 + h]
        output_df.at[bus_idx, ghi_col] = row['ghi']

def fetch_pvgis_data(coord_key, lon, lat, max_retries=3):
    """
    Fetch TMY data from PVGIS API
    Uses the TMY (Typical Meteorological Year) endpoint
    """
    base_url = "https://re.jrc.ec.europa.eu/api/v5_2/tmy"
    
    params = {
        'lat': lat,
        'lon': lon,
        'outputformat': 'json'
    }
    
    for attempt in range(max_retries):
        try:
            print(f"Fetching PVGIS TMY data for {coord_key} (attempt {attempt + 1})")
            
            response = requests.get(base_url, params=params, timeout=30)
            response.raise_for_status()
            
            data = response.json()
            
            # Extract TMY data
            if 'outputs' not in data or 'tmy_hourly' not in data['outputs']:
                print(f"No TMY data available for {coord_key}")
                return coord_key, None
            
            tmy_data = data['outputs']['tmy_hourly']
            
            # Convert to DataFrame
            df = pd.DataFrame(tmy_data)
            
            # Parse datetime and convert from UTC to Pacific time
            df['datetime'] = pd.to_datetime(df['time(UTC)'], format='%Y%m%d:%H%M', utc=True).dt.tz_convert('US/Pacific')
            df['month'] = df['datetime'].dt.month
            df['hour'] = df['datetime'].dt.hour
            
            # Extract GHI (Global Horizontal Irradiance) in W/m²
            df['ghi'] = df['G(h)']  # G(h) is Global irradiance on horizontal plane
            
            # Group by month and hour to get average values
            avg_data = df.groupby(['month', 'hour'])[['ghi']].mean().reset_index()
            
            return coord_key, avg_data
            
        except requests.exceptions.RequestException as e:
            print(f"Request failed for {coord_key} (attempt {attempt + 1}): {e}")
            if attempt < max_retries - 1:
                time.sleep(2 ** attempt)  # Exponential backoff
            else:
                print(f"Failed to fetch data for {coord_key} after {max_retries} attempts")
                return coord_key, None
        except Exception as e:
            print(f"Error processing data for {coord_key}: {e}")
            return coord_key, None

def print_remaining_calls(processed, total):
    remaining = total - processed
    print(f"API calls remaining: {remaining}")
    return remaining

def main():
    bus_df = pd.read_csv('bus_data_missing.csv')
    if TEST_MODE:
        bus_df = bus_df.head(1)

    coord_groups = analyze_unique_coordinates(bus_df)
    total_calls = len(coord_groups)
    processed_calls = 0

    output_df, ghi_cols = create_output_dataframe(bus_df)
    solar_cache = {}

    # Use fewer workers to be respectful to the API
    max_workers = min(8, os.cpu_count() or 4)
    
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = {
            pool.submit(fetch_pvgis_data, key, grp[0]['longitude'], grp[0]['latitude']): key
            for key, grp in coord_groups.items()
        }

        for fut in as_completed(futures):
            key = futures[fut]
            _, avg_df = fut.result()
            solar_cache[key] = avg_df
            processed_calls += 1
            print(f"Processed PVGIS data for {key}")
            print_remaining_calls(processed_calls, total_calls)
            
            # Small delay to be respectful to the API
            time.sleep(0.1)

    # Populate the output dataframe with solar irradiance data
    for key, buses in coord_groups.items():
        avg = solar_cache.get(key)
        for b in buses:
            populate_solar_data(output_df, ghi_cols, b['index'], avg)

    out_name = 'bus_pvgis_ghi_data_missing.csv' if not TEST_MODE else 'bus_pvgis_ghi_test.csv'
    output_df.to_csv(out_name, index=False)
    print(f"Saved PVGIS GHI results to {out_name}")

if __name__ == "__main__":
    main()