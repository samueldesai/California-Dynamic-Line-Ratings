import cdsapi
import pandas as pd
import numpy as np
import os
import tempfile
import zipfile
import math
from datetime import datetime
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

# Test mode flag - set to True to only process the first bus
TEST_MODE = False

def extract_csv_from_zip(zip_filepath):
    try:
        with zipfile.ZipFile(zip_filepath, 'r') as zip_ref:
            csv_files = [f for f in zip_ref.namelist() if f.endswith('.csv')]
            if not csv_files:
                raise ValueError("No CSV found in ZIP")
            csv_name = csv_files[0]
            with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as tmp:
                tmp.write(zip_ref.read(csv_name))
                return tmp.name
    except Exception as e:
        print(f"Error extracting ZIP: {e}")
        return None

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
    u_cols = [f"{m}_U_H{h:02d}" for m in months for h in range(24)]
    v_cols = [f"{m}_V_H{h:02d}" for m in months for h in range(24)]
    t_cols = [f"{m}_Temp_H{h:02d}" for m in months for h in range(24)]
    p_cols = [f"{m}_Pressure_H{h:02d}" for m in months for h in range(24)]
    new_df = pd.DataFrame(np.nan, index=bus_df.index, columns=u_cols + v_cols + t_cols + p_cols)
    return pd.concat([bus_df, new_df], axis=1), u_cols, v_cols, t_cols, p_cols

def populate_weather_data(output_df, u_cols, v_cols, t_cols, p_cols, bus_idx, weather_avg):
    if weather_avg is None:
        return
    for _, row in weather_avg.iterrows():
        m = int(row['month']) - 1
        h = int(row['hour'])
        ucol = u_cols[m*24 + h]
        vcol = v_cols[m*24 + h]
        tcol = t_cols[m*24 + h]
        pcol = p_cols[m*24 + h]
        output_df.at[bus_idx, ucol] = row['u']
        output_df.at[bus_idx, vcol] = row['v']
        output_df.at[bus_idx, tcol] = row['temp']
        output_df.at[bus_idx, pcol] = row['pressure']

def fetch_and_process(coord_key, lon, lat):
    client = cdsapi.Client(url="https://cds.climate.copernicus.eu/api",
                           key="a08eb70e-b147-43de-8a30-2fa79f1a3558")
    req = {
        "variable": ["2m_temperature", "10m_u_component_of_wind", "10m_v_component_of_wind", "surface_pressure"],
        "location": {"longitude": lon, "latitude": lat},
        "date":    ["2020-01-01/2024-12-31"],
        "data_format": "csv"
    }
    with tempfile.NamedTemporaryFile(suffix=".zip", delete=False) as tmp:
        zip_path = tmp.name
    try:
        print(f"Calling API for {coord_key}")
        client.retrieve("reanalysis-era5-single-levels-timeseries", req).download(zip_path)
        csv_path = extract_csv_from_zip(zip_path)
        os.remove(zip_path)
        if csv_path is None:
            return coord_key, None
        df = pd.read_csv(csv_path)
        os.remove(csv_path)
        df['datetime'] = pd.to_datetime(df.iloc[:,0], utc=True).dt.tz_convert('US/Pacific')
        df['month'] = df['datetime'].dt.month
        df['hour']  = df['datetime'].dt.hour
        df['u'] = df.iloc[:,1]
        df['v'] = df.iloc[:,2]
        df['temp'] = df.iloc[:,3]
        df['pressure'] = df.iloc[:,4]
        avg = df.groupby(['month','hour'])[['u','v','temp','pressure']].mean().reset_index()
        return coord_key, avg
    except Exception as e:
        print(f"Failed to fetch/process {coord_key}: {e}")
        return coord_key, None

def print_remaining_calls(processed, total):
    remaining = total - processed
    print(f"API calls remaining: {remaining}")
    return remaining

def main():
    bus_df = pd.read_csv('bus_data_new.csv')
    if TEST_MODE:
        bus_df = bus_df.head(1)

    coord_groups = analyze_unique_coordinates(bus_df)
    total_calls = len(coord_groups)
    processed_calls = 0

    output_df, u_cols, v_cols, t_cols, p_cols = create_output_dataframe(bus_df)
    weather_cache = {}

    max_workers = min(32, os.cpu_count() or 4)
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = {
            pool.submit(fetch_and_process, key, grp[0]['longitude'], grp[0]['latitude']): key
            for key, grp in coord_groups.items()
        }

        for fut in as_completed(futures):
            key = futures[fut]
            _, avg_df = fut.result()
            weather_cache[key] = avg_df
            processed_calls += 1
            print(f"Processed {key}")
            print_remaining_calls(processed_calls, total_calls)

    for key, buses in coord_groups.items():
        avg = weather_cache.get(key)
        for b in buses:
            populate_weather_data(output_df, u_cols, v_cols, t_cols, p_cols, b['index'], avg)

    out_name = 'bus_weather_data_parallel.csv' if not TEST_MODE else 'bus_weather_data_test.csv'
    output_df.to_csv(out_name, index=False)
    print(f"Saved results to {out_name}")

if __name__ == "__main__":
    main()