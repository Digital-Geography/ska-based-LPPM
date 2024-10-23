# coding: utf-8
# -----------------
# created: 20.05.2024
# author: Robby Heusequin
# purpose: ska-based privacy protection method developed for master thesis at the University of Vienna
# ----------------------------------------------

# Import the necessary libraries - note might have to install the libraries
import gpxpy
import requests
import random
import openrouteservice
import os
from geopy.distance import geodesic
from geopy.geocoders import Nominatim
from shapely.geometry import Point, MultiPoint
from shapely.ops import nearest_points


'''
FUNCTIONS
'''


# Function to import gpx file, modify (delete first and last x points) and save it
def modify_gpx_file(input_gpx_file, output_gpx_file):
    # Load the GPX file
    with open(input_gpx_file, 'r') as gpx_file:
        gpx = gpxpy.parse(gpx_file)

    # Access the segment in the gpx file
    segment = gpx.tracks[0].segments[0]

    # Extract the original start and end locations
    original_start_location = (segment.points[0].latitude, segment.points[0].longitude)
    original_end_location = (segment.points[-1].latitude, segment.points[-1].longitude)

    # Print the original start and end locations
    print(f"Original Start Location: {original_start_location}")
    print(f"Original End Location: {original_end_location}")

    # Remove the first 5 and last 5 points
    segment.points = segment.points[5:-5]

    # Extract the new start and end locations
    new_start_location = (segment.points[0].latitude, segment.points[0].longitude)
    new_end_location = (segment.points[-1].latitude, segment.points[-1].longitude)

    # Print the new start and end locations
    print(f"New Start Location: {new_start_location}")
    print(f"New End Location: {new_end_location}")

    # Save the modified GPX file
    with open(output_gpx_file, 'w') as output_gpx_file:
        output_gpx_file.write(gpx.to_xml())

    # Return the original and new locations to use in the algorithm
    locations = {
        'original_start_location': original_start_location,
        'original_end_location': original_end_location,
        'new_start_location': new_start_location,
        'new_end_location': new_end_location
    }

    return locations

# Function to import a gpx file
def import_gpx_file(file_path):
    try:
        with open(file_path, 'r') as gpx_file:
            gpx = gpxpy.parse(gpx_file)
            return gpx
    except FileNotFoundError:
        print("File not found.")


# Function to get the latitude and longitude from an address
def get_lat_lon_from_address(address):
    geolocator = Nominatim(user_agent="ENTER EMAIL ADDRESS")
    lat_lon_address = geolocator.geocode(address)

    # Return the coordinates
    return lat_lon_address.latitude, lat_lon_address.longitude

# Get the nearest neighbour addresses from an input point and get the k nearest addresses
def k_nearest_neighbour_addresses(lat, lon, num_addresses):
    # Create the request to ask all addresses within 1 km
    overpass_url = "http://overpass-api.de/api/interpreter"
    overpass_query = f"[out:json];(node[\"addr:housenumber\"](around:1000,{lat},{lon}););out body;"
    overpass_response = requests.get(overpass_url, params={'data': overpass_query})
    addresses_1000m = overpass_response.json()

    if 'elements' in addresses_1000m:
        addresses = []
        for element in addresses_1000m['elements']:
            if 'tags' in element and 'addr:housenumber' in element['tags']:
                address = element['tags'].get('addr:street', '') + ' ' + element['tags'].get('addr:housenumber', '') + ', ' + element['tags'].get('addr:postcode', '') + ', ' + str(element['lat']) + ', ' + str(element['lon'])
                address_lat = float(element['lat'])
                address_lon = float(element['lon'])
                distance = geodesic((lat, lon), (address_lat, address_lon)).kilometers

                # exclude the original address out of the list with a tolerance of 0,0001 (sometimes the same address is multiple times in the dataset with a slightly different lat and lon)
                tolerance = 0.0001
                if not (abs(address_lat - lat) <= tolerance and abs(address_lon - lon) <= tolerance):
                    addresses.append((address, distance))

                # Create an empty set to store unique addresses - do this because a set does not allow duplicates
                addresses_found = set()
                addresses_no_duplicates = []

                for address, _ in addresses:
                    # Split the address string by comma and take the first part
                    address_street_number = address.split(',')[0].strip()
                    if address_street_number not in addresses_found:
                        addresses_no_duplicates.append((address, _))
                        addresses_found.add(address_street_number)

        # Sort addresses based on distance from the given lat and lon
        # Sorts the list based on the second element of each tuple
        addresses_no_duplicates.sort(key=lambda x: x[1])

        # Return the nearest k addresses
        return [address[0] for address in addresses_no_duplicates[:num_addresses]]
    else:
        return []

# Function to calculate the bounding box as this is needed for the query to get the closest street intersections
def calculate_bbox(gpx_file):
    # Read the GPX file
    with open(gpx_file, 'r') as f:
        gpx = gpxpy.parse(f)

    # Initialize variables to store min and max latitudes and longitudes
    min_lat, max_lat = float('inf'), float('-inf')
    min_lon, max_lon = float('inf'), float('-inf')

    # Iterate over track points
    for track in gpx.tracks:
        for segment in track.segments:
            for point in segment.points:
                # Update min and max latitudes and longitudes
                min_lat = min(min_lat, point.latitude)
                max_lat = max(max_lat, point.latitude)
                min_lon = min(min_lon, point.longitude)
                max_lon = max(max_lon, point.longitude)

    # Add buffer to the bounding box
    buffer = 0.01
    min_lat -= buffer
    max_lat += buffer
    min_lon -= buffer
    max_lon += buffer

    # Define the bounding box (bbox)
    bbox = (min_lat, min_lon, max_lat, max_lon)
    return bbox

# Function to get the closest street intersections
def get_closest_intersections(bbox):
    overpass_url = "http://overpass-api.de/api/interpreter"

    # Overpass Turbo overpass_query template
    overpass_query = f"""
    [out:json][timeout:60][bbox:{bbox[0]},{bbox[1]},{bbox[2]},{bbox[3]}];
    way["highway"~"^(trunk|primary|secondary|tertiary|unclassified|residential|pedestrian|living_street)$"]->.streets;
    node(way_link.streets:3-)->.connections;
    foreach .connections->.connection(
      way(bn.connection);
      if (u(t["name"]) == "< multiple values found >") {{
        (.connection;.intersections;)->.intersections;
      }}
    );
    .intersections out geom;
    """

    overpass_response = requests.post(overpass_url, data={'data': overpass_query})

    # Check if the request was successful
    if overpass_response.status_code == 200:
        # Parse the JSON overpass_response
        intersection_data = overpass_response.json()

        # Extract intersections
        intersections = [
            {'lat': element['lat'], 'lon': element['lon']}
            for element in intersection_data['elements'] if element['type'] == 'node'
        ]

        return intersections
    else:
        print(f"Error: {overpass_response.status_code}")
        print(overpass_response.text)
        return None

# Create a gpx track
def create_gpx_segment(route_coordinates):
    gpx_segment = gpxpy.gpx.GPXTrackSegment()
    for coordinates in route_coordinates:
        lon, lat = coordinates  # OpenRouteService returns coordinates in [lon, lat] format
        gpx_segment.points.append(gpxpy.gpx.GPXTrackPoint(lat, lon))
    return gpx_segment

# Function to calculate the total distance of a track in a gpx file
def calculate_total_distance(file_path):
    # Parse the GPX file
    with open(file_path, 'r') as gpx_file:
        gpx = gpxpy.parse(gpx_file)

    total_distance = 0.0
    for track in gpx.tracks:
        for segment in track.segments:
            total_distance += segment.length_3d()

    return total_distance / 1000.0  # Convert meters to kilometers


'''
DEFINE VARIABLES
'''


# Define the path where the gpx file is saved and a modified version will be saved
input_gpx_file = "DEFINE INPUT GPX FILE"
output_gpx_file = "DEFINE MODIFIED GPX FILE"

# API key for the OpenRouteServices client to request a route
client = openrouteservice.Client(key='ENTER ORS API KEY')

# Import the gpx file and extract the original and temporary start and end point coordinates from the input gpx file
gpx_data_info = modify_gpx_file(input_gpx_file, output_gpx_file)

# Save the coordinates of the original and temporary start and end locations in a variable
original_start = gpx_data_info['original_start_location']
original_end = gpx_data_info['original_end_location']
temporary_start = gpx_data_info['new_start_location']
temporary_end = gpx_data_info['new_end_location']

# Import the modified gpx file with the new temporary start and end locations
gpx_data = import_gpx_file(output_gpx_file)

# Define the k-min and k-max for the nearest addresses
kMin_num_addresses = 5
kMax_num_addresses = 20


'''
START LOCATION
'''


# Get the k-min nearest addresses from the starting location
kMin_nearest_addresses_start = k_nearest_neighbour_addresses(original_start[0], original_start[1], kMin_num_addresses)

# Shifted address k-min from starting address
if kMin_nearest_addresses_start:
    print(f"Nearest {kMin_num_addresses} addresses:")
    for kMin_address_start in kMin_nearest_addresses_start:
        print(kMin_address_start)
else:
    print("No addresses found within the specified range.")

# Select a random address from the list -> k-min shifted address
kMin_shifted_start_address = random.choice(kMin_nearest_addresses_start)
print(f"the shifted start address (k-min) is: {kMin_shifted_start_address}.")

# Save the coordinates from the k-min shifted start address
kMin_shifted_address_start_coordinates = get_lat_lon_from_address(kMin_shifted_start_address)

# Get the k-max nearest addresses from the shifted k-min start location
kMax_nearest_addresses_start = k_nearest_neighbour_addresses(kMin_shifted_address_start_coordinates[0], kMin_shifted_address_start_coordinates[1], kMax_num_addresses)

# Shifted address k-max from the shifted k-min address
if kMax_nearest_addresses_start:
    print(f"Nearest {kMax_num_addresses} addresses:")
    for kMax_address_start in kMax_nearest_addresses_start:
        print(kMax_address_start)
else:
    print("No addresses found within the specified range.")

# Select a random address from the list -> k-max shifted address
kMax_shifted_address_start = random.choice(kMax_nearest_addresses_start)
print(f"the shifted start address (k-max) is: {kMax_shifted_address_start}.")

# Save the coordinates from the k-max shifted start address
kMax_shifted_start_address_coordinates = get_lat_lon_from_address(kMax_shifted_address_start)

# Get the bounding box for querying the street intersections
# bbox (south, west, north, east)
bbox = calculate_bbox(input_gpx_file)
print(f"the bounding box is: {bbox}.")

# Find the nearest intersection point
intersections_start = get_closest_intersections(bbox)

# Create a Shapely Point object for the k-max shifted start address
original_point_start = Point(kMax_shifted_start_address_coordinates[0], kMax_shifted_start_address_coordinates[1])
print(original_point_start)

# Create a list to store Shapely Point objects for the intersections found inside the bounding box
intersection_points = [Point(intersection['lat'], intersection['lon']) for intersection in intersections_start]

# Create a MultiPoint object to store all the intersections
multi_points = MultiPoint(intersection_points)

# Find the nearest intersection point
nearest_intersection_start = nearest_points(original_point_start, multi_points)[1]

# Extract the latitude and longitude of the closest street intersection
nearest_intersection_start_lat = nearest_intersection_start.x
nearest_intersection_start_lon = nearest_intersection_start.y

print(f"The nearest intersection is at {nearest_intersection_start_lat}, {nearest_intersection_start_lon}")

# Save the coordinates in a variable
nearest_intersection_start_lat_lon = nearest_intersection_start_lat, nearest_intersection_start_lon

# Get the coordinates of the closest street intersection and the temporary start location
start_coordinates_lat_lon = ((nearest_intersection_start_lat_lon), (temporary_start))
print(start_coordinates_lat_lon)
# Convert to (longitude, latitude) format - OpenRouteServices needs the coordinates in lon, lat!
start_coordinates_lon_lat = [(lon, lat) for lat, lon in start_coordinates_lat_lon]
print(start_coordinates_lon_lat)

# Request the new walking route from the OpenRouteServices API
shifted_start_route = client.directions(
    coordinates=start_coordinates_lon_lat,
    profile='foot-walking',
    format='geojson'
)

# Extract the coordinates of the start route
start_route_geometry = shifted_start_route['features'][0]['geometry']
start_route_coordinates = start_route_geometry['coordinates']

# Create a new GPX track from the route coordinates for the start location
start_gpx_segment = create_gpx_segment(start_route_coordinates)

# Insert the new segment into the GPX track at the first position
gpx_data.tracks[0].segments.insert(0, start_gpx_segment)


'''
END LOCATION
'''


# Get the k-min nearest addresses from the end location
kMin_nearest_addresses_end = k_nearest_neighbour_addresses(original_end[0], original_end[1], kMin_num_addresses)

# Shifted address k-min from end address
if kMin_nearest_addresses_end:
    print(f"Nearest {kMin_num_addresses} addresses:")
    for kMin_address_end in kMin_nearest_addresses_end:
        print(kMin_address_end)
else:
    print("No addresses found within the specified range.")

# Select a random address from the list -> k-min shifted address
kMin_shifted_address_end = random.choice(kMin_nearest_addresses_end)
print(f"the shifted end address (kMin) is: {kMin_shifted_address_end}.")

# Save the coordinates from the k-min shifted end address
kMin_shifted_address_end_coordinates = get_lat_lon_from_address(kMin_shifted_address_end)

# Get the k-max nearest addresses from the shifted k-min end location
kMax_nearest_addresses_end = k_nearest_neighbour_addresses(kMin_shifted_address_end_coordinates[0], kMin_shifted_address_end_coordinates[1], kMax_num_addresses)

# Shifted address k-max from the shifted k-min address
if kMax_nearest_addresses_end:
    print(f"Nearest {kMax_num_addresses} addresses:")
    for kMax_address_end in kMax_nearest_addresses_end:
        print(kMax_address_end)
else:
    print("No addresses found within the specified range.")

# Select a random address from the list -> k-max shifted address
kMax_shifted_address_end = random.choice(kMax_nearest_addresses_end)
print(f"the shifted end address (kMax) is: {kMax_shifted_address_end}.")

# Save the coordinates from the k-max shifted end address
kMax_shifted_address_end_coordinates = get_lat_lon_from_address(kMax_shifted_address_end)

# Find the nearest intersection point
intersections_end = get_closest_intersections(bbox)

# Create a Shapely Point object for the k-max shifted end address
kMax_shifted_point_end = Point(kMax_shifted_address_end_coordinates[0], kMax_shifted_address_end_coordinates[1])

# Find the nearest intersection
nearest_intersection_end = nearest_points(kMax_shifted_point_end, multi_points)[1]

# Extract the latitude and longitude of the closest street intersection
nearest_intersection_end_lat = nearest_intersection_end.x
nearest_intersection_end_lon = nearest_intersection_end.y

print(f"The nearest intersection is at {nearest_intersection_end_lat}, {nearest_intersection_end_lon}")

# Save the coordinates in a variable
nearest_intersection_end_lat_lon = nearest_intersection_end_lat, nearest_intersection_end_lon

# Get the coordinates of the temporary end location and the closest street intersection
end_coordinates_lat_lon = ((temporary_end), (nearest_intersection_end_lat_lon))
print(end_coordinates_lat_lon)
# Convert to (longitude, latitude) format - OpenRouteServices needs the coordinates in lon, lat!
end_coordinates_lon_lat = [(lon, lat) for lat, lon in end_coordinates_lat_lon]
print(end_coordinates_lon_lat)

# Request the new walking route from the OpenRouteServices API
shifted_end_route = client.directions(
    coordinates=end_coordinates_lon_lat,
    profile='foot-walking',
    format='geojson'
)

# Extract the coordinates of the end route
end_route_geometry = shifted_end_route['features'][0]['geometry']
end_route_coordinates = end_route_geometry['coordinates']

# Create a new GPX track from the route coordinates for the end location
end_gpx_segment = create_gpx_segment(end_route_coordinates)

# Insert the new segment into the GPX track at the second position
gpx_data.tracks[0].segments.insert(2, end_gpx_segment)


'''
SAVING THE MODIFIED GPX FILE
'''


# Define the base path and the masked GPX file name
dir_name = "DEFINE BASE PATH"
masked_name = "DEFINE MASKED GPX FILE NAME"
masked_gpx_file = os.path.join(dir_name + masked_name)

# Save the new GPX file
with open(masked_gpx_file, 'w') as f:
    f.write(gpx_data.to_xml())

print(f"GPX file has been saved as {masked_name}")

# Calculate the original and masked total distance travelled in the GPX files
total_distance_original_gpx = calculate_total_distance(input_gpx_file)
total_distance_masked_gpx = calculate_total_distance(masked_gpx_file)

print(f"The total distance travelled of the original gpx file: {total_distance_original_gpx:.2f} km.")
print(f"The total distance travelled of the masked gpx file: {total_distance_masked_gpx:.2f} km.")