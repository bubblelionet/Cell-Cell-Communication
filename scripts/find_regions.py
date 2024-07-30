import scanpy as sc
import numpy as np
import csv
import pickle
import math
import json
from itertools import combinations
import pandas as pd
import matplotlib.pyplot as plt
from pyvis.network import Network
import networkx as nx
from networkx.algorithms import bipartite
from cdlib import algorithms
from cdlib import NodeClustering
import altair as alt
import sys
import os


sys.path.append("/Users/victoriagao/Documents/MSc/Schwartz_lab/altair-themes/")
if True:  # In order to bypass isort when saving
    import altairThemes

spot_diameter = 89.43  # pixels


def load_gene_ids(file_path):
    gene_ids = []
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            gene_ids.append(line)
    return gene_ids


def load_coordinates(file_path):
    return np.load(file_path)


def load_pickle(file_path):
    with open(file_path, 'rb') as file:
        return pickle.load(file)


def load_pathologist_labels(file_path):
    pathologist_label = []
    with open(file_path) as file:
        csv_file = csv.reader(file, delimiter=",")
        for line in csv_file:
            pathologist_label.append(line)
    return pathologist_label


def create_barcode_serial(cell_barcodes):
    barcode_serial = {cell_code: i for i, cell_code in enumerate(cell_barcodes)}
    return barcode_serial


def create_barcode_info(cell_barcodes, coordinates):
    barcode_info = [["cell_code", "x-coordinate", "y-coordinate", "component", "color", "barcode_number"]]
    for i, cell_code in enumerate(cell_barcodes):
        barcode_info.append([cell_code, coordinates[i, 0], coordinates[i, 1], 0, 0, 0])
    return barcode_info


def combine_barcodes(barcode_type, barcode_coordinate_dict): # record the type (annotation) of each spot (barcode)
    combined_barcode = {key: [barcode_type[key]] + barcode_coordinate_dict[key] for key in barcode_type if key in barcode_coordinate_dict}
    return combined_barcode


def create_combined_barcode_list(combined_barcode):
    combined_barcode_list = [{'barcode': key, 'type': values[0], 'x': values[1], 'y': values[2]} for key, values in combined_barcode.items()]
    return combined_barcode_list


def euclidean_distance(p1, p2):
    return math.sqrt((p1['x'] - p2['x'])**2 + (p1['y'] - p2['y'])**2)


def find_extreme_distances(combined_barcode_list):
    max_distance = 0
    max_distance_barcodes = ()
    min_distance = float('inf')
    closest_barcodes = ()
    """
    Calculates euclidean distance between each one spot and every neighboring spots. Outputs max and min distance, useful for find_connected_components function
    Input: each barcode with their coordinates.
    Returns: max, min distances and barcodes.
    """

    for i, dict1 in enumerate(combined_barcode_list): # Iterate through all pairs of dictionaries and calculate distances - Max and Min
        for j, dict2 in enumerate(combined_barcode_list[i + 1:], start=i + 1):
            distance = euclidean_distance(dict1, dict2)
            if distance > max_distance:
                max_distance = distance
                max_distance_barcodes = (dict1['barcode'], dict2['barcode'])
            if distance < min_distance:
                min_distance = distance
                closest_barcodes = (dict1['barcode'], dict2['barcode'])
    return max_distance, max_distance_barcodes, min_distance, closest_barcodes


def find_connected_components(combined_barcode_list, proximity_threshold):
    G = nx.Graph()
    """
    Given a proximity_threshold, which is the distance between two neighboring spots, create edges between them.
    Output:
        filtered_connected_components: list of sets, each set is a tumor 'region'
    """
    for component in combined_barcode_list:
        G.add_node(component['barcode'], x=component['x'], y=component['y'], type=component['type'])

    for i, component1 in enumerate(combined_barcode_list):
        for j, component2 in enumerate(combined_barcode_list[i + 1:], start=i + 1):
            distance = ((component1['x'] - component2['x'])**2 + (component1['y'] - component2['y'])**2)**0.5
            if distance <= proximity_threshold and component1['type'] == component2['type']:
                G.add_edge(component1['barcode'], component2['barcode'])

    connected_components = list(nx.connected_components(G))
    filtered_connected_components = [s for s in connected_components if len(s) > 1]
    return filtered_connected_components


def save_connected_components(file_path, components):
    with open(file_path, 'wb') as file:
        pickle.dump(components, file)


def main():
    path_dir = '/Users/victoriagao/local_docs/NEST/stored_variables'
    gene_ids_file = os.path.join(path_dir, 'PDAC_experiments/exp1_C1_gene_ids.txt')
    coordinates_file = os.path.join(path_dir, 'PDAC_experiments/exp1_C1_coordinates.npy')
    cell_barcode_file = os.path.join(path_dir, 'PDAC_experiments/exp1_C1_cell_barcode.pkl')
    pathologist_label_file = '/Users/victoriagao/local_docs/NEST/IX_annotation_artifacts_PDAC64630.csv'
    output_file = os.path.join(path_dir, 'exp1_C1_filtered_connected_components.pkl')

    gene_ids = load_gene_ids(gene_ids_file)
    coordinates = load_coordinates(coordinates_file)
    cell_barcodes = load_pickle(cell_barcode_file)
    pathologist_labels = load_pathologist_labels(pathologist_label_file)

    barcode_serial = create_barcode_serial(cell_barcodes)
    barcode_info = create_barcode_info(cell_barcodes, coordinates)

    barcode_type = {line[0]: line[1] for line in pathologist_labels[1:]}
    barcode_coordinate_dict = {item[0]: [item[1], item[2]] for item in barcode_info[1:]}
    combined_barcode = combine_barcodes(barcode_type, barcode_coordinate_dict)
    combined_barcode_list = create_combined_barcode_list(combined_barcode)

    # max_distance, max_distance_barcodes, min_distance, closest_barcodes = find_extreme_distances(combined_barcode_list)

    filtered_connected_components = find_connected_components(combined_barcode_list, proximity_threshold=200) # change proximity threshold if distance between spots are different

    for idx, component_set in enumerate(filtered_connected_components, start=1):
        component_types = {component['type'] for barcode in component_set for component in combined_barcode_list if component['barcode'] == barcode}
        print(f"Connected Component {idx} (Type: {component_types}): {component_set}")

    save_connected_components(output_file, filtered_connected_components)


if __name__ == "__main__":
    main()
