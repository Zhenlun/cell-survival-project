import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import random
import networkx as nx

class Voxel:
    def __init__(self, id, x, y, z):
        self.id = id
        self.x = x
        self.y = y
        self.z = z
        self.dna_id = None
        self.segment_recorded = 0
        self.arm_label = None

class DNA:
    def __init__(self, id, length, centromere_location):
        self.id = id
        self.length = length
        self.centromere_location = centromere_location
        self.current_location = None

def create_sphere(radius):
    sphere = []
    id = 1
    for x in range(-radius, radius+1):
        for y in range(-radius, radius+1):
            for z in range(-radius, radius+1):
                if x**2 + y**2 + z**2 <= radius**2:
                    sphere.append(Voxel(id, x, y, z))
                    id += 1
    return sphere

def create_dna_chains(num_chains, length, centromere_location):
    dna_chains = []
    for i in range(1, num_chains+1):
        dna_chains.append(DNA(i, length, centromere_location))
    return dna_chains

def place_dna_in_sphere(sphere, dna_chains):
    # Create a graph where each voxel is a node and each edge represents a neighboring relationship
    graph = nx.Graph()
    for voxel in sphere:
        graph.add_node(voxel.id, voxel=voxel)
    for i in range(len(sphere)):
        for j in range(i+1, len(sphere)):
            if is_neighbor(sphere[i], sphere[j]):
                graph.add_edge(sphere[i].id, sphere[j].id)

    # Randomly shuffle the DNA chains
    random.shuffle(dna_chains)

    # Assign each DNA chain to a path in the graph
    for dna in dna_chains:
        # Find a random start node that has not been assigned a DNA id
        start = random.choice([node for node in graph.nodes if graph.nodes[node]['voxel'].dna_id is None])
        # Use DFS to find a path of length equal to the DNA length
        path = dfs(graph, start, dna.length)
        if path is None:
            continue  # If no path of sufficient length is found, try with the next DNA chain

        # Assign the voxels in the path to the DNA chain
        for i, node in enumerate(path):
            voxel = graph.nodes[node]['voxel']
            voxel.dna_id = dna.id
            voxel.segment_recorded = i + 1
            voxel.arm_label = 0 if i < dna.centromere_location else 1

        # Remove the used voxels from the graph
        for node in path:
            graph.remove_node(node)

def is_neighbor(voxel1, voxel2):
    # Check if two voxels are neighbors (i.e., their coordinates differ by at most 1 in each dimension)
    return max(abs(voxel1.x - voxel2.x), abs(voxel1.y - voxel2.y), abs(voxel1.z - voxel2.z)) == 1

def dfs(graph, start, length):
    # Use DFS to find a path of a given length in the graph
    stack = [(start, [start])]
    while stack:
        (node, path) = stack.pop()
        for next in set(graph.neighbors(node)) - set(path):
            if len(path) + 1 == length:
                return path + [next]
            else:
                stack.append((next, path + [next]))
    return None

def print_voxels(sphere):
    for voxel in sphere:
        id = voxel.id if voxel.id is not None else 0
        dna_id = voxel.dna_id if voxel.dna_id is not None else 0
        segment_recorded = voxel.segment_recorded if voxel.segment_recorded is not None else 0
        arm_label = voxel.arm_label if voxel.arm_label is not None else 0
        #print(f'{id}'+ ' ' + f'{dna_id}' + ' ' + f'{segment_recorded}' + ' ' + f'{arm_label}')


def visualize_dna_voxels(sphere, dna_chains, view='global', cross_section_z=0, angle=30):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Create a colormap with a unique color for each DNA id
    colormap = cm.get_cmap('rainbow', len(dna_chains))

    for voxel in sphere:
        if voxel.dna_id is not None:
            if view == 'cross_section' and voxel.z != cross_section_z:
                continue  # Skip this voxel if we're in cross-section view and the voxel is not at the cross-section
            # Use the DNA id to get a unique color from the colormap
            color = colormap(voxel.dna_id)
            ax.scatter(voxel.x, voxel.y, voxel.z, c=color)

    # Set the viewing angle
    ax.view_init(elev=10., azim=angle)

    plt.show()

def check_all_voxels_used(sphere):
    for voxel in sphere:
        if voxel.dna_id is None:
            return False
    return True

radius = 10
cube_side_length = 1
num_chains = 46
dna_length = 100
centromere_location = 60

sphere = create_sphere(radius)
dna_chains = create_dna_chains(num_chains, dna_length, centromere_location)
place_dna_in_sphere(sphere, dna_chains)

vid = []
dnaid = []
for voxel in sphere:
    vid.append(voxel.id)
    dnaid.append(voxel.dna_id)
    #print(f'Voxel ID: {voxel.id}, DNA ID: {voxel.dna_id}')

print_voxels(sphere)
visualize_dna_voxels(sphere, dna_chains)
print("All voxels used:", check_all_voxels_used(sphere))
