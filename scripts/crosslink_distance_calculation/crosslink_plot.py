import matplotlib.pyplot as plt
import os
import sys

def plot_crosslink_distances(file_path, label):
    distances = []
    with open(file_path, 'r') as file:
        for line in file:
            distance = float(line.split(':')[1].strip())
            distances.append(distance)

    plt.hist(distances, bins=20, alpha=0.5, label=label)

    plt.title('Crosslink Distances')
    plt.xlabel('Distance (Å)')
    plt.ylabel('Frame')
    plt.legend()

# input is txt file with crosslink distances. eg: 1clv_DSSO_2_distances.txt
file1_path = os.path.join("/home/muskaan/easal_imp/crosslink_distances/", sys.argv[1])
file2_path = os.path.join("/home/muskaan/easal_output_muskaan/crosslink_distances/", sys.argv[1])

plot_crosslink_distances(file1_path, 'IMP')
plot_crosslink_distances(file2_path, 'EASAL')

# Save the figure
plt.savefig(f'crosslink_distances_plot_{sys.argv[1]}.png')
plt.show()
