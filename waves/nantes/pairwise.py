import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load your linked trajectory data into a pandas DataFrame
# Replace 'your_linked_data.csv' with your actual data file
linked_df = pd.read_csv('run2_small_track.csv')

# Function to calculate mean squared inter-particle distance at each time frame
def calculate_mean_squared_distance(dataframe):
    mean_squared_distances = []

    # Group the DataFrame by frame number
    grouped = dataframe.groupby('frame')

    for frame, group in grouped:
        x = group['x'].values
        y = group['y'].values
        N = len(x)

        if N > 1:
            dx = np.subtract.outer(x, x)
            dy = np.subtract.outer(y, y)
            squared_distance_matrix = dx**2 + dy**2
            np.fill_diagonal(squared_distance_matrix, 0)
            mean_squared_distance = np.mean(squared_distance_matrix)
            mean_squared_distances.append(mean_squared_distance)

    return mean_squared_distances

# Calculate mean squared inter-particle distances over time
mean_squared_distances = calculate_mean_squared_distance(linked_df)

# Create a time array corresponding to frame numbers
time = np.arange(len(mean_squared_distances))

# Plot the mean squared inter-particle distance as a function of time
plt.figure(figsize=(10, 6))
plt.plot(time, mean_squared_distances, marker='o', linestyle='-')
plt.xlabel('Time (Frame Number)')
plt.ylabel('Mean Squared Inter-particle Distance')
plt.title('Mean Squared Inter-particle Distance vs. Time')
plt.grid(True)
plt.show()
