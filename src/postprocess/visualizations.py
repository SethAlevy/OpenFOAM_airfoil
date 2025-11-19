import matplotlib.pyplot as plt
import numpy as np


def plot_airfoil(upper_coords, lower_coords, mean_camber_coords, thickness_coords, title="NACA Airfoil", save_path=None):
    """Plot the airfoil geometry including upper surface, lower surface, and mean camber line."""
    plt.figure(figsize=(10, 5))
    plt.plot(upper_coords[0], upper_coords[1], label='Upper Surface', color='b')
    plt.plot(lower_coords[0], lower_coords[1], label='Lower Surface', color='r')
    plt.plot(mean_camber_coords[0], mean_camber_coords[1], label='Mean Camber Line', color='g', linestyle='--')
    plt.plot(thickness_coords[0], thickness_coords[1], label='Thickness Distribution', color='m', linestyle=':')
    plt.title(title)
    plt.xlabel('Chordwise Location')
    plt.ylabel('Thickness / Camber')
    plt.axis('equal')
    plt.grid(True)
    plt.legend()
    if save_path:
        plt.savefig(save_path)
    plt.show()