import matplotlib.pyplot as plt
import numpy as np


def plot_airfoil(
        upper_coords: np.ndarray,
        lower_coords: np.ndarray,
        title: str = "Airfoil",
        save_path: str = None,
        show: bool = True,
        mean_camber_coords: np.ndarray = None,
        thickness_coords: np.ndarray = None,
        chord: float = None,
        alpha: float = 0,
) -> None:
    """
    Plot the airfoil geometry including upper surface, lower surface, mean camber
    line, thickness distribution, and chord.
    
    Args:
        upper_coords (np.ndarray): Coordinates of the upper surface.
        lower_coords (np.ndarray): Coordinates of the lower surface.
        title (str): Title of the plot.
        save_path (str): Path to save the plot image. If None, the plot is shown.
        show (bool): Whether to display the plot.
        mean_camber_coords (np.ndarray): Coordinates of the mean camber line.
        thickness_coords (np.ndarray): Coordinates of the thickness distribution.
        chord (float): Length of the chord line.
        alpha (float): Angle of attack in degrees. Required only to plot the chord line.
    """
    plt.figure(figsize=(10, 5))
    plt.plot(
        upper_coords[0],
        upper_coords[1],
        label='Upper Surface',
        color='b'
    )
    plt.plot(
        lower_coords[0],
        lower_coords[1],
        label='Lower Surface',
        color='r'
    )
    if mean_camber_coords is not None:
        plt.plot(
            mean_camber_coords[0],
            mean_camber_coords[1],
            label='Mean Camber Line',
            color='g',
            linestyle='--'
        )
    if thickness_coords is not None:
        plt.plot(
            thickness_coords[0],
            thickness_coords[1],
            label='Thickness Distribution',
            color='m',
            linestyle=':'
        )
    if chord is not None:
        plt.plot(
            [0, np.cos(np.radians(alpha)) * chord],
            [0, np.sin(np.radians(alpha)) * chord],
            label='Chord Line',
            color='k',
            linestyle='-'
        )

    plt.title(title)
    plt.xlabel('Chordwise Location')
    plt.ylabel('Thickness / Camber')
    plt.axis('equal')
    plt.grid(True)
    plt.legend()
    if save_path:
        plt.savefig(save_path) 
    if show:
        plt.show()
