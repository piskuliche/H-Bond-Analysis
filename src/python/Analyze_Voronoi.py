import numpy as np
import MDAnalysis as mda

def Voronoi_Plot(box, polytopes, occ=None, ax=None, color_by_sides=True, cmap=None):
    """
    This function is taken directly from the Freud github, and modified to change
    how coloring is done. 

    Helper function to draw 2D Voronoi diagram.
    Args:
        box (:class:`freud.box.Box`):
            Simulation box.
        polytopes (:class:`numpy.ndarray`):
            Array containing Voronoi polytope vertices.
        ax (:class:`matplotlib.axes.Axes`): Axes object to plot.
            If :code:`None`, make a new axes and figure object.
            (Default value = :code:`None`).
        color_by_sides (bool):
            If :code:`True`, color cells by the number of sides.
            If :code:`False`, random colors are used for each cell.
            (Default value = :code:`True`).
        cmap (str):
            Colormap name to use (Default value = :code:`None`).
    Returns:
        :class:`matplotlib.axes.Axes`: Axes object with the diagram.
    """
    from matplotlib import cm
    from matplotlib.collections import PatchCollection
    from matplotlib.colorbar import Colorbar
    from matplotlib.patches import Polygon
    from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

    if ax is None:
        fig = plt.figure()
        ax = fig.subplots()

    # Draw Voronoi polytopes
    patches = [Polygon(poly[:, :2]) for poly in polytopes]
    patch_collection = PatchCollection(patches, edgecolors="black", alpha=0.4)


    # Ensure we have enough colors to uniquely identify the cells
    colors = []
    for o in occ:
        if o ==0: 
            colors.append('white')
        else:
            colors.append('red')
        
    patch_collection.set_facecolor(colors)
    ax.add_collection(patch_collection)

    # Draw box
    corners = [[0, 0, 0], [0, 1, 0], [1, 1, 0], [1, 0, 0]]
    # Need to copy the last point so that the box is closed.
    corners.append(corners[0])
    corners = box.make_absolute(corners)[:, :2]
    ax.plot(corners[:, 0], corners[:, 1], color="k")

    # Set title, limits, aspect
    ax.set_title("Voronoi Diagram")
    ax.set_xlim((np.min(corners[:, 0]), np.max(corners[:, 0])))
    ax.set_ylim((np.min(corners[:, 1]), np.max(corners[:, 1])))
    ax.set_aspect("equal", "datalim")

    # Add colorbar for number of sides
    return ax

def Make_Voronoi_Plots(topfile, trjfile):
    """
    """
    

    

if __name__ == "__main__":
    # Do Something