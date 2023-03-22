import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.leaflet import LeafletFinder
from MDAnalysis.lib.distances import calc_bonds
from MDAnalysis.lib.NeighborSearch import AtomNeighborSearch as Neighbors

from numpy.linalg import norm
from mpl_toolkits import mplot3d

import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
import freud
import os

def User_Input():
    import argparse
    parser = argparse.ArgumentParser(description='Code to check Voronoi occupancy')
    parser.add_argument('-top', default = "../step8_5.tpr", type=str, help = 'Topology file name')
    parser.add_argument('-trj', default = "../step8_5.xtc", type=str, help = 'Trajectory file name')
    input_args = parser.parse_args()
    
    return input_args

def voronoi_plot(box, polytopes, occ=None, ax=None, color_by_sides=True, cmap=None):
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


if __name__ == "__main__":
    
    input_args = User_Input()
    u = mda.Universe(input_args.top,input_args.trj)
    patoms = u.select_atoms('name P*')
    watoms = u.select_atoms('name OH2')
    ps     = u.select_atoms('resname *STYRR')
    print(len(np.unique(ps.resids)))
    L = LeafletFinder(u, patoms,pbc=True)
    leafs = [L.groups(0),L.groups(1)]

    try:
        os.mkdir("leaf1_vor")
    except:
        pass
    try:
        os.mkdir("leaf2_vor")
    except:
        pass
    
    count=0
    for ts in u.trajectory[::100]:
        wpoints = watoms.positions
        box=ts.dimensions
        fbox = freud.box.Box(Lx=box[0],Ly=box[1],is2D=True)
        ps_points = ps.center_of_mass(pbc=True, compound='residues')
        ps_points[:,2] = 0.0
        ps_points = fbox.wrap(ps_points)
        print(np.shape(ps_points))
        lfs=["leaf1_vor","leaf2_vor"]
        for i,leaf in enumerate(leafs):
            points = leaf.positions
            points[:,2] = 0.0
            points = fbox.wrap(points)
            dist = fbox.compute_all_distances(ps_points,points)
            cell_occ = np.zeros(np.shape(dist)[1])
            for d in dist:
                loc=np.where(d==np.min(d))[0][0]
                cell_occ[loc]=1
            #vor = Voronoi(points)
            vor = freud.locality.Voronoi()
            cells = vor.compute((fbox,points))
            fig = plt.figure()
            ax = plt.gca()
            #vor.plot(ax=ax)
            voronoi_plot(fbox,cells.polytopes,occ=cell_occ,ax=ax)
            ax.set_xlim((-box[0]/2-10,box[0]/2+10))
            ax.set_ylim((-box[0]/2-10,box[0]/2+10))
            ax.scatter(points[:,0],points[:,1], s=10, c='k')
            ax.scatter(ps_points[:,0],ps_points[:,1],s=20, c='blue')
            plt.savefig("%s/vorplot_%d.png"%(lfs[i],count))
            plt.close()
        count += 1




        


