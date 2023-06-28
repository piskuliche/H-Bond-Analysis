#!/usr/bin/env python
import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt


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
    patch_collection = PatchCollection(patches, edgecolors="black", alpha=0.6)

    # Ensure we have enough colors to uniquely identify the cells
    colors = []
    if occ is not None:
        for o in occ:
            if o == 0:
                colors.append('white')
            elif o == 1:
                colors.append('cyan')
            elif o == 2:
                colors.append('pink')
            else:
                colors.append('purple')

    patch_collection.set_facecolor(colors)
    ax.add_collection(patch_collection)

    # Draw box
    corners = [[0, 0, 0], [0, 1, 0], [1, 1, 0], [1, 0, 0]]
    # Need to copy the last point so that the box is closed.
    corners.append(corners[0])
    corners = box.make_absolute(corners)[:, :2]
    ax.plot(corners[:, 0], corners[:, 1], color="k")

    # Set title, limits, aspect
    ax.set_xlim((np.min(corners[:, 0]), np.max(corners[:, 0])))
    ax.set_ylim((np.min(corners[:, 1]), np.max(corners[:, 1])))
    ax.set_aspect("equal", "datalim")

    # Add colorbar for number of sides
    return ax


def Analyze_Leaflets(mda_U, selection="(resname POPC and name P*)",
                     extra_sel=None):
    """ This function takes a MDAnalysis universe and returns the leaflets of the membrane.

    Args:
        mda_U (MDAnalysis universe): The universe to analyze
        first_leaf (list): The first leaflet
        second_leaf (list): The second leaflet
        selection (list): The selections to use to define the leaflets
        extra_sel (list): Any extra selections to use for other leaflet components

    Returns:
        leaflets (list): A list of the leaflets
    """

    from MDAnalysis.analysis.leaflet import LeafletFinder
    from sklearn.cluster import KMeans

    def _get_leaflet(mda_U, selection, periodicity=True):
        # This function takes a MDAnalysis universe and returns the leaflets of the membrane.
        atom_selection = mda_U.select_atoms(selection)
        L = LeafletFinder(mda_U, atom_selection, pbc=periodicity, cutoff=12.6)
        return L.groups(0), L.groups(1)

    def _cluster_laur(mda_U, leafgroup1, leafgroup2, extra_sel):
        latom_select = mda_U.select_atoms(extra_sel)
        # Grab Z Positions
        z_pos = latom_select.positions[:, 2]
        # Cluster using KMEANS algorithm
        labels = KMeans(n_clusters=2, n_init=10).fit(
            z_pos.reshape(-1, 1)).labels_
        group1, group2 = labels == 0, labels == 1
        # Get the leaflet groups
        laur1, laur2 = latom_select[group1], latom_select[group2]
        # Compare the average z positions of the four groups
        z1, z2 = np.average(leafgroup1.positions[:, 2]), np.average(
            leafgroup2.positions[:, 2])
        lz1, lz2 = np.average(laur1.positions[:, 2]), np.average(
            laur2.positions[:, 2])
        if np.abs(z1-lz1) < np.abs(z2-lz1):
            leafgroup1 = leafgroup1.union(laur1)
            leafgroup2 = leafgroup2.union(laur2)
        else:
            leafgroup1 = leafgroup1.union(laur2)
            leafgroup2 = leafgroup2.union(laur1)
        return leafgroup1, leafgroup2

    first_leaf, second_leaf = [], []
    for ts in mda_U.trajectory:
        # Grabs the leaflets
        l1, l2 = _get_leaflet(mda_U, selection)
        # Adds Laurdan - if necessary
        if extra_sel is not None:
            l1, l2 = _cluster_laur(mda_U, l1, l2, extra_sel)
        # Appends the leaflets to the list
        first_leaf.append(l1)
        second_leaf.append(l2)
    return first_leaf, second_leaf


def Generate_Voronoi_Diagrams(mda_U, first_leaf, second_leaf, ps_selection=None, ifile=1, plot_every=50):
    """This takes the leaflets and generates Voronoi diagrams for each frame.

    Args:
        mda_U (MDAnalysis universe): The universe to analyze
        first_leaf (list): The first leaflet - must have same # frames
        second_leaf (list): The second leaflet - must have same # frames
        ps_selection (str): The selection for the PS points
        ifile (int): The file number being analyzed - for indexing the output dictionary
        plot_every (int): The number of frames to skip between plotting and saving voronoi diagrams

    Returns:
        voronoi_data (dict): A dictionary containing the voronoi diagrams for each frame

    """
    import freud

    def _grab_ps_points(mda_U, ps_selection, fbox):
        # Grabs the PS xy coordinates and sets z to 0
        ps_atoms = mda_U.select_atoms(ps_selection)
        ps_points = ps_atoms.center_of_mass(wrap=True, compound='residues')
        ps_points[:, 2] = 0.0
        ps_points = fbox.wrap(ps_points)
        return ps_points

    def _grab_leaf_points(leaf, fbox):
        # Grabs the leaflet points and sets z to 0
        points = leaf.positions
        points[:, 2] = 0.0
        points = fbox.wrap(points)
        return points

    def _check_occupancy(points, points2, fbox):
        # Checks the occupancy of the voronoi cells
        if points2 is None:
            return None
        dist = fbox.compute_all_distances(points2, points)
        cell_occupancy = np.zeros(np.shape(dist)[1])
        for d in dist:
            loc = np.where(d == np.min(d))[0][0]
            cell_occupancy[loc] = 1
        return cell_occupancy

    def _compare_occupancy(psoccupancy, resnames):
        #TODO Generalize this function
        # Determines the occupancy of the voronoi cells
        popc_molecs = resnames == 'POPC'
        laur_molecs = resnames == 'LAUR'
        occupancy = np.zeros(len(resnames))
        occupancy = occupancy + psoccupancy + laur_molecs*2
        return occupancy

    def _check_array_periodicity(x, y, Lx, Ly):
        # Checks the periodicity of the array
        # Corrects the periodicity if there is a problem (aka crosses the boundary)
        # This makes it so matplotlib can plot the line across the border rather than across the box.
        splits, dxs, dys = [], [], []
        for i in range(1, len(x)):
            dx = x[i]-x[i-1]
            dy = y[i]-y[i-1]
            if np.abs(dx) > Lx/2.0:
                if i not in splits:
                    splits.append(i)
                    dxs.append(dx-Lx*np.round(dx/Lx))
                    dys.append(dy-Ly*np.round(dy/Ly))
            if np.abs(dy) > Ly/2.0:
                if i not in splits:
                    splits.append(i)
                    dxs.append(dx-Lx*np.round(dx/Lx))
                    dys.append(dy-Ly*np.round(dy/Ly))
        return np.array(splits), np.array(dxs), np.array(dys)

    def _reshape_2d_array_by_splits(xarr, yarr, Lx, Ly, scale=1.0):
        outx, outy = {}, {}
        count = 0
        for x, y in zip(xarr, yarr):
            splits, dx, dy = _check_array_periodicity(x, y, Lx, Ly)
            dx = dx * scale
            dy = dy * scale
            if len(splits) > 0:
                for i in range(len(splits)):
                    if i == 0:
                        outx[count] = x[:splits[i]]
                        outy[count] = y[:splits[i]]
                        outx[count] = np.append(
                            outx[count], x[splits[i]-1]+dx[i])
                        outy[count] = np.append(
                            outy[count], y[splits[i]-1]+dy[i])
                    else:
                        outx[count] = np.array(x[splits[i-1]]-dx[i-1])
                        outy[count] = np.array(y[splits[i-1]]-dy[i-1])
                        outx[count] = np.append(
                            outx[count], x[splits[i-1]:splits[i]])
                        outy[count] = np.append(
                            outy[count], y[splits[i-1]:splits[i]])
                        outx[count] = np.append(
                            outx[count], x[splits[i]-1]+dx[i])
                        outy[count] = np.append(
                            outy[count], y[splits[i]-1]+dy[i])
                    count += 1
                outx[count] = np.array(x[splits[-1]]-dx[-1])
                outy[count] = np.array(y[splits[-1]]-dy[-1])
                outx[count] = np.append(outx[count], x[splits[-1]:])
                outy[count] = np.append(outy[count], y[splits[-1]:])
                count += 1
            else:
                outx[count] = x
                outy[count] = y
                count += 1
        return outx, outy

    def _plot_voronoi_diagram(filename, fbox, box, cells, occupancy, lf_xy,
                              laur_xy=None, ps_xy=None, upper=50, lower=-50,
                              figsize=(4, 4), dpi=300):
        """ Plots the voronoi diagram for a given frame

        Args:
            filename (str): The name of the file to save the plot to
            box (tuple): The box dimensions
            cells (freud voronoi object): The voronoi object for the frame
            occupancy (np.array): The occupancy of each cell
            lf_xy (np.array): The leaflet points
            laur_xy (np.array): The laurdan points
            ps_xy (np.array): The ps points
            upper (float): The upper dimension
            lower (float): The lower dimension
            **kwargs: Keyword arguments for the matplotlib figure

        Returns:
            None

        """
        fig = plt.figure(figsize=figsize, dpi=dpi)
        ax = plt.gca()
        if ps_xy is None:
            Voronoi_Plot(fbox, cells.polytopes, occ=None, ax=ax)
        else:
            Voronoi_Plot(fbox, cells.polytopes, occ=occupancy, ax=ax)
            e_x = ps_xy[:, 0]
            e_y = ps_xy[:, 1]
            ax.scatter(e_x, e_y, s=20, c='blue')
            xs = np.reshape(e_x, (-1, 10))
            ys = np.reshape(e_y, (-1, 10))
            # Modify xs and ys to account for periodicity
            plt_x, plt_y = _reshape_2d_array_by_splits(xs, ys, box[0], box[1])
            for xi, yi in zip(plt_x, plt_y):
                ax.plot(plt_x[xi], plt_y[yi], c='blue')
        # Plot the leaflet points
        ax.scatter(lf_xy[:, 0], lf_xy[:, 1], s=10, c='k')
        if laur_xy is not None:
            # Plot the laurdan points
            ax.scatter(laur_xy[:, 0], laur_xy[:, 1], s=15, c='red')
        # Set the axis limits
        # lower, upper = -box[0]/2-10, box[0]/2+10
        # lower, upper = np.round(lower/10)*10, np.round(upper/10)*10
        ax.set_xlim((lower, upper))
        ax.set_ylim((lower, upper))
        # Set the axis labels
        ax.set_xlabel("x (nm)")
        ax.set_ylabel("y (nm)")
        # Set axis ticks
        ax.set_xticks(np.arange(lower, upper, 10))
        ax.set_yticks(np.arange(lower, upper, 10))
        # Set the axis tick labels
        ax.set_xticklabels(np.arange(lower, upper, 10).astype(int))
        ax.set_yticklabels(np.arange(lower, upper, 10).astype(int))
        # Save the figure
        plt.tight_layout()
        plt.savefig("%s" % (filename), dpi=dpi)
        plt.close()
        return

    def _CHECK_OVERLAPS(r, iter=0, frame_index=0):
        # Check for overlaps and correct them
        import sys
        overlap_found = False
        for j, coord1 in enumerate(r):
            for k, coord2 in enumerate(r):
                if j != k:
                    a = (coord1[0], coord1[1])
                    b = (coord2[0], coord2[1])
                    if a == b:
                        print("Overlap found on frame %d" % (frame_index))
                        print("Correcting...")
                        print("Note: This is a hacky fix, but it works for now")
                        r[k] = [coord2[0]+0.0001, coord2[1]+0.0001, 0]
        if overlap_found == True:
            r = _CHECK_OVERLAPS(r, iter=iter+1, frame_index=frame_index)
        if iter > 50:
            print("Too many overlaps found. Exiting...")
            sys.exit()
        return r

    voronoi_data = {}
    ps_xy = None
    nframes = len(mda_U.trajectory)
    # Loop over the trajectory
    for frame, ts in enumerate(mda_U.trajectory):
        frame_index = (ifile-1)*nframes + frame
        box = ts.dimensions
        # Create the box (freud)
        freud_box = freud.box.Box(Lx=box[0], Ly=box[1], is2D=True)
        # Grab PS points if they exist
        if ps_selection is not None:
            ps_xy = _grab_ps_points(mda_U, ps_selection, freud_box)
        lfdirs = ["voronoi_plots/first_leaf/", "voronoi_plots/second_leaf/"]

        voronoi_data[frame_index] = {0: {}, 1: {}}
        # Work on both leaflets
        for i, leaf in enumerate([first_leaf[frame], second_leaf[frame]]):
            # Grab the leaflet points
            lf_xy = _grab_leaf_points(leaf, freud_box)
            # Grab LAUR points
            laur_atoms = mda_U.select_atoms(
                "resname LAUR and group lf", lf=leaf)
            laur_xy = _grab_leaf_points(laur_atoms, freud_box)
            # Check PS occupancy
            ps_occupancy = _check_occupancy(lf_xy, ps_xy, freud_box)
            occupancy = _compare_occupancy(ps_occupancy, leaf.resnames)
            """
            This was used to test the occupancy issues I was getting
            I think the new function is performing better.
            if frame_index == 22357:
                print(ps_occupancy, la_occupancy)
                print(len(ps_occupancy), np.shape(lf_xy))
                print(occupancy)
            """
            lf_xy = _CHECK_OVERLAPS(lf_xy, frame_index=frame_index)
            # Compute the diagram
            vor = freud.locality.Voronoi(freud_box, lf_xy)
            cells = vor.compute((freud_box, lf_xy))
            areas = cells.volumes
            voronoi_data[frame_index][i]['areas'] = areas
            voronoi_data[frame_index][i]['occupancy'] = occupancy
            voronoi_data[frame_index][i]['atomindex'] = leaf.indices
            voronoi_data[frame_index][i]['resnames'] = leaf.resnames
            voronoi_data[frame_index][i]['molids'] = leaf.resids


            # Plot the diagram on the plot_everyth frame
            if frame_index % plot_every == 0:
                pngname = "%sframe_%d.png" % (lfdirs[i], frame_index)
                _plot_voronoi_diagram(pngname, freud_box, box, cells, occupancy, lf_xy,
                                      laur_xy, ps_xy, upper=50, lower=-50, figsize=(4, 4), dpi=300)

    return voronoi_data


def Do_Files(toploc="gro/", trjloc="xtc/", trjprefix='step7_', fstart=1, fstop=10,
             leafsel="(resname POPC and name P*)", laursel="(resname LAUR and name O*)"):
    """ This function loops over the files and does the voronoi/leaflet analysis

    Args:
        toploc (str): The location of the topology files
        trjloc (str): The location of the trajectory files
        trjprefix (str): The prefix of the trajectory files
        fstart (int): The starting file number
        fstop (int): The ending file number
        leafsel (str): The selection string for the leaflet atoms
        laursel (str): The selection string for the laurdan atoms

    Returns:
        None

    Generates:
        voronoi_plots/first_leaf/frame_*.png
        voronoi_plots/second_leaf/frame_*.png
        voronoi_plots/data/vor_*.pkl

    """
    import pickle
    import traceback
    Setup_Safe_Directory("voronoi_plots/")
    Setup_Safe_Directory("voronoi_plots/data/")
    Setup_Safe_Directory("voronoi_plots/first_leaf/")
    Setup_Safe_Directory("voronoi_plots/second_leaf/")
    # Loop Over Files
    for ifile in range(fstart, fstop+1):
        print("Doing file: ", ifile)
        topfile = toploc + trjprefix + str(ifile) + ".gro"
        trjfile = trjloc + trjprefix + str(ifile) + ".xtc"
        try:
            # Load Universe
            u = mda.Universe(topfile, trjfile)
            # Analyze Leaflets
            first_leaf, second_leaf = Analyze_Leaflets(
                u, selection=leafsel, extra_sel=laursel)
            # Voronoi Tesselation
            voronoi_data = Generate_Voronoi_Diagrams(
                u, first_leaf, second_leaf, ps_selection="resname STYRR and name C1", ifile=ifile)
            pickle.dump(voronoi_data, open(
                "voronoi_plots/data/vor_%d.pckl" % ifile, 'wb'))
        except Exception as e:
            print("Failed to do file: ", ifile)
            traceback.print_exc()
            continue

    return


def Setup_Safe_Directory(dirname):
    # Makes an analysis directory if it doesn't exist
    import os
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    return


if __name__ == "__main__":
    Do_Files(leafsel="(resname POPC and name P*)", fstart=1, fstop=100)
