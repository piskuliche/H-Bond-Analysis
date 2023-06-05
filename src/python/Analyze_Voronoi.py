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
            if o ==0: 
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


def Analyze_Leaflets(mda_U, first_leaf=[], second_leaf=[],
                     selection="(resname POPC and name P*)",
                     extra_sel = None):
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
    def _get_leaflet(mda_U, selection, periodicity=True):
        #This function takes a MDAnalysis universe and returns the leaflets of the membrane.
        atom_selection=mda_U.select_atoms(selection)
        L = LeafletFinder(mda_U, atom_selection, pbc=periodicity,cutoff=12.6)
        return L.groups(0), L.groups(1)
    
    def _cluster_laur(mda_U, leafgroup1, leafgroup2, extra_sel):
        latom_select = mda_U.select_atoms(extra_sel)
        # Grab Z Positions
        z_pos = latom_select.positions[:,2]
        # Cluster using KMEANS algorithm
        labels = KMeans(n_clusters=2).fit(z_pos.reshape(-1,1)).labels_
        group1, group2 = labels==0, labels==1
        # Get the leaflet groups
        laur1, laur2 = latom_select[group1], latom_select[group2]
        # Compare the average z positions of the four groups
        z1, z2 = np.average(leafgroup1.positions[:,2]), np.average(leafgroup2.positions[:,2])
        lz1, lz2 = np.average(laur1.positions[:,2]), np.average(laur2.positions[:,2])
        if np.abs(z1-lz1) < np.abs(z2-lz1):
            leafgroup1 = leafgroup1.union(laur1)
            leafgroup2 = leafgroup2.union(laur2)
        else:
            leafgroup1 = leafgroup1.union(laur2)
            leafgroup2 = leafgroup2.union(laur1)
        return leafgroup1, leafgroup2
    

    from MDAnalysis.analysis.leaflet import LeafletFinder
    from sklearn.cluster import KMeans

    for ts in mda_U.trajectory:
        # Grabs the leaflets
        l1, l2 = _get_leaflet(mda_U, selection)
        # Adds Laurdan - if necessary
        if extra_sel is not None:
            l1, l2 = _cluster_laur(mda_U, l1, l2, extra_sel)
        # Appends the leaflets to the list
        first_leaf.append(l1); second_leaf.append(l2)
    return first_leaf, second_leaf
        

def Generate_Voronoi_Diagrams(mda_U, first_leaf, second_leaf, ps_selection=None, ifile=1):
    """
    """
    import freud

    def _grab_ps_points(mda_U, ps_selection, fbox):
        ps_atoms = mda_U.select_atoms(ps_selection)
        ps_points = ps_atoms.center_of_mass(wrap=True, compound='residues')
        ps_points[:,2] = 0.0
        ps_points = fbox.wrap(ps_points)
        return ps_points
    
    def _grab_leaf_points(leaf, fbox):
        points = leaf.positions
        points[:,2] = 0.0
        points = fbox.wrap(points)
        return points
    
    def _check_occupancy(points, points2, fbox):
        if points2 is None: return None
        dist = fbox.compute_all_distances(points2, points)
        cell_occupancy = np.zeros(np.shape(dist)[1])
        for d in dist:
            loc = np.where(d == np.min(d))[0][0]
            cell_occupancy[loc]=1
        return cell_occupancy
    
    def _compare_occupancy(occ1, occ2):
        if occ1 is None or occ2 is None:
            if occ1 is None:
                return occ2*2
            if occ2 is None:
                return occ1
        o1_nz = (occ1 != 0)*1
        o2_nz = (occ2 != 0)*2
        both = o1_nz+o2_nz
        return both

    def _check_array_periodicity(x, y, Lx, Ly):
        splits, dxs, dys = [], [], []
        for i in range(1,len(x)):
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
        for x, y in zip(xarr,yarr):
            splits, dx, dy =_check_array_periodicity(x, y, Lx, Ly)
            dx = dx * scale
            dy = dy * scale
            if len(splits) > 0:
                for i in range(len(splits)):
                    if i == 0:
                        outx[count] = x[:splits[i]]
                        outy[count] = y[:splits[i]]
                        outx[count]=np.append(outx[count],x[splits[i]-1]+dx[i])
                        outy[count]=np.append(outy[count],y[splits[i]-1]+dy[i])
                    else:
                        outx[count] = np.array(x[splits[i-1]-1]+dx[i-1])
                        outy[count] = np.array(y[splits[i-1]-1]+dy[i-1])
                        outx[count]=np.append(outx[count],x[splits[i-1]:splits[i]])
                        outy[count]=np.append(outy[count],y[splits[i-1]:splits[i]])
                        outx[count]=np.append(outx[count],x[splits[i]-1]+dx[i])
                        outy[count]=np.append(outy[count],y[splits[i]-1]+dy[i])
                    count += 1
                outx[count] = np.array(x[splits[-1]]-dx[-1])
                outy[count] = np.array(y[splits[-1]]-dy[-1])
                outx[count] = np.append(outx[count],x[splits[-1]:])
                outy[count] = np.append(outy[count],y[splits[-1]:])
                count += 1
            else:
                outx[count] = x
                outy[count] = y
                count += 1
        return outx, outy



    ps_xy = None
    nframes = len(mda_U.trajectory)
    # Loop over the trajectory
    for frame, ts in enumerate(mda_U.trajectory):
        frame_index = (ifile-1)*nframes + frame
        box=ts.dimensions
        # Create the box (freud)
        freud_box = freud.box.Box(Lx=box[0], Ly=box[1], is2D=True)
        # Grab PS points if they exist
        if ps_selection is not None: ps_xy = _grab_ps_points(mda_U, ps_selection, freud_box)
        lfdirs = ["voronoi_plots/first_leaf/","voronoi_plots/second_leaf/"]
        
        # Work on both leaflets
        for i, leaf in enumerate([first_leaf[frame_index], second_leaf[frame_index]]):
            # Grab the leaflet points
            lf_xy = _grab_leaf_points(leaf,freud_box)
            # Grab LAUR points
            laur_atoms = mda_U.select_atoms("resname LAUR and group lf", lf=leaf)
            laur_xy = _grab_leaf_points(laur_atoms,freud_box)
            # Check PS occupancy
            ps_occupancy = _check_occupancy(lf_xy, ps_xy, freud_box)
            la_occupancy = _check_occupancy(lf_xy, laur_xy, freud_box)
            occupancy = _compare_occupancy(ps_occupancy, la_occupancy)
            # Compute the diagram
            vor = freud.locality.Voronoi(freud_box, lf_xy)
            cells = vor.compute((freud_box,lf_xy))
            # Plot the diagram
            fig = plt.figure(figsize=(4,4), dpi=300)
            ax = plt.gca()
            if ps_xy is None:
                Voronoi_Plot(freud_box, cells.polytopes, occ=None, ax=ax)
            else:
                Voronoi_Plot(freud_box, cells.polytopes, occ=occupancy, ax=ax)
                ps_x=ps_xy[:,0]; ps_y=ps_xy[:,1]
                ax.scatter(ps_x, ps_y, s=20, c='blue')
                # Reshape the xs and ys
                xs = np.reshape(ps_x, (-1, 10)); ys = np.reshape(ps_y, (-1, 10))
                # Modify xs and ys to account for periodicity
                plt_x, plt_y = _reshape_2d_array_by_splits(xs, ys, box[0], box[1])
                for xi,yi in zip(plt_x, plt_y):
                    ax.plot(plt_x[xi],plt_y[yi], c='blue')
            # Plot the leaflet points
            ax.scatter(lf_xy[:,0],lf_xy[:,1], s=10, c='k')
            # Plot the laurdan points
            ax.scatter(laur_xy[:,0],laur_xy[:,1], s=15, c='red')
            # Set the axis limits
            #lower, upper = -box[0]/2-10, box[0]/2+10
            #lower, upper = np.round(lower/10)*10, np.round(upper/10)*10
            lower, upper = -50, 50
            ax.set_xlim((lower,upper))
            ax.set_ylim((lower,upper))
            # Set the axis labels
            ax.set_xlabel("x (nm)")
            ax.set_ylabel("y (nm)")
            # Set axis ticks
            ax.set_xticks(np.arange(lower,upper,10))
            ax.set_yticks(np.arange(lower,upper,10))
            # Set the axis tick labels
            ax.set_xticklabels(np.arange(lower,upper,10).astype(int))
            ax.set_yticklabels(np.arange(lower,upper,10).astype(int))
            # Save the figure
            plt.tight_layout()
            plt.savefig("%sframe_%d.png"%(lfdirs[i],frame_index), dpi=300)
            plt.close()
    return

def Do_Files(toploc="gro/", trjloc="xtc/", trjprefix='step7_', fstart=1, fstop=5,
             leafsel="(resname POPC and name P*)", laursel="(resname LAUR and name O*)"):
    """
    """
    Setup_Safe_Directory("voronoi_plots/")
    Setup_Safe_Directory("voronoi_plots/first_leaf/")
    Setup_Safe_Directory("voronoi_plots/second_leaf/")
    first_leaf, second_leaf = [], []
    # Loop Over Files
    for ifile in range(fstart, fstop+1):
        print("Doing file: ", ifile)
        topfile = toploc + trjprefix + str(ifile) + ".gro"
        trjfile = trjloc + trjprefix + str(ifile) + ".xtc"
        # Load Universe
        u = mda.Universe(topfile, trjfile)
        # Analyze Leaflets
        first_leaf, second_leaf = Analyze_Leaflets(u, first_leaf, second_leaf,
                                                    selection=leafsel, extra_sel=laursel)
        # Voronoi Tesselation
        Generate_Voronoi_Diagrams(u, first_leaf, second_leaf, ps_selection="resname STYRR and name C1", ifile=ifile)

        
    return

def Setup_Safe_Directory(dirname):
    # Makes an analysis directory if it doesn't exist
    import os
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    return

    

if __name__ == "__main__":
    Do_Files(leafsel="(resname POPC and name P*)")