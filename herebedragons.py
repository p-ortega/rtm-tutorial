import os
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import flopy
import pyemu
import geopandas as gpd
from mf6rtm import mup3d, utils
from flopy.utils.gridintersect import GridIntersect
from shapely.geometry import LineString
from pypestutils.pestutilslib import PestUtilsLib
lib = PestUtilsLib()

def get_wel_coords(gwf, name  = "wellin"):
    mg = gwf.modelgrid
    ix = GridIntersect(mg)
    wells = pd.read_csv(os.path.join("data", "wells.csv"))
    wells = gpd.GeoDataFrame(wells, geometry=gpd.points_from_xy(wells.x, wells.y))

    assert name in wells.name.values, f"{name} not in well"
    wells = wells[wells.name == name]
    geom = wells.geometry[wells.name==name].values

    assert len(geom)==1, f"more than one well with name {name} in wells.csv"
    cellid = ix.intersect(geom[0], 'point').cellids
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    mg.plot(ax=ax)
    ix.plot_point(ix.intersect(geom[0], 'point'), ax=ax)
    print(cellid[0])
    return cellid[0]

def make_stress_period_data(coords, rates, add_conc=True):
        if len(coords) != len(rates):
            raise ValueError("Coordinate and rate lists must be the same length")
        return [
            ([cell, q] if add_conc else [cell, q])
            for cell, q in zip(coords, rates)
        ]
    
def make_wel_out(gwf, conservative_tracer = None, mup3d_m=None, wellname = "welout"):

    nper = gwf.tdis.nper.get_data()
    layers = [1,3,5]
    cellid = get_wel_coords(gwf, name  = wellname)
    coords_out = [(lay, cellid[0], cellid[1]) for lay in layers]
    print(coords_out)

    init_rates_out  = [-300,  -30,  -30]                # 3 negatives
    fini_rates_out  = [-400,  -40,  -40]

    init_sp = range(0, 35)   # stress periods 0 – 35
    fini_sp = range(35, nper)  # stress periods 36 – 38
    all_sp  = (*init_sp, *fini_sp)

    # Time‑invariant blocks for each phase
    wellout_init = {sp: make_stress_period_data(coords_out, init_rates_out) for sp in all_sp}
    wellout_fini = {sp: make_stress_period_data(coords_out, fini_rates_out, add_conc=True) 
                    for sp in all_sp}
    wellout_sp_data = {sp: (wellout_init[sp] if sp in init_sp else wellout_fini[sp])
                for sp in all_sp}

    if conservative_tracer is not None:
        wellout_sp_data = append_values_to_inner_lists(wellout_sp_data, 0.0)
        wel_out = flopy.mf6.ModflowGwfwel(gwf, 
                                            stress_period_data=wellout_sp_data, 
                                            auxiliary=conservative_tracer,
                                            pname = 'welout' ,
                                            filename=f'{gwf.name}.welout')
        wel_out.set_all_data_external()
    else:
        wellout_sp_data = append_values_to_inner_lists(wellout_sp_data, [0.0]*len(mup3d_m.components))
        wel_out = flopy.mf6.ModflowGwfwel(gwf, 
                                            stress_period_data=wellout_sp_data, 
                                            auxiliary=mup3d_m.components,
                                            pname = 'welout',
                                            filename=f'{gwf.name}.welout')
        
        wel_out.set_all_data_external()
    return wel_out

def make_chd(gwf, conservative_tracer = None, mup3d_m=None):
    
    l_hd= 0
    domain = gpd.read_file(Path('data', 'domain.gpkg'))
    geom = domain.dissolve().geometry[0]
    minx, miny, maxx, maxy = geom.bounds
    left_boundary = LineString([(minx, miny), (minx+0.1, maxy)])
    right_boundary = LineString([(maxx, miny), (maxx-0.1, maxy)])

    ix = GridIntersect(gwf.modelgrid)
    left_cells = ix.intersect(left_boundary, 'line').cellids.tolist()
    right_cells = ix.intersect(right_boundary, 'line').cellids.tolist()
    left_cells.extend(right_cells)
    boundary_cells = left_cells

    nlay = gwf.dis.nlay.get_data()
    ncpl = gwf.dis.ncpl.get_data()

    if conservative_tracer is not None:
        df_inj = pd.read_csv(os.path.join("data", "ic_aq_chem.csv"), index_col=0)
        assert conservative_tracer in df_inj.index, f"compound {conservative_tracer} not in ic_aq_chem csv"
        c_list = [df_inj.loc[conservative_tracer, 'value']]
        aux = conservative_tracer
        # print(c_list)
    else:
        chdchem = mup3d.ChemStress('chdchem')
        sol_spd = [1]
        chdchem.set_spd(sol_spd)
        mup3d_m.set_chem_stress(chdchem)
        c_list = mup3d_m.chdchem.data[0]
        aux=mup3d_m.components

    chdspd = []
    for i in range(nlay):          # layers
        for icpl in boundary_cells:      # rows
            chdspd.append([(i, icpl), l_hd])           # left boundary
            # chdspd.append([(i, j, ncol - 1), l_hd])    # right boundary
    # print(chdspd)
    for i in range(len(chdspd)):
        chdspd[i].extend(c_list)

    left_line = LineString([(minx, miny), (minx, maxy)])
    right_line = LineString([(maxx, miny), (maxx, maxy)])

    # print(chdspd)
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        maxbound=len(chdspd),
        stress_period_data=chdspd,
        save_flows=True,
        auxiliary=aux,
        pname="CHD",
        filename=f"{gwf.name}.chd",
    )
    chd.set_all_data_external()
    return chd

def get_avg_distance(points, npoints=10):
    """
    Calculate the average distance to the nearest n points for each point in a set of points.

    Parameters
    ----------
    points : numpy array
        Array of points.
    npoints : int
        Number of nearest points to calculate the average distance to.
    Returns
    -------
    average_distances : numpy array
        Array of average distances to the nearest n points for each point in the input array.
    """

    from scipy.spatial import distance
    distances = distance.cdist(points, points, 'euclidean')
    np.fill_diagonal(distances, np.inf)
    nearest_n = np.partition(distances, npoints, axis=1)[:, :npoints]
    average_distances = np.mean(nearest_n, axis=1)
    return average_distances

def get_botms(gwf, ws):
    """
    Get botms using kriging from borehole points.
    
    Parameters
    ----------
    gwf : flopy.mf6.ModflowGwf
        The groundwater flow model object.
    ws : str
        The workspace directory where temporary files will be stored.
    Returns
    -------
    botms : list
        List of bottom elevations for each layer.
    """

    bps = gpd.read_file(os.path.join("data", 'botm.gpkg'))
    bps['x'] = bps.geometry.centroid.x
    bps['y'] = bps.geometry.centroid.y
    bps

    ppeasting = bps.x.values
    ppnorthing = bps.y.values
    anis = 1
    bearing= 0.0
    aa = 1.5 * get_avg_distance(bps[['x','y']].values, 2).max()

    ib = gwf.dis.idomain.get_data()
    # cellids = df.loc[df.layer==layer+1].icpl.values - 1 # zero-based
    easting = gwf.modelgrid.xcellcenters.flatten()
    northing = gwf.modelgrid.ycellcenters.flatten()

    max_pts = 50 # pp are same as cell centers, so kind of irrelevant
    min_pts = 1
    search_dist = 1.e+10
    aa_pp = aa #?
    zone_pp = np.ones_like(ppeasting,dtype=int)
    fac_file = os.path.join(ws,f"factors.bin")

    ib = np.ones_like(easting,dtype=int)
    ipts = lib.calc_kriging_factors_2d(ppeasting,
                                    ppnorthing,
                                    zone_pp,
                                    easting,
                                    northing,
                                    ib.flatten(),
                                    "exp","ordinary",
                                    aa_pp,anis,bearing,search_dist,max_pts,min_pts,fac_file,"binary")

    botms = []
    icpls = gwf.dis.ncpl.get_data()
    for layer in range(1, 13):
        # get COND multiplier
        ppval = bps[f"botm_{layer}"].values
        result = lib.krige_using_file(os.path.join(ws,f"factors.bin"),
                                        "binary",
                                        icpls,
                                        "ordinary",
                                        "none",
                                        np.array(ppval),
                                        np.zeros_like(icpls),
                                        0)
        botms.append(np.round(result['targval'], 1))
    return botms