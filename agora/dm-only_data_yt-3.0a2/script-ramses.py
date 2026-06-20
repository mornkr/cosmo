#####################################################################
#
#  YT SCRIPT TO PRODUCE PROJECTED DARK MATTER DENSITY 
#  USING ONLY THE FINEST DM PARTICLES ON YT-3.0 alpha release 2
#  
#  PLEASE SEE:  https://hub.yt-project.org/nb/5zw7qn
#
#  FUNCTIONALITY AVAILABLE GREATLY THANKS TO:  Matthew Turk et al.
#
#  SCRIPT WRITTEN BY:  Ji-hoon Kim on June 6, 2013
#
#####################################################################

import sys
sys.path.insert(0, "/home/mornkr/yt-3.0")
from yt.config import ytcfg; ytcfg["yt","loglevel"] = "20"
from yt.mods import *
import yt.utilities.lib as au
import numpy as np
import copy
import pylab
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from yt.analysis_modules.halo_finding.api import *

def finest_DM_density(field, data):
    filter = np.where(data["ParticleMassMsun"] <= 340000)

    pos = np.column_stack([data["particle_position_%s" % ax][filter] for ax in 'xyz'])
    d = data.deposit(pos, [data["ParticleMass"][filter]], method = "cic")
    d /= data["CellVolume"]
    return d

RAMSESFieldInfo.add_field(("deposit", "finest_DM_density"),
         function = finest_DM_density,
         validators = [ValidateSpatial()],
         display_name = "\\mathrm{Finest DM Density}",
         units = r"\mathrm{g}/\mathrm{cm}^{3}",
         projected_units = r"\mathrm{g}/\mathrm{cm}^{2}",
         projection_conversion = 'cm')


center = np.array([0.486123, 0.526438, 0.490137])
ds = load("output_00101/info_00101.txt")
print ds.parameters["nstep_coarse"], ds.parameters["levelmax"]
print ds.units["mpchcm"]
print ds.units["mpch"]
print ds.units["cm"]

sp = ds.h.sphere(center, (1.0, 'mpc'))
total_particle_mass = sp.quantities["TotalQuantity"]( [("all","ParticleMassMsun")] )[0]
print "Total particle mass within a radius of 1 Mpc of the center: %0.3e Msun" % total_particle_mass
print ds.h.derived_field_list

pw = ProjectionPlot(ds, "z", ("deposit", "all_density"), weight_field=None, center=center, width=(1.0,'mpch')).save()

w = (6.0, "mpch")
source = ds.h.region(center, center - (w[0]/ds[w[1]])/2.0, center + (w[0]/ds[w[1]])/2.0)
proj = ds.h.proj( ("deposit", "all_density"), 2, weight_field = ("deposit", "all_density"), data_source = source)
pw = proj.to_pw(fields = [("deposit", "all_density")], center = center, width = w)
pw.set_zlim(("deposit","all_density"), 1e-32, 1e-25)
pw.save("info_00101_Projection_z_all_density_subset.png")

proj = ds.h.proj( ("deposit", "finest_DM_density"), 2, weight_field = ("deposit", "finest_DM_density"), data_source = source)
pw = proj.to_pw(fields = [("deposit", "finest_DM_density")], center = center, width = w)
pw.set_zlim(("deposit","finest_DM_density"), 1e-32, 1e-25)
pw.save("info_00101_Projection_z_finest_DM_density_subset.png")

