#####################################################################
#
#  YT SCRIPT TO PRODUCE PROJECTED DARK MATTER DENSITY 
#  USING ONLY THE FINEST DM PARTICLES ON YT-3.0 alpha release 2
#  
#  PLEASE SEE:  https://hub.yt-project.org/nb/tzeizy
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
    filter = np.where((data["particle_type"] == 1) & (data["ParticleMassMsun"] <= 340000))

    pos = np.column_stack([data["particle_position_%s" % ax][filter] for ax in 'xyz'])
    d = data.deposit(pos, [data["ParticleMass"][filter]], method = "cic")
    d /= data["CellVolume"]
    return d

EnzoFieldInfo.add_field(("deposit", "finest_DM_density"),
                        function = finest_DM_density,
                        validators = [ValidateSpatial()],
                        display_name = "\\mathrm{Finest DM Density}",
                        units = r"\mathrm{g}/\mathrm{cm}^{3}",
                        projected_units = r"\mathrm{g}/\mathrm{cm}^{2}",
                        projection_conversion = 'cm')

center = np.array([0.493002, 0.507927, 0.507289])
ds = load("DD0040/data0040")
print ds.parameters["HydroMethod"], ds.parameters["TimeLastDataDump"]
print ds.units["mpchcm"]
print ds.units["mpch"]
print ds.units["cm"]

sp = ds.h.sphere(center, (1.0, 'mpc'))
total_particle_mass = sp.quantities["TotalQuantity"]( [("all","ParticleMassMsun")] )[0]
#total_particle_mass = sp.quantities["TotalQuantity"]("ParticleMassMsun")[0]
print "Total particle mass within a radius of 1 Mpc of the center: %0.3e Msun" % total_particle_mass
print ds.h.derived_field_list

#pw = ProjectionPlot(ds, "x", ("deposit", "all_density"), weight_field=None, center=center)
#pw.zoom(1.01) # Some edge effects because periodicity isn't quite as expected yet
#pw.show()
pw = ProjectionPlot(ds, "z", ("deposit", "all_density"), weight_field=None, center=center, width=(1.0,'mpch')).save()
#pw = ProjectionPlot(ds, "z", ("deposit", "finest_DM_density"), weight_field=None, center=center, width=1./60.).save()
#pw = ProjectionPlot(ds, "z", "particle_density", weight_field=None, center=center, width=1./60.).save()

w = (6.0, "mpch")
source = ds.h.region(center, center - (w[0]/ds[w[1]])/2.0, center + (w[0]/ds[w[1]])/2.0)
proj = ds.h.proj( ("deposit", "all_density"), 2, weight_field = ("deposit", "all_density"), data_source = source)
pw = proj.to_pw(fields = [("deposit", "all_density")], center = center, width = w)
pw.set_zlim(("deposit","all_density"), 1e-32, 1e-25)
#pw.show()
pw.save("data0040_Projection_z_all_density_subset.png")

proj = ds.h.proj( ("deposit", "finest_DM_density"), 2, weight_field = ("deposit", "finest_DM_density"), data_source = source)
pw = proj.to_pw(fields = [("deposit", "finest_DM_density")], center = center, width = w)
pw.set_zlim(("deposit","finest_DM_density"), 1e-32, 1e-25)
pw.save("data0040_Projection_z_finest_DM_density_subset.png")

# pc = PlotCollection(ds, center=center)
# proj = pc.add_projection("particle_density", 2, center=center, weight_field=None, data_source = source)
# proj.set_width(w[0]/ds[w[1]], '1')
# proj.set_autoscale(False)
# proj.set_zlim(1e-6, 1e-1)
# pc.save("data0040_Projection_z_particle_density_subset.png")

# proj2 = ds.h.proj(2, field = "particle_density", weight_field = None, data_source = source)
# frb = FixedResolutionBuffer(proj2, (center[0] - (w[0]/ds[w[1]])/2.0, center[0] + (w[0]/ds[w[1]])/2.0,
#                                     center[1] - (w[0]/ds[w[1]])/2.0, center[1] + (w[0]/ds[w[1]])/2.0), (800, 800))
# pylab.imshow(np.log10(frb["particle_density"]), interpolation='nearest', vmin=-6, vmax=-1, origin='lower')
# pylab.clim(-6, -1)
# pylab.colorbar()
# pylab.savefig("data0040_Projection_z_particle_density_subset.png")
# pylab.clf()
