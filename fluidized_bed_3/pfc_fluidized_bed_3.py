import itasca as it
from pfc_cfd_coupler.pfc_coupler import pfc_coupler

coupler = pfc_coupler()
it.command("call ../fluidized_bed_3/make_ini.p3dat")

coupler.max_dt = 0.001
coupler.smallest_size = 0.00095
coupler.bandwidth = 2*coupler.smallest_size
coupler.initialize()
coupler.solve(4.0)
coupler.plotFluidVel()
coupler.plotPorosity()
coupler.stopSolve()
coupler.close()

it.command("history write 1,2 file 'fluidized_bed_3.txt' truncate")