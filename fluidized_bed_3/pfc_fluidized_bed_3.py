import itasca as it
from pfc_cfd_coupler.pfc_coupler import pfc_coupler

coupler = pfc_coupler()
it.command("call ../fluidized_bed_3/make_ini.p3dat")

coupler.max_dt = 0.00001
coupler.bandwidth = 0
coupler.smallest_size = 0.00095
coupler.initialize()
coupler.solve(2.0)
coupler.plotFluidVel()
coupler.plotPorosity()
coupler.stopSolve()
coupler.close()

#it.command("history write 1,2,3 file 'fluidized_bed_2.txt' truncate")