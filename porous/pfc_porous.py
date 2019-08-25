import itasca as it
from pfc_cfd_coupler.pfc_coupler import pfc_coupler

coupler = pfc_coupler()
it.command("call ../porous/particles.p3dat")

coupler.max_dt = 0.005
coupler.bandwidth = 0
coupler.solve(0.1)
coupler.plotFluidVel()
coupler.plotPorosity()
coupler.stopSolve()
coupler.close()