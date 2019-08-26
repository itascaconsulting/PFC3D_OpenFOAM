import itasca as it
from pfc_cfd_coupler.pfc_coupler import pfc_coupler

coupler = pfc_coupler()
it.command("call ../porous/particles.p3dat")

coupler.max_dt = 0.005
coupler.bandwidth = 0.01
coupler.smallest_size = 0.005
coupler.pow1 = 7.0
coupler.pow2 = 10.0
coupler.pow3 = 10.0
coupler.initialize()
coupler.solve(coupler.max_dt)
#time = 0.0
#while time < 0.01:
#    dt = coupler.max_dt
#    coupler.solve(dt)
#    it.fish.set("cfd_vel",coupler.elements_vel[0][1])
#    time += dt
coupler.plotFluidVel()
coupler.plotPorosity()
coupler.stopSolve()
coupler.close()