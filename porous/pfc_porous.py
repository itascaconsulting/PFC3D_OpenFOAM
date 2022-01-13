import itasca as it
it.command("python-reset-state false")
from pfc_cfd_coupler.pfc_coupler import pfc_coupler

coupler = pfc_coupler()
it.command("call particles.p3dat")

coupler.max_dt = 0.005
coupler.bandwidth = 0.02
coupler.smallest_size = 0.01
1/0
coupler.initialize()
time = 0.0
while time < 0.1:
    dt = coupler.max_dt
    coupler.solve(dt)
    it.fish.set("cfd_vel",coupler.elements_vel[0][1])
    time += dt
coupler.plotFluidVel()
coupler.plotPorosity()
coupler.stopSolve()
coupler.close()