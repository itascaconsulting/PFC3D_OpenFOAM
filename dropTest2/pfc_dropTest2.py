import itasca as it
from pfc_cfd_coupler.pfc_coupler import pfc_coupler

coupler = pfc_coupler()
it.command("""
ball create rad 0.005 x 0.5 y 0.5 z 0.5
ball ini dens 2500
ball prop kn 1e2 ks 1e2 fric 0.25
set gravity 0 0 -9.81
def fluid_time
  global fluid_time = mech.age
end
def appforce
    global appforce = ball.force.app.z(ball.find(1))
end
history add id 1 fish @fluid_time
ball history id 2 zvelocity id 1
ball history id 3 zunbalforce id 1
history add id 4 fish @appforce
plot clear
plot add hist 2 vs 1
plot add ball shape arrow
plot add axes
plot add domain
plot add udvector
""")

coupler.dt = 0.005
coupler.bandwidth = 0.3
coupler.solve(100)
coupler.plotFluidUnitVel()
coupler.close()

print "ball z velocity", it.ball.find(1).vel_z()
it.command("history write 1,2,3,4 file 'droptest2.txt' truncate")