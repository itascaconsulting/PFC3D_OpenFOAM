import itasca as it
from pfc_cfd_coupler.pfc_coupler import pfc_coupler

coupler = pfc_coupler()
it.command("""
cfd porosity poly
cfd buoy on
ball create rad 0.005 x 0.5 y 0.5 z 0.5
ball ini dens 2500
ball prop kn 1e2 ks 1e2 fric 0.25
set gravity 0 0 -9.81
def fluid_time
  global fluid_time = mech.age
end
ball history id 1 zvelocity id 1
history add id 2 fish @fluid_time
plot clear
plot add hist 1 vs 2
plot add cfdelement shape arrow colorby vectorattribute "velocity"
""")

coupler.solve()

print "ball z velocity", it.ball.find(1).vel_z()