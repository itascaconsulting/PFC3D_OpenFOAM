import itasca as it
from pfc_cfd_coupler.pfc_coupler import pfc_coupler

coupler = pfc_coupler()
it.command("""
cfd porosity poly
cfd buoy on
ball generate rad 0.005 number 100 box 0 1 0 0.25 0 1
ball ini dens 2500
ball prop kn 1e2 ks 1e2 fric 0.25
set gravity 0 -9.81 0
def ball_height
  local max = 0
  loop foreach local b ball.list
    if ball.pos.y(b) > max then
      max = ball.pos.y(b)
    endif
  endloop
  ball_height = max
end
def fluid_time
  global fluid_time = mech.age
end
history add id 1 fish @ball_height
history add id 2 fish @fluid_time
plot clear
plot add hist 1 vs 2
plot add cfdelement shape arrow colorby vectorattribute "velocity"
""")

coupler.dt = 0.005
coupler.solve(200)

#print "ball y velocity", it.ball.find(1).vel_y()