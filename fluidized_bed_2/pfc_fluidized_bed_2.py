import itasca as it
from pfc_cfd_coupler.pfc_coupler import pfc_coupler

coupler = pfc_coupler()
it.command("""
set timestep max 1e-4
set deterministic on
set random 1
ball generate rad 0.00075 number 2160 box 0 0.04 0 0.06 0 0.0075
wall generate box 0 0.04 0 0.2 0 0.0075
ball ini dens 2000
cmat default model linear property kn 5000 ks 5000 fric 0.15 dp_nratio 0.1
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
define pressure_drop
    pressure_drop = cfd_pressure_drop 
end
def fluid_time
  global fluid_time = mech.age
end
history add id 1 fish @fluid_time
history add id 2 fish @ball_height
plot clear
plot add ball
plot add axes
plot add udvector
plot add wall transparency 70
cy 20000
ball ini damp 0.7
cy 20000
set mechanical age 0.0
history add id 3 fish @pressure_drop
plot add hist 3 vs 1
""")

coupler.dt = 0.0001
#coupler.bandwidth = 0
coupler.solve(100)
coupler.plotFluidVel()
coupler.close()

it.command("history write 1,2,3 file 'fluidized_bed_2.txt' truncate")