import itasca as it
from pfc_cfd_coupler.pfc_coupler import pfc_coupler

coupler = pfc_coupler()
it.command("""
ball generate rad 0.00075 number 100 box 0 0.04 0 0.005 0 0.0075
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
history add id 1 fish @fluid_time
history add id 2 fish @ball_height
plot clear
plot add hist 2 vs 1
plot add ball shape arrow
plot add axes
plot add domain
plot add udvector
""")

coupler.dt = 0.0001
#coupler.bandwidth = 0
coupler.solve(100)
coupler.plotFluidVel()
coupler.close()

it.command("history write 1,2 file 'fluidized_bed_2.txt' truncate")