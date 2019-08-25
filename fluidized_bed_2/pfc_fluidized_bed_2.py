import itasca as it
from pfc_cfd_coupler.pfc_coupler import pfc_coupler

coupler = pfc_coupler()
it.command("""
restore ../fluidized_bed_2/ini.sav
set timestep max 1e-4
set mechanical age 0.0
history add id 3 fish @pressure_drop
plot add hist 3 vs 1
""")

coupler.dt = 0.0001
coupler.bandwidth = 0.015
coupler.pressureMeasureCell1 = 12
coupler.pressureMeasureCell2 = 1887
coupler.solve(200)
coupler.plotFluidVel()
coupler.stopSolve()
coupler.close()

it.command("history write 1,2,3 file 'fluidized_bed_2.txt' truncate")