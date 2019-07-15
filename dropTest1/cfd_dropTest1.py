from itasca import p2pLinkClient
import numpy as np
from pyDemFoam import pyDemIcoFoam

solver = pyDemIcoFoam()

pfc_link = p2pLinkClient()
#pfc_link.connect("10.0.2.2") # for VirtualBox
pfc_link.connect("")        # for WSL

pfc_link.send_data(solver.nodes())
pfc_link.send_data(solver.elements())
pfc_link.send_data(solver.rho())
pfc_link.send_data(solver.mu())

while True:
    print "waiting for run time"
    deltat = pfc_link.read_data()
    if deltat == 0.0:
        print "solve finished"
        break
    print "got run time", deltat
    solver.n(pfc_link.read_data())
    solver.f(pfc_link.read_data())
    print "got runtime and data"

    solver.set_dt(deltat)
    solver.solve(deltat)

    print "sending data to pfc"
    pfc_link.send_data(solver.p()*solver.rho())
    pfc_link.send_data(solver.gradp()*solver.rho())
    pfc_link.send_data(solver.U())
    print "send finished"

pfc_link.close()
del pfc_link