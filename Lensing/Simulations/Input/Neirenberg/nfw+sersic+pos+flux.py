#!/usr/bin/env python
import glafic
import os

os.chdir('/Users/ainsleylewis/Documents/Astronomy/HPC/Lensing/Simulations/Input/Neirenberg')

glafic.init(0.3, 0.7, -1.0, 0.7, 'image_plane/NFW+SERS+POSFLUX2', -3.5, -3.5, 3.5, 3.5, 0.0025, 0.0025, 1, verb = 0)
glafic.set_secondary('chi2_splane 0', verb = 0)
glafic.set_secondary('chi2_checknimg 0', verb = 0)
glafic.set_secondary('chi2_restart   -1', verb = 0)
glafic.set_secondary('chi2_usemag    0', verb = 0)
glafic.set_secondary('hvary          0', verb = 0)
glafic.set_secondary('ran_seed -122000', verb = 0)

glafic.startup_setnum(2, 0, 1)
glafic.set_lens(1, 'anfw', 0.2900, 8.336294e+10,  1.317348e-02, -5.360185e-02,  3.384265e-01, -4.311449e+01,  7.915958e+01,  0.000000e+00)
glafic.set_lens(2, 'sers', 0.2900, 6.590989e+10,  2.576273e-02,  5.363812e-03,  8.150890e-02,  9.536450e+00,  6.702000e-01,  2.224560e+00)
glafic.set_point(1, 1.7130, 3.126457e-02, -4.735896e-03)

glafic.setopt_lens(1, 0, 0, 0, 0, 0, 0, 0, 0)
glafic.setopt_lens(2, 0, 0, 0, 0, 0, 0, 0, 0)
glafic.setopt_point(1, 0, 0, 0)

glafic.model_init(verb = 0)

# glafic.readobs_point('Eobs_fluxpoint.dat')
# glafic.parprior('prior_NFW_gau_G.dat')
# glafic.optimize()
# glafic.findimg()
# glafic.writecrit(1.7130)
glafic.writelens(1.7130)

glafic.quit()
