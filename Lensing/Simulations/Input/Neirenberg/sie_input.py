#!/usr/bin/env python
import glafic

glafic.init(0.3, 0.7, -1.0, 0.7, '/home/rommulus/Projects/itng_lensing/Simulations/Input/Neirenberg/source_plane/NFW+SHEAR+POSFLUX', -3.5, -3.5, 3.5, 3.5, 0.0025, 0.0025, 1, verb = 0)
glafic.set_secondary('chi2_splane 1', verb = 0)
glafic.set_secondary('chi2_checknimg 0', verb = 0)
glafic.set_secondary('chi2_restart   -1', verb = 0)
glafic.set_secondary('chi2_usemag    0', verb = 0)
glafic.set_secondary('hvary          0', verb = 0)
glafic.set_secondary('ran_seed -122000', verb = 0)

glafic.startup_setnum(3, 0, 1)
glafic.set_lens(1, 'anfw', 0.2900, 3.693054e+11, 2.350164e-02, -4.047035e-02, 2.813979e-01, -4.662619e+01, 4.655859e+01, 0.000000e+00)
glafic.set_lens(2, 'sers', 0.2900, 3.849356e+10, 6.692701e-03, 6.302904e-03, 8.150890e-02, 9.536450e+00, 6.702000e-01, 2.224560e+00)
glafic.set_lens(3, 'pert', 0.2900, 1.7130, 0.0, 0.0, 3.076876e-02, 2.170097e+01, 0.000000e+00, 0.000000e+00)
glafic.set_point(1, 1.7130, 3.015122e-02, -1.076508e-02)

glafic.setopt_lens(1, 0, 0, 0, 0, 0, 0, 0, 0)
glafic.setopt_lens(2, 0, 0, 0, 0, 0, 0, 0, 0)
glafic.setopt_lens(3, 0, 0, 0, 0, 0, 0, 0, 0)
glafic.setopt_point(1, 0, 0, 0)

glafic.model_init(verb = 0)

# glafic.readobs_point('/home/rommulus/Projects/itng_lensing/Simulations/Input/Neirenberg/Eobs_fluxpoint.dat')
# glafic.parprior('/home/rommulus/Projects/itng_lensing/Simulations/Input/Neirenberg/prior_NFW_gau_G.dat')
# glafic.optimize()
# glafic.findimg()
# glafic.writecrit(1.7130)
glafic.writelens(1.7130)

glafic.quit()
