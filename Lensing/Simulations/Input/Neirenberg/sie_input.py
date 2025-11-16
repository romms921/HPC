#!/usr/bin/env python
import glafic

glafic.init(0.3, 0.7, -1.0, 0.7, '/home/rommulus/Projects/itng_lensing/Simulations/Input/Neirenberg/image_plane/NFW+SERSIC+POS+FLUX', -3.5, -3.5, 3.5, 3.5, 0.0025, 0.0025, 1, verb = 0)
glafic.set_secondary('chi2_splane 1', verb = 0)
glafic.set_secondary('chi2_checknimg 0', verb = 0)
glafic.set_secondary('chi2_restart   -1', verb = 0)
glafic.set_secondary('chi2_usemag    0', verb = 0)
glafic.set_secondary('hvary          0', verb = 0)
glafic.set_secondary('ran_seed -122000', verb = 0)

glafic.startup_setnum(2, 0, 1)
glafic.set_lens(1, 'anfw', 0.2900, 7.627179e+11, 3.409167e-02, -1.322709e-02, 9.034902e-03, -5.274117e+01, 4.464198e+01, 0.000000e+00)
glafic.set_lens(2, 'sers', 0.2900, 9.055939e+09, -1.639984e-03,  3.946928e-03,  8.150890e-02,  9.536450e+00,  6.702000e-01,  2.224560e+00)
# glafic.set_lens(3, 'pert', 0.2900, 1.7130, 0.0, 0.0, 4.298076e-02, 1.149450e+01, 0.000000e+00, 0.000000e+00)
glafic.set_point(1, 1.7130, 3.124038e-02, -1.107056e-02)

glafic.setopt_lens(1, 0, 0, 0, 0, 0, 0, 0, 0)
glafic.setopt_lens(2, 0, 0, 0, 0, 0, 0, 0, 0)
# glafic.setopt_lens(3, 0, 0, 0, 0, 0, 0, 0, 0)
glafic.setopt_point(1, 0, 0, 0)

glafic.model_init(verb = 0)

# glafic.readobs_point('/home/rommulus/Projects/itng_lensing/Simulations/Input/Neirenberg/Eobs_fluxpoint.dat')
# glafic.parprior('/home/rommulus/Projects/itng_lensing/Simulations/Input/Neirenberg/prior_NFW_gau_G.dat')
# glafic.optimize()
# glafic.findimg()
# glafic.writecrit(1.7130)
glafic.writelens(1.7130)

glafic.quit()
