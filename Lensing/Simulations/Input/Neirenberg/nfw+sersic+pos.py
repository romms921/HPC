#!/usr/bin/env python
import glafic

glafic.init(0.3, 0.7, -1.0, 0.7, '/home/rommulus/Projects/itng_lensing/Simulations/Input/Neirenberg/source_plane/NFW+SERSIC+POS', -3.5, -3.5, 3.5, 3.5, 0.01, 0.01, 1, verb = 0)
glafic.set_secondary('chi2_splane 0', verb = 0)
glafic.set_secondary('chi2_checknimg 0', verb = 0)
glafic.set_secondary('chi2_restart   -1', verb = 0)
glafic.set_secondary('chi2_usemag    0', verb = 0)
glafic.set_secondary('hvary          0', verb = 0)
glafic.set_secondary('ran_seed -122000', verb = 0)

glafic.startup_setnum(2, 0, 1)
glafic.set_lens(1, 'anfw', 0.2900, 7.143805e+11, 4.139083e-02, -8.519295e-03, 9.405853e-03, -6.266181e+01, 4.481794e+01, 0.000000e+00)
glafic.set_lens(2, 'sers', 0.2900, 1.175461e+10, 7.484725e-04, 5.519739e-03, 8.150890e-02, 9.536450e+00, 6.702000e-01, 2.224560e+00)
glafic.set_point(1, 1.7130, 3.683559e-02, -6.453691e-03)

glafic.setopt_lens(1, 0, 1, 1, 1, 1, 1, 1, 0)
glafic.setopt_lens(2, 0, 1, 1, 1, 0, 0, 0, 0)
glafic.setopt_point(1, 0, 1, 1)

glafic.model_init(verb = 0)

glafic.readobs_point('/home/rommulus/Projects/itng_lensing/Simulations/Input/Neirenberg/Eobs_point.dat')
glafic.parprior('/home/rommulus/Projects/itng_lensing/Simulations/Input/Neirenberg/prior_NFW_gau_G.dat')
glafic.optimize()
glafic.findimg()
glafic.writecrit(1.7130)
glafic.writelens(1.7130)

glafic.quit()
