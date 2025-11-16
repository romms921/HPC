#!/usr/bin/env python
import glafic

glafic.init(0.3, 0.7, -1.0, 0.7, '/home/rommulus/Projects/itng_lensing/Simulations/Input/Neirenberg/image_plane/NFW+SHEAR+POSFLUX', -3.5, -3.5, 3.5, 3.5, 0.01, 0.01, 1, verb = 0)
glafic.set_secondary('chi2_splane 1', verb = 0)
glafic.set_secondary('chi2_checknimg 0', verb = 0)
glafic.set_secondary('chi2_restart   -1', verb = 0)
glafic.set_secondary('chi2_usemag    0', verb = 0)
glafic.set_secondary('hvary          0', verb = 0)
glafic.set_secondary('ran_seed -122000', verb = 0)

glafic.startup_setnum(3, 0, 1)
glafic.set_lens(1, 'anfw', 0.2900, 7.197018e+11,  -1.210343e-02, 5.317669e-03, 7.949442e-02, -3.229676e+01,  4.503662e+01,  0.000000e+00)
glafic.set_lens(2, 'sers', 0.2900, 1.164145e+10, 0.0, 0.0, 0.0815089, 9.53645, 0.6702, 2.22456)
glafic.set_lens(3, 'pert', 0.2900, 1.7130, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0)
glafic.set_point(1, 1.7130, 0.0, 0.0)

glafic.setopt_lens(1, 0, 1, 1, 1, 1, 1, 1, 0)
glafic.setopt_lens(2, 0, 1, 1, 1, 0, 0, 0, 0)
glafic.setopt_lens(3, 0, 0, 0, 0, 1, 1, 0, 0)
glafic.setopt_point(1, 0, 1, 1)

glafic.model_init(verb = 0)

glafic.readobs_point('/home/rommulus/Projects/itng_lensing/Simulations/Input/Neirenberg/Eobs_fluxpoint.dat')
glafic.parprior('/home/rommulus/Projects/itng_lensing/Simulations/Input/Neirenberg/prior_NFW_gau_G.dat')
glafic.optimize()
glafic.findimg()
glafic.writecrit(1.7130)
glafic.writelens(1.7130)

glafic.quit()
