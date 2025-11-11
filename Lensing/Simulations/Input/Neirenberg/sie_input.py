#!/usr/bin/env python
import glafic

path = '/Volumes/T7 Shield/Simulations/Output/SIE_SHEAR_0.001_360.0_1.0'
constraint_file = '/Volumes/T7 Shield/Simulations/Input/pos_point.dat'
prior_file = None
critical_curve = False

glafic.init(0.3, 0.7, -1.0, 0.7, path, -3.5, -3.5, 3.5, 3.5, 0.01, 0.01, 1, verb = 0)
glafic.set_secondary('chi2_splane 1', verb = 0)
glafic.set_secondary('chi2_checknimg 0', verb = 0)
glafic.set_secondary('chi2_restart   -1', verb = 0)
glafic.set_secondary('chi2_usemag    0', verb = 0)
glafic.set_secondary('hvary          0', verb = 0)
glafic.set_secondary('ran_seed -122000', verb = 0)
glafic.startup_setnum(1, 0, 1)
glafic.set_lens(1, 'sie', 0.2300, 2.368933e+02, 0.0, 0.0, 0.1, 20, 0.0, 0.0)
glafic.set_point(1, 0.777, 0.0, 0.0)
glafic.setopt_lens(1, 0, 0, 1, 1, 0, 0, 0, 0)
glafic.setopt_point(1, 0, 1, 1)
glafic.model_init(verb = 0)
glafic.readobs_point(constraint_file)
if prior_file:
    glafic.parprior(prior_file)
glafic.optimize()
glafic.findimg()
if critical_curve:
    glafic.writecrit(0.777)
# glafic.writelens(1.0)
glafic.quit()
