import glafic

path = '/Users/ainsleylewis/Documents/Astronomy/HPC/Lensing/Simulations/Input/System_4/test'
constraint_file = '/Users/ainsleylewis/Documents/Astronomy/HPC/Lensing/Simulations/Input/System_4/pos_point.dat'
prior_file = None
critical_curve = False

glafic.init(0.3089901684739047, 0.6910098315260953, -1.0, 0.6736, path, 0.0, 0.0, 2.12, 2.12, 0.01, 0.01, 1, verb=0)
glafic.set_secondary('chi2_splane 1', verb=0)
glafic.set_secondary('chi2_checknimg 0', verb=0)
glafic.set_secondary('chi2_restart   -1', verb=0)
glafic.set_secondary('chi2_usemag    1', verb=0)
glafic.set_secondary('hvary          0', verb=0)
glafic.set_secondary('ran_seed -122000', verb=0)
glafic.startup_setnum(1, 0, 1)
glafic.set_lens(1, 'sie', 0.297717684517447, 160, 1.06, 1.06, 0.1, 45, 0.0, 0.0)
# glafic.set_lens(2, 'pert', 0.297717684517447, 1.0, 1.06, 1.06, 0.001, 360.0, 0.0, 0.0)
glafic.set_point(1, 1.0, 1.06, 1.06)
glafic.setopt_lens(1, 0, 1, 1, 1, 1, 1, 0, 0)
# glafic.setopt_lens(2, 0, 0, 0, 0, 1, 1, 0, 1)
glafic.setopt_point(1, 0, 1, 1)
glafic.model_init(verb=0)
glafic.readobs_point(constraint_file)
if prior_file:
    glafic.parprior(prior_file)
glafic.optimize()
glafic.findimg()
if critical_curve:
    glafic.writecrit(1.0)
# glafic.writelens(1.0)
glafic.quit()