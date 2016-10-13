import numpy as np
import pylab as pl
pl.style.use('thesisstyle')


err_inf = []
err_inf.append(np.load('./dt001/n60/errfile_inf_sllb_interp.npy'))
#err_inf.append(np.load('./dt00125/n60/errfile_inf_sllb_interp.npy'))
err_inf.append(np.load('./dt002/n60/errfile_inf_sllb_interp.npy'))
#err_inf.append(np.load('./dt0025/n60/errfile_inf_sllb_interp.npy'))
err_inf.append(np.load('./dt004/n60/errfile_inf_sllb_interp.npy'))
#err_inf.append(np.load('./dt005/n60/errfile_inf_sllb_interp.npy'))
err_inf.append(np.load('./dt008/n60/errfile_inf_sllb_interp.npy'))
err_inf.append(np.load('./dt01/n60/errfile_inf_sllb_interp.npy'))
# err_inf.append(np.load('./dt02/n60/errfile_inf_sllb_interp.npy'))
# err_inf.append(np.load('./dt04/n60/errfile_inf_sllb_interp.npy'))
# err_inf.append(np.load('./dt08/n60/errfile_inf_sllb_interp.npy'))
# err_inf.append(np.load('./dt1/n60/errfile_inf_sllb_interp.npy'))

err_l2 = []
err_l2.append(np.load('./dt001/n60/errfile_l2_sllb_interp.npy'))
#err_l2.append(np.load('./dt00125/n60/errfile_l2_sllb_interp.npy'))
err_l2.append(np.load('./dt002/n60/errfile_l2_sllb_interp.npy'))
#err_l2.append(np.load('./dt0025/n60/errfile_l2_sllb_interp.npy'))
err_l2.append(np.load('./dt004/n60/errfile_l2_sllb_interp.npy'))
#err_l2.append(np.load('./dt005/n60/errfile_l2_sllb_interp.npy'))
err_l2.append(np.load('./dt008/n60/errfile_l2_sllb_interp.npy'))
err_l2.append(np.load('./dt01/n60/errfile_l2_sllb_interp.npy'))
# err_l2.append(np.load('./dt02/n60/errfile_l2_sllb_interp.npy'))
# err_l2.append(np.load('./dt04/n60/errfile_l2_sllb_interp.npy'))
# err_l2.append(np.load('./dt08/n60/errfile_l2_sllb_interp.npy'))
# err_l2.append(np.load('./dt1/n60/errfile_l2_sllb_interp.npy'))


index_t = -1
time_step  = [0.01, 0.02, 0.04, 0.08, 0.1]#[0.01, 0.0125, 0.02, 0.025, 0.04, 0.05, 0.08, 0.1]#, 0.2, 0.4, 0.8, 1.0]
time_step = np.asarray(time_step)
err = np.zeros_like(time_step) * 0.0
err2 = np.zeros_like(time_step) * 0.0
x2 = np.zeros_like(time_step) * 0.0

for ind_n in range(np.size(time_step)):
    err[ind_n] = err_inf[ind_n][index_t]
    err2[ind_n] = err_l2[ind_n][index_t]
    x2[ind_n]  = (time_step[ind_n])**-2/10000000.


fig, ax = pl.subplots(1)
pl.plot(time_step, err, '--')
pl.plot(time_step, err2)
pl.plot(time_step, x2, '--')
pl.legend(['$L_\infty$ error', \
           '$L_2$ error', \
           '$x_2$'], loc='best')
pl.xlabel("Time step $\Delta t$")
#pl.xlim(0,1)
ax.set_xlim
ax.set_yscale('log')
ax.set_xscale('log')

pl.xticks(time_step,[str(e) for e in time_step])
pl.show(block=True)
pl.close()
