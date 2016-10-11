import numpy as np
import pylab as pl
pl.style.use('thesisstyle')


number_cells = [28, 40, 50, 60, 70, 80, 100]

err_inf = []
err_l2 = []

for index in range(np.size(number_cells)):
    err_inf.append(np.load('./n'
                           +str(number_cells[index])
                           +'/errfile_inf_sllb_interp.npy'))
    err_l2.append(np.load('./n'
                          +str(number_cells[index])
                          +'/errfile_l2_sllb_interp.npy'))

# err_inf.append(np.load('./n28/errfile_inf_sllb_interp.npy'))
# err_inf.append(np.load('./n40/errfile_inf_sllb_interp.npy'))
# err_inf.append(np.load('./n50/errfile_inf_sllb_interp.npy'))
# err_inf.append(np.load('./n60/errfile_inf_sllb_interp.npy'))
# err_inf.append(np.load('./n70/errfile_inf_sllb_interp.npy'))
# err_inf.append(np.load('./n80/errfile_inf_sllb_interp.npy'))
# err_inf.append(np.load('./n100/errfile_inf_sllb_interp.npy'))


# err_l2.append(np.load('./n28/errfile_l2_sllb_interp.npy'))
# err_l2.append(np.load('./n40/errfile_l2_sllb_interp.npy'))
# err_l2.append(np.load('./n50/errfile_l2_sllb_interp.npy'))
# err_l2.append(np.load('./n60/errfile_l2_sllb_interp.npy'))
# err_l2.append(np.load('./n70/errfile_l2_sllb_interp.npy'))
# err_l2.append(np.load('./n80/errfile_l2_sllb_interp.npy'))
# err_l2.append(np.load('./n100/errfile_l2_sllb_interp.npy'))

index_t = -1
number_cells = np.asarray(number_cells)
err = np.zeros_like(number_cells) * 0.0
err2 = np.zeros_like(number_cells) * 0.0
x3 = np.zeros_like(number_cells) * 0.0
x4 = np.zeros_like(number_cells) * 0.0
x1 = np.zeros_like(number_cells) * 0.0

for ind_n in range(np.size(number_cells)):
    err[ind_n] = err_inf[ind_n][index_t]
    err2[ind_n] = err_l2[ind_n][index_t]
    x3[ind_n]  = (number_cells[ind_n])**-3 * 100
    x1[ind_n]  = (number_cells[ind_n])**-1/10


fig, ax = pl.subplots(1)
pl.plot(number_cells, err, '--')
pl.plot(number_cells, err2)
pl.plot(number_cells, x3,  ':', linewidth=5)
pl.plot(number_cells, x1, ':', linewidth=5)
#pl.plot(number_cells, x5, ':', linewidth=5)
pl.legend(['$L_\infty$ error', \
           '$L_2$ error', \
#           'x3', \
#           'x4', \
           '$x^{-3}$', '$x^{-1}$'], loc='best')
pl.xlabel("Number of points $N=N_1=N_2$")
pl.xlim(number_cells[0]-2,number_cells[-1]+5)
ax.set_xlim
ax.set_yscale('log')
ax.set_xscale('log')

pl.xticks(number_cells,[str(e) for e in number_cells])
pl.show(block=True)
pl.close()
