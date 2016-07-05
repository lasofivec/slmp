import pylab as pl
import prettyplotlib as ppl
from prettyplotlib import brewer2mpl
import numpy as np
import os, sys
from globals_variables import *
# Plotting variables:
levels = np.linspace(PLOT_VAL_MIN, PLOT_VAL_MAX, 25)
pl.style.use('thesisstyle')


def plot_nrb_dens(X_mat, Y_mat, zn, title = '', show = True, save = False, tstep = nstep) :
    listz = []
    pl.clf()
    for npat,nrb in enumerate(geo) :
        Z = zn[npat]
        listz.append(Z)
        # #TO check the mesh :
        # pl.scatter(Y_mat[npat],X_mat[npat])
        # pl.show(block = True)
        # pl.clf()

        pl.contourf(X_mat[npat], Y_mat[npat], Z, levels=levels, cmap=pl.cm.get_cmap(comap))
    pl.title(title)
    pl.colorbar(orientation='horizontal')
    pl.grid()
    pl.xlabel('x')
    pl.ylabel('y')
    pl.axis('image')
    if (show) :
        pl.show(block=True)
    if (save) :
        if ((tstep == 0) or (tstep == -1)) :
            os.system("rm results/results_figures/*.eps")
            os.system("rm results/results_figures/exact_values/*.eps")
        nzeros = ""
        for i in range(len(str(nstep))-len(str(tstep))) :
            nzeros = nzeros + "0"
        pl.savefig("results/results_figures/Density_T_eq_" + nzeros + str(tstep) + ".eps", format='eps', dpi=1000, facecolor='w', edgecolor='none')
    pl.close()
    return listz


def comp_err_time(geo, X, Y, dx, dy, func_init, list_zi, list_errs, \
                  show = False, plot = False, save = False, ax = True, block = False, tval = nstep) :

    nzeros = ""
    for i in range(len(str(nstep))-len(str(tval))) :
        nzeros = nzeros + "0"
    tstep = nzeros + str(tval)

    if (which_f == 5) :
        vi = 0.
        vm = 1.1

    if (plot) :
        pl.clf()

    list_err_inf = []
    list_err_l2  = []
    list_min = []
    list_mass = []
    for npat,nrb in enumerate(geo) :
        # Calculating the corresponding values of knots on physical space
        Xm = dx[npat].reshape((NPTS1*NPTS2))#X - dx[npat]
        Ym = dy[npat].reshape((NPTS1*NPTS2))#Y - dy[npat]

        # Calculation the value on these points :
        zm = func_init(Xm, Ym)
        z  = zm.reshape((NPTS1, NPTS2))

        # Compute error:
        from numpy import linalg as LA
        err_inf = 1./NPTS1/NPTS2 * LA.norm(list_zi[npat]-z, ord=np.inf)
        err_l2  = 1./NPTS1/NPTS2 * LA.norm(list_zi[npat]-z)
        minval = np.min(list_zi[npat])
        sumval = 1./NPTS1/NPTS2 * np.sum(list_zi[npat])
        list_err_inf.append(err_inf)
        list_err_l2.append(err_l2)
        list_min.append(minval)
        list_mass.append(sumval)
        if (show) :
            diff = list_zi[npat]-z
            diff_res = np.amax(np.abs(diff))
            print " * Erreur (L_inf) patch ", npat, " = ", err_inf
            print " * Erreur (L_2)   patch ", npat, " = ", err_l2

        if (plot) :
            pl.contourf(X[npat], Y[npat], z, levels = levels, cmap=pl.cm.get_cmap(comap))

    if (plot) :
        pl.colorbar(orientation='horizontal')
        pl.grid()
        time = tval*dt
        pl.title('Analytical solution of the advection equation at $t =$ ' + str(time)+'\n\nwith '+func_formula)
        pl.xlabel('x')
        pl.ylabel('y')
        pl.axis('image')
        if (save) :
            pl.savefig("results/results_figures/exact_values/exact_density_t"+tstep+".eps", \
                       format='eps', dpi=1000, facecolor='w', edgecolor='none')
        if (block) :
            pl.show(block=True)

    list_errs.append(list_err_inf)
    list_errs.append(list_err_l2)
    list_errs.append(list_min)
    list_errs.append(list_mass)


def plot_errors(list_errs):

    list_err_inf = list_errs[0]
    list_err_l2  = list_errs[1]
    list_minval  = list_errs[2]
    list_mass    = list_errs[3]

    ntime = np.shape(list_err_inf)[0]
    npats = np.shape(list_err_inf)[1]

    list_emt_inf = []
    list_emt_l2  = []
    list_emt_min = []
    list_emt_mass = []
    list_tmp = []

    for tstep in range(ntime) :
        list_tmp.append(tstep*dt)
        max_err_inf = 0.
        max_err_l2  = 0.
        min_val = 0.
        mass_val = 0.
        for n in range(npats) :
            max_err_inf = max(list_err_inf[tstep][n], max_err_inf)
            max_err_l2  = max(list_err_l2[tstep][n], max_err_l2)
            min_val = min(list_minval[tstep][n], min_val)
            mass_val += list_mass[tstep][n]
        list_emt_inf.append(max_err_inf)
        list_emt_l2.append(max_err_l2)
        list_emt_min.append(min_val)
        list_emt_mass.append(mass_val)
    list_emt_mass = [e - list_emt_mass[0] for e in list_emt_mass]

    fig, ax = pl.subplots(1)
    pl.title("$L_2$ and $L_\infty$ errors over time")
    pl.xlabel("Time")
    #    ppl.legend(ax, loc='upper left')
    pl.plot(list_tmp, list_emt_inf, '-')
    pl.plot(list_tmp, list_emt_l2, '--')
    pl.legend(["$L^\infty$ error", "$L^2$ error"], loc='best')
    # *** Saving image :
    fig.savefig("results/results_figures/Error_over_time.eps", \
                format='eps', dpi=1000, facecolor='w', edgecolor='none')
    # *** Showing image :
    pl.show(block = True)
    pl.close()

    fig, ax = pl.subplots(1)
    pl.title("Time evolution of minimal value of $f(t,x,y)$")
    pl.xlabel("Time")
    #    ppl.legend(ax, loc='upper left')
    pl.plot(list_tmp, list_emt_min, '-')
    # *** Saving image :
    fig.savefig("results/results_figures/Minval_over_time.eps", \
                format='eps', dpi=1000, facecolor='w', edgecolor='none')
    # *** Showing image :
    pl.show(block = True)
    pl.close()

    fig, ax = pl.subplots(1)
    pl.title("Time evolution of $\sum_{i,j}f(t,x_i,y_j) - \sum_{i,j}f(0,x_i,y_j)$")
    pl.xlabel("Time")
    #    ppl.legend(ax, loc='upper left')
    pl.plot(list_tmp, list_emt_mass, '-')
    # *** Saving image :
    fig.savefig("results/results_figures/Mass_over_time.eps", \
                format='eps', dpi=1000, facecolor='w', edgecolor='none')
    # *** Showing image :
    pl.show(block = True)
    pl.close()
    ############################
    # *** Save error curve :

    # #if comparing number of points :
    # str_num = str(NPTS1)+'_'+str(NPTS2) 
    # path = "results/errors/space_disc"
    # #if comparing derivative approximation
    # str_num = str(which_deriv)
    # path = "results/errors"
    # #if comparing interpolation method
    if (which_interp == 0) :
        str_num = "_sllb_interp"
    else :
        str_num = "_scipy_"+str(which_interp)
    path = "results/errors/interp_app"

    np.save(path + '/timesteps'   + str_num, list_tmp)
    np.save(path + '/errfile_inf' + str_num, list_emt_inf)
    np.save(path + '/errfile_l2'  + str_num, list_emt_l2)
    write_globals(path, str_num)

    return np.max(list_err_inf)
