##################################################################
# Flageul et al.
# http://dx.doi.org/10.1016/j.ijheatfluidflow.2015.07.009
# Reads xls(x) files and generates figures
##################################################################
from numpy import *
from pylab import *
import xlrd
matplotlib.rc('figure', figsize=(1.3*5.83,1.3*4.13))
matplotlib.rc('text', usetex = True)
size=16
size_legend=14
size_label=20
linewidth=1.5
markersize=10
matplotlib.rc('lines', linewidth=linewidth,markersize=markersize)
matplotlib.rc('font', size=size)
matplotlib.rc('axes', labelsize=size_label, titlesize=size)
matplotlib.rc('legend', fontsize=size_legend)


inc =  xlrd.open_workbook('autocorr.xls')

#
def purge(L):
	return [x for x in L if x <> '']

# open sheet in excel
sh1 = inc.sheet_by_name('solide')
# y-coord
y=purge(sh1.col_values(29)[2:130])
# look for data
t0=purge(sh1.col_values(30)[2:130])
t1=purge(sh1.col_values(31)[2:130])
t2=purge(sh1.col_values(32)[2:130])
t3=purge(sh1.col_values(33)[2:130])

# plot velocity
plot(y,t0,'--',color='r',label=r'$y^+=0$')

plot(y,t1,'-+',color='k',label=r'$y^+=-5.3$')

plot(y,t2,'-x',color='b',label=r'$y^+=-15$')

plot(y,t3,'-o',color='g',markerfacecolor='none',markeredgecolor='g',label=r'$y^+=-77$')

#Graph settings
xscale('linear')
yscale('linear')
axis([-0.1,650,-0.4,1.1])
xlabel(r"$z^+$")
ylabel(r"Spanwise autocorrelation ($T'$)")
legend(bbox_to_anchor=(0.7,0.7),numpoints=1)

savefig("all_autocorr_z_sol.png",bbox_inches='tight')
savefig("all_autocorr_z_sol.pdf",bbox_inches='tight')
