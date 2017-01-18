##################################################################
# Flageul et al.
# http://dx.doi.org/10.1016/j.ijheatfluidflow.2015.07.009
# Reads xls(x) files and generates figures
##################################################################
from numpy import *
from pylab import *
import xlrd
matplotlib.rc('figure', figsize=(5.83,4.13))
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


ref =  xlrd.open_workbook('ref_scal_neumann.xls')
inc =  xlrd.open_workbook('neumann.xls')

#
def purge(L):
	return [x for x in L if x <> '']

# open sheet in excel
sh1 = inc.sheet_by_name('fluctuations')
sh2 = ref.sheet_by_name('diric')
# y-coord
y_inc_bot=purge(sh1.col_values(9)[5:101])
y_inc_top=purge(sh1.col_values(17)[5:101])
y_ref=purge(sh2.col_values(1)[1:49])
# look for data
ut_bot=purge(sh1.col_values(14)[5:101])
vt_bot=purge(sh1.col_values(15)[5:101])
tt_bot=purge(sh1.col_values(16)[5:101])
ut_ref=purge(sh2.col_values(4)[1:49])
vt_ref=purge(sh2.col_values(5)[1:49])
tt_ref=purge(array(sh2.col_values(3)[1:49])**2)
# plot velocity
plot(y_inc_bot,ut_bot,'-',color='r',label=r"$<u'T'>$")
plot(y_inc_bot,vt_bot,'--',color='g',label=r"$<v'T'>$")
plot(y_inc_bot,tt_bot,'-.',color='b',label=r"$<T'^2>$")
plot(y_ref,ut_ref,'+',color='r')
plot(y_ref,vt_ref,'+',color='g')
plot(y_ref,tt_ref,'+',color='b')

#Graph settings
xscale('log')
yscale('log')
axis([0.3,150,10**-4,7])
xlabel(r"$y^+$")
ylabel(r"$<u'T'>$, $<v'T'>$, $<T'^2>$")
legend(bbox_to_anchor=(0.9,0.7),numpoints=1)

savefig("isoq_ut_vt_tt_log.png",bbox_inches='tight')
savefig("isoq_ut_vt_tt_log.pdf",bbox_inches='tight')

xscale('linear')
yscale('linear')
savefig("isoq_ut_vt_tt_lin.png",bbox_inches='tight')
savefig("isoq_ut_vt_tt_lin.pdf",bbox_inches='tight')
