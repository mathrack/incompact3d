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


refd =  xlrd.open_workbook('dirichlet.xls')
refn =  xlrd.open_workbook('neumann.xls')
inc =  xlrd.open_workbook('robin.xls')
cht =  xlrd.open_workbook('conju_ref.xls')

#
def purge(L):
	return [x for x in L if x <> '']

# open sheet in excel
sh1 = inc.sheet_by_name('fluctuations')
sh2d= refd.sheet_by_name('fluctuations')
sh2 = refn.sheet_by_name('fluctuations')
sh3 = cht.sheet_by_name('fluctuations')
# y-coord
y_inc_bot=purge(sh1.col_values(9)[5:101])
y_inc_top=purge(sh1.col_values(17)[5:101])
y_ref=purge(sh2.col_values(9)[5:101])
y_refd=purge(sh2d.col_values(9)[5:101])
y_cht=purge(sh3.col_values(9)[6:101])
# look for data
ut_bot=purge(sh1.col_values(14)[5:101])
vt_bot=purge(sh1.col_values(15)[5:101])
tt_bot=purge(sh1.col_values(16)[5:101])
ut_ref=purge(sh2.col_values(14)[5:101])
vt_ref=purge(sh2.col_values(15)[5:101])
tt_ref=purge(sh2.col_values(16)[5:101])
ut_refd=purge(sh2d.col_values(14)[5:101])
vt_refd=purge(sh2d.col_values(15)[5:101])
tt_refd=purge(sh2d.col_values(16)[5:101])
ut_cht=purge(sh3.col_values(14)[6:101])
vt_cht=purge(sh3.col_values(15)[6:101])
tt_cht=purge(sh3.col_values(16)[6:101])
# plot velocity
#plot(y_inc_bot,ut_bot,'-',color='r',label=r'$uT$')
plot(y_inc_bot,vt_bot,'o',color='r',label=r'$Robin$',markerfacecolor='none',markeredgecolor='r')
#plot(y_inc_bot,tt_bot,'-.',color='r',label=r'$TT$')
#plot(y_ref,ut_ref,'-x',color='g')
plot(y_ref,vt_ref,'x',color='g',label=r'$isoQ$')
#plot(y_ref,tt_ref,'-.x',color='g')
#plot(y_refd,ut_refd,'-+',color='b')
plot(y_refd,vt_refd,'+',color='b',label=r'$isoT$')
#plot(y_refd,tt_refd,'-.+',color='b')
plot(y_cht,vt_cht,'-',color='k',label=r'$Conjug$')

#Graph settings
xscale('log')
yscale('log')
axis([0.3,150,0.2*(10**-4),1])
xlabel(r"$y^+$")
ylabel(r"$<v'T'>$")
legend(bbox_to_anchor=(0.9,0.7),numpoints=1)

savefig("all_vtcht_log.png",bbox_inches='tight')
savefig("all_vtcht_log.pdf",bbox_inches='tight')

xscale('linear')
yscale('linear')
savefig("all_vtcht_lin.png",bbox_inches='tight')
savefig("all_vtcht_lin.pdf",bbox_inches='tight')

xscale('log')
yscale('linear')
savefig("all_vtcht_linlog.png",bbox_inches='tight')
savefig("all_vtcht_linlog.pdf",bbox_inches='tight')
