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
sh2d = refd.sheet_by_name('fluctuations')
sh2n = refn.sheet_by_name('fluctuations')
sh3 = cht.sheet_by_name('fluctuations')
# y-coord
y_inc_bot=purge(sh1.col_values(9)[6:101])
y_inc_top=purge(sh1.col_values(17)[6:101])
y_refd=purge(sh2d.col_values(9)[6:101])
y_refn=purge(sh2n.col_values(9)[6:101])
y_cht=purge(sh3.col_values(9)[6:101])
# look for data
vv_bot=purge(sh1.col_values(11)[6:101])
vt_bot=purge(sh1.col_values(15)[6:101])
tt_bot=purge(sh1.col_values(16)[6:101])
vv_refd=purge(sh2d.col_values(11)[6:101])
vt_refd=purge(sh2d.col_values(15)[6:101])
tt_refd=purge(sh2d.col_values(16)[6:101])
vv_refn=purge(sh2n.col_values(11)[6:101])
vt_refn=purge(sh2n.col_values(15)[6:101])
tt_refn=purge(sh2n.col_values(16)[6:101])
vv_cht=purge(sh3.col_values(11)[6:101])
vt_cht=purge(sh3.col_values(15)[6:101])
tt_cht=purge(sh3.col_values(16)[6:101])
#
yy_kasa=purge(sh2d.col_values(26)[6:77])
vt_kasa=purge(sh2d.col_values(28)[6:77])
# plot velocity
plot(y_inc_bot,vt_bot/(sqrt(vv_bot)*sqrt(tt_bot)),'o',color='r',label=r'$Robin$',markerfacecolor='none',markeredgecolor='r')
plot(y_refn,vt_refn/(sqrt(vv_refn)*sqrt(tt_refn)),'x',color='g',label=r'$isoQ$')
plot(y_refd,vt_refd/(sqrt(vv_refd)*sqrt(tt_refd)),'+',color='b',label=r'$isoT$')
plot(y_cht,vt_cht/(sqrt(vv_cht)*sqrt(tt_cht)),'-',color='k',label=r'$Conjug$')
#plot(yy_kasa,vt_kasa,'+',color='b')

#Graph settings
xscale('log')
yscale('log')
axis([0.3,150,0.1,0.55])
xlabel(r"$y^+$")
ylabel(r"$<v'T'> / v_{RMS} T_{RMS}$")
legend(bbox_to_anchor=(0.8,0.5),numpoints=1)

savefig("corr_vtcht_log.png",bbox_inches='tight')
savefig("corr_vtcht_log.pdf",bbox_inches='tight')

xscale('linear')
yscale('linear')
axis([0.3,20,0.1,0.5])
ax = gca()
draw()
savefig("corr_vtcht.png",bbox_inches='tight')
savefig("corr_vtcht.pdf",bbox_inches='tight')
