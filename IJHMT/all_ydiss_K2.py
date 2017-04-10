##################################################################
# Flageul et al.
# http://doi.org/10.1016/j.ijheatmasstransfer.2017.04.005
# Reads xls file and generates figures
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


refd =  xlrd.open_workbook('./diric.xls')
refn =  xlrd.open_workbook('./neuma.xls')

#
def purge(L):
	return array([x for x in L if x <> ''])

# open sheet in excel
sh2d= refd.sheet_by_name('budget_tt')
sh2 = refn.sheet_by_name('budget_tt')
# y-coord
#y_inc_bot=purge(sh1.col_values(9)[3:101])
y_ref=purge(sh2.col_values(0)[3:101])
y_refd=purge(sh2d.col_values(0)[3:101])
# look for data
#ut_bot=purge(sh1.col_values(14)[3:101])
#vt_bot=purge(sh1.col_values(15)[3:101])
#tt_bot=purge(sh1.col_values(16)[3:101])
tt_ref=-purge(sh2.col_values(1)[3:101])
tt_refd=-purge(sh2d.col_values(1)[3:101])
ydiss_ref=-purge(sh2.col_values(6)[3:101])
ydiss_refd=-purge(sh2d.col_values(6)[3:101])
# plot velocity
#plot(y_inc_bot,ut_bot,'-+',color='r',label=r'$uT$')
#plot(y_inc_bot,vt_bot,'--',color='r',label=r'$vT$')
#plot(y_inc_bot,tt_bot,'o',color='r',label=r'$Robin$',markerfacecolor='none',markeredgecolor='r')
#plot(y_ref,ut_ref,'-x',color='g')
#plot(y_ref,vt_ref,'--x',color='g')
plot(y_ref,ydiss_ref/tt_ref,'s',color='k',label=r'$isoQ$',markerfacecolor='none',markeredgecolor='k')
#plot(y_refd,ut_refd,'-+',color='b')
#plot(y_refd,vt_refd,'--+',color='b')

#Graph settings
#xscale('symlog',linthreshx=1.0)
xscale('linear')
yscale('linear')
axis([0.4,10,0.7,1.02])
xlabel(r"$y^+$")
ylabel(r"$\overline{\partial_n T'_f \partial_n T'_f} / \overline{\nabla T'_f . \nabla T'_f}$")

cht =  xlrd.open_workbook('./g05a1.xls')
sh3 = cht.sheet_by_name('budget_tt')
y_cht=purge(sh3.col_values(0)[3:101])
y_sol=purge(sh3.col_values(8)[3:135])
tt_cht=purge(sh3.col_values(1)[3:101])
ydiss_cht=purge(sh3.col_values(6)[3:101])
plot(y_cht,ydiss_cht/tt_cht,'-+',color='g')

cht =  xlrd.open_workbook('./g2a2.xls')
sh3 = cht.sheet_by_name('budget_tt')
y_cht=purge(sh3.col_values(0)[3:101])
y_sol=purge(sh3.col_values(8)[3:135])
tt_cht=purge(sh3.col_values(1)[3:101])
ydiss_cht=purge(sh3.col_values(6)[3:101])
plot(y_cht,ydiss_cht/tt_cht,'o',color='k',markerfacecolor='none',markeredgecolor='k',label=r'$G=2$')
plot(y_cht,ydiss_cht/tt_cht,'o',color='g',markerfacecolor='none',markeredgecolor='r')#,label=r'$G=2$')
plot(y_cht,ydiss_cht/tt_cht,'--',color='r',label=r'$G_2=\frac{1}{2}$')
plot(y_cht,ydiss_cht/tt_cht,'--o',color='r',markerfacecolor='none',markeredgecolor='r')

#annotate('$Increasing\mbox{ }K$', xy=(0.6,4.3), xytext=(0.6,0.25),arrowprops=dict(facecolor='black', shrink=0.05),horizontalalignment='center',verticalalignment='center')
plot(y_refd,ydiss_refd/tt_refd,'D',color='k',label=r'$isoT$',markerfacecolor='none',markeredgecolor='k')

#text(-0.4, 0.12, r'$Solid$', horizontalalignment='center',verticalalignment='center',color='k')
#text(0.4, 0.12, r'$Fluid$', horizontalalignment='center',verticalalignment='center',color='k')

#legend(bbox_to_anchor=(1.45,0.85),numpoints=1)
savefig("all_ydiss_K2.png",bbox_inches='tight')
savefig("all_ydiss_K2.pdf",bbox_inches='tight')
