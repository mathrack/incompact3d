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
	return [x for x in L if x <> '']

# open sheet in excel
sh2d= refd.sheet_by_name('fluctuations')
sh2 = refn.sheet_by_name('fluctuations')
# y-coord
#y_inc_bot=purge(sh1.col_values(9)[5:101])
y_ref=purge(sh2.col_values(9)[5:101])
y_refd=purge(sh2d.col_values(9)[5:101])
# look for data
#ut_bot=purge(sh1.col_values(14)[5:101])
#vt_bot=purge(sh1.col_values(15)[5:101])
#tt_bot=purge(sh1.col_values(16)[5:101])
ut_ref=purge(sh2.col_values(14)[5:101])
vt_ref=purge(sh2.col_values(15)[5:101])
tt_ref=purge(sh2.col_values(16)[5:101])
ut_refd=purge(sh2d.col_values(14)[5:101])
vt_refd=purge(sh2d.col_values(15)[5:101])
tt_refd=purge(sh2d.col_values(16)[5:101])
# plot velocity
#plot(y_inc_bot,ut_bot,'-+',color='r',label=r'$uT$')
#plot(y_inc_bot,vt_bot,'--',color='r',label=r'$vT$')
#plot(y_inc_bot,tt_bot,'o',color='r',label=r'$Robin$',markerfacecolor='none',markeredgecolor='r')
#plot(y_ref,ut_ref,'-x',color='g')
#plot(y_ref,vt_ref,'--x',color='g')
plot(y_ref,tt_ref,'s',color='k',label=r'$isoQ$',markerfacecolor='none',markeredgecolor='k')
#plot(y_refd,ut_refd,'-+',color='b')
#plot(y_refd,vt_refd,'--+',color='b')

#Graph settings
xscale('symlog',linthreshx=1.0)
#xscale('linear')
yscale('log')
axis([-10,10,0.1,10])
xlabel(r"$y^+$")
ylabel("$\overline{T'^2}$")

cht =  xlrd.open_workbook('./g05a2.xls')
sh3 = cht.sheet_by_name('fluctuations')
y_cht=purge(sh3.col_values(9)[5:101])
y_sol=purge(sh3.col_values(28)[5:135])
ut_cht=purge(sh3.col_values(14)[5:101])
vt_cht=purge(sh3.col_values(15)[5:101])
tt_cht=purge(sh3.col_values(16)[5:101])
tt_sol=purge(sh3.col_values(29)[5:135])
plot(y_cht,tt_cht,'+',color='k',label=r'$G=\frac{1}{2}$')
plot(y_cht,tt_cht,'+',color='r')#,label=r'$G=\frac{1}{2}$')
plot(y_cht,tt_cht,'--+',color='r')
plot(y_sol,tt_sol,'--+',color='r')
cht =  xlrd.open_workbook('./g05a1.xls')
sh3 = cht.sheet_by_name('fluctuations')
y_cht=purge(sh3.col_values(9)[5:101])
y_sol=purge(sh3.col_values(28)[5:135])
ut_cht=purge(sh3.col_values(14)[5:101])
vt_cht=purge(sh3.col_values(15)[5:101])
tt_cht=purge(sh3.col_values(16)[5:101])
tt_sol=purge(sh3.col_values(29)[5:135])
plot(y_cht,tt_cht,'-+',color='g')
plot(y_sol,tt_sol,'-+',color='g')
cht =  xlrd.open_workbook('./g05a05.xls')
sh3 = cht.sheet_by_name('fluctuations')
y_cht=purge(sh3.col_values(9)[5:101])
y_sol=purge(sh3.col_values(28)[5:135])
ut_cht=purge(sh3.col_values(14)[5:101])
vt_cht=purge(sh3.col_values(15)[5:101])
tt_cht=purge(sh3.col_values(16)[5:101])
tt_sol=purge(sh3.col_values(29)[5:135])
plot(y_cht,tt_cht,':+',color='b')
plot(y_sol,tt_sol,':+',color='b')

cht =  xlrd.open_workbook('./g1a2.xls')
sh3 = cht.sheet_by_name('fluctuations')
y_cht=purge(sh3.col_values(9)[5:101])
y_sol=purge(sh3.col_values(28)[5:135])
ut_cht=purge(sh3.col_values(14)[5:101])
vt_cht=purge(sh3.col_values(15)[5:101])
tt_cht=purge(sh3.col_values(16)[5:101])
tt_sol=purge(sh3.col_values(29)[5:135])
plot(y_cht,tt_cht,'x',color='k',label=r'$G=1$')
plot(y_cht,tt_cht,'x',color='r')#,label=r'$G=1$')
plot(y_cht,tt_cht,'--x',color='r')
plot(y_sol,tt_sol,'--x',color='r')
cht =  xlrd.open_workbook('./g1a1_od4.xls')
sh3 = cht.sheet_by_name('fluctuations')
y_cht=purge(sh3.col_values(9)[5:101])
y_sol=purge(sh3.col_values(28)[5:135])
ut_cht=purge(sh3.col_values(14)[5:101])
vt_cht=purge(sh3.col_values(15)[5:101])
tt_cht=purge(sh3.col_values(16)[5:101])
tt_sol=purge(sh3.col_values(29)[5:135])
plot(y_cht,tt_cht,'-x',color='g')
plot(y_sol,tt_sol,'-x',color='g')
cht =  xlrd.open_workbook('./g1a05.xls')
sh3 = cht.sheet_by_name('fluctuations')
y_cht=purge(sh3.col_values(9)[5:101])
y_sol=purge(sh3.col_values(28)[5:135])
ut_cht=purge(sh3.col_values(14)[5:101])
vt_cht=purge(sh3.col_values(15)[5:101])
tt_cht=purge(sh3.col_values(16)[5:101])
tt_sol=purge(sh3.col_values(29)[5:135])
plot(y_cht,tt_cht,':x',color='b')
plot(y_sol,tt_sol,':x',color='b')

cht =  xlrd.open_workbook('./g2a2.xls')
sh3 = cht.sheet_by_name('fluctuations')
y_cht=purge(sh3.col_values(9)[5:101])
y_sol=purge(sh3.col_values(28)[5:135])
ut_cht=purge(sh3.col_values(14)[5:101])
vt_cht=purge(sh3.col_values(15)[5:101])
tt_cht=purge(sh3.col_values(16)[5:101])
tt_sol=purge(sh3.col_values(29)[5:135])
plot(y_cht,tt_cht,'o',color='k',markerfacecolor='none',markeredgecolor='k',label=r'$G=2$')
plot(y_cht,tt_cht,'o',color='g',markerfacecolor='none',markeredgecolor='r')#,label=r'$G=2$')
plot(y_cht,tt_cht,'--',color='r',label=r'$G_2=\frac{1}{2}$')
plot(y_cht,tt_cht,'--o',color='r',markerfacecolor='none',markeredgecolor='r')
plot(y_sol,tt_sol,'--o',color='r',markerfacecolor='none',markeredgecolor='r')
cht =  xlrd.open_workbook('./g2a1.xls')
sh3 = cht.sheet_by_name('fluctuations')
y_cht=purge(sh3.col_values(9)[5:101])
y_sol=purge(sh3.col_values(28)[5:135])
ut_cht=purge(sh3.col_values(14)[5:101])
vt_cht=purge(sh3.col_values(15)[5:101])
tt_cht=purge(sh3.col_values(16)[5:101])
tt_sol=purge(sh3.col_values(29)[5:135])
plot(y_cht,tt_cht,'-',color='g',label=r'$G_2=1$')
plot(y_cht,tt_cht,'-o',color='g',markerfacecolor='none',markeredgecolor='g')
plot(y_sol,tt_sol,'-o',color='g',markerfacecolor='none',markeredgecolor='g')
cht =  xlrd.open_workbook('./g2a05.xls')
sh3 = cht.sheet_by_name('fluctuations')
y_cht=purge(sh3.col_values(9)[5:101])
y_sol=purge(sh3.col_values(28)[5:135])
ut_cht=purge(sh3.col_values(14)[5:101])
vt_cht=purge(sh3.col_values(15)[5:101])
tt_cht=purge(sh3.col_values(16)[5:101])
tt_sol=purge(sh3.col_values(29)[5:135])
plot(y_cht,tt_cht,':',color='b',label=r'$G_2=2$')
plot(y_cht,tt_cht,':o',color='b',markerfacecolor='none',markeredgecolor='b')
plot(y_sol,tt_sol,':o',color='b',markerfacecolor='none',markeredgecolor='b')

#annotate('$Increasing\mbox{ }K$', xy=(0.6,4.3), xytext=(0.6,0.25),arrowprops=dict(facecolor='black', shrink=0.05),horizontalalignment='center',verticalalignment='center')
plot(y_refd,tt_refd,'D',color='k',label=r'$isoT$',markerfacecolor='none',markeredgecolor='k')

text(-0.4, 0.12, r'$Solid$', horizontalalignment='center',verticalalignment='center',color='k')
text(0.4, 0.12, r'$Fluid$', horizontalalignment='center',verticalalignment='center',color='k')

legend(bbox_to_anchor=(1.45,0.85),numpoints=1)
savefig("all_tt.png",bbox_inches='tight')
savefig("all_tt.pdf",bbox_inches='tight')