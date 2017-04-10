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
sh2d= refd.sheet_by_name('moyenne')
sh2 = refn.sheet_by_name('moyenne')

dudy_ref =purge( sh2.col_values(13)[6:101])
dudy_refd=purge(sh2d.col_values(13)[6:101])

dtdy_ref =purge( sh2.col_values(14)[6:101])
dtdy_refd=purge(sh2d.col_values(14)[6:101])
# open sheet in excel
sh2d= refd.sheet_by_name('fluctuations')
sh2 = refn.sheet_by_name('fluctuations')
# y-coord
y_ref=purge(sh2.col_values(9)[6:101])
y_refd=purge(sh2d.col_values(9)[6:101])
# look for data
uv_ref=purge(sh2.col_values(13)[6:101])
vt_ref=purge(sh2.col_values(15)[6:101])
uv_refd=purge(sh2d.col_values(13)[6:101])
vt_refd=purge(sh2d.col_values(15)[6:101])
# plot velocity
#plot(y_inc_bot,ut_bot,'-+',color='r',label=r'$uT$')
#plot(y_inc_bot,vt_bot,'--',color='r',label=r'$vT$')
#plot(y_inc_bot,tt_bot,'o',color='r',label=r'$Robin$',markerfacecolor='none',markeredgecolor='r')
#plot(y_ref,ut_ref,'-x',color='g')
#plot(y_ref,vt_ref,'--x',color='g')
plot(y_ref,dtdy_ref*uv_ref/(dudy_ref*vt_ref),'s',color='k',label=r'$isoQ$',markerfacecolor='none',markeredgecolor='k')
#plot(y_refd,ut_refd,'-+',color='b')
#plot(y_refd,vt_refd,'--+',color='b')

#Graph settings
#xscale('symlog',linthreshx=1.0)
xscale('linear')
yscale('linear')
axis([0.,10,0.1,1.1])
xlabel(r"$y^+$")
ylabel("$Pr_t$")

cht =  xlrd.open_workbook('./g05a2.xls')
sh3 = cht.sheet_by_name('moyenne')
dudy_cht=purge(sh3.col_values(13)[6:101])
dtdy_cht=purge(sh3.col_values(14)[6:101])
sh3 = cht.sheet_by_name('fluctuations')
y_cht=purge(sh3.col_values(9)[6:101])
uv_cht=purge(sh3.col_values(13)[6:101])
vt_cht=purge(sh3.col_values(15)[6:101])
plot(y_cht,dtdy_cht*uv_cht/(dudy_cht*vt_cht),'+',color='k',label=r'$G_1=\frac{1}{2}$')
plot(y_cht,dtdy_cht*uv_cht/(dudy_cht*vt_cht),'+',color='r')#,label=r'$G_1=\frac{1}{2}$')
plot(y_cht,dtdy_cht*uv_cht/(dudy_cht*vt_cht),'--+',color='r')
cht =  xlrd.open_workbook('./g05a1.xls')
sh3 = cht.sheet_by_name('moyenne')
dudy_cht=purge(sh3.col_values(13)[6:101])
dtdy_cht=purge(sh3.col_values(14)[6:101])
sh3 = cht.sheet_by_name('fluctuations')
y_cht=purge(sh3.col_values(9)[6:101])
uv_cht=purge(sh3.col_values(13)[6:101])
vt_cht=purge(sh3.col_values(15)[6:101])
plot(y_cht,dtdy_cht*uv_cht/(dudy_cht*vt_cht),'-+',color='g')
cht =  xlrd.open_workbook('./g05a05.xls')
sh3 = cht.sheet_by_name('moyenne')
dudy_cht=purge(sh3.col_values(13)[6:101])
dtdy_cht=purge(sh3.col_values(14)[6:101])
sh3 = cht.sheet_by_name('fluctuations')
y_cht=purge(sh3.col_values(9)[6:101])
uv_cht=purge(sh3.col_values(13)[6:101])
vt_cht=purge(sh3.col_values(15)[6:101])
plot(y_cht,dtdy_cht*uv_cht/(dudy_cht*vt_cht),':+',color='b')

cht =  xlrd.open_workbook('./g1a2.xls')
sh3 = cht.sheet_by_name('moyenne')
dudy_cht=purge(sh3.col_values(13)[6:101])
dtdy_cht=purge(sh3.col_values(14)[6:101])
sh3 = cht.sheet_by_name('fluctuations')
y_cht=purge(sh3.col_values(9)[6:101])
uv_cht=purge(sh3.col_values(13)[6:101])
vt_cht=purge(sh3.col_values(15)[6:101])
plot(y_cht,dtdy_cht*uv_cht/(dudy_cht*vt_cht),'x',color='k',label=r'$G_1=1$')
plot(y_cht,dtdy_cht*uv_cht/(dudy_cht*vt_cht),'x',color='r')#,label=r'$G_1=1$')
plot(y_cht,dtdy_cht*uv_cht/(dudy_cht*vt_cht),'--x',color='r')
cht =  xlrd.open_workbook('./g1a1_od4.xls')
sh3 = cht.sheet_by_name('moyenne')
dudy_cht=purge(sh3.col_values(13)[6:101])
dtdy_cht=purge(sh3.col_values(14)[6:101])
sh3 = cht.sheet_by_name('fluctuations')
y_cht=purge(sh3.col_values(9)[6:101])
uv_cht=purge(sh3.col_values(13)[6:101])
vt_cht=purge(sh3.col_values(15)[6:101])
plot(y_cht,dtdy_cht*uv_cht/(dudy_cht*vt_cht),'-x',color='g')
cht =  xlrd.open_workbook('./g1a05.xls')
sh3 = cht.sheet_by_name('moyenne')
dudy_cht=purge(sh3.col_values(13)[6:101])
dtdy_cht=purge(sh3.col_values(14)[6:101])
sh3 = cht.sheet_by_name('fluctuations')
y_cht=purge(sh3.col_values(9)[6:101])
uv_cht=purge(sh3.col_values(13)[6:101])
vt_cht=purge(sh3.col_values(15)[6:101])
plot(y_cht,dtdy_cht*uv_cht/(dudy_cht*vt_cht),':x',color='b')

cht =  xlrd.open_workbook('./g2a2.xls')
sh3 = cht.sheet_by_name('moyenne')
dudy_cht=purge(sh3.col_values(13)[6:101])
dtdy_cht=purge(sh3.col_values(14)[6:101])
sh3 = cht.sheet_by_name('fluctuations')
y_cht=purge(sh3.col_values(9)[6:101])
uv_cht=purge(sh3.col_values(13)[6:101])
vt_cht=purge(sh3.col_values(15)[6:101])
plot(y_cht,dtdy_cht*uv_cht/(dudy_cht*vt_cht),'o',color='k',markerfacecolor='none',markeredgecolor='k',label=r'$G_1=2$')
plot(y_cht,dtdy_cht*uv_cht/(dudy_cht*vt_cht),'o',color='g',markerfacecolor='none',markeredgecolor='r')#,label=r'$G_1=2$')
plot(y_cht,dtdy_cht*uv_cht/(dudy_cht*vt_cht),'--',color='r',label=r'$G_2=\frac{1}{2}$')
plot(y_cht,dtdy_cht*uv_cht/(dudy_cht*vt_cht),'--o',color='r',markerfacecolor='none',markeredgecolor='r')
cht =  xlrd.open_workbook('./g2a1.xls')
sh3 = cht.sheet_by_name('moyenne')
dudy_cht=purge(sh3.col_values(13)[6:101])
dtdy_cht=purge(sh3.col_values(14)[6:101])
sh3 = cht.sheet_by_name('fluctuations')
y_cht=purge(sh3.col_values(9)[6:101])
uv_cht=purge(sh3.col_values(13)[6:101])
vt_cht=purge(sh3.col_values(15)[6:101])
plot(y_cht,dtdy_cht*uv_cht/(dudy_cht*vt_cht),'-',color='g',label=r'$G_2=1$')
plot(y_cht,dtdy_cht*uv_cht/(dudy_cht*vt_cht),'-o',color='g',markerfacecolor='none',markeredgecolor='g')
cht =  xlrd.open_workbook('./g2a05.xls')
sh3 = cht.sheet_by_name('moyenne')
dudy_cht=purge(sh3.col_values(13)[6:101])
dtdy_cht=purge(sh3.col_values(14)[6:101])
sh3 = cht.sheet_by_name('fluctuations')
y_cht=purge(sh3.col_values(9)[6:101])
uv_cht=purge(sh3.col_values(13)[6:101])
vt_cht=purge(sh3.col_values(15)[6:101])
plot(y_cht,dtdy_cht*uv_cht/(dudy_cht*vt_cht),':',color='b',label=r'$G_2=2$')
plot(y_cht,dtdy_cht*uv_cht/(dudy_cht*vt_cht),':o',color='b',markerfacecolor='none',markeredgecolor='b')

#annotate('$Increasing\mbox{ }K$', xy=(0.6,4.3), xytext=(0.6,0.25),arrowprops=dict(facecolor='black', shrink=0.05),horizontalalignment='center',verticalalignment='center')
plot(y_refd,dtdy_refd*uv_refd/(dudy_refd*vt_refd),'D',color='k',label=r'$isoT$',markerfacecolor='none',markeredgecolor='k')

#text(-0.4, 0.12, r'$Solid$', horizontalalignment='center',verticalalignment='center',color='k')
#text(0.4, 0.12, r'$Fluid$', horizontalalignment='center',verticalalignment='center',color='k')

#legend(bbox_to_anchor=(1.45,0.85),numpoints=1)
savefig("all_prt.png",bbox_inches='tight')
savefig("all_prt.pdf",bbox_inches='tight')
