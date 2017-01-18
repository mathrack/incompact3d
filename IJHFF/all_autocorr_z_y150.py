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
sh1 = inc.sheet_by_name('robin')
sh2d= inc.sheet_by_name('diric')
sh2n= inc.sheet_by_name('neuma')
sh3 = inc.sheet_by_name('od6_flu')
# y-coord
y_inc_bot=purge(sh1.col_values(20)[2:132])
y_ref=purge(sh2d.col_values(20)[2:132])
y_refn=purge(sh2n.col_values(20)[2:132])
y_cht=purge(sh3.col_values(20)[2:132])
# look for data
diss_bot=purge(sh1.col_values(33)[2:132])

diss_ref=purge(sh2d.col_values(33)[2:132])

diss_refn=purge(sh2n.col_values(33)[2:132])

diss_cht=purge(sh3.col_values(33)[2:132])

# plot velocity
plot(y_cht,diss_cht,'-',color='r',label=r'$Conju$')

plot(y_inc_bot,diss_bot,'o',color='r',markerfacecolor='none',markeredgecolor='r',label=r'$Robin$')

plot(y_ref,diss_ref,'+',color='b',label=r'$isoT$')

plot(y_refn,diss_refn,'x',color='g',label=r'$isoQ$')

#Graph settings
xscale('linear')
yscale('linear')
axis([-0.1,200,-0.1,1.1])
xlabel(r"$z^+$")
ylabel(r"Spanwise autocorrelation ($T'$)")
legend(bbox_to_anchor=(1,1),numpoints=1)

savefig("all_autocorr_z_y150.png",bbox_inches='tight')
savefig("all_autocorr_z_y150.pdf",bbox_inches='tight')
