##################################################################
# Flageul et al.
# http://doi.org/10.1016/j.ijheatmasstransfer.2017.04.005
# Reads xls file and generates figures
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


inc =  xlrd.open_workbook('./autocorr_x.xls')

#
def purge(L):
	return [x for x in L if x <> '']

# open sheet in excel
sh1 = inc.sheet_by_name('Q5')
# y-coord
x=purge(sh1.col_values(0)[2:132])
# look for data
diss_refn=purge(sh1.col_values(2)[2:132])
diss_ref=purge(sh1.col_values(1)[2:132])

# plot velocity
plot(x,diss_refn,'s',color='k',label=r'$isoQ$',markerfacecolor='none',markeredgecolor='k')
plot(x,diss_ref,'D',color='k',label=r'$isoT$',markerfacecolor='none',markeredgecolor='k')

diss_cht=purge(sh1.col_values(6)[2:132])
plot(x,diss_cht,':+',color='b')
diss_cht=purge(sh1.col_values(7)[2:132])
plot(x,diss_cht,'-+',color='g')
diss_cht=purge(sh1.col_values(8)[2:132])
plot(x,diss_cht,'--+',color='r')

diss_cht=purge(sh1.col_values(9)[2:132])
plot(x,diss_cht,':x',color='b')
diss_cht=purge(sh1.col_values(10)[2:132])
plot(x,diss_cht,'-x',color='g')
diss_cht=purge(sh1.col_values(11)[2:132])
plot(x,diss_cht,'--x',color='r')

diss_cht=purge(sh1.col_values(12)[2:132])
plot(x,diss_cht,':o',color='b',markerfacecolor='none',markeredgecolor='b')
diss_cht=purge(sh1.col_values(13)[2:132])
plot(x,diss_cht,'-o',color='g',markerfacecolor='none',markeredgecolor='g')
diss_cht=purge(sh1.col_values(14)[2:132])
plot(x,diss_cht,'--o',color='r',markerfacecolor='none',markeredgecolor='r')

#Graph settings
xscale('linear')
yscale('linear')
axis([-0.1,2000,-0.1,1.1])
xlabel(r"$x^+$")
ylabel(r"Streamwise autocorrelation ($\partial_y T'$)")

savefig("all_autocorrQ_x_y5.png",bbox_inches='tight')
savefig("all_autocorrQ_x_y5.pdf",bbox_inches='tight')
