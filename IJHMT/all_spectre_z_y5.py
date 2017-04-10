##################################################################
# Flageul et al.
# http://doi.org/10.1016/j.ijheatmasstransfer.2017.04.005
# Reads xls file and generates figures
##################################################################
from numpy import *
from pylab import *
from scipy import fftpack
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


inc =  xlrd.open_workbook('./autocorr_z.xls')

#
def purge(L):
	return [x for x in L if x <> '']

# open sheet in excel
sh1 = inc.sheet_by_name('T5')
# y-coord
x=array(purge(sh1.col_values(0)[2:132]))
x=x*2*math.acos(-1)/(x[128]*x[1])
# look for data
diss_refn=purge(sh1.col_values(2)[3:132])
spec_refn=fftpack.dct(diss_refn)
diss_ref=purge(sh1.col_values(1)[3:132])
spec_ref=fftpack.dct(diss_ref)

# plot velocity
plot(x[1:],pow(x[1:],5/3)*spec_refn,'s',color='k',markerfacecolor='none',markeredgecolor='k')
plot(x[1:],pow(x[1:],5/3)*spec_ref,'D',color='k',markerfacecolor='none',markeredgecolor='k')

diss_cht=purge(sh1.col_values(6)[3:132])
spec_cht=fftpack.dct(diss_cht)
plot(x[1:],pow(x[1:],5/3)*spec_cht,':+',color='b')
diss_cht=purge(sh1.col_values(7)[3:132])
spec_cht=fftpack.dct(diss_cht)
plot(x[1:],pow(x[1:],5/3)*spec_cht,'-+',color='g')
diss_cht=purge(sh1.col_values(8)[3:132])
spec_cht=fftpack.dct(diss_cht)
plot(x[1:],pow(x[1:],5/3)*spec_cht,'--+',color='r')

diss_cht=purge(sh1.col_values(9)[3:132])
spec_cht=fftpack.dct(diss_cht)
plot(x[1:],pow(x[1:],5/3)*spec_cht,':x',color='b')
diss_cht=purge(sh1.col_values(10)[3:132])
spec_cht=fftpack.dct(diss_cht)
plot(x[1:],pow(x[1:],5/3)*spec_cht,'-x',color='g')
diss_cht=purge(sh1.col_values(11)[3:132])
spec_cht=fftpack.dct(diss_cht)
plot(x[1:],pow(x[1:],5/3)*spec_cht,'--x',color='r')

diss_cht=purge(sh1.col_values(12)[3:132])
spec_cht=fftpack.dct(diss_cht)
plot(x[1:],pow(x[1:],5/3)*spec_cht,':o',color='b',markerfacecolor='none',markeredgecolor='b')
diss_cht=purge(sh1.col_values(13)[3:132])
spec_cht=fftpack.dct(diss_cht)
plot(x[1:],pow(x[1:],5/3)*spec_cht,'-o',color='g',markerfacecolor='none',markeredgecolor='g')
diss_cht=purge(sh1.col_values(14)[3:132])
spec_cht=fftpack.dct(diss_cht)
plot(x[1:],pow(x[1:],5/3)*spec_cht,'--o',color='r',markerfacecolor='none',markeredgecolor='r')

#Graph settings
xscale('linear')
yscale('log')
axis([0.002,1.4,0.0002,2])
xlabel(r"${k_z}^+$")
ylabel(r"Compensated spanwise spectra ($T'$)")

savefig("all_spectre_z_y5.png",bbox_inches='tight')
savefig("all_spectre_z_y5.pdf",bbox_inches='tight')
