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


diric =  xlrd.open_workbook('./g1a1_convergence_solid_domain.xlsx')

#
def purge(L):
	return [x for x in L if x <> '']

# open sheet in excel
diri = diric.sheet_by_name('Sheet1')
# y-coord
x=array(purge(diri.col_values(1)[0:0]))#:132]))

y = [100000*0.002*(iii-1) for iii in range(1, 15)]
t = [diri.col_values(2+3*i)[0] for i in range(0, 14)]

# plot velocity
plot(y,t,'-+',color='k')

#Graph settings
xscale('linear')
yscale('linear')
ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
axis([0.,2700,0.00000064,0.00000101])##axis([-1,-0.8,0.00000064,0.00000113])
xlabel(r"$t$")#xlabel(r"$y$")
ylabel("$\overline{T'^2}$")

savefig("convergence_t2_solide.png",bbox_inches='tight')
savefig("convergence_t2_solide.pdf",bbox_inches='tight')
