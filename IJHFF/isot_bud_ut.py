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


ref =  xlrd.open_workbook('ref_scal_diric.xls')
inc =  xlrd.open_workbook('dirichlet.xls')

#
def purge(L):
	return [x for x in L if x <> '']

# open sheet in excel
sh1 = inc.sheet_by_name('budget_ut')
sh2 = ref.sheet_by_name('diric')
# y-coord
y_inc_bot=purge(sh1.col_values(0)[3:99])
y_ref=purge(sh2.col_values(1)[1:73])
# look for data
diss_bot=purge(sh1.col_values(1)[3:99])
prod_bot=purge(sh1.col_values(2)[3:99])
corp_bot=purge(sh1.col_values(3)[3:99])
tbdf_bot=purge(sh1.col_values(4)[3:99])
vsdf_bot=purge(sh1.col_values(5)[3:99])
prod_ref=purge(sh2.col_values(2)[533:605])
diss_ref=purge(sh2.col_values(3)[533:605])
corp_ref=purge(sh2.col_values(4)[533:605])
tbdf_ref=purge(sh2.col_values(5)[533:605])
vsdf_ref=purge(sh2.col_values(6)[533:605])
# plot velocity
plot(y_inc_bot,diss_bot,'-',color='r',label=r'$Diss$')
plot(y_inc_bot,prod_bot,'-',color='g',label=r'$Prod$')
plot(y_inc_bot,corp_bot,'-.',color='b',label=r'$P-Cor$')
plot(y_inc_bot,tbdf_bot,':',color='k',label=r'$Tb.Df.$')
plot(y_inc_bot,vsdf_bot,'--',color='c',label=r'$Vs.Df.$')
plot(y_ref,diss_ref,'+',color='r')
plot(y_ref,prod_ref,'+',color='g')
plot(y_ref,corp_ref,'+',color='b')
plot(y_ref,tbdf_ref,'+',color='k')
plot(y_ref,vsdf_ref,'+',color='c')

#Graph settings
xscale('log')
yscale('linear')
axis([0.3,150,-0.3,0.5])
xlabel(r"$y^+$")
ylabel(r"Budget $<u'T'>$")
legend(bbox_to_anchor=(1,1),numpoints=1)

savefig("isot_bud_ut.png",bbox_inches='tight')
savefig("isot_bud_ut.pdf",bbox_inches='tight')
