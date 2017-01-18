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


#ref =  xlrd.open_workbook('ref_scal_diric.xls')
ref = xlrd.open_workbook('dirichlet.xls')
refn= xlrd.open_workbook('neumann.xls')
inc =  xlrd.open_workbook('robin.xls')
cht =  xlrd.open_workbook('conju_ref.xls')

#
def purge(L):
	return [x for x in L if x <> '']

# open sheet in excel
sh1 = inc.sheet_by_name('budget_ut')
#sh2 = ref.sheet_by_name('diric')
sh2 = ref.sheet_by_name('budget_ut')
sh2n= refn.sheet_by_name('budget_ut')
sh3 = cht.sheet_by_name('budget_ut')
# y-coord
y_inc_bot=purge(sh1.col_values(0)[3:99])
#y_ref=purge(sh2.col_values(1)[1:73])
y_ref=purge(sh2.col_values(0)[3:99])
y_refn=purge(sh2n.col_values(0)[3:99])
y_cht=purge(sh3.col_values(0)[3:99])
# look for data
diss_bot=purge(sh1.col_values(1)[3:99])
prod_bot=purge(sh1.col_values(2)[3:99])
corp_bot=purge(sh1.col_values(3)[3:99])
tbdf_bot=purge(sh1.col_values(4)[3:99])
vsdf_bot=purge(sh1.col_values(5)[3:99])
#prod_ref=purge(sh2.col_values(2)[533:605])
#diss_ref=purge(sh2.col_values(3)[533:605])
#corp_ref=purge(sh2.col_values(4)[533:605])
#tbdf_ref=purge(sh2.col_values(5)[533:605])
#vsdf_ref=purge(sh2.col_values(6)[533:605])
diss_ref=purge(sh2.col_values(1)[3:99])
prod_ref=purge(sh2.col_values(2)[3:99])
corp_ref=purge(sh2.col_values(3)[3:99])
tbdf_ref=purge(sh2.col_values(4)[3:99])
vsdf_ref=purge(sh2.col_values(5)[3:99])
diss_refn=purge(sh2n.col_values(1)[3:99])
prod_refn=purge(sh2n.col_values(2)[3:99])
corp_refn=purge(sh2n.col_values(3)[3:99])
tbdf_refn=purge(sh2n.col_values(4)[3:99])
vsdf_refn=purge(sh2n.col_values(5)[3:99])
diss_cht=purge(sh3.col_values(1)[3:99])
prod_cht=purge(sh3.col_values(2)[3:99])
corp_cht=purge(sh3.col_values(3)[3:99])
tbdf_cht=purge(sh3.col_values(4)[3:99])
vsdf_cht=purge(sh3.col_values(5)[3:99])
# plot velocity
plot(y_cht,diss_cht,'-',color='r',label=r'$Diss$')
plot(y_cht,prod_cht,'-',color='g',label=r'$Prod$')
plot(y_cht,corp_cht,'-.',color='b',label=r'$P-Cor$')
plot(y_cht,tbdf_cht,':',color='k',label=r'$Tb.Df.$')
plot(y_cht,vsdf_cht,'--',color='c',label=r'$Vs.Df.$')
plot(y_inc_bot,diss_bot,'-o',color='r',markerfacecolor='none',markeredgecolor='r')
plot(y_inc_bot,prod_bot,'-o',color='g',markerfacecolor='none',markeredgecolor='g')
plot(y_inc_bot,corp_bot,'-.o',color='b',markerfacecolor='none',markeredgecolor='b')
plot(y_inc_bot,tbdf_bot,':o',color='k',markerfacecolor='none',markeredgecolor='k')
plot(y_inc_bot,vsdf_bot,'--o',color='c',markerfacecolor='none',markeredgecolor='c')
plot(y_ref,diss_ref,'-+',color='r')
plot(y_ref,prod_ref,'-+',color='g')
plot(y_ref,corp_ref,'-.+',color='b')
plot(y_ref,tbdf_ref,':+',color='k')
plot(y_ref,vsdf_ref,'--+',color='c')
plot(y_refn,diss_refn,'-x',color='r')
plot(y_refn,prod_refn,'-x',color='g')
plot(y_refn,corp_refn,'-.x',color='b')
plot(y_refn,tbdf_refn,':x',color='k')
plot(y_refn,vsdf_refn,'--x',color='c')

#Graph settings
xscale('log')
yscale('linear')
axis([0.3,150,-0.3,0.5])
xlabel(r"$y^+$")
ylabel(r"Budget $<u'T'>$")
legend(bbox_to_anchor=(1,1),numpoints=1)

savefig("all_bud_utcht.png",bbox_inches='tight')
savefig("all_bud_utcht.pdf",bbox_inches='tight')
