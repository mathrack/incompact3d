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

#
def purge(L):
	return [x for x in L if x <> '']

# open sheet in excel
sh1a = inc.sheet_by_name('moyenne')
sh2da = refd.sheet_by_name('moyenne')
sh2na = refn.sheet_by_name('moyenne')
# look for data
dudy_bot=array(purge(sh1a.col_values(13)[6:101]))
dtdy_bot=array(purge(sh1a.col_values(14)[6:101]))
dudy_refd=array(purge(sh2da.col_values(13)[6:101]))
dtdy_refd=array(purge(sh2da.col_values(14)[6:101]))
dudy_refn=array(purge(sh2na.col_values(13)[6:101]))
dtdy_refn=array(purge(sh2na.col_values(14)[6:101]))

# open sheet in excel
sh1 = inc.sheet_by_name('fluctuations')
sh2d = refd.sheet_by_name('fluctuations')
sh2n = refn.sheet_by_name('fluctuations')
# y-coord
y_inc_bot=array(purge(sh1.col_values(9)[6:101]))
y_inc_top=array(purge(sh1.col_values(17)[6:101]))
y_refd=array(purge(sh2d.col_values(9)[6:101]))
y_refn=array(purge(sh2n.col_values(9)[6:101]))
# look for data
uv_bot=array(purge(sh1.col_values(13)[6:101]))
vt_bot=array(purge(sh1.col_values(15)[6:101]))
uv_refd=array(purge(sh2d.col_values(13)[6:101]))
vt_refd=array(purge(sh2d.col_values(15)[6:101]))
uv_refn=array(purge(sh2n.col_values(13)[6:101]))
vt_refn=array(purge(sh2n.col_values(15)[6:101]))
# plot velocity
plot(y_inc_bot,uv_bot*dtdy_bot/(vt_bot*dudy_bot),'-',color='r',label=r'$Robin$')
plot(y_refn,uv_refn*dtdy_refn/(vt_refn*dudy_refn),'-x',color='g',label=r'$isoQ$')
plot(y_refd,uv_refd*dtdy_refd/(vt_refd*dudy_refd),'-+',color='b',label=r'$isoT$')

#plot(y_inc_bot,0.71*(1+(y_inc_bot**2)/(2*20)-exp(-y_inc_bot/sqrt(20)))/(1+(0.71*y_inc_bot**2)/(2*20)-(1)*exp(-y_inc_bot*sqrt(0.71/20))/(1+1)),'-.',color='r')
#plot(y_refd,0.71*(1+(y_inc_bot**2)/(2*20)-exp(-y_inc_bot/sqrt(20)))/(1+(0.71*y_inc_bot**2)/(2*20)-exp(-y_inc_bot*sqrt(0.71/20))),'-.',color='g')
#plot(y_refn,0.71*(1+(y_inc_bot**2)/(2*20)-exp(-y_inc_bot/sqrt(20)))/(1+(0.71*y_inc_bot**2)/(2*20)),'-.',color='b')

#Graph settings
xscale('log')
yscale('log')
axis([0.3,25,0.095,1.1])
xlabel(r"$y^+$")
#ylabel("$Pr_t$")
legend(bbox_to_anchor=(0.75,0.4),numpoints=1)

savefig("all_prandtl_turb_log.png",bbox_inches='tight')
savefig("all_prandtl_turb_log.pdf",bbox_inches='tight')

xscale('linear')
yscale('linear')
savefig("all_prandtl_turb_lin.png",bbox_inches='tight')
savefig("all_prandtl_turb_lin.pdf",bbox_inches='tight')

xscale('log')
yscale('linear')
savefig("all_prandtl_turb_linlog.png",bbox_inches='tight')
savefig("all_prandtl_turb_linlog.pdf",bbox_inches='tight')
