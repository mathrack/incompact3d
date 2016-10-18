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


diric =  xlrd.open_workbook('./diric_flu_jj_2_nx_256_nz_256_nt_000075000.xls')
neuma =  xlrd.open_workbook('./neuma_flu_jj_2_nx_256_nz_256_nt_000075000.xls')
robin =  xlrd.open_workbook('./robin_flu_jj_2_nx_256_nz_256_nt_000075000.xls')
conju =  xlrd.open_workbook('./conju_flu_jj_2_nx_512_nz_512_nt_000075000.xlsx')
#
def purge(L):
	return [x for x in L if x <> '']
# open sheet in excel
diri = diric.sheet_by_name('Feuille1')
neum = neuma.sheet_by_name('Feuille1')
robi = robin.sheet_by_name('Feuille1')
conj = conju.sheet_by_name('Feuille1')
# y-coord
x=array(purge(diri.col_values(0)[1:129]))
z=array(purge(diri.row_values(0)[1:129]))
X, Z = meshgrid(x, z)
xx=array(purge(conj.col_values(0)[1:257]))
zz=array(purge(conj.row_values(0)[1:257]))
XX, ZZ = meshgrid(xx, zz)
# look for data
h_diri=zeros((128,128))
for ii in range(1,128):
	h_diri[ii-1]=array(purge(diri.col_values(ii)[1:129]))
h_neum=zeros((128,128))
for ii in range(1,128):
	h_neum[ii-1]=array(purge(neum.col_values(ii)[1:129]))
h_robi=zeros((128,128))
for ii in range(1,128):
	h_robi[ii-1]=array(purge(robi.col_values(ii)[1:129]))
h_conj=zeros((256,256))
for ii in range(1,256):
	h_conj[ii-1]=array(purge(conj.col_values(ii)[1:257]))

# plot velocity
axis([0.,1000,0.,90])
ylabel(r"$z^+$")
xlabel(r"$x^+$")
CS = contour(X,Z,h_diri,levels = [0.1,-0.1],colors='r')
text(390,10,r"$isoT$",color='r')
text(180,70,r"$isoT$",color='r')
#CS = contour(X,Z,h_neum,levels = [0.1,-0.1],colors='g')
#text(85,7,r"$isoQ$",color='g')
#text(5,61,r"$isoQ$",color='g')
CS = contour(X,Z,h_robi,levels = [0.1,-0.1],colors='b')
text(650,20,r"$Robin$",color='b')
text(315,70,r"$Robin$",color='b')
CS = contour(XX,ZZ,h_conj,levels = [0.1,-0.1],colors='k')
text(220,10,r"$Conju.$",color='k')
text(10,65,r"$Conju.$",color='k')

savefig("autocorr_Q_y0.png",bbox_inches='tight')
savefig("autocorr_Q_y0.pdf",bbox_inches='tight')
