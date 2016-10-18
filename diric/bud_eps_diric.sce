nb_point=193;
nb_dns=73;
ly=2;
beta=0.225;
re=2280;
pr=0.71;

maxfiles(100);
format('e',20);
exec('/home/C91471/Téléchargements/plotlib/loader.sce')

y=zeros(nb_point,1);
ydns=zeros(nb_dns,1);

ume =zeros(y);
udns=zeros(ydns);

phime =zeros(y);
phidns=zeros(ydns);

uume =zeros(y);
uudns=zeros(ydns);
vvme =zeros(y);
vvdns=zeros(ydns);
wwme =zeros(y);
wwdns=zeros(ydns);
uvme =zeros(y);
uvdns=zeros(ydns);

uphime   =zeros(y);
uphidns  =zeros(ydns);
vphime   =zeros(y);
vphidns  =zeros(ydns);
phiphime =zeros(y);
phiphidns=zeros(ydns);

uuDissdns=zeros(ydns);
uuProddns=zeros(ydns);
uuPstrdns=zeros(ydns);
uuTbDfdns=zeros(ydns);
uuVsDfdns=zeros(ydns);
uuSumdns =zeros(ydns);

vvDissdns=zeros(ydns);
vvPstrdns=zeros(ydns);
vvPDifdns=zeros(ydns);
vvTbDfdns=zeros(ydns);
vvVsDfdns=zeros(ydns);
vvSumdns =zeros(ydns);

wwDissdns=zeros(ydns);
wwPstrdns=zeros(ydns);
wwTbDfdns=zeros(ydns);
wwVsDfdns=zeros(ydns);
wwSumdns =zeros(ydns);

uvDissdns=zeros(ydns);
uvProddns=zeros(ydns);
uvPstrdns=zeros(ydns);
uvPrDfdns=zeros(ydns);
uvTbDfdns=zeros(ydns);
uvVsDfdns=zeros(ydns);
uvSumdns =zeros(ydns);

utproddns=zeros(ydns);
utdissdns=zeros(ydns);
utpstrdns=zeros(ydns);
uttdifdns=zeros(ydns);
utmdifdns=zeros(ydns);
utsumdns =zeros(ydns);

vtproddns=zeros(ydns);
vtdissdns=zeros(ydns);
vtpstrdns=zeros(ydns);
vttdifdns=zeros(ydns);
vtmdifdns=zeros(ydns);
vtpdifdns=zeros(ydns);
vtsumdns =zeros(ydns);

ttproddns=zeros(ydns);
ttdissdns=zeros(ydns);
tttdifdns=zeros(ydns);
ttmdifdns=zeros(ydns);
ttsumdns =zeros(ydns);

cd /home/C91471/Documents/Bibliographie/Rapport3A/images/diric/raw

function [yy,var] = read_my_stat_1d(filename,nb)
  [fd_um, ier]=mopen(filename,'r');
  for i=1:nb
    tmp_m=mfscanf(fd_um,'%f %f'); 
    yy(i)=tmp_m(1);
    var(i)=tmp_m(2);
  end
  mclose(fd_um);
endfunction;

[y,ume]      = read_my_stat_1d('um1d.dat',nb_point);
[y,vme]      = read_my_stat_1d('vm1d.dat',nb_point);
[y,wme]      = read_my_stat_1d('wm1d.dat',nb_point);
[y,uume]     = read_my_stat_1d('uum1d.dat',nb_point);
[y,vvme]     = read_my_stat_1d('vvm1d.dat',nb_point);
[y,wwme]     = read_my_stat_1d('wwm1d.dat',nb_point);
[y,uvme]     = read_my_stat_1d('uvm1d.dat'    ,nb_point);
[y,phime]    = read_my_stat_1d('phim1d.dat'   ,nb_point);
[y,uphime]   = read_my_stat_1d('uphim1d.dat'  ,nb_point);
[y,vphime]   = read_my_stat_1d('vphim1d.dat'  ,nb_point);
[y,phiphime] = read_my_stat_1d('phiphim1d.dat',nb_point);

cd /home/C91471/Documents/Rapports/2014_10_27_IJHFF/images/reference
[fd_dns, ier]=mopen('dynamique.dat','r')
// thermique de Kasagi
//[fd_dns_t, ier]=mopen('thermique_flux.dat','r') // Flux thermique imposé
[fd_dns_t, ier]=mopen('thermique_temp.dat','r') // Température imposée
// thermique de Tiselj
[fd_xl, sst, sh_names, sh_pos] = xls_open('Re150Pr0.7.xls');
[numeric, text] = xls_read(fd_xl,sh_pos(3));

header_dns=mgetl(fd_dns,133);
for i=1:nb_dns
  tmp_dns=mfscanf(fd_dns,'%f %f %f %f %f %f %f');
  ydns(i)=tmp_dns(2);
  udns(i)=tmp_dns(3);
  uudns(i)=tmp_dns(4);
  vvdns(i)=tmp_dns(5);
  wwdns(i)=tmp_dns(6);
  ppdns(i)=tmp_dns(7);
end
uudns=uudns.*uudns;
vvdns=vvdns.*vvdns;
wwdns=wwdns.*wwdns;
ppdns=ppdns.*ppdns;
tmp_dns=mgetl(fd_dns,3);
for i=1:nb_dns
  tmp_dns=mfscanf(fd_dns,'%f %f %f %f %f');
  uvdns(i)=tmp_dns(3);
end
tmp_dns=mgetl(fd_dns,229);
for i=1:nb_dns
  tmp_dns=mfscanf(fd_dns,'%f %f %f %f %f %f %f');
  uuDissdns(i)=tmp_dns(3);
  uuProddns(i)=tmp_dns(4);
  uuPstrdns(i)=tmp_dns(5);
  uuTbDfdns(i)=tmp_dns(6);
  uuVsDfdns(i)=tmp_dns(7);
end
uuSumdns=uuDissdns+uuProddns+uuPstrdns+uuTbDfdns+uuVsDfdns;
tmp_dns=mgetl(fd_dns,4);
for i=1:nb_dns
  tmp_dns=mfscanf(fd_dns,'%f %f %f %f %f %f %f');
  vvDissdns(i)=tmp_dns(3);
  vvPstrdns(i)=tmp_dns(4);
  vvPDifdns(i)=tmp_dns(5);
  vvTbDfdns(i)=tmp_dns(6);
  vvVsDfdns(i)=tmp_dns(7);
end
vvSumdns=vvDissdns+vvPstrdns+vvPDifdns+vvTbDfdns+vvVsDfdns;
tmp_dns=mgetl(fd_dns,4);
for i=1:nb_dns
  tmp_dns=mfscanf(fd_dns,'%f %f %f %f %f %f');
  wwDissdns(i)=tmp_dns(3);
  wwPstrdns(i)=tmp_dns(4);
  wwTbDfdns(i)=tmp_dns(5);
  wwVsDfdns(i)=tmp_dns(6);
end
wwSumdns=wwDissdns+wwPstrdns+wwTbDfdns+wwVsDfdns;
tmp_dns=mgetl(fd_dns,4);
for i=1:nb_dns
  tmp_dns=mfscanf(fd_dns,'%f %f %f %f %f %f %f %f');
  uvDissdns(i)=tmp_dns(3);
  uvProddns(i)=tmp_dns(4);
  uvPstrdns(i)=tmp_dns(5);
  uvPrDfdns(i)=tmp_dns(6);
  uvTbDfdns(i)=tmp_dns(7);
  uvVsDfdns(i)=tmp_dns(8);
end
uvSumdns=uvDissdns+uvProddns+uvPstrdns+uvPrDfdns+uvTbDfdns+uvVsDfdns;
tmp_dns=mgetl(fd_dns,80);
for i=1:nb_dns
  tmp_dns=mfscanf(fd_dns,'%f %f %f %f %f %f %f %f %f %f');
  EPSprod1(i)=tmp_dns(3);
  EPSprod2(i)=tmp_dns(4);
  EPSprod3(i)=tmp_dns(5);
  EPSprod4(i)=tmp_dns(6);
  EPSpress(i)=tmp_dns(7);
  EPSturbu(i)=tmp_dns(8);
  EPSvisco(i)=tmp_dns(9);
  EPSdissi(i)=tmp_dns(10);
end


// température imposée
header_dns=mgetl(fd_dns_t,167);
for i=1:nb_dns
  tmp_dns=mfscanf(fd_dns_t,'%f %f %f %f %f %f %f %f');
  phidns(i)=tmp_dns(3);
  phiphidns(i)=tmp_dns(7);
end
tmp_dns=mgetl(fd_dns_t,3);
for i=1:nb_dns
  tmp_dns=mfscanf(fd_dns_t,'%f %f %f %f');
  uphidns(i)=tmp_dns(3);
  vphidns(i)=tmp_dns(4);
end
phiphidns=phiphidns.*phiphidns;
vphidns=-vphidns;
tmp_dns=mgetl(fd_dns_t,797-315+1);
for i=1:nb_dns
  tmp_dns=mfscanf(fd_dns_t,'%f %f %f %f %f %f');
  utproddns(i)=tmp_dns(2);
  utdissdns(i)=tmp_dns(3);
  utpstrdns(i)=tmp_dns(4);
  uttdifdns(i)=tmp_dns(5);
  utmdifdns(i)=tmp_dns(6);
end
utsumdns=utproddns+utdissdns+utpstrdns+uttdifdns+utmdifdns;
tmp_dns=mgetl(fd_dns_t,6);
for i=1:nb_dns
  tmp_dns=mfscanf(fd_dns_t,'%f %f %f %f %f %f %f');
  vtproddns(i)=tmp_dns(2);
  vtdissdns(i)=tmp_dns(3);
  vtpstrdns(i)=tmp_dns(4);
  vttdifdns(i)=tmp_dns(5);
  vtmdifdns(i)=tmp_dns(6);
  vtpdifdns(i)=tmp_dns(7);
end
vtsumdns=vtproddns+vtdissdns+vtpstrdns+vttdifdns+vtmdifdns+vtpdifdns;
tmp_dns=mgetl(fd_dns_t,5);
for i=1:nb_dns
  tmp_dns=mfscanf(fd_dns_t,'%f %f %f %f %f');
  ttproddns(i)=tmp_dns(2);
  ttdissdns(i)=tmp_dns(3);
  tttdifdns(i)=tmp_dns(4);
  ttmdifdns(i)=tmp_dns(5);
end
ttsumdns=ttproddns+ttdissdns+tttdifdns+ttmdifdns;
tmp_dns=mgetl(fd_dns_t,5);
for i=1:nb_dns
  tmp_dns=mfscanf(fd_dns_t,'%f %f %f %f %f %f %f'); // Mg_Prod        T_Prod         Diss          T_Trans 	G_Prod         M_Diff
  EPST_Pr12(i)=tmp_dns(2);
  EPST_Pro4(i)=tmp_dns(3);
  EPST_Diss(i)=tmp_dns(4);
  EPST_TbDf(i)=tmp_dns(5);
  EPST_P3(i)=tmp_dns(6);
  EPST_Diff(i)=tmp_dns(7);
end

// flux imposé
//header_dns=mgetl(fd_dns_t,144);
//for i=1:nb_dns
//  tmp_dns=mfscanf(fd_dns_t,'%f %f %f');
//  phidns(i)=tmp_dns(3);
//end
//tmp_dns=mgetl(fd_dns_t,3);
//for i=1:nb_dns
//  tmp_dns=mfscanf(fd_dns_t,'%f %f %f');
//  phiphidns(i)=tmp_dns(3);
//end
//tmp_dns=mgetl(fd_dns_t,444-292+1);
//for i=1:nb_dns
//  tmp_dns=mfscanf(fd_dns_t,'%f %f %f %f');
//  uphidns(i)=tmp_dns(3);
//end
//tmp_dns=mgetl(fd_dns_t,3);
//for i=1:nb_dns
//  tmp_dns=mfscanf(fd_dns_t,'%f %f %f %f %f %f %f');
//  vphidns(i)=tmp_dns(3);
//end
//phiphidns=phiphidns.*phiphidns;

ny=nb_point/2+1;
ym=zeros(ny,1);
yp=zeros(ny,1);
um=zeros(ny,1);
up=zeros(ny,1);
uum=zeros(ny,1);
uup=zeros(ny,1);
vvm=zeros(ny,1);
vvp=zeros(ny,1);
wwm=zeros(ny,1);
wwp=zeros(ny,1);
uvm=zeros(ny,1);
uvp=zeros(ny,1);
phim=zeros(ny,1);
phip=zeros(ny,1);
uphim=zeros(ny,1);
uphip=zeros(ny,1);
vphim=zeros(ny,1);
vphip=zeros(ny,1);
phiphim=zeros(ny,1);
phiphip=zeros(ny,1);

for i=1:ny
   um(i)=ume(i);
  uum(i)=uume(i);
  vvm(i)=vvme(i);
  wwm(i)=wwme(i);
  uvm(i)=uvme(i);

  phim(i)   = phime(i);
  uphim(i)  = uphime(i);
  vphim(i)  = vphime(i);
  phiphim(i)= phiphime(i);

   up(i)= ume(nb_point-i+1);
  uup(i)= uume(nb_point-i+1);
  vvp(i)= vvme(nb_point-i+1);
  wwp(i)= wwme(nb_point-i+1);
  uvp(i)= uvme(nb_point-i+1);

  phip(i)   = phime(nb_point-i+1);
  uphip(i)  = uphime(nb_point-i+1);
  vphip(i)  = vphime(nb_point-i+1);
  phiphip(i)= phiphime(nb_point-i+1);
end

function dudy=dery(ux,my_point,ly,beta)

pi=acos(-1);
yinf=-ly/2;
den=2*beta*yinf;
xnum=-yinf-sqrt(pi*pi*beta*beta+yinf*yinf);
alpha=abs(xnum/den);

alfa1y = 2;
alfa2y = 1/4;
alfajy = 1/3;
alfamy = 1/4;
alfany = 2;

dy=ly/(my_point-1);

af1y=-(5/2)/dy;
bf1y=     2/dy;
cf1y= (1/2)/dy;
af2y= (3/4)/dy;
afjy= (7/9)/dy;
bfjy=(1/36)/dy;
afmy= (3/4)/dy;
afny=-(5/2)/dy;
bfny=     2/dy;
cfny= (1/2)/dy;

ffy(1)         =alfa1y;
ffy(2)         =alfa2y;
ffy(my_point-2)=alfajy;
ffy(my_point-1)=alfamy;
ffy(my_point)  =0;
fcy(1)         =1;
fcy(2)         =1;
fcy(my_point-2)=1;
fcy(my_point-1)=1;
fcy(my_point  )=1;
fby(1)         =alfa2y;
fby(2)         =alfajy;
fby(my_point-2)=alfamy;
fby(my_point-1)=alfany;
fby(my_point  )=0;
for i=3:my_point-3
   ffy(i)=alfajy;
   fcy(i)=1;
   fby(i)=alfajy;
end
for i=1:my_point
   fwy(i)=fcy(i);
end
for i=2:my_point
   if (fby(i-1)==0) then
      fsy(i)=0;
   else
      fsy(i)=fby(i-1)/fwy(i-1);
   end
   fwy(i)=fwy(i)-ffy(i-1)*fsy(i);
end
for i=1:my_point
   fwy(i)=1/fwy(i);
end

yeta(1)=-1/2;
for i=2:my_point
   yeta(i)=(i-1)*(1/(my_point-1))-0.5;
end
for i=1:my_point
  ppy(i) = ly*( alpha/pi + (1/pi/beta) * sin(pi*yeta(i)) * sin(pi*yeta(i)) );
end

yp(1)=0;
for i=2:my_point
   den=2*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi);
   den1=sqrt(alpha*beta+1);
   xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi);
   den3=((sin(pi*yeta(i)))*(sin(pi*yeta(i)))/beta/pi)+alpha/pi;
   den4=2*alpha*beta-cos(2*pi*yeta(i))+1;
   xnum1=(atan(xnum*tan(pi*yeta(i))))*den4/den1/den3/den;
   cst=sqrt(beta)*pi/(2*sqrt(alpha)*sqrt(alpha*beta+1));
   if yeta(i)<0.5 then
      yp(i)=xnum1-cst+ly;
   end
   if yeta(i)==0.5 then
      yp(i)=ly;
   end
   if yeta(i)>0.5 then
      yp(i)=xnum1+cst+ly;
   end
end

dudy(1)=af1y*ux(1)+bf1y*ux(2)+cf1y*ux(3);
dudy(2)=af2y*(ux(3)-ux(1));
for i=3:my_point-2
   dudy(i)=afjy*(ux(i+1)-ux(i-1)) + bfjy*(ux(i+2)-ux(i-2));
end
dudy(my_point-1)=afmy*(ux(my_point)-ux(my_point-2));
dudy(my_point  )=-afny*ux(my_point)-bfny*ux(my_point-1)-cfny*ux(my_point-2);

for i=2:my_point
   dudy(i)=dudy(i)-dudy(i-1)*fsy(i);
end

dudy(my_point)=dudy(my_point)*fwy(my_point);
for i=my_point-1:-1:1
   dudy(i)=(dudy(i)-ffy(i)*dudy(i+1))*fwy(i);
end

for i=1:my_point
   dudy(i)=dudy(i)*ppy(i);
end

endfunction;

dudy=dery(ume,nb_point,ly,beta);

utau_m=sqrt( (abs(dudy(1)))/re );
utau_p=sqrt( (abs(dudy(nb_point)))/re );
utau=sqrt(0.5*(utau_m*utau_m+utau_p*utau_p))

for i=1:ny
  ym(i)=re*utau_m*y(i);
  yp(i)=re*utau_p*(2-y(nb_point-i+1));

  um(i)=um(i)/utau_m;
  up(i)=up(i)/utau_p;

  uum(i)=uum(i)/(utau_m*utau_m);
  uup(i)=uup(i)/(utau_p*utau_p);
  vvm(i)=vvm(i)/(utau_m*utau_m);
  vvp(i)=vvp(i)/(utau_p*utau_p);
  wwm(i)=wwm(i)/(utau_m*utau_m);
  wwp(i)=wwp(i)/(utau_p*utau_p);

  uvm(i)=uvm(i)/(utau_m*utau_m);
  uvp(i)=uvp(i)/(utau_p*utau_p);
end
dudy(1:ny)=dudy(1:ny)/abs(dudy(1));

dTdy=dery(phime,nb_point,ly,beta);
Ttau_m=(abs(dTdy(1)))/(re*pr*utau_m);
Ttau_p=(abs(dTdy(nb_point)))/(re*pr*utau_p);

for i=1:ny
  phim(i)=phim(i)/Ttau_m;
  phip(i)=phip(i)/Ttau_p;

  uphim(i)=uphim(i)/(Ttau_m*utau_m);
  uphip(i)=uphip(i)/(Ttau_p*utau_p);
  vphim(i)=vphim(i)/(Ttau_m*utau_m);
  vphip(i)=vphip(i)/(Ttau_p*utau_p);

  phiphim(i)=phiphim(i)/(Ttau_m*Ttau_m);
  phiphip(i)=phiphip(i)/(Ttau_p*Ttau_p);
end
dTdy(1:ny)=dTdy(1:ny)*pr/abs(dTdy(1));

cd /home/C91471/Documents/Bibliographie/Rapport3A/images/diric/raw

dudy=dery(ume,nb_point,ly,beta);
dTdy=dery(phime,nb_point,ly,beta);

re*utau

function dudy2=deryy(uy,my_point,ly,beta)

pi=acos(-1);
ordre6=%T;

// Ordre 6
  fpi2=4.;
// SVV : Ordre 4
  c1=exp( -( (1.-2./3.)/(0.3-2./3.) )**2 );
  nu0snu= 3.;
  kcdx2 = (nu0snu+1.)*pi*pi;
  kmdx2 = (c1*nu0snu+1.)*(4./9.)*pi*pi;

if ordre6 == %T then
  myalf=(45.*fpi2*pi*pi-272.)/(2.*(45.*fpi2*pi*pi-208.));
  mya  =((6.-9.*myalf)/4.);
  myb  =((-3.+24*myalf)/5.);
  myc  =((2.-11.*myalf)/20.);
else
  myalf = (64.*kmdx2-27.*kcdx2-96)                    / (64*kmdx2-54.*kcdx2+48);
  mya   = (54.*kcdx2-15.*kcdx2*kmdx2+12.)             / (64*kmdx2-54.*kcdx2+48);
  myb   = (192.*kmdx2-216.*kcdx2+24.*kcdx2*kmdx2-48.) / (64*kmdx2-54.*kcdx2+48);
  myc   = (54.*kcdx2-9.*kcdx2*kmdx2-108.)             / (64*kmdx2-54.*kcdx2+48);
end

dy=ly/(my_point-1);
dy2=dy*dy;

yinf=-ly/2;
den=2*beta*yinf;
xnum=-yinf-sqrt(pi*pi*beta*beta+yinf*yinf);
alpha=abs(xnum/den);

yeta(1)=-1/2;
for i=2:my_point
   yeta(i)=(i-1)*(1/(my_point-1))-0.5;
end
for i=1:my_point
  ppy(i) = ly*( alpha/pi + (1/pi/beta) * sin(pi*yeta(i)) * sin(pi*yeta(i)) );
  pp2y(i)=ppy(i)*ppy(i);
  pp4y(i)=( -2/beta * cos(pi*yeta(i)) * sin(pi*yeta(i)) );
end

alsa1y= 11.;
alsa2y= 1./10.;
alsa3y= 2./11.;
alsajy= myalf;
alsaty= 2./11.;
alsamy= 1./10.;
alsany= 11.;

as1y  = (13.    )/dy2;
bs1y  =-(27.    )/dy2;
cs1y  = (15.    )/dy2;
ds1y  =-(1.     )/dy2;
as2y  = (6./5.  )/dy2;
as3y  = (12./11.)/dy2;
bs3y  = (3./44. )/dy2;
asjy  = mya/dy2;
bsjy  = myb/(4.*dy2);
csjy  = myc/(9.*dy2);
asty  = (12./11.)/dy2;
bsty  = (3./44. )/dy2;
asmy  = (6./5.  )/dy2;
asny  = (13.    )/dy2;
bsny  =-(27.    )/dy2;
csny  = (15.    )/dy2;
dsny  =-(1.     )/dy2;

sfy(1)   =alsa1y;
sfy(2)   =alsa2y;
sfy(3)   =alsa3y;
sfy(my_point-3)=alsajy;
sfy(my_point-2)=alsaty;
sfy(my_point-1)=alsamy;
sfy(my_point)  =0.;
scy(1)   =1.;
scy(2)   =1.;
scy(3)   =1.;
scy(my_point-3)=1.;
scy(my_point-2)=1.;
scy(my_point-1)=1.;
scy(my_point  )=1.;
sby(1)   =alsa2y;
sby(2)   =alsa3y;
sby(3)   =alsajy;
sby(my_point-3)=alsaty;
sby(my_point-2)=alsamy;
sby(my_point-1)=alsany;
sby(my_point  )=0.;
for j=4:my_point-4
   sfy(j)=alsajy;
   scy(j)=1.;
   sby(j)=alsajy;
end

for i=1:my_point
   swy(i)=scy(i);
end
for i=2:my_point
   if (sby(i-1)==0) then
      ssy(i)=0;
   else
      ssy(i)=sby(i-1)/swy(i-1);
   end
   swy(i)=swy(i)-sfy(i-1)*ssy(i);
end
for i=1:my_point
   swy(i)=1./swy(i);
end

ty(1)=as1y*uy(1)+bs1y*uy(2)+cs1y*uy(3)+ds1y*uy(4);
ty(2)=as2y*(uy(3)-uy(2)-uy(2)+uy(1));
ty(3)=as3y*(uy(4)-uy(3)-uy(3)+uy(2))+bs3y*(uy(5)-uy(3)-uy(3)+uy(1));
for j=4:my_point-3
   ty(j)=asjy*(uy(j+1)-uy(j)-uy(j)+uy(j-1))+bsjy*(uy(j+2)-uy(j)-uy(j)+uy(j-2))+csjy*(uy(j+3)-uy(j)-uy(j)+uy(j-3));
end
ty(my_point-2)=asty*(uy(my_point-1)-uy(my_point-2)-uy(my_point-2)+uy(my_point-3))+bsty*(uy(my_point)-uy(my_point-2)-uy(my_point-2)+uy(my_point-4));
ty(my_point-1)=asmy*(uy(my_point  )-uy(my_point-1)-uy(my_point-1)+uy(my_point-2));
ty(my_point)=asny*uy(my_point)+bsny*uy(my_point-1)+csny*uy(my_point-2)+dsny*uy(my_point-3);

for j=2:my_point
   ty(j)=ty(j)-ty(j-1)*ssy(j);
end
ty(my_point)=ty(my_point)*swy(my_point);
for j=my_point-1:-1:1
   ty(j)=(ty(j)-sfy(j)*ty(j+1))*swy(j);
end

dudy=dery(uy,my_point,ly,beta);

dudy2=ty.*pp2y-pp4y.*dudy;

endfunction;

// P1
[y,dudyme] = read_my_stat_1d('dudym1d.dat',nb_point);
dudyme=dudyme/dudyme(1);
myepsxy=0.;
[y,tmp]=read_my_stat_1d('dudvdxm1d.dat',nb_point); myepsxy=myepsxy+tmp;
[y,tmp]=read_my_stat_1d('dudvdym1d.dat',nb_point); myepsxy=myepsxy+tmp;
[y,tmp]=read_my_stat_1d('dudvdzm1d.dat',nb_point); myepsxy=myepsxy+tmp;
myepsxy_m=myepsxy(1:ny)*(-2/((re**2)*(utau_m**4)));
scf(); plot(ydns,EPSprod1,'k',re*utau_m*y(1:ny),dudyme(1:ny).*myepsxy_m(1:ny),'b');
epsdyn_P1=dudyme(1:ny).*myepsxy_m(1:ny);


// P2
[y,dudyme] = read_my_stat_1d('dudym1d.dat',nb_point);
dudyme=dudyme/dudyme(1);
mytmp=0.;
[y,tmp]=read_my_stat_1d('dudxdudym1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dvdxdvdym1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dwdxdwdym1d.dat',nb_point); mytmp=mytmp+tmp;
mytmp_m=mytmp(1:ny)*(-2/((re**2)*(utau_m**4)));
scf(); plot(ydns,EPSprod2,'k',re*utau_m*y(1:ny),dudyme(1:ny).*mytmp_m(1:ny),'b');
epsdyn_P2=dudyme(1:ny).*mytmp_m(1:ny);

// P3
[y,dudyyme] = read_my_stat_1d('dudyym1d.dat',nb_point);
dudyyme=dudyyme/((re**2)*(utau_m**4));
[y,vdudyme] = read_my_stat_1d('vdudym1d.dat',nb_point);
vdudyme=vdudyme/(re*(utau_m**2));
scf(); plot(ydns,EPSprod3,'k',re*utau_m*y(1:ny),-2*dudyyme(1:ny).*vdudyme(1:ny),'b');
epsdyn_P3=-2*dudyyme(1:ny).*vdudyme(1:ny);

// P4
mytmp=0.;
[y,tmp]=read_my_stat_1d('dudx_dududxm1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dudx_dududym1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dudx_dududzm1d.dat',nb_point); mytmp=mytmp+tmp;

[y,tmp]=read_my_stat_1d('dudy_dudvdxm1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dudy_dudvdym1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dudy_dudvdzm1d.dat',nb_point); mytmp=mytmp+tmp;

[y,tmp]=read_my_stat_1d('dudz_dudwdxm1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dudz_dudwdym1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dudz_dudwdzm1d.dat',nb_point); mytmp=mytmp+tmp;

[y,tmp]=read_my_stat_1d('dvdx_dvdudxm1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dvdx_dvdudym1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dvdx_dvdudzm1d.dat',nb_point); mytmp=mytmp+tmp;

[y,tmp]=read_my_stat_1d('dvdy_dvdvdxm1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dvdy_dvdvdym1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dvdy_dvdvdzm1d.dat',nb_point); mytmp=mytmp+tmp;

[y,tmp]=read_my_stat_1d('dvdz_dvdwdxm1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dvdz_dvdwdym1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dvdz_dvdwdzm1d.dat',nb_point); mytmp=mytmp+tmp;

[y,tmp]=read_my_stat_1d('dwdx_dwdudxm1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dwdx_dwdudym1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dwdx_dwdudzm1d.dat',nb_point); mytmp=mytmp+tmp;

[y,tmp]=read_my_stat_1d('dwdy_dwdvdxm1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dwdy_dwdvdym1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dwdy_dwdvdzm1d.dat',nb_point); mytmp=mytmp+tmp;

[y,tmp]=read_my_stat_1d('dwdz_dwdwdxm1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dwdz_dwdwdym1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dwdz_dwdwdzm1d.dat',nb_point); mytmp=mytmp+tmp;

mytmp=-2*mytmp/((re**3)*(utau_m**6));
scf(); plot(ydns,EPSprod4,'k',re*utau_m*y(1:ny),mytmp(1:ny),'b');
epsdyn_P4=mytmp(1:ny);

// Turbulent transport
mytmp=0.;
[y,tmp]=read_my_stat_1d('vdudx2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('vdudy2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('vdudz2m1d.dat',nb_point); mytmp=mytmp+tmp;

[y,tmp]=read_my_stat_1d('vdvdx2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('vdvdy2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('vdvdz2m1d.dat',nb_point); mytmp=mytmp+tmp;

[y,tmp]=read_my_stat_1d('vdwdx2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('vdwdy2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('vdwdz2m1d.dat',nb_point); mytmp=mytmp+tmp;

mytmp2=dery(mytmp,nb_point,ly,beta);
mytmp=-mytmp2/((re**3)*(utau_m**6));
scf(); plot(ydns,EPSturbu,'k',re*utau_m*y(1:ny),mytmp(1:ny),'b');
epsdyn_TU=mytmp(1:ny);

// Pression
mytmp=0.;
[y,tmp]=read_my_stat_1d('dpdvdxm1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dpdvdym1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dpdvdzm1d.dat',nb_point); mytmp=mytmp+tmp;
mytmp2=dery(mytmp,nb_point,ly,beta);
mytmp=-mytmp2;
scf(); plot(ydns,EPSpress,'k',re*utau_m*y(1:ny),mytmp(1:ny),'b');
scf(); semilogx(ydns(2:nb_dns),EPSpress(2:nb_dns),'k',re*utau_m*y(2:ny),mytmp(2:ny),'b');
epsdyn_PR=mytmp(1:ny);

// Diffusion
mytmp=0.;
[y,tmp]=read_my_stat_1d('dudx2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dudy2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dudz2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dvdx2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dvdy2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dvdz2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dwdx2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dwdy2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dwdz2m1d.dat',nb_point); mytmp=mytmp+tmp;
mytmp2=deryy(mytmp,nb_point,ly,beta);
mytmp=mytmp2/((re**4)*(utau**6));
scf(); plot(ydns,EPSvisco,'k',re*utau_m*y(1:ny),mytmp(1:ny),'b');
epsdyn_DF=mytmp(1:ny);

// Dissipation
mytmp=0.;
[y,tmp]=read_my_stat_1d('dudxx2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dudyy2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dudzz2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dudxy2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dudxz2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dudyz2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dvdxx2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dvdyy2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dvdzz2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dvdxy2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dvdxz2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dvdyz2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dwdxx2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dwdyy2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dwdzz2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dwdxy2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dwdxz2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dwdyz2m1d.dat',nb_point); mytmp=mytmp+tmp;
mytmp=-2*mytmp/((re**4)*(utau**6));
scf(); plot(ydns,EPSdissi,'k',re*utau_m*y(1:ny),mytmp(1:ny),'b');
epsdyn_DS=mytmp(1:ny);

ym=re*utau_m*y(1:ny);
scf();
semilogx(ym(2:ny),epsdyn_P1(2:ny),'-r',ym(2:ny),epsdyn_P2(2:ny),'-g',ym(2:ny),epsdyn_P3(2:ny),'-b',ym(2:ny),epsdyn_P4(2:ny),'-k',ym(2:ny),epsdyn_TU(2:ny),'--r',ym(2:ny),epsdyn_PR(2:ny),'--g',ym(2:ny),epsdyn_DF(2:ny),'--b',ym(2:ny),epsdyn_DS(2:ny),'--k');
legend('$P1$','$P2$','$P3$','$P4$','$Turb.$','$Pres.$','$Diff.$','$Diss.$',4);
xlabel('$y+$');
ylabel('$Budget Epsilon$');

scf();
semilogx(ym(2:ny),epsdyn_P1(2:ny).*ym(2:ny),'-r',ym(2:ny),epsdyn_P2(2:ny).*ym(2:ny),'-g',ym(2:ny),epsdyn_P3(2:ny).*ym(2:ny),'-b',ym(2:ny),epsdyn_P4(2:ny).*ym(2:ny),'-k',ym(2:ny),epsdyn_TU(2:ny).*ym(2:ny),'--r',ym(2:ny),epsdyn_PR(2:ny).*ym(2:ny),'--g',ym(2:ny),epsdyn_DF(2:ny).*ym(2:ny),'--b',ym(2:ny),epsdyn_DS(2:ny).*ym(2:ny),'--k');
xlabel('$y+$');
ylabel('$Budget Epsilon$');

scf();
semilogx(ym(2:ny),epsdyn_P1(2:ny).*ym(2:ny).*ym(2:ny),'-r',ym(2:ny),epsdyn_P2(2:ny).*ym(2:ny).*ym(2:ny),'-g',ym(2:ny),epsdyn_P3(2:ny).*ym(2:ny).*ym(2:ny),'-b',ym(2:ny),epsdyn_P4(2:ny).*ym(2:ny).*ym(2:ny),'-k',ym(2:ny),epsdyn_TU(2:ny).*ym(2:ny).*ym(2:ny),'--r',ym(2:ny),epsdyn_PR(2:ny).*ym(2:ny).*ym(2:ny),'--g',ym(2:ny),epsdyn_DF(2:ny).*ym(2:ny).*ym(2:ny),'--b',ym(2:ny),epsdyn_DS(2:ny).*ym(2:ny).*ym(2:ny),'--k');
xlabel('$y+$');
ylabel('$Budget Epsilon$');

scf();
semilogx(ym(2:ny),epsdyn_P1(2:ny),'-r',ydns(2:nb_dns),EPSprod1(2:nb_dns),'o-r',ym(2:ny),epsdyn_P2(2:ny),'-g',ydns(2:nb_dns),EPSprod2(2:nb_dns),'o-g',ym(2:ny),epsdyn_P3(2:ny),'-b',ydns(2:nb_dns),EPSprod3(2:nb_dns),'o-b',ym(2:ny),epsdyn_P4(2:ny),'-k',ydns(2:nb_dns),EPSprod4(2:nb_dns),'o-k',ym(2:ny),epsdyn_TU(2:ny),'--r',ydns(2:nb_dns),EPSturbu(2:nb_dns),'o--r',ym(2:ny),epsdyn_PR(2:ny),'--g',ydns(2:nb_dns),EPSpress(2:nb_dns),'o--g',ym(2:ny),epsdyn_DF(2:ny),'--b',ydns(2:nb_dns),EPSvisco(2:nb_dns),'o--b',ym(2:ny),epsdyn_DS(2:ny),'--k',ydns(2:nb_dns),EPSdissi(2:nb_dns),'o--k');

write_csv([ym,epsdyn_P1,epsdyn_P2,epsdyn_P3,epsdyn_P4,epsdyn_TU,epsdyn_PR,epsdyn_DF,epsdyn_DS], "../csv/dynamique.csv"," ");

// Thermique
// Kasagi
// scf(); plot(ydns,EPST_Pr12,'-r',ydns,EPST_Pro4,'-g',ydns,EPST_Diss,'-b',ydns,EPST_TbDf,'-k',ydns,EPST_P3,'--r',ydns,EPST_Diff,'--g');

// Diffusion
mytmp=0.;
[y,tmp]=read_my_stat_1d('dphidx2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dphidy2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dphidz2m1d.dat',nb_point); mytmp=mytmp+tmp;
mytmp2=deryy(mytmp,nb_point,ly,beta);
mytmp=mytmp2/((pr**2)*(re**4)*(utau_m**4)*(Ttau_m**2));
scf(); plot(ydns,EPST_Diff,'-k',ym(1:ny),mytmp(1:ny),'-b');
scf(); semilogx(ydns(2:nb_dns),EPST_Diff(2:nb_dns),'-k',ym(2:ny),mytmp(2:ny),'-b');
epsteta_diff=mytmp(1:ny);

// Dissipation
mytmp=0.;
[y,tmp]=read_my_stat_1d('dtdxx2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dtdyy2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dtdzz2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dtdxy2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dtdxz2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dtdyz2m1d.dat',nb_point); mytmp=mytmp+tmp;
mytmp=-2*mytmp/((pr**2)*(re**4)*(utau_m**4)*(Ttau_m**2));
scf(); plot(ydns,EPST_Diss,'-k',ym(1:ny),mytmp(1:ny),'-b');
epsteta_diss=mytmp(1:ny);

// Turbulent transport
mytmp=0.;
[y,tmp]=read_my_stat_1d('vdtdx2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('vdtdy2m1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('vdtdz2m1d.dat',nb_point); mytmp=mytmp+tmp;
mytmp2=dery(mytmp,nb_point,ly,beta);
mytmp=-mytmp2/((pr**1)*(re**3)*(utau_m**4)*(Ttau_m**2));
scf(); plot(ydns,EPST_TbDf,'-k',ym(1:ny),mytmp(1:ny),'-b');
epsteta_turbu=mytmp(1:ny);

// P4
mytmp=0.;
[y,tmp]=read_my_stat_1d('dtdx_dudtdxm1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dtdx_dudtdym1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dtdx_dudtdzm1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dtdy_dvdtdxm1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dtdy_dvdtdym1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dtdy_dvdtdzm1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dtdz_dwdtdxm1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dtdz_dwdtdym1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dtdz_dwdtdzm1d.dat',nb_point); mytmp=mytmp+tmp;
mytmp=-2*mytmp/((pr**1)*(re**3)*(utau_m**4)*(Ttau_m**2));
scf(); plot(ydns,EPST_Pro4,'-k',ym(1:ny),mytmp(1:ny),'-b');
epsteta_p4=mytmp(1:ny);

// P3
mytmp=0.;
[y,tmp]=read_my_stat_1d('dphidyym1d.dat',nb_point); mytmp=mytmp+tmp; mytmp2=mytmp;
mytmp=0.;
[y,tmp]=read_my_stat_1d('vdtdym1d.dat',nb_point); mytmp=mytmp+tmp;
mytmp=mytmp.*mytmp2;
mytmp=-2*mytmp/((pr**1)*(re**3)*(utau_m**4)*(Ttau_m**2));
scf(); plot(ydns,EPST_P3,'-k',ym(1:ny),mytmp(1:ny),'-b');
epsteta_p3=mytmp(1:ny);

// P12
mytmp=0.;
[y,tmp]=read_my_stat_1d('dphidym1d.dat',nb_point); mytmp=mytmp+tmp; mytmp2=mytmp;
mytmp=0.;
[y,tmp]=read_my_stat_1d('dvdphidxm1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dvdphidym1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dvdphidzm1d.dat',nb_point); mytmp=mytmp+tmp;
mytmpP1=mytmp.*mytmp2;
mytmp=0.;
[y,tmp]=read_my_stat_1d('dudym1d.dat',nb_point); mytmp=mytmp+tmp; mytmp2=mytmp;
mytmp=0.;
[y,tmp]=read_my_stat_1d('dtdxdtdym1d.dat',nb_point); mytmp=mytmp+tmp;
mytmpP2=mytmp.*mytmp2;
mytmp=mytmpP1+mytmpP2;
mytmp=-2*mytmp/((pr**1)*(re**3)*(utau_m**4)*(Ttau_m**2));
scf(); plot(ydns,EPST_Pr12,'-k',ym(1:ny),mytmp(1:ny),'-b');
epsteta_p12=mytmp(1:ny);
epsteta_p1=-2*mytmpP1(1:ny)/((pr**1)*(re**3)*(utau_m**4)*(Ttau_m**2));
epsteta_p2=-2*mytmpP2(1:ny)/((pr**1)*(re**3)*(utau_m**4)*(Ttau_m**2));

// Source
mytmp=0.;
[y,tmp]=read_my_stat_1d('dudphidxm1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dudphidym1d.dat',nb_point); mytmp=mytmp+tmp;
[y,tmp]=read_my_stat_1d('dudphidzm1d.dat',nb_point); mytmp=mytmp+tmp;
mytmp=-2*mytmp/((pr**2)*(re**2)*(utau_m**3)*(Ttau_m**1));
epsteta_source=mytmp(1:ny);

scf(); plot(ym,epsteta_p12,'-k',ym,epsteta_p3,'-g',ym,epsteta_p4,'-b',ym,epsteta_turbu,'--k',ym,epsteta_diff,'--g',ym,epsteta_diss,'--b');
legend('$P1+P2$','$P3$','$P4$','$Turb.$','$Diff.$','$Diss.$',4);
xlabel('$y+$');
ylabel('$Budget Epsilon_\theta$');

scf(); semilogx(ym(2:ny),epsteta_p12(2:ny),'-k',ym(2:ny),epsteta_p3(2:ny),'-g',ym(2:ny),epsteta_p4(2:ny),'-b',ym(2:ny),epsteta_turbu(2:ny),'--k',ym(2:ny),epsteta_diff(2:ny),'--g',ym(2:ny),epsteta_diss(2:ny),'--b');
legend('$P1+P2$','$P3$','$P4$','$Turb.$','$Diff.$','$Diss.$',4);
xlabel('$y+$');
ylabel('$Budget Epsilon_\theta$');

scf(); semilogx(ym(2:ny),epsteta_p12(2:ny).*ym(2:ny),'-k',ym(2:ny),epsteta_p3(2:ny).*ym(2:ny),'-g',ym(2:ny),epsteta_p4(2:ny).*ym(2:ny),'-b',ym(2:ny),epsteta_turbu(2:ny).*ym(2:ny),'--k',ym(2:ny),epsteta_diff(2:ny).*ym(2:ny),'--g',ym(2:ny),epsteta_diss(2:ny).*ym(2:ny),'--b');
xlabel('$y+$');
ylabel('$Budget Epsilon_\theta$');

scf(); semilogx(ym(2:ny),epsteta_p12(2:ny).*ym(2:ny).*ym(2:ny),'-k',ym(2:ny),epsteta_p3(2:ny).*ym(2:ny).*ym(2:ny),'-g',ym(2:ny),epsteta_p4(2:ny).*ym(2:ny).*ym(2:ny),'-b',ym(2:ny),epsteta_turbu(2:ny).*ym(2:ny).*ym(2:ny),'--k',ym(2:ny),epsteta_diff(2:ny).*ym(2:ny).*ym(2:ny),'--g',ym(2:ny),epsteta_diss(2:ny).*ym(2:ny).*ym(2:ny),'--b');
xlabel('$y+$');
ylabel('$Budget Epsilon_\theta$');

scf(); semilogx(ym(2:ny),epsteta_p12(2:ny),'-k',ym(2:ny),epsteta_p3(2:ny),'-g',ym(2:ny),epsteta_p4(2:ny),'-b',ym(2:ny),epsteta_turbu(2:ny),'--k',ym(2:ny),epsteta_diff(2:ny),'--g',ym(2:ny),epsteta_diss(2:ny),'--b',ydns(2:nb_dns),EPST_Pr12(2:nb_dns),'-ok',ydns(2:nb_dns),EPST_P3(2:nb_dns),'-og',ydns(2:nb_dns),EPST_Pro4(2:nb_dns),'-ob',ydns(2:nb_dns),EPST_TbDf(2:nb_dns),'--ok',ydns(2:nb_dns),EPST_Diff(2:nb_dns),'--og',ydns(2:nb_dns),EPST_Diss(2:nb_dns),'--ob');

write_csv([ym,epsteta_p1,epsteta_p2,epsteta_p3,epsteta_p4,epsteta_turbu,epsteta_diff,epsteta_diss], "../csv/thermique.csv"," ");
