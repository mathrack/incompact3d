nb_point=193;
nb_dns=73;
ly=2;
beta=0.225;
re=2280;
pr=0.71;
dt=0.002;

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

cd /home/C91471/Documents/Bibliographie/Rapport3A/images/robin_27/raw

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

write_csv([y,ume,vme,wme,phime], "../csv/raw_1.csv"," ");
write_csv([y,uume,vvme,wwme,uvme,uphime,vphime,phiphime], "../csv/raw_2.csv"," ");

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

cd /home/C91471/Documents/Bibliographie/Rapport3A/images/robin_27/raw

write_csv([ym,um,phim-phim(1),dudy(1:ny),dTdy(1:ny)], "../csv/moy1.csv"," ");
write_csv([yp,up,phip-phip(1)], "../csv/moy2.csv"," ");
write_csv([ym,uum,vvm,wwm,-uvm,uphim,-vphim,phiphim], "../csv/fluct1.csv"," ");
write_csv([yp,uup,vvp,wwp,uvp,uphip,vphip,phiphip], "../csv/fluct2.csv"," ");

dudy=dery(ume,nb_point,ly,beta);
dTdy=dery(phime,nb_point,ly,beta);

scf();
semilogx(ym(2:ny),um(2:ny),'+-r',yp(2:ny),up(2:ny),'+--g',ydns(2:nb_dns),udns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$<U_x>$');

scf();
semilogx(ym(2:ny),uum(2:ny),'+-r',yp(2:ny),uup(2:ny),'+--g',ydns(2:nb_dns),uudns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$<Rxx>$');

scf();
semilogx(ym(2:ny),vvm(2:ny),'+-r',yp(2:ny),vvp(2:ny),'+--g',ydns(2:nb_dns),vvdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$<Ryy>$');

scf();
semilogx(ym(2:ny),wwm(2:ny),'+-r',yp(2:ny),wwp(2:ny),'+--g',ydns(2:nb_dns),wwdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$<Rzz>$');

scf();
semilogx(ym(2:ny),-uvm(2:ny),'+-r',yp(2:ny),uvp(2:ny),'+--g',ydns(2:nb_dns),uvdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$-<Rxy>$');

scf();
//semilogx(ym(2:ny),phim(2:ny),'+-r',yp(2:ny),phip(2:ny),'+--g',numeric(3:50,1),numeric(3:50,7),'+k');
semilogx(ym(2:ny),phim(2:ny),'+-r',yp(2:ny),phip(2:ny)-phip(1),'+--g',ydns(2:nb_dns),phidns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',2);
xlabel('$y+$');
ylabel('$<phi>$');

scf();
//semilogx(ym(2:ny),uphim(2:ny),'+-r',yp(2:ny),uphip(2:ny),'+--g',numeric(3:50,1),numeric(3:50,9),'+k');
semilogx(ym(2:ny),uphim(2:ny),'+-r',yp(2:ny),uphip(2:ny),'+--g',ydns(2:nb_dns),uphidns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',2);
xlabel('$y+$');
ylabel('$<uphi>$');

scf();
//semilogx(ym(2:ny),-vphim(2:ny),'+-r',yp(2:ny),vphip(2:ny),'+--g',numeric(3:50,1),numeric(3:50,10),'+k');
semilogx(ym(2:ny),-vphim(2:ny),'+-r',yp(2:ny),vphip(2:ny),'+--g',ydns(2:nb_dns),vphidns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',2);
xlabel('$y+$');
ylabel('$<vphi>$');

scf();
//semilogx(ym(2:ny),phiphim(2:ny),'+-r',yp(2:ny),phiphip(2:ny),'+--g',numeric(3:50,1),numeric(3:50,8)**2,'+k');
semilogx(ym(2:ny),phiphim(2:ny),'+-r',yp(2:ny),phiphip(2:ny),'+--g',ydns(2:nb_dns),phiphidns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',2);
xlabel('$y+$');
ylabel('$<phiphi>$');

//scf(); 
//semilogx(ydns(2:nb_dns),uuDissdns(2:nb_dns),'+-r',ydns(2:nb_dns),uuProddns(2:nb_dns),'+-g',ydns(2:nb_dns),uuPstrdns(2:nb_dns),'+-k',ydns(2:nb_dns),uuTbDfdns(2:nb_dns),'+r',ydns(2:nb_dns),uuVsDfdns(2:nb_dns),'+g',ydns(2:nb_dns),uuSumdns(2:nb_dns),'ok'); 
//legend('$Diss$','$Prod$','$P-strain$','$Turb Diff$','$Visc Diff$','$Sum$',4);
//xlabel('$y+$');
//ylabel('$Budget <Rxx>$');

//scf(); 
//semilogx(ydns(2:nb_dns),vvDissdns(2:nb_dns),'+-r',ydns(2:nb_dns),vvPDifdns(2:nb_dns),'+-g',ydns(2:nb_dns),vvPstrdns(2:nb_dns),'+-k',ydns(2:nb_dns),vvTbDfdns(2:nb_dns),'+r',ydns(2:nb_dns),vvVsDfdns(2:nb_dns),'+g',ydns(2:nb_dns),vvSumdns(2:nb_dns),'ok'); 
//legend('$Diss$','$P-diff$','$P-strain$','$Turb Diff$','$Visc Diff$','$Sum$',4);
//xlabel('$y+$');
//ylabel('$Budget <Ryy>$');

//scf(); 
//semilogx(ydns(2:nb_dns),wwDissdns(2:nb_dns),'+-r',ydns(2:nb_dns),wwPstrdns(2:nb_dns),'+-k',ydns(2:nb_dns),wwTbDfdns(2:nb_dns),'+r',ydns(2:nb_dns),wwVsDfdns(2:nb_dns),'+g',ydns(2:nb_dns),wwSumdns(2:nb_dns),'ok'); 
//legend('$Diss$','$P-strain$','$Turb Diff$','$Visc Diff$','$Sum$',4);
//xlabel('$y+$');
//ylabel('$Budget <Rzz>$');

//scf(); 
//semilogx(ydns(2:nb_dns),uvDissdns(2:nb_dns),'+-r',ydns(2:nb_dns),uvProddns(2:nb_dns),'+-g',ydns(2:nb_dns),uvPstrdns(2:nb_dns),'+-k',ydns(2:nb_dns),uvPrDfdns(2:nb_dns),ydns(2:nb_dns),uvTbDfdns(2:nb_dns),'+r',ydns(2:nb_dns),uvVsDfdns(2:nb_dns),'+g',ydns(2:nb_dns),uvSumdns(2:nb_dns),'ok'); 
//legend('$Diss$','$Prod$','$P-strain$','$Pres Diff$','$Turb Diff$','$Visc Diff$','$Sum$',4);
//xlabel('$y+$');
//ylabel('$Budget <Rxy>$');

//scf(); 
//semilogx(ydns(2:nb_dns),utdissdns(2:nb_dns),'+-r',ydns(2:nb_dns),utproddns(2:nb_dns),'+-g',ydns(2:nb_dns),utpstrdns(2:nb_dns),'+-k',ydns(2:nb_dns),utmdifdns(2:nb_dns),'+r',ydns(2:nb_dns),uttdifdns(2:nb_dns),'+g',ydns(2:nb_dns),utsumdns(2:nb_dns),'ok'); 
//legend('$Diss$','$Prod$','$P-strain$','$M Diff$','$Turb Diff$','$Sum$',4);
//xlabel('$y+$');
//ylabel('$Budget <uT>$');

//scf(); 
//semilogx(ydns(2:nb_dns),vtdissdns(2:nb_dns),'+-r',ydns(2:nb_dns),vtproddns(2:nb_dns),'+-g',ydns(2:nb_dns),vtpstrdns(2:nb_dns),'+-k',ydns(2:nb_dns),vtmdifdns(2:nb_dns),'+r',ydns(2:nb_dns),vttdifdns(2:nb_dns),'+g',ydns(2:nb_dns),vtpdifdns(2:nb_dns),'+k',ydns(2:nb_dns),utsumdns(2:nb_dns),'ok'); 
//legend('$Diss$','$Prod$','$P-strain$','$M Diff$','$Turb Diff$','$Pres Diff$','$Sum$',4);
//xlabel('$y+$');
//ylabel('$Budget <vT>$');

//scf(); 
//semilogx(ydns(2:nb_dns),ttdissdns(2:nb_dns),'+-r',ydns(2:nb_dns),ttproddns(2:nb_dns),'+-g',ydns(2:nb_dns),ttmdifdns(2:nb_dns),'+r',ydns(2:nb_dns),tttdifdns(2:nb_dns),'+g',ydns(2:nb_dns),ttsumdns(2:nb_dns),'ok'); 
//legend('$Diss$','$Prod$','$M Diff$','$Turb Diff$','$Sum$',4);
//xlabel('$y+$');
//ylabel('$Budget <TT>$');

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

// Budget Rxx
// Production
uuprodm=-2*uvm.*dudy(1:ny)/(re*(utau_m**2));
for i=1:ny uuprodp(i)=-2*uvp(i).*dudy(nb_point-i+1)/(re*(utau_p**2)); end
scf();
semilogx(ym(2:ny),uuprodm(2:ny),'+-r',yp(2:ny),uuprodp(2:ny),'+--g',ydns(2:nb_dns),uuProddns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$Prod - Rxx$');

// Turbulent diffusion
[y,vuume] = read_my_stat_1d('vuum1d.dat',nb_point);
dvuudy=dery(vuume,nb_point,ly,beta);
uutbdfm=-dvuudy(1:ny)/(re*(utau_m**4));
for i=1:ny uutbdfp(i)=-dvuudy(nb_point-i+1)/(re*(utau_p**4)); end
scf();
semilogx(ym(2:ny),uutbdfm(2:ny),'+-r',yp(2:ny),uutbdfp(2:ny),'+--g',ydns(2:nb_dns),uuTbDfdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$TurbDiff - Rxx$');

// Viscous diffusion
drxxdyy=deryy(uume,nb_point,ly,beta);
uuvsdfm=drxxdyy(1:ny)/((re**2)*(utau_m**4));
for i=1:ny uuvsdfp(i)=drxxdyy(nb_point-i+1)/((re**2)*(utau_p**4)); end
scf();
semilogx(ym(2:ny),uuvsdfm(2:ny),'+-r',yp(2:ny),uuvsdfp(2:ny),'+--g',ydns(2:nb_dns),uuVsDfdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$ViscDiff - Rxx$');

// Corrélation vitesse - grad(p)
[y,udpdxme] = read_my_stat_1d('udpdxm1d.dat',nb_point);
uupstrm=-2*udpdxme(1:ny)/(dt*re*(utau_m**4));
for i=1:ny uupstrp(i)=-2*udpdxme(nb_point-i+1)/(dt*re*(utau_p**4)); end
scf();
semilogx(ym(2:ny),uupstrm(2:ny),'+-r',yp(2:ny),uupstrp(2:ny),'+--g',ydns(2:nb_dns),uuPstrdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$Pres-Str - Rxx$');

// Dissipation
[y,dudx2me] = read_my_stat_1d('dudx2m1d.dat',nb_point);
[y,dudy2me] = read_my_stat_1d('dudy2m1d.dat',nb_point);
[y,dudz2me] = read_my_stat_1d('dudz2m1d.dat',nb_point);
d2ume=dudx2me+dudy2me+dudz2me;
uudissm=-2*d2ume(1:ny)/((re**2)*(utau_m**4));
for i=1:ny uudissp(i)=-2*d2ume(nb_point-i+1)/((re**2)*(utau_p**4)); end
scf();
semilogx(ym(2:ny),uudissm(2:ny),'+-r',yp(2:ny),uudissp(2:ny),'+--g',ydns(2:nb_dns),uuDissdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$Dissipation - Rxx$');

uusumm=uudissm+uuprodm+uupstrm+uutbdfm+uuvsdfm;
write_csv([ym,uudissm,uuprodm,uupstrm,uutbdfm,uuvsdfm,uusumm], "../csv/uu1.csv"," ");
uusump=uudissp+uuprodp+uupstrp+uutbdfp+uuvsdfp;
write_csv([yp,uudissp,uuprodp,uupstrp,uutbdfp,uuvsdfp,uusump], "../csv/uu2.csv"," ");

scf();
semilogx(ym(2:ny),uudissm(2:ny),'+-r',ym(2:ny),uuprodm(2:ny),'+-g',ym(2:ny),uupstrm(2:ny),'+-k',ym(2:ny),uutbdfm(2:ny),'+r',ym(2:ny),uuvsdfm(2:ny),'+g',ym(2:ny),uusumm(2:ny),'ok');
legend('$Diss$','$Prod$','$P-strain$','$Turb Diff$','$Visc Diff$','$Sum$',4);
xlabel('$y+$');
ylabel('$Budget <Rxx>, bottom$');
scf();
semilogx(ym(2:ny),uudissm(2:ny),'+-r',ym(2:ny),uuprodm(2:ny),'+-g',ym(2:ny),uupstrm(2:ny),'+-k',ym(2:ny),uutbdfm(2:ny),'+r',ym(2:ny),uuvsdfm(2:ny),'+g',ydns(2:nb_dns),uuDissdns(2:nb_dns),'o-r',ydns(2:nb_dns),uuProddns(2:nb_dns),'o-g',ydns(2:nb_dns),uuPstrdns(2:nb_dns),'o-k',ydns(2:nb_dns),uuTbDfdns(2:nb_dns),'or',ydns(2:nb_dns),uuVsDfdns(2:nb_dns),'og');
xlabel('$y+$');
ylabel('$Budget <Rxx>, bottom$');
scf();
semilogx(yp(2:ny),uudissp(2:ny),'+-r',yp(2:ny),uuprodp(2:ny),'+-g',yp(2:ny),uupstrp(2:ny),'+-k',yp(2:ny),uutbdfp(2:ny),'+r',yp(2:ny),uuvsdfp(2:ny),'+g',ydns(2:nb_dns),uuDissdns(2:nb_dns),'o-r',ydns(2:nb_dns),uuProddns(2:nb_dns),'o-g',ydns(2:nb_dns),uuPstrdns(2:nb_dns),'o-k',ydns(2:nb_dns),uuTbDfdns(2:nb_dns),'or',ydns(2:nb_dns),uuVsDfdns(2:nb_dns),'og');
xlabel('$y+$');
ylabel('$Budget <Rxx>, top$');

// Budget Ryy
// Turbulent diffusion
[y,vvvme] = read_my_stat_1d('vvvm1d.dat',nb_point);
dvvvdy=dery(vvvme,nb_point,ly,beta);
vvtbdfm=-dvvvdy(1:ny)/(re*(utau_m**4));
for i=1:ny vvtbdfp(i)=-dvvvdy(nb_point-i+1)/(re*(utau_p**4)); end
scf();
semilogx(ym(2:ny),vvtbdfm(2:ny),'+-r',yp(2:ny),vvtbdfp(2:ny),'+--g',ydns(2:nb_dns),vvTbDfdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$TurbDiff - Ryy$');

// Viscous diffusion
dryydyy=deryy(vvme,nb_point,ly,beta);
//dryydy=dery(vvme,nb_point,ly,beta);
//dryydyy=dery(dryydy,nb_point,ly,beta);
vvvsdfm=dryydyy(1:ny)/((re**2)*(utau_m**4));
for i=1:ny vvvsdfp(i)=dryydyy(nb_point-i+1)/((re**2)*(utau_p**4)); end
scf();
semilogx(ym(2:ny),vvvsdfm(2:ny),'+-r',yp(2:ny),vvvsdfp(2:ny),'+--g',ydns(2:nb_dns),vvVsDfdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$ViscDiff - Ryy$');

// Corrélation vitesse - grad(p)
[y,vdpdyme] = read_my_stat_1d('vdpdym1d.dat',nb_point);
vvpstrm=-2*vdpdyme(1:ny)/(dt*re*(utau_m**4));
for i=1:ny vvpstrp(i)=-2*vdpdyme(nb_point-i+1)/(dt*re*(utau_p**4)); end
scf();
semilogx(ym(2:ny),vvpstrm(2:ny),'+-r',yp(2:ny),vvpstrp(2:ny),'+--g',ydns(2:nb_dns),vvPstrdns(2:nb_dns)+vvPDifdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$Pres-Str - Ryy$');

// Dissipation
[y,dvdx2me] = read_my_stat_1d('dvdx2m1d.dat',nb_point);
[y,dvdy2me] = read_my_stat_1d('dvdy2m1d.dat',nb_point);
[y,dvdz2me] = read_my_stat_1d('dvdz2m1d.dat',nb_point);
d2vme=dvdx2me+dvdy2me+dvdz2me;
vvdissm=-2*d2vme(1:ny)/((re**2)*(utau_m**4));
for i=1:ny vvdissp(i)=-2*d2vme(nb_point-i+1)/((re**2)*(utau_p**4)); end
scf();
semilogx(ym(2:ny),vvdissm(2:ny),'+-r',yp(2:ny),vvdissp(2:ny),'+--g',ydns(2:nb_dns),vvDissdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$Dissipation - Ryy$');

vvsumm=vvdissm+vvpstrm+vvtbdfm+vvvsdfm;
write_csv([ym,vvdissm,vvpstrm,vvtbdfm,vvvsdfm,vvsumm], "../csv/vv1.csv"," ");
vvsump=vvdissp+vvpstrp+vvtbdfp+vvvsdfp;
write_csv([yp,vvdissp,vvpstrp,vvtbdfp,vvvsdfp,vvsump], "../csv/vv2.csv"," ");

scf();
semilogx(ym(2:ny),vvdissm(2:ny),'+-r',ym(2:ny),vvpstrm(2:ny),'+-k',ym(2:ny),vvtbdfm(2:ny),'+r',ym(2:ny),vvvsdfm(2:ny),'+g',ym(2:ny),vvsumm(2:ny),'ok');
legend('$Diss$','$P-strain$','$Turb Diff$','$Visc Diff$','$Sum$',4);
xlabel('$y+$');
ylabel('$Budget <Ryy>, bottom$');
scf();
semilogx(ym(2:ny),vvdissm(2:ny),'+-r',ym(2:ny),vvpstrm(2:ny),'+-k',ym(2:ny),vvtbdfm(2:ny),'+r',ym(2:ny),vvvsdfm(2:ny),'+g',ydns(2:nb_dns),vvDissdns(2:nb_dns),'o-r',ydns(2:nb_dns),vvPDifdns(2:nb_dns)+vvPstrdns(2:nb_dns),'o-k',ydns(2:nb_dns),vvTbDfdns(2:nb_dns),'or',ydns(2:nb_dns),vvVsDfdns(2:nb_dns),'og')
xlabel('$y+$');
ylabel('$Budget <Ryy>, bottom$');

// Budget Rzz
// Turbulent diffusion
[y,vwwme] = read_my_stat_1d('vwwm1d.dat',nb_point);
dvwwdy=dery(vwwme,nb_point,ly,beta);
wwtbdfm=-dvwwdy(1:ny)/(re*(utau_m**4));
for i=1:ny wwtbdfp(i)=-dvwwdy(nb_point-i+1)/(re*(utau_p**4)); end
scf();
semilogx(ym(2:ny),wwtbdfm(2:ny),'+-r',yp(2:ny),wwtbdfp(2:ny),'+--g',ydns(2:nb_dns),wwTbDfdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$TurbDiff - Rzz$');

// Viscous diffusion
drzzdyy=deryy(wwme,nb_point,ly,beta);
wwvsdfm=drzzdyy(1:ny)/((re**2)*(utau_m**4));
for i=1:ny wwvsdfp(i)=drzzdyy(nb_point-i+1)/((re**2)*(utau_p**4)); end
scf();
semilogx(ym(2:ny),wwvsdfm(2:ny),'+-r',yp(2:ny),wwvsdfp(2:ny),'+--g',ydns(2:nb_dns),wwVsDfdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$ViscDiff - Rzz$');

// Corrélation vitesse - grad(p)
[y,wdpdzme] = read_my_stat_1d('wdpdzm1d.dat',nb_point);
wwpstrm=-2*wdpdzme(1:ny)/(dt*re*(utau_m**4));
for i=1:ny wwpstrp(i)=-2*wdpdzme(nb_point-i+1)/(dt*re*(utau_p**4)); end
scf();
semilogx(ym(2:ny),wwpstrm(2:ny),'+-r',yp(2:ny),wwpstrp(2:ny),'+--g',ydns(2:nb_dns),wwPstrdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$Pres-Str - Rzz$');

// Dissipation
[y,dwdx2me] = read_my_stat_1d('dwdx2m1d.dat',nb_point);
[y,dwdy2me] = read_my_stat_1d('dwdy2m1d.dat',nb_point);
[y,dwdz2me] = read_my_stat_1d('dwdz2m1d.dat',nb_point);
d2wme=dwdx2me+dwdy2me+dwdz2me;
wwdissm=-2*d2wme(1:ny)/((re**2)*(utau_m**4));
for i=1:ny wwdissp(i)=-2*d2wme(nb_point-i+1)/((re**2)*(utau_p**4)); end
scf();
semilogx(ym(2:ny),wwdissm(2:ny),'+-r',yp(2:ny),wwdissp(2:ny),'+--g',ydns(2:nb_dns),wwDissdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$Dissipation - Rzz$');

wwsumm=wwdissm+wwpstrm+wwtbdfm+wwvsdfm;
write_csv([ym,wwdissm,wwpstrm,wwtbdfm,wwvsdfm,wwsumm], "../csv/ww1.csv"," ");
wwsump=wwdissp+wwpstrp+wwtbdfp+wwvsdfp;
write_csv([yp,wwdissp,wwpstrp,wwtbdfp,wwvsdfp,wwsump], "../csv/ww2.csv"," ");
scf();
semilogx(ym(2:ny),wwdissm(2:ny),'+-r',ym(2:ny),wwpstrm(2:ny),'+-k',ym(2:ny),wwtbdfm(2:ny),'+r',ym(2:ny),wwvsdfm(2:ny),'+g',ym(2:ny),wwsumm(2:ny),'ok');
legend('$Diss$','$P-strain$','$Turb Diff$','$Visc Diff$','$Sum$',4);
xlabel('$y+$');
ylabel('$Budget <Rzz>, bottom$');
scf();
semilogx(ym(2:ny),wwdissm(2:ny),'+-r',ym(2:ny),wwpstrm(2:ny),'+-k',ym(2:ny),wwtbdfm(2:ny),'+r',ym(2:ny),wwvsdfm(2:ny),'+g',ydns(2:nb_dns),wwDissdns(2:nb_dns),'o-r',ydns(2:nb_dns),wwPstrdns(2:nb_dns),'o-k',ydns(2:nb_dns),wwTbDfdns(2:nb_dns),'or',ydns(2:nb_dns),wwVsDfdns(2:nb_dns),'og');
xlabel('$y+$');
ylabel('$Budget <Rzz>, bottom$');

// Budget Rxy
// Production
uvprodm=-vvm.*dudy(1:ny)/(re*(utau_m**2));
for i=1:ny uvprodp(i)=-vvp(i).*dudy(nb_point-i+1)/(re*(utau_p**2)); end
scf();
semilogx(ym(2:ny),uvprodm(2:ny),'+-r',yp(2:ny),-uvprodp(2:ny),'+--g',ydns(2:nb_dns),uvProddns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$Prod - Rxy$');

// Turbulent diffusion
[y,uvvme] = read_my_stat_1d('uvvm1d.dat',nb_point);
duvvdy=dery(uvvme,nb_point,ly,beta);
uvtbdfm=-duvvdy(1:ny)/(re*(utau_m**4));
for i=1:ny uvtbdfp(i)=-duvvdy(nb_point-i+1)/(re*(utau_p**4)); end
scf();
semilogx(ym(2:ny),uvtbdfm(2:ny),'+-r',yp(2:ny),-uvtbdfp(2:ny),'+--g',ydns(2:nb_dns),uvTbDfdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$TurbDiff - Rxy$');

// Viscous diffusion
drxydyy=deryy(uvme,nb_point,ly,beta);
//drxydy=dery(uvme,nb_point,ly,beta);
//drxydyy=dery(drxydy,nb_point,ly,beta);
uvvsdfm=drxydyy(1:ny)/((re**2)*(utau_m**4));
for i=1:ny uvvsdfp(i)=drxydyy(nb_point-i+1)/((re**2)*(utau_p**4)); end
scf();
semilogx(ym(2:ny),uvvsdfm(2:ny),'+-r',yp(2:ny),-uvvsdfp(2:ny),'+--g',ydns(2:nb_dns),uvVsDfdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$ViscDiff - Rxy$');

// Corrélation vitesse - grad(p)
[y,udpdyme] = read_my_stat_1d('udpdym1d.dat',nb_point);
[y,vdpdxme] = read_my_stat_1d('vdpdxm1d.dat',nb_point);
uvpstrm=-(udpdyme(1:ny)+vdpdxme(1:ny))/(dt*re*(utau_m**4));
for i=1:ny uvpstrp(i)=-(udpdyme(nb_point-i+1)+vdpdxme(nb_point-i+1))/(dt*re*(utau_p**4)); end
scf();
semilogx(ym(2:ny),uvpstrm(2:ny),'+-r',yp(2:ny),-uvpstrp(2:ny),'+--g',ydns(2:nb_dns),uvPstrdns(2:nb_dns)+uvPrDfdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$Pres-Str - Rxy$');

// Dissipation
[y,dudvdxme] = read_my_stat_1d('dudvdxm1d.dat',nb_point);
[y,dudvdyme] = read_my_stat_1d('dudvdym1d.dat',nb_point);
[y,dudvdzme] = read_my_stat_1d('dudvdzm1d.dat',nb_point);
d2uvme=dudvdxme+dudvdyme+dudvdzme;
uvdissm=-2*d2uvme(1:ny)/((re**2)*(utau_m**4));
for i=1:ny uvdissp(i)=-2*d2uvme(nb_point-i+1)/((re**2)*(utau_p**4)); end
scf();
semilogx(ym(2:ny),uvdissm(2:ny),'+-r',yp(2:ny),-uvdissp(2:ny),'+--g',ydns(2:nb_dns),uvDissdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$Dissipation - Rxy$');

uvsumm=uvdissm+uvprodm+uvpstrm+uvtbdfm+uvvsdfm;
write_csv([ym,uvdissm,uvprodm,uvpstrm,uvtbdfm,uvvsdfm,uvsumm], "../csv/uv1.csv"," ");
uvsump=-uvdissp-uvprodp-uvpstrp-uvtbdfp-uvvsdfp;
write_csv([yp,-uvdissp,-uvprodp,-uvpstrp,-uvtbdfp,-uvvsdfp,uvsump], "../csv/uv2.csv"," ");
scf();
semilogx(ym(2:ny),uvdissm(2:ny),'+-r',ym(2:ny),uvprodm(2:ny),'+-g',ym(2:ny),uvpstrm(2:ny),'+-k',ym(2:ny),uvtbdfm(2:ny),'+r',ym(2:ny),uvvsdfm(2:ny),'+g',ym(2:ny),uvsumm(2:ny),'ok');
legend('$Diss$','$Prod$','$P-strain$','$Turb Diff$','$Visc Diff$','$Sum$',4);
xlabel('$y+$');
ylabel('$Budget <Rxy>, bottom$');
scf();
semilogx(ym(2:ny),uvdissm(2:ny),'+-r',ym(2:ny),uvprodm(2:ny),'+-g',ym(2:ny),uvpstrm(2:ny),'+-k',ym(2:ny),uvtbdfm(2:ny),'+r',ym(2:ny),uvvsdfm(2:ny),'+g',ydns(2:nb_dns),uvDissdns(2:nb_dns),'o-r',ydns(2:nb_dns),uvProddns(2:nb_dns),'o-g',ydns(2:nb_dns),uvPstrdns(2:nb_dns)+uvPrDfdns(2:nb_dns),'o-k',ydns(2:nb_dns),uvTbDfdns(2:nb_dns),'or',ydns(2:nb_dns),uvVsDfdns(2:nb_dns),'og'); 
xlabel('$y+$');
ylabel('$Budget <Rxy>, bottom$');

// Budget U * T
// Production
[y,vphime] = read_my_stat_1d('vphim1d.dat',nb_point);
utprod=-vphime.*dudy-uvme.*dTdy;
utprodm=utprod(1:ny)/(utau_m*utau_m*Ttau_m*re*utau_m);
for i=1:ny utprodp(i)=utprod(nb_point-i+1)/(utau_p*utau_p*Ttau_p*re*utau_p); end
scf();
semilogx(ym(2:ny),utprodm(2:ny),'+-r',yp(2:ny),utprodp(2:ny),'+--g',ydns(2:nb_dns),utproddns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$Production - U_x T$');

// Dissipation
[y,dudtdxme] = read_my_stat_1d('dudphidxm1d.dat',nb_point);
[y,dudtdyme] = read_my_stat_1d('dudphidym1d.dat',nb_point);
[y,dudtdzme] = read_my_stat_1d('dudphidzm1d.dat',nb_point);
d2utme=dudtdxme+dudtdyme+dudtdzme;
utdissm=-(1+1/pr)*d2utme(1:ny)/((re**2)*(utau_m**3)*Ttau_m);//utdissm=-2*d2utme(1:ny)/((re**2)*(utau_m**3)*Ttau_m);
for i=1:ny utdissp(i)=-(1+1/pr)*d2utme(nb_point-i+1)/((re**2)*(utau_p**3)*Ttau_p); end //utdissp(i)=-2*d2utme(nb_point-i+1)/((re**2)*(utau_p**3)*Ttau_p); end
scf();
semilogx(ym(2:ny),utdissm(2:ny),'+-r',yp(2:ny),utdissp(2:ny),'+--g',ydns(2:nb_dns),utdissdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$Dissipation - U_x T$');

// Corrélation température / grad(p)
[y,tdpdxme] = read_my_stat_1d('phidpdxm1d.dat',nb_point);
utpstrm=-tdpdxme(1:ny)/(dt*re*(utau_m**3)*Ttau_m);
for i=1:ny utpstrp(i)=-tdpdxme(nb_point-i+1)/(dt*re*(utau_p**3)*Ttau_p); end
scf();
semilogx(ym(2:ny),utpstrm(2:ny),'+-r',yp(2:ny),utpstrp(2:ny),'+--g',ydns(2:nb_dns),utpstrdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',3);
xlabel('$y+$');
ylabel('$Pres-Str - U_x T$');

// Turbulent diffusion
[y,tuvme] = read_my_stat_1d('phiuvm1d.dat',nb_point);
dtuvdy=dery(tuvme,nb_point,ly,beta);
uttbdfm=-dtuvdy(1:ny)/(re*(utau_m**3)*Ttau_m);
for i=1:ny uttbdfp(i)=-dtuvdy(nb_point-i+1)/(re*(utau_p**3)*Ttau_p); end
scf();
semilogx(ym(2:ny),uttbdfm(2:ny),'+-r',yp(2:ny),uttbdfp(2:ny),'+--g',ydns(2:nb_dns),uttdifdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$TurbDiff - U_x T$');

// Viscous diffusion
[y,utme] = read_my_stat_1d('uphim1d.dat',nb_point);
[y,udtdxxme] = read_my_stat_1d('udtdxxm1d.dat',nb_point);
[y,udtdyyme] = read_my_stat_1d('udtdyym1d.dat',nb_point);
[y,udtdzzme] = read_my_stat_1d('udtdzzm1d.dat',nb_point);
udeltatme=udtdxxme+udtdyyme+udtdzzme;
dutdyy=deryy(utme,nb_point,ly,beta)+(1/pr-1)*udeltatme+(1/pr-1)*d2utme;
utvsdfm=dutdyy(1:ny)/((re**2)*(utau_m**3)*Ttau_m);
for i=1:ny utvsdfp(i)=dutdyy(nb_point-i+1)/((re**2)*(utau_p**3)*Ttau_p); end
scf();
semilogx(ym(2:ny),utvsdfm(2:ny),'+-r',yp(2:ny),utvsdfp(2:ny),'+--g',ydns(2:nb_dns),utmdifdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$ViscDiff - U_x T$');

utsumm=utdissm+utprodm+utpstrm+uttbdfm+utvsdfm;
write_csv([ym,utdissm,utprodm,utpstrm,uttbdfm,utvsdfm,utsumm], "../csv/ut1.csv"," ");
utsump=utdissp+utprodp+utpstrp+uttbdfp+utvsdfp;
write_csv([yp,utdissp,utprodp,utpstrp,uttbdfp,utvsdfp,utsump], "../csv/ut2.csv"," ");
scf();
semilogx(ym(2:ny),utdissm(2:ny),'+-r',ym(2:ny),utprodm(2:ny),'+-g',ym(2:ny),utpstrm(2:ny),'+-k',ym(2:ny),uttbdfm(2:ny),'+r',ym(2:ny),utvsdfm(2:ny),'+g',ym(2:ny),utsumm(2:ny),'ok');
legend('$Diss$','$Prod$','$P-strain$','$Turb Diff$','$Visc Diff$','$Sum$',4);
xlabel('$y+$');
ylabel('$Budget <U_xT>, bottom$');

// Budget V * T
// Production
vtprod=-vvme.*dTdy;
vtprodm=vtprod(1:ny)/(utau_m*utau_m*Ttau_m*re*utau_m);
for i=1:ny vtprodp(i)=vtprod(nb_point-i+1)/(utau_p*utau_p*Ttau_p*re*utau_p); end
scf();
semilogx(ym(2:ny),vtprodm(2:ny),'+-r',yp(2:ny),-vtprodp(2:ny),'+--g',ydns(2:nb_dns),vtproddns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$Production - U_y T$');

// Dissipation
[y,dvdtdxme] = read_my_stat_1d('dvdphidxm1d.dat',nb_point);
[y,dvdtdyme] = read_my_stat_1d('dvdphidym1d.dat',nb_point);
[y,dvdtdzme] = read_my_stat_1d('dvdphidzm1d.dat',nb_point);
d2vtme=dvdtdxme+dvdtdyme+dvdtdzme;
vtdissm=-(1+1/pr)*d2vtme(1:ny)/((re**2)*(utau_m**3)*Ttau_m);//-2*d2vtme(1:ny)/((re**2)*(utau_m**3)*Ttau_m);
for i=1:ny vtdissp(i)=-(1+1/pr)*d2vtme(nb_point-i+1)/((re**2)*(utau_p**3)*Ttau_p); end//-2*d2vtme(nb_point-i+1)/((re**2)*(utau_p**3)*Ttau_p); end
scf();
semilogx(ym(2:ny),vtdissm(2:ny),'+-r',yp(2:ny),-vtdissp(2:ny),'+--g',ydns(2:nb_dns),vtdissdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$Dissipation - U_y T$');

// Corrélation température / grad(p)
[y,tdpdyme] = read_my_stat_1d('phidpdym1d.dat',nb_point);
vtpstrm=-tdpdyme(1:ny)/(dt*re*(utau_m**3)*Ttau_m);
for i=1:ny vtpstrp(i)=-tdpdyme(nb_point-i+1)/(dt*re*(utau_p**3)*Ttau_p); end
scf();
//semilogx(ym(2:ny),vtpstrm(2:ny),'+-r',yp(2:ny),-vtpstrp(2:ny),'+--g',ydns(2:nb_dns),vtpstrdns(2:nb_dns)+vtpdifdns(2:nb_dns),'+k');
semilogx(ym(2:ny),vtpstrm(2:ny),'+-r',yp(2:ny),-vtpstrp(2:ny),'+--g',ydns(2:nb_dns),vtpstrdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',3);
xlabel('$y+$');
ylabel('$Pres-Str - U_y T$');

// Turbulent diffusion
[y,tvvme] = read_my_stat_1d('phivvm1d.dat',nb_point);
dtvvdy=dery(tvvme,nb_point,ly,beta);
vttbdfm=-dtvvdy(1:ny)/(re*(utau_m**3)*Ttau_m);
for i=1:ny vttbdfp(i)=-dtvvdy(nb_point-i+1)/(re*(utau_p**3)*Ttau_p); end
scf();
semilogx(ym(2:ny),vttbdfm(2:ny),'+-r',yp(2:ny),-vttbdfp(2:ny),'+--g',ydns(2:nb_dns),vttdifdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$TurbDiff - U_y T$');

// Viscous diffusion
[y,vtme] = read_my_stat_1d('vphim1d.dat',nb_point);
[y,vdtdxxme] = read_my_stat_1d('vdtdxxm1d.dat',nb_point);
[y,vdtdyyme] = read_my_stat_1d('vdtdyym1d.dat',nb_point);
[y,vdtdzzme] = read_my_stat_1d('vdtdzzm1d.dat',nb_point);
vdeltatme=vdtdxxme+vdtdyyme+vdtdzzme;
dvtdyy=deryy(vtme,nb_point,ly,beta)+(1/pr-1)*vdeltatme+(1/pr-1)*d2vtme;
//dvtdy=dery(vtme,nb_point,ly,beta);
//dvtdyy=dery(dvtdy,nb_point,ly,beta);
vtvsdfm=dvtdyy(1:ny)/((re**2)*(utau_m**3)*Ttau_m);
for i=1:ny vtvsdfp(i)=dvtdyy(nb_point-i+1)/((re**2)*(utau_p**3)*Ttau_p); end
scf();
semilogx(ym(2:ny),vtvsdfm(2:ny),'+-r',yp(2:ny),-vtvsdfp(2:ny),'+--g',ydns(2:nb_dns),vtmdifdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$ViscDiff - U_y T$');

vtsumm=vtdissm+vtprodm+vtpstrm+vttbdfm+vtvsdfm;
write_csv([ym,vtdissm,vtprodm,vtpstrm,vttbdfm,vtvsdfm,vtsumm], "../csv/vt1.csv"," ");
vtsump=-vtdissp-vtprodp-vtpstrp-vttbdfp-vtvsdfp;
write_csv([yp,-vtdissp,-vtprodp,-vtpstrp,-vttbdfp,-vtvsdfp,vtsump], "../csv/vt2.csv"," ");
scf();
semilogx(ym(2:ny),vtdissm(2:ny),'+-r',ym(2:ny),vtprodm(2:ny),'+-g',ym(2:ny),vtpstrm(2:ny),'+-k',ym(2:ny),vttbdfm(2:ny),'+r',ym(2:ny),vtvsdfm(2:ny),'+g',ym(2:ny),vtsumm(2:ny),'ok');
legend('$Diss$','$Prod$','$P-strain$','$Turb Diff$','$Visc Diff$','$Sum$',4);
xlabel('$y+$');
ylabel('$Budget <U_yT>, bottom$');

// Budget T * T
// Production
ttprod=-vphime.*dTdy;
ttprodm=ttprod(1:ny)/(utau_m*Ttau_m*Ttau_m*re*utau_m);
for i=1:ny ttprodp(i)=ttprod(nb_point-i+1)/(utau_p*Ttau_p*Ttau_p*re*utau_p); end
scf();
semilogx(ym(2:ny),ttprodm(2:ny),'+-r',yp(2:ny),ttprodp(2:ny),'+--g',ydns(2:nb_dns),ttproddns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$Production - T T$');

// Dissipation
[y,dtdx2me] = read_my_stat_1d('dphidx2m1d.dat',nb_point);
[y,dtdy2me] = read_my_stat_1d('dphidy2m1d.dat',nb_point);
[y,dtdz2me] = read_my_stat_1d('dphidz2m1d.dat',nb_point);
d2ttme=(dtdx2me+dtdy2me+dtdz2me)/(2*pr);
ttdissm=-2*d2ttme(1:ny)/((re**2)*(utau_m**2)*(Ttau_m**2));
ttdissym=-2*(dtdy2me(1:ny)/(2*pr))/((re**2)*(utau_m**2)*(Ttau_m**2));
for i=1:ny ttdissp(i)=-2*d2ttme(nb_point-i+1)/((re**2)*(utau_p**2)*(Ttau_p**2)); end
scf();
semilogx(ym(2:ny),ttdissm(2:ny),'+-r',yp(2:ny),ttdissp(2:ny),'+--g',ydns(2:nb_dns),ttdissdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$Dissipation - T T$');

// Turbulent diffusion
[y,ttvme] = read_my_stat_1d('vphi2m1d.dat',nb_point);
dttvdy=dery(ttvme,nb_point,ly,beta)/2;
tttbdfm=-dttvdy(1:ny)/(re*(utau_m**2)*(Ttau_m**2));
for i=1:ny tttbdfp(i)=-dttvdy(nb_point-i+1)/(re*(utau_p**2)*(Ttau_p**2)); end
scf();
semilogx(ym(2:ny),tttbdfm(2:ny),'+-r',yp(2:ny),tttbdfp(2:ny),'+--g',ydns(2:nb_dns),tttdifdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$TurbDiff - T T$');

// Viscous diffusion
[y,ttme] = read_my_stat_1d('phiphim1d.dat',nb_point);
dttdyy=deryy(ttme,nb_point,ly,beta)/(2*pr);
//dvtdy=dery(vtme,nb_point,ly,beta);
//dvtdyy=dery(dvtdy,nb_point,ly,beta);
ttvsdfm=dttdyy(1:ny)/((re**2)*(utau_m**2)*(Ttau_m**2));
for i=1:ny ttvsdfp(i)=dttdyy(nb_point-i+1)/((re**2)*(utau_p**2)*(Ttau_p**2)); end
scf();
semilogx(ym(2:ny),ttvsdfm(2:ny),'+-r',yp(2:ny),ttvsdfp(2:ny),'+--g',ydns(2:nb_dns),ttmdifdns(2:nb_dns),'+k');
legend('$m$','$p$','$DNS$',4);
xlabel('$y+$');
ylabel('$ViscDiff - T T$');

ttsumm=ttdissm+ttprodm+tttbdfm+ttvsdfm;
write_csv([ym,ttdissm,ttprodm,tttbdfm,ttvsdfm,ttsumm,ttdissym], "../csv/tt1.csv"," ");
ttsump=ttdissp+ttprodp+tttbdfp+ttvsdfp;
write_csv([yp,ttdissp,ttprodp,tttbdfp,ttvsdfp,ttsump], "../csv/tt2.csv"," ");
scf();
semilogx(ym(2:ny),ttdissm(2:ny),'+-r',ym(2:ny),ttprodm(2:ny),'+-g',ym(2:ny),tttbdfm(2:ny),'+r',ym(2:ny),ttvsdfm(2:ny),'+g',ym(2:ny),ttsumm(2:ny),'ok');
legend('$Diss$','$Prod$','$Turb Diff$','$Visc Diff$','$Sum$',4);
xlabel('$y+$');
ylabel('$Budget <TT>, bottom$');
scf();
semilogx(ym(2:ny),ttdissm(2:ny),'+-r',ym(2:ny),ttprodm(2:ny),'+-g',ym(2:ny),tttbdfm(2:ny),'+r',ym(2:ny),ttvsdfm(2:ny),'+g',ydns(2:nb_dns),ttdissdns(2:nb_dns),'o-r',ydns(2:nb_dns),ttproddns(2:nb_dns),'o-g',ydns(2:nb_dns),ttmdifdns(2:nb_dns),'or',ydns(2:nb_dns),tttdifdns(2:nb_dns),'og');
xlabel('$y+$');
ylabel('$Budget <TT>, bottom$');
