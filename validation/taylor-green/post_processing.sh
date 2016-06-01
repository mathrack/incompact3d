for i in ab2*; do cd $i; tail -n20 */log | grep Liw | cut -c21- > tmp; paste ../ncell.dat tmp > liw; rm tmp; cd ..; done
for i in ab2*; do cd $i; tail -n20 */log | grep L2w | cut -c21- > tmp; paste ../ncell.dat tmp > l2w; rm tmp; cd ..; done
for i in ab2*; do cd $i; tail -n20 */log | grep T | grep inf | cut -c21- > tmp; paste ../ncell.dat tmp > tinf; rm tmp; cd ..; done
for i in ab2*; do cd $i; tail -n20 */log | grep T2 | cut -c21- > tmp; paste ../ncell.dat tmp > t2; rm tmp; cd ..; done

for i in ab2*; do cd $i; tail -n20 */log_imp | grep Liw | cut -c21- > tmp; paste ../ncell.dat tmp > liw_imp; rm tmp; cd ..; done
for i in ab2*; do cd $i; tail -n20 */log_imp | grep L2w | cut -c21- > tmp; paste ../ncell.dat tmp > l2w_imp; rm tmp; cd ..; done
for i in ab2*; do cd $i; tail -n20 */log_imp | grep T | grep inf | cut -c21- > tmp; paste ../ncell.dat tmp > tinf_imp; rm tmp; cd ..; done
for i in ab2*; do cd $i; tail -n20 */log_imp | grep T2 | cut -c21- > tmp; paste ../ncell.dat tmp > t2_imp; rm tmp; cd ..; done
