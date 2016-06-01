for i in ab2*; do cd $i; for j in 0* 128; do cd $j; pwd; make; ./incompact3d > NAME_OF_YOUR_LOG_FILE_HERE; make clean; rm *dat; cd ..; done; cd ..; done
