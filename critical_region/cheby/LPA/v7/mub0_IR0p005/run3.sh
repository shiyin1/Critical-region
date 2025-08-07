#!/bin/bash
rm exam/exe
rm exam/*.out
rm exam/*.o
rm adata/*.dat
rm bdata/*.dat
rm dVd1rdata/*.dat
cd exam
make
rm *.o
cd ..
rm -rf Tem*
for a in {1..10} 
do
mkdir Tem$a
cp -r exam/exe Tem$a
mkdir Tem$a/buffer
echo $a >Tem$a/m1.dat
echo -e "#!/bin/bash\ncd Tem$a\n./exe" >Tem$a/run.sh
chmod 755 run.sh
done
sh Tem1/run.sh &
sh Tem2/run.sh &
sh Tem3/run.sh &
sh Tem4/run.sh &
sh Tem5/run.sh &
sh Tem6/run.sh &
sh Tem7/run.sh &
sh Tem8/run.sh &
sh Tem9/run.sh &
sh Tem10/run.sh &
wait
for a in {1..10} 
do
mv Tem$a/buffer/a.dat adata/a$a.dat
mv Tem$a/buffer/b.dat bdata/b$a.dat
mv Tem$a/buffer/dVd1rho.dat dVd1rdata/dvdr$a.dat
done
rm -rf Tem*