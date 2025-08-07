#!/bin/bash
rm exam/exe
rm exam/*.out
rm exam/*.o
cd exam
make
rm *.o
cd ..
rm -rf data/*
rm -rf Tem*
for a in {1..58} 
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
sh Tem11/run.sh &
sh Tem12/run.sh &
sh Tem13/run.sh &
sh Tem14/run.sh &
sh Tem15/run.sh &
sh Tem16/run.sh &
sh Tem17/run.sh &
sh Tem18/run.sh &
sh Tem19/run.sh &
sh Tem20/run.sh &
sh Tem21/run.sh &
sh Tem22/run.sh &
sh Tem23/run.sh &
sh Tem24/run.sh &
sh Tem25/run.sh &
sh Tem26/run.sh &
sh Tem27/run.sh &
sh Tem28/run.sh &
sh Tem29/run.sh &
sh Tem30/run.sh &
sh Tem31/run.sh &
sh Tem32/run.sh &
sh Tem33/run.sh &
sh Tem34/run.sh &
sh Tem35/run.sh &
sh Tem36/run.sh &
sh Tem37/run.sh &
sh Tem38/run.sh &
sh Tem39/run.sh &
sh Tem40/run.sh &
sh Tem41/run.sh &
sh Tem42/run.sh &
sh Tem43/run.sh &
sh Tem44/run.sh &
sh Tem45/run.sh &
sh Tem46/run.sh &
sh Tem47/run.sh &
sh Tem48/run.sh &
sh Tem49/run.sh &
sh Tem50/run.sh &
sh Tem51/run.sh &
sh Tem52/run.sh &
sh Tem53/run.sh &
sh Tem54/run.sh &
sh Tem55/run.sh &
sh Tem56/run.sh &
sh Tem57/run.sh &
sh Tem58/run.sh &
wait
for a in {1..58} 
do
mv Tem$a/buffer/mf.dat mfdata/mf$a.dat
mv Tem$a/buffer/Zs.dat Zphidata/Z$a.dat
done