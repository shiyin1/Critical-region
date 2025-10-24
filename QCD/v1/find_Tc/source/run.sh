#!/bin/bash
rm exam/exe
rm exam/*.out
rm exam/*.o
cd exam
make
rm *.o
rm *.mod
cd ..
rm -rf data/*
rm -rf Tem*
for a in {1..10} 
do
mkdir Tem$a
cp -r exam/eqcd Tem$a
mkdir Tem$a/BUFFER
echo $a >Tem$a/m1.dat
echo -e "#!/bin/bash\ncd Tem$a\n./eqcd" >Tem$a/run.sh
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