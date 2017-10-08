for x in 100
do
#for y in 1.6 3.2 6.3 12.6 17.8 35.5 70.8 
for y in 1.6
do
for z in 0.25 # 5.0
do
qsub -v N=$x,T=$y,s=$z -N Cat_N${x}_T${y}_s${z} submit_g.sh
done
done
done
