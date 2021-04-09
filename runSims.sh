

g++ -o outf eff_pop_2d_rec_lastrow.cpp -lgsl

for b in  1.5 1.75 2 2.25 2.5 2.75 3 3.25 3.5 4 4.5 5 5.5 6 6.5 7 
do
	for s in 1 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 
	do
		#for (( i = 0; i < 50 ; i++ )) ### Outer for loop ###
		#do

		./outf -B $b -S $s #-I $[i*u]
	  	cnt=$[cnt+1]
	  	echo $cnt "out of" $tot "\n"

		#done
	done
done