

g++ -o outf eff_pop_2d_rec_lastrow.cpp -lgsl

for b in  1.5 1.75 2 2.25 2.5 2.75 3 3.25 3.5 4
do
	for s in 1 2 4 6 8 16 32  
	do
		for (( i = 0; i < 10 ; i++ )) ### Outer for loop ###
		do
			./outf -B $b -S $s -I $[i*u]
		  	cnt=$[cnt+1]
		  	echo $cnt "out of" $tot "\n"

		done
	done
done