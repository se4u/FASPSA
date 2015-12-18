for i in 15 30 60 120 240 480 960 ;
do
    tail -n 5 ../res/timing_${i}.txt;
done
