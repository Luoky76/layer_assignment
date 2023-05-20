for index in 2 3 6 7 9 11 12 14 16 19
do 
    for i in `seq 1 5`
    do
        echo $index | ./main > solution/threads_adaption/data${index}_out${i}.txt
    done
done