for value in {1..10}
do

echo "RUNNING 1 THREADS"
./main -threads 1

done

for value in {1..10}
do

echo "RUNNING 2 THREADS"
./main -threads 2

done

for value in {1..10}
do

echo "RUNNING 4 THREADS"
./main -threads 4

done

for value in {1..10}
do

echo "RUNNING 8 THREADS"
./main -threads 8

done

for value in {1..10}
do

echo "RUNNING 16 THREADS"
./main -threads 16

done