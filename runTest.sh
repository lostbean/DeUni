cd dist/build/tester-DeUni

for i in {1..9}
do	
	rm -Rf test$i
	mkdir ./test$i
	cd test$i
	../tester-DeUni 2> /dev/null > testLog.txt &
	echo $! >> ../pidList
	cd ..
done
