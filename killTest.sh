cd dist/build/tester-DeUni

for i in $(cat ./pidList)
do	
	echo "Killing $i"
	kill $i
done
	rm pidList
