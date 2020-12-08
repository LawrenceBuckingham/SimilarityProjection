for f in *.hpp *.h; do
	if [ "$f" != "QutBio.h" ]; then
		echo '#include "'$f'"'
	fi 
done >QutBio.h

echo done.