if (( $# != 2 ))
then
	echo 'Expected 2 command line arguments:
	swissprot-fastaFile  The name of a file containing the test set.
	homolog-file         The name of the homolog file.
	'
	exit 1
fi

fastaFile=$1
homologFile=$2

awk '
BEGIN {
	fileNum = 0;
}

FNR == 1 {
	fileNum ++;
}

fileNum == 1 {
	split( $0, parts ,"|" );
	queryIds[parts[3]] = 1;
}

fileNum == 2 {
	if ( length($1) > 0 && $1 in queryIds ) {
		printf( "%s", $1 );
		
		for ( i = 3; i < NF; i++ ) {
			if ( ! ($i in queryIds ) ) {
				printf(" %s", $i);
			}
		}
		
		printf("\n");
	}
}

' $fastaFile $homologFile 
