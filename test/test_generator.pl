#!/usr/bin/perl

do{
	($m,$n)=(int(rand()*14)+2,int(rand()*14)+2);
}while( $m < $n );

print "$m $n\n";

for( $i = 0; $i < $m; $i++ ){
	for( $j = 0; $j < $n + 1; $j++ ){
		$x = (rand()*100);
		$x = (rand() > 0.5 ? $x : 0 );
		print sprintf("%.4f ",$x);
	}
print "\n";
}