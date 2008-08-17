#!/usr/bin/perl

$l=<STDIN>;

while(<STDIN>){
	$c.=$_;
}

while($c =~ /((\d|\.)+)/g){
	push @x, $1;
}

#print join(" ",@x)."\n";
($m,undef,$n) = ($l =~ /(\d+)(.*?)(\d+)/);

#print "m=$m, n=$n\n";

open(MAC,">/tmp/maxima.mac");

print MAC 'ttyoff : true;'."\n";
print MAC 'load(lsquares);'."\n";
print MAC 'numer : true;'."\n";
print MAC 'M : matrix (';

print MAC '[';
for($i=0;$i < $#x+1;$i++){
	print MAC $x[$i];
	if( ($i+1) % 3 == 0 ){
		print MAC ']';
		if( $i+2 < $#x+1 ){
			print MAC ', [';
		}
	}else{
		print MAC ', ';
	}
}
print MAC ' );'."\n";
print MAC 'r : lsquares_estimates ( M, [';
for($i=1;$i < $n+1;$i++){
	print MAC 'a'.$i;
	if( $i+1 < $n+1 ){
		print MAC ', ';
	}
}
print MAC ', b], ';
for($i=1;$i < $n+1;$i++){
	print MAC 'a'.$i.'*x'.$i;
	if( $i+1 < $n+1 ){
		print MAC '+';
	}
}
print MAC '=b, [';
for($i=1;$i < $n+1;$i++){
	print MAC 'x'.$i;
	if( $i+1 < $n+1 ){
		print MAC ', ';
	}
}
print MAC '] );'."\n";
print MAC 'r : r[1];'."\n";
print MAC 'for i:1 step 1 thru length(r) do r[i] : part(r[i],2);'."\n";
print MAC 'printf( true, "~{~,4f ~}~%", r );'."\n";

close(MAC);

system('maxima -b /tmp/maxima.mac | tail -n 1');
