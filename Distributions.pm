package Statistics::Distributions;

use strict;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);
use vars qw($n $p $x $m $y $i $c $d $e $q $delta $round $z);
use constant PI => 3.1415926536;
require Exporter;

@ISA = qw(Exporter AutoLoader);
# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.
@EXPORT_OK = qw(chisqrdistr tdistr fdistr udistr);
$VERSION = '0.02';


# Preloaded methods go here.
   
sub chisqrdistr { # Percentage points  X’(x’,þ)
    $n=shift;

    
    $p=shift;
    if (($n<=0) || ((abs($n)-(abs(int($n))))!=0)) {
	die "Invalid n: $n\n"; # degree of freedom
    }
    if (($p<=0) || ($p>1)) {
	die "Invalid p: $p\n"; 
    }
    &_subchisqr;
    if ($x) {
	my $precision=abs int(log10(abs $x)-5);
	my $x=sprintf"%.".$precision."f",$x;
	return $x;
    }
}

sub udistr { # Percentage points   N(0,1’)
    $p=shift;
    if (($p>1) || ($p<=0)) {
	die "Invalid p: $p\n";
    }
    &_subu;
    if ($x) {
	$x=sprintf"%.".abs(int(log10(abs $x)-5))."f",$x;
    }
    return $x;
}

sub tdistr { # Percentage points   t(x,þ)
    $n=shift;
    $p=shift;
    if (($n<=0) || ((abs($n)-(abs(int($n))))!=0)) {
	die "Invalid n: $n\n";
    }
    if (($p<=0) || ($p>=1)) {
	die "Invalid p: $p\n";
    }
    &_subt;
    if ($x) {
	$x=sprintf "%.".abs(int(log10(abs $x)-5))."f",$x;
    }
    return $x;
}

sub fdistr { # Percentage points  F(x,þü,þý)
    $n=shift;
    $m=shift;
    $p=shift;
    if (($n<=0) || ((abs($n)-(abs(int($n))))!=0)) {
	die "Invalid n: $n\n"; # first degree of freedom
    }
    if (($m<=0) || ((abs($m)-(abs(int($m))))!=0)) {
	die "Invalid m: $m\n"; # second degree of freedom
    }
    if (($p<=0) || ($p>1)) {
	die "Invalid p: $p\n";
    }
    &_subf;
    if ($x) {
	$x=sprintf"%.".abs(int(log10(abs $x)-5))."f",$x;
    }  
    return $x;
}

sub _subu {
    $y=-log(4*$p*(1-$p));
    $x=.5824238515E-5+$y*(-.104527497E-5+$y*(.8360937017E-7+$y*(-.3231081277E-8+$y*(.3657763036E-10+$y*.6936233982E-12))));
    $x=sqrt($y*(1.570796288+$y*(.03706987906+$y*(-.8364353589E-3+$y*(-.2250947176E-3+$y*(.6841218299E-5+$y*$x))))));
    $x=-$x if ($p>.5);
    $y=$x;
}

sub _subchisqr2 {
    $y=abs($x);
    $p=0;
    if ($y>100) {
	$p=1-$p if ($x<0);
	$y=$p;
    }
    elsif ($y<1.9) {
	$p=(1+$y*(.049867347+$y*(.0211410061+$y*(.0032776263+$y*(.0000380036+$y*(.0000488906+$y*.000005383))))))**-16/2;
	$p=1-$p if ($x<0);
	$y=$p;
    }
    else {
	for ($i=18;$i>=1;$i--) {
	    $p=$i/($y+$p);
	}
	$p=exp(-.5*$y*$y)/sqrt(2*PI)/($y+$p);
	$p=1-$p if ($x<0);
	$y=$p;
    }
}

sub log10 {
my $n = shift;
return log($n)/log(10);
}
    
sub _subt {
    if (($p>=1) || ($p<=0)) {
	die "Invalid p: $p\n";
    }
    &_subu;
    $x=$y;
    $y=$x**2;
    $a=($y+1)/4;
    $b=((5*$y+16)*$y+3)/96;
    $c=(((3*$y+19)*$y+17)*$y-15)/384;
    $d=((((79*$y+776)*$y+1482)*$y-1920)*$y-945)/92160;
    $e=(((((27*$y+339)*$y+930)*$y-1782)*$y-765)*$y+17955)/368640;
    $x=$x*(1+($a+($b+($c+($d+$e/$n)/$n)/$n)/$n)/$n);
    if ($n>(log10($p)**2+3)) {
	$y=$x;
    }
    else {
	do {	
	    $q=$p;
	    &_subt2;
	    $p=$y;
	    $b=$n+1;
	    $a=exp(($b*log($b/($n+$x*$x))+log($n/$b/2/PI)-1+(1/$b-1/$n)/6)/2);
	    $y=$x;
	    $x=$x+($p-$q)/$a;
	    $p=$q;
	    $delta=$x-$y;
	    $round=sprintf("%.".abs(int(log10(abs $x)-4))."f",$delta);
	} while (($x) && ($round!=0));
	$y=$x;
    }
}

sub _subt2 {
    $y=atan2($x/sqrt($n),1);
    $z=cos $y**2;
    if ($n%2==0) {
	$a=sin($y)/2;
	$b=.5;
    }
    else {
	$b=.5+$y/PI;
	if ($n==1) {
	    $a=0;
	}
	else {   
	    $a=sin($y)*cos($y)/PI;
	}
    }
    $y=1;
    for ($i=$n-2;$i>=2;$i-=2) {
	$y=1+($i-1)/$i*$z*$y;
    } 
    $p=1-($b+$a*$y);
    if ($p<0) {
	$p=0;
    }
    $y=$p;
}

sub _subf {
    if (($p>=1) || ($p<=0)) {
	die "Invalid p: $p\n";
    }
    if ($p==1) {
	$x=0;
	$y=$x;
    }
    elsif ($m==1) {
	$m=$p;
	$p=.5-$p/2;
	&_subt;
	$p=$m;
	$m=1;
	$x=1/$y**2;
	$y=$x;
    }
    elsif ($n==1) {
	$n=$m;
	$p=$p/2;
	&_subt;
	$n=1;
	$p=$p*2;
	$x=$y**2;
	$y=$x;
    }
    elsif ($m==2) {
	$p=1-$p;
	$m=$n;
	$n=2;
	&_subchisqr;
	$x=$y;
	$a=$n-2;
	$x=$x/$n*(1+(($x-$a)/2+(((4*$x-11*$a)*$x+$a*(7*$n-10))/24+(((2*$x-10*$a)*$x+$a*(17*$n-26))*$x-$a*$a*(9*$n-6))/48/$m)/$m)/$m);
	$p=1-$p;
	$n=$m;
	$m=2;
	$x=1/$x;
	$y=$x;
    }
    elsif ($n>$m) {
	$p=1-$p;
	$d=$n;
	$n=$m;
	$m=$d;
	&_subf2;
	$x=1/$x;
	$d=$m;
	$m=$n;
	$n=$d;
	$p=1-$p;
	$y=$x;
    }
    else {
	&_subf2;
	$y=$x;
    }
}

sub _subf2 {
    &_subchisqr;
    $x=$y;
    $a=$n-2;
    $x=$x/$n*(1+(($x-$a)/2+(((4*$x-11*$a)*$x+$a*(7*$n-10))/24+(((2*$x-10*$a)*$x+$a*(17*$n-26))*$x-$a*$a*(9*$n-6))/48/$m)/$m)/$m);
    do {
	$d=$x;
	$c=$p;
	&S6240;
	$p=$c;
	$z=$n+$m;
	$z=exp(($z*log($z/($n*$x+$m))+($n-2)*log $x+log($n*$m/$z)-log(4*PI)-(1/$n+1/$m-1/$z)/6)/2);
	$x=$x+($y-$p)/$z;
    } while (abs($d-$x)>3e-4);
}

sub _subchisqr {
    if (($p>1) || ($p<=0)) {
	die "Invalid p: $p\n";
    }
    elsif ($p==1){
	$x=0;
	$y=$x;
    }
    elsif ($n==1) {
	$q=$p;
	$p=$q/2;
	&_subu;
	$x=$y*$y;
	$p=$q;
	$y=$x;
    }
    elsif ($n==2) {
	$x=-2*log($p);
	$y=$x;
    }
    else {
	&_subu;
	$x=$y;
	$y=$x*$x;
	$x=$n+sqrt(2*$n)*$x+2/3*($y-1)+$x*($y-7)/9/sqrt(2*$n)-2/405/$n*($y*(3*$y+7)-16);
	$x=0 if ($x<0);
	if ($n>100) {
	    $y=$x;
	}
	else {
	    do {
		$b=$x;
		$q=$p;
		if ($x<0) {
		    $p=1;
		    $y=$p;
		}
		elsif ($n>100) {
		    $z=$x;
		    $x=(($x/$n)**(1/3)-(1-2/9/$n))/sqrt(2/9/$n);
		    &_subchisqr2;
		    $p=$y;
		    $x=$z;
		    $y=$p;
		}
		elsif ($x>400) {
		    $p=0;
		    $y=$p;
		}
		else {
		    $a=exp(-$x/2);
		    $p=$a;
		    $y=2;
		    if (($n % 2)!=0) {
			$z=$x;
			$x=sqrt($x);
			&_subchisqr2;
			$p=2*$y;
			$a=sqrt(2/PI)*$a/$x;
			$x=$z;
			$y=1;
		    }
		    for ($i=$y;$i<=($n-2);$i+=2) {
			$a=$a*$x/$i;
			$p+=$a;
		    }
		    $y=$p;
		}
		$p=$y;
		$z=exp((($n-1)*log($x/$n)-log(4*PI*$x)+$n-$x-1/$n/6)/2);
		$x=$x+($p-$q)/$z;
		$p=$q;
		$x=sprintf("%.5f",$x);
	    } while (($n<31) && (abs($b-$x)>1e-4));
	    $y=$x;
	}
    }
}

# Autoload methods go after =cut, and are processed by the autosplit program.

1;
__END__
# Below is the stub of documentation for your module. You better edit it!

=head1 NAME

Statistics::Distributions - Perl module for calculating critical values of common statistical distributions

=head1 SYNOPSIS

  use Statistics::Distributions;

  $chis=chisqrdistr (2,.05);
  print "Chi-squared-crit (2 degrees of freedom, 95th percentile = 0.05 level) = $chis\n";
  $u=udistr (.05);
  print "u-crit (95th percentile = 0.05 level) = $u\n";
  $t=tdistr (1,.005);
  print "t-crit (1 degree of freedom, 99.5th percentile = 0.005 level) =$t\n";
  $f=fdistr (1,3,.01);
  print "F-crit (1 degree of freedom in numerator, 3 degrees of freedom in denominator, 99th percentile = 0.01 level) = $f\n";


=head1 DESCRIPTION

This Perl module calulates percentage points (5 significant digits) of the u (standard normal) distribution, the student's t distribution, the chi-square distribution and the F distribution.
These critical values are needed to perform statistical tests, like the u test, the t test, the F test and the chi-squared test, and to calculate confidence intervals.

=head1 AUTHOR

Michael Kospach, mike.perl@gmx.at

=head1 SEE ALSO

Statistics::ChiSquare, Statistics::Table::t, Statistics::Table::F
    perl(1).

=cut
