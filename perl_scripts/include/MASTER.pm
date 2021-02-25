# By Gevorg Grigoryan

package MASTER;
$DEF_VER_CONTROL->{date}->{MASTER} = "05/22/2015";

use strict;
use GENERAL;
use PDB;

sub rmsdCut {
  my $pdb = shift;
  my $max = shift;
  my $L0 = shift;
  GENERAL::requireArgs($pdb, $max, $L0);

  # figure out how many segments and how long each
  my @N;
  my @res = $pdb->ConRes();
  my $n = 0;
  for (my $i = 1; $i < scalar(@res); $i++) {
    if (($i == scalar(@res) - 1) || (!PDB::isBondedTo($res[$i-1], $res[$i]))) {
      push(@N, $n+1); $n = 0;
    } else {
      $n++;
    }
  }

  return rmsdCutL(\@N, $max, $L0);
}

sub rmsdCutL {
  my $N = shift;
  my $max = shift;
  my $L0 = shift;
  GENERAL::requireArgs($N, $max, $L0);

  return $max/sqrt(GENERAL::sum($N)/df($N, $L0));
}

sub df {
  my $N = shift;
  my $L0 = shift;
  GENERAL::requireArgs($N, $L0);

  my $a = exp(-1/$L0);
  # disjoint segments are counted as independent, so their correlation with respect to each other is zero
  my $c = 0;
  for (my $i = 0; $i < scalar(@$N); $i++) {
    my $n = $N->[$i];
    $c += ($a/(1-$a))*($n-1) - ((($a/(1-$a)))**2) * (1 - $a**($n-1));
  }
  my $n = GENERAL::sum($N);
  return $n*(1 - (2/($n*($n-1)))*$c);
}

1;
