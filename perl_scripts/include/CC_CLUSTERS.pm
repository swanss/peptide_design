package CC_CLUSTERS;

# The purpose of this module is to implement functions related to operating with point, pair (and maybe eneventually higher-order) clusters in 
# parallel dimeric coiled coils (again, maybe eventually this can be extended for any oligomerization state and orientation or even for any generic 
# symmetrical or non-symmetrical system.
# In this module it is assumed that coiled-coil positions are named a-f using a digit to designate heptad number and a 'p' to designate the 
# opposing strand. For example: a d e g a1 d1 e1 g1 ap dp ep gp ap1 dp1 ep1 gp1. A heptad is always assumed to start with an 'a' position, even if 
# the sequence does not (this is just a convenient convention). So for example, if a sequence starts with f, the pattern would be something like 
# this: f g a1 b1 d1 e1 f1... (or g a1 d1 e1 g1 if we only have adeg). So 'g a d e g1 a1 d1 e1' would be an invalid pattern.

use GENERAL;

=head2

 Title   :  reducePoint
 Usage   :  CC_CLUSTERS::reducePoint($pi);
 Function:  Reduces a point cluster to its most basic representation
 Returns :
 Args    :  1. heptad assignment of the position (e.g. g1)

=cut

sub reducePoint {
  my $pi = shift;
  GENERAL::requireArgs($pi);

  $pi =~ /^(.)(p?)(\d*)$/;
  my $hi = $1;

  return $hi;
}


=head2

 Title   :  reducePair
 Usage   :  CC_CLUSTERS::reducePair($pi, $pj);
 Function:  Reduces a pair cluster to its most basic representation
 Returns :
 Args    :  1. heptad assignment of the first position (e.g. g1)
            2. heptad assignment of the second position (e.g. ap1)

=cut

sub reducePair {
  my $pi = shift;
  my $pj = shift;
  GENERAL::requireArgs($pi, $pj);
  my $L = shift;
  $L = 2 if (!defined($L));

  my $rpi = "";
  my $rpj = "";

  $pi =~ /^(.)(p?)(\d*)$/;
  my $hi = $1;
  if (defined($2) && ($2 ne "")) { $pi = 1; }
  else { $pi = 0; }
  if (defined($3) && ($3 ne "")) { $ci = $3; }
  else { $ci = 0; }

  $pj =~ /^(.)(p?)(\d*)$/;
  my $hj = $1;
  if (defined($2) && ($2 ne "")) { $pj = 1; }
  else { $pj = 0; }
  if (defined($3) && ($3 ne "")) { $cj = $3; }
  else { $cj = 0; }

  # things too far away do not interact
  if ((abs($ci - $cj) > $L) || (($cj == $ci + $L) && ($hj gt $hi)) || (($ci == $cj + $L) && ($hi gt $hj))) {
    return ($rpi, $rpj, 0, 0);
  }

  # if both primed, go to the other chain
  if ($pi && $pj) { $pi = 0; $pj = 0; }

  # go to first heptad pair
  my $d = GENERAL::min($ci, $cj);
  $ci -= $d; $cj -= $d;

  # determine order
  if (($ci < $cj) || (($ci == $cj) && ($hi le $hj))) {
    # i before j
    $rpi .= $hi; $rpj .= $hj;
    $rpj .= "p" if ($pi || $pj);
    $rpi .= $ci; $rpj .= $cj;
    return ($rpi, $rpj, 0, (($ci == $cj) && ($hi eq $hj)));
  } elsif (($ci > $cj) || (($ci == $cj) && ($hi gt $hj))) {
    # j before i
    $rpi .= $hj; $rpj .= $hi;
    $rpj .= "p" if ($pi || $pj);
    $rpi .= $cj; $rpj .= $ci;
    return ($rpi, $rpj, 1, (($ci == $cj) && ($hi eq $hj)));
  }
}

=head2

 Title   :  getPointCF
 Usage   :  CC_CLUSTERS::getPointCF($S, $cl, $aa);
 Function:  Returns the value of a point cluster function
 Returns :
 Args    :  1. Reference to a cluster function hash, where all CF's are hashed by name
            2. Cluster definition are returned by CC_CLUSTERS::reducePoint
            3. Name of amino acid

=cut

sub getPointCF {
  my $S = shift;
  my $cl = shift;
  my $aa = shift;
  
  my $cfn = "$cl $aa";
  if(defined($S->{$cfn})) { return $S->{$cfn}; }
  else { return 0; }
}

=head2

 Title   :  getPairCF
 Usage   :  CC_CLUSTERS::getPairCF($S, \@cl, $aai, $aaj);
 Function:  Returns the value of a pair cluster function
 Returns :
 Args    :  1. Reference to a cluster function hash, where all CF's are hashed by name
            2. Cluster definition are returned by CC_CLUSTERS::reducePair
            3. Name of first amino acid
            4. Name of second amino acid (in the order of the positions passed to CC_CLUSTERS::reducePair)

=cut

sub getPairCF {
  my $S = shift;
  my $cl = shift;
  my $aai = shift;
  my $aaj = shift;
  
  return 0 if ($cl->[0] eq "");
  # if cluster defined in the order we are visiting the pair
  my $cfnr = "$cl->[0]-$cl->[1] $aai $aaj";
  # cluster defined opposite to the order we are visiting the pair
  my $cfnf = "$cl->[0]-$cl->[1] $aaj $aai";
  # if homo-typic interaction, try both orientations (they are the same, but CF is defined only one way)
  if ($cl->[3]) {
    if (defined($S->{$cfnr})) { return $S->{$cfnr}; }
    elsif (defined($S->{$cfnf})) { return $S->{$cfnf}; }
    else { return 0; }
  # if hetero-typic interaction, check order in which the cluster is defined
  } elsif ($cl->[2]) {
    if (defined($S->{$cfnf})) { return $S->{$cfnf} }
    else { return 0; }
  } else {
    if (defined($S->{$cfnr})) { return $S->{$cfnr} }
    else { return 0; }
  }
}


1;
