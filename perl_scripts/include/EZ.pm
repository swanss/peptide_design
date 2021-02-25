# module for calculating EZ energies
# author: Gevorg Grigoryan
# start date: March 2, 2010
# info: initial data for EZ potential taken from the DeGrado web server scripts

package EZ;

use GENERAL;
use SEQUENCE;
use PDB;

=head2

 Title   :  new
 Usage   :  my $ez = EZ::new("version")
 Function:  EZ object constructor.
 Returns :  new EZ object.
 Args    :  1. version of potential to use

=cut

sub new {
  my $self = {};
  my $ver = shift;
  $ver = "default" if (!defined($ver));

  bless $self;
  $self->{version} = $ver;

  if ($ver eq "default") {
    $self->{mwidth} = 28.5;

    my %potential;
    $potential{"A"} = "sig";
    $potential{"D"} = "sig";
    $potential{"E"} = "sig";
    $potential{"F"} = "sig";
    $potential{"G"} = "sig";
    $potential{"H"} = "sig";
    $potential{"I"} = "sig";
    $potential{"K"} = "sig";
    $potential{"L"} = "sig";
    $potential{"M"} = "sig";
    $potential{"N"} = "sig";
    $potential{"P"} = "sig";
    $potential{"Q"} = "sig";
    $potential{"R"} = "sig";
    $potential{"S"} = "sig";
    $potential{"T"} = "sig";
    $potential{"V"} = "sig";
    $potential{"C"} = "flat";
    $potential{"W"} = "gauss";
    $potential{"Y"} = "gauss";
    $self->{pars}->{potential} = \%potential;

    my %e0;
    $e0{"A"} = -0.29;
    $e0{"D"} = 1.19;
    $e0{"E"} = 1.30;
    $e0{"F"} = -0.8;
    $e0{"G"} = -0.01;
    $e0{"H"} = 0.75;
    $e0{"I"} = -0.56;
    $e0{"K"} = 1.66;
    $e0{"L"} = -0.64;
    $e0{"M"} = -0.28;
    $e0{"N"} = 0.89;
    $e0{"P"} = 0.83;
    $e0{"Q"} = 1.21;
    $e0{"R"} = 1.55;
    $e0{"S"} = 0.10;
    $e0{"T"} = 0.01;
    $e0{"V"} = -0.47;
    $self->{pars}->{e0} = \%e0;

    my %zmid;
    $zmid{"A"} = 10.22;
    $zmid{"D"} = 14.25;
    $zmid{"E"} = 14.66;
    $zmid{"F"} = 19.67;
    $zmid{"G"} = 13.86;
    $zmid{"H"} = 12.26;
    $zmid{"I"} = 14.34;
    $zmid{"K"} = 11.11;
    $zmid{"L"} = 17.34;
    $zmid{"M"} = 18.04;
    $zmid{"N"} = 12.78;
    $zmid{"P"} = 18.09;
    $zmid{"Q"} = 10.46;
    $zmid{"R"} = 9.34;
    $zmid{"S"} = 13.86;
    $zmid{"T"} = 13.86;
    $zmid{"V"} = 11.35;
    $self->{pars}->{zmid} = \%zmid;

    my %n;
    $n{"A"} = 4.67;
    $n{"D"} = 8.98;
    $n{"E"} = 4.16;
    $n{"F"} = 7.12;
    $n{"G"} = 6.00;
    $n{"H"} = 2.77;
    $n{"I"} = 10.69;
    $n{"K"} = 2.09;
    $n{"L"} = 8.61;
    $n{"M"} = 7.13;
    $n{"N"} = 6.28;
    $n{"P"} = 3.53;
    $n{"Q"} = 2.59;
    $n{"R"} = 4.68;
    $n{"S"} = 6.00;
    $n{"T"} = 6.00;
    $n{"V"} = 4.97;
    $self->{pars}->{n} = \%n;

    my %emin;
    $emin{"W"} = -0.85;
    $emin{"Y"} = -0.42;
    $self->{pars}->{emin} = \%emin;

    my %zmin;
    $zmin{"W"} = 11.65;
    $zmin{"Y"} = 13.04;
    $self->{pars}->{zmin} = \%zmin;

    my %sig2;
    $sig2{"W"} = 51.84;
    $sig2{"Y"} = 38.44;
    $self->{pars}->{sig2} = \%sig2;

  } else {
    GENERAL::error("Uknown EZ version '$ver' requested!");
  }

  return $self;
}

=head2

 Title   :  scoreHelicalSequenceVertical
 Usage   :  my $score = $ez->scoreHelicalSequenceVertical("ACDEFGRTYKL...", 0.0)
            my $score = $ez->scoreHelicalSequenceVertical(\@seq, 0.0)
 Function:  EZ score for the sequence when inserted as a helix into the middle of a membrane, with a given dZ offset from the middle
 Returns :  EZ score
 Args    :  1. sequence string (single-letter code) or sequence array reference (either single or triple-letter code)
            2. optional: dZ offset from the middle of the membrane

=cut

sub scoreHelicalSequenceVertical {
  my $self = shift;
  my $seq = shift;
  GENERAL::requireArgs($self, $seq);
  my $dzoff = shift; $dzoff = 0 if (!defined($dzoff));

  my @seq;
  if (ref($seq) =~ /ARRAY/) {
    @seq = @$seq;
  } else {
    @seq = split("", $seq);
  }

  my $par = $self->{pars};
  my $ez = 0;
  my $d = 1.51; # helical rise per residue
  for (my $i = 0; $i < scalar(@seq); $i++) {
    my $z = abs(($i - (scalar(@seq) - 1)/2)*$d + $dzoff);
    my $aa = $seq[$i]; $aa = SEQUENCE::t2s($aa) if (length($aa) == 3);
    GENERAL::assert(defined($par->{potential}{$aa}), "potential not defined for amino acid '$aa'");
    if ($par->{potential}{$aa} eq "sig") {
      $ez += $par->{e0}{$aa}/(($z/$par->{zmid}{$aa})**($par->{n}{$aa}) + 1);
    } elsif ($par->{potential}{$aa} eq "gauss") {
      $ez += $par->{emin}{$aa} * exp(-(($z - $par->{zmin}{$aa})**2)/(2*$par->{sig2}{$aa}));
    }
  }

  return $ez;
}

=head2

 Title   :  insertSequenceVertical
 Usage   :  my $dz = $ez->insertSequenceVertical("ACDEFGRTYKL...", $dzmin, $dzmax)
 Function:  Finds the optimal delta Z offset for insering a sequence vertically
 Returns :  Optimal dZ offset
 Args    :  1. sequence string (single-letter code) or sequence array reference (either single or triple-letter code)
            2. optional: dZ lower bound value (default is -20)
            3. optimal: dZ upper bound value (default is +20)

=cut

sub insertSequenceVertical {
  my $self = shift;
  my $seq = shift;
  GENERAL::requireArgs($self, $seq);
  my $dzmi = shift; $dzmi = -20 if (!defined($dzmi));
  my $dzma = shift; $dzma = 20 if (!defined($dzma));

  my $N = 200; # number of interations
  my $nf = 4; # number of focusing runs (zooming ins)
  my $fwin = 20; # size of the focusing window
  my $bsc = undef; # best score
  my $bdz = undef; # best deltaZ
  for (my $fi = 1; $fi <= $nf; $fi++) {
    my @scores;
    for (my $i = 0; $i < $N; $i++) {
      push(@scores, $self->scoreHelicalSequenceVertical($seq, $dzmi + $i*($dzma - $dzmi)/($N-1)));
    }
    my $mi = -1; $bsc = GENERAL::min(\@scores, \$mi);
    $bdz = $dzmi + $mi*($dzma - $dzmi)/($N-1);
    $dzmi_n = GENERAL::max($dzmi + ($mi - $fwin/2)*($dzma - $dzmi)/($N-1), $dzmi);
    $dzma_n = GENERAL::min($dzmi + ($mi + $fwin/2)*($dzma - $dzmi)/($N-1), $dzma);
    $dzmi = $dzmi_n; $dzma = $dzma_n;
  }
  return $bdz;
}

1;
