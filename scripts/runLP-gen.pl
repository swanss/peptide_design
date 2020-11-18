use GENERAL;
use ROTAMER;
use DEFINITIONS;
use SEQUENCE;
use Getopt::Long;
use File::Spec;

my %opts = GetInput();
my $isspec = (scalar(@{$opts{ce}}) > 0);

# -- move to scratch space
GENERAL::createLocalSpace();
my $odir = GENERAL::GetDir();
GENERAL::cchdir($LOCAL_DEF);

# -- write data file for LP
my $ofh = GENERAL::GetOutFH("energy.dat");
print "Writing data file for LP model...\n";

# read/write energy table
my ($ns, $n, $L, $I, $sites) = formatEnergyData($opts{e}, $ofh, \%opts);

# read/write energy tables for any additional constraints
for (my $i = 0; $i < scalar(@{$opts{ce}}); $i++) {
  $opts{preI} = $I;
  my ($nsi, $ni) = (formatEnergyData($opts{ce}->[$i], $ofh, \%opts, 0, "$i"))[0, 1];
  GENERAL::assert(($nsi == $ns) && ($ni == $n), "the number of sites and self-energy terms in the main energy table and any constraint energy tables must be the same!");
}

# write number of solutions required
$ofh->printf("param iterations := 1 ;\n"); # number of solutions required
$ofh->printf("end;\n\n");
close($ofh);

# -- parse sequence constraints
my @sc; # the minimal/maximal number of matching amino acids and the constraining sequences
if (defined($opts{d})) {
  my $diffSeqs = GENERAL::file2array($opts{d});
  foreach my $diffSeq (@$diffSeqs) {
    my @arr = split(" ", $diffSeq);
    GENERAL::assert(scalar(@arr) == scalar(@$sites) + 2, "based on the energy file, expected the sequence '$diffSeq' to have " . scalar(@$sites) . " amino acids.");
    GENERAL::assert(GENERAL::isInteger($arr[0]) && GENERAL::isInteger($arr[1]), "expected two entries of '$diffSeq to have minimal and maximal number of identities required.");
    my @const = (shift @arr, shift @arr);
    my @seq;
    for (my $i = 0; $i < scalar(@$sites); $i++) {
      my $sn = $sites->[$i];
      my $aa = $arr[$i];
      next if ($aa eq "***");
      GENERAL::assert(defined($I->{$sn}{$aa}), "in parsing sequence constraint file '$opts{d}', did not find amino acid '$aa' at site '$sn' as allowed.");
      push(@const, $I->{$sn}{$aa});
    }
    push(@sc, \@const);
  }
}

# -- parse equivalence constraints
my @eqC;
if (defined($opts{'q'})) {
  my $ifh = GENERAL::GetInFH($opts{'q'});
  foreach my $cline (<$ifh>) {
    $cline = GENERAL::Trim($cline);
    next if ($cline eq "");
    my @cline = split(" ", $cline);
    GENERAL::assert(scalar(@cline) == 2, "expected two entries in equivalence constraint '$cline'");
    GENERAL::assert(defined($I->{$cline[0]}) && defined($I->{$cline[1]}), "could not find at least one site from equivalence constraint '$cline'");
    my $si = $I->{$cline[0]};
    my $sj = $I->{$cline[1]};
    GENERAL::assert(scalar(keys(%$si)) == scalar(keys(%$sj)), "the two sites from equivalence constraint '$cline' have different number of allowed amino acids");
    foreach my $aa (keys(%$si)) {
      # make sure the two sites have equivalent amino-acid choices
      GENERAL::assert(defined($sj->{$aa}), "the two sites from equivalence constraint '$cline' have different sets of allowed amino acids");
      my @eq = ($si->{$aa}, $sj->{$aa});
      push(@eqC, \@eq);
    }
    
  }
  close($ifh);
}

# -- parse generalized sequence constraints
my @gc; # the maximum number of matching amino acids and the constraining sequences
if (defined($opts{g})) {
  my $file = GENERAL::file2array($opts{g});
  foreach my $consLine (@$file) {
    $consLine =~ s/#.*$//g;
    $consLine = GENERAL::Trim($consLine); next if ($consLine eq "");
    # parse constraint
    my @consLine = split("\t", $consLine);
    GENERAL::assert(scalar(@consLine) == 4, "each generalized constraint in file $opts{g} should have exactly four space-separated fields!");
    my @cSites = split(" ", $consLine[0]);
    my @cAAs = split(" ", $consLine[1]);
    my $cmin = $consLine[2]; GENERAL::assert(GENERAL::isInteger($cmin), "lower bound in constraint '$consLine' must be integer!");
    my $cmax = $consLine[3]; GENERAL::assert(GENERAL::isInteger($cmax), "upper bound in constraint '$consLine' must be integer!");

    # encode constraint
    my @aaChoices;
    foreach my $sn (@cSites) {
      GENERAL::assert(GENERAL::findin($sn, @$sites) >= 0, "site $sn in constraint '$consLine' is not a valid design site!");
      foreach my $aa (@cAAs) {
        next if (!defined($I->{$sn}{$aa})); # if this amino acid is not allowed at this site to begin with, it will never contribute to the total count
        push(@aaChoices, $I->{$sn}{$aa});
      }
    }
    my %const;
    $const{aaChoices} = \@aaChoices;
    $const{min} = $cmin;
    $const{max} = $cmax;
    push(@gc, \%const);
  }
}

# -- run optimizer
# first try LP
my ($mathFile, $addConstNames) = getMathFile($isspec ? "cks_lp_gnu_spec.mod" : "cks_lp_gnu.mod", \@sc, \@gc, \@eqC, $opts{cb});
my ($E, $X) = optimize($mathFile, $n, $addConstNames, \%opts);
if ($X eq 0) { print "NO SOLUTIONS!\n"; die("\n"); }

# is the solution integral?
my $isint = 1;
foreach my $x (@$X) {
  if (!GENERAL::isInteger($x)) {
    $isint = 0;
  }
}

# if not, try ILP
if (! $isint) {
  print "LP solution not integral...\nTrying ILP...\n";
  ($mathFile, $addConstNames) = getMathFile($isspec ? "cks_ilp_gnu_spec.mod" : "cks_ilp_gnu.mod", \@sc, \@gc, \@eqC, $opts{cb});
#  GENERAL::csystem("$BIN_DEF/glpsol --math $mathFile --data energy.dat --output out");
#  ($E, $X) = parseOptimizerOutput("out", $n, @$addConstNames);
  ($E, $X) = optimize($mathFile, $n, $addConstNames, \%opts);
  if ($X eq 0) { print "NO SOLUTIONS!\n"; die("\n"); }
} else {
  print "LP solution is integral!\n";
}

# -- express solution as a usual string
my @sol = GENERAL::ones($ns, -1);
for (my $i = 1; $i <= scalar(@$X); $i++) {
  if ($X->[$i-1]) {
    my $p = $L->{$i}{pos};
    my $c = $L->{$i}{choice};
    $sol[$p] = $c;
  }
}
GENERAL::assert(GENERAL::findin(-1, @sol) == -1, "not all positions were assigned!");
my $solSeq = "";
if (defined($opts{s})) {
  foreach my $aa (@sol) { $solSeq .= SEQUENCE::t2s($aa); }
} else {
  $solSeq = join(" ", @sol);
}

# first print additional/custom constraints
print "\n";
for (my $i = 1; $i < scalar(@$E); $i++) {
  printf("CONST: $addConstNames->[$i-1]\t$E->[$i]\n");
}

# then the final solution and sequence
print "$E->[0] $solSeq\n";
print "$solSeq\n";

# -- clean up
GENERAL::cchdir($odir);
GENERAL::destroyLocalSpace();



sub optimize {
  my $mathFile = shift;
  my $n = shift;
  my $addVarNames = shift;
  my $opts = shift;
  my $solf = "out";

  if ($opts->{slv} ne "CPLEX") {
    GENERAL::csystem("$BIN_DEF/glpsol --math $mathFile --data energy.dat --output $solf");
    ($E, $X) = parseGLPSolOutput($solf, $n, @$addVarNames);
  } else {
    GENERAL::csystem("$BIN_DEF/glpsol --math $mathFile --data energy.dat --wlp problem.cplex --check");
    GENERAL::csystem("echo \"read problem.cplex lp\noptimize\nwrite $solf sol\nquit\n\" | /home/grigoryanlab/library/CPLEX_Studio1263/cplex/bin/x86-64_linux/cplex");
    ($E, $X) = parseCPLEXOutput($solf, $n, @$addVarNames);
  }
  GENERAL::crm($solf);
  return ($E, $X);
}

=head2

 Title   :  formatEnergyData
 Usage   :  formatEnergyData("energy.dat", $ofh, 1, "")
 Function:  Reads a flat-text formatted energy table and writes the data into an open data file for ILP.
 Returns :  A tuple: number of sites idetified in the erngy table, number of self energy terms identified in the energy table,
            and a hash reference with mappings between node variable indices, sites, and amino-acid choices.
 Args    :  1. Input flat-text energy table.
            2. File handle for the output data file.
            3. Hash reference with various parameters
            4. optional: whether to define the various size parameters needed for defining the problem. Default is yes.
            5. optional: name of the energy table; will name the associated vertex and edge costs accordingly. Default is empty string.
=cut

sub formatEnergyData {
  my $efile = shift;
  my $ofh = shift;
  my $opts = shift;
  my $writeSizeInfo = shift; $writeSizeInfo = 1 if (!defined($writeSizeInfo));
  my $name = shift; $name = "" if (!defined($name));
  my $ifh = GENERAL::GetInFH($efile);

  # read self energies
  my (%S, %I, %L); my $ns = 0; my $n = 0; my @sites;
  my $line = undef;
  while (1) {
    $line = <$ifh>; last if (!$line);
    $line = myTrim($line); next if ($line eq "");
    my @line = split(" ", $line);
    last if (scalar(@line) > 3);
    GENERAL::assert(scalar(@line) == 3, "expected line '$line' to be a self-energy line with 3 entries");
    my $s = $line[0];
    my $aa = $line[1];
    next if (defined($opts->{preI}) && !defined($opts->{preI}->{$s}{$aa}));
    if (!defined($S{$s})) {
      push(@sites, $s);
      $ns++;
    }
    my $e = pseudoEnergy($line[2], $opts->{l}, $opts->{z});

    if (!defined($opts->{preI}) && defined($opts->{n}) && defined($S{$s}) && (scalar(keys(%{$S{$s}})) >= $opts->{n})) {
      # if too many values at this position, find previous worst
      my ($me, $maa) = maxHashValue($S{$s});
      # if the new value is worse than that, skip
      next if ($me < $e);
      # otherwise, remove the old amino-acid choice, and add the new one in its place
      my $no = $I{$s}{$maa};
      delete($S{$s}{$maa}); delete($I{$s}{$maa});
      $S{$s}{$aa} = $e;
      $I{$s}{$aa} = $no;
      $L{$no}{choice} = $aa;
    } else {
      $n++;
      $S{$s}{$aa} = $e;
      $I{$s}{$aa} = defined($opts->{preI}) ? $opts->{preI}->{$s}{$aa} : $n;
      $L{$n}{choice} = $aa; $L{$n}{pos} = GENERAL::findin($s, @sites);
    }
  }

  # Write size info
  if ($writeSizeInfo) {
    $ofh->printf("data;\n\n");
    $ofh->printf("param num_posn := %d ;\n", $ns);
    $ofh->printf("param num_nodes := %d ;\n", $n);
    $ofh->printf("param posn_size :=\n");
    for (my $i = 0; $i < scalar(@sites); $i++) {
      $ofh->printf("%d %d\n", $i+1, scalar(keys(%{$S{$sites[$i]}})));
    }
    $ofh->printf(";\n\n");
  }

  # Write self energies
  $ofh->printf("param costV$name :=\n");
  foreach my $i (@sites) {
    # apply sort so that amino acids are always visited in the same order (even across different energy tables)
    foreach my $c (sort(keys(%{$S{$i}}))) {
      $ofh->printf("%d %f\n", $I{$i}{$c}, $S{$i}{$c});
    }
  }
  $ofh->printf(";\n\n");

  # Write pair energies
  $ofh->printf("param costE$name :=\n");
  while (defined($line)) {
    $line = myTrim($line);
    if ($line ne "") {
      my @line = split(" ", GENERAL::Trim($line));
      GENERAL::assert(scalar(@line) == 5, "expected line '$line' to be a pair-energy line with 5 entries");
      my $i = $line[0]; my $ci = $line[2];
      my $j = $line[1]; my $cj = $line[3];
      my $eij = pseudoEnergy($line[4], $opts->{l}, $opts->{z});
      if (!defined($opts->{c}) || (abs($eij) > $opts->{c})) {
        if (defined($I{$i}{$ci}) && defined($I{$j}{$cj})) {
          $ofh->printf("%d %d %f\n", $I{$i}{$ci}, $I{$j}{$cj}, $eij);
        } else {
          GENERAL::error("pair term ($i, $ci) - ($j, $cj) involves choices not encountered in the self part of the energy table") if (!defined($opts->{n}));
        }
      }
    }
    $line = <$ifh>;
  } 
  $ofh->printf(";\n\n");
  close($ifh);
  return ($ns, $n, \%L, \%I, \@sites);
}

sub GetInput {
  my $usage = GENERAL::usage("Sets up and runs Mona Singh's LP code for a general \"design\" optimization problem.",
                            "Required options:", "",
                            "-e", "energy file with self and pair terms. IMPORTANT: since nothing but this table is required, it must explicitly list ".
                                  "all self-energy terms, even those that are zero. Zero pair terms can (and should) be skipped. One term per line. The ".
                                  "format for self terms is 'S ABC ener', where S is some string uniquely identifying the site (can be an integer, but does not ".
                                  "have to be), ABC is the amino acid (need not be a three-letter code, again anything is fine), and ener is the corresponding ".
                                  "energy. Similarly, the format for pair terms is 'S1 S2 ABC DEF ener', meaning that the pair energy for amino acid ABC at site ".
                                  "S1 and amino acid DEF at site S2 is ener (obviously, the same naming convention for sites and amino acids must be used in ".
                                  "self and pair terms).",
                            "Optional:", "", "-c", "cutoff to apply to the absolute value of pair interactions (those below the cutoff are not considered).",
                            "-l", "flag: if specified, will apply -log() to the scores.",
                            "-d", "a file specifying sequences for constraining the optimal sequence. Each line should look like ".
                                  "'N1 N2 AA1 AA2 ... AAn', where AA1 through AAn are the amino acids of a sequence and N1 and N2 are integers designating ".
                                  "the minimal and maximal number of identities to this sequence that the optiaml sequence is allowed to have, respectively. Naturally, the ".
                                  "length of sequences in this file, and the amino acid names should be consistent with the energy table above. String '***' ".
                                  "in place of amino-acid name means that the corresponding position is skipped and not considered when counting identities.",
                            "-g", "a file with generalized sequence constraints. Each line, corresponding to one constraint, should consist of four tab-separated fields. ".
                                  "The first should be a space-separated list of positions, the second a space-separated list of amino acids, and the next two integers, corresponding to ".
                                  "the lower and upper bound. The interpretation is that the number of times the specified amino acids are used in the specified positions (in total), ".
                                  "should not exceed the upper bound and should not be below the lower bound.",
                            "-q", "a file specifying equivalent residues. Each line should have two sites (separated by space) to be constrained as ".
                                  "equivalent (i.e., to be occupied with the same amino acid). Sites should be the same as in the energy file.",
                            "--ce", "energy table for an additional variable to constrain (can be specified multiple times for multiple constraints)",
                            "--cb", "constraint bounds, should be specified as --cb 'min max', where min and max are the upper and lower bounds, respectively. ".
                                    "There should be as many --cb optionis specified as there are --ce options, with constraint boudns and energy tables ".
                                    "corresponding to each other in the order specified. If either the lower or upper bound are not needed, specify NA ".
                                    "instead of a numerical value.",
                            "-n", "limit the number of amino acids per position to this number. Will choose the best based on their self energy.",
                            "-z", "if -l is specified, and a number <= 0 is encountered, will set the result to the value given here. Otherwise, will throw an error.",
                            "-s", "if specified, will print the solution sequence in single-letter code.",
                            "--slv", "solver to use. Defaults to CPLEX studio, but specify glpk to use GLPK glpsol.");
  my %opts;
  $opts{ce} = newArrayRef(); $opts{cb} = newArrayRef();
  GetOptions (\%opts, "e=s", "c=s", "l", "d=s", "q=s", "n=s", "z=s", "ce=s" => $opts{ce}, "cb=s" => $opts{cb}, "g=s", "slv=s", "s");
  if (!(defined($opts{'e'}))) {
    die($usage);
  }
  GENERAL::assert(scalar(@{$opts{ce}}) == scalar(@{$opts{cb}}), "must have the same number of constraint energy tables as constraint bounds!");

  $opts{e} = File::Spec->rel2abs($opts{e});
  $opts{d} = File::Spec->rel2abs($opts{d}) if (defined($opts{d}));
  $opts{g} = File::Spec->rel2abs($opts{g}) if (defined($opts{g}));
  $opts{'q'} = File::Spec->rel2abs($opts{'q'}) if (defined($opts{'q'}));
  GENERAL::assert(GENERAL::isInteger($opts{n}) && ($opts{n} > 0), "-n must be a positive integer!") if (defined($opts{n}));
  GENERAL::assert(GENERAL::isNumeric($opts{z}), "-z must be numeric!") if (defined($opts{z}));
  $opts{slv} = "CPLEX" if (!defined($opts{slv}));
  $opts{slv} = uc($opts{slv});
  GENERAL::assert(($opts{slv} eq "GLPK") || ($opts{slv} eq "CPLEX"), "unknown solver '$opts{slv}'");

  for (my $i = 0; $i < scalar(@{$opts{ce}}); $i++) {
    GENERAL::assert(-e $opts{ce}->[$i] ? 1 : 0, "no such file '$opts{ce}->[$i]'");
    $opts{ce}->[$i] = File::Spec->rel2abs($opts{ce}->[$i]);
    my @b = split(" ", GENERAL::Trim($opts{cb}->[$i]));
    GENERAL::assert(scalar(@b) == 2, "bound specifications are expected to have two entries, but got '$opts{cb}->[$i]' instead");
    for (my $k = 0; $k < scalar(@b); $k++) {
      GENERAL::assert(GENERAL::isNumeric($b[$k]) || ($b[$k] eq "NA"), "bound values should be either numeric or NA, but got '$opts{cb}->[$i]' instead");
    }
    my %b; $b{min} = $b[0]; $b{max} = $b[1];
    $opts{cb}->[$i] = \%b;
  }
  return %opts;
}

sub myTrim {
  my $line = shift;
  $line =~ s/#.*$//g;
  $line = GENERAL::Trim($line);
  return $line;
}

sub pseudoEnergy {
  my $e = shift;
  my $takeLog = shift;
  my $zeroValue = shift;
  GENERAL::assert(GENERAL::isNumeric($e), "value '$e' is not numeric!");

  if (defined($takeLog)) {
    if ($e <= 0) {
      if (defined($zeroValue)) {
        return $zeroValue;
      } else {
        GENERAL::error("cannot take log of '$e'");
      }
    } else {
      return -log($e);
    }
  } else {
    return $e;
  }
}

# Parses the output of the optimizer and returns the values of node variables and total energy
sub parseCPLEXOutput {
  my $ofile = shift;
  my $N = shift;
  GENERAL::requireArgs($ofile, $N);
  my @varNames = @_; # names of additional variables to extract values for
  my $fh = GENERAL::GetInFH($ofile);
  my @varE = GENERAL::ones(scalar(@varNames), "N/A"); # for backwards compatibility; it's a pain to extra constraint activities from CPLEX

  # first, check if there is no solution
  my $st; GENERAL::assert(GENERAL::skipTo($fh, "solutionStatusString=", \$st), "could not find status line in file '$ofile'");
  return (0, 0) if ($st !~ /feasible solution|optimal/);
  $fh->seek(0, 0); GENERAL::assert(GENERAL::skipTo($fh, "objectiveValue=", \$st), "could not find objective value in file '$ofile'");
  GENERAL::assert($st =~ /objectiveValue="(.+)"/ ? 1 : 0, "could not parse out objective value from line '$st'");
  unshift(@varE, $1); # make sure to place energy as the first variable value

  # if there appears to be, parse it out
  my @X = GENERAL::ones($N, 0); # space for node variables
  GENERAL::assert(GENERAL::skipTo($fh, "variables"), "could not find the variables in file '$ofile'");
  foreach my $line (<$fh>) {
    if ($line =~ /variable name="X\(([^"]+)\)" .+ value="([^"]+)"/) {
      my $i = $1; my $val = $2;
      GENERAL::assert(GENERAL::isInteger($i) && GENERAL::isNumeric($val) && ($i <= $N), "value or index not numeric or out of range ($N) on line '$line'");
      $X[$i-1] = $val + 0.0;
    }
  }
  return (\@varE, \@X);
}


# Parses the output of the optimizer and returns the values of node variables and total energy
sub parseGLPSolOutput {
  my $ofile = shift;
  my $N = shift;
  GENERAL::requireArgs($ofile, $N);
  my @varNames = @_; # names of additional variables to extract values for
  unshift(@varNames, "energy"); # make sure energy is the first name
  my $fh = GENERAL::GetInFH($ofile);

  # first, check if there is no solution
  my $st; GENERAL::assert(GENERAL::skipTo($fh, "Status:", \$st), "could not find Status line in file '$ofile'");
  return (0, 0) if ($st =~ /UNDEFINED|INTEGER EMPTY/);

  # if there appears to be, parse it out
  my @X = GENERAL::ones($N, 0); # space for node variables
  my @varE = GENERAL::ones(scalar(@varNames), undef);
  GENERAL::assert(GENERAL::skipTo($fh, "-----"), "could not find header line in file '$ofile'");
  my $i = 0;
  foreach my $line (<$fh>) {
    for (my $k = 0; $k < scalar(@varNames); $k++) {
      if (!defined($varE[$k]) && ($line =~ /$varNames[$k]/)) {
        GENERAL::assert(length($line) >= 36, "could not parse activity info for variable $varNames[$k], '$line'");
        my $val = substr($line, 23, 13); GENERAL::assert(GENERAL::isNumeric($val), "could not parse a numeric activity value from line '$line'");
        $varE[$k] = $val + 0.0;
        next;
      }
    }
    next if ($line !~ /X\[/);
    my @arr = split(" ", $line);
    if ($i >= $N) { GENERAL::error("Too many node variables (only $N expected)!"); }
    $X[$i] = $arr[3] + 0.0;
    $i++;
  }
  return (\@varE, \@X);
}

sub getMathFile {
  my $origMathFile = shift;
  my $sc = shift;    # sequence constraints
  my $gc = shift;    # generalized sequence constraints
  my $ec = shift;    # equivalence constraints
  my $etabc = shift; # additional energy table constraints
  $sc = newArrayRef() if (!defined($sc));
  $gc = newArrayRef() if (!defined($gc));
  $ec = newArrayRef() if (!defined($ec));
  $etabc = newArrayRef() if (!defined($etabc));
  my @addConstNames; # will record names of the additional constraints to be explicitly displayed at the end
  my $ifh = GENERAL::GetInFH(DEFINITIONS::getParmFile($origMathFile));
  my $tmpf = "_math_file_tmp.mod";
  my $ofh = GENERAL::GetOutFH($tmpf);

  while (<$ifh>) {
    if (/# additional parameters/) {
      for (my $i = 0; $i < scalar(@$etabc); $i++) {
        if (((defined($etabc->[$i]->{min}) && ($etabc->[$i]->{min} ne "NA"))) || ((defined($etabc->[$i]->{max}) && ($etabc->[$i]->{max} ne "NA")))) {
          $ofh->printf("\nparam costV$i {V} default 0.0;\nparam costE$i {Efull} default 0.0;\n\n");
          push(@addConstNames, "etab$i\_c");
        }
      }
    } elsif (/# additional constraints/) {
      for (my $i = 0; $i < scalar(@$sc); $i++) {
        my $c = $sc->[$i];
        my $Nmin = $c->[0];
        my $Nmax = $c->[1];
        my @ii = @$c[2 .. scalar(@$c)-1];
        $ofh->printf("subject to notTooSimilar$i:\n");
        $ofh->printf("    X[" . join("] + X[", @ii) . "] <= $Nmax;\n\n");
        $ofh->printf("subject to notTooDisimilar$i:\n");
        $ofh->printf("    X[" . join("] + X[", @ii) . "] >= $Nmin;\n\n");
      }
      for (my $i = 0; $i < scalar(@$gc); $i++) {
        my $c = $gc->[$i];
        $ofh->printf("subject to genSeqConstUB$i:\n");
        $ofh->printf("    X[" . join("] + X[", @{$c->{aaChoices}}) . "] <= $c->{max};\n\n");
        $ofh->printf("subject to genSeqConstLB$i:\n");
        $ofh->printf("    X[" . join("] + X[", @{$c->{aaChoices}}) . "] >= $c->{min};\n\n");
      }
      for (my $i = 0; $i < scalar(@$ec); $i++) {
        my $si = $ec->[$i]->[0];
        my $sj = $ec->[$i]->[1];
        $ofh->printf("subject to equivalentSites$i:\n");
        $ofh->printf("    X[$si] = X[$sj];\n\n");
      }
      for (my $i = 0; $i < scalar(@$etabc); $i++) {
        if ((defined($etabc->[$i]->{min}) && ($etabc->[$i]->{min} ne "NA"))) {
          $ofh->printf("subject to etab$i\_cLB:\n    (sum {v in V} costV$i\[v] * X[v]) + (sum {(u,v) in E} costE$i\[u,v] * Y[u,v]) >= $etabc->[$i]->{min};\n\n");
        }
        if ((defined($etabc->[$i]->{max}) && ($etabc->[$i]->{max} ne "NA"))) {
          $ofh->printf("subject to etab$i\_cUB:\n    (sum {v in V} costV$i\[v] * X[v]) + (sum {(u,v) in E} costE$i\[u,v] * Y[u,v]) <= $etabc->[$i]->{max};\n\n");
        }
      }
    } elsif (/POSITIVE_PAIRS_FILL/) {
      my $expr = "";
      for (my $i = 0; $i < scalar(@$etabc); $i++) {
        if (((defined($etabc->[$i]->{min}) && ($etabc->[$i]->{min} ne "NA"))) || ((defined($etabc->[$i]->{max}) && ($etabc->[$i]->{max} ne "NA")))) {
          $expr .= " or (exists {u in C[i], v in C[j]} costE$i\[u,v] > 0)";
        }
      }
      $_ =~ s/POSITIVE_PAIRS_FILL/$expr/;
      $ofh->printf($_);
    } else {
      $ofh->printf($_);
    }
  }
  return ($tmpf, \@addConstNames);
}

sub maxHashValue {
  my $h = shift;
  my $mv = undef;
  my $mk = undef;
  foreach my $k (keys(%$h)) {
    if (!defined($mv) || ($mv < $h->{$k})) {
      $mv = $h->{$k};
      $mk = $k;
    }
  }
  return ($mv, $mk);
}

sub newArrayRef {
  my @arr;
  return \@arr;
}
