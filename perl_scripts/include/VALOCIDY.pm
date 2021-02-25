package VALOCIDY;

use DEFINITIONS;
use CHARMM;
use GENERAL;
use NAMD;
use TieArrayC;
use Math::Libm ':all';
use Inline C => 'DATA' => CCFLAGS => '-O3';
#use Math::Random::Secure;

# Implements various functions based on the VALOCIDY approach (VAluation of LOcal COnfiguration Integral by DYnamics),
# includeing the MODISTE approach for global partition function estimation
# (MOlecular Dynamics-based Integration with the assistance of Sampling with Temperature Elevation)

=head2

 Title   :  integrateMinimum
 Usage   :
 Function:  Given a dihedral-based definition of a region of phase space, estimates its total partition function using
            molecular dynamics and utilizing the "inverse" importance sampling approach. Implements the VALOCIDY algorithm.
 Returns :  Returns an estimate of the total partition function as well as the variance of the estimate.
 Args    :  1. A pre-setup CHARMM structure with the system loaded, ideally already within the desired region of phase space.
            2. A reference to an array of lower-bound IC values.
            3. A reference to an array of upper-bound IC values.
            4. A reference to array of IC definitions. For each IC, the entry is an array of four strings which are interpreted as CHARMM selections.
               corresponding to the four atoms defining the IC.
            5. Hash reference defining various parameters of the dynamics. Currently, the following parameters are supported (hash key - explanation):
               ne   - number of equilibration steps in femtoseconds (before integration begins)
               nc   - number of data-collections steps in femtoseconds
               ns   - structure sampling interval in femtoseconds
               T    - simulation temperature
               seed - optional: value of random seed for repeating realizations (picked randomly from time and pid by default)
               minN - optional: number of minimization steps before beginning MD (default 10000)

=cut

sub integrateMinimum {
  my $chrm = shift;
  my $Cdef = shift;
  my $parMD = shift;
  GENERAL::requireArgs($chrm, $Cdef, $parMD);
  if (defined($parMD->{dbg}) && defined($parMD->{startPDB})) {
    GENERAL::SaveStruct($Cdef, GENERAL::GetBase($parMD->{startPDB}) . "_cdef.out");
  }
  my $time = time();
  my $measureTime = (defined($parMD->{dbg}) && ($parMD->{dbg}));

  # setup some default parameters
  $parMD->{seed} = $$ + int(rand()*100000) if (!defined($parMD->{seed}));
  $parMD->{T} = 298.15 if (!defined($parMD->{T}));
  my $waitC = 0.001*scalar(@$Cdef)/10; $waitC = $parMD->{waitC} if (defined($parMD->{waitC}));
  my $R = 1.9858775/1000;
  my $echeck = ""; $echeck = "ECHECK $parMD->{echeck}" if (defined($parMD->{echeck}));

  # Limit the dynamics of the system to within the desired area of phase space by setting up steep bounding potentials
  # just outside of the allowed range of phase-space coordinates (will call this r hereafter).
  setupCoordinateRestraints($chrm, $Cdef, $parMD);

  $chrm->send("prnlev 5\nwrnlev 5\nic fill\nprint ic\nenergy\n")  if (defined($parMD->{dbg}) && $parMD->{dbg});
  if (defined($parMD->{dbg}) && defined($parMD->{startPDB})) {
    $chrm->writeCRD("premin.crd"); $chrm->wait('./', 1, 300, 0); GENERAL::csystem("cp premin.crd " . GENERAL::GetBase($parMD->{startPDB}) . "_premin.crd");
    my $pdb = PDB::new("premin.crd"); $pdb->writePDB(GENERAL::GetBase($parMD->{startPDB}) . "_premin.pdb", "");
  }

  my $outdcd = "_datacoll_run.dcd";
  # Run long dynamics
  if (!defined($parMD->{traj})) {
    $chrm->send(sprintf("mini sd nstep %d\n", (defined($parMD->{minN}) ? $parMD->{minN} : 10000))); # start with a locally-minimized low-energy structure
    if (defined($parMD->{startPDB})) {
      $chrm->writeCRD("startsnap.crd"); $chrm->wait('./', 1, 300, 0);
      my $pdb = PDB::new("startsnap.crd"); $pdb->writePDB($parMD->{startPDB}, "");
    }

    # equilibration phase
    my $rstf = "_equil_restart.rst";
    $chrm->send("open unit 49 write card name \"$rstf\"\n");
    $chrm->send(constantTemperatureDynamicsCHARMM('temp', $parMD->{T}, 'steps', $parMD->{ne}, 'echeck', $echeck, 'seed', $parMD->{seed}, 'IUNWRI', 49, 'step', $parMD->{step}));
    $parMD->{seed}++;

    # data collection phase
    $chrm->send("open unit 49 read card name \"$rstf\"\n");
    my $outene = "_datacoll_run.ener";
    my $nrstf = "_datacoll_1.rst";
    $chrm->send("open unit 50 write file name \"$outdcd\"\n");
    $chrm->send("open unit 51 write card name \"$outene\"\n");
    $chrm->send("open unit 52 write card name \"$nrstf\"\n");
    $chrm->send(constantTemperatureDynamicsCHARMM('temp', $parMD->{T}, 'steps', $parMD->{nc}, 'echeck', $echeck, 'seed', $parMD->{seed}, 'restart', 1,
                                                  'NSAVC', $parMD->{ns}, 'IUNREA', 49, 'IUNWRI', 52, 'IUNCRD', 50, 'KUNIT', 51, 'step', $parMD->{step}));
    $chrm->closeCard(50);
    $chrm->closeCard(51);
    $chrm->closeCard(52);
    $chrm->wait('./', 1, 1000, 0); # wait until the buffer is clear (data-collection run is complete)
    $chrm->send("prnlev 5\nwrnlev 5\nic fill\nprint ic\nenergy\n")  if (defined($parMD->{dbg}) && $parMD->{dbg});
  } else {
    GENERAL::csystem("cp $parMD->{traj} $outdcd");
  }
  if ($measureTime) { printf("--> TIME: setup and trajectory generation took %d seconds\n", time() - $time); $time = time(); }

  # Analyze run: extract energies and phase-space coordinates
  printf("Before analyzeRun: " . `grep VmPeak /proc/$$/status`) if (defined($parMD->{dbg}) && ($parMD->{dbg} > 1));
  my $ret = analyzeRun($chrm, $outdcd, $parMD->{dbg}, "?GEO", $Cdef);
  printf("After analyzeRun: " . `grep VmPeak /proc/$$/status`) if (defined($parMD->{dbg}) && ($parMD->{dbg} > 1));;
  if ($measureTime) { printf("--> TIME: analyzeRun took %d seconds\n", time() - $time); $time = time(); }

  my $ener = $ret->{ener}; my $C = $ret->{C};
  $chrm->send("MMFP\nGEO RESET\nEND\n\n");

  # determine the visited range for those degrees of freedom that were not fixed, but are to be integrated over
  for (my $ci = 0; $ci < scalar(@$Cdef); $ci++) {
    my $cdef = $Cdef->[$ci];
    if (defined($cdef->{determineRange}) && $cdef->{determineRange}) {
      my $mi = $ret->{C}->[$ci]->[0]; my $ma = $mi;
      if (($cdef->{type} =~ /^(DIHE|IMPR|ANGL)$/) && (!defined($cdef->{equilVal}))) {
        $cdef->{equilVal} = VALOCIDY::angularMean($ret->{C}->[$ci]);
      }
      foreach my $snv (@{$ret->{C}->[$ci]}) {
        if ($cdef->{type} =~ /^(DIHE|IMPR|ANGL)$/) {
          $mi = VALOCIDY::angleDiff(GENERAL::min(VALOCIDY::angleDiff($mi, $cdef->{equilVal}), VALOCIDY::angleDiff($snv, $cdef->{equilVal})), -$cdef->{equilVal});
          $ma = VALOCIDY::angleDiff(GENERAL::max(VALOCIDY::angleDiff($ma, $cdef->{equilVal}), VALOCIDY::angleDiff($snv, $cdef->{equilVal})), -$cdef->{equilVal});
        } else {
          $mi = GENERAL::min($mi, $snv);
          $ma = GENERAL::max($ma, $snv);
        }
      }
      $Cdef->[$ci]->{min} = $mi; $Cdef->[$ci]->{max} = $ma;
    }
  }
  printf("Before integrateIIS: " . `grep VmPeak /proc/$$/status`) if (defined($parMD->{dbg}) && ($parMD->{dbg} > 1));
  if ($measureTime) { printf("--> TIME: pre-integration setup took %d seconds\n", time() - $time); $time = time(); }

  # integrate
  my $ans = integrateIIS($Cdef, $ret->{C}, $ret->{ener}, %{$parMD});
  printf("After integrateIIS: " . `grep VmPeak /proc/$$/status`) if (defined($parMD->{dbg}) && ($parMD->{dbg} > 1));
  if ($measureTime) { printf("--> TIME: integration took %d seconds\n", time() - $time); $time = time(); }
  if (wantarray()) {
    return ($ans->{G}, $ans->{G}, $ans->{Q}, $ans->{varQlog}, $ans->{logQ});
  } else {
    return $ans;
  }

}


=head2

 Title   :  globalSampling
 Usage   :
 Function:  Implements the MODISTE approach in dihedral phase space.
 Returns :
 Args    :  1. A pre-setup CHARMM structure with the system loaded.
            2. A reference to array of IC definitions used to defined phase space for wide sampling. For each IC, the entry is an array of four strings which are
               interpreted as CHARMM selections corresponding to the four atoms defining the IC.
            3. Hash reference defining various parameters of the dynamics. Currently, the following parameters are supported (hash key - explanation):
               ne       - number of equilibration steps in femtoseconds (before sampling begins; during this time the energy/temperature equilibrate)
               nc       - number of data-collections steps in femtoseconds (the amount of actual high-temperature sampling)
               ns       - structure sampling interval in femtoseconds
               T        - simulation temperature
               seed     - optional: value of random seed for repeating realizations (picked randomly from time and pid by default)
               addcoors - optional: additional degrees of freedom, in the same format as before, for integrating locally. Each additional degree of freedom
                          contains hash keys "min" and "max" that define the range to allow these to vary during local integration.

=cut

sub globalSampling {
  my $chrm = shift;
  my $ICdef = shift;
  my $parSYM = shift;

  # setup some default parameters
  $parSYM->{seed} = $$ + int(rand()*100000) if (!defined($parSYM->{seed}));
  my $R = 1.9858775/1000;
  my $PI = M_PI;

  # if min/max are specified for global integration coordinates, then these need to be restrained (we are really doing a more focused global integration)
  setupCoordinateRestraints($chrm, $ICdef, $parSYM);

  # equilibration phase
  my $rstf = "_ws_equil_restart.rst";
  $chrm->send("open unit 49 write card name \"$rstf\"\n");
  $chrm->send("
      DYNA LEAP VERLET STRT NSTEP $parSYM->{nhe} TIMESTEP 0.001 -
      IPRFRQ 2000 IHTFRQ 0 IEQFRQ 0 NTRFRQ 100  -
      IUNREA -1 IUNWRI 49 IUNCRD -1 IUNVEL -1 KUNIT -1 -
      NPRINT 1000 NSAVC 100 NSAVV 0 IHBFRQ 0 INBFRQ 25  -
      TCONst  TCOUpling 5.0  TREFerence $parSYM->{Th} ISEED $parSYM->{seed} -
      FIRSTT $parSYM->{Th} FINALT $parSYM->{Th}  -
      IASORS 0 IASVEL 1 ISCVEL 0 ICHECW 0 TWINDH 0.0 TWINDL 0.0\n");
  $parSYM->{seed}++;

  # data collection phase
  $chrm->send("open unit 49 read card name \"$rstf\"\n");
  my $outdcd = "_ws_datacoll_run.dcd";
  my $outene = "_ws_datacoll_run.ener";
  my $nrstf = "_ws_datacoll_1.rst";
  $chrm->send("open unit 50 write file name \"$outdcd\"\n");
  $chrm->send("open unit 51 write card name \"$outene\"\n");
  $chrm->send("open unit 52 write card name \"$nrstf\"\n");
  $chrm->send("
    DYNA LEAP VERLET REST NSTEP $parSYM->{nhc} TIMESTEP 0.001 -
    IPRFRQ 2000 IHTFRQ 0 IEQFRQ 0 NTRFRQ 100  -
    IUNREA 49 IUNWRI 52 IUNCRD 50 IUNVEL -1 KUNIT 51 -
    NPRINT 1000 NSAVC $parSYM->{nhs} NSAVV 0 IHBFRQ 0 INBFRQ 25  -
    TCONst  TCOUpling 5.0  TREFerence $parSYM->{Th} ISEED $parSYM->{seed} -
    FIRSTT $parSYM->{Th} FINALT $parSYM->{Th}  -
    IASORS 0 IASVEL 1 ISCVEL 0 ICHECW 0 TWINDH 0.0 TWINDL 0.0\n");
  $chrm->closeCard(50);
  $chrm->closeCard(51);
  $chrm->closeCard(52);
  $chrm->wait('./', 1, 1000, 0); # wait until the buffer is clear (data-collection run is complete)
  $chrm->send("MMFP\nGEO RESET\nEND\n\n");

  # Analyze run: extract energies and ICs
  my $ret = analyzeRun($chrm, $outdcd, undef, undef, $ICdef);
  my $ener = $ret->{ener}; my $IC = $ret->{C};

  my $N = scalar(@$ener); # number of snapshots
  my $m = scalar(@$IC); # number of phase-space defining coordinates
  my $V = 360**$m; my $logV = $m*log(360); # total volume of phase space

  # Decide on delta's for each phase-space coordinate for expanding sampled structures
  my $delA = 30;
  $delA = $parSYM->{delA} if (defined($parSYM->{delA}));
  printf("Phase space around sampled structures will be expanded by +/- %f degrees\n", $delA);
  my $Vi = $delA**$m; my $logVi = $m*log($delA); # volume of phase-space elements that will be locally integrated

  # Generalized Monte-Carlo Integration
  my (@A, @Qa, @Hloc); my $Qa = 0; my $logQa = undef; my $nL = $parSYM->{N};
  # -- calculate the normalization constant for our high temp trajectory-derived PDF
  my $sig = 20; # steepsness (in degrees) of the PDF's drop off away from high-temperature sampled structures
  $sig = $parSYM->{sk} if (defined($parSYM->{sk}));
  my $Ck = 1/($sig * sqrt(2*$PI) * erf(180/($sig*sqrt(2))));

  # add any additional degrees of freedom for local integration
  if (defined($parSYM->{addcoors})) {
    push(@$ICdef, @{$parSYM->{addcoors}});
  }
  # save the original min/max, if defined
  for (my $ci = 0; $ci < scalar(@$ICdef); $ci++) {
    $ICdef->[$ci]->{gmin} = $ICdef->[$ci]->{min} if (defined($ICdef->[$ci]->{min}));
    $ICdef->[$ci]->{gmax} = $ICdef->[$ci]->{max} if (defined($ICdef->[$ci]->{max}));
  }

  # -- sample with this PDF
  my $i = 0;
  while ($i < $nL) {
    # pick a random snapshot
    my $sni = GENERAL::irand(0, $N-1);

    # pick a normally-distributed random displacement vector
    my @p = normalRandomVector($m, 0, $sig, 180); my @dp;
    for (my $ci = 0; $ci < $m; $ci++) {
      $p[$ci] = VALOCIDY::angleDiff($IC->[$ci]->[$sni], -$p[$ci]);
      push(@dp, VALOCIDY::angleDiff($p[$ci], $IC->[$ci]->[$sni]));
    }

    # compute the PDF of this pick
    my $pi = 0;
    for (my $snk = 0; $snk < $N; $snk++) {
      my $tmp_p = 1;
      for (my $ci = 0; $ci < $m; $ci++) {
        $tmp_p *= $Ck * exp(-((VALOCIDY::angleDiff($p[$ci], $IC->[$ci]->[$snk]))**2)/(2*($sig**2)));
      }
      $pi += $tmp_p;
    }
    $pi = $pi/$N;

    # create limits of phase-space defining variables (for the additional coordinates,
    # their limits were either defined outside of this function or were left undefined to integrate over the observed range)
    my $ok = 1;
    for (my $ci = 0; $ci < $m; $ci++) {
      $ICdef->[$ci]->{min} = VALOCIDY::angleDiff($p[$ci], $delA);
      $ICdef->[$ci]->{max} = VALOCIDY::angleDiff($p[$ci], -$delA);
      # if doing a focused global sampling, make sure the local area of phase space is fully contained within the total area we are interested in
      if (defined($ICdef->[$ci]->{gmin}) && (!VALOCIDY::angleWithin($ICdef->[$ci]->{min}, $ICdef->[$ci]->{gmin}, $ICdef->[$ci]->{gmax}) || !VALOCIDY::angleWithin($ICdef->[$ci]->{min}, $ICdef->[$ci]->{gmin}, $ICdef->[$ci]->{gmax}))) { $ok = 0; last; }
    }
    next if ($ok == 0);

    $parSYM->{startPDB} = sprintf("%s\_V%04d.pdb", $parSYM->{o}, $i+1) if (defined($parSYM->{o}));

    # if picked, compute the local low-temperature partition function corresponding to this snapshot
    printf("> On iteration %d picked snapshot %d with displacement vector norm %f (degrees), probability density %e...\n", $i+1, $sni+1, norm(@dp), $pi);
    windToSnapshot($chrm, $outdcd, $sni+1);                                                               # wind to the picked snapshot
    my $lcrdf = "local$i.crd"; $chrm->writeCRD($lcrdf); $chrm->wait();                                    # dump snahshop to a CRD
    my $lchrm = VALOCIDY::setupCHARMM($parSYM, "icrdf", $lcrdf, "chcmdlog" => "local$i.inp", "rnd" => 0); # start a clean CHARMM instance, load the local the CRD and integrate locally
    my $a = integrateMinimum($lchrm, $ICdef, $parSYM);                                                    # local integration
    if (defined($parSYM->{o}) && defined($parSYM->{dbg})) {
      GENERAL::csystem("cp _datacoll_run.dcd \"$parSYM->{o}.$i.dcd\"");
    }
    $lchrm->finish();
    #GENERAL::csystem("cp local$i.inp " . GENERAL::GetBase($parSYM->{startPDB}) . ".inp") if ((defined($parSYM->{o})) && defined($parSYM->{dbg}));

    # store visited point
    push(@A, $a);
    $Qa += $a->{Q}/$pi/$Vi; # quantity to average
    push(@Qa, $a->{Q}/$pi/$Vi);
    if (defined($logQa)) { $logQa = logSumExp($logQa, $a->{logQ} - log($pi) - $logVi); }
    else { $logQa = $a->{logQ} - log($pi) - $logVi; }
    push(@logQa, $a->{logQ} - log($pi) - $logVi);
    push(@Hloc, $a->{H}); # local enthalpy
    push(@Hw, $a->{logQ} - $logVi - log($pi)); # log of numerators of weights for local enthalpy for calculating global enthalpy (the denomenator is global partition function)
    printf(">>> Current estimate of Qall = %e (-RTlogQ = %f kcal/mol), Hall = %f kcal/mol\n\n", $Qa/($i+1), -$R*($parSYM->{T})*($logQa - log($i+1)), weighByLogFrac(\@Hw, $logQa-log($i+1), \@Hloc));
    $i++;
  }
  $Qa /= $nL;
  $logQa -= log($nL);
  my $H = weighByLogFrac(\@Hw, $logQa, \@Hloc);

  if (wantarray()) {
    return ($Qa, \@Qa, $logQa, \@logQa);
  } else {
    my %ret;
    $ret{Q} = $Qa;
    $ret{logQ} = $logQa;
    $ret{Qarr} = \@Qa;
    $ret{logQarr} = \@logQa;
    $ret{H} = $H;
    return \%ret;
  }

}


=head2

 Title   :  integrateIIS
 Usage   :
 Function:  Computes the local configurational integral structures sampled with the Boltzmann PDF, using Inverse Importance Sampling
 Returns :
 Args    :  1. Array of definition of phase-space variables
            2. Low bounds of coordinates, defining the lower bound of phase-space region to integrate
            3. Upper bounds of coordinates, defining the upper bound of phase-space region to integrate
            4. Table representing the trajectory of coordinate values. Array of arrays, first by coordinate, then by snapshot.
               This table will be modified, in place, to remove snapshots that are outside of the desired area of phase space.
            5. Table representing the energy trajectory. Each row has two elements - total energy of the system, and restraint energy portion
            2. Hash reference defining various parameters. Currently, the following parameters are supported (hash key - explanation):
               T    - REQUIRED: simulation temperature
               dbg  - optional: set the debug flag

=cut

sub integrateIIS {
  my $Cdef = shift;
  my $C = shift;
  my $ener = shift;
  my %par = @_;
  GENERAL::requireArgs($Cdef, $C, $ener);
  GENERAL::assert(defined($par{T}), "Temperature must be defined for integration!");
  my $R = 1.9858775/1000;
  my $PI = M_PI; # pi in radians
  if (!defined($par{refDensityType})) { $par{refDensityType} = 'singleGaussian'; }

  # As a pretty aggressive debugging step, replace the actually observed trajectory of structures with a uniformly
  # distributed one (in terms of the degrees of freedom integrated over) with a flag energy. In this case, if
  # everything is working well, the resultant Q should just correspond to the volume of phase space integrated over.
  if (defined($par{dbg}) && ($par{dbg} > 15)) {
    printf("!!!!!!! Aggressive debugging -- will replace the observed trajectory with a uniformly distributed one !!!!!!!\n");
    uniformlyDistributedTrajectory($Cdef, $C, $ener, 1);
  }

  # Go through all snapshots and remove those that are not strictly within the specified area of phase space.
  my $Nsn = scalar(@$ener);  # number of snapshots
  my $Nco = scalar(@$Cdef); # number of coordinates
  my (@ener, @rener, @sni);
  for (my $sni = 0; $sni < $Nsn; $sni++) {
    my $f = 1; # flag of whether this snapshot is within the IC range prescribed
    for (my $ci = 0; $ci < $Nco; $ci++) {
      if (!coordinateWithinLimits($Cdef->[$ci], $C->[$ci]->[$sni], $Cdef->[$ci]->{min}, $Cdef->[$ci]->{max})) {
        # the snapshot failed, so set the flag and remove this snapshot from the other arrays
        $f = "$ci: $Cdef->[$ci]->{type} $C->[$ci]->[$sni] ought to be between $Cdef->[$ci]->{min} and $Cdef->[$ci]->{max}"; last;
      }
    }
    # if this snapshot passed, push the IC from it onto the matrix of ICs that is by IC index first and snapshot index second
    if ($f eq 1) {
      push(@ener, $ener->[$sni]->[0]);
      push(@rener, $ener->[$sni]->[1]);
      push(@sni, $sni);
    } else {
      GENERAL::assert(GENERAL::isNumeric($ener->[$sni]->[1]) && ($ener->[$sni]->[1] gt 0), "Snapshot $sni is outside of desired phase-space area ($f), but bounding potential is not above zero ($ener->[$sni]->[1])!");
    }
  }
  if (defined($par{dbg}) && $par{dbg} && ((scalar(@ener) == 0) || $par{dbg} > 4)) {
    for (my $sni = 0; $sni < $Nsn; $sni++) {
      for (my $ci = 0; $ci < $Nco; $ci++) {
        if (!coordinateWithinLimits($Cdef->[$ci], $C->[$ci]->[$sni], $Cdef->[$ci]->{min}, $Cdef->[$ci]->{max}) || ($par{dbg} > 4)) {
          # the snapshot failed, so set the flag and remove this snapshot from the other arrays
          printf("]]] $sni: [$Cdef->[$ci]->{min} $Cdef->[$ci]->{max}] $C->[$ci]->[$sni] '$Cdef->[$ci]->{type}' : ");
          print(join(" x ", @{$Cdef->[$ci]->{atoms}})) if (uc($Cdef->[$ci]->{type}) =~ /^(DIHE|DIST|AXISDIST|PLANEDIST|ANGL|IMPR)$/);
          printf(", axis [%s] - [%s]", join(" ", @{$Cdef->[$ci]->{axisA}}), join(" ", @{$Cdef->[$ci]->{axisB}})) if (uc($Cdef->[$ci]->{type}) =~ /^AXISDIST$/);
          printf(", axis [%s] - [%s]", join(" ", @{$Cdef->[$ci]->{planePoint}}), join(" ", @{$Cdef->[$ci]->{planeNorm}})) if (uc($Cdef->[$ci]->{type}) =~ /^PLANEDIST$/);
          printf("\n");
        }
      }
    }
  }
  GENERAL::assert(scalar(@ener), "For some reason, after dynamics, no structures were found within the desired area of phase space!");
  $Nsn = scalar(@ener);
  for (my $ci = 0; $ci < $Nco; $ci++) {
    my @newarr; tie @newarr, TieArrayC => "double"; $#newarr = scalar(@sni)-1;
    for (my $i = 0; $i < scalar(@sni); $i++) {
      $newarr[$i] = $C->[$ci]->[$sni[$i]];
    }
    $C->[$ci] = \@newarr;
  }

  # Build a function f to integrate based on the observed evolution of r.
  # It depends only on r_i the lowest-energy structure (to place the minimum correctly) and the standard deviation of r_i (to make it vary in a relevant way)
  my $li; my @Cstd = GENERAL::ones($Nco, 0);
  GENERAL::min(\@ener, \$li); # lowest-energy structure index
  printf("REPORT: %d/%d structures within requested region, out of which $ener[$li] is the lowest energy...\n", scalar(@ener), scalar(@$ener));
  my $minE = GENERAL::max(\@ener); for (my $sni = 0; $sni < scalar(@ener); $sni++) { $ener[$sni] -= $minE; } # for fewer problems with large/small numbers, factor out the lowest energy
  printf("Will subtract reference energy %f\n", $minE) if (defined($par{dbg}) && ($par{dbg} > 1));
  my @Cm = GENERAL::ones($Nco, 0); # will hold the means of each Gaussian
  for (my $ci = 0; $ci < $Nco; $ci++) {
    next if (defined($Cdef->[$ci]->{integrate}) && ($Cdef->[$ci]->{integrate} == 0)); # certain variables are only for defining the relevant phase space and not for integrating
    # by default, integrate around the mean
    if (uc($Cdef->[$ci]->{type}) =~ /^(DIHE|IMPR)$/) {
      $Cm[$ci] = VALOCIDY::angularMean($C->[$ci]);
#      $Cm[$ci] = $C->[$ci]->[$li];
    } else {
      $Cm[$ci] = GENERAL::mean($C->[$ci]);
#      $Cm[$ci] = $C->[$ci]->[$li];
    }
    
    # for dihedral degrees of freedom, map all values to [mu - 180; mu + 180], where mu is the mean value around which integration happens for the given df
    # doing this makes it such that angleDiff(x, mu) (where x is some value from the trajectory) just becomes x - mu
    if (uc($Cdef->[$ci]->{type}) =~ /^(DIHE|IMPR)$/) {
      for (my $si = 0; $si < scalar(@{$C->[$ci]}); $si++) {
        $C->[$ci]->[$si] = VALOCIDY::angleDiff($C->[$ci]->[$si], $Cm[$ci]) + $Cm[$ci];
      }
      $Cdef->[$ci]->{max} = VALOCIDY::angleDiff($Cdef->[$ci]->{max}, $Cm[$ci]) + $Cm[$ci];
      $Cdef->[$ci]->{min} = VALOCIDY::angleDiff($Cdef->[$ci]->{min}, $Cm[$ci]) + $Cm[$ci];
    }
    if (!defined($Cdef->[$ci]->{sigma})) {
      $Cstd[$ci] = findOptimalSigma($Cdef->[$ci], $Cm[$ci], $C->[$ci]);
      printf("OPTIMAL SIGMA $Cstd[$ci] for $ci (mean $Cm[$ci], range [$Cdef->[$ci]->{min} $Cdef->[$ci]->{max}]), $Cdef->[$ci]->{type}\n") if (defined($par{dbg}));
    } else {
      $Cstd[$ci] = $Cdef->[$ci]->{sigma};
    }
    if ($Cstd[$ci] eq 0) {
      GENERAL::warning("Coordinate $ci, of type $Cdef->[$ci]->{type} did not detectably vary (sigma = $Cstd[$ci]), so will not integrate over it: " . join(" x ", @{$Cdef->[$ci]->{atoms}}));
      $Cdef->[$ci]->{integrate} = 0;
    }
  }

  # pick some analytcal function F whose integral will be computed in two ways - analytically or using the points we've sampled via the Boltzmann probability density 
  my $n = scalar(@ener);
  my $I = 1; # the value of the integral of function F we are "aiming" to compute using the Boltzmann probability as the biasing function
  my $logI = 0; # the log of the integral
  my @f = GENERAL::ones($n); # array of values of function f at the sampled structures
  my @logf = GENERAL::ones($n, 0); # array of logs of values of function f at the sampled structures

  if ($par{refDensityType} eq 'singleGaussian') {
    # Compute the analytical integral of f over all of the region of phase space
    for (my $ci = 0; $ci < $Nco; $ci++) {
      next if (defined($Cdef->[$ci]->{integrate}) && ($Cdef->[$ci]->{integrate} == 0));
      my $s = $Cstd[$ci]; my $f;
      if (uc($Cdef->[$ci]->{type}) =~ /^(DIHE|IMPR)$/) {
        # The function we are integrating is a Gaussian of angles, with the difference defined via angleDiff() (that's the function we use below for computing
        # function values at specific points). By defining the difference to be the angular difference (which makes sense from the perspective of the potential),
        # we essentially cause the number axis to wrap around. If we measure differences with some reference angle a0, this wrapping around happens at a0+pi. So,
        # if our integration ranges are defined in such a way that when going from the min to the max angle we span a0+pi, then the integral is really split into
        # two fragments that have to be computed separately.
        if (isCircularGaussianSplit($Cdef->[$ci], $Cm[$ci])) {
          $f = ($s*sqrt($PI/2)) * (erf(180/(sqrt(2)*$s)) + erf(($Cm[$ci] - $Cdef->[$ci]->{min})/(sqrt(2)*$s))) +
               ($s*sqrt($PI/2)) * (erf(($Cdef->[$ci]->{max} - $Cm[$ci])/(sqrt(2)*$s)) + erf(180/(sqrt(2)*$s)));
        } else {
          $f = ($s*sqrt($PI/2)) * (erf(($Cdef->[$ci]->{max} - $Cm[$ci])/(sqrt(2)*$s)) + erf(($Cm[$ci] - $Cdef->[$ci]->{min})/(sqrt(2)*$s)));
        }
      } elsif (uc($Cdef->[$ci]->{type}) =~ /^(DIST|AXISDIST|PLANEDIST|INPLACEDIST|ANGL)$/) {
        $f = ($s*sqrt($PI/2)) * (erf(($Cdef->[$ci]->{max} - $Cm[$ci])/(sqrt(2)*$s)) - erf(($Cdef->[$ci]->{min} - $Cm[$ci])/(sqrt(2)*$s)));
      } else {
        GENERAL::error("Unrecognized coordinate type '$Cdef->[$ci]->{type}'");
      }
      $I = $I * $f;
if ($f <= 0) { printf("ERROR: f = $f\n"); printf("ci = $ci\n"); printf("type = $Cdef->[$ci]->{type}\n"); GENERAL::SaveStruct($Cdef, "$$.ttt.save"); printf("see info in $$.ttt.save\n"); die; }
      $logI += log($f);
      printf("f = %e, log(f) = %f, sig = %e, mi = %f, ma = %f, x0 = %f, refE = %f\n", $f, log($f), $s, $Cdef->[$ci]->{min}, $Cdef->[$ci]->{max}, $Cm[$ci], $minE) if (defined($par{dbg}) && ($par{dbg} > 1));
    }
    print "I = $I, log(I) = $logI\n";

    # Compute the numerical version of the same integral to within the factor of the partition function,
    # using importance-sampling MC with the importance function being the boltzman distribution.
    for (my $ci = 0; $ci < $Nco; $ci++) {
      next if (defined($Cdef->[$ci]->{integrate}) && ($Cdef->[$ci]->{integrate} == 0));
      my $s = $Cstd[$ci];
      for (my $i = 0; $i < $Nsn; $i++) {
        if (uc($Cdef->[$ci]->{type}) =~ /^(DIHE|IMPR|DIST|AXISDIST|PLANEDIST|INPLACEDIST|ANGL)$/) {
          $f[$i] *= exp(-(($C->[$ci]->[$i] - $Cm[$ci])**2)/(2*$s**2));
          $logf[$i] += -(($C->[$ci]->[$i] - $Cm[$ci])**2)/(2*$s**2);
        } else {
          GENERAL::error("Unrecognized coordinate type '$Cdef->[$ci]->{type}'");
        }
        printf("%05d %05d ($Cdef->[$ci]->{type}) %e %f %f %f\n", $i, $ci, $f[$i], $logf[$i], $ener[$i], $C->[$ci]->[$i]) if (defined($par{dbg}) && ($par{dbg} > 1));
      }
    }
  } elsif ($par{refDensityType}  eq 'GaussianKernelDensityEstimate') {
    my (@type, @Cmi, @Cma);
    for (my $ci = 0; $ci < $Nco; $ci++) {
      if (defined($Cdef->[$ci]->{integrate}) && ($Cdef->[$ci]->{integrate} == 0)) { push(@type, -1); }
      elsif (uc($Cdef->[$ci]->{type}) =~ /^(DIHE|IMPR)$/) { push(@type, 2); }
      elsif (uc($Cdef->[$ci]->{type}) =~ /^(DIST|AXISDIST|PLANEDIST|INPLACEDIST|ANGL)$/) { push(@type, 1) }
      else { GENERAL::error("Unrecognized coordinate type '$Cdef->[$ci]->{type}'"); }
      if (defined($Cdef->[$ci]->{integrate}) && ($Cdef->[$ci]->{integrate} == 0)) { push(@Cmi, 0); push(@Cma, 0); }
      else { push(@Cmi, $Cdef->[$ci]->{min}); push(@Cma, $Cdef->[$ci]->{max}); }
    }
    $logI = GaussianKernelDensity($C, \@type, \@Cmi, \@Cma, \@logf);
    $I = exp($logI);
    for (my $sni = 0; $sni < scalar(@logf); $sni++) { $f[$sni] = exp($logf[$sni]); }
  } else {
    GENERAL::error("Unknown reference PDF type '$par{refDensityType}'");
  }

  my $Eh = 0; # expectation in the MC-based approximation of the integral of F computed above
  my $H = 0; # local enthalpy (to be divided by the umbrella correction factor)
  my $logEh = undef; # log of this expectation
  my $uc = 0; my @uc; # the umbrella correction factor (the expectation of exp[beta*U], where U is the biasing umbrella potential)
  for (my $i = 0; $i < $n; $i++) {
    # calculate the Jacobian (for transforming from integrating coordinates to Cartesian) of this snapshot
    my $J = 1; my $logJ = 0;
    for (my $ci = 0; $ci < $Nco; $ci++) {
      next if (defined($Cdef->[$ci]->{integrate}) && ($Cdef->[$ci]->{integrate} == 0));
      next if (defined($Cdef->[$ci]->{jacobian}) && ($Cdef->[$ci]->{jacobian} == 0));
      if (uc($Cdef->[$ci]->{type}) =~ /^(DIST)$/) {
        $J = $J * (($C->[$ci]->[$i])**2);
        $logJ = $logJ + log(($C->[$ci]->[$i])**2);
      } elsif (uc($Cdef->[$ci]->{type}) =~ /^(ANGL)$/) {
        $J = $J * sin(($PI) * $C->[$ci]->[$i] / 180);
        $logJ = $logJ + log(sin(($PI) * $C->[$ci]->[$i] / 180));
      } else {
        GENERAL::error("Unrecognized coordinate type '$Cdef->[$ci]->{type}' for Jacobian");
      }
    }
    $Eh += exp($logf[$i] + $ener[$i]/$R/$par{T} - $logJ)/$n;
    if (defined($logEh)) {
      $logEh = logSumExp($logEh, $logf[$i] + $ener[$i]/$R/$par{T} - $logJ);
    } else { $logEh = $logf[$i] + $ener[$i]/$R/$par{T} - $logJ; }
    $uc += exp($R*$par{T}*$rener[$i])/$n;
    push(@uc, exp($R*$par{T}*$rener[$i]));
    $H += ($ener[$i] - $rener[$i]) * exp($R*$par{T}*$rener[$i])/$n;
    printf("log(f) = %f; E/RT = %f, Er/RT = %f, logEh = %f, Eh = %f, logJ = %f, J = %f\n", $logf[$i], $ener[$i]/$R/$par{T}, $rener[$i]/$R/$par{T}, $logEh, $Eh, $logJ, $J) if (defined($par{dbg}) && ($par{dbg} > 2));
  }
  $logEh -= log($n);
  $H = $minE + $H/$uc; # because mean({H[r] + Hoff} * exp(U[r]))/mean(exp(U[r])) = mean(H[r] * exp(U[r]))/mean(exp(U[r])) + Hoff

  # Deduce the partition function of interest and compute its variance
  my $varEh = 0; # compute the variance of the estimated expectation $Eh
  my $varEhlog = -inf; # log of this variance
  for (my $i = 0; $i < $n; $i++) {
    $varEh += (($f[$i]/exp(-$ener[$i]/$R/$par{T}) - $Eh)**2)/$n/($n-1);
    $varEhlog = logSumExp($varEhlog, 2*($logEh + log(abs(exp($logf[$i] + $ener[$i]/$R/$par{T} - $logEh) - 1))) - log($n) - log($n-1));
#    $varEhlog = logSumExp($varEhlog, 2*log(abs($f[$i]/exp(-$ener[$i]/$R/$par{T}) - $Eh)) - log($n) - log($n-1));
  }
  my $logQ = $logI - $logEh; # log of the estimate of the configuration integral
  my $Q = exp($logQ); # estimate of the configuration integral
  my $varQlog = 2*($logI  - 2*$logEh) + $varEhlog; # log of the variance of this estimate
  if ($uc != 1) {
    if (GENERAL::stdev(\@uc) != 0) { # sometimes, there is only one snapshot with a biasing potential, so don't worry about such cases
      $Q = $Q*$uc; $logQ = $logQ + log($uc);
      print "---- Old varQlog = $varQlog...\n" if ($par{dbg});
      $varQlog = logSumExp($varQlog + 2*log($uc), 2*$logQ + 2*log(GENERAL::stdev(\@uc))); # varQ_new = varQ_old*uc^2 + Q_old^2 * var(uc)
      printf("---- New varQlog = $varQlog (uc = $uc +/- %f), logQ = %e...\n", GENERAL::stdev(\@uc), $logQ) if ($par{dbg});
    }
  }
  my $Gerr = 10**2; # some large value
  if ($varQlog/2 < $logQ) { $Gerr = $R*$par{T}*(logSumExp($logQ, $varQlog/2) - ($logQ + log(1 - exp($varQlog/2 - $logQ)))); } # error in G

  # correct for factoring out of the smallest energy
  $Q *= exp(-$minE/($R*$par{T}));
  $logQ -= $minE/($R*$par{T});
  $varQlog -= 2*$minE/($R*$par{T});

  printf("Configuration integral for this minimum Q = e^{%e} +/- e^{%e} (-RTlog(Q) = %f +/- %f); H = %f; umb. corr. = %f +/- %f\n", $logQ, $varQlog/2, -$R*$par{T}*$logQ, $Gerr, $H, $uc, GENERAL::stdev(\@uc));

  printf("At the end of integrateIIS: " . `grep VmPeak /proc/$$/status`) if (defined($par{dbg}) && ($par{dbg} > 1));
  my %ret;
  $ret{G} = -$R*$par{T}*$logQ;
  $ret{Gerr} = $Gerr;
  $ret{Q} = $Q;
  $ret{logQ} = $logQ;
  $ret{varQlog} = $varQlog;
  $ret{H} = $H;
  return \%ret;
}


# ------------------------------ Helper functions --------------------------------------

=head2

 Title   :  findOptimalSigma
 Usage   :  findOptimalSigma($Cdef->[$ci], $C->[$ci]->[$li], $C->[$ci])
 Function:  Employs Newton's method to find the standard deviation for the given degree of freedom
            that results in the lowest cross entropy.
 Returns :  Optimal value for the standard deviation of the given degree of freedom that minimizes cross entropy
 Args    :  1. coordinate definition
            2. mean value around which integration will take place
            3. trajectory vector
            

=cut

sub findOptimalSigma {
  my $cdef = shift;
  my $mu = shift;
  my $C = shift;
  my $etol = 10**(-10); # error tolerance for root finding
  my $case = 1;
  $case = 2 if ((uc($cdef->{type}) =~ /^(DIHE|IMPR)$/) && isCircularGaussianSplit($cdef, $mu));
  my $P = 180;   # pi in degrees
  my $PI = M_PI; # pi in radians

  # first, compute the sample standard deviation that will be needed many times
  my $sum = 0; my $sums = 0;
  for (my $i = 0; $i < scalar(@$C); $i++) {
    $sum += $C->[$i] - $mu;
    $sums += ($C->[$i] - $mu)**2;
  }
  my $var = $sums/scalar(@$C) - ($sum/scalar(@$C))**2;

  # next, apply Newton's method
  my $s = sqrt($var); # a good starting point is the observed standard deviation
  if ($s eq 0) {
    GENERAL::error("coordinate $cdef->{type} [" . join(", ",@{ $cdef->{atoms}}) . "] seems to have remained constant in the simulation!");
  }
  for (my $i = 1; $i <= 1000; $i++) {
    my ($f, $fp) = crossEntropyDerivatives($var, $cdef->{min}, $cdef->{max}, $mu, $s, $case);
    return $s if (abs($f) < $etol);   # found root
    last if (abs($fp) < $etol);          # the derivative is too shallow - value will not improve much, so check limit
    $s -= $f/$fp;
  }

  # looks like Newton's method is not converging, which may be because there is no root,
  # so compare current cross entropy with that at infinite sigma
  my $n = scalar(@$C);
  my ($cce, $ice); # current value of cross entropy and its limit as sigma approaches infinity
  if ($case == 2) {
    $cce = $var/(2*$s**2) + log(sqrt($PI/2)) + log($s) + log(erf(($cdef->{max} - $mu)/(sqrt(2)*$s)) + erf(($mu - $cdef->{min})/(sqrt(2)*$s)) + 2*erf($P/(sqrt(2)*$s)));
    $ice = log(sqrt($PI/2)) + log(sqrt(2/$PI)*($cdef->{max} - $cdef->{min} + 2*$P));
  } else {
    $cce = $var/(2*$s**2) + log(sqrt($PI/2)) + log($s) + log(erf(($cdef->{max} - $mu)/(sqrt(2)*$s)) + erf(($mu - $cdef->{min})/(sqrt(2)*$s)));
    $ice = log(sqrt($PI/2)) + log(sqrt(2/$PI)*($cdef->{max} - $cdef->{min}));
  }
  if ($ice < $cce) {
    $s = 10**10;
  }
  return $s;
}

=head2

 Title   :  crossEntropyDerivatives
 Usage   :  crossEntropyDerivatives($var, $cdef->{min}, $cdef->{max}, $mu, $s, $case)
 Function:  Computes the first and second partial derivatives of cross entropy with respect to the reference-PDF 
            standard deviation in the direction of the given degree of freedom. For simplicity, actually computes
            the first derivative of cross entropy (H) times sigma^2 (e.g. (dH/ds)*s^2) and the first derivative of
            that function (e.g. d((dH/ds)*s^2))/ds). This gives fewer terms and is easier to compute. This does not
            loose any roots, because sigma = inf is still a "root" for the modified function (e.g. the modified first and
            second derivatives are also zero at sigma -> inf).
 Returns :  First and second derivatives.
 Args    :  1. coordinate definition
            2. mean value around which integration will take place
            3. trajectory vector
            

=cut

sub crossEntropyDerivatives {
  my $S = shift;
  my $xmin = shift;
  my $xmax = shift;
  my $xm = shift;
  my $s = shift;
  my $case = shift;
  my $PI = M_PI;

  my $A = $xmax - $xm;
  my $B = $xm - $xmin;
  my $a = $A/$s; my $a2 = $a**2; my $a3 = $a**3;
  my $b = $B/$s; my $b2 = $b**2; my $b3 = $b**3;
  my $P = 180;
  my $c = $P/$s; my $c2 = $c**2; my $c3 = $c**3;
  my $ts2 = 2 * $s**2;
  my $sq2 = sqrt(2);
  my ($d1, $d2);
  if ($case == 2) {
    $d1 = -$S/$s + $s - $sq2*($A*exp(-$A**2/$ts2) + $B*exp(-$B**2/$ts2) + 2*$P*exp(-$P**2/$ts2))/(sqrt($PI) * (erf($A/sqrt($ts2)) + erf($B/sqrt($ts2)) + 2*erf($P/sqrt($ts2))));
    $d2 = $S/$s**2 + 1 - ($a3*exp(-$a2/2) + $b3*exp(-$b2/2) + 2*$c3*exp(-$c2/2))/(sqrt($PI/2)* (erf($a/$sq2) + erf($b/$sq2) + 2*erf($c/$sq2))) - ($a*exp(-$a2/2) + $b*exp(-$b2/2) + 2*$c*exp(-$c2/2))**2/(($PI/2)*(erf($a/$sq2) + erf($b/$sq2) + 2*erf($c/$sq2))**2);
  } else {
    $d1 = -$S/$s + $s - $sq2*($A*exp(-$A**2/$ts2) + $B*exp(-$B**2/$ts2))/(sqrt($PI) * (erf($A/sqrt($ts2)) + erf($B/sqrt($ts2))));
    $d2 = $S/$s**2 + 1 - ($a3*exp(-$a2/2) + $b3*exp(-$b2/2))/(sqrt($PI/2)* (erf($a/$sq2) + erf($b/$sq2))) - ($a*exp(-$a2/2) + $b*exp(-$b2/2))**2/(($PI/2)*(erf($a/$sq2) + erf($b/$sq2))**2);
  }
  return ($d1, $d2);
}

=head2

 Title   :  isCircularGaussianSplit
 Usage   :  isCircularGaussianSplit($Cdef->[$ci], $C->[$ci]->[$li])
 Function:  Given the definition of a dihedral degree of freedom, and a mean integration value,
            determines whether when the Gaussian defined on a circle, around the mean value, is
            unwrapped onto a line, the overall integral becomes split into integrals over two segments.
 Notes   :  The function we are integrating is a Gaussian of angles, with the difference defined via angleDiff(). By defining the difference to be the angular
            difference (which makes sense from the perspective of the potential), we essentially cause the number axis to wrap around. If we measure differences
            with some reference angle a0, this wrapping around happens at a0+pi. So, if our integration ranges are defined in such a way that when going from 
            the min to the max angle we span a0+pi, then the integral is really split into two fragments that have to be computed separately.
 Returns :
 Args    :  1. definition of a dihedral degree of freedom
            2. mean integration angle in degrees

=cut

sub isCircularGaussianSplit {
  my $cdef = shift;
  my $mu = shift;
  if (defined($cdef->{fullCircle}) || (VALOCIDY::angleDiff($cdef->{min}, $cdef->{max}) == 0) || (VALOCIDY::angleDiffCCW($mu + 180, $cdef->{min}) < VALOCIDY::angleDiffCCW($cdef->{max}, $cdef->{min}))) {
    return 1;
  }
  return 0;
}

=head2

 Title   :  prepareForBuilding
 Usage   :  prepareForBuilding($pdb, $psf, CHARMM::readCHARMMParameters($opts{toppar}->{par}, 1), $build)
 Function:  Prepares the build structure returned by NAMD::generateBAT for building/sampling of the given
            PDB structure. Fetches atoms and force-field parameters corresponding to the build instructions.
 Returns :  Nothing
 Args    :  1. PDB object
            2. PSF object
            3. CHARMM parameter hash
            4. build structure

=cut

sub prepareForBuilding {
  my $pdb = shift;
  my $psf = shift;
  my $B = shift;
  my $A = shift;
  my $T = shift;
  my $placement = shift;
  my $P = shift;
  my $PI = M_PI;
  my $af = (1/180) * $PI;

  # create some hashes in PSF for convenience
  $psf->{bondH} = (); $psf->{anglH} = (); $psf->{diheH} = ();
  foreach my $p (@{$psf->{bonds}}) {
    $psf->{bondH}{$p->[0]}{$p->[1]} = 1;
    $psf->{bondH}{$p->[1]}{$p->[0]} = 1;
  }
  foreach my $p (@{$psf->{angles}}) {
    $psf->{anglH}{$p->[0]}{$p->[1]}{$p->[2]} = 1;
    $psf->{anglH}{$p->[2]}{$p->[1]}{$p->[0]} = 1;
  }
  foreach my $p (@{$psf->{dihedrals}}) {
    $psf->{diheH}{$p->[0]}{$p->[1]}{$p->[2]}{$p->[3]} = 1;
    $psf->{diheH}{$p->[3]}{$p->[2]}{$p->[1]}{$p->[0]} = 1;
  }
  foreach my $p (@{$psf->{impropers}}) {
    $psf->{imprH}{$p->[0]}{$p->[1]}{$p->[2]}{$p->[3]} = 1;
    $psf->{imprH}{$p->[3]}{$p->[2]}{$p->[1]}{$p->[0]} = 1;
  }

  # find atoms to be placed in the PDB object, annotate BAT with corresponding force-field parameters, where possible
  my (@place, @pB, @pA, @pT);
  for (my $i = 0; $i < scalar(@$placement); $i++) {
    my @placeOp;
    foreach my $ai (@{$placement->[$i]}) {
      $psf->{atoms}->[$ai]->{pdbatom} = PDB::getAtomInRes($pdb->getResByInd($psf->{atoms}->[$ai]->{segname}, $psf->{atoms}->[$ai]->{resnum}), $psf->{atoms}->[$ai]->{name});
      push(@placeOp, $psf->{atoms}->[$ai]->{pdbatom});
    }
    push(@place, \@placeOp);
  }
  foreach my $b (@$B) {
    push(@pB, defined($psf->{bondH}{$b->[0]}{$b->[1]}) ? CHARMM::getBondPar($P, $psf->{atoms}->[$b->[0]]->{type}, $psf->{atoms}->[$b->[1]]->{type}, -1) : undef);
  }
  foreach my $a (@$A) {
    push(@pA, defined($psf->{anglH}{$a->[0]}{$a->[1]}{$a->[2]}) ? CHARMM::getAnglPar($P, $psf->{atoms}->[$a->[0]]->{type}, $psf->{atoms}->[$a->[1]]->{type}, $psf->{atoms}->[$a->[2]]->{type}, -1) : undef);
  }
  foreach my $t (@$T) {
    my @at = ($psf->{atoms}->[$t->[0]]->{type}, $psf->{atoms}->[$t->[1]]->{type}, $psf->{atoms}->[$t->[2]]->{type}, $psf->{atoms}->[$t->[3]]->{type});
    if (defined($psf->{diheH}{$t->[0]}{$t->[1]}{$t->[2]}{$t->[3]})) {
      push(@pT, CHARMM::getDihePar($P, $at[0], $at[1], $at[2], $at[3], -1));
    } elsif (defined($psf->{imprH}{$t->[0]}{$t->[1]}{$t->[2]}{$t->[3]})) {
      push(@pT, CHARMM::getImprPar($P, $at[0], $at[1], $at[2], $at[3], -1));
    } else {
      # is there a different improper defined with the same four atoms, which effectively restricts this improper?
      my $suc = 0;
      my %D = ($t->[0], 0, $t->[1], 1, $t->[2], 2, $t->[3], 3);
      my %A = ($t->[0], 0, $t->[1], 1, $t->[2], 2); my %B = ($t->[1], 1, $t->[2], 2, $t->[3], 3); # the two planes defining our dihedral
      foreach my $p (@{$psf->{impropers}}) {
        # --- condition O: the defined improper has the same atoms as the desired dihedral
        next if (!(defined($D{$p->[0]}) && defined($D{$p->[1]}) && defined($D{$p->[2]}) && defined($D{$p->[3]})));

        # --- condition 1: the defined improper shares one plane with the needed improper
        my @sharedPlane;
        if (defined($A{$p->[0]}) && defined($A{$p->[1]}) && defined($A{$p->[2]}) || defined($B{$p->[0]}) && defined($B{$p->[1]}) && defined($B{$p->[2]})) {
          @sharedPlane = ($p->[0], $p->[1], $p->[2]);
        } elsif (defined($A{$p->[1]}) && defined($A{$p->[2]}) && defined($A{$p->[3]}) || defined($B{$p->[1]}) && defined($B{$p->[2]}) && defined($B{$p->[3]})) {
          @sharedPlane = ($p->[1], $p->[2], $p->[3]);
        } else { next; }
        # annotate all atoms: the defined dihedral used will be N-S-O-O1, and the desired dihedral will be S-O-N-O1
        my @pc = ($p->[1], $p->[2]); my @tc = ($t->[1], $t->[2]); # central atoms in the two impropers
        my @O = GENERAL::intersect(\@pc, \@tc); GENERAL::assert(scalar(@O) == 1, "two impropers sharing only one plane are expected to have exactly one atom in common between the two central atoms!");
        $O = $O[0]; # atom that is shared etween all three planes - the shared plane plus the two different planes, one for each improper
        my @O1 = GENERAL::setdiff($p, \@sharedPlane); GENERAL::assert(scalar(@O) == 1, "all improper atoms minus three shared between two impropers should equal to one atom!");
        $O1 = $O1[0]; # O-O1 is the side shared by the two planes unique to each dihedral
        $S = ($p->[1] == $O) ? $p->[2] : $p->[1]; # atom shared by the two planes of the defined dihedral, unique to the defined dihedral
        $N = ($t->[1] == $O) ? $t->[2] : $t->[1]; # atom shared by the two planes of the desired dihedral, unique to the desired dihedral

        # condition 2: all of the necessary (pseudo) bonds and angles are restricted from the PSF/force-field
        my ($par, $sign); # parameter corresponding to defined dihedral, and sign given our specific interpretation of order
        if ((($p->[0] == $N) && ($p->[1] == $S) && ($p->[2] == $O) && ($p->[3] == $O1)) || (($p->[3] == $N) && ($p->[2] == $S) && ($p->[1] == $O) && ($p->[0] == $O1))) {
          $sign = 1;
          $par = CHARMM::getImprPar($P, $psf->{atoms}->[$p->[0]]->{type}, $psf->{atoms}->[$p->[1]]->{type}, $psf->{atoms}->[$p->[2]]->{type}, $psf->{atoms}->[$p->[3]]->{type}, -1);
        } elsif ((($p->[0] == $N) && ($p->[1] == $O) && ($p->[2] == $S) && ($p->[3] == $O1)) || (($p->[3] == $N) && ($p->[2] == $O) && ($p->[1] == $S) && ($p->[0] == $O1))) {
          $sign = -1;
          $par = CHARMM::getImprPar($P, $psf->{atoms}->[$p->[0]]->{type}, $psf->{atoms}->[$p->[2]]->{type}, $psf->{atoms}->[$p->[1]]->{type}, $psf->{atoms}->[$p->[3]]->{type}, -1);
        } else {
          GENERAL::error("Something strange happened - after annotating the atoms of defined improper $p->[0] $p->[1] $p->[2] $p->[3], I can't find it anymore!");
        }
        if (!defined($par)) { GENERAL::warning("improper $t->[0] $t->[1] $t->[2] $t->[3] in PSF, but no parameter for it found!"); next; }
        my $phi1 = $af*$sign*$par->{psi0}; # defined improper equilibrium value, given our prefered of central atoms
        my $al1 = VALOCIDY::getPseudoAngleEquil($P, $psf, $S, $O, $O1); next if (!defined($al1));
        my $al2 = VALOCIDY::getPseudoAngleEquil($P, $psf, $N, $O, $O1); next if (!defined($al2));
        my $beta = VALOCIDY::getPseudoAngleEquil($P, $psf, $S, $O, $N); next if (!defined($beta));

        # compute the equilibrium value and force constant scale factor of desired dihedral
        my $phi2 = asin(sin($phi1) * sin($al1)/sin($al2));
        my $ksf = abs(cos($phi1)*(sin($al1)/sin($al2))/(sqrt(1 - (sin($phi1) * sin($al1)/sin($al2))**2)));
        # check if it is the complimentary angle that we need
        if (cos($al1)*tan($beta) < sin($al1)*cos($phi1)) {
          $phi2 = $PI - $phi2;
        }

        # create a new improper parameter
        my %npar;
        $npar{Kpsi} = $par->{Kpsi} * $ksf;
        $npar{n} = $par->{n};
        if ((($t->[0] == $S) && ($t->[1] == $O) && ($t->[2] == $N) && ($t->[3] == $O1)) || (($t->[3] == $S) && ($t->[2] == $O) && ($t->[1] == $N) && ($t->[0] == $O1))) {
          $npar{psi0} = $phi2/$af;
        } elsif ((($t->[0] == $S) && ($t->[1] == $N) && ($t->[2] == $O) && ($t->[3] == $O1)) || (($t->[3] == $S) && ($t->[2] == $N) && ($t->[1] == $O) && ($t->[0] == $O1))) {
          $npar{psi0} = -$phi2/$af;
        } else {
          GENERAL::error("Something strange happened - after annotating the atoms of desired improper $t->[0] $t->[1] $t->[2] $t->[3], I can't find it anymore!");
        }
        push(@pT, \%npar);
        $suc = 1;
        last;
      }
      if (!$suc) {
        push(@pT, undef);
        GENERAL::warning("could not derive equilibrium value/force constant for torsion $psf->{atoms}->[$t->[0]]->{tclSel} - $psf->{atoms}->[$t->[1]]->{tclSel} - $psf->{atoms}->[$t->[2]]->{tclSel} - $psf->{atoms}->[$t->[3]]->{tclSel}");
      }
    }
  }
  return (\@place, \@pB, \@pA, \@pT);
}


=head2

 Title   :  getPseudoAngleEquil
 Usage   :  getPseudoAngleEquil($P, $psf, $A, $B, $C)
 Function:  Given atoms $A-$B-$C, if the corresponding angle is defined in PSF/force-field parameters, then returns its
            equilibrium value, and otherwise tries to compute the equilibrium value from other bonds/angles defined in PSF/pars.
 Returns :  Equilibrium value or undef if fails
 Args    :  1. CHARMM parameters object
            2. PSF object
            3. Atom 1 index
            4. Atom 2 index
            5. Atom 3 index

=cut

sub getPseudoAngleEquil {
  my $P = shift;
  my $psf = shift;
  my $A = shift;
  my $B = shift;
  my $C = shift;
  my $eval = undef;

  # is the angle defined in PSF/parames?
  if (defined($psf->{angleH}{$A}{$B}{$C})) {
    my $par = CHARMM::getAnglPar($P, $psf->{atoms}->[$A]->{type}, $psf->{atoms}->[$B]->{type}, $psf->{atoms}->[$C]->{type}, -1);
    $eval = $par->{Theta0} if (defined($par));
  }
  return $eval if (defined($eval));

  # if not, try to calculate its equilibrium value
  my $AB = getPseudoBondEquil($P, $psf, $A, $B);
  my $BC = getPseudoBondEquil($P, $psf, $B, $C);
  my $AC = getPseudoBondEquil($P, $psf, $A, $C);
  return undef if (!defined($AB) || !defined($BC) || !defined($AC));
  my $cos = ($AB**2 + $BC**2 - $AC**2)/(2*$AB*$BC);
  return acos($cos);
}

=head2

 Title   :  getPseudoBondEquil
 Usage   :  getPseudoBondEquil($P, $psf, $A, $B)
 Function:  Given atoms $A-$B, if the corresponding bond is defined in PSF/force-field parameters, then returns its
            equilibrium value, and otherwise tries to compute the equilibrium value from other bonds/angles defined in PSF/pars.
 Returns :  Equilibrium value or undef if fails
 Args    :  1. CHARMM parameters object
            2. PSF object
            3. Atom 1 index
            4. Atom 2 index

=cut

sub getPseudoBondEquil {
  my $P = shift;
  my $psf = shift;
  my $A = shift;
  my $B = shift;
  my $eval = undef;
  my $PI = M_PI;
  my $af = (1/180) * $PI;

  # is the angle defined in PSF/parames?
  if (defined($psf->{bondH}{$A}{$B})) {
    my $par = CHARMM::getBondPar($P, $psf->{atoms}->[$A]->{type}, $psf->{atoms}->[$B]->{type}, -1);
    $eval = $par->{b0} if (defined($par));
  }
  return $eval if (defined($eval));

  # if not, try to calculate its equilibrium value - try to find some atom X, to which both atoms are bonded
  foreach my $angle (@{$psf->{angles}}) {
    next if (!((($angle->[0] == $A) && ($angle->[2] == $B)) || (($angle->[0] == $B) && ($angle->[2] == $A))));
    my $X = $angle->[1];
    my $AX = CHARMM::getBondPar($P, $psf->{atoms}->[$A]->{type}, $psf->{atoms}->[$X]->{type}, -1);
    my $BX = CHARMM::getBondPar($P, $psf->{atoms}->[$B]->{type}, $psf->{atoms}->[$X]->{type}, -1);
    my $alpha = CHARMM::getAnglPar($P, $psf->{atoms}->[$A]->{type}, $psf->{atoms}->[$X]->{type}, $psf->{atoms}->[$B]->{type}, -1);
    next if (!defined($AX) || !defined($BX) || !defined($alpha));
    return sqrt($AX->{b0}**2 + $BX->{b0}**2 - 2*$AX->{b0}*$BX->{b0}*cos($af*$alpha->{Theta0}));
  }
  return undef;
}

=head2

 Title   :  buildCoordinates
 Usage   :  buildCoordinates($pdb, $build)
 Function:  Rebuilds the coordinates of the given PDB from the build instructions and corresponding coordinate value specified.
 Returns :  False if for some reason the build fails (i.e. one of the important angles is 180), true otherwise
 Args    :  1. PDB object
            2. Build structure definition array
            

=cut

sub buildCoordinates {
  my $pdb = shift;
  my $cB = shift;
  my $cA = shift;
  my $cT = shift;
  my $place = shift;
  my $PI = M_PI;
  my $af = (1/180) * $PI;

  my $bi = 0; $ai = 0; $ti = 0;
  for (my $i = 0; $i < scalar(@$place); $i++) {
    my $a = $place->[$i]->[0];
    if ($i > 2) {
      if (!PDB::buildAtom($a, $place->[$i]->[1], $place->[$i]->[2], $place->[$i]->[3], $cB->[$bi], $cA->[$ai], $cT->[$ti])) { return 0; }
      $bi++; $ai++; $ti++;
    } elsif ($i == 0) {
      $a->{xcoor} = 0.0; $a->{ycoor} = 0.0; $a->{zcoor} = 0.0;
    } elsif ($i == 1) {
      $a->{xcoor} = $cB->[$bi];
      $a->{ycoor} = 0.0;
      $a->{zcoor} = 0.0;
      $bi++;
    } elsif ($i == 2) {
      $a->{xcoor} = $place->[$i]->[1]->{xcoor} + $cB->[$bi]*cos($PI - $af*$cA->[$ai]);
      $a->{ycoor} = $cB->[$bi]*sin($PI - $af*$cA->[$ai]);
      $a->{zcoor} = 0.0;
      $bi++; $ai++;
    }
  }
  return 1;
}

=head2

 Title   :  setupCoordinateRestraints
 Usage   :  setupCoordinateRestraints($chrm, );
 Function:  Given a CHARMM structure, setups up coordinate restraints for those coordinates that need to be limited.
 Returns :
 Args    :  1. CHARMM object
            2. 

=cut

sub setupCoordinateRestraints {
  my $chrm = shift;
  my $Cdef = shift;
  my $opts = shift;
  GENERAL::requireArgs($chrm, $Cdef, $opts);

  my $FC = {}; $FC = $opts->{FC} if (defined($opts->{FC}));
  $FC->{dihe} = 1000 if (!defined($FC->{dihe}));
  $FC->{dist} = 100 if (!defined($FC->{dist}));
  $FC->{inplacedist} = 100 if (!defined($FC->{inplacedist}));
  $FC->{axisdist} = 100 if (!defined($FC->{axisdist}));
  $FC->{planedist} = 100 if (!defined($FC->{planedist}));
  my $dA = 3; $dA = $opts->{mA} if (defined($opts->{mA})); # margin (in degrees) into the desired area where the restraining potential picks up - SHOULD BE REPLACED WITH COORDINATE-SPECIFI MARGINS!

  $chrm->send("skipe exclude GEO\n"); # append GEO to the list of calculated terms (in case a previous skipe was run)
  $chrm->send("MMFP\n");
  my $Nr = 0; # number of restraints
  for (my $ci = 0; $ci < scalar(@$Cdef); $ci++) {
    next if ((defined($Cdef->[$ci]->{determineRange}) && $Cdef->[$ci]->{determineRange}) || (!defined($Cdef->[$ci]->{min})) || (!defined($Cdef->[$ci]->{max})) || (defined($Cdef->[$ci]->{doNotConstrain}) && $Cdef->[$ci]->{doNotConstrain}));
    if (uc($Cdef->[$ci]->{type}) eq "DIHE") {
      $Nr += 8;
    } elsif (uc($Cdef->[$ci]->{type}) eq "DIST") {
      $Nr += 4;
    } elsif (uc($Cdef->[$ci]->{type}) eq "INPLACEDIST") {
      $Nr += 2;
    } elsif (uc($Cdef->[$ci]->{type}) eq "AXISDIST") {
      $Nr += 2;
    } elsif (uc($Cdef->[$ci]->{type}) eq "PLANEDIST") {
      $Nr += 2;
    } else {
      GENERAL::error("Unknown phase-space coordinate type '$Cdef->[$ci]->{type}'");
    }
  }
  my $c = 0;
  for (my $ci = 0; $ci < scalar(@$Cdef); $ci++) {
    # NOTE: keyword 'INSIDE' means keep the vlaue below the reference value (no potential when value is below reference value)
    #       keyword 'OUTSIDE' means keep the vlaue above the reference value (no potential when value is above reference value)
    next if ((defined($Cdef->[$ci]->{determineRange}) && $Cdef->[$ci]->{determineRange}) || (!defined($Cdef->[$ci]->{min})) || (!defined($Cdef->[$ci]->{max})) || (defined($Cdef->[$ci]->{doNotConstrain}) && $Cdef->[$ci]->{doNotConstrain}));
    if (uc($Cdef->[$ci]->{type}) eq "DIHE") {
      # the following makes sure that the constraint is smooth upon the dihedral changing sign (between +/-180)
      my $off = -($Cdef->[$ci]->{min} + $Cdef->[$ci]->{max})/2; # the amount by which to move both lower and upper bounds to make sure that the center is at zero
      my $mi = $Cdef->[$ci]->{min} + $off; my $ma = $Cdef->[$ci]->{max} + $off; # the lower and upper bounds of the moved constraints (these are the new constrains if the offset is applied to the coordinate)
      if ($ma < $mi) {
        $off = 180 + $off; $off = $off - 360 if ($off > 180);
        $ma = 180 + $ma; $mi = $mi - 180; # keep in mind that here ma is -mi (they are symmetric about +/- PI), so doing this will result in both beeing between -PI and PI and symmetric about zero
      }
      my ($off_in, $off_out); # CHARMM always wants a positive offset, which is subtracted for INSIDE and added for OUTSIDE
      if ($off < 0) {
        $off_in = abs($off); $off_out = $off + 360;
      } else {
        $off_in = 360 - $off; $off_out = $off;
      }
      $ma -= $dA if ($ma > $dA); $mi += $dA if ($mi < $dA); # make the potential start slightly before the edge of the desired area of phase space
      $Cdef->[$ci]->{pmin} = $mi; $Cdef->[$ci]->{pmax} = $ma; # min/max where potential was applied
      printf("[$Cdef->[$ci]->{min} $Cdef->[$ci]->{max}] will offset by $off and apply potential below $ma and above $mi...\n")  if (defined($opts->{dbg}) && $opts->{dbg});
      $chrm->send("GEO " . ($c == 0 ? sprintf("MAXGEO %d", $Nr) : "") . " SPHERE DIHEDRAL TREF $mi DTOFF $off_out HARMONIC OUTSIDE FORCE $FC->{dihe} P1 -1 -\n");
      $chrm->send("select $Cdef->[$ci]->{atoms}->[0] end select $Cdef->[$ci]->{atoms}->[1] end select $Cdef->[$ci]->{atoms}->[2] end select $Cdef->[$ci]->{atoms}->[3] end\n");
      $chrm->send("GEO SPHERE DIHEDRAL TREF $ma DTOFF $off_in HARMONIC INSIDE FORCE $FC->{dihe} P1 -1 -\n");
      $chrm->send("select $Cdef->[$ci]->{atoms}->[0] end select $Cdef->[$ci]->{atoms}->[1] end select $Cdef->[$ci]->{atoms}->[2] end select $Cdef->[$ci]->{atoms}->[3] end\n");
    } elsif (uc($Cdef->[$ci]->{type}) eq "DIST") {
      my $mi = $Cdef->[$ci]->{min} + $Cdef->[$ci]->{margin}; my $ma = $Cdef->[$ci]->{max} - $Cdef->[$ci]->{margin};
      $Cdef->[$ci]->{pmin} = $mi; $Cdef->[$ci]->{pmax} = $ma; # min/max where potential was applied
      $chrm->send("GEO " . ($c == 0 ? sprintf("MAXGEO %d", $Nr) : "") . " SPHERE DISTANCE DROFF $mi HARMONIC OUTSIDE FORCE $FC->{dist} -\n");
      $chrm->send("select $Cdef->[$ci]->{atoms}->[0] end select $Cdef->[$ci]->{atoms}->[1] end\n");
      $chrm->send("GEO SPHERE DISTANCE DROFF $ma HARMONIC INSIDE FORCE $FC->{dist} -\n");
      $chrm->send("select $Cdef->[$ci]->{atoms}->[0] end select $Cdef->[$ci]->{atoms}->[1] end\n");
    } elsif (uc($Cdef->[$ci]->{type}) eq "INPLACEDIST") {
      my $ma = $Cdef->[$ci]->{max} - $Cdef->[$ci]->{margin};
      $Cdef->[$ci]->{pmax} = $ma; # max where potential was applied
      $chrm->send("GEO " . ($c == 0 ? sprintf("MAXGEO %d", $Nr) : "") . " SPHERE DROFF $ma HARMONIC INSIDE FORCE $FC->{inplacedist} -\n");
      $chrm->send("XREF $Cdef->[$ci]->{origX} YREF $Cdef->[$ci]->{origY} ZREF $Cdef->[$ci]->{origZ} select $Cdef->[$ci]->{atoms}->[0] end\n");
    } elsif (uc($Cdef->[$ci]->{type}) eq "AXISDIST") {
      my $mi = $Cdef->[$ci]->{min} + $Cdef->[$ci]->{margin}; my $ma = $Cdef->[$ci]->{max} - $Cdef->[$ci]->{margin};
      $Cdef->[$ci]->{pmin} = $mi; $Cdef->[$ci]->{pmax} = $ma; # min/max where potential was applied
      my $xref = $Cdef->[$ci]->{axisA}->[0]; my $yref = $Cdef->[$ci]->{axisA}->[1]; my $zref = $Cdef->[$ci]->{axisA}->[2];
      my $xdir = $Cdef->[$ci]->{axisB}->[0] - $xref; my $ydir = $Cdef->[$ci]->{axisB}->[1] - $yref; my $zdir = $Cdef->[$ci]->{axisB}->[2] - $zref;
      my $L = sqrt($xdir**2 + $ydir**2 + $zdir**2); $xdir /= $L; $ydir /= $L; $zdir /= $L;
      $chrm->send("GEO " . ($c == 0 ? sprintf("MAXGEO %d", $Nr) : "") . " CYLINDER DROFF $mi HARMONIC OUTSIDE FORCE $FC->{axisdist} -\n");
      $chrm->send("XREF $xref YREF $yref ZREF $zref XDIR $xdir YDIR $ydir ZDIR $zdir select $Cdef->[$ci]->{atoms}->[0] end\n");
      $chrm->send("GEO CYLINDER DROFF $ma HARMONIC INSIDE FORCE $FC->{axisdist} -\n");
      $chrm->send("XREF $xref YREF $yref ZREF $zref XDIR $xdir YDIR $ydir ZDIR $zdir select $Cdef->[$ci]->{atoms}->[0] end\n");
    } elsif (uc($Cdef->[$ci]->{type}) eq "PLANEDIST") {
      my $mi = $Cdef->[$ci]->{min} + $Cdef->[$ci]->{margin}; my $ma = $Cdef->[$ci]->{max} - $Cdef->[$ci]->{margin};
      $Cdef->[$ci]->{pmin} = $mi; $Cdef->[$ci]->{pmax} = $ma; # min/max where potential was applied
      my $xref = $Cdef->[$ci]->{planePoint}->[0]; my $yref = $Cdef->[$ci]->{planePoint}->[1]; my $zref = $Cdef->[$ci]->{planePoint}->[2];
      my $xdir = $Cdef->[$ci]->{planeNorm}->[0]; my $ydir = $Cdef->[$ci]->{planeNorm}->[1]; my $zdir = $Cdef->[$ci]->{planeNorm}->[2];
      my $L = sqrt($xdir**2 + $ydir**2 + $zdir**2); $xdir /= $L; $ydir /= $L; $zdir /= $L;
      my $fc = defined($Cdef->[$ci]->{forceConst}) ? $Cdef->[$ci]->{forceConst} : $FC->{planedist};
      $chrm->send("GEO " . ($c == 0 ? sprintf("MAXGEO %d", $Nr) : "") . " PLANE DROFF $mi HARMONIC OUTSIDE FORCE $fc -\n");
      $chrm->send("XREF $xref YREF $yref ZREF $zref XDIR $xdir YDIR $ydir ZDIR $zdir select $Cdef->[$ci]->{atoms}->[0] end\n");
      $chrm->send("GEO PLANE DROFF $ma HARMONIC INSIDE FORCE $fc -\n");
      $chrm->send("XREF $xref YREF $yref ZREF $zref XDIR $xdir YDIR $ydir ZDIR $zdir select $Cdef->[$ci]->{atoms}->[0] end\n");
    } else {
      GENERAL::error("Unrecognized coordinate type '$Cdef->[$ci]->{type}'");
    }
    $c++;
  }
  $chrm->send("END\n");
}

=head2

 Title   :  analyzeRun
 Usage   :  my $ret = analyzeRun($chrm, $outdcd, $parMD->{dbg}, "?GEO", $Cdef);
 Function:  Given a trajectory file, returns various energies and geometric coordinates. Uses VMD to compute
            the geometric coordinates, which is faster.
 Returns :
 Args    :  1. CHARMM object
            2. Trajectory file
            3. Debug level (integer)
            4. Coordinate definition hash reference
            5. Flag that specifies whether the average RMS fluctuations should be computed

=cut

sub analyzeRun {
  my $chrm = shift;
  my $trajf = shift;
  my $dbg = shift;
  my $terms = shift; # additional terms to save from each snapshot
  $terms = "" if (!defined($terms));
  my $Cdef = shift; # definitions of phase space-defining coordinates
  if (!defined($Cdef)) { my @tmp; $Cdef = \@tmp; }
  my $avrmsf = shift; $avrmsf = 0 if (!defined($avrmsf)); # flag - if set, will return average structure and RMSF for each atom
  my $sels = shift; $sels = GENERAL::arrayRef() if (!defined($sels)); # will also print total energy for each of these selections, if specified
  my $findMinEnerCoor = shift; $findMinEnerCoor = 0 if (!defined($findMinEnerCoor)); # if specified, will copy the lowest-energy snapshot into comparison coordinates
  my $time = time();
  my $measureTime = $dbg;

  my $tag = $$ . "-" . lc(GENERAL::GetMachine());
  my $acrdf = "_avg.$tag.crd";
  my $cref = "_ef.$tag.dat";
  my $psff = "_tmp.$tag.psf";

  # first get RMSF and energies with CHARMM
  $chrm->send("OPEN UNIT 21 READ FILE NAME $trajf\n");
  if ($avrmsf) {
    $chrm->send("COOR DYNA COMP FIRSTU 21 NUNIT 1 SELE all END ORIENT SELE all END\n");
    $chrm->send("WRITE COOR COMP CARD NAME $acrdf\n* avg coords with rms fluct in weight array\n*\n\n");
    $chrm->send("REWIND UNIT 21\n");
  }
  $chrm->send("TRAJ FIRSTU 21 NUNIT 1\n");
  $chrm->send("SET K = 1\n");
  $chrm->send("SET MINENER = 9999999\n");
  $chrm->send("OPEN UNIT 11 WRITE CARD NAME $cref\n");
  # loops can not be in an interactive stream, so switch over
  my $tmpsf = "_strea.tmp.$tag.inp"; my $ofh = GENERAL::GetOutFH($tmpsf);
  $ofh->print("LABEL LOOP\n");
  $ofh->print("TRAJ READ\n");
  $ofh->print("energy\n");
  $ofh->print("WRITE TITLE UNIT 11\n");
  $ofh->print("* ?ENER $terms\n*\n\n");
  foreach my $sele (@$sels) {
    $ofh->print("inte select $sele end\n");
    $ofh->print("WRITE TITLE UNIT 11\n");
    $ofh->print("* ?ENER $terms\n*\n\n");
  }
  if ($findMinEnerCoor) {
    $ofh->print("IF ?ENER LT \@MINENER COOR COPY COMP sele all end\n");
    $ofh->print("IF ?ENER LT \@MINENER SET MINENER = ?ENER\n");
  }
  $ofh->print("INCR K BY 1\n");
  $ofh->print("IF K LE ?NFILE GOTO LOOP\n");
  close($ofh);
  $chrm->send("stream $tmpsf\n");
  $chrm->closeCard(21); $chrm->closeCard(11);
  $chrm->wait('./', 1, 300, 0); # wait until the buffer is clear (window run is complete)
  if ($measureTime) { printf("--> TIME (analyzeRun): extracting energies via CHARMM took %d seconds\n", time() - $time); $time = time(); }

  # then get coordinates with VMD (much faster)
  $chrm->writePSF($psff);
  $chrm->wait("./");
  my $C = NAMD::getCoordinatesVMD(undef, 'CDEF', $Cdef, 'TRAJ', $trajf, 'PSF', $psff);
  if ($measureTime) { printf("--> TIME (analyzeRun): extracting coordinates via VMD took %d seconds\n", time() - $time); $time = time(); }

  # read average structure and RMSF information and energies
  my $avg; $avg = PDB::new($acrdf) if ($avrmsf);
  my $ifh = GENERAL::GetInFH($cref);
  my @ener;
  while (1) {
    my $line = <$ifh>; last if (!defined($line));
    # -- total energy/other terms
    my @line = split(" ", $line);
    GENERAL::assert(scalar(@line) >= 1, "Expected at least one element in each energy line of output file $cref...");
    if (scalar(@line) == 1) {
      push(@ener, $line[0]);
    } else {
      GENERAL::assert(scalar(@{$ener[0]}) == scalar(@line), "Inconsistent number of energy terms inenergy line of output file $cref...") if (scalar(@ener) > 0);
      push(@ener, \@line);
    }
  }
  close($ifh);

  GENERAL::crm($acrdf) if ($avrmsf && !defined($dbg));
  GENERAL::crm($cref, $tmpsf, $psff) if (!defined($dbg));

  # error check
  GENERAL::assert(scalar(@ener) == scalar(@{$C->[0]}), sprintf("Read %d snapshot energies and %d snapshot coordinate sets!", scalar(@ener), scalar(@{$C->[0]}))) if (scalar(@$Cdef) > 0);

  # compose the return hash
  my %ret;
  $ret{ener} = \@ener;
  $ret{avrmsf} = $avg if ($avrmsf);
  $ret{C} = $C;
  if ($measureTime) { printf("--> TIME (analyzeRun): post-processing took %d seconds\n", time() - $time); $time = time(); }
  return \%ret;
}


=head2

 Title   :  uniformlyDistributedTrajectory
 Usage   :  uniformlyDistributedTrajectory($ener, $Cdef);
 Function:  This function is primarily for debugging purposes. Given a previously read trajectory and
            corresponding coordinate definitions, it "flattens" the trajectory (i.e. replaces it with
            one uniformly sampled from the allowed region and with zero energies everywhere).
 Returns :
 Args    :  1. Coordinate definitions
            2. Coordinate array
            3. Energy array
            4. Flag to specify that certain degrees of freedom should be thrown away

=cut

sub uniformlyDistributedTrajectory {
  my $Cdef = shift;
  my $C = shift;
  my $ener = shift;
  my $removeSomeDF = shift;
  $removeSomeDF = 0 if (!defined($removeSomeDF));

  # eliminate non-dihedral degrees of freedom
  my @toRemove;
  if ($removeSomeDF) {
    for (my $ci = scalar(@$Cdef)-1; $ci >= 0 ; $ci--) {
      if (uc($Cdef->[$ci]->{type}) !~ /^(DIHE|IMPR|DIST|ANGL)$/) {
        splice(@$Cdef, $ci, 1);
        splice(@$C, $ci, 1);
      }
    }
  }

  # make sure allowed ranges are defined for all coordinates
  for (my $ci = 0; $ci < scalar(@$Cdef); $ci++) {
    GENERAL::assert(defined($Cdef->[$ci]->{min}) && defined($Cdef->[$ci]->{max}) && GENERAL::isNumeric($Cdef->[$ci]->{min}) && GENERAL::isNumeric($Cdef->[$ci]->{max}), 
    "valid ranges for all coordinates must be specified (not specified for coordinate $ci of type $Cdef->[$ci]->{type}");
    # since we will be generating a uniform distribution not in Cartesian space, but in whatever the degrees of freedom are defined to be, the Jacobian is not
    # necessary (of course, the final integral will correspond to the integral in terms of these given degrees of freedom, and not Cartesian coordinates like
    # we'd want if were actually calculating free energies)
    $Cdef->[$ci]->{jacobian} = 0;
  }

  # flatten the trajectory
  for (my $sni = 0; $sni < scalar(@$ener); $sni++) {
    for (my $ci = 0; $ci < scalar(@$Cdef); $ci++) {
      if (uc($Cdef->[$ci]->{type}) =~ /^(DIHE|IMPR)$/) {
        if (defined($Cdef->[$ci]->{fullCircle})) {
#          $C->[$ci]->[$sni] = -180 + Math::Random::Secure::rand(360);
          $C->[$ci]->[$sni] = -180 + rand(360);
        } else {
#          $C->[$ci]->[$sni] = VALOCIDY::angleDiff(Math::Random::Secure::rand(1)*VALOCIDY::angleDiffCCW($Cdef->[$ci]->{max}, $Cdef->[$ci]->{min}) + $Cdef->[$ci]->{min}, 0);
          $C->[$ci]->[$sni] = VALOCIDY::angleDiff(rand(1)*VALOCIDY::angleDiffCCW($Cdef->[$ci]->{max}, $Cdef->[$ci]->{min}) + $Cdef->[$ci]->{min}, 0);
        }
      } elsif (uc($Cdef->[$ci]->{type}) =~ /^(DIST|AXISDIST|PLANEDIST|INPLACEDIST|ANGL)$/) {
#        $C->[$ci]->[$sni] = Math::Random::Secure::rand(1)*($Cdef->[$ci]->{max} - $Cdef->[$ci]->{min}) + $Cdef->[$ci]->{min};
        $C->[$ci]->[$sni] = rand(1)*($Cdef->[$ci]->{max} - $Cdef->[$ci]->{min}) + $Cdef->[$ci]->{min};
      } else {
        GENERAL::error("Unrecognized coordinate type '$Cdef->[$ci]->{type}'");
      }
    }
    for (my $i = 0; $i < scalar(@{$ener->[$sni]}); $i++) { $ener->[$sni]->[$i] = 0; }
  }
}

# loads a run trajectory and extracts information about structural variations
sub analyzeRunSlowerAllCHARMM {
  my $chrm = shift;
  my $trajf = shift;
  my $dbg = shift;
  my $terms = shift; # additional terms to save from each snapshot
  $terms = "" if (!defined($terms));
  my $Cdef = shift; # definitions of phase space-defining coordinates
  if (!defined($Cdef)) { my @tmp; $Cdef = \@tmp; }
  my $avrmsf = shift; $avrmsf = 0 if (!defined($avrmsf)); # flag - if set, will return average structure and RMSF for each atom

  my $acrdf = "_avg.crd";
  my $cref = "_ef.dat";

  $chrm->send("OPEN UNIT 21 READ FILE NAME $trajf\n");
  if ($avrmsf) {
    $chrm->send("COOR DYNA COMP FIRSTU 21 NUNIT 1 SELE all END ORIENT SELE all END\n");
    $chrm->send("WRITE COOR COMP CARD NAME $acrdf\n* avg coords with rms fluct in weight array\n*\n\n");
    $chrm->send("REWIND UNIT 21\n");
  }
  $chrm->send("TRAJ FIRSTU 21 NUNIT 1\n");
  $chrm->send("SET K = 1\n");
  $chrm->send("OPEN UNIT 11 WRITE CARD NAME $cref\n");
  # loops can not be in an interactive stream, so switch over
  my $tmpsf = "_strea.tmp.inp"  . GENERAL::irand(0, 100000); my $ofh = GENERAL::GetOutFH($tmpsf);

  # first, pre-define all selections (much faster)
  my $k = 1; my %selNames;
  foreach my $cdef (@$Cdef) {
    next if (!defined($cdef->{atoms}));
    foreach my $a (@{$cdef->{atoms}}) {
      next if (defined($selNames{$a}));
      my $sname = "atom$k";
      $selNames{$a} = $sname;
      $ofh->print("define $sname select $a end\n");
      $k++;
    }
  }
  # CHARMM has a limit for the number of definitions allowed (I've currently recompiled my version to have this set to 20000)
  if (scalar(keys(%selNames)) > 20000) {
    GENERAL::error("Number of unique atoms to select for measuring integration coordinates is above CHARMM's limit of 20,000. Recompile CHARMM to increase MAXSKY (in fcm/selcta.fcm) and change the limit here too");
  }
  $ofh->print("LABEL LOOP\n");
  $ofh->print("TRAJ READ\n");
  $ofh->print("energy\n");
  $ofh->print("WRITE TITLE UNIT 11\n");
  $ofh->print("* ?ENER $terms\n*\n\n");
  foreach my $cdef (@$Cdef) {
    if (uc($cdef->{type}) eq "DIHE") {
      $ofh->print("QUICK select $selNames{$cdef->{atoms}->[0]} end select $selNames{$cdef->{atoms}->[1]} end select $selNames{$cdef->{atoms}->[2]} end select $selNames{$cdef->{atoms}->[3]} end\n");
      $ofh->print("WRITE TITLE UNIT 11\n");
      $ofh->print("* ?PHI\n*\n\n");
    } elsif (uc($cdef->{type}) eq "DIST") {
      $ofh->print("QUICK select $selNames{$cdef->{atoms}->[0]} end select $selNames{$cdef->{atoms}->[1]} end\n");
      $ofh->print("WRITE TITLE UNIT 11\n");
      $ofh->print("* ?DIST\n*\n\n");
    } elsif (uc($cdef->{type}) eq "INPLACEDIST") {
      $ofh->print("QUICK select $selNames{$cdef->{atoms}->[0]} end\n");
      $ofh->print("WRITE TITLE UNIT 11\n");
      $ofh->print("* ?XVAL ?YVAL ?ZVAL\n*\n\n");
    } elsif (uc($cdef->{type}) eq "ANGL") {
      $ofh->print("QUICK select $selNames{$cdef->{atoms}->[0]} end select $selNames{$cdef->{atoms}->[1]} end select $selNames{$cdef->{atoms}->[2]} end\n");
      $ofh->print("WRITE TITLE UNIT 11\n");
      $ofh->print("* ?THET\n*\n\n");
    } elsif (uc($cdef->{type}) eq "AXISDIST") {
      # distance from an atom to an axis will be measure by measuring the distance between the ends of the axis and the atom
      $ofh->print("QUICK select $selNames{$cdef->{atoms}->[0]} end\n");
      $ofh->print("WRITE TITLE UNIT 11\n");
      $ofh->print("* ?XVAL ?YVAL ?ZVAL\n*\n\n");
    } elsif (uc($cdef->{type}) eq "PLANEDIST") {
      # distance from an atom to an axis will be measure by measuring the distance between the ends of the axis and the atom
      $ofh->print("QUICK select $selNames{$cdef->{atoms}->[0]} end\n");
      $ofh->print("WRITE TITLE UNIT 11\n");
      $ofh->print("* ?XVAL ?YVAL ?ZVAL\n*\n\n");
    } else {
      GENERAL::error("Unrecognized phase-space coordinate type '$cdef->{type}'");
    }
  }
  $ofh->print("INCR K BY 1\n");
  $ofh->print("IF K LE ?NFILE GOTO LOOP\n");
  close($ofh);
  $chrm->send("stream $tmpsf\n");
  $chrm->closeCard(21); $chrm->closeCard(11);
  $chrm->wait('./', 1, 300, 0); # wait until the buffer is clear (window run is complete)

  # array of arrays that will hold all coordinate value/snapshots, first by coordinate then by snapshot
  my @C = GENERAL::ones(scalar(@$Cdef), 0);
  for (my $i = 0; $i < scalar(@$Cdef); $i++) {
    my @arr; tie @arr, TieArrayC => "double";
    $C[$i] = \@arr;
  }
  # read average structure and RMSF information
  my $avg; $avg = PDB::new($acrdf) if ($avrmsf);
  my $ifh = GENERAL::GetInFH($cref);
  my @ener;
  while (1) {
    my $line = <$ifh>; last if (!defined($line));
    # -- total energy/other terms
    my @line = split(" ", $line);
    GENERAL::assert(scalar(@line) >= 1, "Expected at least one element in each energy line of output file $cref...");
    if (scalar(@line) == 1) {
      push(@ener, $line[0]);
    } else {
      GENERAL::assert(scalar(@{$ener[0]}) == scalar(@line), "Inconsistent number of energy terms inenergy line of output file $cref...") if (scalar(@ener) > 0);
      push(@ener, \@line);
    }

    # -- phase-space coordinates
    for (my $ci = 0; $ci < scalar(@$Cdef); $ci++) {
      my $cdef = $Cdef->[$ci];
      my $line = <$ifh>;
      GENERAL::assert(defined($line), "Inconsistent number of phase-space coordinates in output file $cref...");
      my @line = split(" ", $line);
      if (uc($cdef->{type}) eq "DIHE") {
        GENERAL::assert(scalar(@line) == 1, "Expected only one value in DIHE line of output file $cref '$line'...");
        push(@{$C[$ci]}, $line[0] + 0.0);
      } elsif (uc($cdef->{type}) eq "DIST") {
        GENERAL::assert(scalar(@line) == 1, "Expected only one value in DIST line of output file $cref '$line'...");
        push(@{$C[$ci]}, $line[0] + 0.0);
      } elsif (uc($cdef->{type}) eq "INPLACEDIST") {
        GENERAL::assert(scalar(@line) == 3, "Expected three values in INPLACEDIST line of output file $cref '$line'...");
        push(@{$C[$ci]}, sqrt(($line[0] - $cdef->{origX})**2 + ($line[1] - $cdef->{origY})**2 + ($line[2] - $cdef->{origZ})**2));
      } elsif (uc($cdef->{type}) eq "ANGL") {
        GENERAL::assert(scalar(@line) == 1, "Expected only one value in ANGL line of output file $cref '$line'...");
        push(@{$C[$ci]}, $line[0] + 0.0);
      } elsif (uc($cdef->{type}) eq "AXISDIST") {
        # distance from an atom to an axis will be measure by measuring the distance between the ends of the axis and the atom
        GENERAL::assert(scalar(@line) == 3, "Expected three values in AXISDIST line of output file $cref '$line'...");
        my $d1 = sqrt(($line[0] - $cdef->{axisA}->[0])**2 + ($line[1] - $cdef->{axisA}->[1])**2 + ($line[2] - $cdef->{axisA}->[2])**2); # distance from atom to one end of axis
        my $d2 = sqrt(($line[0] - $cdef->{axisB}->[0])**2 + ($line[1] - $cdef->{axisB}->[1])**2 + ($line[2] - $cdef->{axisB}->[2])**2); # distance from atom to other end of axis
        my $L = sqrt(($cdef->{axisA}->[0] - $cdef->{axisB}->[0])**2 + ($cdef->{axisA}->[1] - $cdef->{axisB}->[1])**2 + ($cdef->{axisA}->[2] - $cdef->{axisB}->[2])**2); # length of axis
        push(@{$C[$ci]}, sqrt(($d1 + $d2 - $L)*($L + $d1 - $d2)*($L + $d2 - $d1)*($L + $d1 + $d2))/(2*$L));
      } elsif (uc($cdef->{type}) eq "PLANEDIST") {
        GENERAL::error("This error is placed here to signal that code for calculating PLANEDIST in VALOCIDY::analyzeRunSlowerAllCHARMM is not yet tested. It is probably correct, but please test before using.");
        # distance from an atom to an axis will be measure by measuring the distance between the ends of the axis and the atom
        GENERAL::assert(scalar(@line) == 3, "Expected three values in PLANEDIST line of output file $cref '$line'...");
        my @v1 = ($line[0] - $cdef->{axisPoint}->[0], $line[2] - $cdef->{axisPoint}->[1], $line[2] - $cdef->{axisPoint}->[2]);
        my @v2 = ($cdef->{axisPoint}->[0], $cdef->{axisPoint}->[1], $cdef->{axisPoint}->[2]);
        my $L = sqrt($v2[0]**2 + $v2[1]**2 + $v2[2]**2);
        $v2[0] /= $L; $v2[1] /= $L; $v2[2] /= $L;
        push(@{$C[$ci]}, $v1[0]*$v2[0] + $v1[1]*$v2[1] + $v1[2]*$v2[2]);
      } else {
        GENERAL::error("Unrecognized phase-space coordinate type '$cdef->{type}'");
      }
    }
  }
  close($ifh);

  GENERAL::crm($acrdf) if ($avrmsf && !defined($dbg));
  GENERAL::crm($cref, $tmpsf) if (!defined($dbg));

  # compose the return hash
  my %ret;
  $ret{ener} = \@ener;
  $ret{avrmsf} = $avg if ($avrmsf);
  $ret{C} = \@C;
  return \%ret;
}

# computes the logarithm of the sum of the exponents of the given arguments - log(sum(exp(x))), where x is the passed array
sub logSumExp {
  my $A = GENERAL::max(@_); # largest exponent, will be used for factorization
  my $B = 0; # log(e^x1 + e^x2 + e^x3 + ...) = log(e^A*[e^(x1-A) + e^(x2-A) + e^(x3-A) + ...]) = A + log(e^(x1-A) + e^(x2-A) + e^(x3-A) + ...)
  foreach my $b (@_) {
    $B += exp($b - $A);
  }
  return $A + log($B);
}

# Given an array of values, computes a weighted average given fractional contributions as logarimths of numerators and denomenators
sub weighByLogFrac {
  my $logn = shift;
  my $logd = shift;
  my $arr = shift;
  my $narrf = (ref($logn) =~ /ARRAY/) ? 1 : 0;
  my $darrf = (ref($logd) =~ /ARRAY/) ? 1 : 0;

  my $av = 0; my ($num, $den);
  for (my $i = 0; $i < scalar(@$arr); $i++) {
    if ($narrf) { $num = $logn->[$i]; }
    else { $num = $logn; }

    if ($darrf) { $den = $logd->[$i]; }
    else { $den = $logd; }

    $av += exp($num - $den) * $arr->[$i] / scalar(@$arr);
  }

  return $av;
}

# Returns 1 if the given angle is between the low and high limits going in the positive (counter-clockwise) direction
sub angleWithin {
  my $a = shift;
  my $mi = shift;
  my $ma = shift;
  my $eps = 10**(-5);
#  GENERAL::assert((abs($a) <= 180) && (abs($mi) <= 180) && (abs($ma) <= 180), "Angles specified must be between -180 and +180 ($a, $mi, and $ma were specified)"); # no need for this since we are doing mods

  # special case: when min and max are the same angle, then the entire circle is allowed, so any angle is good
  return 1 if (abs(mod($mi, 360) - mod($ma, 360)) < $eps);

  # special case: angle is very close to either border, within machine error - count as within
  return 1 if ((abs(mod($a, 360) - mod($mi, 360)) < $eps) || (abs(mod($a, 360) - mod($ma, 360)) < $eps));

  # if counter-clockwise distance between min and max is larger than between min and the angle, then we are within
  if (mod((mod($a, 360) - mod($mi, 360)), 360) <= mod((mod($ma, 360) - mod($mi, 360)), 360)) {
    return 1;
  } else {
    return 0;
  }

}

# Check whether the specified coordinate satisfies the limits of phase space imposed
sub coordinateWithinLimits {
  my $cdef = shift;
  my $cval = shift;
  my $cmin = shift;
  my $cmax = shift;

  my $f = 1;
  if (uc($cdef->{type}) =~ /^(DIHE|IMPR|ANGL)$/) {
    $f = 0 if (!VALOCIDY::angleWithin($cval, $cmin, $cmax));
  } elsif (uc($cdef->{type}) =~ /^(DIST|AXISDIST|PLANEDIST)$/) {
    $f = 0 if (($cval < $cmin) || ($cval > $cmax));
  } elsif (uc($cdef->{type}) eq "INPLACEDIST") {
    $f = 0 if ($cval > $cmax);
  } else {
    GENERAL::error("Unrecognized coordinate type '$cdef->{type}'");
  }
  return $f;
}

# Returns 1 with probability equal to the number given and 0 otherwise.
sub pickWithProbability {
  my $p = shift;
  GENERAL::assert(GENERAL::isNumeric($p) && ($p >= 0) && ($p <= 1), "Probability must be [0; 1], but got '$p'");

  return 1 if ($p == 1);
  return 0 if ($p == 0);

  my $n = 1; # decimal place currently visited (in logs of 10)
  while (1) {
    my $d = GENERAL::irand(0, 9); # digit at current decimal place in random number
    my $di = int($p*(10**($n))) - 10*int($p*(10**($n-1))); # digit at current decimal place of probability cutoff
    if ($di < $d) {
      return 0;
    } elsif ($di > $d) {
      return 1;
    }
    $n++;
  }
}

# reads a trajectory files and winds it up to the given snapshot such that the coordinates of that snapshot are in the main set
sub windToSnapshot {
  my $chrm = shift;
  my $trajf = shift;
  my $sni = shift;

#$chrm->send("prnlev 5\nwrnlev 5\nprint coor\n");
  $chrm->send("SET K = 1\n");
  $chrm->send("OPEN UNIT 21 READ FILE NAME $trajf\n");
  $chrm->send("TRAJECTORY FIRSTU 21\n");
  # loops can not be in an interactive stream, so switch over
  my $tmpsf = "_strea.tmp.inp"  . GENERAL::irand(0, 100000); my $ofh = GENERAL::GetOutFH($tmpsf);
  $ofh->print("LABEL LOOP\n");
  $ofh->print("TRAJ READ\n");
  $ofh->print("INCR K BY 1\n");
  $ofh->print("IF K LE $sni GOTO LOOP\n");
  close($ofh);
  $chrm->send("stream $tmpsf\n");
  $chrm->closeCard(21);
  $chrm->print("ic fill\n");
#$chrm->send("prnlev 5\nwrnlev 5\nprint coor\nprint ic\n");
  $chrm->wait('./', 1, 300); # wait until the buffer is clear (window run is complete)
  GENERAL::crm($tmpsf);

}


# Randomizes the phi/psi angles of the structure until one is encountered with a negative total energy
sub randomizeStructure {
  my $chrm = shift;
  my $seed = int(rand()*100000);

  # Copyright (c) 1994
  # Molecular Simulations Inc.
  # All Rights Reserved
  # This stream file randomly generates values for the dihedral 
  # angles phi and psi in all residues of a protein
  my $tmpsf = "_strea.tmp.inp"  . GENERAL::irand(0, 100000); my $ofh = GENERAL::GetOutFH($tmpsf);
  $ofh->printf("! Main loop\nset SEED $seed\nLABEL START\n");
  # Temporarily save internal coordinate table; keep only those
  # dihedrals that correspond to phi and psi and assign random
  # values to these internal coordinates; print results to
  # output file
  $ofh->printf("IC SAVE\n");
  $ofh->printf("IC KEEP FIRST SELE ATOM * * N .OR. ATOM * * C END\n");
  $ofh->printf("IC KEEP SECO SELE ATOM * * CA .OR. ATOM * * N END\n");
  $ofh->printf("IC KEEP THIRD SELE ATOM * * C .OR. ATOM * * CA END\n");
  $ofh->printf("IC KEEP FOURTH SELE ATOM * * N .OR. ATOM * * C END\n");
  $ofh->printf("IC RAND ISEED \@SEED\n");
  # Restore full internal coordinate table preserving randomly
  # assigned values of phi and psi; initialize and rebuild all
  # Cartesian coordinates
  $ofh->printf("IC REST PRES\n");
  $ofh->printf("COOR INIT\n");
  $ofh->printf("IC SEED 1 N 1 CA 1 C\n");
  $ofh->printf("IC BUIL\n");
  # Compute an energy
  $ofh->printf("mini sd nstep 100\n");
  $ofh->printf("GETE\n");
  # Check value of energy; if less than 0.0, accept conformation; if not, generate another one
  $ofh->printf("SET 2 ?ENER\n");
  $ofh->printf("INCR SEED BY 1\n");
  $ofh->printf("IF 2 GT 0.0 GOTO START\n");
  close($ofh);
  $chrm->send("stream $tmpsf\n");
  $chrm->wait('./', 1, 300); # wait until the structure is successfully randomized
  GENERAL::crm($tmpsf);
}


# Calculates the difference between two angles.
# The difference is guaranteed to be between -180 and 180, where the negative sign means "clockwise" difference and positive "counter-clockwise".
sub angleDiff {
  my $a = shift;
  my $b = shift;
  GENERAL::requireArgs($a, $b);

  my $da = mod((mod($a, 360) - mod($b, 360)), 360);
  $da -= 360 if ($da > 180);

  return $da;
}

# Calculates the difference between two angles, in the counter-clockwise direction. The difference is guaranteed to be between 0 and 360.
sub angleDiffCCW {
  my $a = shift;
  my $b = shift;
  GENERAL::requireArgs($a, $b);

  my $da = mod((mod($a, 360) - mod($b, 360)), 360);

  return $da;
}


=head2

 Title   :  angularStdev
 Usage   :  my $std = angularStdev(\@angles, $min, $max);
 Function:  Computes the standard deviation of an array of angles. Uses a specified reference
            angle (interpreted as the lower bound) to disambiguate angle averaging. See angularMean().
 Returns :
 Args    :  1. reference to an array of angles
            2. lower bound of angles (in the CCW-is-positive sense)

=cut

sub angularStdev {
  my $arr = shift;
  my $mi = shift;
  my $n = scalar(@$arr);

  # compute mean
  my $av = VALOCIDY::angularMean($arr, $mi);

  # compute squares of deviations from mean
  my @d2;
  foreach my $a (@$arr) {
    push(@d2, (VALOCIDY::angleDiff($av, $a))**2);
  }
  return sqrt(GENERAL::mean(\@d2));
}

=head2

 Title   :  angularMean
 Usage   :  my $mean = angularMean(\@angles, $min, $max);
 Function:  Computes the mean of an array of angles. In order to resolve the ambiguities involved
            in averaging angles (what is the average between -PI and PI, is it zero? Or is it PI?),
            it can use the second parameter as the reference lower-bound angle and average the distance
            from that angle to all others in the counter-clockwise direction. If this reference is not
            specified, it will try to find a lower-bound reference that will minimize the average deviation
            from the eventually found mean.
 Returns :
 Args    :  1. reference to an array of angles
            2. optional: lower bound of angles (in the CCW-is-positive sense)

=cut

sub angularMean {
  my $arr = shift;
  my $mi = shift;

  # make a reasonable reference if one is not given
  if (!defined($mi)) {
    # map all angles from 0 to 360 and sort them in ascending order
    my @arr = @$arr;
    for (my $i = 0; $i < scalar(@arr); $i++) { $arr[$i] = mod($arr[$i], 360); }
    @arr = sort {$a <=> $b} @arr;
    # find the biggest gap between consecutive angles
    my @gap = ($arr[-1], $arr[0]); # last gap
    for (my $i = 0; $i < scalar(@arr)-1; $i++) {
      if (VALOCIDY::angleDiffCCW($arr[$i+1], $arr[$i]) > VALOCIDY::angleDiffCCW($gap[1], $gap[0])) {
        @gap = ($arr[$i], $arr[$i+1]);
      }
    }
    $mi = $gap[0] + VALOCIDY::angleDiffCCW($gap[1], $gap[0])/2;
  }

  my $av = 0;
  foreach my $a (@$arr) {
    $av += VALOCIDY::angleDiffCCW($a, $mi);
  }
  $av /= scalar(@$arr);       # now this is the average distance in counter-clockwise direction (always positive) from the lower-bound to the set of angles given
  $av = mod($mi + $av, 360);  # and this is then the average angle
  $av -= 360 if ($av > 180);

  return $av;
}

# Computes the standard deviation of an array of angles. Angles are expected to be in the range of -PI to PI.
# In order to resolve the ambiguities involved in averaging angles (what is the average between -PI and PI, is it zero? Or is it PI?),
# it uses the second parameter as the reference lower-bound angle and averages the distance from that angle to all others in the counter-clockwise
# direction (since the reference angle is interpreted as the lower bound).
sub angularStdevOld {
  my $arr = shift;
  my $mi = shift;
  my $n = scalar(@$arr);

  # compute mean
  my $av = 0;
  foreach my $a (@$arr) {
    $av += VALOCIDY::angleDiffCCW($a, $mi);
  }
  $av /= $n; # now this is the average distance in counter-clockwise direction (always positive) from the lower-bound to the set of angles given
  $av = $mi + $av; # and this is then the average angle
  $av -= 360 if ($av > 180);

  # compute squares of deviations from mean
  my @d2;
  foreach my $a (@$arr) {
    push(@d2, (VALOCIDY::angleDiff($av, $a))**2);
  }
  return sqrt(GENERAL::mean(\@d2));
}


# proper mod
sub mod {
  return $_[0] - floor($_[0] / $_[1]) * $_[1];
}


# Moves structure to the given point in phase space. Currently, supports only dihedral angles as coordinates.
sub moveToPoint {
  my $chrm = shift;
  my $cdef = shift; # phase-space coordinate definitions
  my $p = shift; # the actual point (values of coordinates)
  GENERAL::assert(scalar(@$p) == scalar(@$cdef), sprintf("%d coordinates and %d coordinate definitions given...", scalar(@$p), scalar(@$cdef)));

  $chrm->send("ic fill\n"); # fill the internal coordinates given current Cartesian coordinates
  $chrm->send("coor init select all end\n"); # remove all coordinates
  # edit the IC table
  for (my $ci = 0; $ci < scalar(@$p); $ci++) {
    GENERAL::assert(uc($cdef->[$ci]->{type}) eq "DIHE", "Currently only DIHE coordinate type is supported!");
    $chrm->send("define ggg select $cdef->[$ci]->{atoms}->[0] end\nset A ?SELATOM\n");
    $chrm->send("define ggg select $cdef->[$ci]->{atoms}->[1] end\nset B ?SELATOM\n");
    $chrm->send("define ggg select $cdef->[$ci]->{atoms}->[2] end\nset C ?SELATOM\n");
    $chrm->send("define ggg select $cdef->[$ci]->{atoms}->[3] end\nset D ?SELATOM\n");
    $chrm->send("ic edit\n\n");
    $chrm->send("dihe BYNUM \@A BYNUM \@B BYNUM \@C BYNUM \@D $p->[$ci]\n");
    $chrm->send("end\n");
  }
  $chrm->send("ic seed 1 n 1 ca 1 c\nic build\n\n"); # place all coordinates back given the modified IC
}


# Returns a vector of random numbers picked from a uniform distribution between the given bounds.
sub uniformRandomVector {
  my $n = shift;
  my $mi = shift;
  my $ma = shift;

  my @v;
  for (my $i = 0; $i < $n; $i++) {
    push(@v, rand($ma-$mi) + $mi);
  }
  return @v;
}

# Returns a vector of random numbers picked from a normal distribution of given mean and standard deviation. Optimally,
# samples inly within +/- some value away from the mean.
sub normalRandomVector {
  my $n = shift;
  my $u = shift;
  my $std = shift;
  my $max = shift;

  my @v;
  for (my $i = 0; $i < $n; $i++) {
    while (1) {
      my $v = GENERAL::gaussianRand()*$std;
      if (!defined($max) || (abs($v) <= $max)) {
        push(@v, $v + $u); last;
      }
    }
  }
  return @v;
}

# Vector magnitude
sub norm {
  my $ss = 0;
  foreach my $v (@_) {
    $ss += $v**2;
  }
  return sqrt($ss);
}

# create a 2D matrix
sub ones2D {
  my $x = shift;
  my $y = shift;
  my $val;

  my @arr2D = GENERAL::ones($x, undef);
  for (my $i = 0; $i < scalar(@arr2D); $i++) {
    my @arr = GENERAL::ones($y, $val);
    $arr2D[$i] = \@arr;
  }
  return \@arr2D;
}

# pick a pair of indices according to a given histogram of counts
sub randomlyPickFromHist2D {
  my $hist2D = shift;
  my $N = shift;
  my $n = 0;
  my $r = rand();
  for (my $ii = 0; $ii < scalar(@$hist2D); $ii++) {
    for (my $jj = 0; $jj < scalar(@{$hist2D->[$ii]}); $jj++) {
      $n += $hist2D->[$ii][$jj]; # by keeping track of counts, this accumulation is exact (no numerical error)
      if ($n/$N >= $r) {
        # picked bin $ii, $jj: sample a pair of angles uniformly from it
        return ($ii, $jj);
      }
    }
  }
}

=head2

 Title   :  getDihedralDefinitions
 Usage   :
 Function:  Returns an array set of dihedral definitions of the molecule currently loaded in CHARMM.
 Returns :
 Args    :  1. CHARMM object
            3. Backbone-only flag - if set, only backbone dihedrals will be considered (set by default)
            4. Phi/Psi-only flag - if set, only phi/psi will be considered (not set by default)
            5. Skip impropers flag - if set, impropers will not be considered (set by default)

=cut

sub getDihedralDefinitions {
  my $chrm = shift;
  my $bbo = shift; $bbo = 1 if (!defined($bbo));
  my $phipsio = shift; $phipsio = 0 if (!defined($phipsio));
  my $simp = shift; $simp = 1 if (!defined($simp));

  # Get backbone dihedral definitions
  my $icf = "_$$\_ic.out";
  $chrm->send("ic fill\n");
  $chrm->send("OPEN UNIT 22 WRITE CARD NAME $icf\n");
  $chrm->send("WRITE IC UNIT 22\n* G\n*\n\n"); # unit closes automatically
  $chrm->wait('./', 1, 10); # wait until the buffer is clear (window run is complete)
  my @ICdef; my $ifh = GENERAL::GetInFH($icf);
  while (<$ifh>) {
    next if (/^\s*\*/);
    my @line = split(" ", GENERAL::Trim($_));
    next if (scalar(@line) != 14); # an IC line has 14 entries
    my @names = ($line[2], $line[4], $line[6], $line[8]);
    my @resid = ($line[1], $line[3], $line[5], $line[7]);
    next if (join(" ", @names) =~ /\?\?/); # skip IC with atoms not in the structure
    next if (join(" ", @resid) =~ /-99/); # skip IC with atoms not in the structure
    next if ($bbo && (join(" ", @names) !~ /^(N|C|CA) (N|C|CA) (N|C|CA) (N|C|CA)$/));
    next if ($phipsio && (join(" ", @names) !~ /^(N|C) (CA|N) (C|CA) (N|C)$/));
    my @sels = ("IRES $line[1] .AND. TYPE $line[2]", "IRES $line[3] .AND. TYPE $line[4]", "IRES $line[5] .AND. TYPE $line[6]", "IRES $line[7] .AND. TYPE $line[8]");
    my $dihe = $line[11];
    my $imp = 0;
    if ($names[2] =~ /\*/) {
      $imp = 1; $names[2] =~ s/\*//;
      next if ($simp);
    }
    my %c; $c{type} = "DIHE"; $c{atoms} = \@sels; $c{dihe} = $dihe;
    push(@ICdef, \%c);
  }
  close($ifh);
  GENERAL::crm($icf);

  return @ICdef;
}


=head2

 Title   :  constantTemperatureDynamicsCHARMM
 Usage   :
 Function:  Returns a string that specifies a constant-temperature dynamics run in CHARMM with proper thermostat.
 Returns :
 Args    :

=cut

sub constantTemperatureDynamicsCHARMM {
  my %opts = @_;
  GENERAL::assert(defined($opts{temp}) && defined($opts{steps}), "not all required options specified");
  setDefault(\%opts, "echeck", "");
  setDefault(\%opts, "seed", $$ + int(rand()*100000));
  setDefault(\%opts, "restart", 0);
  setDefault(\%opts, "NSAVC", 100);
  setDefault(\%opts, "IUNREA", -1);
  setDefault(\%opts, "IUNWRI", -1);
  setDefault(\%opts, "ISVFRQ", 100);
  setDefault(\%opts, "IUNCRD", -1);
  setDefault(\%opts, "KUNIT", -1);
  setDefault(\%opts, "QREF", 10);
  setDefault(\%opts, "step", 0.001);
  setDefault(\%opts, "FIRSTT", $opts{temp});
  setDefault(\%opts, "integrator", "VVER");
  setDefault(\%opts, "FBETA", 0.1);
  my $start = "STRT"; $start = "REST" if ($opts{restart});

  my $cmd = "";
  if ($opts{integrator} eq "VVER") {
    # Nose-Hoover thermostat (is supposed to produce the Boltzmann ensemble)
    # VVER is a newer integrator that is more efficient
    $cmd = "
    DYNA VVER NOSE QREF $opts{QREF} TREF $opts{temp} NCYC 5 FIRSTT $opts{FIRSTT} FINALT $opts{temp} -
    $start NSTEP $opts{steps} TIMESTP $opts{step} $opts{echeck} -
    IPRFRQ 2000 IHTFRQ 0 IEQFRQ 0 NTRFRQ 100  -
    IUNREA $opts{IUNREA} IUNWRI $opts{IUNWRI} ISVFRQ $opts{ISVFRQ} IUNCRD $opts{IUNCRD} IUNVEL -1 KUNIT $opts{KUNIT} -
    NPRINT 1000 NSAVC $opts{NSAVC} NSAVV 0 IHBFRQ 0 INBFRQ 25 IMGFrq 25 -
    ISEED $opts{seed} -
    IASORS 0 IASVEL 1 ISCVEL 0 ICHECW 0 TWINDH 0.0 TWINDL 0.0";
  } elsif ($opts{integrator} eq "LEAP") {
    # Langevin dynamics
    $cmd = "SCALAR FBETA SET $opts{FBETA} select .not. hydrogen end
    DYNA LEAP LANGEVIN $start NSTEP $opts{steps} TIMESTEP $opts{step} $opts{echeck} -
    IPRFRQ 2000 IHTFRQ 0 IEQFRQ 0 NTRFRQ 0  -
    IUNREA $opts{IUNREA} IUNWRI $opts{IUNWRI} IUNCRD $opts{IUNCRD} IUNVEL -1 KUNIT $opts{KUNIT} -
    NPRINT 1000 NSAVC $opts{NSAVC} NSAVV 0 IHBFRQ 0 INBFRQ 25 IMGFrq 25 -
    TBATH $opts{temp} ISEED $opts{seed} -
    FIRSTT $opts{FIRSTT} FINALT $opts{temp}  -
    IASORS 0 IASVEL 1 ISCVEL 0 ICHECW 0 TWINDH 0.0 TWINDL 0.0";
    if (defined($opts{selfGuided})) {
      $cmd .= " -\n    SGLD TSGAVG 0.2 SGFT 1.0";
      $cmd .= " SGSTOPT SGSTOPR" if (!defined($opts{selfGuidedMobileCOM}) || !$opts{selfGuidedMobileCOM});
    }
  } else {
    GENERAL::error("unknown integrator requested '$opts{integrator}'");
  }

  return $cmd . "\n";

  # Simple Verlet integrator with constant temperature -- the simple way of keeping temperature constant; not sure if it gives a proper ensemble
#  return "
#      DYNA LEAP VERLET $start NSTEP $opts{steps} TIMESTEP 0.001 $opts{echeck} -
#      IPRFRQ 2000 IHTFRQ 0 IEQFRQ 0 NTRFRQ 1000  -
#      IUNREA $opts{IUNREA} IUNWRI $opts{IUNWRI} IUNCRD $opts{IUNCRD} IUNVEL -1 KUNIT $opts{KUNIT} -
#      NPRINT 1000 NSAVC $opts{NSAVC} NSAVV 0 IHBFRQ 0 INBFRQ 25  -
#      TCONst  TCOUpling 5.0  TREFerence $opts{temp} ISEED $opts{seed} -
#      FIRSTT $opts{temp} FINALT $opts{temp}  -
#      IASORS 0 IASVEL 1 ISCVEL 0 ICHECW 0 TWINDH 0.0 TWINDL 0.0\n";
}


=head2

 Title   :  setupCHARMM
 Usage   :
 Function:  A generic CHARMM setup function for use with VALOCIDY/CHARMM calculations for uniformity.
 Returns :  An initialized CHARMM object.
 Args    :

=cut

sub setupCHARMM {
  my $opts = shift;
  my %opts = %$opts;
  GENERAL::assert(scalar(@_) % 2 == 0, "expected 'key' -> 'value' pairs");
  for (my $i = 0; $i < scalar(@_); $i += 2) {
    $opts{$_[$i]} = $_[$i+1];
  }

  # decide on interaction cutoffs
  my $cutnb = 7.0; my $ctonnb = 6.0; my $ctofnb = 6.5;
  if (defined($opts{d})) {
    $CUSTOM_PARAM_DEF = $opts{d};
    my $par = GENERAL::LoadStruct(DEFINITIONS::getParmFile("CHARMM.param"), 1);
    $cutnb = $par->{cutnb};
    $ctonnb = $par->{ctonnb};
    $ctofnb = $par->{ctofnb};
  }

  # Setup system
  my $chrm = CHARMM::new(DEFINITIONS::getExec("charmm.33b1.Gevorg.MMFP"), $opts{chcmdlog}, $opts{co}, (defined($opts{co}) || defined($opts{dbg})), 1);
  $chrm->{par}->{param} = $opts{par}; my %parfs;
  if ($opts{e} eq 2) { %parfs = $chrm->loadParm('EEF1'); }
  elsif ($opts{e} eq 4) { %parfs = $chrm->loadParm('IMM1'); }
  else { %parfs = $chrm->loadParm(''); }
  $opts->{toppar} = \%parfs;

  # Set up peptide in given conformation
  $chrm->setupFromCRD($opts{icrdf}, $opts{ter});
  $chrm->send("hbuild\n");

  # Set up energy function
  if (!GENERAL::isInteger($opts{e})) {
    $opts{e} =~ s/\\n/\n/g;
    $chrm->verbose($opts{e});
  } elsif ($opts{e} == 1) {
    $chrm->verbose("update atom rdiel switch eps 1.0 e14fac 0.4 cutnb $cutnb ctonnb $ctonnb ctofnb $ctofnb nbxmod 5 vswitch wmin 0.0 vatom vdistance");
    $chrm->verbose("skip all exclude BOND ANGLE DIHE IMPR VDW ELEC");
  } elsif ($opts{e} == 2) {
    $chrm->{par}->{eef_temp} = $opts{eeft};
    $chrm->setupEEF1();
    $chrm->verbose("skip all exclude BOND ANGLE DIHE IMPR VDW ELEC ASP");
  } elsif ($opts{e} == 3) {
    $chrm->verbose("skip all exclude BOND ANGLE DIHE IMPR VDW");
  } elsif ($opts{e} == 4) {
    $chrm->{par}->{eef_temp} = $opts{eeft};
    $chrm->setupIMM1(30.0);
  } else {
    die "Error: unknown energy model picked - $opts{e}\n";
  }

  # Fix things if necessary
  if (defined($opts{f})) {
    my $bb = CHARMM::backbone();
    $opts{f} =~ s/BACKBONE/\($bb\)/;
    $chrm->verbose("cons fix sele $opts{f} end");
    if (defined($opts{saveFixed})) {
      $chrm->writePDB($opts{saveFixed}, "sele $opts{f} end");
      $chrm->wait();
    }
  }

  # SHAKE
  $chrm->send("shake bond param\n") if ($opts{shake} eq "all");
  $chrm->send("shake bonh param\n") if ($opts{shake} eq "hyd");

  # Randomize starting structure, if necessary
  if ($opts{rnd}) {
    printf("Randomizing starting structure...\n");
    randomizeStructure($chrm);
  }

  return $chrm;
}

sub setDefault {
  my $hash = shift;
  my $key = shift;
  my $defval = shift;

  $hash->{$key} = $defval if (!defined($hash->{$key}));
}

# --------------------------------------- Old functions ----------------------------------------------
=head2

 Title   :  globalSampling_old
 Usage   :
 Function:  Implements the MODISTE approach.
 Returns :
 Args    :  1. A pre-setup CHARMM structure with the system loaded, ideally already within the desired region of phase space.
            3. A reference to array of IC definitions used to defined phase space. For each IC, the entry is an array of four strings which are
               interpreted as CHARMM selections corresponding to the four atoms defining the IC.
            2. Hash reference defining various parameters of the dynamics. Currently, the following parameters are supported (hash key - explanation):
               ne   - number of equilibration steps in femtoseconds (before sampling begins; during this time the energy/temperature equilibrate)
               nc   - number of data-collections steps in femtoseconds (the amount of actual high-temperature sampling)
               ns   - structure sampling interval in femtoseconds
               T    - simulation temperature
               seed - optional: value of random seed for repeating realizations (picked randomly from time and pid by default)

=cut
# 
# sub globalSampling_old {
#   my $chrm = shift;
#   my $ICdef = shift;
#   my $parSYM = shift;
# 
#   # setup some default parameters
#   $parSYM->{seed} = $$ + int(rand()*100000) if (!defined($parSYM->{seed}));
#   my $R = 1.9858775/1000;
# 
#   my %CIdef; # hash IC definitions such that we can quickly check whether a certain IC is within the set of phase-space coordinates
#   for (my $ici = 0; $ici < scalar(@$ICdef); $ici++) {
#     $ICdef{"$ICdef->[$ici]->[0] - $ICdef->[$ici]->[1] - $ICdef->[$ici]->[2] - $ICdef->[$ici]->[3]"} = $ici;
#     $ICdef{"$ICdef->[$ici]->[3] - $ICdef->[$ici]->[2] - $ICdef->[$ici]->[1] - $ICdef->[$ici]->[0]"} = $ici;
#   }
# 
#   # equilibration phase
#   my $rstf = "_ws_equil_restart.rst";
#   $chrm->send("open unit 49 write card name \"$rstf\"\n");
#   $chrm->send("
#       DYNA LEAP VERLET STRT NSTEP $parSYM->{nhe} TIMESTEP 0.001 -
#       IPRFRQ 2000 IHTFRQ 0 IEQFRQ 0 NTRFRQ 1000  -
#       IUNREA -1 IUNWRI 49 IUNCRD -1 IUNVEL -1 KUNIT -1 -
#       NPRINT 1000 NSAVC 100 NSAVV 0 IHBFRQ 0 INBFRQ 25  -
#       TCONst  TCOUpling 5.0  TREFerence $parSYM->{Th} ISEED $parSYM->{seed} -
#       FIRSTT $parSYM->{Th} FINALT $parSYM->{Th}  -
#       IASORS 0 IASVEL 1 ISCVEL 0 ICHECW 0 TWINDH 0.0 TWINDL 0.0\n");
#   $parSYM->{seed}++;
# 
#   # data collection phase
#   $chrm->send("open unit 49 read card name \"$rstf\"\n");
#   my $outdcd = "_wsWS_datacoll_run.dcd";
#   my $outene = "_ws_datacoll_run.ener";
#   my $nrstf = "_ws_datacoll_1.rst";
#   $chrm->send("open unit 50 write file name \"$outdcd\"\n");
#   $chrm->send("open unit 51 write card name \"$outene\"\n");
#   $chrm->send("open unit 52 write card name \"$nrstf\"\n");
#   $chrm->send("
#     DYNA LEAP VERLET REST NSTEP $parSYM->{nhc} TIMESTEP 0.001 -
#     IPRFRQ 2000 IHTFRQ 0 IEQFRQ 0 NTRFRQ 1000  -
#     IUNREA 49 IUNWRI 52 IUNCRD 50 IUNVEL -1 KUNIT 51 -
#     NPRINT 1000 NSAVC $parSYM->{nhs} NSAVV 0 IHBFRQ 0 INBFRQ 25  -
#     TCONst  TCOUpling 5.0  TREFerence $parSYM->{Th} ISEED $parSYM->{seed} -
#     FIRSTT $parSYM->{Th} FINALT $parSYM->{Th}  -
#     IASORS 0 IASVEL 1 ISCVEL 0 ICHECW 0 TWINDH 0.0 TWINDL 0.0\n");
#   $chrm->closeCard(50);
#   $chrm->closeCard(51);
#   $chrm->closeCard(52);
#   $chrm->wait('./', 1, ($parSYM->{nhe} + $parSYM->{nhc})*0.001); # wait until the buffer is clear (data-collection run is complete)
# 
#   # Analyze run: extract energies and ICs
#   my ($avg, $rmsf, $ener, $IC) = analyzeRun($chrm, $outdcd);
# 
#   # Create a rearranged IC matrix so that it is by IC index first, and snapshot next. Also, in this matrix, preserve only those ICs that phase space is defined by
#   # (and in the same order in which bounds on these ICs are specified).
#   my @IC;
#   for (my $sni = 0; $sni < scalar(@{$IC->{dihe}}); $sni++) {
#     my @snapIC = GENERAL::ones(scalar(@$ICdef), -999); # array of ICs (only those defining phase space) in this snapshot
#     for (my $ici = 0; $ici < scalar(@{$IC->{dihe}->[$sni]}); $ici++) {
#       my $icstr = join(" - ", @{$IC->{sels}->[$ici]});
#       next if (!defined($ICdef{$icstr})); # skip this IC if it is not one with which phase space is defined
#       $snapIC[$ICdef{$icstr}] = $IC->{dihe}->[$sni]->[$ici];
#     }
#     foreach my $ic (@snapIC) { GENERAL::assert($ic ne -999, "After analyzing snapshot $sni, not all phase-space defining ICs were discovered!"); }
#     GENERAL::assert(scalar(@snapIC) == scalar(@$ICdef), "After analyzing snapshot $sni, found a different number of phase-space defining ICs than what was expected!");
#     for (my $i = 0; $i < scalar(@snapIC); $i++) {
#       push(@{$IC[$i]}, $snapIC[$i]);
#     }
#   }
#   my $N = scalar(@{$IC[0]}); # number of snapshots
#   my $m = scalar(@IC); # number of phase-space defining coordinates
#   my $V = 360**$m; my $logV = $m*log(360); # total volume of phase space
# 
#   # Try to compute Q at the high temperature directly
#   my (@logQhf, @Qhf); # factors, the sum of which is Q
#   foreach my $en (@$ener) {
#     push(@Qhf, exp($en/$R/$parSYM->{Th}));
#     push(@logQhf, $en/$R/$parSYM->{Th});
#   }
#   my $Qh = $V/GENERAL::mean(\@Qhf); my $stdQhf = GENERAL::stdev(\@Qhf);
#   my $Qherr = $V/(GENERAL::mean(\@Qhf) - $stdQhf) - $V/(GENERAL::mean(\@Qhf) + $stdQhf);
#   printf("At high temperature, Q is estimated to be %e +/- %e (<exp[bE]>_Th = %e +/- %e, V = %e)\n", $Qh, $Qherr, GENERAL::mean(\@Qhf), $stdQhf, $V);
# 
#   # Decide on delta's for each phase-space coordinate for expanding sampled structures, based on how much they vary in sampling
#   # (roughly compute the range each dihedral angle needs to have so that the trajectory covers all of phase space)
# #  my $delA = 180/($N**(1/$m));
#   my $delA = 30;
#   printf("Phase space around sampled structures will be expanded by +/- %f degrees\n", $delA);
#   my $Vi = $delA**$m; my $logVi = $m*log($delA); # volume of phase-space elements that will be locally integrated
# 
#   # Monte-Carlo Integration
#   my $sni = undef; my $pi = undef; my @A; my $Qa = 0; my $logQa = 0; my $nL = $parSYM->{N};
#   for (my $i = 0; $i < $nL; $i++) {
#     # pick an area of phase space out of the space covered by the high-temp run +/- alpha expansion
#     while (1) {
#       # pick a random snapshot
#       $sni = GENERAL::irand(0, $N-1);
# 
#       # find all snapshots in the trajectory within the deltaA voxel of this snapshot and compute p_i
#       $pi = 0;
#       for (my $snj = 0; $snj < $N; $snj++) {
#         my $f = 1;
#         for (my $ici = 0; $ici < $m; $ici++) {
#           if (abs(angleDiff($IC[$ici]->[$snj], $IC[$ici]->[$sni])) > $delA) {
#             $f = 0; last;
#           }
#         }
#         $pi += $f;
#       }
#       $pi /= $N;
# 
#       # pick this snapshot with probability p_i
# #      last if (pickWithProbability($pi));
#       last if (pickWithProbability(1));
#     }
# 
#     # create limits of phase-space defining variables
#     my (@ICmi, @ICma);
#     for (my $ici = 0; $ici < $m; $ici++) {
#       push(@ICmi, angleDiff($IC[$ici]->[$sni], $delA));
#       push(@ICma, angleDiff($IC[$ici]->[$sni], -$delA));
#     }
# 
#     # if picked, compute the local low-temperature partition function corresponding to this snapshot
#     printf("> On iteration %d picked snapshot %d with probability %f (E = %f)...\n", $i+1, $sni+1, $pi, $ener->[$sni]);
#     windToSnapshot($chrm, $outdcd, $sni+1); # wind to the picked snapshot so that we start within the desired phase-space region
#     my @a = integrateMinimum($chrm, \@ICmi, \@ICma, $ICdef, $parSYM);
#     push(@A, @a);
#     $Qa += $V*$a[2]/$pi/$Vi;
#     if (defined($logQa)) { $logQa = logSumExp($logQa, $logV + $a[4] - log($pi) - $logVi); }
#     else { $logQa = $a[4] - log($pi); }
#     printf(">>> Current estimate of Qall = %e (-RTlogQ = %f kcal/mol)\n\n", $Qa/($i+1), -$R*($parSYM->{T})*($logQa - log($i+1)));
#   }
#   $Qa /= $nL;
#   $logQa -= log($nL);
# 
#   return $Qa;
# 
# }

# -------------- OLD CODE -------------------------------

#   # -- calculate the density expectation analytically (needs to be around 1/Vi)
#   my $pE = 0;
#   for (my $sni = 0; $sni < $N; $sni++) {
# print "$sni\n";
#     my $tmp = 1;
#     # squared terms
#     for (my $ci = 0; $ci < $m; $ci++) {
#       $tmp *= erf(180/$sig)/(2 * $sig * sqrt(M_PI) * ((erf(180/($sig * sqrt(2))))**2));
#     }
#     $pE += $tmp;
#     for (my $snj = $sni+1; $snj < $N; $snj++) {
#       $tmp = 1;
#       for (my $ci = 0; $ci < $m; $ci++) {
#         my $a = abs(angleDiff($IC[$ci]->[$sni], $IC[$ci]->[$snj]));
#         $tmp *= (exp(-($a**2)/(4*$sig**2)) * erf((360 - $a)/(2*$sig)) + exp(-((360 - $a)**2)/(4*$sig**2)) * erf($a/(2*$sig))) / (2 * $sig * sqrt(M_PI) * ((erf(180/($sig * sqrt(2))))**2));
#       }
#       $pE += 2*$tmp;
#     }
#   }
#   $pE = $pE / ($N**2);
#   printf("The expectation of probability density is %e, whereas 1/Vi = %e\n", $pE, 1/$Vi);
# 
#   # -- optimize value of sigma such that the sampling distribution is of appropriate "wideness" (density expectation is roughly 1/Vi)
#   # use bisection method
#   my $csig = $sig; my (@a, @b); # initial trial and solution bracket (upper and lower bounds of solution)
#   while (1) {
#     printf("Tring sigma = %f degrees\n", $csig);
#     # -- calculate the density expectation numerically
#     my @pE; my $nspe = 1000;
#     for (my $i = 0; $i < $nspe; $i++) {
#       # pick a random snapshot
#       my $sni = GENERAL::irand(0, $N-1);
# 
#       # pick a normally-distributed random displacement vector
#       my @p = normalRandomVector($m, 0, $csig, 180); my @dp;
#       for (my $ci = 0; $ci < $m; $ci++) {
#         $p[$ci] = angleDiff($IC[$ci]->[$sni], -$p[$ci]);
#         push(@dp, angleDiff($p[$ci], $IC[$ci]->[$sni]));
#       }
# 
#       # compute the PDF of this pick
#       my $pi = 0;
#       for (my $snk = 0; $snk < $N; $snk++) {
#         my $tmp_p = 1;
#         for (my $ci = 0; $ci < $m; $ci++) {
#           $tmp_p *= $Ck * exp(-((angleDiff($p[$ci], $IC[$ci]->[$snk]))**2)/(2*($csig**2)));
#         }
#         $pi += $tmp_p;
#       }
#       $pi = $pi/$N;
#       push(@pE, $pi);
#     }
#     my $pE = GENERAL::mean(\@pE);
#     my $dmpE = GENERAL::stdev(\@pE)/sqrt($nspe); # deviation of the mean
#     printf("After %d sample points, average density expectation is %e +/- %e\n", $nspe, $pE, $dmpE);
#     if (abs($pE - 1/$Vi) < $dmpE) {
#       # solution found
#       $sig = $csig;
#       printf("Solution found, sig = %f leads to density expectation is %e +/- %e and 1/Vi = %e\n", $sig, $pE, $dmpE, 1/$Vi);
#       last;
#     } elsif ($pE > 1/$Vi) {
#       # overshot
#       @b = ($pE, $csig);
#       if (scalar(@a) != 0) {
#         $csig = ($a[1] + $b[1])/2; # mid point
# printf("New upper bound found, going for mid-point...\n");
#       } else {
#         $csig = 1.5*$csig; # sigma is too small
# printf("New upper bound found, looking for lower bound...\n");
#       }
#     } elsif ($pE < 1/$Vi) {
#       # undershot
#       @a = ($pE, $csig);
#       if (scalar(@b) != 0) {
#         $csig = ($a[1] + $b[1])/2; # mid point
# printf("New lower bound found, going for mid-point...\n");
#       } else {
#         $csig = $csig/1.5; # sigma is too large
# printf("Lower bound found, looking for upper bound...\n");
#       }
#     }
#   }

1;



__DATA__
__C__

double logSumExpC (double logI, double logJ) {
  double max = (logI > logJ) ? logI : logJ;
  double del = (logI > logJ) ? logJ - max : logI - max;
  return max + log(1 + exp(del));
}

double mymod (double a, double b) {
  return a - floor(a / b) * b;
}

double angleDiff_c (double a, double b) {
  double da = mymod((mymod(a, 360) - mymod(b, 360)), 360);
  if (da > 180) { da -= 360; }
  return da;
}

double angleDiffCCW_c (double a, double b) {
  return mymod((mymod(a, 360) - mymod(b, 360)), 360);
}

int angleWithin_c (double a, double mi, double ma) {
  double eps = pow(10, -5);

  // special case: when min and max are the same angle, then the entire circle is allowed, so any angle is good
  if (fabs(mymod(mi, 360) - mymod(ma, 360)) < eps) { return 1; }

  // special case: angle is very close to either border, within machine error - count as within
  if ((fabs(mymod(a, 360) - mymod(mi, 360)) < eps) || (fabs(mymod(a, 360) - mymod(ma, 360)) < eps)) { return 1; }

  // if counter-clockwise distance between min and max is larger than between min and the angle, then we are within
  if (mymod((mymod(a, 360) - mymod(mi, 360)), 360) <= mymod((mymod(ma, 360) - mymod(mi, 360)), 360)) {
    return 1;
  } else {
    return 0;
  }
}


// _traj - trajectory double array: by coordinate index first, snapshot index second
// _type - array of coordinate types. 1 - distance-like (bond length, bond angle, ...), 2 - angle-like (dihedral or improper diheral angle)
// _min  - array of coordinate lower bounds of integration
// _max  - array of coordinate upper bounds of integration
// _logD - array of length equal to the number of snapshots which will be populated with logs of functionl values at each snapshot (when divided by the overall integral, probability densities)
// Returns the value of the overall integral (normalization constant)
double GaussianKernelDensity (AV* _traj, AV* _type, AV* _min, AV* _max, AV* _logD) {
  double **traj, *std, *logD, *max, *min, logI, logIi, logVij, logpref, d;
  int *type, ns, nc, nic, i, j, c;
  AV *arr;

  // get dimensions
  nc = av_len(_traj) + 1;                         // number of coordinates
  if (nc == 0) { return 0; }
  ns = av_len(SvRV(*av_fetch(_traj, 0, 0))) + 1;  // number of snapshots

  // allocate arrays
  traj = (double**) malloc(nc*sizeof(double*));
  for (i = 0; i < nc; i++) { traj[i] = (double*) malloc(ns*sizeof(double)); }
  type = (int*) malloc(nc*sizeof(double));
  min = (double*) malloc(nc*sizeof(double));
  max = (double*) malloc(nc*sizeof(double));
  std = (double*) malloc(nc*sizeof(double));
  logD = (double*) malloc(ns*sizeof(double));

  // copy data
  nic = nc; // number of coordinates over which will integrate
  for (i = 0; i < nc; i++) {
    type[i] = (int) SvIV(*av_fetch(_type, i, 0));
    min[i] = (double) SvNV(*av_fetch(_min, i, 0));
    max[i] = (double) SvNV(*av_fetch(_max, i, 0));
    switch (type[i]) {
      case -1:                   // special "do not integrate" type
        nic--; continue; break;
      case 1:
        d = max[i] - min[i];
        break;
      case 2:
        d = fabs(angleDiffCCW_c(max[i], min[i]));
        if (d == 0) { d = 360; } // since we are on a circle, max == min means the whole circle is the range
        break;
      default:
        printf("Error in GaussianKernelDensity: don't know coordinate type %d\n", type[c]);
        exit(-1);
    }
    if (d <= 0) { printf("Error in GaussianKernelDensity: for coordinate %d, type %d, min and max don't make sense: [%f; %f]\n", c+1, type[c], min[c], max[c]); exit(-1); }
    std[i] = d/ns; // standard deviation is set "dynamically", based on the range of this df and how many sample points we have
    arr = SvRV(*av_fetch(_traj, i, 0));
    for (j = 0; j < ns; j++) {
      traj[i][j] = (double) SvNV(*av_fetch(arr, j, 0));
    }
  }

  // Gaussian pre-factor (the product of 1/sqrt(2*pi)/sigma_i for all coordinates)
  logpref = -(nic/2.0)*log(2*M_PI);
  for (i = 0; i < nc; i++) {
    if (type[c] == -1) { continue; }
    logpref -= log(std[i]);
  }

  // compute function value at each snapshot
  for (i = 0; i < ns; i++) {
    logD[i] = 0;
    for (j = 0; j < ns; j++) {
      logVij = logpref;  // the log of the value of the function at point i due to the kernel at j
      for (c = 0; c < nc; c++) {
        switch (type[c]) {
          case -1: continue; break; // special "do not integrate" type
          case 1:
            d = traj[c][i] - traj[c][j];
            break;
          case 2:
            d = angleDiff_c(traj[c][i], traj[c][j]);
            break;
          default:
            printf("Error: don't know coordinate type %d\n", type[c]);
            exit(-1);
        }
        logVij -= d*d/(2*std[c]*std[c]);
      }
      if (j == 0) { logD[i] = logVij; }
      else { logD[i] = logSumExpC(logD[i], logVij); }
    }
  }

  // compute integral
  logI = 0;
  for (i = 0; i < ns; i++) {
    // integral of multi-variable Gaussian due to snapshot i (integral of products of independent Gaussians = product of integrals of Gaussians)
    logIi = 0;
    for (c = 0; c < nc; c++) {
      switch (type[c]) {
        case -1: continue; break; // special "do not integrate" type
        case 1:
          if ((max[c] < traj[c][i]) || (min[c] > traj[c][i])) { printf("Error in GaussianKernelDensity: coordinate %d, type %d has value %f, which is outside of [%f; %f]\n", c+1, type[c], traj[c][i], min[c], max[c]); exit(-1); }
          logIi += log(0.5*erf((max[c] - traj[c][i])/sqrt(2)/std[c]) - 0.5*erf((min[c] - traj[c][i])/sqrt(2)/std[c]));
          break;
        case 2:
          if (!angleWithin_c(traj[c][i], min[c], max[c])) { printf("Error in GaussianKernelDensity: coordinate %d, type %d has value %f, which is outside of [%f; %f]\n", c+1, type[c], traj[c][i], min[c], max[c]); exit(-1); }
          logIi += log(0.5*erf(angleDiffCCW_c(max[c], traj[c][i])/sqrt(2)/std[c]) + 0.5*erf(angleDiffCCW_c(traj[c][i], min[c])/sqrt(2)/std[c]));
          break;
        default:
          printf("Error in GaussianKernelDensity: don't know coordinate type %d\n", type[c]);
          exit(-1);
      }
    }
    if (i == 0) { logI = logIi; }
    else { logI = logSumExpC(logI, logIi); }
  }
  
  // copy data back
  for (i = 0; i < ns; i++) {
    av_store(_logD, i, newSVnv(logD[i]));
  }

  // free arrays
  for (i = 0; i < nc; i++) { free(traj[i]); }
  free(traj);
  free(type);
  free(min);
  free(max);
  free(std);
  free(logD);

  return logI;
}

