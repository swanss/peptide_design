#
# module ENERGY for Design
#
# POD documentation - main docs before the code

=head1 NAME

ENERGY module:  calculate all the eneregy terms.

=head1 SYNOPSIS

=head1 DESCRIPTION

    Interface between several Modules. DO all the necesseay 
    calculations in the reference state, apo state, complexed state.
    Also including the function for post design refinement calculation 
    and energy analysis.

    (1) pairwise vdw, electrostatic, hbond
    (2) self: vdw, electrostatic, solvation, hbond 
    (3) a variety of solvation implementation
    (4) knowledge based potential: a hybrid approach


=head1 EXAMPLES

=head2 Reporting Bugs

=head1 AUTHOR Jiangang Chen

chen2001@mit.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a "_".

=cut



package ENERGY;
#$DEF_VER_CONTROL->{ver}->{ENERGY} = "1.0";
$DEF_VER_CONTROL->{date}->{ENERGY} = "06/17/04";

use strict;
use Parallel::MPI::Simple;
use CHARMM;
use PDB;
use GENERAL;
#use STRIDE;
use ROTAMER;
use POSIX;
use DEFINITIONS;
use Exporter ();

# Set up the energy tracker data structure - it will be a global
# variable so that we don't have to keep reading it from a file
# and only have to write it on updates. Any script/module importing
# this module will be able to initialize this structure. By default
# it will start out empty.
our @ISA = qw(Exporter);
our @EXPORT = qw($TERM_TRACKER);
our $TERM_TRACKER = ();
$TERM_TRACKER->{pair} = {};
$TERM_TRACKER->{self} = {};
$TERM_TRACKER->{unfold} = {};



=head2 

 Title   :  StatisticalPotential
 Usage   :  ENERGY::StatistaicalPotential(rotamer object)
 Function:  convert dunbrack's statistical potential into pseudo-free energy
            using deltaG = -RTln(P) where P is the probability of given rotamer
 To-do   :  (1) add a valid check. don't calculate those with flag of large energy (clash)  
            (2) support subrotamer??? 
                the expansion of subrotamer make this potential slightly complicated. 
                The solution is: do it after subrotamer mean field average calculation. 
            (3) is there a better way to do it????
            (4) may be i should including entropy term. by using 
                (Pi)*ln(Pi) term. But I am not so convincing to add it. 
            (5) how to support modified rotamer library????
 Returns :
 Args    :

=cut

sub StatisticalPotential {
  my $rotamer=shift;

  my $energy = GENERAL::GetOutFH($FS_DEF->{stat_pot_tabf});

  my $sitelist= $rotamer->{dslist};
  foreach my $site (@{$sitelist})   {
    foreach my $res ( @{$site->{reslist}})  {
      foreach my $rot ( @{$res->{rotamerlist}}) {
          if ($rot->{valid} == 1) {
          # note: here we are using log10 function from
          # POSIX Module. The default built-in log function
          # is natural logarithm.
          # we only want to collect data for those without
          # being eliminaed by self energy calculation
          my $G= -1.37*log10($rot->{probability});
          $energy->printf("%7.4f\n", $G);
        }
      }
    }
  }
}


=head2 

 Title   :  writeScrwlSeq  
 Usage   : 
 Function:  write the sequence file the SCRWL. 
            basically, it defined the sites where side placement will be done
            with low case letter.  what this function do is: 
            for the protein with N design sites, creat N sequence input files, 
            each input file contains: 
                 (A) fix one deign site    (upper case) 
                 (B) define the rest of N design sites  (low case)
 Example : The problem: 

           The Solution: see the example code below

 Returns : 
 Args    : (1) Rotamer object 
           (2) chain id for the rotamer to be fixed
           (3) indexed residue number( renumbered residue number) for the rotamer to be fixed
           (4) pdb file fr origin protein structure  

=cut


sub writeScrwlSeq {
    my $rotamer=shift;
    my $chain=shift;
    my $iresnum=shift;
    my $pdbfile =shift;
    my $fname=$chain.$iresnum.'.seq';
    
    my $tempfasta='temp.fasta';
    system("stride $pdbfile -r$chain -q$tempfasta 1>& /dev/null");

    my $temp = &GENERAL::GetInFH('temp.fasta');
    my $out= GENERAL::GetOutFH($fname);

    my $seq;
    while (<$temp>) {
    next if ( /^$/ );
    next if ( /^>/ );
    chomp($_);
    $seq .=$_;
    }

    my @seq=split("", $seq);

    # define the position need to re-packed, we include the place we want to 
    # fix here also 

    for (my $j=0; $j <= @{$rotamer->{dslist}}; $j++ )  {
      # get the interface residue index
      my $resno=$rotamer->{dslist}->[$j]->{iresnum}-1;
      # low case for interface residues
      $seq[$resno]=lc($seq[$resno]);
    }

    # now define the rotamer position I want to keep fixed

    my $resname;
    my $k=0;
    for (my $i=0; $i <=$#seq; $i++) {
      $k++;
      if ( $i == ($iresnum-1)) {
        $resname=uc($seq[$i]);
      } else {
        $resname=$seq[$i];
      }
      if ($k < 60) {
        $out->print($resname);
      } elsif ($k == 60) {
        # one line for 60 residues
        $out->print($resname);
        $out->print("\n");
        # reset counter
        $k=0;
      }
    }

    &GENERAL::Remove('temp.fasta');
}


=head2

 Title   :  addTrackerTerm
 Usage   :  my $term = ENERGY::addTrackerTerm($struct_ptr, $category, $term_tag, $term_description, @colomn_names);
 Function:  Creates a new energy term in the term tracker structure. This function is just created for
            convenience since this action has to be taken by every function that creates an energy term.
 Returns :  A pointer to the newly inserted term.
 Args    :  (1) pointer to the term tracker structure
            (2) category of the term (e.g. unfold, self, pair)
            (3) term tag - a unique term id, using which it will be hashed into the table
            (4) description of the term (a string)
            (5) array of column names
=cut
sub addTrackerTerm {
  my $tt = shift;
  my $category = shift;
  my $tag = shift;
  my $desc = shift;
  my @col_names = @_;

  my $en;
  if (defined($tt->{$category})) {
    $en = $tt->{$category};
  } else {
    $tt->{$category} = {};
    $en = $tt->{$category};
  }
  $en->{$tag} = {};
  $en = $en->{$tag};
  $en->{name} = $desc;
  $en->{'eschemes'} = {};
  $en->{col_names} = \@col_names;
  $en->{files} = ();

  return $en;
}

=head2

 Title   :  addSchemeToTerm
 Usage   :  my $scheme = ENERGY::addSchemeToTerm($term, $scheme_name, @columns_used, @colomn_scale_factors);
 Function:  Creates a new energy scheme scenario in the given term of the term tracker structure.
            This function is just created for convenience since this action has to be taken by every
            function that creates an energy term for every energy scheme the term is used under.
 Returns :  A pointer to the newly inserted scheme.
 Args    :  (1) pointer to the term
            (2) name of the energy scheme (string - will be used as the hash key)
            (3) array of column names
            (4) array of column scale factors
=cut
sub addSchemeToTerm {
  my $term = shift;
  my $name = shift;
  my @arr = @_;
  if (scalar(@arr)%2) {
    die "Error in addSchemeToTerm: there must be as many column indices as scale factors:\n". join(" ", @arr) . "\n";
  }
  my @cols_used = @arr[0..scalar(@arr)/2-1];
  my @sf = @arr[scalar(@arr)/2..scalar(@arr)-1];

  $term->{'eschemes'}->{$name} = {};
  my $esch = $term->{'eschemes'}->{$name};
  $esch->{cols_used} = \@cols_used;
  $esch->{sf} = \@sf;

  return $esch;
}



=head2

 Title   : getAtomicSA
 Usage   : ENERGY::getAtomicSA($pdb);
 Function: Returns the solvent accessible areas of the atoms in the specified
           PDB structure. Uses naccess to do this.
 Returns : array of atomic accessible surface areas 
 Args    : (1) PDB structure
           (2) Probe radius
           (3) Optional: a reference to an array of atoms. If specified,
               only areas of these atoms will be returned. Otherwise, all
               atoms will be considered.
           (4) Optional: a double hash table of all atomic radii hashed by
               residue name first and atom name second. If this parameter is
               not passed and is needed, it is created by calling PDB::sizeLookup().
               It is needed when the previous parameter is passed (to determine
               distances between vdW surfaces of atoms). It is a good idea to
               pass this parameter if this function is called many times in a row
               to calculate the SAA of a subset of atoms in the structure
               (so that the size file is not read in and parsed every time).
           (5) Optional: a pointer to a hash table. If specified, it will be populated
               with accessible surface areas of atoms other than the ones
               considered according to parameter 3 - but only those, which
               occlude the atoms specified in parameter 3 from water. 
               The hashing is by atom string consistent with PDB::atomStr().
               If parameter 3 is not specified then all atoms are considered,
               in which case this hash table, naturally, will not be populated.
           (6) Optional: value for the z-slice parameter. The default is 0.1 A.

=cut

sub getAtomicSA {
  my $pdb = shift;
  my $pr = shift;
  GENERAL::requireArgs($pdb, $pr);

  # Just in case we are running in parallel, append rank number to all files
  my $cdir = GENERAL::GetDir();
  my $base = "atomic_sa_tmp_$PRANK_DEF";
  my $tpdbfile = "$base\.pdb";

  my $atomsp = shift;
  my $R = shift;
  my $Oa = shift; # reference to "occluding atoms" hash table
  if (defined($Oa) && (ref($Oa) !~ /HASH/)) {
    GENERAL::error("Parameter \"other atoms\" is not a valid hash pointer!");
  }
  my $zslice = shift; if (!defined($zslice)) { $zslice = 0.1; }
  my @atoms; # source atoms
  my @oatoms; # source atoms plus occluded atoms
  if (!defined($atomsp)) {
    # If need to do all atoms, then print all
    @atoms = PDB::conAtoms($pdb);
    $pdb->writePDB($tpdbfile, "renumber");
  } elsif (scalar(@$atomsp) == 0) {
    return ();
  } elsif (defined($Oa)) {
    # We'll need to know the vdW radii of all atoms
    if (!defined($R)) { $R = PDB::sizeLookup(); }
    
    # If need to return a subset as well as the rest, then print source atoms,
    # all atoms, which occlude source atoms and all atoms, which occlude the atoms occluding source atoms
    @atoms = @{$atomsp};
    
    @oatoms = @{$pdb->atomsWithin($atomsp, 2*$pr, $R)};
    my $patoms = $pdb->atomsWithin(\@oatoms, 2*$pr, $R);
    my $ofh = GENERAL::GetOutFH($tpdbfile);
    foreach my $a (@$patoms) { $ofh->print(PDB::_renumberpdbLine($a) . "\n"); }
    close($ofh);
  } else {
    # We'll need to know the vdW radii of all atoms
    if (!defined($R)) { $R = PDB::sizeLookup(); }
    
    # If we are only interested in a subset of source atoms, we need not consider all the atoms
    # in the molecule - just those, which occlude the source atoms from water.
    @atoms = @{$atomsp};

    # Print only those atoms in the structure less than two probe radii away from any source atom
    my $patoms = $pdb->atomsWithin($atomsp, 2*$pr, $R); # all atoms within 2 probe radii of the source atoms
    my $ofh = GENERAL::GetOutFH($tpdbfile);
    foreach my $a (@$patoms) { $ofh->print(PDB::_renumberpdbLine($a) . "\n"); }
    close($ofh);
  }

  # Run naccess
  my $vdwf = DEFINITIONS::getParmFile('naccess_vdw.dat');
  GENERAL::csystem("$NACCESS_DEF $tpdbfile -r $vdwf -a -y -p $pr -z $zslice > /dev/null");

  # Check for errors
  my $log_size = GENERAL::GetFileInfo("$base\.log", "lines");
  if ($log_size != 17) {
    GENERAL::error("There was a problem with the last naccess run in directory " . GENERAL::GetDir() . " (vdW radius file used $vdwf).");
  }

  # Read in results (parse output file)
  my @SA;
  my $ofh = GENERAL::GetInFH("$base\.asa");
  my $S = {}; # hash table of results
  while (<$ofh>) {
    next if ($_ !~ /^ATOM/);
    chomp($_);
    my @strs = split(" ", $_);
    $S->{"$strs[4]_$strs[3]_$strs[5]_$strs[2]"} = $strs[9];
  }
  close($ofh);

  # Find areas corresponding to requrested source atoms
  my $sent = "source";
  foreach my $atom (@atoms) {
    my $as = PDB::atomStr($atom);
    if (!defined($S->{$as})) {
      GENERAL::error("No area found for $as!");
    }
    push(@SA, $S->{$as});
    $S->{$as} = $sent;
  }
  # If necessary, find areas corresponding to occluding atoms
  if (defined($Oa)) {
    foreach my $atom (@oatoms) {
      my $as = PDB::atomStr($atom);
      if (!defined($S->{$as})) {
        GENERAL::error("No area found for occluding atom $as!");
      }
      unless ($S->{$as} eq $sent) { $Oa->{$as} = $S->{$as}; }
    }
  }

  # Clean up output
  GENERAL::crm("$base\.pdb", "$base\.log", "$base\.asa"); # remove all file

  return @SA;
}


=head2

 Title   : pepAtomicSolvation
 Usage   : ENERGY::pepAtomicSolvation($pdb, $eps_in, $output_format, $atoms);
 Function: Calculates the atomic solvation energies of the specified
           atoms in the specified structure using pep.
 Returns : Either a list (or hash) of energies or a list (or hash) of
           Born radii.
 Args    : (1) PDB of the structure
           (2) Internal dielectric constant
           (3) Optional: a reference to a double hash table of atomic radii.
               If this parameter is not specified, this table will be read in
               using PDB::sizeLookup().
           (4) Optional: a reference to an array of atoms of interest. If
               this parameter is not specified, all atoms will be considered.
           (5) Optional: output format. 1 means array of energies, 2 means hash of
               energies (by atom string), 3 means array of Born radii, 4
               means hash of Born radii. The default is 4.
           (6) Optional: version of pep to use (speed/accuracy trade-off).
               Currently available ones are 5, 10, and 21.
=cut

sub pepAtomicSolvation_old {
  my $pdb = shift;
  my $eps_in = shift;
  GENERAL::requireArgs($pdb, $eps_in);
  my $R = shift;
  if (!defined($R)) {
    $R = PDB::sizeLookup();
  }
  my $atoms = shift;
  if (!defined($atoms)) {
    $atoms = $pdb->conAtoms(undef, 1);
  }
  my $out = shift;
  if (!defined $out) { $out = 4; }
  if (($out != 1) && ($out != 2) && ($out != 3) && ($out != 4)) {
    GENERAL::error("Unknown output type $out!");
  }
  my $ver = shift;
  if (!defined($ver)) { $ver = 5; }
  if (($ver != 5) && ($ver != 10) && ($ver != 21)) {
    GENERAL::error("Unrecognized version of pep requested - $ver!");
  }
  my $conv_f = 332;

  # Since extensive network I/O seems to cause errors, I'll try to do much of I/O locally
#  my $cdir = GENERAL::GetDir();
#  GENERAL::cchdir($LOCAL_DEF);
  # Just in case we are running in parallel, append process rank to all file names
#  my $wdir = "pep_tmp_$PRANK_DEF";
#  GENERAL::cmkdir($wdir);
#  GENERAL::cchdir($wdir);
  my $mol = "mol.pqr_$PRANK_DEF";
  my $src = "src.pqr_$PRANK_DEF";
  my $surf = "mol.pt_$PRANK_DEF";
  my $alog = "asurf.log_$PRANK_DEF";
  my $plog = "pep.log_$PRANK_DEF";
  # Write molecule PQR file
  PDB::writePQR($pdb, $mol, $R);
  # Write sources PQR file
  PDB::writePQR($atoms, $src, $R);
  # Write asurf to get surface file
  GENERAL::csystem("echo -e \"$mol\n$surf\n1.4\nb\n\" | $ASURF_DEF > $alog");
  my $err = `grep -iE \"(error|warning)\" $alog`; chomp($err);
  if (!($err eq "")) {
    GENERAL::error("There was a problem with the previous asurf call in directory " . GENERAL::GetDir());
  }

  # PEP inputs
  #pep <<EOF
  #tst.pqr                 <- complete molecule
  #sources.pqr             <- list of atoms to calc. GF's for
  #NULL                    <- optional xyz file (see below)
  #tst.surf                <- output from asurf (a ".surf" file)
  #Pot/                    <- directory to write output data
  #1 20                    <- range of atoms in sources.pqr to calculate
  #1.4                     <- solvent radius (Ang)
  #2.0                     <- Stern layer thickness (Ang)
  #0.500                   <- Salt concentration (M)
  #4.0                     <- internal dielectric constant of solute
  #80.0                    <- solvent dielectric constant
  #2                       <- Boundary conditions for initial grid ( see below)
  #0                       <- center of first grid (see below)
  #1                       <- initial grid spacing is integer
  #4.0                     <- coarsest allowed grid spacing (Ang)
  #0.25                    <- focusing factor (see below)
  #0                       <- autospace option
  #4                       <- maximum number of focusing iterations
  #1.0e-5                  <- dlap matrix solver error tolerance
  #1                       <- reaction field flag
  #b                       <- output file format (b=binary a=ascii)
  #EOF

  # Call PEP to calculate Green functions in water and
  # parse through the output files to get energy in water
  my @enw;
  for (my $i=0; $i < scalar(@{$atoms}); $i++) {
#    print "Running pep on atom $i in water...\n";
    my $a = $atoms->[$i];
    if ($a->{charge} != 0) {
      my $j = $i+1;
      my $cmd;
      if ($ver == 5) {
        my $pepwin = "$mol\n$src\nNULL\n$surf\n./\n$j $j\n1.4\n2.0\n0.0\n$eps_in\n80.0\n2\n0\n1\n18.4\n0.25\n0\n20\n1.0e-5\n1\na\n";
        $cmd = "echo -e \"$pepwin\" | $PEP5_DEF > $plog";
      } elsif ($ver == 10) {
        my $pepwin = "$mol\n$src\nNULL\n$surf\n./\n$j $j\n1.4\n2.0\n0.0\n$eps_in\n80.0\n2\n0\n1\n9.2\n0.25\n0\n8\n1.0e-5\n1\na\n";
        $cmd = "echo -e \"$pepwin\" | $PEP11_DEF > $plog";
      } elsif ($ver == 21) {
        my $pepwin = "$mol\n$src\nNULL\n$surf\n./\n$j $j\n1.4\n2.0\n0.0\n$eps_in\n80.0\n2\n0\n1\n4.0\n0.25\n0\n8\n1.0e-5\n1\na\n";
        $cmd = "echo -e \"$pepwin\" | $PEP21_DEF > $plog";
      }
      GENERAL::csystem($cmd);
      my $err = `grep -iE \"(warning|error)\" $plog`; chomp($err);
      if (!($err eq "")) {
        GENERAL::error("There was a problem with the previous pep call in directory " . GENERAL::GetDir());
      }
      # Open output file corresponding the the atom
      my $ofile = sprintf("%-5s%3s%6d", $a->{atomname}, $a->{residue}->{resname}, $a->{residue}->{iresnum});
      $ofile =~ s/ /_/g;
      my $fh = GENERAL::GetInFH($ofile);
      # Read in the index of the atom
      <$fh>; my $ind = <$fh>; chomp($ind);
      # Skip the necessary number of lines to get to the potential at the atom
      for (my $j=1; $j <= $ind-1; $j++) { <$fh>; }
      # Get the potential at the atom
      my $pot = <$fh>; chomp($pot);
      # This potential is in electron/Angstrem under the assumption of 1e at the source atom
      if (($out == 3) || ($out == 4)) {
        $enw[$i] = $conv_f*$pot;
      } else {
        $enw[$i] = $conv_f*$pot*$a->{charge}*$a->{charge}*0.5;
      }
      close($fh);
      # Delete the output file (otherwise won't do the calculation for this atom again for gas)
      GENERAL::crm("$ofile");
    } else {
      $enw[$i] = 0;
    }
  }

  # Call PEP to calculate Green function in gas and
  # parse through the output file to get energy in gas
  my @sen;
  for (my $i=0; $i < scalar(@{$atoms}); $i++) {
#    print "Running pep on atom $i in gas...\n";
    my $a = $atoms->[$i];
    if ($a->{charge} != 0) {
      my $j = $i+1;
      my $cmd;
      if ($ver == 5) {
        my $pepgin = "$mol\n$src\nNULL\n$surf\n./\n$j $j\n1.4\n2.0\n0.0\n$eps_in\n$eps_in\n2\n0\n1\n18.4\n0.25\n0\n20\n1.0e-5\n1\na\n";
        $cmd = "echo -e \"$pepgin\" | $PEP5_DEF > $plog";
      } elsif ($ver == 10) {
        my $pepgin = "$mol\n$src\nNULL\n$surf\n./\n$j $j\n1.4\n2.0\n0.0\n$eps_in\n$eps_in\n2\n0\n1\n9.2\n0.25\n0\n8\n1.0e-5\n1\na\n";
        $cmd = "echo -e \"$pepgin\" | $PEP11_DEF > $plog";
      } elsif ($ver == 21) {
        my $pepgin = "$mol\n$src\nNULL\n$surf\n./\n$j $j\n1.4\n2.0\n0.0\n$eps_in\n$eps_in\n2\n0\n1\n4.0\n0.25\n0\n8\n1.0e-5\n1\na\n";
        $cmd = "echo -e \"$pepgin\" | $PEP21_DEF > $plog";
      }
      GENERAL::csystem($cmd);
      my $err = `grep -iE \"(warning|error)\" $plog`; chomp($err);
      if (!($err eq "")) {
        GENERAL::error("There was a problem with the previous pep call in direcotry " . GENERAL::GetDir());
      }
      my $ofile = sprintf("%-5s%3s%6d", $a->{atomname}, $a->{residue}->{resname}, $a->{residue}->{iresnum});
      $ofile =~ s/ /_/g;
      my $fh = GENERAL::GetInFH($ofile);
      # Read in the index of the atom
      <$fh>; my $ind = <$fh>; chomp($ind);
      # Skip the necessary number of lines to get to the potential at the atom
      for (my $j=1; $j <= $ind-1; $j++) { <$fh>; }
      # Get the potential at the atom
      my $pot = <$fh>; chomp($pot);
      # This potential is in electron/Angstrem under the assumption of 1e at the source atom
      if (($out == 3) || ($out == 4)) {
        $sen[$i] = $enw[$i] - $conv_f*$pot;
      } else {
        $sen[$i] = $enw[$i] - $conv_f*$pot*$a->{charge}*$a->{charge}*0.5;
      }
      close($fh);
      # Delete the output file (otherwise won't do the calculation for this atom again for gas)
      GENERAL::crm("$ofile");
    } else {
      $sen[$i] = 0;
    }
  }
  GENERAL::crm($mol, $src, $surf, $alog, $plog); # remove all files
#  GENERAL::crm(glob("*")); # remove all files
#  GENERAL::cchdir("..");
#  GENERAL::crmdir("$wdir");
#  GENERAL::cchdir($cdir); # go back to global directory

  # Return solvation energy
  if ($out == 1) {
    return @sen;
  } elsif ($out == 2) {
    my %sen;
    my $i=0;
    foreach my $a (@{$atoms}) {
      my $astr = PDB::atomStr($a);
      $sen{$astr} = $sen[$i];
      $i++;
    }
    return %sen;
  } elsif ($out == 3) {
    my @br;
    my $i=0;
    foreach my $a (@{$atoms}) {
      if ($a->{charge} != 0) {
        $br[$i] = -332*(1/$eps_in-1/80)/$sen[$i];
      } else {
        $br[$i] = 1;
      }
      $i++;
    }
    return @br;
  } elsif ($out == 4) {
    my %br;
    my $i=0;
    foreach my $a (@{$atoms}) {
      my $astr = PDB::atomStr($a);
      if ($a->{charge} != 0) {
        $br{$astr} = -332*(1/$eps_in-1/80)/$sen[$i];
      } else {
        $br{$astr} = 1;
      }
      $i++;
    }
    return %br;
  }
}

=head2

 Title   : pepAtomicSolvation
 Usage   : ENERGY::pepAtomicSolvation($pdb, $eps_in, $output_format, $atoms);
 Function: Calculates the atomic solvation energies of the specified
           atoms in the specified structure using pep. Calculates as many atoms
           as possible per each pep run to save time.
 Returns : Either a list (or hash) of energies or a list (or hash) of
           Born radii.
 Args    : (1) PDB of the structure
           (2) Internal dielectric constant
           (3) Optional: a reference to a double hash table of atomic radii.
               If this parameter is not specified, this table will be read in
               using PDB::sizeLookup().
           (4) Optional: a reference to an array of atoms of interest. If
               this parameter is not specified, all atoms will be considered.
           (5) Optional: output format. 1 means array of energies, 2 means hash of
               energies (by atom string), 3 means array of Born radii, 4
               means hash of Born radii. The default is 4.
           (6) Optional: version of pep to use (speed/accuracy trade-off).
               Currently available ones are 5, 10, and 21.
=cut

sub pepAtomicSolvation {
  my $pdb = shift;
  my $eps_in = shift;
  GENERAL::requireArgs($pdb, $eps_in);
  my $R = shift;
  if (!defined($R)) {
    $R = PDB::sizeLookup();
  }
  my $atoms = shift; # array of atoms the potential at which is of interest
  if (!defined($atoms)) {
    $atoms = $pdb->conAtoms(undef, 1);
  }
  my @catoms; # charged and neutral atoms
  foreach my $a (@{$atoms}) {
    if ($a->{charge} != 0) { push(@catoms, $a); }
  }
  my $out = shift;
  if (!defined $out) { $out = 4; }
  if (($out != 1) && ($out != 2) && ($out != 3) && ($out != 4)) {
    GENERAL::error("Unknown output type $out!");
  }
  my $ver = shift;
  if (!defined($ver)) { $ver = 5; }
  # Version-dependent parameters
  my ($parw, $parg, $gs, $pepbin);
  if ($ver == 5) {
    $gs = 5; # grid size (hard coded in pep version) - grid extends +/- grid size
    my $ms = 18.4; # default maximum grid spacing
    my $nfi = 20; # number of focusing iteractions
    # make sure the molecule fits entirely in the coarsest grid (adjust grid spacing if not)
    # the grid extends over a total length of 2*$gs (from -$sg to +$gs), but adjustments of up to 1 grid point 
    # may be done by the program to place the central point onto a grid point, so we have to conservatively
    # estimate the grid's size as 2*$gs-1 (read PEP documentation about centering around the geometric center)
    my @box = PDB::box($atoms, $R);
    if ((2*$gs-1)*$ms < GENERAL::max($box[1]-$box[0], $box[3]-$box[2], $box[5]-$box[4])) {
      my $nms = GENERAL::max($box[1]-$box[0], $box[3]-$box[2], $box[5]-$box[4])/(2*$gs-1);
#printf(GENERAL::max($box[1]-$box[0], $box[3]-$box[2], $box[5]-$box[4]) . " changing $ms x $nfi to ");
      $nfi = GENERAL::ceil($nfi*log($ms)/log($nms));
      $ms = $nms;
#printf("$ms x $nfi\n");
    }
    $parw = "1.4\n2.0\n0.0\n$eps_in\n80.0\n2\n0\n1\n$ms\n0.25\n0\n$nfi\n1.0e-5\n1\na\n";
    $parg = "1.4\n2.0\n0.0\n$eps_in\n$eps_in\n2\n0\n1\n$ms\n0.25\n0\n$nfi\n1.0e-5\n1\na\n";
    $pepbin = $PEP5_DEF;
  } elsif ($ver == 10) {
    $gs = 11; # grid size (hard coded in pep version) - grid extends +/- grid size
    my $ms = 9.2; # default maximum grid spacing
    my $nfi = 8; # number of focusing iteractions
    my @box = PDB::box($atoms);
    if ((2*$gs-1)*$ms < GENERAL::max($box[1]-$box[0], $box[3]-$box[2], $box[5]-$box[4])) {
      my $nms = GENERAL::max($box[1]-$box[0], $box[3]-$box[2], $box[5]-$box[4])/(2*$gs-1);
      $nfi = GENERAL::cei($nfi*log($ms)/log($nms));
      $ms = $nms;
    }
    $parw = "1.4\n2.0\n0.0\n$eps_in\n80.0\n2\n0\n1\n$ms\n0.25\n0\n$nfi\n1.0e-5\n1\na\n";
    $parg = "1.4\n2.0\n0.0\n$eps_in\n$eps_in\n2\n0\n1\n$ms\n0.25\n0\n$nfi\n1.0e-5\n1\na\n";
    $pepbin = $PEP11_DEF;
  } elsif ($ver == 21) {
    $gs = 21; # grid size (hard coded in pep version) - grid extends +/- grid size
    my $ms = 4; # default maximum grid spacing
    my $nfi = 8; # number of focusing iteractions
    my @box = PDB::box($atoms);
    if ((2*$gs-1)*$ms < GENERAL::max($box[1]-$box[0], $box[3]-$box[2], $box[5]-$box[4])) {
      my $nms = GENERAL::max($box[1]-$box[0], $box[3]-$box[2], $box[5]-$box[4])/(2*$gs-1);
      $nfi = GENERAL::ceil($nfi*log($ms)/log($nms));
      $ms = $nms;
    }
    $parw = "1.4\n2.0\n0.0\n$eps_in\n80.0\n2\n0\n1\n$ms\n0.25\n0\n$nfi\n1.0e-5\n1\na\n";
    $parg = "1.4\n2.0\n0.0\n$eps_in\n$eps_in\n2\n0\n1\n$ms\n0.25\n0\n$nfi\n1.0e-5\n1\na\n";
    $pepbin = $PEP21_DEF;
  } else {
    GENERAL::error("Unrecognized version of pep requested - $ver!");
  }
  my $conv_f = 332;

  my $mol = "mol.pqr_$PRANK_DEF";
  my $src = "src.pqr_$PRANK_DEF";
  my $surf = "mol.pt_$PRANK_DEF";
  my $alog = "asurf.log_$PRANK_DEF";
  my $plog = "pep.log_$PRANK_DEF";
  if (scalar(@catoms) != 0) {
    # Note, when writing PQR mol and src files for PEP make sure NOT to renumber the structures - i.e. the atom numbers in the
    # src file MUST be the same as in the mol file (PEP looks for source atom coordinates in the mol file given the atom string derived
    # from the source file, this atom strings must match for PEP to charge the right atom.
    # Write molecule PQR file
    PDB::writePQR($pdb, $mol, $R, 0);
    # Write sources PQR file
    PDB::writePQR(\@catoms, $src, $R, 0); # write only charged atoms
    # Write asurf to get surface file
    GENERAL::csystem("echo -e \"$mol\n$surf\n1.4\nb\n\" | $ASURF_DEF > $alog");
    my $err = `grep -iE \"(error|warning)\" $alog`; chomp($err);
    if (!($err eq "")) {
      GENERAL::error("There was a problem with the previous asurf call in directory " . GENERAL::GetDir());
    }
  }

  # PEP inputs
  #pep <<EOF
  #tst.pqr                 <- complete molecule
  #sources.pqr             <- list of atoms to calc. GF's for
  #NULL                    <- optional xyz file (see below)
  #tst.surf                <- output from asurf (a ".surf" file)
  #Pot/                    <- directory to write output data
  #1 20                    <- range of atoms in sources.pqr to calculate
  #1.4                     <- solvent radius (Ang)
  #2.0                     <- Stern layer thickness (Ang)
  #0.500                   <- Salt concentration (M)
  #4.0                     <- internal dielectric constant of solute
  #80.0                    <- solvent dielectric constant
  #2                       <- Boundary conditions for initial grid ( see below)
  #0                       <- center of first grid (see below)
  #1                       <- initial grid spacing is integer
  #4.0                     <- coarsest allowed grid spacing (Ang)
  #0.25                    <- focusing factor (see below)
  #0                       <- autospace option
  #4                       <- maximum number of focusing iterations
  #1.0e-5                  <- dlap matrix solver error tolerance
  #1                       <- reaction field flag
  #b                       <- output file format (b=binary a=ascii)
  #EOF

  # Call PEP to calculate Green functions in water
  my @enw;
  if (scalar(@catoms) != 0) {
    my $cmd = "$mol\n$src\nNULL\n$surf\n./\n1 " . scalar(@catoms) . "\n$parw";
    GENERAL::csystem("echo -e \"$cmd\" | $pepbin > $plog");
    my $err = `grep -iE \"(warning|error)\" $plog`; chomp($err);
    if (!($err eq "")) {
      GENERAL::error("There was a problem with the previous pep call in directory " . GENERAL::GetDir());
    }
  }
  # Parse through the output files to get energy in water
  for (my $i = 0; $i < scalar(@{$atoms}); $i++) {
    my $a = $atoms->[$i];
    next if ($a->{charge} == 0);
    # Open output file corresponding the the atom
    my $ofile = sprintf("%-5s%3s%6d\_%d", $a->{atomname}, $a->{residue}->{resname}, $a->{residue}->{iresnum}, $a->{atomunx});
    $ofile =~ s/ /_/g;
    my $arr = GENERAL::file2array($ofile);
    my $pot = $arr->[$arr->[1]+1] + 0.0;
    # This potential is in electron/Angstrem under the assumption of 1e at the source atom
    if (($out == 3) || ($out == 4)) {
      $enw[$i] = $conv_f*$pot;
    } else {
      $enw[$i] = $conv_f*$pot*$a->{charge}*$a->{charge}*0.5;
    }
    # Delete the output file (otherwise won't do the calculation for this atom again for gas)
    GENERAL::crm("$ofile");
  }

  # Call PEP to calculate Green function in gas
  my @sen; # solvation energies
  if (scalar(@catoms) != 0) {
    my $cmd = "$mol\n$src\nNULL\n$surf\n./\n1 " . scalar(@catoms) . "\n$parg";
    GENERAL::csystem("echo -e \"$cmd\" | $pepbin > $plog");
    my $err = `grep -iE \"(warning|error)\" $plog`; chomp($err);
    if (!($err eq "")) {
      GENERAL::error("There was a problem with the previous pep call in direcotry " . GENERAL::GetDir());
    }
  }
  # Parse through the output file to get energy in gas
  for (my $i=0; $i < scalar(@{$atoms}); $i++) {
    my $a = $atoms->[$i];
    if ($a->{charge} == 0) { $sen[$i] = 0; next; }
    my $ofile = sprintf("%-5s%3s%6d\_%d", $a->{atomname}, $a->{residue}->{resname}, $a->{residue}->{iresnum}, $a->{atomunx});
    $ofile =~ s/ /_/g;
    my $arr = GENERAL::file2array($ofile);
    my $pot = $arr->[$arr->[1]+1] + 0.0;
    if (($out == 3) || ($out == 4)) {
      $sen[$i] = $enw[$i] - $conv_f*$pot;
    } else {
      $sen[$i] = $enw[$i] - $conv_f*$pot*$a->{charge}*$a->{charge}*0.5;
    }
    # Delete the output file (otherwise won't do the calculation for this atom again for gas)
    GENERAL::crm("$ofile");
  }
  if (scalar(@catoms) != 0) { GENERAL::crm($mol, $src, $surf, $alog, $plog); } # remove all files

  # Return solvation energy
  if ($out == 1) {
    return @sen;
  } elsif ($out == 2) {
    my %sen;
    my $i = -1;
    foreach my $a (@{$atoms}) {
      $i++;
      my $astr = PDB::atomStr($a);
      $sen{$astr} = $sen[$i];
    }
    return %sen;
  } elsif ($out == 3) {
    my @br;
    my $i = -1;
    foreach my $a (@{$atoms}) {
      $i++;
      if ($a->{charge} == 0) { $br[$i] = -1; next; }
      $br[$i] = -332*(1/$eps_in-1/80)/$sen[$i];
    }
    return @br;
  } elsif ($out == 4) {
    my %br;
    my $i = -1;
    foreach my $a (@{$atoms}) {
      $i++;
      my $astr = PDB::atomStr($a);
      if ($a->{charge} == 0) { $br{$astr} = -1; next; }
      $br{$astr} = -332*(1/$eps_in-1/80)/$sen[$i];
    }
    return %br;
  }
}

=head2

 Title   :  getGBSelf
 Usage   :  my $Gpol = ENERGY::getGBSelf(\%born_radii, $eps_in, \@atoms);
 Function:  Calculates the sum of atomic polarization energies of the given atomic collection.
 Returns :  Energy (scalar)
 Args    :  (1) Reference to a hash table of atomic Born radii (hashed by atom string)
            (2) Internal dielectric constant.
            (3) Reference to an array of atoms that comprize the structure.
=cut

sub getGBSelf {
  my $br = shift;
  my $eps_in = shift;
  my $a = shift;
  GENERAL::requireArgs($br, $eps_in, $a);

  my $Gpol = 0;
  for (my $i = 0; $i < scalar(@$a); $i++) {
    $Gpol += -166*(1/$eps_in - 1/80)*($a->[$i]->{charge})**2/$br->{PDB::atomStr($a->[$i])} if ($a->[$i]->{charge} != 0);
  }
  return $Gpol;
}


=head2

 Title   :  getGBInterAA
 Usage   :  my $Gpol = ENERGY::getGBInterAA(\%born_radii, $eps_in, \@atoms1, \@atoms1);
 Function:  Calculates the polarization energy associated with screening of interactions 
            between atoms within the specified group groups.
 Returns :  Energy (scalar)
 Args    :  (1) Reference to a hash table of atomic Born radii (hashed by atom string)
            (2) Internal dielectric constant.
            (3) Reference to an array of atoms.
            (4) Optional: cutoff distance for an interaction.
=cut

sub getGBInterAA {
  my $br = shift;
  my $eps_in = shift;
  my $a = shift;
  GENERAL::requireArgs($br, $eps_in, $a);
  my $cutoff = shift;

  my $Gpol = 0; my $d = 0;
  for (my $i = 0; $i < scalar(@$a) - 1; $i++) {
    next if ($a->[$i]->{charge} == 0);
    for (my $j = $i+1; $j < scalar(@$a); $j++) {
      $d = PDB::atomDist($a->[$i], $a->[$j]);
      next if (defined($cutoff) && ($d > $cutoff));
      $Gpol += -332*(1/$eps_in - 1/80)*($a->[$i]->{charge})*($a->[$j]->{charge})/fGB($br->{PDB::atomStr($a->[$i])}, $br->{PDB::atomStr($a->[$j])}, $d) if ($a->[$j]->{charge} != 0);
    }
  }
  return $Gpol;
}


=head2

 Title   :  getGBInterAB
 Usage   :  my $Gpol = ENERGY::getGBInterAB(\%born_radii, $eps_in, \@atoms);
 Function:  Calculates the polarization energy associated with screening of interactions 
            between the two (separate) specified atom groups.
 Returns :  Energy (scalar)
 Args    :  (1) Reference to a hash table of atomic Born radii (hashed by atom string)
            (2) Internal dielectric constant.
            (3) Reference to an array of atoms in the first group.
            (3) Reference to an array of atoms in the second group.
            (4) Optional: cutoff distance for an interaction.
=cut

sub getGBInterAB {
  my $br = shift;
  my $eps_in = shift;
  my $atoms1 = shift;
  my $atoms2 = shift;
  GENERAL::requireArgs($br, $eps_in, $atoms1, $atoms2);
  my $cutoff = shift;

  my $Gpol = 0; my $d = 0;
  foreach my $a1 (@$atoms1) {
    next if ($a1->{charge} == 0);
    foreach my $a2 (@$atoms2) {
      $d = PDB::atomDist($a1, $a2);
      next if (defined($cutoff) && ($d > $cutoff));
      $Gpol += -332*(1/$eps_in - 1/80)*($a1->{charge})*($a2->{charge})/fGB($br->{PDB::atomStr($a1)}, $br->{PDB::atomStr($a2)}, $d) if ($a2->{charge} != 0);
    }
  }
  return $Gpol;
}


=head2

 Title   :  fGB
 Usage   :  my $fgb = ENERGY::fGB($born_radius1, $born_radius2, $distance);
 Function:  For use with the Generalized Born equation.
 Returns :  fgb (scalar)
 Args    :  1. Born radius of atom 1
            2. Born radius of atom 2
            3. Distance between atoms 1 and 2
=cut

sub fGB {
  my $ai = shift;
  my $aj = shift;
  my $rij = shift;
  return sqrt($rij**2 + $ai*$aj*exp(-$rij**2/(4*$ai*$aj)));
}


=head2

 Title   :  getCoulombTotal
 Usage   :  my $G = ENERGY::getCoulombTotal($eps_in, \@atoms);
            my $G = ENERGY::getCoulombTotal($eps_in, \@atoms, \@target_atoms);
 Function:  Calculates the total Coulombic energy of the given atomic collection in
            a uniform dielectric medium.
 Returns :  Energy (scalar)
 Args    :  (1) Dielectric constant.
            (2) Reference to an array of atoms that comprize the structure.
            (3) Optional: cutoff distance for an interaction.
=cut

sub getCoulombTotal {
  my $eps = shift;
  my $atoms = shift;
  GENERAL::requireArgs($eps, $atoms);
  my $cutoff = shift;

  my $Ge = 0; my $d = 0;
  for (my $i = 0; $i < scalar(@$atoms)-1; $i++) {
    next if ($atoms->[$i]->{charge} == 0);
    for (my $j = $i+1; $j < scalar(@$atoms); $j++) {
      $d = PDB::atomDist($atoms->[$i], $atoms->[$j]);
      next if (defined($cutoff) && ($d > $cutoff));
      $Ge += 332*(1/$eps)*($atoms->[$i]->{charge})*($atoms->[$j]->{charge})/$d if ($atoms->[$j]->{charge} != 0);
    }
  }
  return $Ge;
}

=head2

 Title   :  getCoulombInter
 Usage   :  my $G = ENERGY::getCoulombInter($eps_in, \@atoms1, \@atoms2);
 Function:  Calculates the Coulombic energy of interaction between two groups
            of atoms.
 Returns :  Energy (scalar)
 Args    :  (1) Dielectric constant.
            (2) Reference to an array of atoms in the first group.
            (3) Reference to an array of atoms in the second group.
            (4) Optional: cutoff distance for an interaction.
=cut

sub getCoulombInter {
  my $eps = shift;
  my $atoms1 = shift;
  my $atoms2 = shift;
  GENERAL::requireArgs($eps, $atoms1, $atoms2);
  my $cutoff = shift;

  my $Ge = 0; my $d = 0;
  foreach my $a1 (@$atoms1) {
    next if ($a1->{charge} == 0);
    foreach my $a2 (@$atoms2) {
      $d = PDB::atomDist($a1, $a2);
      next if (defined($cutoff) && ($d > $cutoff));
      $Ge += (1/$eps)*($a1->{charge})*($a2->{charge})*332/$d if ($a2->{charge} != 0);
    }
  }
  return $Ge;
}


=head2

 Title   :  runDelphiSolvation
 Usage   :  my $G = ENERGY::runDelphiSolvation($pdb, ...);
 Function:  Runs delphi to find the polarization energy and returns the output (contents of fort.16).
 Returns :  Either the total energy (kcal/mol) of the specified atoms or an array of potentials (kt/e)
            at the specified atoms.
 Args    :  (1) pdb struct of the molecule.
            (2) internal dielectric constant
            (3) Pointer to array of atoms to charge.
            (4) Pointer to array of atoms to output potential at.
            (5) Optinal: grid size (integer). Default is 65.
            (6) optional: flag specifying whether the total energy or the array of potentials at
                atoms specified in parameter 2 is to be returned. 1 means total energy, 2 means
                array of atomic potentials. 3 means that a hash of atomic potentials will be returned
                (hash keys correspond to PDB::atomStr). 1 is the default.
            (7) optional: flag specifying whether the 2 stage mechanism for calculating solvation
                energy is used (water - gas) or the direct method (asking Delphi to output the reaciton
                field energy). 1 means the direct method, 2 means 2 stage method (2 is default).

=cut

sub runDelphiSolvation {
  my $pdb = shift;
  my $eps_in = shift;
  my $catoms = shift; # charged atoms
  my $satoms = shift; # site atoms (only potential in these atoms will be output)
  GENERAL::requireArgs($pdb, $eps_in, $catoms, $satoms);
  my $gsize = shift; if (!defined($gsize)) { $gsize = 65; }
  my $out = shift;
  if (!defined $out) { $out = 1; }
  my $method = shift;
  if (!defined $method) { $method = 2; }
  if (($method == 1) && ($out == 2)) {
    GENERAL::error("Output $out and method $method not compatible!");
  }
  my (@arr, @arr1, $en, $ans);
  my $conv_f = 0.592; # kt/e -> kcal/(mol*e) conversion factor

  $pdb->PDB::writePDB("fort.13", "");
  GENERAL::csystem("cp " . DEFINITIONS::getParmFile("fort.11") . " .");

  # create charge file
  my $fh = GENERAL::GetOutFH("fort.12");
  print $fh "atom__resnumbc_charge_\n";
  foreach my $a (@{$catoms}) {
    printf($fh "%-4s%5s%4d%1s%8.5f\n", $a->{atomname}, $a->{residue}->{resname}, $a->{residue}->{resnum}, $a->{residue}->{chain}->{id}, $a->{charge});
  }
  close($fh);

  # create site file
  $fh = GENERAL::GetOutFH("fort.15");
  foreach my $a (@{$satoms}) {
    # we can't charge site atoms, so temporarily set charge to zero
    my $crg = $a->{charge};
    $a->{charge} = 0;
    print $fh PDB::_pdbLine($a) . "\n";
    $a->{charge} = $crg;
  }
  close($fh);

  if ($method == 1) {
    GENERAL::csystem("cat " . DEFINITIONS::getParmFile("waterlow_sol.prm") . " | grep -Ev \"gsize\" > fort.10");
    GENERAL::csystem("echo \"gsize=$gsize\nindi=$eps_in\" >> fort.10");
    GENERAL::csystem("/mit/gevorg/bin/gdelphi > /dev/null");
    GENERAL::csystem("cat " . DEFINITIONS::getParmFile("waterhigh_sol.prm") . " | grep -Ev \"gsize\" > fort.10");
    GENERAL::csystem("echo \"gsize=$gsize\nindi=$eps_in\" >> fort.10");
    GENERAL::csystem("/mit/gevorg/bin/gdelphi > /dev/null");
    $ans = `grep -E "^corrected reaction field" energy.dat`;
    $ans =~ /corrected reaction field\s+en\s+(\S+)\s+kt/;
    $en = $1*$conv_f;
    GENERAL::crm("energy.dat");
  } elsif ($method == 2) {
    # run delphi to get solvation energy
    GENERAL::csystem("cat " . DEFINITIONS::getParmFile("gaslow_t.prm") . " | grep -Ev \"gsize\" > fort.10");
    GENERAL::csystem(sprintf("echo \"gsize=$gsize\nindi=%f\nexdi=$eps_in\" >> fort.10", $eps_in + 0.000001));
    GENERAL::csystem("$DELPHI_DEF > /dev/null");
    GENERAL::csystem("cat " . DEFINITIONS::getParmFile("gashigh_t.prm") . " | grep -Ev \"gsize\" > fort.10");
    GENERAL::csystem(sprintf("echo \"gsize=$gsize\nindi=%f\nexdi=$eps_in\" >> fort.10", $eps_in + 0.000001));
    GENERAL::csystem("$DELPHI_DEF > /dev/null");
    if ($out == 1) {
      $ans = `grep -E "^ total energy" fort.16`;
      $ans =~ /total energy =\s+(\S+)\s+kt/;
      $en = $1*$conv_f;
    } elsif (($out == 2) || ($out == 3)) {
      @arr = ENERGY::readDelphiFRC("fort.16", 5);
    }

    GENERAL::csystem("cat " . DEFINITIONS::getParmFile("waterlow_t.prm") . " | grep -Ev \"gsize\" > fort.10");
    GENERAL::csystem("echo \"gsize=$gsize\nindi=$eps_in\" >> fort.10");
    GENERAL::csystem("$DELPHI_DEF > /dev/null");
    GENERAL::csystem("cat " . DEFINITIONS::getParmFile("waterhigh_t.prm") . " | grep -Ev \"gsize\" > fort.10");
    GENERAL::csystem("echo \"gsize=$gsize\nindi=$eps_in\" >> fort.10");
    GENERAL::csystem("$DELPHI_DEF > /dev/null");
    if ($out == 1) {
      $ans = `grep -E "^ total energy" fort.16`;
      $ans =~ /total energy =\s+(\S+)\s+kt/;
      $en = $1*$conv_f - $en;
    } elsif (($out == 2) || ($out == 3)) {
      @arr1 = ENERGY::readDelphiFRC("fort.16", 5);
      for (my $i=0; $i < scalar(@arr); $i++) {
        $arr[$i] = ($arr1[$i] - $arr[$i])*$conv_f;
      }
    }
    GENERAL::crm("fort.16");
  } else {
    die "Error in RunDelphISolvation: unknown method $method!\n";
  }
  GENERAL::crm("fort.10", "fort.11", "fort.12", "fort.13", "fort.15");
  if (-e "ARCDAT") { GENERAL::crm("ARCDAT"); }
  if (-e "temp.phi") { GENERAL::crm("temp.phi"); }
  if (-e "perNew") { GENERAL::crm("perNew"); }
  if ($out == 1) {
    return $en;
  } elsif ($out == 2) {
    return @arr;
  } elsif ($out == 3) {
    my %H;
    for (my $i = 0; $i < scalar(@$satoms); $i++) {
      $H{PDB::atomStr($satoms->[$i])} = $arr[$i];
    }
    return \%H;
  } else {
    GENERAL::error("Output type $out unknown");
  }
}


=head2

 Title   :  readDelphiFRC
 Usage   :  ENERGY::readDelphiFRC($frc_file, $col_number);
 Function:  Returns an array containing the specified column of the specified
            delphi output formatted FRC file.
 Returns :  Array.
 Args    :  (1) FRC file name
            (2) column number to return

=cut

sub readDelphiFRC {
  my $frc = shift;
  my $col = shift;

  my $frch = GENERAL::GetInFH($frc);
  my $i = 0;
  my @col;
  while (<$frch>) {
    $i++;
    # skip formatted FRC header
    next if ($i <= 12);
    last if (/total energy/);

    # DelPhi writes all entries of the formatter FRC file in f10.4 format (see src/wrtsit4.f line 626-724)
    push(@col, GENERAL::Trim(substr($_, 10*($col-1), 10)));
  }
  close($frch);
  return @col;
}


=head2

 Title   :  interfaceEPIC
 Usage   :  ENERGY::interfaceEPIC($ofile, $terms, 'self'/'pair', \%T, \%P, \%parm, \@ai, $inpfi, [\@aj, $inpfj], $ex12, $ex13, $ex14);
 Function:  Interfaces with EPIC to calculate the necessary term and write it to the specified file.
 Returns :  Nothing
 Args    :  (0)  Output file name
            (1)  Names of terms to calculate
            (2)  Type of calculation (self or pair). In case of 'self', expects one atom array and one file with coordinates.
                 With pair, expects two array atoms and two coordinate files.
            (3)  Topology hash reference.
            (4)  Atomic parameter hash reference.
            (5)  Calculation parameter hash reference (things like nonbond cutoffs, temperature, dielectric constant - whatever is
                 needed for the calculation. Will through error if a necessary value is missing in the hash.
            (6)  First atom array reference.
            (7)  Coordinate file name for the first set of atoms.
            (8)  Optional: second atom array reference (only for 'pair')
            (9)  Optional: coordinate file name for the second set of atoms (only for 'pair')
            (10) Optional: 1-2 exclusion list
            (11) Optional: 1-3 exclusion list
            (12) Optional: 1-4 exclusion list

=cut

sub interfaceEPIC {
  my $ofile = shift;
  my $terms = GENERAL::Trim(shift);
  my $ps = uc(shift);
  my $T = shift;
  my $P = shift;
  my $parm = shift;
  my $atomsi = shift;
  my $inpi = shift;
  GENERAL::requireArgs($ofile, $terms, $ps, $T, $P, $parm, $atomsi, $inpi);
  if (($ps ne "PAIR") && ($ps ne "SELF") && ($ps ne "CORR")) { GENERAL::error("Expected parameter ps to be either 'PAIR' or 'SELF', while it is '$ps'."); }
  my %par = @_;
  my $atomsj = defined($par{atomsj}) ? $par{atomsj} : undef;
  my $inpj = defined($par{inpj}) ? $par{inpj} : undef;
  if ((($ps eq "PAIR") || ($ps eq "CORR")) && ((!defined($par{inpj})) || (!defined($par{atomsj})))) { GENERAL::error("Two atom sets required with 'PAIR' or 'CORR' options.") }
  my $sf = defined($par{sf}) ? $par{sf} : 1.0;
  my $ex12 = defined($par{ex12}) ? $par{ex12} : undef;
  my $ex13 = defined($par{ex13}) ? $par{ex13} : undef;
  my $ex14 = defined($par{ex14}) ? $par{ex14} : undef;
  my $prog = defined($par{prog}) ? $par{prog} : undef;

  # Since I am going to be calling epic.exe a lot, copy it locally the first time
  GENERAL::csystem("cp $BIN_DEF/epic.exe $LOCAL_DEF/epic.exe") if (! -e "$LOCAL_DEF/epic.exe");

  # Call EPIC
  if (!defined($prog)) {
    $prog = new IO::File;
    $prog->open("| $LOCAL_DEF/epic.exe");
  }
  if ($ps eq "SELF") {
    $prog->printf("s $ofile $terms $inpi $sf\n");
  } elsif ($ps eq "PAIR") {
    $prog->printf("p $ofile $terms $inpi $inpj $sf\n");
  } else {
    $prog->printf("c $ofile $terms $inpi $inpj $sf\n");
  }

  my @terms = split(" ", $terms);
  foreach my $term (@terms) {
    if ($term =~ /EEF/) {
      # Certain calculation parameters are expected
      if (!defined($parm->{eef_cutnb})) { GENERAL::error("Expected 'eef_cutnb' in parameter hash for EEF1 calculation"); }
      if (!defined($parm->{eef_temp})) { GENERAL::error("Expected 'eef_temp' in parameter hash for EEF1 calculation"); }
      # Prepare parameter string
      my $param = "";
      my @aas; push(@aas, $atomsi); if (defined($atomsj)) { push(@aas, $atomsj); }
      foreach my $aa (@aas) {
        my $gi = 1; my %G;
        foreach my $ai (@$aa) {
          # fish out atom type and group info
          my $rni = $ai->{residue}->{resname};
          if (defined($ai->{residue}->{patch})) { # see if residue is patched
            $rni .= "_patch_" . $ai->{residue}->{patch};
            PDB::createPatchedResidue($T, $ai->{residue}->{resname}, $ai->{residue}->{patch}, $rni, 0);
          }
          my $rs = PDB::resStr($ai->{residue});
          my $ani = $ai->{atomname};
          my $ati = $T->{residues}->{$rni}->{atoms}->{$ani}->{type};
          if (!defined($ati)) { GENERAL::error("No atom type defined for atom $ani in residue $rni"); }
          my $gn = $T->{residues}->{$rni}->{atoms}->{$ani}->{group};
          if (!defined($gn)) { GENERAL::error("No group defined for atom $ani in residue $rni"); }
          # make sure groups are numbered 1 ... N
          if (!defined($G{"$rs$gn"})) { $G{"$rs$gn"} = $gi; $gi++; }
          # print parameters
          if (defined($P->{eef1}{$ati}{V})) { $param .= sprintf("%f", $P->{eef1}{$ati}{V}); }
          else { GENERAL::error("EEF1 parameter V not defined for atom type $ati"); }
          if (defined($P->{eef1}{$ati}{Gref})) { $param .= sprintf(" %f", $P->{eef1}{$ati}{Gref}); }
          else { GENERAL::error("EEF1 parameter Gref not defined for atom type $ati"); }
          if (defined($P->{eef1}{$ati}{Gfree})) { $param .= sprintf(" %f", $P->{eef1}{$ati}{Gfree}); }
          else { GENERAL::error("EEF1 parameter Gfree not defined for atom type $ati"); }
          if (defined($P->{eef1}{$ati}{Href})) { $param .= sprintf(" %f", $P->{eef1}{$ati}{Href}); }
          else { GENERAL::error("EEF1 parameter Href not defined for atom type $ati"); }
          if (defined($P->{eef1}{$ati}{CPref})) { $param .= sprintf(" %f", $P->{eef1}{$ati}{CPref}); }
          else { GENERAL::error("EEF1 parameter CPref not defined for atom type $ati"); }
          if (defined($P->{eef1}{$ati}{SigW})) { $param .= sprintf(" %f", $P->{eef1}{$ati}{SigW}); }
          else { GENERAL::error("EEF1 parameter SigW not defined for atom type $ati"); }
          $param .= sprintf(" %f", $G{"$rs$gn"});
          if (defined($P->{vdW}{$ati}{rmin})) { $param .= sprintf(" %f\n", $P->{vdW}{$ati}{rmin}); }
          else { GENERAL::error("van der Waals radius not defined for atom type $ati"); }
        }
      }

      if ($parm->{potm} eq "linrep") {
        $prog->printf("$parm->{eef_cutnb} $parm->{eef_temp} $parm->{potm}\n");
      } else {
        $prog->printf("$parm->{eef_cutnb} $parm->{eef_temp}\n");
      }
      $prog->printf("exclude\n");
      if (defined($ex12)) { foreach my $ex (@$ex12) { $prog->printf("%d %d\n", $ex->[0]->{ai}, $ex->[1]->{ai}); } }
      if (defined($ex13)) { foreach my $ex (@$ex13) { $prog->printf("%d %d\n", $ex->[0]->{ai}, $ex->[1]->{ai}); } }
      $prog->printf("endexclude\n");
      $prog->printf("alt\n");
      $prog->printf("endalt\n");
      $prog->printf("$param");
    } elsif ($term =~ /VDW/) {
      # Certain calculation parameters are expected
      if (!defined($parm->{cutnb})) { GENERAL::error("Expected 'cutnb' in parameter hash for vdW calculation"); }
      if (!defined($parm->{ctonnb})) { GENERAL::error("Expected 'ctonnb' in parameter hash for vdW calculation"); }
      if (!defined($parm->{ctofnb})) { GENERAL::error("Expected 'ctofnb' in parameter hash for vdW calculation"); }
      # Prepare parameter string
      my $param = "";
      my @aas; push(@aas, $atomsi); if (defined($atomsj)) { push(@aas, $atomsj); }
      foreach my $aa (@aas) {
        foreach my $ai (@$aa) {
          # fish out atom type info
          my $rni = $ai->{residue}->{resname};
          if (defined($ai->{residue}->{patch})) { # see if residue is patched
            $rni .= "_patch_" . $ai->{residue}->{patch};
            PDB::createPatchedResidue($T, $ai->{residue}->{resname}, $ai->{residue}->{patch}, $rni, 0);
          }
          my $ati = $T->{residues}->{$rni}->{atoms}->{$ai->{atomname}}->{type};
          if (!defined($ati)) { GENERAL::error("No atom type defined for atom $ai->{atomname} in residue " . PDB::resStr($ai->{residue})); }
          # print parameters
          if (defined($P->{vdW}{$ati}{emin})) { $param .= sprintf("%f", $P->{vdW}{$ati}{emin}); }
          else { GENERAL::error("No well depth defined for atom type $ati"); }
          if (defined($P->{vdW}{$ati}{rmin})) { $param .= sprintf(" %f", $P->{vdW}{$ati}{rmin}); }
          else { GENERAL::error("No van der Waals radius defined for atom type $ati"); }
          if (defined($P->{vdW}{$ati}{emin_alt})) {
            $param .= sprintf(" %f %f\n", $P->{vdW}{$ati}{emin_alt}, $P->{vdW}{$ati}{rmin_alt});
          } else {
            $param .= sprintf(" %f %f\n", $P->{vdW}{$ati}{emin}, $P->{vdW}{$ati}{rmin});
          }
        }
      }

      $prog->printf("$parm->{cutnb} $parm->{ctonnb} $parm->{ctofnb} $parm->{potm}\n");
      $prog->printf("exclude\n");
      if (defined($ex12)) { foreach my $ex (@$ex12) { $prog->printf("%d %d\n", $ex->[0]->{ai}, $ex->[1]->{ai}); } }
      if (defined($ex13)) { foreach my $ex (@$ex13) { $prog->printf("%d %d\n", $ex->[0]->{ai}, $ex->[1]->{ai}); } }
      $prog->printf("endexclude\n");
      $prog->printf("alt\n");
      if (defined($ex14)) { foreach my $ex (@$ex14) { $prog->printf("%d %d\n", $ex->[0]->{ai}, $ex->[1]->{ai}); } }
      $prog->printf("endalt\n");
      $prog->printf("$param");
    } elsif ($term =~ /DDE/) {
      # Certain calculation parameters are expected
      if (!defined($parm->{cutnb})) { GENERAL::error("Expected 'cutnb' in parameter hash for DDE calculation"); }
      if (!defined($parm->{ctonnb})) { GENERAL::error("Expected 'ctonnb' in parameter hash for DDE calculation"); }
      if (!defined($parm->{ctofnb})) { GENERAL::error("Expected 'ctofnb' in parameter hash for DDE calculation"); }
      if (!defined($parm->{rdie})) { GENERAL::error("Expected 'rdie' in parameter hash for DDE calculation"); }
      if (!defined($parm->{e14fac})) { GENERAL::error("Expected 'e14fac' in parameter hash for DDE calculation"); }
      # Prepare parameter string
      my $param = "";
      my @aas; push(@aas, $atomsi); if (defined($atomsj)) { push(@aas, $atomsj); }
      foreach my $aa (@aas) {
        foreach my $ai (@$aa) {
          my $rni = $ai->{residue}->{resname};
          if (defined($ai->{residue}->{patch})) { # see if residue is patched
            $rni .= "_patch_" . $ai->{residue}->{patch};
            PDB::createPatchedResidue($T, $ai->{residue}->{resname}, $ai->{residue}->{patch}, $rni, 0);
          }
          my $aci = $T->{residues}->{$rni}->{atoms}->{$ai->{atomname}}->{charge};
          if (!defined($aci)) { GENERAL::error("No charge defined for atom $ai->{atomname} in residue " . PDB::resStr($ai->{residue})); }
          $param .= sprintf("%f ", $aci);

          # only for use with 'linrep', print van der Waals radius information for each atom
          my $ati = $T->{residues}->{$rni}->{atoms}->{$ai->{atomname}}->{type};
          if (!defined($ati)) { GENERAL::error("No atom type defined for atom $ai->{atomname} in residue " . PDB::resStr($ai->{residue})); }
          if (defined($P->{vdW}{$ati}{rmin})) { $param .= sprintf(" %f", $P->{vdW}{$ati}{rmin}); }
          else { GENERAL::error("No van der Waals radius defined for atom type $ati"); }
          if (defined($P->{vdW}{$ati}{rmin_alt})) { $param .= sprintf(" %f\n", $P->{vdW}{$ati}{rmin_alt}); }
          else { $param .= sprintf(" %f\n", $P->{vdW}{$ati}{rmin}); }
        }
      }

      # Call ECPI
      if ($parm->{potm} eq "vdwonlylinrep") {
        $prog->printf("$parm->{rdie} $parm->{e14fac} $parm->{cutnb} $parm->{ctonnb} $parm->{ctofnb} switch\n");
      } else {
        $prog->printf("$parm->{rdie} $parm->{e14fac} $parm->{cutnb} $parm->{ctonnb} $parm->{ctofnb} $parm->{potm}\n");
      }
      $prog->printf("exclude\n");
      if (defined($ex12)) { foreach my $ex (@$ex12) { $prog->printf("%d %d\n", $ex->[0]->{ai}, $ex->[1]->{ai}); } }
      if (defined($ex13)) { foreach my $ex (@$ex13) { $prog->printf("%d %d\n", $ex->[0]->{ai}, $ex->[1]->{ai}); } }
      $prog->printf("endexclude\n");
      $prog->printf("alt\n");
      if (defined($ex14)) { foreach my $ex (@$ex14) { $prog->printf("%d %d\n", $ex->[0]->{ai}, $ex->[1]->{ai}); } }
      $prog->printf("endalt\n");
      $prog->printf("$param");
    } else {
      GENERAL::error("Do not know how to treat term '$terms'");
    }
  }

  if (exists($par{prog})) { return $prog; }
  $prog->close();
}


=head2

 Title   :  writeUnfoldEPIC
 Usage   :  
 Function:  Calculates and write all oligo-peptide unfolded state energy terms required by the given energy scheme that EPIC can handle.
 Returns :
 Args    : 1. ROTAMER object corresponding to the design problem
           2. Energy scheme
           3. Unfoldmer length
           4. EPIC calculation parameter hash
           5. Resume flag - optional.

=cut

sub writeUnfoldEPIC {
  my $rotamer = shift;
  my $escheme = shift;
  my $um = shift;
  my $parm = shift;
  my $resume = shift;
  my $terpatch = shift;
  if (!defined($resume)) { $resume = 0; }
  GENERAL::requireArgs($rotamer, $escheme, $um, $parm, $resume);

  # -------------------------------------------------------
  # Get the term tracker structure and setup the new entry
  my %en;
  if ($PRANK_DEF == 0) {
    # all terms EPIC can calculate
    $en{VDW} = ENERGY::addTrackerTerm($TERM_TRACKER, 'unfold', 'unfold_vdw_epic', 'Unfolded vdW using EPIC', "$um-peptide i-t + i-i vdW interaction energy");
    ENERGY::addSchemeToTerm($en{VDW}, '1p', (0), ('unfold_sf'));
    ENERGY::addSchemeToTerm($en{VDW}, 'gb-eef', (0), ('unfold_sf')); # test energy function
    ENERGY::addSchemeToTerm($en{VDW}, 'vdwp', (0), ('unfold_sf'));
    $en{EEF1} = ENERGY::addTrackerTerm($TERM_TRACKER, 'unfold', 'unfold_eef1_epic', 'Unfolded EEF1 using EPIC', "$um-peptide i-t EEF1 desolvation energy");
    ENERGY::addSchemeToTerm($en{EEF1}, '1p', (0), ('unfold_sf'));
    $en{DDE} = ENERGY::addTrackerTerm($TERM_TRACKER, 'unfold', 'unfold_dde_epic', 'Unfolded DDE using EPIC', "$um-peptide i-t + i-i DDE interaction energy");
    ENERGY::addSchemeToTerm($en{DDE}, '1p', (0), ('unfold_sf'));
  }
  # -------------------------------------------------------

  # Do we need to run this term under the current escheme?
  my $time = time(); my @terms;
  if (($escheme eq '1p') || ($escheme eq 'gb-eef') || ($escheme eq 'vdwp')) {
    if ($escheme eq '1p') { @terms = ('VDW', 'EEF1', 'DDE'); }
    if ($escheme eq 'gb-eef') { @terms = ('VDW'); }
    if ($escheme eq 'vdwp') { @terms = ('VDW'); }
    my $epic = undef;
    $rotamer->numify();
    # Divide the task among processes
    my $Nt = $rotamer->countResTotal();
    my $a = GENERAL::splitTasks($Nt, $NPROC_DEF);
    my $N = $a->[$PRANK_DEF]->{end} - $a->[$PRANK_DEF]->{beg} + 1;
    # Calculate each term
    foreach my $term (@terms) {
      print "Process $PRANK_DEF: Calculating unfold '$term' using EPIC...\n";
      # Have we done this before?
      my $tablef;
      if ($term eq "VDW") { $tablef = "$FS_DEF->{dunfold}/vdw_epic.tab"; }
      elsif ($term eq "DDE") { $tablef = "$FS_DEF->{dunfold}/dde_epic.tab"; }
      elsif ($term eq "EEF1") { $tablef = "$FS_DEF->{dunfold}/eef1_epic.tab"; }
      else { GENERAL::error("Do not know what to call final table for term '$term'"); }
      # If resuming a previous run, see if the total table has already been created
      if ($resume && GENERAL::fileCheck($tablef, "e") && (GENERAL::GetFileInfo($tablef, 'lines') == $rotamer->countSelfTotal())) {
        GENERAL::GLog("Skipping unfold EPIC '$term' calculation - has been completed already...");
        print "Skipping unfold EPIC '$term' calculation - has been completed already...\n";
        next;
      }
      # Read appropriate topology and parameter tables
      my ($P, $T);
      if ($term =~ /EEF/) {
        my $Peef = CHARMM::readEEFParameters(DEFINITIONS::getParmFile("solvpar.inp"));
        $T = CHARMM::readCHARMMTopology(DEFINITIONS::getParmFile("toph19_eef1.inp"));
        $T = CHARMM::readCHARMMTopology(DEFINITIONS::getParmFile("patchtop19_eef1.inp"), $T);
        $P = CHARMM::readCHARMMParameters(DEFINITIONS::getParmFile("param19_eef1.inp"));
        CHARMM::expandCHARMMParameters($P, $T);
        foreach my $at (keys(%$Peef)) {
          foreach my $val (keys(%{$Peef->{$at}})) {
            $P->{eef1}{$at}{$val} = $Peef->{$at}{$val};
          }
        }
      } elsif ($term =~ /(VDW|DDE)/) {
        $T = CHARMM::readCHARMMTopology(DEFINITIONS::getParmFile("toph19.inp"));
        $T = CHARMM::readCHARMMTopology(DEFINITIONS::getParmFile("patchtop19.inp"), $T);
        $P = CHARMM::readCHARMMParameters(DEFINITIONS::getParmFile("param19_90.inp"));
        CHARMM::expandCHARMMParameters($P, $T);
      }
      # Do the assigned part of the work
      my $num = 0;
      if ($N > 0) {
        # loop over site and residue
        for (my $si = 0; $si < scalar(@{$rotamer->{dslist}}); $si++) {
          my $ds = $rotamer->{dslist}->[$si];
          # Create template, number template atoms
          my $tfile = "$FS_DEF->{dunfold}/$ds->{chain}$ds->{iresnum}\.out";
          my $tpdb = PDB::new("$FS_DEF->{dunfold}/$ds->{chain}$ds->{iresnum}\.pdb");
          # patch the template if terminal patches are specified and this fragment includes the termini.
          # because the designed sites are ultimately derived from the template, this will properly patch everything
          if (defined($terpatch) && ($ds->{iresnum} <= 3)) {
            my $fres = $tpdb->{chain}->[0]->{res}->[0]; $fres->{patch} = $terpatch->{$fres->{chain}->{id}}{fir};
            PDB::createPatchedResidue($T, $fres->{resname}, $fres->{patch}, undef, 0);
          }
          if (defined($terpatch) && ($ds->{iresnum} >= $terpatch->{$tpdb->{chain}->[0]->{id}}{len}-2)) {
            my $lres = $tpdb->{chain}->[0]->{res}->[-1]; $lres->{patch} = $terpatch->{$lres->{chain}->{id}}{lst};
            PDB::createPatchedResidue($T, $lres->{resname}, $lres->{patch}, undef, 0);
          }
          my @ta = $tpdb->conAtoms();
          for (my $i = 0; $i < scalar(@ta); $i++) { $ta[$i]->{ai} = $i+1; }
          # Find central residue
          my $cr;
          if ($ds->{iresnum} < int($um/2)+1) { $cr = $ds->{iresnum}; }
          else { $cr = int($um/2) + 1; }
          my $res = $tpdb->{chain}->[0]->{res}->[$cr-1];
          for (my $ri = 0; $ri < scalar(@{$ds->{reslist}}); $ri++) {
            $num++;
            next if (($num < $a->[$PRANK_DEF]->{beg}) || ($num > $a->[$PRANK_DEF]->{end}));
            my $dresi = $ds->{reslist}->[$ri];
            my $rpdbi = $FS_DEF->{rot_pdbf}; $rpdbi =~ s/%%/$si/; $rpdbi =~ s/%/$ri/;
            my $inpi = $FS_DEF->{rot_outf}; $inpi =~ s/%%/$si/; $inpi =~ s/%/$ri/;
            # read one rotamer for this residue
            my $rotsi = PDB::new($rpdbi, "PDB", "CHARMM", 1);
            my $roti = $rotsi->{chain}->[0]->{res}->[0];
            # paste the sidechain into the template and number sidechain
            PDB::replaceSidechain($roti, $res);
            my @sa = PDB::sidechain($res);
            for (my $i = 0; $i < scalar(@sa); $i++) { $sa[$i]->{ai} = $i+1; }
            # get exclusion lists for rot-template and rot-rot interactions
            my ($ex12, $ex13, $ex14, $ex12s, $ex13s, $ex14s) = PDB::exclusionList(\@sa, "1-2 1-3 1-4", $T, 3);
#printf("12:\n"); PDB::printExclusionList($ex12);
#printf("12s:\n"); PDB::printExclusionList($ex12s);
#printf("13:\n"); PDB::printExclusionList($ex13);
#printf("13s:\n"); PDB::printExclusionList($ex13s);
#printf("14:\n"); PDB::printExclusionList($ex14);
#printf("14s:\n"); PDB::printExclusionList($ex14s);
            # PRO hack (to avoid interactions with the amide H, which is supposed to be absent from the template)
            if ($roti->{resname} eq "PRO") {
              if ($parm->{param} ne "19") {
                GENERAL::error("Don't know how to apply a PRO-specific exclusion list modification in parameter set $parm->{param}");
              }
              # 1. add (everything in side chain - H) into the exclusion list (hack)
              my $H = PDB::getAtomInRes($res, "H", -1);
              if ($H ne -1) {
                foreach my $sa (@sa) { my @arr = ($sa, $H); push(@$ex12, \@arr); }
                # 2. temporarily add H into the topology of PRO (put into a group by itself)
                $T->{residues}->{PRO}->{atoms}->{H}->{type} = "H";
                $T->{residues}->{PRO}->{atoms}->{H}->{group} = scalar(keys(%{$T->{residues}->{PRO}->{groups}})) + 1;
                $T->{residues}->{PRO}->{atoms}->{H}->{charge} = 0.0;
              }
            } elsif ($roti->{resname} =~ /PHE|TYR/) {
              # PHE and TYR hacks (to avoid in-ring 1-4 interactions - no interactions for these in CHARMM)
              if ($parm->{param} ne "19") {
                GENERAL::error("Don't know how to apply a PHE/TYR-specific exclusion list modification in parameter set $parm->{param}");
              }
              my @nex = ("CG", "CZ", "CD1", "CE2", "CE1", "CD2");
              for (my $i = 0; $i < scalar(@nex); $i += 2) {
                my @arr = (PDB::getAtomInRes($res, $nex[$i]), PDB::getAtomInRes($res, $nex[$i+1])); push(@$ex12s, \@arr);
              }
            } elsif ($roti->{resname} eq "TRP") {
              # TRP hack (to avoid in-ring/across-rings interactions, 1-4 or other - no interactions for these in CHARMM)
              my @nex = ("CG","CZ3", "CG","CZ2", "CD2","CH2", "CE2","CZ3", "CE3","CD1", "CE3","CZ2", "CE3","NE1", "CD1","CZ2", "NE1","CH2",
                         "CG","CH2", "CD1","CZ3", "CD1","CH2", "CZ3","NE1");
              for (my $i = 0; $i < scalar(@nex); $i += 2) {
                my @arr = (PDB::getAtomInRes($res, $nex[$i]), PDB::getAtomInRes($res, $nex[$i+1])); push(@$ex12s, \@arr);
              }
            }
            # calculate rot-template intereaction
            $epic = ENERGY::interfaceEPIC("$FS_DEF->{dunfold}/$si-$ri-rt.tab", $term, "pair", $T, $P, $parm, \@sa, $inpi, 'atomsj', \@ta, 'inpj', $tfile, 'sf', 1.0, 'ex12', $ex12, 'ex13', $ex13, 'ex14', $ex14, 'prog', $epic);
            # calculate rot-rot interaction
            if ($term ne "EEF1") {
              $epic = ENERGY::interfaceEPIC("$FS_DEF->{dunfold}/$si-$ri-rr.tab", $term, "self", $T, $P, $parm, \@sa, $inpi, 'sf', 1.0, 'ex12', $ex12s, 'ex13', $ex13s, 'ex14', $ex14s, 'prog', $epic);
            }
            if (time() - $time > 30*$LOCAL_TIMEOUT) { GENERAL::touchLocalSpace(); $time = time(); }
            # undo PRO hack
            if ($roti->{resname} eq "PRO") {
              delete($T->{residues}->{PRO}->{atoms}->{H});
              delete($T->{residues}->{PRO}->{atoms}->{H}->{group});
            }
          }
        }
      }
      if (defined($epic)) { $epic->close(); $epic = undef; } # close to cause a flush

      # Wait until all processes are done before putting the data together
      GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF, $NFS_DELAY);
      # Combine all tables into one sorted by rotamer
      if ($PRANK_DEF == 0) {
        my $ofh = GENERAL::GetOutFH($tablef);
        for (my $si = 0; $si < scalar(@{$rotamer->{dslist}}); $si++) {
          my $ds = $rotamer->{dslist}->[$si];
          for (my $ri = 0; $ri < scalar(@{$ds->{reslist}}); $ri++) {
            if ($term eq "EEF1") {
              my $rt = GENERAL::file2array("$FS_DEF->{dunfold}/$si-$ri-rt.tab");
              for (my $i = 0; $i < scalar(@$rt); $i++) {
                $ofh->printf("%e\n", $rt->[$i]);
              }
              GENERAL::crm("$FS_DEF->{dunfold}/$si-$ri-rt.tab");
            } else {
              my $rr = GENERAL::file2array("$FS_DEF->{dunfold}/$si-$ri-rr.tab");
              my $rt = GENERAL::file2array("$FS_DEF->{dunfold}/$si-$ri-rt.tab");
              if (scalar(@$rr) != scalar(@$rt)) { GENERAL::error(sprintf("%d r-r energies, but %d r-t energies!", scalar(@$rr), scalar(@$rt))); }
              for (my $i = 0; $i < scalar(@$rr); $i++) {
                $ofh->printf("%e\n", $rr->[$i]+$rt->[$i]);
              }
              GENERAL::crm("$FS_DEF->{dunfold}/$si-$ri-rr.tab");
              GENERAL::crm("$FS_DEF->{dunfold}/$si-$ri-rt.tab");
            }
          }
        }
        close($ofh);
        # Document the created energy table
        push(@{$en{$term}->{files}}, File::Spec->abs2rel($tablef, $FS_DEF->{dhome}));
      }
      # Wait until all old files are erased before moving on
      GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF);
    }
  }

  # Save term tracker
  if ($PRANK_DEF == 0) {
    GENERAL::SaveStruct($TERM_TRACKER, $FS_DEF->{term_trackerf});
  }
}


=head2

 Title   :  writeUnfoldSA
 Usage   :  ENERGY::writeUnfoldSA($rotamer, $escheme, $unfoldmer);
 Function:  Calculates the solvent accessible area in the oligopeptide unfolded
            state. Uses naccess for this purpose.
 Returns :  Nothing
 Args    :  (1) rotamer object
            (2) energy scheme id
            (3) length of the unfoldmer
            (4) optional: an array of atom type regular expressions (will be applied to
                strings returned by PDB::atomStr())
            (5) optional: resume flag (default 0). If set, will not overwrite a previously
                computed complete table.

=cut

sub writeUnfoldSA {

  my $rotamer = shift;
  my $escheme = shift;
  my $unfoldmer = shift;
  my $resume = shift;
  if (!defined($resume)) { $resume = 0; }
  GENERAL::requireArgs($rotamer, $escheme, $unfoldmer);

  # If resuming a previous run, see if the total table has already been created
  if ($resume && GENERAL::fileCheck($FS_DEF->{unfold_SA_tabf}, "e") && (GENERAL::GetFileInfo($FS_DEF->{unfold_SA_tabf}, 'lines') == $rotamer->countSelfTotal())) {
    GENERAL::GLog("Skipping unfold SA calculation - has been completed already...");
    return;
  }
  
  my @types = ('_C[^_]*$', '_N[^_]*$', '_O[^_]*$', '_S[^_]*$', '_H[^_]*$');
  my $pr = 1.4;

  # -------------------------------------------------------
  # Get the term tracker structure and setup the new entry
  my $en;
  if ($PRANK_DEF == 0) {
    if (scalar(@_)) {
      @types = @_;
      # If atom types specified, we don't know how to treat them in other energy schemes.
      # For the current one, we assume uniform scaling by sa_sf.
      my @names;
      foreach my $t (@types) { push(@names, "SASA of $t"); }
      $en = ENERGY::addTrackerTerm($TERM_TRACKER, 'unfold', 'unfold_sa', 'Oligopeptide unfolded state solvent accessible area (scaled).', @names);
      ENERGY::addSchemeToTerm($en, $escheme, GENERAL::linspace(0, scalar(@types)-1, 1), GENERAL::ones(scalar(@types), "sa_sf"));
    } else {
      $en = ENERGY::addTrackerTerm($TERM_TRACKER, 'unfold', 'unfold_sa', 'Oligopeptide unfolded state solvent accessible area (scaled).',
            ("SASA of C in unfolded",
              "SASA of N in unfolded",
              "SASA of O in unfolded",
              "SASA of S in unfolded",
              "SASA of H in unfolded"));
      # Do the following for every energy scheme this term belongs to
      ENERGY::addSchemeToTerm($en, '2', (0, 1, 2, 3, 4), ('C_sa_sf', 'N_sa_sf', 'O_sa_sf', 'S_sa_sf', 'H_sa_sf'));
      ENERGY::addSchemeToTerm($en, '2b', (0, 1, 2, 3, 4), ('C_sa_sf', 'N_sa_sf', 'O_sa_sf', 'S_sa_sf', 'H_sa_sf'));
      ENERGY::addSchemeToTerm($en, '2c', (0, 1, 2, 3, 4), ('C_sa_sf', 'N_sa_sf', 'O_sa_sf', 'S_sa_sf', 'H_sa_sf'));
      ENERGY::addSchemeToTerm($en, 'c1', (0, 1, 2, 3, 4), ('C_sa_sf', 'N_sa_sf', 'O_sa_sf', 'S_sa_sf', 'H_sa_sf'));
      ENERGY::addSchemeToTerm($en, '3', (0, 1, 2, 3, 4), ('C_sa_sf', 'N_sa_sf', 'O_sa_sf', 'S_sa_sf', 'H_sa_sf'));
      ENERGY::addSchemeToTerm($en, '4', (0, 1, 2, 3, 4), ('C_sa_sf', 'N_sa_sf', 'O_sa_sf', 'S_sa_sf', 'H_sa_sf'));
    }
  }
  # -------------------------------------------------------

  # Do we need to run this term under the current escheme?
  if (($escheme eq '2') || ($escheme eq '2b') || ($escheme eq '2c') || ($escheme eq 'c1')  || ($escheme eq '3')  || ($escheme eq '4')) {
    my $cdir = GENERAL::GetDir();
    my $wdir = "$LOCAL_DEF/sa_unfold_$PRANK_DEF";
    GENERAL::cmkdir($wdir);
    GENERAL::cchdir($wdir);

    # Divide the task among processes
    my $a = GENERAL::splitTasks($rotamer->countSelfTotal(), $NPROC_DEF);
    my $N = $a->[$PRANK_DEF]->{end} - $a->[$PRANK_DEF]->{beg} + 1;
    if ($N > 0) {
      $rotamer->markRotamers($a->[$PRANK_DEF]->{beg}, $a->[$PRANK_DEF]->{end}, "visit");


      # Read in atomic vdW radii so that ENERGY::getAtomicSA does not need
      # to do this every time it is called.
      my $R = PDB::sizeLookup();
      print "Running unfolded SA calculation (process $PRANK_DEF)...\n";

      # Open output file
      my $energytab = "$FS_DEF->{dunfold}/tmp_$PRANK_DEF";
      my $efh = GENERAL::GetOutFH($energytab);

      # Visit each site
      my $i = -1;
      foreach my $site (@{$rotamer->{dslist}}) {
        $i++;
        next if (!$site->{visit});
        my $chain = $site->{chain};
        my $iresnum = $site->{iresnum};
        # load the corresponding unfoldmer
        my $ufname = "$FS_DEF->{dunfold}/$chain$iresnum\.pdb";
        my $updb = PDB::new($ufname);
        # get the residue corresponding to this site
        my $cr;
        if ($iresnum < int($unfoldmer/2)+1) { $cr = $iresnum; }
        else { $cr = int($unfoldmer/2) + 1; }
        my $res = $updb->{chain}->[0]->{res}->[$cr-1];

        my $j = -1;
        foreach my $r ( @{$site->{reslist}}) {
          $j++;
          next if (!$r->{visit});
          print "Unfold SA for site $i residue $j (process $PRANK_DEF)...\n";
          # Read in all the rotamers of the current aa at the current site
          my $rotpdb = $FS_DEF->{rot_pdbf};
          $rotpdb =~ s/%%/$i/;
          $rotpdb =~ s/%/$j/;
          my $rots = PDB::new($rotpdb, "PDB", "CHARMM");
          my @rots = $rots->ConRes();

          # Place each rotamer onto the unfoldmer
          foreach (my $ri=0; $ri < scalar(@{$r->{rotamerlist}}); $ri++) {
            next if (!$r->{rotamerlist}->[$ri]->{visit});
            my $rot = $rots[$ri];
            PDB::copyResidue($rot, $res);
            # calculate the solvent accessible area (by atom) of the rotamer
            my @sa = ENERGY::getAtomicSA($updb, $pr, $res->{atom}, $R);
            if (scalar(@sa) != scalar(@{$res->{atom}})) {
              die "Error in writeUnfoldSA: " . scalar(@{$res->{atom}}) . " atoms in\n".
                  "residue " . $updb->resStr($res) . ", " .
                  "but only " . scalar(@sa) . " atom areas returned by naccess!\n";
            }
            # distribute the area by atom type
            my %ta;
            foreach my $type (@types) { $ta{$type} = 0; }
            for (my $ai = 0; $ai < scalar(@{$res->{atom}}); $ai++) {
              my $atom = $res->{atom}->[$ai];
              foreach my $type (@types) {
                if (PDB::atomStr($atom) =~ /$type/) {
                  $ta{$type} += $sa[$ai];
                  last;
                }
              }
            }
            # Save the areas
            foreach my $type (@types) {
              $efh->printf("%.3f ", $ta{$type});
            }
            $efh->printf("\n");
          }
        }
      }
      close($efh);
      $rotamer->unmarkRotamers("visit");
    }

    # Wait untill all processes get to this point
    GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF, $NFS_DELAY);

    # Process 0 is to collect all data together
    if ($PRANK_DEF == 0) {
#      my $i = -1; # site number
      my $k = 1; # rotamer number
      my $pr = 0; # process rank of the currently examined process
      my $pfh; # file handle of file made by the currently examined process
      my $energytab = $FS_DEF->{unfold_SA_tabf};
      my $efh = GENERAL::GetOutFH($energytab); # output file handle
      foreach my $site (@{$rotamer->{dslist}}) {
#        $i++;
#        # Open output file
#        my $energytab = $FS_DEF->{unfold_SA_tabf};
#        $energytab =~ s/%%/$site->{chain}/; $energytab =~ s/%/$site->{iresnum}/;
#        my $efh = GENERAL::GetOutFH($energytab);

        foreach my $res (@{$site->{reslist}}) {
          foreach my $rot (@{$res->{rotamerlist}}) {
            if ($k == $a->[$pr]->{beg}) {
              $pfh = GENERAL::GetInFH("$FS_DEF->{dunfold}/tmp_$pr");
            }
            my $line = <$pfh>;
            if (!defined($line)) {
              GENERAL::error("Process $pr has fewer lines than expected!");
            }
            print $efh $line;
            if ($k == $a->[$pr]->{end}) {
              my $line = <$pfh>;
              if (defined($line)) {
                GENERAL::error("Process $pr has more lines than expected!");
              }
              close($pfh);
              GENERAL::crm("$FS_DEF->{dunfold}/tmp_$pr");
              $pr++;
            }
            $k++;
          }
        }
      }
      close($efh);
      # Document the created energy table
      push(@{$en->{files}}, File::Spec->abs2rel($energytab, $FS_DEF->{dhome}));
    }
    GENERAL::cchdir($cdir);
    GENERAL::crmdir($wdir);
  }

  # Save term tracker
  if ($PRANK_DEF == 0) {
    GENERAL::SaveStruct($TERM_TRACKER, $FS_DEF->{term_trackerf});
  }
#  GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF, $NFS_DELAY);
}




=head2

 Title   :  writeUnfoldGB
 Usage   :  ENERGY::writeUnfoldGB($rotamer, $escheme, $unfoldmer);
 Function:  Calculates the electrostatic energy in the oligopeptide unfolded
            state as the sum of the polarization energy and the electrostatic
            energy in the reference state ("vacuum").
 Returns :  Nothing
 Args    :  (1) rotamer object
            (2) energy scheme id
            (3) length of the unfoldmer
            (4) internal dielectric constant
            (5) optional: resume flag (default 0). If set, will not overwrite a previously
                computed complete table.
=cut

sub writeUnfoldGB {

  my $rotamer = shift;
  my $escheme = shift;
  my $unfoldmer = shift;
  my $eps_in = shift;
  my $resume = shift;
  if (!defined($resume)) { $resume = 0; }
  GENERAL::requireArgs($rotamer, $escheme, $unfoldmer, $eps_in);

  # If resuming a previous run, see if the total table has already been created
  if ($resume && GENERAL::fileCheck($FS_DEF->{unfold_GB_tabf}, "e") && (GENERAL::GetFileInfo($FS_DEF->{unfold_GB_tabf}, 'lines') == $rotamer->countSelfTotal())) {
    GENERAL::GLog("Skipping unfold GB calculation - has been completed already...");
    return;
  }
  
  # -------------------------------------------------------
  # Get the term tracker structure and setup the new entry
  my $en;
  if ($PRANK_DEF == 0) {
    $en = ENERGY::addTrackerTerm($TERM_TRACKER, 'unfold', 'unfold_gb', "Oligopeptide unfolded state electrostatic energy using \'smart template GB\'.",
             ("Unfold atomic self polarization energies (qi^2/a)",
              "Unfold screening of sc-template interactions",
              "Unfold screening of sc-sc interactions",
              "Unfold Coulombic sc-template interactions",
              "Unfold Coulombic sc-sc interactions"));
    # Do the following for every energy scheme this term belongs to
    ENERGY::addSchemeToTerm($en, '2', (0, 1, 2, 3), ('1.0', '1.0', '1.0', '1.0'));
    ENERGY::addSchemeToTerm($en, '2b', (0, 1, 3), ('1.0', '1.0', '1.0'));
    ENERGY::addSchemeToTerm($en, '2c', (0, 1, 2, 3, 4), ('1.0', '1.0', '1.0', '1.0', '1.0'));
    ENERGY::addSchemeToTerm($en, 'c1', (0, 1, 2, 3, 4), ('1.0', '1.0', '1.0', '1.0', '1.0'));
  }
  # -------------------------------------------------------

  # Do we need to run this term under the current escheme?
  if (($escheme eq '2') || ($escheme eq '2b') || ($escheme eq '2c') || ($escheme eq 'c1')) {
    # Go to the appropriate directory
    my $cdir = GENERAL::GetDir();
    my $wdir = "$LOCAL_DEF/gb_unfold_$PRANK_DEF";
    GENERAL::cmkdir($wdir);
    GENERAL::cchdir($wdir);

    # Divide the task among processes
    my $a = GENERAL::splitTasks($rotamer->countSelfTotal(), $NPROC_DEF);
    my $N = $a->[$PRANK_DEF]->{end} - $a->[$PRANK_DEF]->{beg} + 1;
    if ($N > 0) {
      $rotamer->markRotamers($a->[$PRANK_DEF]->{beg}, $a->[$PRANK_DEF]->{end}, "visit");

      # Open output file for the current process
      my $efh = GENERAL::GetOutFH("$FS_DEF->{dunfold}/tmp_$PRANK_DEF");

      # Read in atomic vdW radii and atomic charges.
      my $R = PDB::sizeLookup();
      my $C = PDB::chargeLookup();

      print "Running unfolded GB calculation (process $PRANK_DEF)...\n";
      my $i = -1; # site number
      foreach my $site (@{$rotamer->{dslist}}) {
        $i++;
        next if (!$site->{visit});
        my $chain = $site->{chain};
        my $iresnum = $site->{iresnum};
        # load the corresponding unfoldmer
        my $ufname = "$FS_DEF->{dunfold}/$chain$iresnum\.pdb";
        my $updb = PDB::new($ufname);
        $updb->charge($C);
        # get the residue corresponding to this site
        my $cr;
        if ($iresnum < int($unfoldmer/2)+1) {
          $cr = $iresnum;
        } else {
          $cr = int($unfoldmer/2) + 1;
        }
        my $res = $updb->{chain}->[0]->{res}->[$cr-1];

        my $j = -1; # residue number
        foreach my $r ( @{$site->{reslist}}) {
          $j++;
          next if (!$r->{visit});
          print "Unfold GB for site $i residue $j (process $PRANK_DEF)...\n";
          # Read in all the rotamers of the current aa at the current site
          my $rotpdb = $FS_DEF->{rot_pdbf};
          $rotpdb =~ s/%%/$i/;
          $rotpdb =~ s/%/$j/;
          my $rots = PDB::new($rotpdb, "PDB", "CHARMM");
          my @rots = $rots->ConRes();
          if (scalar(@rots) != scalar(@{$r->{rotamerlist}})) {
            GENERAL::error("Inconsistent number of rotamers in the ROTAMER object and $rotpdb!");
          }

          # Place each rotamer onto the unfoldmer
          foreach (my $ri=0; $ri < scalar(@{$r->{rotamerlist}}); $ri++) {
            next if (!$r->{rotamerlist}->[$ri]->{visit});
            my $rot = $rots[$ri];
            PDB::copyResidue($rot, $res);
            # charge the structure
            PDB::charge($res->{atom}, $C);
            # Calculate the polarization energy
            # first, get the Born radii of all atoms in the unfoldmer
            my @ta = $updb->conAtoms("b_$updb->{chain}->[0]->{id}_$cr");
            my @a = @ta;
            my @sa = PDB::sidechain($res);
            push(@a, @sa);
            my %br = ENERGY::pepAtomicSolvation($updb, $eps_in, $R, \@a);
            my $Gs = ENERGY::getGBSelf(\%br, $eps_in, \@sa);
            my $Gst = ENERGY::getGBInterAB(\%br, $eps_in, \@ta, \@sa);
            my $Gss = ENERGY::getGBInterAA(\%br, $eps_in, \@sa);
            # calculate the reference state (vacuum) electrostatic energy
            my $Gvs = ENERGY::getCoulombInter($eps_in, \@ta, \@sa);
            my $Gvi = ENERGY::getCoulombTotal($eps_in, \@sa);
            # Save the values
            $efh->printf("%.3f %.3f %.3f %.3f %.3f\n", $Gs, $Gst, $Gss, $Gvs, $Gvi);
          }
        }
      }
      close($efh);
      $rotamer->unmarkRotamers("visit");
    }

    # Wait untill all processes get to this point
    GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF, $NFS_DELAY);

    # Process 0 is to collect all data together
    if ($PRANK_DEF == 0) {
#      my $i = -1; # site number
      my $k = 1; # rotamer number
      my $pr = 0; # process rank of the currently examined process
      my $pfh; # file handle of file made by the currently examined process
      my $energytab = $FS_DEF->{unfold_GB_tabf};
      my $efh = GENERAL::GetOutFH($energytab); # output file handle
      foreach my $site (@{$rotamer->{dslist}}) {
#        $i++;
#        # Open output file
#        my $energytab = $FS_DEF->{unfold_GB_tabf};
#        $energytab =~ s/%%/$site->{chain}/; $energytab =~ s/%/$site->{iresnum}/;
#        my $efh = GENERAL::GetOutFH($energytab);

        foreach my $res (@{$site->{reslist}}) {
          foreach my $rot (@{$res->{rotamerlist}}) {
            if ($k == $a->[$pr]->{beg}) {
              $pfh = GENERAL::GetInFH("$FS_DEF->{dunfold}/tmp_$pr");
            }
            my $line = <$pfh>;
            if (!defined($line)) {
              GENERAL::error("Process $pr has fewer lines than expected!");
            }
            print $efh $line;
            if ($k == $a->[$pr]->{end}) {
              my $line = <$pfh>;
              if (defined($line)) {
                GENERAL::error("Process $pr has more lines than expected!");
              }
              close($pfh);
              GENERAL::crm("$FS_DEF->{dunfold}/tmp_$pr");
              $pr++;
            }
            $k++;
          }
        }
      }
      close($efh);
      # Document the created energy table
      push(@{$en->{files}}, File::Spec->abs2rel($energytab, $FS_DEF->{dhome}));
    }
    GENERAL::cchdir($cdir);
    GENERAL::crmdir($wdir);
  }

  # Save term tracker
  if ($PRANK_DEF == 0) {
    GENERAL::SaveStruct($TERM_TRACKER, $FS_DEF->{term_trackerf});
  }
#  GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF, $NFS_DELAY);
}



=head2

 Title   :  writeUnfoldEEF
 Usage   :  ENERGY::writeUnfoldEEF($rotamer,'vdw', 'elec','hbond','dihe','asp');
 Function:  Generates and runs CHARMM scripts for unfolded reference state calculations.
            One script is generated per design site.
 Returns :  A string concatination of script names (separated by a space)
 Args    :  (1) rotamer object
            (2) energy scheme id
            NOTE:  in charmm, eisenberg's asp value and
            eef1 value are stored in same variable asp!
            (3) unfoldemer length

=cut

sub writeUnfoldEEF {

  my $rotamer = shift;
  my $escheme = shift;
  my $unfoldmer = shift;

  GENERAL::requireArgs($rotamer, $escheme, $unfoldmer);

  # -------------------------------------------------------
  # Get the term tracker structure and setup the new entry
  my $en;
  if ($PRANK_DEF == 0) {
#    $th = GENERAL::LoadStruct($FS_DEF->{term_trackerf}, 1);
    $en = ENERGY::addTrackerTerm($TERM_TRACKER, 'unfold', 'unfold_eef', 'Solvation energy of unfolded state calculated using CHARMM EEF',
           ("Desovlation energy from EEF",
           "Van Der Waals energy from EEF",
           "Electrostatic energy from EEF with reduced charge"));
    # Do the following for every energy scheme this term belongs to
    ENERGY::addSchemeToTerm($en, '1', (0), ('unfold_sf'));
    ENERGY::addSchemeToTerm($en, '1h', (0), ('unfold_sf'));
    ENERGY::addSchemeToTerm($en, 'c1', (0), ('unfold_sf'));
    ENERGY::addSchemeToTerm($en, '1e', (0, 2), ('unfold_sf', 'unfold_sf'));
    ENERGY::addSchemeToTerm($en, '1m', (0, 2), ('unfold_sf', 'unfold_sf'));
  }
  # -------------------------------------------------------

  # Do we need to run this term under the current escheme?
  if (($escheme eq '1') || ($escheme eq '1h') || ($escheme eq 'c1') || ($escheme eq '1e') || ($escheme eq '1m')) {
    # Divide the task among processes
    my $a = GENERAL::splitTasks(scalar(@{$rotamer->{dslist}}), $NPROC_DEF);
    my $N = $a->[$PRANK_DEF]->{end} - $a->[$PRANK_DEF]->{beg} + 1;
    if ($N > 0) {
      # Go to the appropriate directory
      my $cdir = GENERAL::GetDir();
      GENERAL::cchdir($FS_DEF->{dunfold});
      print "Running unfolded EEF calculation (process $PRANK_DEF)...\n";
      my $chidef=$rotamer->{chidef};
      my @inp = ();

      # Create scripts
      for (my $si = 0; $si < scalar(@{$rotamer->{dslist}}); $si++) {
        next if (($si+1 < $a->[$PRANK_DEF]->{beg}) || ($si+1 > $a->[$PRANK_DEF]->{end}));
        my $site = $rotamer->{dslist}->[$si];

        # for each designed site, get chain id and residue number
        my $j=0;
        my $chain = $site->{chain};
        my $lchain = lc($chain);
        my $iresnum = $site->{iresnum};
        my $charmminp = "$LOCAL_DEF/$chain$iresnum\_eef.inp";
        push(@inp, $charmminp);

        # now write charmm input script
        # we are using charmm instead of multe.
        # since we want to turn on eef1 potential
        my $charmm = CHARMM::new("$CHARMM_DEF", $charmminp);

        # using eef1 top and parm file
        $charmm->loadParm('EEF1');
        my $fname = $chain.$iresnum.'.crd';

        # given a crd file, setup sequence read in coordinates
        # build missing atoms and hydrogens.
        $charmm->setupFromCRD($fname);
        # find the index of the design rotamer
        my $cr;
        if ($iresnum < int($unfoldmer/2)+1) {
          $cr = $iresnum;
        } else {
          $cr = int($unfoldmer/2) + 1;
        }

        # Document the created energy table
        my $energytab = "site\_$si\_eef.tab";
        MPI_Send(File::Spec->rel2abs($energytab), 0, 0, $MPI_COMM_WORLD_DEF);
        $charmm->openCard($energytab, 18);
        $charmm->defineUnfoldTemplate($chain, $cr);

        foreach my $r ( @{$site->{reslist}}) {
          my $resname = $r->{resname};
          $charmm->deleteDesignSite();
          $charmm->verbose("bomlev -2");
          $charmm->patchRes($resname, $chain, $cr);
          $charmm->defineUnfoldTemplate($chain, $cr);
          $charmm->defineGroup($chain, $cr);

          my $k=0;
          foreach my $rot ( @{$r->{rotamerlist}}) {
            $charmm->placeRotamer($chain, $resname, $cr, $chidef, $rot);
            if ($k == 0) {
              if ($escheme eq "1m") { $charmm->setupIMM1(); }
              else { $charmm->setupEEF1(); }
            }
            $charmm->inte('sele group end', 'sele template end', 18, "asp vdw elec");
            $k++;
          }
          # If we just patched a proline, the amide hydrogen was removed, so re-initialize the structure
          if ($resname =~ /PRO/) {
            $charmm->verbose("delete atom select all end");
            $charmm->setupFromCRD($fname);
            $charmm->defineUnfoldTemplate($chain, $cr);
          }

          $j++;
          # increase the index $j for $Reslist
        }
        $charmm->closeCard(18);
        $charmm->endofstory();
      }
      my $unf_eef_inp = join(" ", @inp);

      # Run the scripts
      CHARMM::runCharmm($CHARMM_DEF, $unf_eef_inp);
      # Do we need to delete the scripts?
      if ($FS_DEF->{del_charmm_inp}) {
        GENERAL::crm(@inp);
      }
      GENERAL::cchdir($cdir);
    }
    
    # Update and save term tracker before the actual energy file is written. This is done in
    # the interest of time (so that only one delay has to be inserted)
    push(@{$en->{files}}, File::Spec->abs2rel($FS_DEF->{unfold_eef_tabf}, $FS_DEF->{dhome}));
    if ($PRANK_DEF == 0) {
      GENERAL::SaveStruct($TERM_TRACKER, $FS_DEF->{term_trackerf});
    }
    # Wait to make sure that the energy tracker file is consistent as well as that individual
    # energy tables are consistent before we concatinate them into one
    GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF, $NFS_DELAY);
    
    # Concatinate into one energy file and document it
    if ($PRANK_DEF == 0) {
      my $ofh = GENERAL::GetOutFH($FS_DEF->{unfold_eef_tabf});
      for (my $pr = 0; $pr < $NPROC_DEF; $pr++) {
        for (my $pi = 1; $pi <= ($a->[$pr]->{end} - $a->[$pr]->{beg} + 1); $pi++) {
          my $energytab = MPI_Recv($pr, 0, $MPI_COMM_WORLD_DEF);
          my $fh = GENERAL::GetInFH($energytab);
          while (<$fh>) {
            s/[ \t]+/ /g; s/^[ \t]//; s/[ \t]$//;
            $ofh->print($_);
          }
          close($fh);
          GENERAL::crm($energytab);
        }
      }
      close($ofh);
    }
  }
}



=head2

 Title   :  writeUnfoldEEF
 Usage   :  ENERGY::writeUnfoldEEF($rotamer, $escheme, $crd_file, $unfoldmer);
 Function:  Generates and runs CHARMM scripts for unfolded EEF calculations.
 Returns :  Nothing
 Args    :  (1) rotamer object
            (2) energy scheme id
            (3) CRD file of the original structure
            (4) unfoldemer length
            (5) optional: resume flag (default 0). If set, will not overwrite a previously
                computed complete table.

=cut

sub writeUnfoldEEF_fast {

  my $rotamer = shift;
  my $escheme = shift;
  my $pdb = shift;
  my $unfoldmer = shift;
  my $resume = shift;
  if (!defined($resume)) { $resume = 0; }
  GENERAL::requireArgs($rotamer, $escheme, $unfoldmer);

  # If resuming a previous run, see if the total table has already been created
  if ($resume && GENERAL::fileCheck($FS_DEF->{unfold_eef_tabf}, "e") && (GENERAL::GetFileInfo($FS_DEF->{unfold_eef_tabf}, 'lines') == $rotamer->countSelfTotal())) {
    GENERAL::GLog("Skipping unfold EEF calculation - has been completed already...");
    return;
  }
  
  # -------------------------------------------------------
  # Get the term tracker structure and setup the new entry
  my $en;
  if ($PRANK_DEF == 0) {
    $en = ENERGY::addTrackerTerm($TERM_TRACKER, 'unfold', 'unfold_eef', 'Solvation energy of unfolded state calculated using CHARMM EEF',
           ("Desovlation energy from EEF",
           "Van Der Waals energy from EEF",
           "Electrostatic energy from EEF with reduced charge"));
    # Do the following for every energy scheme this term belongs to
    ENERGY::addSchemeToTerm($en, '1', (0), ('unfold_sf'));
    ENERGY::addSchemeToTerm($en, '1h', (0), ('unfold_sf'));
    ENERGY::addSchemeToTerm($en, '1s', (0), ('unfold_sf'));
    ENERGY::addSchemeToTerm($en, 'c1', (0), ('unfold_sf'));
    ENERGY::addSchemeToTerm($en, '1e', (0, 2), ('unfold_sf', 'unfold_sf'));
    ENERGY::addSchemeToTerm($en, '1m', (0, 2), ('unfold_sf', 'unfold_sf'));
#    ENERGY::addSchemeToTerm($en, '1p', (0), ('unfold_sf'));
  }
  # -------------------------------------------------------

  # Do we need to run this term under the current escheme?
  if (($escheme eq '1') || ($escheme eq '1h') || ($escheme eq 'c1') || ($escheme eq '1s') || ($escheme eq '1e') || ($escheme eq '1m')) {
    # Go to the appropriate directory
    my $cdir = GENERAL::GetDir();
    GENERAL::cchdir($FS_DEF->{dunfold});
    $rotamer->numify();
    # Divide the task among processes
    my $a = GENERAL::splitTasks(scalar(@{$rotamer->{dslist}}), $NPROC_DEF);
    my $N = $a->[$PRANK_DEF]->{end} - $a->[$PRANK_DEF]->{beg} + 1;
    # Make a poly GLY - EEF is sensitive to this (because it is defined through group interactions)
    my $crdf = "$LOCAL_DEF/polyGLY.crd";
    my $cpdb = $pdb->clone();
    foreach my $chain (@{$cpdb->{chain}}) {
      foreach my $res (@{$chain->{res}}) {
        PDB::sMutate($res, "GLY");
      }
    }
    $cpdb->writeCRD($crdf, "");

    if ($N > 0) {
      my $charmminp = "$LOCAL_DEF/unfold_eef$PRANK_DEF.inp";
      my $charmm = CHARMM::new("$CHARMM_DEF", $charmminp);
      $charmm->loadParm('EEF1');
      $charmm->setupFromCRD($crdf);
      if ($escheme eq "1m") { $charmm->setupIMM1(); }
      else { $charmm->setupEEF1(); }


      print "Running unfolded EEF calculation (process $PRANK_DEF)...\n";
      my $chidef = $rotamer->{chidef};
      my @inp = ();

      # Create the script
      for (my $si = 0; $si < scalar(@{$rotamer->{dslist}}); $si++) {
        next if (($si+1 < $a->[$PRANK_DEF]->{beg}) || ($si+1 > $a->[$PRANK_DEF]->{end}));
        my $site = $rotamer->{dslist}->[$si];
        my $energytab = "site\_$si\_eef.tab";

        # if resuming a previous run, skip if this site has been completely treated
        if ($resume && GENERAL::fileCheck($energytab, "e") && (GENERAL::GetFileInfo($energytab, 'lines') == $site->{numvalidrot})) {
          print "Skipping site $si - has already been calculated in a previous run!\n";
          next;
        }

        my $chain = $site->{chain};
        my $iresnum = $site->{iresnum};
        my @sel;
        for (my $i = $iresnum - ($unfoldmer-1)/2; $i <= $iresnum + ($unfoldmer-1)/2; $i++) {
          next if ($i <= 0);
          push(@sel, ".byres. atom $chain $i N");
        }

        # Document the created energy table
#        MPI_Send(File::Spec->rel2abs($energytab), 0, 0, $MPI_COMM_WORLD_DEF);
        $charmm->openCard($energytab, 18);

        foreach my $r ( @{$site->{reslist}}) {
          my $resname = $r->{resname};
          $charmm->patchRes($resname, $chain, $iresnum);
          # Need to define things after patching - things change upon patching
          $charmm->verbose("define template select (" . join(" .or. ", @sel) . ")-\n.and. (" . CHARMM::backbone() . ") end\n");

          my $k=0;
          foreach my $rot ( @{$r->{rotamerlist}}) {
            $charmm->placeRotamer($chain, $resname, $iresnum, $chidef, $rot);
            if ($k == 0) {
              if ($escheme eq "1m") { $charmm->setupIMM1(); }
              else { $charmm->setupEEF1(); }
            }
            $charmm->inte("sele .byres. (atom $chain $iresnum N) .and. .not. template end", "sele template end", 18, "asp vdw elec");
            $k++;
          }
          # If we just patched a proline, the amide hydrogen was removed, so re-initialize the structure
          if ($resname =~ /PRO/) {
            $charmm->verbose("delete atom select all end");
            $charmm->setupFromCRD($crdf);
          }
        }
        $charmm->closeCard(18);
        # Put GLY back
        $charmm->patchRes("GLY", $chain, $iresnum);
      }
      $charmm->endofstory();
      # Do we need to delete the scripts?
      if ($FS_DEF->{del_charmm_inp}) {
        GENERAL::crm($charmminp);
      }
    }
    GENERAL::crm($crdf);
    
    # Wait to make sure that the energy tracker file is consistent as well as that individual
    # energy tables are consistent before we concatinate them into one
    GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF, $NFS_DELAY);
    
    # Concatinate into one energy file and document it
    if ($PRANK_DEF == 0) {
      my $ofh = GENERAL::GetOutFH($FS_DEF->{unfold_eef_tabf});
      for (my $si = 0; $si < scalar(@{$rotamer->{dslist}}); $si++) {
        my $energytab = "site\_$si\_eef.tab";
#      for (my $pr = 0; $pr < $NPROC_DEF; $pr++) {
#        for (my $pi = 1; $pi <= ($a->[$pr]->{end} - $a->[$pr]->{beg} + 1); $pi++) {
#          my $energytab = MPI_Recv($pr, 0, $MPI_COMM_WORLD_DEF);
        my $fh = GENERAL::GetInFH($energytab);
        while (<$fh>) {
          s/[ \t]+/ /g; s/^[ \t]//; s/[ \t]$//;
          $ofh->print($_);
        }
        close($fh);
        GENERAL::crm($energytab);
#        }
#      }
      }
      close($ofh);
      push(@{$en->{files}}, File::Spec->abs2rel($FS_DEF->{unfold_eef_tabf}, $FS_DEF->{dhome}));
    }
    GENERAL::cchdir($cdir);
  }
  # Update and save term tracker
  if ($PRANK_DEF == 0) {
    GENERAL::SaveStruct($TERM_TRACKER, $FS_DEF->{term_trackerf});
  }
}


=head2

 Title   :  writeUnfoldDDE_old
 Usage   :  ENERGY::writeUnfoldDDE_old($rotamer, $escheme);
 Function:  Generate and run a CHARMM scripts for unfold reference state calculations.
            One script is generated for each design site.
 Returns :  Nothing.
 Args    :  (1) rotamer object
            (2) energy scheme id
            (3) unfoldmer length

=cut

sub writeUnfoldDDE_old {
  GENERAL::error("This function is deprecated - not guaranteed to handle all parameters or behave correctly!");

  my $rotamer = shift;
  my $escheme = shift;
  my $unfoldmer = shift;

  GENERAL::requireArgs($rotamer, $escheme, $unfoldmer);

  # -------------------------------------------------------
  # Get the term tracker structure and setup the new entry
  my $en;
  if ($PRANK_DEF == 0) {
#    $th = GENERAL::LoadStruct($FS_DEF->{term_trackerf}, 1);
    $en = ENERGY::addTrackerTerm($TERM_TRACKER, 'unfold', 'unfold_dde', 'CHARMM (oligopeptide unfolded state): electrostatic energy (dde), VdW energy, internal torsion, H-bond energy.',
           ("van der Waals energy from MultE calculation",
            "Electrostatic energy from MultE with reduced charge",
            "Internal torsion energy of the rotamer",
            "H-bond energy from MulteE calculation"));
    # Do the following for every energy scheme this term belongs to
    ENERGY::addSchemeToTerm($en, '1', (0, 1, 2), ('unfold_sf', 'unfold_sf', 'unfold_sf'));
    ENERGY::addSchemeToTerm($en, '1h', (0, 1, 2, 3), ('unfold_sf', 'unfold_sf', 'unfold_sf', 'unfold_sf'));
    ENERGY::addSchemeToTerm($en, 'vdw', (0), ('unfold_sf'));
    ENERGY::addSchemeToTerm($en, '2', (0, 2), ('unfold_sf', 'unfold_sf'));
    ENERGY::addSchemeToTerm($en, '2b', (0, 2), ('unfold_sf', 'unfold_sf'));
    ENERGY::addSchemeToTerm($en, '2c', (0, 2), ('unfold_sf', 'unfold_sf'));
    ENERGY::addSchemeToTerm($en, 'c1', (0, 1, 2, 3), ('unfold_sf', 'unfold_sf', 'unfold_sf', 'unfold_sf'));
    ENERGY::addSchemeToTerm($en, '3', (0, 1, 2), ('unfold_sf', 'unfold_sf', 'unfold_sf'));
    ENERGY::addSchemeToTerm($en, '4', (1, 2), ('unfold_sf', 'unfold_sf'));
    ENERGY::addSchemeToTerm($en, '1e', (0, 2), ('unfold_sf', 'unfold_sf'));
    ENERGY::addSchemeToTerm($en, '1m', (0, 2), ('unfold_sf', 'unfold_sf'));
  }
  # -------------------------------------------------------

  # Do we need to run this term under the current escheme?
  if (($escheme eq '1') || ($escheme eq '1h') || ($escheme eq 'vdw') || ($escheme eq '2') || ($escheme eq '2b') ||
      ($escheme eq '2c') || ($escheme eq 'c1') || ($escheme eq '3') || ($escheme eq '4') || ($escheme eq '1e') || ($escheme eq '1m')) {
    # Divide the task among processes
    my $a = GENERAL::splitTasks(scalar(@{$rotamer->{dslist}}), $NPROC_DEF);
    my $N = $a->[$PRANK_DEF]->{end} - $a->[$PRANK_DEF]->{beg} + 1;
    if ($N > 0) {
      my $cdir = GENERAL::GetDir();
      GENERAL::cchdir($FS_DEF->{dunfold});
      my @inp = (); # will hold the names of created charmm scripts
      my $chidef = $rotamer->{chidef};

      print "Running unfolded DDE, VdW, torsion calculation (process $PRANK_DEF)...\n";
      # Create scripts
      for (my $si = 0; $si < scalar(@{$rotamer->{dslist}}); $si++) {
        next if (($si+1 < $a->[$PRANK_DEF]->{beg}) || ($si+1 > $a->[$PRANK_DEF]->{end}));
        # for each designed site, get chain id and residue number
        my $site = $rotamer->{dslist}->[$si];
        my $chain = $site->{chain};
        my $lchain = lc($chain);
        my $iresnum = $site->{iresnum};
        my $charmminp = "$LOCAL_DEF/$chain$iresnum\_multe.inp";
        push(@inp, $charmminp);

        # now write charmm input script
        # we are using charmm instead of multe.
        # since we want to turn on eef1 potential

        my $charmm = CHARMM::new($CHARMM_DEF, $charmminp);

        # using eef1 top and parm file
        $charmm->loadParm("scaled hbond");
        my $fname = $chain.$iresnum.'.crd';

        # given a crd file, setup sequence read in coordinates
        # build missing atoms and hydrogens.
        $charmm->setupFromCRD($fname);
        # find the index of the design rotamer
        my $cr;
        if ($iresnum < int($unfoldmer/2)+1) {
          $cr = $iresnum;
        } else {
          $cr = int($unfoldmer/2) + 1;
        }

        # Document the created energy table
        my $energytab = "site\_$si\_eef.tab";
        MPI_Send(File::Spec->rel2abs($energytab), 0, 0, $MPI_COMM_WORLD_DEF);
        $energytab = File::Spec->abs2rel($energytab);
        $charmm->openCard($energytab, 18);
        $charmm->defineUnfoldTemplate($chain, $cr);
        foreach my $r ( @{$site->{reslist}}) {
          my $resname = $r->{resname};
          $charmm->deleteDesignSite();
          $charmm->patchRes($resname, $chain, $cr);
          $charmm->defineUnfoldTemplate($chain, $cr);
          $charmm->defineGroup($chain, $cr);

          my $k=0;
          foreach my $rot ( @{$r->{rotamerlist}}) {
            $charmm->placeRotamer($chain, $resname, $cr, $chidef, $rot);
            $charmm->setupMulteNonBonded(1) if ($k==0);
            $charmm->inte("sele group end", "sele group .or. template end", 18, "vdw elec dihe hbon");
            $k++;
          }
          # If we just patched a proline, the amide hydrogen was removed, so re-initialize the structure (if not the last 
          if ($resname =~ /PRO/) {
            $charmm->verbose("delete atom select all end");
            $charmm->setupFromCRD($fname);
            $charmm->defineUnfoldTemplate($chain, $cr);
          }
        }
        $charmm->closeCard(18);
        $charmm->endofstory();
      }
      my $unf_multe_inp = join(" ", @inp);

      # Run the scripts
      CHARMM::runCharmm($CHARMM_DEF, $unf_multe_inp);
      
      # Do we need to delete the scripts?
      if ($FS_DEF->{del_charmm_inp}) {
        GENERAL::crm(@inp);
      }
      GENERAL::cchdir($cdir);
    }
    
    # Update and save term tracker before the actual energy file is written. This is done in
    # the interest of time (so that only one delay has to be inserted)
    push(@{$en->{files}}, File::Spec->abs2rel($FS_DEF->{unfold_multe_tabf}, $FS_DEF->{dhome}));
    if ($PRANK_DEF == 0) {
      GENERAL::SaveStruct($TERM_TRACKER, $FS_DEF->{term_trackerf});
    }
    # Wait to make sure that the energy tracker file is consistent as well as that individual
    # energy tables are consistent before we concatinate them into one
    GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF, $NFS_DELAY);
    
    # Concatinate into one energy file and document it
    if ($PRANK_DEF == 0) {
      my $ofh = GENERAL::GetOutFH($FS_DEF->{unfold_multe_tabf});
      for (my $pr = 0; $pr < $NPROC_DEF; $pr++) {
        for (my $pi = 1; $pi <= ($a->[$pr]->{end} - $a->[$pr]->{beg} + 1); $pi++) {
          my $energytab = MPI_Recv($pr, 0, $MPI_COMM_WORLD_DEF);
          my $fh = GENERAL::GetInFH($energytab);
          while (<$fh>) {
            s/[ \t]+/ /g; s/^[ \t]//; s/[ \t]$//;
            $ofh->print($_);
          }
          close($fh);
          GENERAL::crm($energytab);
        }
      }
      close($ofh);
    }
  }
}


=head2

 Title   :  writeUnfoldDDE
 Usage   :  ENERGY::writeUnfoldDDE($rotamer, $escheme);
 Function:  Generate and run a CHARMM scripts for unfold reference state calculations.
            One script is generated for each design site.
 Returns :  Nothing.
 Args    :  (1) rotamer object
            (2) energy scheme id
            (3) unfoldmer length
            (4) resume flag (optional)
            (5) custom rdie to use with DDE

=cut

sub writeUnfoldDDE {

  my $rotamer = shift;
  my $escheme = shift;
  my $pdb = shift;
  my $unfoldmer = shift;
  my $resume = shift;
  my $rdie = shift;
  if (!defined($resume)) { $resume = 0; }
  GENERAL::requireArgs($rotamer, $escheme, $unfoldmer);

  # If resuming a previous run, see if the total table has already been created
  if ($resume && GENERAL::fileCheck($FS_DEF->{unfold_multe_tabf}, "e") && (GENERAL::GetFileInfo($FS_DEF->{unfold_multe_tabf}, 'lines') == $rotamer->countSelfTotal())) {
    GENERAL::GLog("Skipping unfold DDE calculation - has been completed already...");
    return;
  }
  
  # -------------------------------------------------------
  # Get the term tracker structure and setup the new entry
  my $en;
  if ($PRANK_DEF == 0) {
    $en = ENERGY::addTrackerTerm($TERM_TRACKER, 'unfold', 'unfold_dde', "CHARMM ($unfoldmer-peptide unfolded state)",
           ("van der Waals energy from MultE calculation",
            "Electrostatic energy from MultE with reduced charge",
            "Internal torsion energy of the rotamer",
            "H-bond energy from MulteE calculation"));
    # Do the following for every energy scheme this term belongs to
    ENERGY::addSchemeToTerm($en, '1', (0, 1, 2), ('unfold_sf', 'unfold_sf', 'unfold_sf'));
    ENERGY::addSchemeToTerm($en, '1h', (0, 1, 2, 3), ('unfold_sf', 'unfold_sf', 'unfold_sf', 'unfold_sf'));
    ENERGY::addSchemeToTerm($en, '1s', (0, 1, 2, 3), ('unfold_sf', 'unfold_sf', 'unfold_sf', 'unfold_sf'));
    ENERGY::addSchemeToTerm($en, 'vdw', (0), ('unfold_sf'));
    ENERGY::addSchemeToTerm($en, 'vdws', (0), ('unfold_sf'));
    ENERGY::addSchemeToTerm($en, '2', (0, 2), ('unfold_sf', 'unfold_sf'));
    ENERGY::addSchemeToTerm($en, '2b', (0, 2), ('unfold_sf', 'unfold_sf'));
    ENERGY::addSchemeToTerm($en, '2c', (0, 2), ('unfold_sf', 'unfold_sf'));
    ENERGY::addSchemeToTerm($en, 'c1', (0, 1, 2, 3), ('unfold_sf', 'unfold_sf', 'unfold_sf', 'unfold_sf'));
    ENERGY::addSchemeToTerm($en, '3', (0, 1, 2), ('unfold_sf', 'unfold_sf', 'unfold_sf'));
    ENERGY::addSchemeToTerm($en, '4', (1, 2), ('unfold_sf', 'unfold_sf'));
    ENERGY::addSchemeToTerm($en, '1p', (2), ('unfold_sf'));
    ENERGY::addSchemeToTerm($en, 'vdwp', (2), ('unfold_sf'));
    ENERGY::addSchemeToTerm($en, '1e', (0, 2), ('unfold_sf', 'unfold_sf'));
    ENERGY::addSchemeToTerm($en, '1m', (0, 2), ('unfold_sf', 'unfold_sf'));
  }
  # -------------------------------------------------------

  # Do we need to run this term under the current escheme?
  if (($escheme eq '1') || ($escheme eq '1h') || ($escheme eq 'vdw') || ($escheme eq '2') || ($escheme eq '2b') ||
      ($escheme eq '2c') || ($escheme eq 'c1') || ($escheme eq '3') || ($escheme eq '4') || ($escheme eq '1s') || 
      ($escheme eq 'vdws') || ($escheme eq '1p') || ($escheme eq 'vdwp') || ($escheme eq '1e') || ($escheme eq '1m')) {
    my $cdir = GENERAL::GetDir();
    GENERAL::cchdir($FS_DEF->{dunfold});
    $rotamer->numify();
    # Divide the task among processes
    my $a = GENERAL::splitTasks(scalar(@{$rotamer->{dslist}}), $NPROC_DEF);
    my $N = $a->[$PRANK_DEF]->{end} - $a->[$PRANK_DEF]->{beg} + 1;
    my $crdf = "$LOCAL_DEF/polyGLY.crd";
    $pdb->writeCRD($crdf, "");

    if ($N > 0) {
      my $charmminp = "$LOCAL_DEF/unfold_multe$PRANK_DEF.inp";
      my $charmm = CHARMM::new($CHARMM_DEF, $charmminp);
      $charmm->loadParm("scaled hbond");
      $charmm->setupFromCRD($crdf);

      my $chidef = $rotamer->{chidef};
      print "Running unfolded DDE, VdW, torsion calculation (process $PRANK_DEF)...\n";
      # Create scripts
      for (my $si = 0; $si < scalar(@{$rotamer->{dslist}}); $si++) {
        next if (($si+1 < $a->[$PRANK_DEF]->{beg}) || ($si+1 > $a->[$PRANK_DEF]->{end}));
        # for each designed site, get chain id and residue number
        my $site = $rotamer->{dslist}->[$si];
        my $chain = $site->{chain};
        my $iresnum = $site->{iresnum};
        my $energytab = "site\_$si\_dde.tab";

        # if resuming a previous run, skip if this site has been completely treated
        if ($resume && GENERAL::fileCheck($energytab, "e") && (GENERAL::GetFileInfo($energytab, 'lines') == $site->{numvalidrot})) {
          print "Skipping site $si - has already been calculated in a previous run!\n";
          next;
        }

        my @sel;
        for (my $i = $iresnum - ($unfoldmer-1)/2; $i <= $iresnum + ($unfoldmer-1)/2; $i++) {
          next if ($i <= 0);
          push(@sel, ".byres. atom $chain $i N");
        }

        # Document the created energy table
#        MPI_Send(File::Spec->rel2abs($energytab), 0, 0, $MPI_COMM_WORLD_DEF);
        $charmm->openCard($energytab, 18);

        foreach my $r (@{$site->{reslist}}) {
          my $resname = $r->{resname};
          $charmm->patchRes($resname, $chain, $iresnum);
          # Need to define things after patching - things change upon patching
          $charmm->verbose("define template select (" . join(" .or. ", @sel) . ")-\n.and. (" . CHARMM::backbone() . ") end\n");

          my $k=0;
          foreach my $rot ( @{$r->{rotamerlist}}) {
            $charmm->placeRotamer($chain, $resname, $iresnum, $chidef, $rot);
            $charmm->setupMulteNonBonded(1, $rdie) if ($k==0);
            $charmm->inte("sele .byres. (atom $chain $iresnum N) .and. .not. template end", "sele .byres. (atom $chain $iresnum N) .or. template end", 18, "vdw elec dihe hbon");
            $k++;
          }
          # If we just patched a proline, the amide hydrogen was removed, so re-initialize the structure (if not the last 
          if ($resname =~ /PRO/) {
            $charmm->verbose("delete atom select all end");
            $charmm->setupFromCRD($crdf);
          }
        }
        $charmm->closeCard(18);
        # Put GLY back
        $charmm->patchRes("GLY", $chain, $iresnum);
      }
      $charmm->endofstory();

      # Do we need to delete the scripts?
      if ($FS_DEF->{del_charmm_inp}) {
        GENERAL::crm($charmminp);
      }
    }
    GENERAL::crm($crdf);

#    # Update and save term tracker before the actual energy file is written. This is done in
#    # the interest of time (so that only one delay has to be inserted)
#    if ($PRANK_DEF == 0) {
#      GENERAL::SaveStruct($TERM_TRACKER, $FS_DEF->{term_trackerf});
#    }
    # Wait to make sure that the energy tracker file is consistent as well as that individual
    # energy tables are consistent before we concatinate them into one
    GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF, $NFS_DELAY);
    
    # Concatinate into one energy file and document it
    if ($PRANK_DEF == 0) {
      my $ofh = GENERAL::GetOutFH($FS_DEF->{unfold_multe_tabf});
      for (my $si = 0; $si < scalar(@{$rotamer->{dslist}}); $si++) {
        my $energytab = "site\_$si\_dde.tab";
#      for (my $pr = 0; $pr < $NPROC_DEF; $pr++) {
#        for (my $pi = 1; $pi <= ($a->[$pr]->{end} - $a->[$pr]->{beg} + 1); $pi++) {
#        my $energytab = MPI_Recv($pr, 0, $MPI_COMM_WORLD_DEF);
        my $fh = GENERAL::GetInFH($energytab);
        while (<$fh>) {
          s/[ \t]+/ /g; s/^[ \t]//; s/[ \t]$//;
          $ofh->print($_);
        }
        close($fh);
        GENERAL::crm($energytab);
#        }
#      }
      }
      close($ofh);
      push(@{$en->{files}}, File::Spec->abs2rel($FS_DEF->{unfold_multe_tabf}, $FS_DEF->{dhome}));
    }
    GENERAL::cchdir($cdir);
  }
  # Save term tracker
  if ($PRANK_DEF == 0) {
    GENERAL::SaveStruct($TERM_TRACKER, $FS_DEF->{term_trackerf});
  }
}


sub writeRefSolv {
  my $rotamer = shift;
  my $escheme = shift;

  # -------------------------------------------------------
  # Get the term tracker structure and setup the new entry
#  my $th = GENERAL::LoadStruct($FS_DEF->{term_trackerf}, 1);
  my $ens = ENERGY::addTrackerTerm($TERM_TRACKER, 'self', 'self_sol_ref', 'Reference solvation energy determined by residue type (in folded state)',
           ("Reference solvation energy"));
  # Do the following for every energy scheme this term belongs to
  ENERGY::addSchemeToTerm($ens, '1', (0), ('1.0'));
  ENERGY::addSchemeToTerm($ens, '1h', (0), ('1.0'));
  ENERGY::addSchemeToTerm($ens, 'c1', (0), ('1.0'));

  # Reference solvation pertains to both self and unfolded energy
  my $enu = ENERGY::addTrackerTerm($TERM_TRACKER, 'unfold', 'unfold_sol_ref', 'Reference solvation energy determined by residue type (in unfolded state)',
           ("Reference solvation energy"));
  # Do the following for every energy scheme this term belongs to
  ENERGY::addSchemeToTerm($enu, '1', (0), ('1.0'));
  ENERGY::addSchemeToTerm($enu, '1h', (0), ('1.0'));
  ENERGY::addSchemeToTerm($enu, 'c1', (0), ('1.0'));
  # -------------------------------------------------------

  # Do we need to run this term under the current escheme?
  if (($escheme eq '1') || ($escheme eq '1h') || ($escheme eq 'c1')) {
    print "Writing the reference solvation energy...\n";
    my $fh = GENERAL::GetInFH(DEFINITIONS::getParmFile("reference_eef1"));
    my %solv;

    # Hash reference solvation energies by aa name
    while (<$fh>) {
      my @line=split(" ", $_);
      $solv{$line[0]}=$line[1];
    }
    close($fh);

    # Write out reference solvation energies for each rotamer
    my $outfh = GENERAL::GetOutFH($FS_DEF->{eef_solv_tabf});
    my $sitelist = $rotamer->{dslist};
    my $counter=0;
    foreach my $s (@{$sitelist}) {
      # loop over all the residues
      foreach my $r (@{$s->{reslist}}) {
        my $resname = $r->{resname};
        # for each residue, get the energies for all the  rotamers
        foreach (@{$r->{rotamerlist}}) {
          # print the energy value into the file
          if (!defined($solv{$resname})) { GENERAL::error("No reference solvation energy defined for residue $resname!"); }
          $outfh->printf("%14.5f\n", $solv{$resname});
        }
      }
    }
    close($outfh);

    # Document the reference solvation table file (for both self and unfolded energy)
    push(@{$ens->{files}}, File::Spec->abs2rel($FS_DEF->{eef_solv_tabf}, $FS_DEF->{dhome}));
    push(@{$enu->{files}}, File::Spec->abs2rel($FS_DEF->{eef_solv_tabf}, $FS_DEF->{dhome}));
  }

  # Save term tracker
  if ($PRANK_DEF == 0) {
    GENERAL::SaveStruct($TERM_TRACKER, $FS_DEF->{term_trackerf});
  }
}



=head2

 Title   :  writeSelfSA
 Usage   :  ENERGY::writeSelfSA($rotamer, $original_pdb, $escheme);
 Function:  Calculates the solvent accessible area of each rotamer in the context of
            the template. Template here is defined as the wild type structure without
            all the design side chains.
 Returns :  Nothing
 Args    :  (1) rotamer object
            (2) original PDB structure
            (3) energy scheme id
            (4) optional: an array of atom type regular expressions (will be applied to
                strings returned by PDB::atomStr())
            (5) optional: resume flag (default 0). If set, will not overwrite a previously
                computed complete table.
=cut

sub writeSelfSA {

  my $rotamer = shift;
  my $pdb = shift;
  my $escheme = shift;
  my $resume = shift;
  if (!defined($resume)) { $resume = 0; }
  GENERAL::requireArgs($rotamer, $pdb, $escheme);

  # If resuming a previous run, see if the total table has already been created
  if ($resume && GENERAL::fileCheck($FS_DEF->{self_SA_tabf}, "e") && (GENERAL::GetFileInfo($FS_DEF->{self_SA_tabf}, 'lines') == $rotamer->countSelfTotal())) {
    GENERAL::GLog("Skipping self SA calculation - has been completed already...");
    return;
  }
  
  my @types = ('_C[^_]*$', '_N[^_]*$', '_O[^_]*$', '_S[^_]*$', '_H[^_]*$');
  my $pr = 1.4;

  # -------------------------------------------------------
  # Get the term tracker structure and setup the new entry
  my $en;
  if ($PRANK_DEF == 0) {
    if (scalar(@_)) {
      @types = @_;
      # If atom types specified, we don't know how to treat them in other energy schemes.
      # For the current one, we assume uniform scaling by sa_sf.
      my @names;
      foreach my $t (@types) { push(@names, "SASA of $t in templ"); }
      foreach my $t (@types) { push(@names, "SASA of templ $t burried by rot"); }
      $en = ENERGY::addTrackerTerm($TERM_TRACKER, 'self', 'self_sa', 'Self solvent accessible surface area (scaled).', @names);
      ENERGY::addSchemeToTerm($en, $escheme, GENERAL::linspace(0, scalar(@types)-1, 1), GENERAL::ones(scalar(@types), "sa_sf"));
    } else {
      $en = ENERGY::addTrackerTerm($TERM_TRACKER, 'self', 'self_sa', 'Self solvent accessible surface area (scaled).',
            ("SASA of rot C in templ",
              "SASA of rot N in templ",
              "SASA of rot O in templ",
              "SASA of rot S in templ",
              "SASA of rot H in templ",
              "SASA of templ C burried by rot",
              "SASA of templ N burried by rot",
              "SASA of templ O burried by rot",
              "SASA of templ S burried by rot",
              "SASA of templ H burried by rot"));
      # Do the following for every energy scheme this term belongs to
      ENERGY::addSchemeToTerm($en, '2', (0,1,2,3,4,5,6,7,8,9),
      ('C_sa_sf','N_sa_sf','O_sa_sf','S_sa_sf','H_sa_sf','-C_sa_sf*psa_sf','-N_sa_sf*psa_sf','-O_sa_sf*psa_sf','-S_sa_sf*psa_sf','-H_sa_sf*psa_sf'));
      ENERGY::addSchemeToTerm($en, '2b', (0,1,2,3,4,5,6,7,8,9),
      ('C_sa_sf','N_sa_sf','O_sa_sf','S_sa_sf','H_sa_sf','-C_sa_sf*psa_sf','-N_sa_sf*psa_sf','-O_sa_sf*psa_sf','-S_sa_sf*psa_sf','-H_sa_sf*psa_sf'));
      ENERGY::addSchemeToTerm($en, '2c', (0,1,2,3,4,5,6,7,8,9),
      ('C_sa_sf','N_sa_sf','O_sa_sf','S_sa_sf','H_sa_sf','-C_sa_sf*psa_sf','-N_sa_sf*psa_sf','-O_sa_sf*psa_sf','-S_sa_sf*psa_sf','-H_sa_sf*psa_sf'));
      ENERGY::addSchemeToTerm($en, 'c1', (0,1,2,3,4,5,6,7,8,9),
      ('C_sa_sf','N_sa_sf','O_sa_sf','S_sa_sf','H_sa_sf','-C_sa_sf*psa_sf','-N_sa_sf*psa_sf','-O_sa_sf*psa_sf','-S_sa_sf*psa_sf','-H_sa_sf*psa_sf'));
      ENERGY::addSchemeToTerm($en, '3', (0,1,2,3,4,5,6,7,8,9),
      ('C_sa_sf','N_sa_sf','O_sa_sf','S_sa_sf','H_sa_sf','-C_sa_sf*psa_sf','-N_sa_sf*psa_sf','-O_sa_sf*psa_sf','-S_sa_sf*psa_sf','-H_sa_sf*psa_sf'));
      ENERGY::addSchemeToTerm($en, '4', (0,1,2,3,4,5,6,7,8,9),
      ('C_sa_sf','N_sa_sf','O_sa_sf','S_sa_sf','H_sa_sf','-C_sa_sf*psa_sf','-N_sa_sf*psa_sf','-O_sa_sf*psa_sf','-S_sa_sf*psa_sf','-H_sa_sf*psa_sf'));
    }
  }
  # -------------------------------------------------------

  # Do we need to run this term under the current escheme?
  if (($escheme eq '2') || ($escheme eq '2b') || ($escheme eq '2c') || ($escheme eq 'c1') || ($escheme eq '3') || ($escheme eq '4')) {
    my $cdir = GENERAL::GetDir();
    my $wdir = "$LOCAL_DEF/sa_self_$PRANK_DEF";
    GENERAL::cmkdir($wdir);
    GENERAL::cchdir($wdir);
    # Divide the task among processes
    my $a = GENERAL::splitTasks($rotamer->countSelfTotal(), $NPROC_DEF);
    my $N = $a->[$PRANK_DEF]->{end} - $a->[$PRANK_DEF]->{beg} + 1;
    if ($N > 0) {
      $rotamer->markRotamers($a->[$PRANK_DEF]->{beg}, $a->[$PRANK_DEF]->{end}, "visit");
      print "Running self SA calculation (process $PRANK_DEF)...\n";

      # Read in atomic vdW radii so that ENERGY::getAtomicSA does not need
      # to do this every time it is called.
      my $R = PDB::sizeLookup();
      # Create template PDB and a list of non-design atoms
      my $tpdb = $pdb->clone();
      my $rules = "";
      foreach my $site (@{$rotamer->{dslist}}) {
        my $res = $tpdb->getResByInd($site->{chain}, $site->{iresnum});
        PDB::sMutate($res, 'none'); # remove side chain
        $rules .= "n\_$site->{chain}\_$site->{iresnum} ";
      }
      $rules = GENERAL::Trim($rules);
      my @nda = $tpdb->conAtoms($rules); # a list of non-design atoms
      
      # Calculate sa's of non-design atoms (make hash table)
      my %ndasa = ();
      my @ndasa = ENERGY::getAtomicSA($tpdb, $pr, \@nda, $R);
      foreach my $atom (@nda) { $ndasa{PDB::atomStr($atom)} = shift(@ndasa); }

      # Open output file
      my $energytab = "$FS_DEF->{dself}/tmp_$PRANK_DEF";
      my $efh = GENERAL::GetOutFH($energytab);

      my $i = -1;
      foreach my $site (@{$rotamer->{dslist}}) {
        $i++;
        next if (!$site->{visit});
        my $chain = $site->{chain};
        my $iresnum = $site->{iresnum};
        my $c = $tpdb->getChain($chain);
        my $res = $c->{res}->[$iresnum-1];
        my $ores = PDB::cloneResidue($res);

        my $j = -1;
        foreach my $r ( @{$site->{reslist}}) {
          $j++;
          next if (!$r->{visit});
          print "Self SA for site $i residue $j (process $PRANK_DEF)...\n";
          # Read in all the rotamers of the current aa at the current site
          my $rotpdb = $FS_DEF->{rot_pdbf};
          $rotpdb =~ s/%%/$i/;
          $rotpdb =~ s/%/$j/;
          my $rots = PDB::new($rotpdb, "PDB", "CHARMM");
          my @rots = $rots->ConRes();

          # Place each rotamer onto the template
          foreach (my $ri=0; $ri < scalar(@{$r->{rotamerlist}}); $ri++) {
            next if (!$r->{rotamerlist}->[$ri]->{visit});
            my $rot = $rots[$ri];
            PDB::copyResidue($rot, $res);
            # calculate the solvent accessible areas (by atom) of the rotamer and of the non-design atoms
            my %bsa = (); # atoms buried by rotamer
            my @sa = ENERGY::getAtomicSA($tpdb, $pr, $res->{atom}, $R, \%bsa);
            # distribute the area by atom type - rotamer accessible area and non-design area burried by the rotamer
            my (%ra, %ndba);
            foreach my $type (@types) { $ra{$type} = 0; $ndba{$type} = 0; }
            for (my $ai = 0; $ai < scalar(@{$res->{atom}}); $ai++) {
              my $atom = $res->{atom}->[$ai];
              foreach my $type (@types) {
                if (PDB::atomStr($atom) =~ /$type/) {
                  $ra{$type} += $sa[$ai];
                  last;
                }
              }
            }
            foreach my $astr (keys(%bsa)) {
              # if this buried atom does not belong to the set of non-design atoms, skip
              next if (!defined($ndasa{$astr}));
              foreach my $type (@types) {
                if ($astr =~ /$type/) {
                  $ndba{$type} += $ndasa{$astr} - $bsa{$astr};
                  last;
                }
              }
            }
            # Save the rotamer accessible areas
            foreach my $type (@types) {
              $efh->printf("%.3f ", $ra{$type});
            }
            # Save the template areas burried by the rotamer
            foreach my $type (@types) {
              $efh->printf("%.3f ", $ndba{$type});
            }
            $efh->printf("\n");
          }
        }
        # Put the template residue back
        PDB::copyResidue($ores, $res);

        # Once we are done with this site, free the clone of the template residue
        PDB::freeResidue($ores);
      }
      close($efh);
      $rotamer->unmarkRotamers("visit");
    }

    # Wait untill all processes get to this point
    GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF, $NFS_DELAY);

    # Process 0 is to collect all data together
    if ($PRANK_DEF == 0) {
      my $k = 1; # rotamer number
      my $pr = 0; # process rank of the currently examined process
      my $pfh; # file handle of file made by the currently examined process
      # Open output file
      my $energytab = $FS_DEF->{self_SA_tabf};
      my $efh = GENERAL::GetOutFH($energytab);
      foreach my $site (@{$rotamer->{dslist}}) {
        foreach my $res (@{$site->{reslist}}) {
          foreach my $rot (@{$res->{rotamerlist}}) {
            if ($k == $a->[$pr]->{beg}) {
              $pfh = GENERAL::GetInFH("$FS_DEF->{dself}/tmp_$pr");
            }
            my $line = <$pfh>;
            if (!defined($line)) {
              GENERAL::error("Process $pr has fewer lines than expected!");
            }
            print $efh $line;
            if ($k == $a->[$pr]->{end}) {
              my $line = <$pfh>;
              if (defined($line)) {
                GENERAL::error("Process $pr has more lines than expected!");
              }
              close($pfh);
              GENERAL::crm("$FS_DEF->{dself}/tmp_$pr");
              $pr++;
            }
            $k++;
          }
        }
      }
      close($efh);
      # Document the created energy table
      push(@{$en->{files}}, File::Spec->abs2rel($energytab, $FS_DEF->{dhome}));
    }
    GENERAL::cchdir($cdir);
    GENERAL::crmdir($wdir);
  }
  if ($PRANK_DEF == 0) {
    GENERAL::SaveStruct($TERM_TRACKER, $FS_DEF->{term_trackerf});
  }
#  GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF, $NFS_DELAY);
}




=head2

 Title   :  writeSelfGB
 Usage   :  ENERGY::writeSelfGB($rotamer, $original_pdb, $escheme, $internal_epsilon);
 Function:  Calculates the polarization energy of every rotamer
            in the context of the template. The template here is
            defined as the wild type structure, except the rotamer
            in question.
 Returns :  Nothing
 Args    :  (1) rotamer object
            (2) original PDB structure
            (3) energy scheme id
            (4) internal dielectric constant
            (5) optional: resume flag (default 0). If set, will not overwrite a previously
                computed complete table.
=cut

sub writeSelfGB {

  my $rotamer = shift;
  my $pdb = shift;
  my $escheme = shift;
  my $eps_in = shift;
  my $resume = shift;
  if (!defined($resume)) { $resume = 0; }
  my %opts = @_;
  GENERAL::requireArgs($rotamer, $pdb, $escheme, $eps_in);

  # If resuming a previous run, see if the total table has already been created
  if ($resume && GENERAL::fileCheck($FS_DEF->{self_GB_tabf}, "e") && (GENERAL::GetFileInfo($FS_DEF->{self_GB_tabf}, 'lines') == $rotamer->countSelfTotal())) {
    GENERAL::GLog("Skipping self GB calculation - has been completed already...");
    return;
  }
  
  # -------------------------------------------------------
  # Get the term tracker structure and setup the new entry
  my $en;
  if ($PRANK_DEF == 0) {
    $en = ENERGY::addTrackerTerm($TERM_TRACKER, 'self', 'self_gb', "Self electrostatic energy using \'smart template GB\'.",
             ("Self atomic self polarization energies (qi^2/a)",
              "Self screening of sc-template interactions",
              "Self screening of sc-sc interactions",
              "Self Coulombic sc-template interactions",
              "Self Coulombic sc-sc interactions"));
    # Do the following for every energy scheme this term belongs to
    ENERGY::addSchemeToTerm($en, '2', (0, 1, 2, 3), ('1.0', '1.0', '1.0', '1.0'));
    ENERGY::addSchemeToTerm($en, '2b', (0, 1, 3), ('1.0', '1.0', '1.0'));
    ENERGY::addSchemeToTerm($en, '2c', (0, 1, 2, 3, 4), ('1.0', '1.0', '1.0', '1.0', '1.0'));
    ENERGY::addSchemeToTerm($en, 'gb-eef', (1, 3), ('1.0', '1.0')); # test energy function
  }
  # -------------------------------------------------------

  # caring cutoff distance (interactions beyond this, as well as the effect of anything beyond this on Born radii will be ignored)
  my $dcut = 20.0;
  # Do we need to run this term under the current escheme?
  if (($escheme eq '2') || ($escheme eq '2b') || ($escheme eq '2c') || ($escheme eq 'c1') || ($escheme eq 'gb-eef')) {
    my %br; # all Born radii for this process
    my $cdir = GENERAL::GetDir();
    my $wdir = "$LOCAL_DEF/pep_tmp_$PRANK_DEF";
    GENERAL::cmkdir($wdir);
    GENERAL::cchdir($wdir);
    
    # Divide the task among processes
    my $a = GENERAL::splitTasks($rotamer->countSelfTotal(), $NPROC_DEF);
    my $N = $a->[$PRANK_DEF]->{end} - $a->[$PRANK_DEF]->{beg} + 1;
    if ($N > 0) {
      $rotamer->numify();
      $rotamer->markRotamers($a->[$PRANK_DEF]->{beg}, $a->[$PRANK_DEF]->{end}, "visit");
      print "Running self GB calculation (process $PRANK_DEF)...\n";

      # Create template PDB we can mess with to calculate Born radii
      # (design sites occupied with WT residues) - we'll call it "shape template".
      my $stpdb;
      if (defined($opts{gb_shape_template})) { $stpdb = PDB::new($opts{gb_shape_template}); }
      else { $stpdb = $pdb->clone(); }

      # Create template PDB we can use to calculate interactions with template
      # (backbone + non-design sites - no design sites) - we'll call it "real template".
      my $rtpdb = $pdb->clone();
      foreach my $site (@{$rotamer->{dslist}}) {
        my $c = $rtpdb->getChain($site->{chain});
        my $res = $c->{res}->[$site->{iresnum}-1];
        # remove side chain
        PDB::sMutate($res, 'none');
      }

      # Read in atomic vdW radii and charges.
      my $R = PDB::sizeLookup();
      my $C = PDB::chargeLookup();
      $stpdb->charge($C);
      $rtpdb->charge($C);

      # Calculate the Born radii of real template atoms in the context of the shape template.
      my %pbr; # all Born radii calculated by this process
      if ($PRANK_DEF == 0) {
        print "Getting template radii (process $PRANK_DEF)...\n";
        %br = ENERGY::pepAtomicSolvation($stpdb, $eps_in, $R, $rtpdb->conAtoms("", 1));
        for (my $pr = 1; $pr < scalar(@$a); $pr++) {
          if ($a->[$pr]->{end} - $a->[$pr]->{beg} + 1 > 0) {
            MPI_Send(\%br, $pr, 0, $MPI_COMM_WORLD_DEF);
          }
        }
      } else {
        %br = %{MPI_Recv(0, 0, $MPI_COMM_WORLD_DEF)};
      }

      my $i = -1;
      foreach my $site (@{$rotamer->{dslist}}) {
        $i++;
        # Open output file
        my $energytab = "$FS_DEF->{dself}/tmp_$PRANK_DEF\_$i";
        if ($resume && GENERAL::fileCheck($energytab, "e") && (GENERAL::GetFileInfo($energytab, 'lines') == $site->{numrot})) {
          print "Skipping $energytab self GB calculation - has been completed already...\n";
          next;
        }
        my $efh = undef;
        next if (!$site->{visit});
        my $chain = $site->{chain};
        my $iresnum = $site->{iresnum};

        # Create working versions of the shape and real templates - we don't need the atoms really far away from the current site
        my $tmps = $stpdb;
        my $tmpr = $rtpdb;
        my $ress = $stpdb->resWithin($stpdb->getResByInd($chain, $iresnum)->{atom}, $dcut);
        my $resr = $rtpdb->resWithin($rtpdb->getResByInd($chain, $iresnum)->{atom}, $dcut);
        $stpdb->resetValidResidues(0);
        $rtpdb->resetValidResidues(0);
        foreach (@$ress) { $_->{valid} = 1; }
        foreach (@$resr) { $_->{valid} = 1; }
        $stpdb = $stpdb->clone(1);
        $rtpdb = $rtpdb->clone(1);
        
        my $sr = $stpdb->getResByInd($chain, $iresnum, 1, 2); # original residue in the shape template
        my $rr = $rtpdb->getResByInd($chain, $iresnum, 1, 2); # original residue in the real template
        
#        # Clone the original residue in the shape and real templates corresponding to this site
#        my $osr = PDB::cloneResidue($sr);
#        my $orr = PDB::cloneResidue($rr);

        my $j = -1;
        foreach my $r ( @{$site->{reslist}}) {
          $j++;
          next if (!$r->{visit});
          print "Self GB for site $i residue $j (process $PRANK_DEF)...\n";
          # Read in all the rotamers of the current aa at the current site
          my $rotpdb = $FS_DEF->{rot_pdbf};
          $rotpdb =~ s/%%/$i/;
          $rotpdb =~ s/%/$j/;
          my $rots = PDB::new($rotpdb, "PDB", "CHARMM");
          my @rots = $rots->ConRes();

          # Place each rotamer onto the template
          foreach (my $ri=0; $ri < scalar(@{$r->{rotamerlist}}); $ri++) {
            next if (!$r->{rotamerlist}->[$ri]->{visit});
            if (!defined($efh)) { $efh = GENERAL::GetOutFH($energytab); }
            my $rot = $rots[$ri];
            PDB::copyResidue($rot, $sr);
            PDB::charge($sr->{atom}, $C);
            PDB::copyResidue($rot, $rr);
            PDB::charge($rr->{atom}, $C);
            # in the rotamer PDB files, the residue number actually indicates the rotamer number
            $sr->{tag} = $sr->{resnum};
            # first, get the Born radii of the rotamer side chain atoms
            my %sbr = ENERGY::pepAtomicSolvation($stpdb, $eps_in, $R, $sr->{atom});
            # then add these radii to the hash table of all radii
            foreach my $k (keys(%sbr)) {
              if (defined($br{$k})) { GENERAL::error("A Born radius for atom $k already defined!"); }
              $br{$k} = $sbr{$k};
              if ($PRANK_DEF != 0) {
                if (defined($pbr{$k})) { GENERAL::error("A Born radius for atom $k already defined!"); }
                $pbr{$k} = $sbr{$k};
              }
            }

            # Create a list of atoms in the (current rotamer sidechain) and
            # (real template + current rotamer sidechain + the backbone corresponding to the current rotamer)
            # For backbone atoms of the current rotamer I want to use the Born radii calculated in the context
            # of this rotamer being in the site - it's probably more accurate than using Born radii derived from the
            # shape template.
            my @sa = PDB::sidechain($sr);
            my @ba = PDB::backbone($sr);
            my @ta = $rtpdb->PDB::conAtoms("n_$chain\_$iresnum");
            push(@ta, @ba);
            # calculate the polarization energy of the rotamer in question
            my $Gs = ENERGY::getGBSelf(\%br, $eps_in, \@sa);
            my $Gst = ENERGY::getGBInterAB(\%br, $eps_in, \@sa, \@ta); # specify sc atoms first for faster calculation
            my $Gss = ENERGY::getGBInterAA(\%br, $eps_in, \@sa);
            # calculate the reference state (vacuum) electrostatic energy
            my $Gvs = ENERGY::getCoulombInter($eps_in, \@ta, \@sa);
            my $Gvi = ENERGY::getCoulombTotal($eps_in, \@sa);
            # Save the values
            $efh->printf("%.3f %.3f %.3f %.3f %.3f\n", $Gs, $Gst, $Gss, $Gvs, $Gvi);
          }
        }
        # Return the original versions of the shape and real templates
        $stpdb = $tmps;
        $rtpdb = $tmpr;
#        # Put the template residue back
#        PDB::copyResidue($osr, $sr);
#        undef($sr->{tag});
#        PDB::copyResidue($orr, $rr);

#        # Once we are done with this site, free the clone of the template residue
#        PDB::freeResidue($osr);
#        PDB::freeResidue($orr);
        
        # If this process visited this site, close the file
        if (defined($efh)) { close($efh); }
      }
      if ($PRANK_DEF != 0) {
        # Send to process 0 Born radii calculated by this process
        MPI_Send(\%pbr, 0, 0, $MPI_COMM_WORLD_DEF);
      }
    }

    # Wait untill all processes get to this point
    GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF, $NFS_DELAY);
    
    # Process 0 is to collect all data together. The completion of the receives will make sure
    # that all the computing processes have finished.
    if ($PRANK_DEF == 0) {
      # First collect and save all Born radii
      for (my $pr = 1; $pr < scalar(@$a); $pr++) {
        if ($a->[$pr]->{end} - $a->[$pr]->{beg} + 1 > 0) {
          my $pbr = MPI_Recv($pr, 0, $MPI_COMM_WORLD_DEF);
          foreach my $k (keys(%$pbr)) {
            if (defined($br{$k})) { GENERAL::error("Parent process already has a Born radius for $k!"); }
            $br{$k} = $pbr->{$k};
          }
        }
      }
      GENERAL::SaveStruct(\%br, $FS_DEF->{born_radf});

      # Then, collect all energies
      my $k = 0; # rotamer number
      my $pr = 0; # process rank of the currently examined process
      my $pfh; # file handle of file made by the currently examined process
      # Open output file
      my $energytab = $FS_DEF->{self_GB_tabf};
      my $efh = GENERAL::GetOutFH($energytab);
      for (my $i = 0; $i < scalar(@{$rotamer->{dslist}}); $i++) {
        $pfh = GENERAL::GetInFH("$FS_DEF->{dself}/tmp_$pr\_$i");
        my $site = $rotamer->{dslist}->[$i];
        foreach my $res (@{$site->{reslist}}) {
          foreach my $rot (@{$res->{rotamerlist}}) {
            if ($k == $a->[$pr]->{end}) {
              my $line = <$pfh>;
              if (defined($line)) {
                GENERAL::error("File $FS_DEF->{dself}/tmp_$pr\_$i has more lines than expected!");
              }
              close($pfh);
              GENERAL::crm("$FS_DEF->{dself}/tmp_$pr\_$i");
              $pr++;
            }
            $k++;
            if ($k == $a->[$pr]->{beg}) {
              $pfh = GENERAL::GetInFH("$FS_DEF->{dself}/tmp_$pr\_$i");
            }
            my $line = <$pfh>;
            if (!defined($line)) {
              GENERAL::error("File $FS_DEF->{dself}/tmp_$pr\_$i has fewer lines than expected!");
            }
            print $efh $line;
          }
        }
        my $line = <$pfh>;
        if (defined($line)) {
          GENERAL::error("File $FS_DEF->{dself}/tmp_$pr\_$i has more lines than expected!");
        }
        close($pfh);
        GENERAL::crm("$FS_DEF->{dself}/tmp_$pr\_$i");
        close($pfh);
      }
      close($efh);
      # Document the created energy table
      push(@{$en->{files}}, File::Spec->abs2rel($energytab, $FS_DEF->{dhome}));
    }
    GENERAL::cchdir($cdir);
    GENERAL::crmdir($wdir);
  }
  if ($PRANK_DEF == 0) {
    GENERAL::SaveStruct($TERM_TRACKER, $FS_DEF->{term_trackerf});
  }
#  GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF, $NFS_DELAY);
}



=head2

 Title   :  writeSelfEEF
 Usage   :
 Function:
 Returns :
 Args    :

=cut


sub writeSelfEEF {
  my $rotamer = shift;
  my $crdf = shift;
  my $escheme = shift;
  my $resume = shift;
  if (!defined($resume)) { $resume = 0; }
  GENERAL::requireArgs($rotamer, $crdf, $escheme, $resume);

  # If resuming a previous run, see if the total table has already been created
  if ($resume && GENERAL::fileCheck($FS_DEF->{self_eef_tabf}, "e") && (GENERAL::GetFileInfo($FS_DEF->{self_eef_tabf}, 'lines') == $rotamer->countSelfTotal())) {
    GENERAL::GLog("Skipping self EEF calculation - has been completed already...");
    return;
  }

  # -------------------------------------------------------
  # Get the term tracker structure and setup the new entry
  my $en = ENERGY::addTrackerTerm($TERM_TRACKER, 'self', 'self_eef', "Self energy using CHARMM EEF",
           ("Desovlation energy from EEF",
           "Van Der Waals energy from EEF",
           "Electrostatic energy from EEF with reduced charge"));
  # Do the following for every energy scheme this term belongs to
  ENERGY::addSchemeToTerm($en, '1', (0), ('eef_sf'));
  ENERGY::addSchemeToTerm($en, '1h', (0), ('eef_sf'));
  ENERGY::addSchemeToTerm($en, '1s', (0), ('eef_sf'));
  ENERGY::addSchemeToTerm($en, 'c1', (0), ('eef_sf'));
  ENERGY::addSchemeToTerm($en, '1e', (0, 2), ('eef_sf', 'eef_sf'));
  ENERGY::addSchemeToTerm($en, '1m', (0, 2), ('eef_sf', 'eef_sf'));
#  ENERGY::addSchemeToTerm($en, '1p', (0), ('eef_sf'));
  # -------------------------------------------------------

  # Do we need to run this term under the current escheme?
  if (($escheme eq '1') || ($escheme eq '1h') || ($escheme eq 'c1') || ($escheme eq '1s') || ($escheme eq '1e') || ($escheme eq '1m')) {
    my $pdb = PDB::new($crdf);
    my $wdir = "$LOCAL_DEF/self_eef_$PRANK_DEF";
    GENERAL::cmkdir($wdir);

    my $cdir = GENERAL::GetDir();
    GENERAL::cchdir($FS_DEF->{dselfeef});

    print "Running self EEF calculation...\n";
    # ------------ First do desolvation by EEF -------------------
    my $sitelist = $rotamer->{dslist};
    my $chidef = $rotamer->{chidef};

    my $self_eef_inp = "$wdir/self_eef1.inp";
    # charmm does not understand absolute paths for output
    my $energytab = File::Spec->abs2rel($FS_DEF->{self_eef_tabf});
    # Document the created energy table
    push(@{$en->{files}}, File::Spec->abs2rel($FS_DEF->{self_eef_tabf}, $FS_DEF->{dhome}));

    my $charmm = CHARMM::new($CHARMM_DEF, $self_eef_inp);
    $charmm->loadParm('EEF1');
    $charmm->setupFromCRD($crdf);
    $charmm->openCard($energytab, 18);
    $charmm->defineTemplate($rotamer);

    foreach my $site (@{$sitelist}) {
      # for each designed site, get chain id and residue number
      my $chain = $site->{chain};
      my $iresnum = $site->{iresnum};
      my $ores = $pdb->getResByInd($chain, $iresnum);
      foreach my $r (@{$site->{reslist}}) {
        my $resname = $r->{resname};
        $charmm->verbose("! -- VISITING RESIDUE $resname...");
        $charmm->deleteDesignSite();
        $charmm->verbose("bomlev -2");
        $charmm->patchRes($resname, $chain, $iresnum);
        $charmm->defineTemplate($rotamer);
        $charmm->defineGroup($chain, $iresnum);

        my $k=0;
        foreach my $rot ( @{$r->{rotamerlist}} ) {
          #  reference $angle is a list of list
          $charmm->placeRotamer($chain, $resname, $iresnum, $chidef, $rot);
          if ($k == 0) {
            if ($escheme eq "1m") { $charmm->setupIMM1(); }
            else { $charmm->setupEEF1(); }
          }
          $charmm->inte('sele group end', 'sele template end', 18, "asp vdw elec");
          $k++;   # keep track of the number of rotamers
        }
        # - if we just patched a Pro, the amide hydrogen may have been removed
        # - if we just patched a Gly, atom type of CA may have changed
        # - if there used to be a GLY, atom type of CA may have changed
        # so re-initialize the template structure
        if (($resname =~ /PRO|GLY/) || ($ores->{resname} =~ /GLY/)) {
          $charmm->verbose("delete atom select all end");
          $charmm->setupFromCRD($crdf);
          $charmm->defineTemplate($rotamer);
        }
      }
    }
    $charmm->closeCard(18);
    $charmm->endofstory();

    # Do we need to delete the script
    if ($FS_DEF->{del_charmm_inp}) {
      GENERAL::crm("$self_eef_inp");
      GENERAL::crmdir($wdir);
    }

    GENERAL::cchdir($cdir);
  }

  # Save term tracker
  GENERAL::SaveStruct($TERM_TRACKER, $FS_DEF->{term_trackerf});
}


sub writeCrystaSelflEEFInput {

    my $rotamer = shift;
    my $fname=shift;
    # my @energy=@_;

    GENERAL::cmkdir("self/crystal");
    GENERAL::csystem ("cp $fname self/crystal/");
    GENERAL::cchdir("self/crystal");

    my $sitelist=$rotamer->{dslist};
    my $chidef=$rotamer->{chidef};

    my $fname1= 'crystalself_eef1.inp';
    my $energytab= 'crystalself_eef1.tab';

    my $charmm=CHARMM::new($CHARMM_DEF, $fname1);
    $charmm->loadParm('EEF1');
    $charmm->setupFromCRD($fname);
    $charmm->openCard($energytab, 18);
    $charmm->defineTemplate($rotamer);


    foreach my $site (@{$sitelist}) {

        # for each designed site, get chain id and residue number

        my $chain=$site->{chain};
        my $iresnum=$site->{iresnum};

        foreach my $r ( @{$site->{reslist}}) {

            my $resname=$r->{resname};
            $charmm->deleteDesignSite();
            $charmm->verbose("bomlev -2");
            $charmm->patchRes($resname, $chain, $iresnum);
            $charmm->defineTemplate($rotamer);
            $charmm->defineGroup($chain,$iresnum);

            my $k=0;
            foreach my $rot ( @{$r->{rotamerlist}} ) {
                #  reference $angle is a list of list
                my $chi =$rot->{chi};

                $charmm->rebuildRotamer($chidef, $chi, $resname, $chain, $iresnum);
                $charmm->setupEEF1() if ($k==0);
                $charmm->getEnergy();
                $charmm->selectINTE('sele group end', 'sele template end');
                #$charmm->calculateTotal();
                $charmm->writeCard(18, '* ?asp ?vdw ?elec');

                $k++;   # keep track of the number of rotamers
            }

        }
    }

    $charmm->closeCard(18);
    $charmm->endofstory();

    GENERAL::cchdir("../");
    GENERAL::cchdir("../");
}


=head2

 Title   :  writeSelfDDE
 Usage   :
 Function:
 Returns :
 Args    : 1. ROTAMER object corresponding to the design problem
           2. CRD file of structure
           3. Energy scheme
           4. Resume flag - optional.
           5. Custom rdie to use with DDE - optional.

=cut

sub writeSelfDDE {

  my $rotamer = shift;
  my $crdf = shift;
  my $escheme = shift;
  my $resume = shift;
  my $rdie = shift;
  if (!defined($resume)) { $resume = 0; }
  GENERAL::requireArgs($rotamer, $crdf, $escheme, $resume);

  # If resuming a previous run, see if the total table has already been created
  if ($resume && GENERAL::fileCheck($FS_DEF->{self_multe_tabf}, "e") && (GENERAL::GetFileInfo($FS_DEF->{self_multe_tabf}, 'lines') == $rotamer->countSelfTotal())) {
    GENERAL::GLog("Skipping self DDE calculation - has been completed already...");
    return;
  }
  
  # -------------------------------------------------------
  # Get the term tracker structure and setup the new entry
  my $en = ENERGY::addTrackerTerm($TERM_TRACKER, 'self', 'self_dde', 'CHARMM (interaction with template + self)',
           ("van der Waals energy from MultE calculation",
            "Electrostatic energy from MultE with reduced charge",
            "Internal torsion energy of the rotamer",
            "H-bond energy from MulteE calculation"));
  # Do the following for every energy scheme this term belongs to
  ENERGY::addSchemeToTerm($en, '1', (0, 1, 2), ('vdw_sf', 'dde_self_sf', 1.0));
  ENERGY::addSchemeToTerm($en, '1h', (0, 1, 2, 3), ('vdw_sf', 'dde_self_sf', 1.0, 'hbond_sf'));
  ENERGY::addSchemeToTerm($en, 'vdw', (0), ('vdw_sf'));
  ENERGY::addSchemeToTerm($en, '2', (0, 2), ('vdw_sf', 1.0));
  ENERGY::addSchemeToTerm($en, '2b', (0, 2), ('vdw_sf', 1.0));
  ENERGY::addSchemeToTerm($en, '2c', (0, 2), ('vdw_sf', 1.0));
  ENERGY::addSchemeToTerm($en, 'c1', (0, 1, 2, 3), ('vdw_sf', 'dde_self_sf', 1.0, 'hbond_sf'));
  ENERGY::addSchemeToTerm($en, '3', (0, 1, 2), ('vdw_sf', 'dde_self_sf', 1.0));
  ENERGY::addSchemeToTerm($en, '4', (1, 2), ('dde_self_sf', 1.0));
  ENERGY::addSchemeToTerm($en, '1p', (2), (1.0));
  ENERGY::addSchemeToTerm($en, 'gb-eef', (2), (1.0)); # test energy function
  ENERGY::addSchemeToTerm($en, 'vdwp', (2), (1.0));
  ENERGY::addSchemeToTerm($en, '1e', (0, 2), ('vdw_sf', 1.0));
  ENERGY::addSchemeToTerm($en, '1m', (0, 2), ('vdw_sf', 1.0));
  # -------------------------------------------------------

  # Do we need to run this term under the current escheme?
  if (($escheme eq '1') || ($escheme eq '1h') || ($escheme eq 'vdw') || ($escheme eq '2') || 
      ($escheme eq '2b') || ($escheme eq '2c') || ($escheme eq 'c1') || ($escheme eq '3') || ($escheme eq '4') || ($escheme eq '1p') ||
      ($escheme eq 'gb-eef') || ($escheme eq 'vdwp') || ($escheme eq '1e') || ($escheme eq '1m')) {
    my $pdb = PDB::new($crdf);
    my $wdir = "$LOCAL_DEF/self_dde_$PRANK_DEF";
    GENERAL::cmkdir($wdir);

    my $cdir = GENERAL::GetDir();
    GENERAL::cchdir($FS_DEF->{dselfmulte});

    print "Running self DDE, VdW, torsion calculation...\n";
    my $sitelist=$rotamer->{dslist};
    my $chidef=$rotamer->{chidef};

    my $self_multe_inp = "$wdir/self_multe.inp";
    # charmm does not understand absolute paths for output
    my $energytab = File::Spec->abs2rel($FS_DEF->{self_multe_tabf});
    # Document the created energy table
    push(@{$en->{files}}, File::Spec->abs2rel($FS_DEF->{self_multe_tabf}, $FS_DEF->{dhome}));

    my $charmm = CHARMM::new($CHARMM_DEF, $self_multe_inp);
    $charmm->loadParm("scaled hbond");
    $charmm->setupFromCRD($crdf);
    $charmm->openCard($energytab, 18);
    $charmm->defineTemplate($rotamer);

    my $i=-1;
    foreach my $site (@{$sitelist}) {

      # for each designed site, get chain id and residue number
      $i++;
      my $chain = $site->{chain};
      my $iresnum = $site->{iresnum};
      my $ores = $pdb->getResByInd($chain, $iresnum);

      my $j=-1;
      foreach my $r ( @{$site->{reslist}}) {

        $j++;
        my $resname = $r->{resname};
        $charmm->verbose("bomlev -2");
        $charmm->patchRes($resname, $chain, $iresnum);
        $charmm->defineTemplate($rotamer);
        $charmm->defineGroup($chain,$iresnum);
        my $k=0;
        foreach my $rot ( @{$r->{rotamerlist}} ) {
          $charmm->placeRotamer($chain, $resname, $iresnum, $chidef, $rot);
          $charmm->setupMulteNonBonded(1, $rdie) if ($k==0);
          $charmm->inte('sele group end', 'sele group .or. template end', 18, "vdw elec dihe hbon");
          $k++;   # keep track of the number of rotamers
        }
        # - if we just patched a Pro, the amide hydrogen may have been removed
        # - if we just patched a Gly, atom type of may have CA changed
        # - if there used to be a GLY, atom type of CA may have changed
        # so re-initialize the template structure
        if (($resname =~ /PRO|GLY/) || ($ores->{resname} =~ /GLY/)) {
          $charmm->verbose("delete atom select all end");
          $charmm->setupFromCRD($crdf);
          $charmm->defineTemplate($rotamer);
        }
      }
    }
    $charmm->closeCard(18);
    $charmm->endofstory();
    # Do we need to delete the script?
    if ($FS_DEF->{del_charmm_inp}) {
      GENERAL::crm("$self_multe_inp");
      GENERAL::crmdir($wdir);
    }

    GENERAL::cchdir($cdir);
  }

  # Save term tracker
  GENERAL::SaveStruct($TERM_TRACKER, $FS_DEF->{term_trackerf});
}



=head2

 Title   :  writeSelfDDE_min
 Usage   :
 Function:  For testing purposes, evaluates selv DDE and vdW energies after a short step of minimization in CHARMM.
 Returns :
 Args    : 1. ROTAMER object corresponding to the design problem
           2. CRD file of structure
           3. Energy scheme
           4. Resume flag - optional.
           5. Custom rdie to use with DDE - optional.

=cut


sub writeSelfDDE_min {

  my $rotamer = shift;
  my $crdf = shift;
  my $escheme = shift;
  my $resume = shift;
  my $rdie = shift;
  if (!defined($resume)) { $resume = 0; }
  GENERAL::requireArgs($rotamer, $crdf, $escheme, $resume);

  # If resuming a previous run, see if the total table has already been created
  if ($resume && GENERAL::fileCheck($FS_DEF->{self_multe_tabf}, "e") && (GENERAL::GetFileInfo($FS_DEF->{self_multe_tabf}, 'lines') == $rotamer->countSelfTotal())) {
    GENERAL::GLog("Skipping self DDE calculation - has been completed already...");
    return;
  }
  
  # -------------------------------------------------------
  # Get the term tracker structure and setup the new entry
  my $en = ENERGY::addTrackerTerm($TERM_TRACKER, 'self', 'self_dde_min', 'CHARMM (interaction with template + self minimized)',
           ("van der Waals energy from MultE calculation",
            "Electrostatic energy from MultE with reduced charge",
            "Internal torsion energy of the rotamer",
            "H-bond energy from MulteE calculation"));
  # Do the following for every energy scheme this term belongs to
  ENERGY::addSchemeToTerm($en, '1s', (0, 1, 2), ('vdw_sf', 'dde_self_sf', 1.0));
  ENERGY::addSchemeToTerm($en, 'vdws', (0), ('vdw_sf'));
  # -------------------------------------------------------

  # Do we need to run this term under the current escheme?
  if (($escheme eq '1s') || ($escheme eq 'vdws')) {
    my $pdb = PDB::new($crdf);
    my $wdir = "$LOCAL_DEF/self_dde_min$PRANK_DEF";
    my @terms = ("vdw", "elec",  "dihe", "hbon");
    GENERAL::cmkdir($wdir);

    my $cdir = GENERAL::GetDir();
    GENERAL::cchdir($FS_DEF->{dselfmulte});

    print "Running self DDE, VdW, torsion calculation...\n";
    my $sitelist = $rotamer->{dslist};
    my $chidef = $rotamer->{chidef};

    my $self_multe_inp = "$wdir/self_multe.inp";
    # charmm does not understand absolute paths for output
    my $energytab = File::Spec->abs2rel(GENERAL::GetBase($FS_DEF->{self_multe_tabf}) . "_min.tab");
    # Document the created energy table
    push(@{$en->{files}}, File::Spec->abs2rel($energytab, $FS_DEF->{dhome}));

    my $charmm = CHARMM::new($CHARMM_DEF, $self_multe_inp);
    $charmm->loadParm("scaled hbond");
    $charmm->setupFromCRD($crdf);
    $charmm->send("skipe all excl bond excl angle excl dihe excl impr " . join(" excld ", @terms) . "\n");
    $charmm->openCard($energytab, 18);
    $charmm->defineTemplate($rotamer);

    my $i=-1;
    foreach my $site (@{$sitelist}) {

      # for each designed site, get chain id and residue number
      $i++;
      my $chain = $site->{chain};
      my $iresnum = $site->{iresnum};
      my $ores = $pdb->getResByInd($chain, $iresnum);

      my $j=-1;
      foreach my $r ( @{$site->{reslist}}) {

        $j++;
        my $resname = $r->{resname};
        $charmm->verbose("bomlev -2");
        $charmm->patchRes($resname, $chain, $iresnum);
        $charmm->defineTemplate($rotamer);
        $charmm->defineGroup($chain,$iresnum);
        my $k=0;
        foreach my $rot ( @{$r->{rotamerlist}} ) {
          $charmm->placeRotamer($chain, $resname, $iresnum, $chidef, $rot);
          $charmm->setupMulteNonBonded(1, $rdie) if ($k==0);
          $charmm->send("cons fix sele template end\n");
          $charmm->send("mini sd nstep 5 step 0.1 ihbfrq 1 inbfrq 1\n");
          $charmm->inte('sele group end', 'sele group .or. template end', 18, join(" ", @terms));
          $k++;   # keep track of the number of rotamers
        }
        # - if we just patched a Pro, the amide hydrogen may have been removed
        # - if we just patched a Gly, atom type of may have CA changed
        # - if there used to be a GLY, atom type of CA may have changed
        # so re-initialize the template structure
        if (($resname =~ /PRO|GLY/) || ($ores->{resname} =~ /GLY/)) {
          $charmm->verbose("delete atom select all end");
          $charmm->setupFromCRD($crdf);
          $charmm->defineTemplate($rotamer);
        }
      }
    }
    $charmm->closeCard(18);
    $charmm->endofstory();

    # Do we need to delete the script?
    if ($FS_DEF->{del_charmm_inp}) {
      GENERAL::crm("$self_multe_inp");
      GENERAL::crmdir($wdir);
    }

    GENERAL::cchdir($cdir);
  }

  # Save term tracker
  GENERAL::SaveStruct($TERM_TRACKER, $FS_DEF->{term_trackerf});
}


=head2

 Title   :  writeSelfEPIC
 Usage   :
 Function:  Calculates and write all energy terms required by the given energy scheme that EPIC can handle.
 Returns :
 Args    : 1. ROTAMER object corresponding to the design problem
           2. PDB object
           3. Energy scheme
           4. Resume flag - optional.

=cut

sub writeSelfEPIC {

  my $rotamer = shift;
  my $pdb = shift;
  my $escheme = shift;
  my $parm = shift;
  my $resume = shift;
  if (!defined($resume)) { $resume = 0; }
  GENERAL::requireArgs($rotamer, $pdb, $escheme, $parm, $resume);
  my $impdbs = shift; $impdbs = () if (!defined($impdbs));
  my $terpatch = shift;

  # -------------------------------------------------------
  # Get the term tracker structure and setup the new entry
  my %en;
  if ($PRANK_DEF == 0) {
    # all terms EPIC can calculate
    $en{VDW} = ENERGY::addTrackerTerm($TERM_TRACKER, 'self', 'self_vdw_epic', 'Self vdW energies using EPIC', "i-t + i-i vdW interaction energy");
    ENERGY::addSchemeToTerm($en{VDW}, '1p', (0), ('vdw_sf'));
    ENERGY::addSchemeToTerm($en{VDW}, 'gb-eef', (0), ('vdw_sf')); # test energy function
    ENERGY::addSchemeToTerm($en{VDW}, 'vdwp', (0), ('vdw_sf'));
    $en{EEF1} = ENERGY::addTrackerTerm($TERM_TRACKER, 'self', 'self_eef1_epic', 'Self EEF1 energies using EPIC', "i-t EEF1 desolvation energy");
    ENERGY::addSchemeToTerm($en{EEF1}, '1p', (0), ('eef_sf'));
    ENERGY::addSchemeToTerm($en{EEF1}, 'gb-eef', (0), ('eef_sf')); # test energy function
    $en{DDE} = ENERGY::addTrackerTerm($TERM_TRACKER, 'self', 'self_dde_epic', 'Self DDE energies using EPIC', "i-t + i-i DDE interaction energy");
    ENERGY::addSchemeToTerm($en{DDE}, '1p', (0), ('dde_pair_sf'));
  }
  # -------------------------------------------------------

  # Do we need to run this term under the current escheme?
  my $time = time(); my @terms;
  if (($escheme eq '1p') || ($escheme eq 'gb-eef') || ($escheme eq 'vdwp')) {
    if ($escheme eq '1p') { @terms = ('VDW', 'EEF1', 'DDE'); }
    if ($escheme eq 'gb-eef') { @terms = ('VDW', 'EEF1'); }
    if ($escheme eq 'vdwp') { @terms = ('VDW'); }
    my $epic = undef;
    $rotamer->numify();
    # create template
    my $tpdb = $pdb->clone();
    my $ctpdb = defined($parm->{ct}) ? PDB::new($parm->{ct}) : undef; # common template, if there is one
    foreach my $ds (@{$rotamer->{dslist}}) {
      my $res = PDB::getResByInd($tpdb, $ds->{chain}, $ds->{iresnum});
      PDB::sMutate($res, "NONE");
    }
    # number template atoms
    my @ta = $tpdb->conAtoms();
    my @cta = defined($parm->{ct}) ? $ctpdb->conAtoms() : (); # common template, if there is one;
    for (my $i = 0; $i < scalar(@ta); $i++) { $ta[$i]->{ai} = $i+1; }
    my $tfile = "$FS_DEF->{drots}/template.out"; my $ctfile = "$FS_DEF->{drots}/commontemplate.out"; my @imtfiles;
    for (my $im = 0; $im < scalar(@$impdbs); $im++) { push(@imtfiles, GENERAL::GetBase($tfile) . ".$im.out"); }
    # write template coordinates
    if ($PRANK_DEF == 0) {
      my $fh = GENERAL::GetOutFH($tfile); my @imfhs;
      for (my $im = 0; $im < scalar(@$impdbs); $im++) { push(@imfhs, GENERAL::GetOutFH($imtfiles[$im])); }
      $fh->printf("%d 1\n", scalar(@ta));
      for (my $im = 0; $im < scalar(@$impdbs); $im++) { $imfhs[$im]->printf("%d 1\n", scalar(@ta)); }
      foreach my $a (@ta) {
        $fh->printf("%f %f %f\n", $a->{xcoor}, $a->{ycoor}, $a->{zcoor});
        for (my $im = 0; $im < scalar(@$impdbs); $im++) {
          my $ia = $impdbs->[$im]->getAtom($a->{residue}->{chain}->{id}, $a->{residue}->{iresnum}, $a->{atomname});
          $imfhs[$im]->printf("%f %f %f\n", $ia->{xcoor}, $ia->{ycoor}, $ia->{zcoor});
        }
      }
      close($fh);
      for (my $im = 0; $im < scalar(@$impdbs); $im++) { close($imfhs[$im]); }
      $tpdb->writePDB("$FS_DEF->{drots}/template.pdb", "");

      # if there is a common template (one for all images), write it too
      $ctpdb->writeDUMMY($ctfile) if (defined($parm->{ct}));
    }
    GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF); # wait until template gets written

    # Divide the task among processes
    my $Nt = $rotamer->countResTotal();
    my $a = GENERAL::splitTasks($Nt, $NPROC_DEF);
    my $N = $a->[$PRANK_DEF]->{end} - $a->[$PRANK_DEF]->{beg} + 1;
    # Calculate each term
    foreach my $term (@terms) {
      print "Process $PRANK_DEF: Calculating self '$term' using EPIC...\n";
      # Have we done this before?
      my $tablef;
      if ($term eq "VDW") { $tablef = "$FS_DEF->{dself}/vdw_epic.tab"; }
      elsif ($term eq "DDE") { $tablef = "$FS_DEF->{dself}/dde_epic.tab"; }
      elsif ($term eq "EEF1") { $tablef = "$FS_DEF->{dself}/eef1_epic.tab"; }
      else { GENERAL::error("Do not know what to call final table for term '$term'"); }
      # If resuming a previous run, see if the total table has already been created
      if ($resume && GENERAL::fileCheck($tablef, "e") && (GENERAL::GetFileInfo($tablef, 'lines') == $rotamer->countSelfTotal())) {
        GENERAL::GLog("Skipping self EPIC '$term' calculation - has been completed already...");
        print "Skipping self EPIC '$term' calculation - has been completed already...\n";
        next;
      }
      # Read appropriate topology and parameter tables
      my ($P, $T);
      if ($term =~ /EEF/) {
        my $Peef = CHARMM::readEEFParameters(DEFINITIONS::getParmFile("solvpar.inp"));
        $T = CHARMM::readCHARMMTopology(DEFINITIONS::getParmFile("toph19_eef1.inp"));
        $T = CHARMM::readCHARMMTopology(DEFINITIONS::getParmFile("patchtop19_eef1.inp"), $T);
        $P = CHARMM::readCHARMMParameters(DEFINITIONS::getParmFile("param19_eef1.inp"));
        CHARMM::expandCHARMMParameters($P, $T);
        foreach my $at (keys(%$Peef)) {
          foreach my $val (keys(%{$Peef->{$at}})) {
            $P->{eef1}{$at}{$val} = $Peef->{$at}{$val};
          }
        }
      } elsif ($term =~ /(VDW|DDE)/) {
        $T = CHARMM::readCHARMMTopology(DEFINITIONS::getParmFile("toph19.inp"));
        $T = CHARMM::readCHARMMTopology(DEFINITIONS::getParmFile("patchtop19.inp"), $T);
        $P = CHARMM::readCHARMMParameters(DEFINITIONS::getParmFile("param19_90.inp"));
        CHARMM::expandCHARMMParameters($P, $T);
      }
      # Do the assigned part of the work
      my $num = 0;
      if ($N > 0) {
        # patch the template if terminal patches are specified
        if (defined($terpatch)) {
          foreach my $chain (@{$tpdb->{chain}}) {
            my $fres = $chain->{res}->[0]; $fres->{patch} = $terpatch->{$chain->{id}}{fir};
            my $lres = $chain->{res}->[-1]; $lres->{patch} = $terpatch->{$chain->{id}}{lst};
            PDB::createPatchedResidue($T, $fres->{resname}, $fres->{patch}, undef, 0);
            PDB::createPatchedResidue($T, $lres->{resname}, $lres->{patch}, undef, 0);
          }
        }
        # loop over site and residue
        for (my $si = 0; $si < scalar(@{$rotamer->{dslist}}); $si++) {
          my $ds = $rotamer->{dslist}->[$si];
          my $res = PDB::getResByInd($tpdb, $ds->{chain}, $ds->{iresnum});
          my $ores = PDB::cloneResidue($res);
          for (my $ri = 0; $ri < scalar(@{$ds->{reslist}}); $ri++) {
            $num++;
            next if (($num < $a->[$PRANK_DEF]->{beg}) || ($num > $a->[$PRANK_DEF]->{end}));
            my $dresi = $ds->{reslist}->[$ri];
            my $rpdbi = $FS_DEF->{rot_pdbf}; $rpdbi =~ s/%%/$si/; $rpdbi =~ s/%/$ri/;
            my $inpi = $FS_DEF->{rot_outf}; $inpi =~ s/%%/$si/; $inpi =~ s/%/$ri/;
            my @iminpi; for (my $im = 0; $im < scalar(@$impdbs); $im++) { my $iminpi = $FS_DEF->{rot_ioutf}; $iminpi =~ s/%%%/$im/; $iminpi =~ s/%%/$si/; $iminpi =~ s/%/$ri/; push(@iminpi, $iminpi); }
            # read one rotamer for this residue
            my $rotsi = PDB::new($rpdbi, "PDB", "CHARMM", 1);
            my $roti = $rotsi->{chain}->[0]->{res}->[0];
            # paste the sidechain into the template and number sidechain
            PDB::replaceSidechain($roti, $res);
            my @sa = PDB::sidechain($res);
            for (my $i = 0; $i < scalar(@sa); $i++) { $sa[$i]->{ai} = $i+1; }
            # get exclusion lists for rot-template and rot-rot interactions
            my ($ex12, $ex13, $ex14, $ex12s, $ex13s, $ex14s) = PDB::exclusionList(\@sa, "1-2 1-3 1-4", $T, 3);
#printf("12:\n"); PDB::printExclusionList($ex12);
#printf("12s:\n"); PDB::printExclusionList($ex12s);
#printf("13:\n"); PDB::printExclusionList($ex13);
#printf("13s:\n"); PDB::printExclusionList($ex13s);
#printf("14:\n"); PDB::printExclusionList($ex14);
#printf("14s:\n"); PDB::printExclusionList($ex14s);
            # PRO hack (to avoid interactions with the amide H, which is supposed to be absent from the template)
            if ($roti->{resname} eq "PRO") {
              if ($parm->{param} ne "19") {
                GENERAL::error("Don't know how to apply a PRO-specific exclusion list modification in parameter set $parm->{param}");
              }
              # 1. add (everything in side chain - H) into the exclusion list (hack)
              my $H = PDB::getAtomInRes($res, "H", -1);
              if ($H ne -1) {
                foreach my $sa (@sa) { my @arr = ($sa, $H); push(@$ex12, \@arr); }
                # 2. temporarily add H into the topology of PRO (put in a group by itself)
                $T->{residues}->{PRO}->{atoms}->{H}->{type} = "H";
                $T->{residues}->{PRO}->{atoms}->{H}->{group} = scalar(keys(%{$T->{residues}->{PRO}->{groups}})) + 1;
                $T->{residues}->{PRO}->{atoms}->{H}->{charge} = 0.0;
              }
            } elsif ($roti->{resname} =~ /PHE|TYR/) {
              # PHE and TYR hacks (to avoid in-ring 1-4 interactions - no interactions for these in CHARMM)
              if ($parm->{param} ne "19") {
                GENERAL::error("Don't know how to apply a PHE/TYR-specific exclusion list modification in parameter set $parm->{param}");
              }
              my @nex = ("CG", "CZ", "CD1", "CE2", "CE1", "CD2");
              for (my $i = 0; $i < scalar(@nex); $i += 2) {
                my @arr = (PDB::getAtomInRes($res, $nex[$i]), PDB::getAtomInRes($res, $nex[$i+1])); push(@$ex12s, \@arr);
              }
            } elsif ($roti->{resname} eq "TRP") {
              # TRP hack (to avoid in-ring/across-rings interactions, 1-4 or other - no interactions for these in CHARMM)
              my @nex = ("CG","CZ3", "CG","CZ2", "CD2","CH2", "CE2","CZ3", "CE3","CD1", "CE3","CZ2", "CE3","NE1", "CD1","CZ2", "NE1","CH2",
                         "CG","CH2", "CD1","CZ3", "CD1","CH2", "CZ3","NE1");
              for (my $i = 0; $i < scalar(@nex); $i += 2) {
                my @arr = (PDB::getAtomInRes($res, $nex[$i]), PDB::getAtomInRes($res, $nex[$i+1])); push(@$ex12s, \@arr);
              }
            }
            # calculate rot-template intereaction
            $epic = ENERGY::interfaceEPIC("$FS_DEF->{dself}/$si-$ri-rt.tab", $term, "pair", $T, $P, $parm, \@sa, $inpi, 'atomsj', \@ta, 'inpj', $tfile, 'sf', 1.0, 'ex12', $ex12, 'ex13', $ex13, 'ex14', $ex14, 'prog', $epic);
            $epic = ENERGY::interfaceEPIC("$FS_DEF->{dself}/$si-$ri-rct.tab", $term, "pair", $T, $P, $parm, \@sa, $inpi, 'atomsj', \@cta, 'inpj', $ctfile, 'sf', 1.0, 'prog', $epic) if (defined($parm->{ct}));
            for (my $im = 0; $im < scalar(@$impdbs); $im++) {
              $epic = ENERGY::interfaceEPIC("$FS_DEF->{dself}/$si-$ri-rt.$im.tab", $term, "pair", $T, $P, $parm, \@sa, $inpi, 'atomsj', \@ta, 'inpj', $imtfiles[$im], 'sf', 1.0, 'prog', $epic);
              $epic = ENERGY::interfaceEPIC("$FS_DEF->{dself}/$si-$ri-rr.$im.tab", $term, "corr", $T, $P, $parm, \@sa, $inpi, 'atomsj', \@sa, 'inpj', $iminpi[$im], 'sf', 1.0, 'prog', $epic);
            }
            # calculate rot-rot interaction
            if ($term ne "EEF1") {
              $epic = ENERGY::interfaceEPIC("$FS_DEF->{dself}/$si-$ri-rr.tab", $term, "self", $T, $P, $parm, \@sa, $inpi, 'sf', 1.0, 'ex12', $ex12s, 'ex13', $ex13s, 'ex14', $ex14s, 'prog', $epic);
            }
            if (time() - $time > 30*$LOCAL_TIMEOUT) { GENERAL::touchLocalSpace(); $time = time(); }
            # undo PRO hack
            if ($roti->{resname} eq "PRO") {
              delete($T->{residues}->{PRO}->{atoms}->{H});
            }
          }
          PDB::replaceSidechain($ores, $res);
        }
      }
      if (defined($epic)) { $epic->close(); $epic = undef; } # close to cause a flush

      # Wait until all processes are done before putting the data together
      GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF, $NFS_DELAY);
      # Combine all tables into one sorted by rotamer
      if ($PRANK_DEF == 0) {
        my $ofh = GENERAL::GetOutFH($tablef);
        for (my $si = 0; $si < scalar(@{$rotamer->{dslist}}); $si++) {
          my $ds = $rotamer->{dslist}->[$si];
          for (my $ri = 0; $ri < scalar(@{$ds->{reslist}}); $ri++) {
            # calculate the contribution from the images
            my @imself;
            for (my $im = 0; $im < scalar(@$impdbs); $im++) {
              my $rt = GENERAL::file2array("$FS_DEF->{dself}/$si-$ri-rt.$im.tab");
              my $rr = GENERAL::file2array("$FS_DEF->{dself}/$si-$ri-rr.$im.tab");
              if (scalar(@$rr) != scalar(@$rt)) { GENERAL::error(sprintf("image self terms: %d r-r energies, but %d r-t energies!", scalar(@$rr), scalar(@$rt))); }
              @imself = GENERAL::ones(scalar(@$rt), 0) if (scalar(@imself) == 0);
              for (my $i = 0; $i < scalar(@$rt); $i++) {
                $imself[$i] += ($rt->[$i] + $rr->[$i])/2;
              }
              GENERAL::crm("$FS_DEF->{dself}/$si-$ri-rt.$im.tab");
              GENERAL::crm("$FS_DEF->{dself}/$si-$ri-rr.$im.tab");
            }
            # primary contribution
            if ($term eq "EEF1") {
              my $rt = GENERAL::file2array("$FS_DEF->{dself}/$si-$ri-rt.tab");
              my $rct = undef; $rct = GENERAL::file2array("$FS_DEF->{dself}/$si-$ri-rct.tab") if (defined($parm->{ct}));
              for (my $i = 0; $i < scalar(@$rt); $i++) {
                $ofh->printf("%e\n", $rt->[$i] + ((scalar(@$impdbs) == 0) ? 0 : $imself[$i]) + ((defined($parm->{ct})) ? $rct->[$i] : 0));
              }
              GENERAL::crm("$FS_DEF->{dself}/$si-$ri-rt.tab");
              GENERAL::crm("$FS_DEF->{dself}/$si-$ri-rct.tab") if (defined($parm->{ct}));
            } else {
              my $rr = GENERAL::file2array("$FS_DEF->{dself}/$si-$ri-rr.tab");
              my $rt = GENERAL::file2array("$FS_DEF->{dself}/$si-$ri-rt.tab");
              my $rct = undef; $rct = GENERAL::file2array("$FS_DEF->{dself}/$si-$ri-rct.tab") if (defined($parm->{ct}));
              if (scalar(@$rr) != scalar(@$rt)) { GENERAL::error(sprintf("%d r-r energies, but %d r-t energies!", scalar(@$rr), scalar(@$rt))); }
              if ((defined($parm->{ct})) && (scalar(@$rct) != scalar(@$rt))) { GENERAL::error(sprintf("%d r-t energies, but %d r-ct energies!", scalar(@$rt), scalar(@$rct))); }
              if ((scalar(@$impdbs) > 0) && (scalar(@imself) != scalar(@$rt))) { GENERAL::error(sprintf("%d r-t energies, but %d image self energies!", scalar(@$rt), scalar(@imself))); }
              for (my $i = 0; $i < scalar(@$rr); $i++) {
                $ofh->printf("%e\n", $rr->[$i] + $rt->[$i] + ((scalar(@$impdbs) == 0) ? 0 : $imself[$i]) + ((defined($parm->{ct})) ? $rct->[$i] : 0));
              }
              GENERAL::crm("$FS_DEF->{dself}/$si-$ri-rr.tab");
              GENERAL::crm("$FS_DEF->{dself}/$si-$ri-rt.tab");
              GENERAL::crm("$FS_DEF->{dself}/$si-$ri-rct.tab") if (defined($parm->{ct}));
            }
          }
        }
        close($ofh);
        # Document the created energy table
        push(@{$en{$term}->{files}}, File::Spec->abs2rel($tablef, $FS_DEF->{dhome}));
      }
      # Wait until all old files are erased before moving on
      GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF);
    }
  }

  # Save term tracker
  if ($PRANK_DEF == 0) {
    GENERAL::SaveStruct($TERM_TRACKER, $FS_DEF->{term_trackerf});
  }
}


=head2

 Title   :  writePairwiseSA
 Usage   :  ENERGY::writePairwiseSA($rotamer, $pdb, $escheme);
 Function:  Calculates the solvent accessible area of every rotamer at
            every site buried by every rotamer at every other site. This
            requires that the self surface area terms has been calculated
            with writeSelfSA (i.e. the solvent accessible area buried by
            the template).
 Returns :  Nothing
 Args    :  (1) rotamer object
            (2) original PDB structure
            (3) energy scheme
            (4) optional: an array of atom type regular expressions (will be applied to
                strings returned by PDB::atomStr())
=cut

sub writePairwiseSA {

  my $rotamer = shift;
  my $pdb = shift;
  my $escheme = shift;
  GENERAL::requireArgs($rotamer, $pdb, $escheme);
  my @types = ('_C[^_]*$', '_N[^_]*$', '_O[^_]*$', '_S[^_]*$', '_H[^_]*$');
  my $pr = 1.4; # probe radius
  # -------------------------------------------------------
  # Get the term tracker structure and setup the new entry
  my $en;
  if ($PRANK_DEF == 0) {
    if (scalar(@_)) {
      @types = @_;
      # If atom types specified, we don't know how to treat them in other energy schemes.
      # For the current one, we assume uniform scaling by sa_sf.
      my @names;
      foreach my $t (@types) { push(@names, "Pairwise SASA of $t"); }
      $en = ENERGY::addTrackerTerm($TERM_TRACKER, 'pair', 'pair_sa', 'Pair solvent accessible surface area (scaled).', @names);
      ENERGY::addSchemeToTerm($en, $escheme, GENERAL::linspace(0, scalar(@types)-1, 1), GENERAL::ones(scalar(@types), "sa_sf"));
    } else {
      $en = ENERGY::addTrackerTerm($TERM_TRACKER, 'pair', 'pair_sa', 'Pair solvent accessible surface area (scaled).',
              ("Pairwise SASA of C",
              "Pairwise SASA of N",
              "Pairwise SASA of O",
              "Pairwise SASA of S",
              "Pairwise SASA of H"));
      # Do the following for every energy scheme this term belongs to
      ENERGY::addSchemeToTerm($en, '2', (0, 1, 2, 3, 4), 
      ('-C_sa_sf*psa_sf', '-N_sa_sf*psa_sf', '-O_sa_sf*psa_sf', '-S_sa_sf*psa_sf', '-H_sa_sf*psa_sf'));
      ENERGY::addSchemeToTerm($en, '2b', (0, 1, 2, 3, 4), 
      ('-C_sa_sf*psa_sf', '-N_sa_sf*psa_sf', '-O_sa_sf*psa_sf', '-S_sa_sf*psa_sf', '-H_sa_sf*psa_sf'));
      ENERGY::addSchemeToTerm($en, '2c', (0, 1, 2, 3, 4), 
      ('-C_sa_sf*psa_sf', '-N_sa_sf*psa_sf', '-O_sa_sf*psa_sf', '-S_sa_sf*psa_sf', '-H_sa_sf*psa_sf'));
      ENERGY::addSchemeToTerm($en, 'c1', (0, 1, 2, 3, 4), 
      ('-C_sa_sf*psa_sf', '-N_sa_sf*psa_sf', '-O_sa_sf*psa_sf', '-S_sa_sf*psa_sf', '-H_sa_sf*psa_sf'));
      ENERGY::addSchemeToTerm($en, '3', (0, 1, 2, 3, 4), 
      ('-C_sa_sf*psa_sf', '-N_sa_sf*psa_sf', '-O_sa_sf*psa_sf', '-S_sa_sf*psa_sf', '-H_sa_sf*psa_sf'));
      ENERGY::addSchemeToTerm($en, '4', (0, 1, 2, 3, 4), 
      ('-C_sa_sf*psa_sf', '-N_sa_sf*psa_sf', '-O_sa_sf*psa_sf', '-S_sa_sf*psa_sf', '-H_sa_sf*psa_sf'));
    }
  }
  # -------------------------------------------------------
  # Do we need to run this term under the current escheme?
  if (($escheme eq '2') || ($escheme eq '2b') || ($escheme eq '2c') || ($escheme eq 'c1') || ($escheme eq '3') || ($escheme eq '4')) {
    my $cdir = GENERAL::GetDir();
    my $wdir = "$LOCAL_DEF/sa_pair_$PRANK_DEF";
    GENERAL::cmkdir($wdir);
    GENERAL::cchdir($wdir);
    
    # Divide the task among processes
    my $a = GENERAL::splitTasks($rotamer->countValidPairTotal(), $NPROC_DEF);
    my $N = $a->[$PRANK_DEF]->{end} - $a->[$PRANK_DEF]->{beg} + 1;
    if ($N > 0) {
      $rotamer->markValidRotamerPairs($a->[$PRANK_DEF]->{beg}, $a->[$PRANK_DEF]->{end}, "beg1", "beg2", "end1", "end2");

      # First, read the self term (will be needed later)
      my $stab = ETABLE2::readEnergyTerm($TERM_TRACKER->{self}->{self_sa});

      # Index the rotamers so that we know which entry in the self table
      # each rotamer corresponds to
      $rotamer->indexRotamers('ind');

      # We'll need atomic radii to for distance calculations
      my $R = PDB::sizeLookup();

      # Then read in all rotamers of all residues at all sites
      my @R;
      for (my $i = 0; $i < scalar(@{$rotamer->{dslist}}); $i++) {
        my $site = $rotamer->{dslist}->[$i];
        my @r;
        for (my $j = 0; $j < scalar(@{$site->{reslist}}); $j++) {
          # Read in all the rotamers of the current aa at the current site
          my $rotpdb = $FS_DEF->{rot_pdbf};
          $rotpdb =~ s/%%/$i/;
          $rotpdb =~ s/%/$j/;
          my $rots = PDB::new($rotpdb, "PDB", "CHARMM");
          push(@r, $rots);
        }
        push(@R, \@r);
      }
      # Create template PDB
      my $tpdb = $pdb->clone();
      foreach my $site (@{$rotamer->{dslist}}) {
        my $c = $tpdb->getChain($site->{chain});
        my $res = $c->{res}->[$site->{iresnum}-1];
        # remove side chain
        PDB::sMutate($res, 'none');
      }

      # Open output file
      my $energytab = "$FS_DEF->{dpair}/tmp_$PRANK_DEF";
      my $efh = GENERAL::GetOutFH($energytab);
      # Go through all pairs of rotamers designated for this process
      my $first = 1; my $last = 0;
      for (my $i1=0; $i1 < scalar(@{$rotamer->{dslist}})-1; $i1++) {
        my $site1 = $rotamer->{dslist}->[$i1];
        next if ($first && !$site1->{beg1});
        my $chain1 = $site1->{chain};
        my $iresnum1 = $site1->{iresnum};
        my $c1 = $tpdb->getChain($chain1);
        my $res1 = $c1->{res}->[$iresnum1-1];
        my $ores1 = PDB::cloneResidue($res1);

        for (my $j1 = 0; $j1 < scalar(@{$site1->{reslist}}); $j1++) {
          my $re1 = $site1->{reslist}->[$j1];
          next if ($first && !$re1->{beg1});
          # Read in all the rotamers of the current aa at the current site
          my $rots1 = $R[$i1][$j1];

          # Place each valid rotamer onto the template
          for (my $k1 = 0; $k1 < scalar(@{$re1->{rotamerlist}}); $k1++) {
            my $ro1 = $re1->{rotamerlist}->[$k1];
            next if ($first && !$ro1->{beg1});
            next if ($ro1->{valid} == 0);
            print "Pair SA for site $i1 residue $j1 rotamer $k1 - ... (process $PRANK_DEF)\n";

            my $rot1 = $rots1->{chain}->[0]->{res}->[$k1];
            PDB::copyResidue($rot1, $res1);

            for (my $i2 = $i1+1; $i2 < scalar(@{$rotamer->{dslist}}); $i2++) {
              my $site2 = $rotamer->{dslist}->[$i2];
              next if ($first && !$site2->{beg2});
              my $chain2 = $site2->{chain};
              my $iresnum2 = $site2->{iresnum};
              my $c2 = $tpdb->getChain($chain2);
              my $res2 = $c2->{res}->[$iresnum2-1];
              my $ores2 = PDB::cloneResidue($res2);

              for (my $j2 = 0; $j2 < scalar(@{$site2->{reslist}}); $j2++) {
                my $re2 = $site2->{reslist}->[$j2];
                next if ($first && !$re2->{beg2});
                # Read in all the rotamers of the current aa at the current site
                my $rots2 = $R[$i2][$j2];

                # Place each valid rotamer onto the template
                for (my $k2 = 0; $k2 < scalar(@{$re2->{rotamerlist}}); $k2++) {
                  my $ro2 = $re2->{rotamerlist}->[$k2];
                  next if ($first && !$ro2->{beg2});
                  next if ($ro2->{valid} == 0);
                  $first = 0;

                  my $rot2 = $rots2->{chain}->[0]->{res}->[$k2];

                  # Are the two rotamers too far apart?
                  if (PDB::closest($rot2, $rot1, $R) >= 2*$pr) {
                    # These two rotamers do not "obstract" one another from solvent
                    foreach my $type (@types) {
                      $efh->printf("0 ");
                    }
                    $efh->printf("\n");
                  } else {
                    # If not, need to do the calculation
                    PDB::copyResidue($rot2, $res2);

                    # calculate the solvent accessible area (by atom) of the two rotamers
                    my @atms = (@{$res1->{atom}}, @{$res2->{atom}});
                    my @sa = ENERGY::getAtomicSA($tpdb, $pr, \@atms, $R);
                    if (scalar(@sa) != scalar(@atms)) {
                      GENERAL::error(scalar(@atms) . " atoms in" . PDB::resStr($res1) . " and " . PDB::resStr($res2) .
                                     ", but " . scalar(@sa) . " atom areas returned by naccess!");
                    }
                    # distribute the area by atom type
                    my %ta;
                    foreach my $type (@types) { $ta{$type} = 0; }
                    for (my $ai = 0; $ai < scalar(@atms); $ai++) {
                      my $atom = $atms[$ai];
                      for (my $ti = 0; $ti < scalar(@types); $ti++) {
                        my $type = $types[$ti];
                        if (PDB::atomStr($atom) =~ /$type/) {
                          $ta{$type} += $sa[$ai];
                          last;
                        }
                      }
                    }
                    # Save the areas
                    for (my $ti = 0; $ti < scalar(@types); $ti++) {
                      my $type = $types[$ti];
                      # what we want is the area mutually buried by the pair
                      my $a = $stab->[$ro1->{ind}-1]->[$ti] + $stab->[$ro2->{ind}-1]->[$ti] - $ta{$type};
                      $efh->printf("%.3f ", $a);
                    }
                    $efh->printf("\n");
                  }
                  if ($site1->{end1} && $re1->{end1} && $ro1->{end1} && $site2->{end2} && $re2->{end2} && $ro2->{end2}) {
                    $last = 1;
                    last;
                  }
                }
                last if ($last);
              }
              # Put the template residue back
              PDB::copyResidue($ores2, $res2);

              # Once we are done with this site, free the clone of the template residue
              PDB::freeResidue($ores2);

              last if ($last);
            }
            last if ($last);
          }
          last if ($last);
        }
        # Put the template residue back
        PDB::copyResidue($ores1, $res1);

        # Once we are done with this site, free the clone of the template residue
        PDB::freeResidue($ores1);
        last if ($last);
      }
      close($efh);
      $rotamer->unmarkRotamers('beg1', 'beg2', 'end1', 'end2', 'ind');
    }

    # Wait untill all processes get to this point
    GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF, $NFS_DELAY);

    # Process 0 is to collect all data together
    if ($PRANK_DEF == 0) {
      my $k = 1; # rotamer pair number
      my $pr = 0; # process rank of the currently examined process
      my $pfh; # file handle of file made by the currently examined process
      # Open output file
      my $energytab = $FS_DEF->{pair_SA_tabf};
      my $efh = GENERAL::GetOutFH($energytab, "bzip2");
      for (my $i1=0; $i1 < scalar(@{$rotamer->{dslist}}) - 1; $i1++) {
        my $site1 = $rotamer->{dslist}->[$i1];
        foreach my $res1 (@{$site1->{reslist}}) {
          foreach my $rot1 (@{$res1->{rotamerlist}}) {
            next if (!$rot1->{valid});
            for (my $i2=$i1+1; $i2 < scalar(@{$rotamer->{dslist}}); $i2++) {
              my $site2 = $rotamer->{dslist}->[$i2];
              foreach my $res2 (@{$site2->{reslist}}) {
                foreach my $rot2 (@{$res2->{rotamerlist}}) {
                  next if (!$rot2->{valid});
                  if ($k == $a->[$pr]->{beg}) {
                    $pfh = GENERAL::GetInFH("$FS_DEF->{dpair}/tmp_$pr");
                  }
                  my $line = <$pfh>;
                  if (!defined($line)) {
                    GENERAL::error("Process $pr has fewer lines than expected ($a->[$pr]->{beg} to $a->[$pr]->{end})!");
                  }
                  print $efh $line;
                  if ($k == $a->[$pr]->{end}) {
                    my $line = <$pfh>;
                    if (defined($line)) {
                      GENERAL::error("Process $pr has more lines than expected ($a->[$pr]->{beg} to $a->[$pr]->{end})!");
                    }
                    close($pfh);
                    GENERAL::crm("$FS_DEF->{dpair}/tmp_$pr");
                    $pr++;
                  }
                  $k++;
                }
              }
            }
          }
        }
      }
      # Document the created energy table
      push(@{$en->{files}}, File::Spec->abs2rel($energytab, $FS_DEF->{dhome}));
    }
    GENERAL::cchdir($cdir);
    GENERAL::crmdir($wdir);
  }
  if ($PRANK_DEF == 0) {
    GENERAL::SaveStruct($TERM_TRACKER, $FS_DEF->{term_trackerf});
  }
#  GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF, $NFS_DELAY);
}


=head2

 Title   :  writePairwiseGB
 Usage   :  ENERGY::writePairwiseGB($rotamer, $escheme, $epsilon);
 Function:  Calculates the polarization energy associated with the interaction
            of every pair of rotamers in the context of the template.
 Returns :  Nothing
 Args    :  (1) rotamer object
            (2) energy scheme
            (3) internal dielectric
=cut

sub writePairwiseGB {

  my $rotamer = shift;
  my $escheme = shift;
  my $eps_in = shift;
  GENERAL::requireArgs($rotamer, $escheme, $eps_in);

  # -------------------------------------------------------
  # Get the term tracker structure and setup the new entry
#  my $th = GENERAL::LoadStruct($FS_DEF->{term_trackerf}, 1);
  my $en = ENERGY::addTrackerTerm($TERM_TRACKER, 'pair', 'pair_gb', "Pair electrostatic energy using \'smart template GB\'.",
           ("Pair GB Polarization Energy",
            "Pair \'quasi-vacuum\' Electrostatic Energy"));
  # Do the following for every energy scheme this term belongs to
  ENERGY::addSchemeToTerm($en, '2', (0, 1), ('1.0', '1.0'));
  ENERGY::addSchemeToTerm($en, '2b', (0, 1), ('1.0', '1.0'));
  ENERGY::addSchemeToTerm($en, '2c', (0, 1), ('1.0', '1.0'));
  ENERGY::addSchemeToTerm($en, 'c1', (0, 1), ('1.0', '1.0'));
  ENERGY::addSchemeToTerm($en, 'gb-eef', (0, 1), ('1.0', '1.0')); # test energy function
  # -------------------------------------------------------

  # Do we need to run this term under the current escheme?
  if (($escheme eq '2') || ($escheme eq '2b') || ($escheme eq '2c') || ($escheme eq 'c1') || ($escheme eq 'gb-eef')) {
    my $time = time();
    my $cdir = GENERAL::GetDir();
    # Divide the task among processes
    my $a = GENERAL::splitTasks($rotamer->countValidPairTotal(), $NPROC_DEF);
    my $N = $a->[$PRANK_DEF]->{end} - $a->[$PRANK_DEF]->{beg} + 1;
    if ($N > 0) {
      $rotamer->markValidRotamerPairs($a->[$PRANK_DEF]->{beg}, $a->[$PRANK_DEF]->{end}, "beg1", "beg2", "end1", "end2");

      print "Running pair GB calculation...\n";

      # First, read in Born radii and charge file
      my %br = GENERAL::LoadStruct($FS_DEF->{born_radf});
      my $C = PDB::chargeLookup();

      # Then read in all rotamers of all residues at all sites and charge them
      my @R;
      for (my $i = 0; $i < scalar(@{$rotamer->{dslist}}); $i++) {
        my $site = $rotamer->{dslist}->[$i];
        my @r;
        for (my $j = 0; $j < scalar(@{$site->{reslist}}); $j++) {
          # Read in all the rotamers of the current aa at the current site
          my $rotpdb = $FS_DEF->{rot_pdbf};
          $rotpdb =~ s/%%/$i/;
          $rotpdb =~ s/%/$j/;
          my $rots = PDB::new($rotpdb, "PDB", "CHARMM");
          PDB::charge($rots, $C);
          push(@r, $rots);
        }
        push(@R, \@r);
      }

      # Open output file
      my $energytab = "$FS_DEF->{dpair}/tmp_$PRANK_DEF";
      my $efh = GENERAL::GetOutFH($energytab);
      # Go through all pairs of rotamers designated for this process
      my $first = 1; my $last = 0;
      for (my $i1=0; $i1 < scalar(@{$rotamer->{dslist}})-1; $i1++) {
        my $site1 = $rotamer->{dslist}->[$i1];
        next if ($first && !$site1->{beg1});

        for (my $j1 = 0; $j1 < scalar(@{$site1->{reslist}}); $j1++) {
          my $re1 = $site1->{reslist}->[$j1];
          next if ($first && !$re1->{beg1});
          my $rots1 = $R[$i1][$j1];

          for (my $k1 = 0; $k1 < scalar(@{$re1->{rotamerlist}}); $k1++) {
            my $ro1 = $re1->{rotamerlist}->[$k1];
            next if ($first && !$ro1->{beg1});
            next if ($ro1->{valid} == 0);

            print "Pair GB for site $i1 residue $j1 rotamer $k1 - ...(process $PRANK_DEF)\n";
            my $rot1 = $rots1->{chain}->[0]->{res}->[$k1];
            # Tag the rotamer with the rotamer number and "place" it onto the correct site
            $rot1->{tag} = $rot1->{resnum};
            $rot1->{iresnum} = $site1->{iresnum};
            my @sa1 = PDB::sidechain($rot1);

            for (my $i2 = $i1+1; $i2 < scalar(@{$rotamer->{dslist}}); $i2++) {
              my $site2 = $rotamer->{dslist}->[$i2];
              next if ($first && !$site2->{beg2});

              for (my $j2 = 0; $j2 < scalar(@{$site2->{reslist}}); $j2++) {
                my $re2 = $site2->{reslist}->[$j2];
                next if ($first && !$re2->{beg2});
                my $rots2 = $R[$i2][$j2];

                for (my $k2 = 0; $k2 < scalar(@{$re2->{rotamerlist}}); $k2++) {
                  my $ro2 = $re2->{rotamerlist}->[$k2];
                  next if ($first && !$ro2->{beg2});
                  next if ($ro2->{valid} == 0);
                  $first = 0;

                  my $rot2 = $rots2->{chain}->[0]->{res}->[$k2];
                  # Tag the rotamer with the rotamer number and "place" it onto the correct site
                  $rot2->{tag} = $rot2->{resnum};
                  $rot2->{iresnum} = $site2->{iresnum};
                  my @sa2 = PDB::sidechain($rot2);

                  # Calculate the solvation energy associated with the interaction between rot1 and rot2
                  my $Gpol = ENERGY::getGBInterAB(\%br, $eps_in, \@sa1, \@sa2);
                  # Calculate the electrostatic energy of interaction between rot1 and rot2 in the reference state
                  my $Gv = ENERGY::getCoulombInter($eps_in, \@sa1, \@sa2);
                  $efh->printf("%.3f %.3f\n", $Gpol, $Gv);

                  if ($site1->{end1} && $re1->{end1} && $ro1->{end1} && $site2->{end2} && $re2->{end2} && $ro2->{end2}) {
                    $last = 1;
                    last;
                  }
                }
                last if ($last);
              }
              last if ($last);
              if (time() - $time > 30*$LOCAL_TIMEOUT) { GENERAL::touchLocalSpace(); $time = time(); }
            }
            last if ($last);
          }
          last if ($last);
        }
        last if ($last);
      }
      close($efh);
      $rotamer->unmarkRotamers('beg1', 'beg2', 'end1', 'end2', 'ind');
    }


    # Wait untill all processes get to this point
    GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF, $NFS_DELAY);

    # Process 0 is to collect all data together
    if ($PRANK_DEF == 0) {
      my $k = 1; # rotamer pair number
      my $pr = 0; # process rank of the currently examined process
      my $pfh; # file handle of file made by the currently examined process
      # Open output file
      my $energytab = $FS_DEF->{pair_GB_tabf};
      my $efh = GENERAL::GetOutFH($energytab, "bzip2");
      for (my $i1=0; $i1 < scalar(@{$rotamer->{dslist}}) - 1; $i1++) {
        my $site1 = $rotamer->{dslist}->[$i1];
        foreach my $res1 (@{$site1->{reslist}}) {
          foreach my $rot1 (@{$res1->{rotamerlist}}) {
            next if (!$rot1->{valid});
            for (my $i2=$i1+1; $i2 < scalar(@{$rotamer->{dslist}}); $i2++) {
              my $site2 = $rotamer->{dslist}->[$i2];
              foreach my $res2 (@{$site2->{reslist}}) {
                foreach my $rot2 (@{$res2->{rotamerlist}}) {
                  next if (!$rot2->{valid});
                  if ($k == $a->[$pr]->{beg}) {
                    $pfh = GENERAL::GetInFH("$FS_DEF->{dpair}/tmp_$pr");
                  }
                  my $line = <$pfh>;
                  if (!defined($line)) {
                    GENERAL::error("Process $pr has fewer lines than expected ($a->[$pr]->{beg} to $a->[$pr]->{end})!");
                  }
                  print $efh $line;
                  if ($k == $a->[$pr]->{end}) {
                    my $line = <$pfh>;
                    if (defined($line)) {
                      GENERAL::error("Process $pr has more lines than expected ($a->[$pr]->{beg} to $a->[$pr]->{end})!");
                    }
                    close($pfh);
                    GENERAL::crm("$FS_DEF->{dpair}/tmp_$pr");
                    $pr++;
                  }
                  $k++;
                }
              }
            }
          }
        }
      }
      # Document the created energy table
      push(@{$en->{files}}, File::Spec->abs2rel($energytab, $FS_DEF->{dhome}));
    }
  }
  if ($PRANK_DEF == 0) {
    GENERAL::SaveStruct($TERM_TRACKER, $FS_DEF->{term_trackerf});
  }
#  GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF, $NFS_DELAY);
}



=head2

 Title   :  writePairwiseEEF
 Usage   :
 Function:
 Returns :
 Args    :

=cut


sub writePairwiseEEF {

  my $rotamer = shift;
  my $escheme = shift;
  GENERAL::requireArgs($rotamer, $escheme);
  my $resume = shift;
  if (!defined($resume)) { $resume = 0; }
  
  # If resuming a previous run, see if the total table has already been created
  if ($resume && GENERAL::fileCheck($FS_DEF->{pair_eef_tot_tabf}, "e") && (GENERAL::GetFileInfo($FS_DEF->{pair_eef_tot_tabf}, 'lines') == $rotamer->countValidPairTotal())) {
    GENERAL::GLog("Skipping pair EEF calculation - has been completed already...");
    return;
  }
  
  # -------------------------------------------------------
  # Get the term tracker structure and setup the new entry
  my $en;
  if ($PRANK_DEF == 0) {
    $en = ENERGY::addTrackerTerm($TERM_TRACKER, 'pair', 'pair_eef', 'Pairwise energy using CHARMM EEF1',
           ("Desovlation energy from EEF",
           "Van Der Waals energy from EEF",
           "Electrostatic energy from EEF with reduced charge"));
    # Do the following for every energy scheme this term belongs to
    ENERGY::addSchemeToTerm($en, '1', (0), ('eef_sf'));
    ENERGY::addSchemeToTerm($en, '1h', (0), ('eef_sf'));
    ENERGY::addSchemeToTerm($en, 'c1', (0), ('eef_sf'));
    ENERGY::addSchemeToTerm($en, '1e', (0, 2), ('eef_sf', 'eef_sf'));
    ENERGY::addSchemeToTerm($en, '1m', (0, 2), ('eef_sf', 'eef_sf'));
  }
  # -------------------------------------------------------

  # Do we need to run this term under the current escheme?
  if (($escheme eq '1') || ($escheme eq '1h') || ($escheme eq 'c1') || ($escheme eq '1e') || ($escheme eq '1m')) {
    my @terms = ("asp", "vdw", "elec");
    my $ctonnb = 7.0; my $ctofnb = 9.0; my $cutnb = 10.0; # CHARMM distance cutoffs
    my $dt = 10**10;
    $dt = 10*60 if ($NPROC_DEF > 1); # time interval, in seconds, for tasks split re-negotiation
    my $T = time(); # current time in seconds
    $rotamer->numify();
    $rotamer->extent(); # calculate the maxmimum extent of each side chain and get the coordinates of the base atom
    my $sitelist = $rotamer->{dslist};
    my $cdir = GENERAL::GetDir();
    my $wdir = "$LOCAL_DEF/eef_pair_$PRANK_DEF";
    GENERAL::cmkdir($wdir);
    # We want to make the number of pair files in a directory to grow linearly with
    # the number of design residues and not quadratically, so we create another layers of directories.
    if ($PRANK_DEF == 0) {
      for (my $s = 0; $s < scalar(@{$rotamer->{dslist}}); $s++) {
        for (my $res = 0; $res < scalar(@{$rotamer->{dslist}->[$s]->{reslist}}); $res++) {
          GENERAL::cmkdir("$FS_DEF->{dpair_eef}/$s\_$res");
        }
      }
    }
    MPI_Barrier($MPI_COMM_WORLD_DEF);
    print "Running pair EEF calculation (process $PRANK_DEF)...\n";
    my $Nt = $rotamer->countResPairTotal();
    # Divide the task among processes
    my $a = GENERAL::splitTasksNonContig($Nt, $NPROC_DEF);
    my @pa = @{$a->[$PRANK_DEF]}; # asignment for this process
    my @comp; # a list of completed tasks
    while (scalar(@pa) > 0) {
      my $pairnum = 1; # index of the pair
      my $sitelist = $rotamer->{dslist};
      my $chidef = $rotamer->{chidef};
      # loop over all the pairwise term
      for (my $s1 = 0; $s1 < scalar(@$sitelist); $s1++) {
        my $site1 = $sitelist->[$s1];
        my $chain1 = $sitelist->[$s1]->{chain};
        my $iresnum1 = $sitelist->[$s1]->{iresnum};
        for (my $res1 = 0; $res1 < scalar(@{$sitelist->[$s1]->{reslist}}); $res1++) {
          my $resname1 = $sitelist->[$s1]->{reslist}->[$res1]->{resname};
          my $re1 = $sitelist->[$s1]->{reslist}->[$res1];
          # Read all rotamers of this residue at this site
          my $rotpdb1 = $FS_DEF->{rot_pdbf};
          $rotpdb1 =~ s/%%/$s1/;
          $rotpdb1 =~ s/%/$res1/;
          my $rots1 = PDB::new($rotpdb1, "PDB", "CHARMM");
          my @rots1 = $rots1->ConRes();

          for (my $s2 = $s1+1; $s2 < scalar(@{$sitelist}); $s2++) {
            my $chain2 = $sitelist->[$s2]->{chain};
            my $iresnum2 = $sitelist->[$s2]->{iresnum};
            for (my $res2 = 0; $res2 < scalar(@{$sitelist->[$s2]->{reslist}}); $res2++) {
              last if ((scalar(@pa) == 0) || (time() - $T > $dt));
              if ($pairnum != $pa[0]) { $pairnum++; next; }
              my $resname2 = $sitelist->[$s2]->{reslist}->[$res2]->{resname};
              my $re2 = $sitelist->[$s2]->{reslist}->[$res2];
              
              my $energytab = "$FS_DEF->{dpair_eef}/$s1\_$res1/paireef$pairnum.tab";
              $energytab =~ s/\%/$pairnum/;
              $energytab = File::Spec->abs2rel($energytab);

              # if resuming a previous run, skip if this pair of residues has been completely treated
              if ($resume && GENERAL::fileCheck($energytab, "e") && (GENERAL::GetFileInfo($energytab, 'lines') == 
                       $sitelist->[$s1]->{reslist}->[$res1]->{numvalidrot}*$sitelist->[$s2]->{reslist}->[$res2]->{numvalidrot})) {
                print "Skipping pair $pairnum - has already been calculated in a previous run!\n";
                $pairnum++; push(@comp, shift @pa);
                next;
              }

              # See if the two residues are close enough to bother calling CHARMM
              my $d = sqrt(($re1->{base_coor}->[0]-$re2->{base_coor}->[0])**2 +
                           ($re1->{base_coor}->[1]-$re2->{base_coor}->[1])**2 +
                           ($re1->{base_coor}->[2]-$re2->{base_coor}->[2])**2);
              if ($d - $re1->{extent} - $re2->{extent} > $cutnb) {
                #print "Process $PRANK_DEF: pair $pairnum too distant - filling with zeros...\n";
                my $fh = GENERAL::GetOutFH($energytab);
                my $np = ($re1->{numvalidrot})*(($re2->{numvalidrot}));
                for (my $i = 1; $i <= $np; $i++) {
                  $fh->printf(("0 "x(scalar(@terms)) . "\n"));
                }
                close($fh);
                $pairnum++; push(@comp, shift @pa);
                next;
              }

              # Read all rotamers of this residue at this site
              my $rotpdb2 = $FS_DEF->{rot_pdbf};
              $rotpdb2 =~ s/%%/$s2/;
              $rotpdb2 =~ s/%/$res2/;
              my $rots2 = PDB::new($rotpdb2, "PDB", "CHARMM");
              my @rots2 = $rots2->ConRes();

              my $charmminp = "$wdir/paireef".$pairnum.".inp";
              my $charmm = CHARMM::new($CHARMM_DEF, $charmminp);
              $charmm->loadPairParm('EEF1');
              $charmm->openCard($energytab, 18);
              for (my $rot1 = 0; $rot1<=$#{$sitelist->[$s1]->{reslist}->[$res1]->{rotamerlist}}; $rot1++) {
                next if ($sitelist->[$s1]->{reslist}->[$res1]->{rotamerlist}->[$rot1]->{valid} == 0);
                my $rotamer1 = $rots1[$rot1];
                $charmm->setPairSequence($resname1, 'A');
                $charmm->setRotCoord("A", 1, $rotamer1);
                for(my $rot2 = 0; $rot2<=$#{$sitelist->[$s2]->{reslist}->[$res2]->{rotamerlist}}; $rot2++) {
                  next if ($sitelist->[$s2]->{reslist}->[$res2]->{rotamerlist}->[$rot2]->{valid} == 0);
                  my $rotamer2 = $rots2[$rot2];
                  $charmm->setPairSequence($resname2, 'B');
                  $charmm->setRotCoord("B", 1, $rotamer2);
                  if ($escheme eq "1m") { $charmm->setupIMM1(); }
                  else { $charmm->setupEEF1(1, $ctonnb, $ctofnb, $cutnb); }
#                  $charmm->getEnergy();
#                  $charmm->selectINTE('sele segid A end', 'sele segid B end');
#                  $charmm->writeCard(18, '* ?asp ?vdw ?elec');
                  $charmm->inte('sele segid A end', 'sele segid B end', 18, join(" ", @terms));
                  $charmm->verbose("\ndele atom select SEGID B end\n");
                }
                $charmm->verbose("\ndele atom select SEGID A end\n");
              }
              $charmm->closeCard(18);
              $charmm->endofstory();
              # Do we need to delete the script?
              if ($FS_DEF->{del_charmm_inp}) { GENERAL::crm("$charmminp"); }
              $pairnum++; push(@comp, shift @pa);
            }
            last if ((scalar(@pa) == 0) || (time() - $T > $dt));
          }
          last if ((scalar(@pa) == 0) || (time() - $T > $dt));
        }
        last if ((scalar(@pa) == 0) || (time() - $T > $dt));
      }
      # If I have reached here, I have either finished my tasks or it is time to re-negotiate tasks splitting
      if (scalar(@pa) == 0) { print "Process $PRANK_DEF finished all its tasks...\n"; }
      print "Process $PRANK_DEF waiting to re-negotiate task splitting...\n";
      MPI_Barrier($MPI_COMM_WORLD_DEF);
      print "Process $PRANK_DEF negotiating...\n";
      my $tcomp = (); # list of completed tasks among all processes
      if ($PRANK_DEF != 0) {
        # Tell root what I have done
        MPI_Send(\@comp, 0, 0, $MPI_COMM_WORLD_DEF);
        # Receive from root a list of all completed tasks
        $tcomp = MPI_Recv(0, 0, $MPI_COMM_WORLD_DEF);
      } else {
        push(@$tcomp, @comp); # what root has done
        for (my $i = 1; $i < $NPROC_DEF; $i++) {
          next if (scalar(@{$a->[$i]}) == 0); # is this process supposed to be doing anything
          my $ar = MPI_Recv($i, 0, $MPI_COMM_WORLD_DEF); # how much proc $i has done
          push(@$tcomp, @$ar);
        }
        for (my $i = 1; $i < $NPROC_DEF; $i++) {
          next if (scalar(@{$a->[$i]}) == 0); # is this process supposed to be doing anything
          MPI_Send($tcomp, $i, 0, $MPI_COMM_WORLD_DEF);
        }
      }
      # If not yet all done, re-negotiate and re-set time
      if (scalar(@$tcomp) != $Nt) {
        $a = GENERAL::splitTasksNonContig($Nt, $NPROC_DEF, $tcomp);
        @pa = @{$a->[$PRANK_DEF]};
        $T = time();
        printf("Process $PRANK_DEF done negotiating and has been assigned %d tasks\n", scalar(@pa));
      } else { last; }
    }
    print "Process $PRANK_DEF done - no more tasks...\n";
    # Wait until all processes are done before putting the data together
    GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF, $NFS_DELAY);
    
    # Combine all tables into one sorted by rotamer
    if ($PRANK_DEF == 0) {
      my $pairnum=0;
      my $tablef = $FS_DEF->{pair_eef_tot_tabf};
      my $tfh = GENERAL::GetOutFH($tablef, "bzip2");
      for (my $s1=0; $s1 <=$#{$sitelist}; $s1++) {
        for (my $ri1=0; $ri1<=$#{$sitelist->[$s1]->{reslist}}; $ri1++) {
          # Collect all files handles for tables involving this residue at this site
          # and all residues with which this residue interacts
          my @fh = ();
          my @f = ();
          my @res2 = ();
          my $res1 = $sitelist->[$s1]->{reslist}->[$ri1];
          for (my $s2=$s1+1; $s2 <=$#{$sitelist}; $s2++) {
            for (my $ri2=0; $ri2<=$#{$sitelist->[$s2]->{reslist}}; $ri2++) {
              push(@res2, $sitelist->[$s2]->{reslist}->[$ri2]);
              $pairnum++;
              my $energytab = "$FS_DEF->{dpair_eef}/$s1\_$ri1/paireef$pairnum.tab";
              $energytab =~ s/\%/$pairnum/;
              push(@f, $energytab);
              push(@fh, GENERAL::GetInFH($energytab));
            }
          }
          # Now merge all tables involving the current residue
          foreach my $rot1 (@{$res1->{rotamerlist}}) {
            next if ($rot1->{valid} == 0);
            for (my $ri2=0; $ri2 < scalar(@res2); $ri2++) {
              my $res2 = $res2[$ri2];
              foreach my $rot2 (@{$res2->{rotamerlist}}) {
                next if ($rot2->{valid} == 0);
                my $fh = $fh[$ri2];
                my $line = <$fh>;
                # Get rid of extra blank spaces - waste of space
                $line =~ s/[ \t]+/ /g; $line =~ s/^[ \t]//; $line =~ s/[ \t]$//;
                print $tfh $line;
              }
            }
          }
          # Close all opened files
          foreach my $fh (@fh) {
            close($fh);
          }
          # Delete smaller table files if needed
          if ($FS_DEF->{del_pair_eef_tabf}) { GENERAL::csystem("rm -r $FS_DEF->{dpair_eef}/$s1\_$ri1/"); }
        }
      }
      close($tfh);
      # Document the created energy table
      push(@{$en->{files}}, File::Spec->abs2rel($tablef, $FS_DEF->{dhome}));
    }
    if ($FS_DEF->{del_charmm_inp}) { GENERAL::crmdir("$wdir"); }
  }

  # Save term tracker
  if ($PRANK_DEF == 0) {
    GENERAL::SaveStruct($TERM_TRACKER, $FS_DEF->{term_trackerf});
  }
#  GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF, $NFS_DELAY);
}


sub writePairwiseEEF_fast {

  my $rotamer = shift;
  my $escheme = shift;
  GENERAL::requireArgs($rotamer, $escheme);
  my $resume = shift;
  if (!defined($resume)) { $resume = 0; }

  # If resuming a previous run, see if the total table has already been created
  if ($resume && GENERAL::fileCheck($FS_DEF->{pair_eef_tot_tabf}, "e") && (GENERAL::GetFileInfo($FS_DEF->{pair_eef_tot_tabf}, 'lines') == $rotamer->countValidPairTotal())) {
    GENERAL::GLog("Skipping pair EEF calculation - has been completed already...");
    return;
  }

  # -------------------------------------------------------
  # Get the term tracker structure and setup the new entry
  my $en;
  if ($PRANK_DEF == 0) {
    $en = ENERGY::addTrackerTerm($TERM_TRACKER, 'pair', 'pair_eef_fast', 'Pairwise energy using CHARMM EEF1',
           ("Desovlation energy from EEF",
           "Van Der Waals energy from EEF",
           "Electrostatic energy from EEF with reduced charge"));
    # Do the following for every energy scheme this term belongs to
    ENERGY::addSchemeToTerm($en, '1m', (0, 2), ('eef_sf', 'eef_sf'));
  }
  # -------------------------------------------------------

  # Do we need to run this term under the current escheme?
  if (($escheme eq '1m')) {
    my $time = time();
    my @terms = ("asp", "vdw", "elec");
    my $ctonnb = 7.0; my $ctofnb = 9.0; my $cutnb = 10.0; # CHARMM distance cutoffs
    my $dt = 10**10;
    $dt = 10*60 if ($NPROC_DEF > 1); # time interval, in seconds, for tasks split re-negotiation
    my $T = time(); # current time in seconds
    $rotamer->numify();
    $rotamer->extent(); # calculate the maxmimum extent of each side chain and get the coordinates of the base atom
    my $sitelist = $rotamer->{dslist};
    my $cdir = GENERAL::GetDir();
    my $wdir = "$LOCAL_DEF/eef_pair_$PRANK_DEF";
    GENERAL::cmkdir($wdir);
    # We want to make the number of pair files in a directory to grow linearly with
    # the number of design residues and not quadratically, so we create another layers of directories.
    if ($PRANK_DEF == 0) {
      for (my $s = 0; $s < scalar(@{$rotamer->{dslist}}); $s++) {
        for (my $res = 0; $res < scalar(@{$rotamer->{dslist}->[$s]->{reslist}}); $res++) {
          GENERAL::cmkdir("$FS_DEF->{dpair_eef}/$s\_$res");
        }
      }
    }
    MPI_Barrier($MPI_COMM_WORLD_DEF);
    my $Nt = $rotamer->countResPairTotal();
    # Divide the task among processes
    my $a = GENERAL::splitTasks($rotamer->countResPairTotal(), $NPROC_DEF);
    my $N = $a->[$PRANK_DEF]->{end} - $a->[$PRANK_DEF]->{beg} + 1;
    if ($N > 0) {
      print "Running pair EEF calculation (process $PRANK_DEF)...\n";
      my $sitelist = $rotamer->{dslist};
      my $chidef = $rotamer->{chidef};
      my $pairnum = 0;

      my $charmminp = "$wdir/paireef$PRANK_DEF.inp";
      my $charmm = CHARMM::new($CHARMM_DEF, $charmminp);
      $charmm->loadPairParm('EEF1');
      # loop over all the pairwise term
      for (my $s1 = 0; $s1 < scalar(@$sitelist); $s1++) {
        my $site1 = $sitelist->[$s1];
        my $chain1 = $sitelist->[$s1]->{chain};
        my $iresnum1 = $sitelist->[$s1]->{iresnum};
        for (my $res1 = 0; $res1 < scalar(@{$sitelist->[$s1]->{reslist}}); $res1++) {
          my $resname1 = $sitelist->[$s1]->{reslist}->[$res1]->{resname};
          my $re1 = $sitelist->[$s1]->{reslist}->[$res1];
          # Read all rotamers of this residue at this site
          my $rotpdb1 = $FS_DEF->{rot_pdbf};
          $rotpdb1 =~ s/%%/$s1/;
          $rotpdb1 =~ s/%/$res1/;
          my $rots1 = PDB::new($rotpdb1, "PDB", "CHARMM");
          my @rots1 = $rots1->ConRes();

          for (my $s2 = $s1+1; $s2 < scalar(@{$sitelist}); $s2++) {
            my $chain2 = $sitelist->[$s2]->{chain};
            my $iresnum2 = $sitelist->[$s2]->{iresnum};
            for (my $res2 = 0; $res2 < scalar(@{$sitelist->[$s2]->{reslist}}); $res2++) {
              $pairnum++;
              next if (($pairnum < $a->[$PRANK_DEF]->{beg}) || ($pairnum > $a->[$PRANK_DEF]->{end}));
              
              my $resname2 = $sitelist->[$s2]->{reslist}->[$res2]->{resname};
              my $re2 = $sitelist->[$s2]->{reslist}->[$res2];
              
              my $energytab = "$FS_DEF->{dpair_eef}/$s1\_$res1/paireef$pairnum.tab";
              $energytab =~ s/\%/$pairnum/;
              $energytab = File::Spec->abs2rel($energytab);

              # if resuming a previous run, skip if this pair of residues has been completely treated
              if ($resume && GENERAL::fileCheck($energytab, "e") && (GENERAL::GetFileInfo($energytab, 'lines') == 
                       $sitelist->[$s1]->{reslist}->[$res1]->{numvalidrot}*$sitelist->[$s2]->{reslist}->[$res2]->{numvalidrot})) {
                print "Skipping pair $pairnum - has already been calculated in a previous run!\n";
                next;
              }

              # See if the two residues are close enough to bother calling CHARMM
              my $d = sqrt(($re1->{base_coor}->[0]-$re2->{base_coor}->[0])**2 +
                           ($re1->{base_coor}->[1]-$re2->{base_coor}->[1])**2 +
                           ($re1->{base_coor}->[2]-$re2->{base_coor}->[2])**2);
              if ($d - $re1->{extent} - $re2->{extent} > $cutnb) {
                #print "Process $PRANK_DEF: pair $pairnum too distant - filling with zeros...\n";
                my $fh = GENERAL::GetOutFH($energytab);
                my $np = ($re1->{numvalidrot})*(($re2->{numvalidrot}));
                for (my $i = 1; $i <= $np; $i++) {
                  $fh->printf(("0 "x(scalar(@terms)) . "\n"));
                }
                close($fh);
                next;
              }

              # Read all rotamers of this residue at this site
              my $rotpdb2 = $FS_DEF->{rot_pdbf};
              $rotpdb2 =~ s/%%/$s2/;
              $rotpdb2 =~ s/%/$res2/;
              my $rots2 = PDB::new($rotpdb2, "PDB", "CHARMM");
              my @rots2 = $rots2->ConRes();

              $charmm->openCard($energytab, 18);
              $charmm->setPairSequence($resname1, 'A');
              for (my $rot1 = 0; $rot1<=$#{$sitelist->[$s1]->{reslist}->[$res1]->{rotamerlist}}; $rot1++) {
                next if ($sitelist->[$s1]->{reslist}->[$res1]->{rotamerlist}->[$rot1]->{valid} == 0);
                my $rotamer1 = $rots1[$rot1];
                $charmm->setPairSequence($resname2, 'B');
                $charmm->setRotCoord("A", 1, $rotamer1);
                my $first = 1;
                for(my $rot2 = 0; $rot2<=$#{$sitelist->[$s2]->{reslist}->[$res2]->{rotamerlist}}; $rot2++) {
                  next if ($sitelist->[$s2]->{reslist}->[$res2]->{rotamerlist}->[$rot2]->{valid} == 0);
                  my $rotamer2 = $rots2[$rot2];
                  $charmm->setRotCoord("B", 1, $rotamer2);
                  if ($escheme eq "1m") { $charmm->setupIMM1(); }
                  else { $charmm->setupEEF1(1, $ctonnb, $ctofnb, $cutnb); }
                  $first = 0;
                  $charmm->inte("sele segid A end", "sele segid B end", 18, join(" ", @terms));
                }
                $charmm->verbose("dele atom select SEGID B end");
              }
              $charmm->verbose("dele atom select SEGID A end");
              $charmm->closeCard(18);
              if (time() - $time > 30*$LOCAL_TIMEOUT) { GENERAL::touchLocalSpace(); $time = time(); }
            }
          }
        }
      }
      $charmm->endofstory();
      # Do we need to delete the script?
      if ($FS_DEF->{del_charmm_inp}) { GENERAL::crm("$charmminp"); }
    }
    # Wait until all processes are done before putting the data together
    GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF, $NFS_DELAY);
    
    # Combine all tables into one sorted by rotamer
    if ($PRANK_DEF == 0) {
      my $pairnum=0;
      my $tablef = $FS_DEF->{pair_eef_tot_tabf};
      my $tfh = GENERAL::GetOutFH($tablef, "bzip2");
      for (my $s1=0; $s1 <=$#{$sitelist}; $s1++) {
        for (my $ri1=0; $ri1<=$#{$sitelist->[$s1]->{reslist}}; $ri1++) {
          # Collect all files handles for tables involving this residue at this site
          # and all residues with which this residue interacts
          my @fh = ();
          my @f = ();
          my @res2 = ();
          my $res1 = $sitelist->[$s1]->{reslist}->[$ri1];
          for (my $s2=$s1+1; $s2 <=$#{$sitelist}; $s2++) {
            for (my $ri2=0; $ri2<=$#{$sitelist->[$s2]->{reslist}}; $ri2++) {
              push(@res2, $sitelist->[$s2]->{reslist}->[$ri2]);
              $pairnum++;
              my $energytab = "$FS_DEF->{dpair_eef}/$s1\_$ri1/paireef$pairnum.tab";
              $energytab =~ s/\%/$pairnum/;
              push(@f, $energytab);
              push(@fh, GENERAL::GetInFH($energytab));
            }
          }
          # Now merge all tables involving the current residue
          foreach my $rot1 (@{$res1->{rotamerlist}}) {
            next if ($rot1->{valid} == 0);
            for (my $ri2=0; $ri2 < scalar(@res2); $ri2++) {
              my $res2 = $res2[$ri2];
              foreach my $rot2 (@{$res2->{rotamerlist}}) {
                next if ($rot2->{valid} == 0);
                my $fh = $fh[$ri2];
                my $line = <$fh>;
                # Get rid of extra blank spaces - waste of space
                $line =~ s/[ \t]+/ /g; $line =~ s/^[ \t]//; $line =~ s/[ \t]$//;
                print $tfh $line;
              }
            }
          }
          # Close all opened files
          foreach my $fh (@fh) {
            close($fh);
          }
          # Delete smaller table files if needed
          if ($FS_DEF->{del_pair_eef_tabf}) { GENERAL::csystem("rm -r $FS_DEF->{dpair_eef}/$s1\_$ri1/"); }
        }
      }
      close($tfh);
      # Document the created energy table
      push(@{$en->{files}}, File::Spec->abs2rel($tablef, $FS_DEF->{dhome}));
    }
    if ($FS_DEF->{del_charmm_inp}) { GENERAL::crmdir("$wdir"); }
  }

  # Save term tracker
  if ($PRANK_DEF == 0) {
    GENERAL::SaveStruct($TERM_TRACKER, $FS_DEF->{term_trackerf});
  }
}


=head2

 Title   :  writePairwiseEEFnoCHARMM
 Usage   :
 Function:
 Returns :
 Args    :

=cut

sub writePairwiseEEFnoCHARMM {

  my $rotamer = shift;
  my $escheme = shift;
  GENERAL::requireArgs($rotamer, $escheme);
  my $resume = shift;
  if (!defined($resume)) { $resume = 0; }
  
  # If resuming a previous run, see if the total table has already been created
  if ($resume && GENERAL::fileCheck($FS_DEF->{pair_eef_tot_tabf}, "e") && (GENERAL::GetFileInfo($FS_DEF->{pair_eef_tot_tabf}, 'lines') == $rotamer->countValidPairTotal())) {
    GENERAL::GLog("Skipping pair EEF calculation - has been completed already...");
    return;
  }

  # -------------------------------------------------------
  # Get the term tracker structure and setup the new entry
  my $en;
  if ($PRANK_DEF == 0) {
    $en = ENERGY::addTrackerTerm($TERM_TRACKER, 'pair', 'pair_eef_nochrm', 'Pairwise energy using standalone EEF1',
           ("Desovlation energy from EEF"));
    # Do the following for every energy scheme this term belongs to
    ENERGY::addSchemeToTerm($en, '1', (0), ('eef_sf'));
    ENERGY::addSchemeToTerm($en, '1h', (0), ('eef_sf'));
    ENERGY::addSchemeToTerm($en, '1s', (0), ('eef_sf'));
    ENERGY::addSchemeToTerm($en, 'c1', (0), ('eef_sf'));
  }
  # -------------------------------------------------------

  # Do we need to run this term under the current escheme?
  if (($escheme eq '1') || ($escheme eq '1h') || ($escheme eq 'c1') || ($escheme eq '1s')) {
    my $time = time();
    $rotamer->numify();
    # Load CHARMM runtime parameters
    my $charmm_par = GENERAL::LoadStruct(DEFINITIONS::getParmFile("CHARMM.param"), 1);
    my $cutnb = $charmm_par->{cutnb}; # CHARMM distance cutoff
    # Load CHARMM parameter and topology files and EEF solvpar file
    my $Peef = CHARMM::readEEFParameters(DEFINITIONS::getParmFile("solvpar.inp"));
    my $T = CHARMM::readCHARMMTopology(DEFINITIONS::getParmFile("toph19_eef1.inp"));
    $T = CHARMM::readCHARMMTopology(DEFINITIONS::getParmFile("patchtop19_eef1.inp"), $T);
    my $P = CHARMM::readCHARMMParameters(DEFINITIONS::getParmFile("param19_eef1.inp"));
    CHARMM::expandCHARMMParameters($P, $T);
    # Calculate the maxmimum extent of each side chain and get the coordinates of the base atom
    $rotamer->extent();
    my $sitelist = $rotamer->{dslist};
    # We want to make the number of pair files in a directory to grow linearly with
    # the number of design residues and not quadratically, so we create another layers of directories.
    if ($PRANK_DEF == 0) {
      for (my $s = 0; $s < scalar(@{$rotamer->{dslist}}); $s++) {
        for (my $res = 0; $res < scalar(@{$rotamer->{dslist}->[$s]->{reslist}}); $res++) {
          GENERAL::cmkdir("$FS_DEF->{dpair_eef}/$s\_$res");
        }
      }
    }
    MPI_Barrier($MPI_COMM_WORLD_DEF);
    my $Nt = $rotamer->countResPairTotal();
    # Divide the tasks among processes
    my $a = GENERAL::splitTasks($rotamer->countResPairTotal(), $NPROC_DEF);
    my $N = $a->[$PRANK_DEF]->{end} - $a->[$PRANK_DEF]->{beg} + 1;
    if ($N > 0) {
      print "Running pair EEF calculation (process $PRANK_DEF)...\n";
      my $sitelist = $rotamer->{dslist};
      my $chidef = $rotamer->{chidef};
      my $pairnum = 0;

      # Loop over pairs of sites and pairs of residues
      for (my $si = 0; $si < scalar(@$sitelist); $si++) {
        my $sitei = $sitelist->[$si];
        for (my $ri = 0; $ri < scalar(@{$sitei->{reslist}}); $ri++) {
          my $resi = $sitei->{reslist}->[$ri];
          # Read all rotamers of this residue at this site
          my $rpdbi = $FS_DEF->{rot_pdbf};
          $rpdbi =~ s/%%/$si/;
          $rpdbi =~ s/%/$ri/;
          my $rotsi = PDB::new($rpdbi, "PDB", "CHARMM", 1);
          my $roti = $rotsi->{chain}->[0]->{res}->[0];
          my @sci = PDB::sidechain($roti);
          my $parami = "";
          foreach my $ai (@sci) {
            my $ati = $T->{residues}->{$roti->{resname}}->{atoms}->{$ai->{atomname}}->{type};
            $parami .= sprintf("%f %f %f %f %f %f %d %f\n", $Peef->{$ati}{V}, $Peef->{$ati}{Gref}, $Peef->{$ati}{Gfree}, $Peef->{$ati}{Href}, $Peef->{$ati}{CPref}, $Peef->{$ati}{SigW}, $T->{residues}->{$resi->{resname}}->{atoms}->{$ai->{atomname}}->{group}, $P->{vdW}{$ati}{rmin});
          }

          for (my $sj = $si+1; $sj < scalar(@{$sitelist}); $sj++) {
            my $sitej = $sitelist->[$sj];
            for (my $rj = 0; $rj < scalar(@{$sitej->{reslist}}); $rj++) {
              $pairnum++;
              next if (($pairnum < $a->[$PRANK_DEF]->{beg}) || ($pairnum > $a->[$PRANK_DEF]->{end}));
              my $resj = $sitej->{reslist}->[$rj];
              
              my $energytab = "$FS_DEF->{dpair_eef}/$si\_$ri/paireef$pairnum.tab";
              $energytab = File::Spec->abs2rel($energytab);

              # if resuming a previous run, skip if this pair of residues has been completely treated
              if ($resume && GENERAL::fileCheck($energytab, "e") && (GENERAL::GetFileInfo($energytab, 'lines') == 
                  $resi->{numvalidrot}*$resj->{numvalidrot})) {
                print "Skipping pair $pairnum - has already been calculated in a previous run!\n";
                next;
              }

              # See if the two residues are close enough to bother calculating
              my $d = sqrt(($resi->{base_coor}->[0]-$resj->{base_coor}->[0])**2 +
                           ($resi->{base_coor}->[1]-$resj->{base_coor}->[1])**2 +
                           ($resi->{base_coor}->[2]-$resj->{base_coor}->[2])**2);
              if ($d - $resi->{extent} - $resj->{extent} > $cutnb) {
                #print "Process $PRANK_DEF: pair $pairnum too distant - filling with zeros...\n";
                my $fh = GENERAL::GetOutFH($energytab);
                my $np = ($resi->{numvalidrot})*(($resj->{numvalidrot}));
                $fh->printf("0\n" x $np);
                close($fh);
                next;
              }

              # Read all rotamers of this residue at this site
              my $rpdbj = $FS_DEF->{rot_pdbf};
              $rpdbj =~ s/%%/$sj/;
              $rpdbj =~ s/%/$rj/;
              my $rotsj = PDB::new($rpdbj, "PDB", "CHARMM", 1);
              my $rotj = $rotsj->{chain}->[0]->{res}->[0];
              my @scj = PDB::sidechain($rotj);
              my $paramj = "";
              foreach my $aj (@scj) {
                my $atj = $T->{residues}->{$rotj->{resname}}->{atoms}->{$aj->{atomname}}->{type};
                $paramj .= sprintf("%f %f %f %f %f %f %d %f\n", $Peef->{$atj}{V}, $Peef->{$atj}{Gref}, $Peef->{$atj}{Gfree}, $Peef->{$atj}{Href}, $Peef->{$atj}{CPref}, $Peef->{$atj}{SigW}, $T->{residues}->{$resj->{resname}}->{atoms}->{$aj->{atomname}}->{group}, $P->{vdW}{$atj}{rmin});
              }
              
              # Call eef1.exe
              my $inpi = $FS_DEF->{rot_outf}; $inpi =~ s/%%/$si/; $inpi =~ s/%/$ri/;
              my $inpj = $FS_DEF->{rot_outf}; $inpj =~ s/%%/$sj/; $inpj =~ s/%/$rj/;
              my $eef = new IO::File;
              $eef->open("| $BIN_DEF/eef1.exe $inpi $inpj 298.15 10 > $energytab");
              $eef->printf("$parami$paramj");
              $eef->close();
              if (time() - $time > 30*$LOCAL_TIMEOUT) { GENERAL::touchLocalSpace(); $time = time(); }
            }
          }
        }
      }
    }
    # Wait until all processes are done before putting the data together
    GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF, $NFS_DELAY);
    
    # Combine all tables into one sorted by rotamer
    if ($PRANK_DEF == 0) {
      my $pairnum=0;
      my $tablef = $FS_DEF->{pair_eef_tot_tabf};
      my $tfh = GENERAL::GetOutFH($tablef, "bzip2");
      for (my $s1=0; $s1 <=$#{$sitelist}; $s1++) {
        for (my $ri1=0; $ri1<=$#{$sitelist->[$s1]->{reslist}}; $ri1++) {
          # Collect all files handles for tables involving this residue at this site
          # and all residues with which this residue interacts
          my @fh = ();
          my @f = ();
          my @res2 = ();
          my $res1 = $sitelist->[$s1]->{reslist}->[$ri1];
          for (my $s2=$s1+1; $s2 <=$#{$sitelist}; $s2++) {
            for (my $ri2=0; $ri2<=$#{$sitelist->[$s2]->{reslist}}; $ri2++) {
              push(@res2, $sitelist->[$s2]->{reslist}->[$ri2]);
              $pairnum++;
              my $energytab = "$FS_DEF->{dpair_eef}/$s1\_$ri1/paireef$pairnum.tab";
              $energytab =~ s/\%/$pairnum/;
              push(@f, $energytab);
              push(@fh, GENERAL::GetInFH($energytab));
            }
          }
          # Now merge all tables involving the current residue
          foreach my $rot1 (@{$res1->{rotamerlist}}) {
            next if ($rot1->{valid} == 0);
            for (my $ri2=0; $ri2 < scalar(@res2); $ri2++) {
              my $res2 = $res2[$ri2];
              foreach my $rot2 (@{$res2->{rotamerlist}}) {
                next if ($rot2->{valid} == 0);
                my $fh = $fh[$ri2];
                my $line = <$fh>;
                print $tfh $line;
              }
            }
          }
          # Close all opened files
          foreach my $fh (@fh) {
            close($fh);
          }
          # Delete smaller table files if needed
          if ($FS_DEF->{del_pair_eef_tabf}) { GENERAL::csystem("rm -r $FS_DEF->{dpair_eef}/$s1\_$ri1/"); }
        }
      }
      close($tfh);
      # Document the created energy table
      push(@{$en->{files}}, File::Spec->abs2rel($tablef, $FS_DEF->{dhome}));
    }
  }

  # Save term tracker
  if ($PRANK_DEF == 0) {
    GENERAL::SaveStruct($TERM_TRACKER, $FS_DEF->{term_trackerf});
  }
}


=head2

 Title   : writePairwiseDDE
 Usage   :
 Function: Writes and runs a series of charmm scripts for generating pairwise energies
           (of the rotamers defiend in the given rotamer object) using multe.
           Returns an array of script names.
 Returns :
 Args    : 1. ROTAMER object corresponding to the design problem.
           2. Energy scheme.
           3. Resume flag - optional.
           4. Custom rdie to use with DDE - optional.

=cut

sub writePairwiseDDE {
  GENERAL::error("This function is depricated - is not guaranteed to handle all parameters or behave correctly!");

  my $rotamer = shift;
  my $escheme = shift;
  GENERAL::requireArgs($rotamer, $escheme);
  my $resume = shift;
  if (!defined($resume)) { $resume = 0; }
  my $rdie = shift;

  # If resuming a previous run, see if the total table has already been created
  if ($resume && GENERAL::fileCheck($FS_DEF->{pair_multe_tot_tabf}, "e") && (GENERAL::GetFileInfo($FS_DEF->{pair_multe_tot_tabf}, 'lines') == $rotamer->countValidPairTotal())) {
    GENERAL::GLog("Skipping pair DDE calculation - has been completed already...");
    return;
  }

  # -------------------------------------------------------
  # Get the term tracker structure and setup the new entry
  my $en;
  if ($PRANK_DEF == 0) {
#    $th = GENERAL::LoadStruct($FS_DEF->{term_trackerf}, 1);
   $en = ENERGY::addTrackerTerm($TERM_TRACKER, 'pair', 'pair_dde', 'CHARMM (rot - rot interaction)',
           ("Van Der energy from MultE calculation",
            "Electrostatic energy from MultE with reduced charge",
            "H-bond energy from MulteE calculation"));
    # Do the following for every energy scheme this term belongs to
    ENERGY::addSchemeToTerm($en, '1', (0, 1), ('vdw_sf', 'dde_pair_sf'));
    ENERGY::addSchemeToTerm($en, '1h', (0, 1, 2), ('vdw_sf', 'dde_pair_sf', 'hbond_sf'));
    ENERGY::addSchemeToTerm($en, 'vdw', (0), ('vdw_sf'));
    ENERGY::addSchemeToTerm($en, '2', (0), ('vdw_sf'));
    ENERGY::addSchemeToTerm($en, '2b', (0), ('vdw_sf'));
    ENERGY::addSchemeToTerm($en, '2c', (0), ('vdw_sf'));
    ENERGY::addSchemeToTerm($en, 'c1', (0, 1, 2), ('vdw_sf', 'dde_pair_sf', 'hbond_sf'));
    ENERGY::addSchemeToTerm($en, '3', (0, 1), ('vdw_sf', 'dde_pair_sf'));
    ENERGY::addSchemeToTerm($en, '4', (1), ('dde_pair_sf'));
    ENERGY::addSchemeToTerm($en, '1e', (0), ('vdw_sf'));
    ENERGY::addSchemeToTerm($en, '1m', (0), ('vdw_sf'));
  }
  # -------------------------------------------------------

  # Do we need to run this term under the current escheme?
  if (($escheme eq '1') || ($escheme eq '1h') || ($escheme eq 'vdw') || ($escheme eq '2') ||
      ($escheme eq '2b') || ($escheme eq '2c') || ($escheme eq 'c1') || ($escheme eq '3') || ($escheme eq '4') ||
      ($escheme eq '1e') || ($escheme eq '1m')) {
    $rotamer->numify();
    my $sitelist = $rotamer->{dslist};
    my $cdir = GENERAL::GetDir();
    my $wdir = "$LOCAL_DEF/dde_pair_$PRANK_DEF";
    GENERAL::cmkdir($wdir);
    # We want to make the number of pair files in a directory to grow linearly with
    # the number of design residues and not quadratically, so we create another layers of directories.
    if ($PRANK_DEF == 0) {
      for (my $s = 0; $s < scalar(@{$rotamer->{dslist}}); $s++) {
        for (my $res = 0; $res < scalar(@{$rotamer->{dslist}->[$s]->{reslist}}); $res++) {
          GENERAL::cmkdir("$FS_DEF->{dpair_multe}/$s\_$res");
        }
      }
    }
    MPI_Barrier($MPI_COMM_WORLD_DEF);
    # Divide the task among processes
    my $a = GENERAL::splitTasks($rotamer->countResPairTotal(), $NPROC_DEF);
    my $N = $a->[$PRANK_DEF]->{end} - $a->[$PRANK_DEF]->{beg} + 1;
    if ($N > 0) {
     my @terms = ("vdw", "elec", "hbon");
      print "Running pair DDE, VdW calculation (process $PRANK_DEF)...\n";

      my $sitelist = $rotamer->{dslist};
      my $chidef = $rotamer->{chidef};
      my $pairnum = 0;
      # loop over all the pairwise terms
      for (my $s1 = 0; $s1 < scalar(@{$rotamer->{dslist}}); $s1++) {
        my $site1 = $rotamer->{dslist}->[$s1];
        my $chain1 = $site1->{chain};
        my $iresnum1 = $site1->{iresnum};
        for (my $r1 = 0; $r1 < scalar(@{$site1->{reslist}}); $r1++) {
          my $res1 = $site1->{reslist}->[$r1];
          my $resname1 = uc($res1->{resname});

          for (my $s2 = $s1+1; $s2 < scalar(@{$rotamer->{dslist}}); $s2++) {
            my $site2 = $rotamer->{dslist}->[$s2];
            my $chain2 = $site2->{chain};
            my $iresnum2 = $site2->{iresnum};
            for (my $r2 = 0; $r2 < scalar(@{$site2->{reslist}}); $r2++) {
              my $res2 = $site2->{reslist}->[$r2];
              my $resname2 = uc($res2->{resname});

              # every residue pair has a pair.inp
              $pairnum++;
              next if (($pairnum < $a->[$PRANK_DEF]->{beg}) || ($pairnum > $a->[$PRANK_DEF]->{end}));

              # --- Make sure multe does not run on empty files -----------------------------------------------
              if (($res1->{numvalidrot})*($res2->{numvalidrot}) == 0) { next; }
              # -----------------------------------------------------------------------------------------------

              my $energytab = "$FS_DEF->{dpair_multe}/$s1\_$r1/pairmulte$pairnum.tab";
              $energytab =~ s/\%/$pairnum/;
              $energytab = File::Spec->abs2rel($energytab);
              # if resuming a previous run, skip if this pair of residues has been completely treated
              if ($resume && GENERAL::fileCheck($energytab, "e") && (GENERAL::GetFileInfo($energytab, 'lines') == $res1->{numvalidrot}*$res2->{numvalidrot})) {
                print "Skipping pair $pairnum - has already been calculated in a previous run!\n";
                next;
              }
              # Don't run multe for GLY because it freaks out for empty sidechains
              if (($resname1 =~ /GLY/) || ($resname2 =~ /GLY/)) {
                my $fh = GENERAL::GetOutFH($energytab);
                my $np = ($res1->{numvalidrot})*(($res2->{numvalidrot}));
                for (my $i = 1; $i <= $np; $i++) {
                  $fh->printf(("0 "x(scalar(@terms)) . "\n"));
                }
                close($fh);
                next;
              }
              my $charmminp = "$wdir/pairmulte$pairnum.inp";
              my $charmm = CHARMM::new($CHARMM_DEF, $charmminp);
              $charmm->loadPairParm("scaled hbond");
              $charmm->setPairSequence($resname1, 'A');
              $charmm->setPairSequence($resname2, 'B');
              $charmm->verbose("ic param all\nhbuild\n\n");
              $charmm->setupMultePairNonBonded($resname1, $resname2, 1, $rdie);
              $charmm->verbose("skipe all exclude " . join(" exclude ", @terms));
              $charmm->selectINTE('sele segid A end', 'sele segid B end', 'keep');
              $charmm->setupMultePair($s1, $r1, $s2, $r2, $energytab, @terms);
              $charmm->endofstory();

              # Do we need to delete the script
              if ($FS_DEF->{del_charmm_inp}) { GENERAL::crm("$charmminp"); }
            }
          }
        }
      }
    }
    if ($FS_DEF->{del_charmm_inp}) { GENERAL::crmdir("$wdir"); }
    # Wait until all processes are done before putting the data together
    GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF, $NFS_DELAY);

    # Combine all tables into one sorted by rotamer
    if ($PRANK_DEF == 0) {
      my $pairnum=0;
      my $tablef = $FS_DEF->{pair_multe_tot_tabf};
      my $tfh = GENERAL::GetOutFH($tablef, "bzip2");
      for (my $s1 = 0; $s1 <= $#{$sitelist}; $s1++) {
        for (my $ri1 = 0; $ri1 <= $#{$sitelist->[$s1]->{reslist}}; $ri1++) {
          # Collect all files handles for tables involving this residue at this site
          # and all residues with which this residue interacts
          my @fh = ();
          my @f = ();
          my @res2 = ();
          my $res1 = $sitelist->[$s1]->{reslist}->[$ri1];
          for (my $s2 = $s1+1; $s2 <= $#{$sitelist}; $s2++) {
            for (my $ri2 = 0; $ri2 <= $#{$sitelist->[$s2]->{reslist}}; $ri2++) {
              push(@res2, $sitelist->[$s2]->{reslist}->[$ri2]);
              $pairnum++;
              my $energytab = "$FS_DEF->{dpair_multe}/$s1\_$ri1/pairmulte$pairnum.tab";
              $energytab =~ s/\%/$pairnum/;
              push(@f, $energytab);
              push(@fh, GENERAL::GetInFH($energytab));
            }
          }
          # Now merge all tables involving the current residue
          foreach my $rot1 (@{$res1->{rotamerlist}}) {
            next if ($rot1->{valid} == 0);
            for (my $ri2 = 0; $ri2 < scalar(@res2); $ri2++) {
              my $res2 = $res2[$ri2];
              foreach my $rot2 (@{$res2->{rotamerlist}}) {
                next if ($rot2->{valid} == 0);
                my $fh = $fh[$ri2];
                my $line = <$fh>;
                # Get rid of extra blank spaces - waste of space
                $line =~ s/[ \t]+/ /g; $line =~ s/^[ \t]//; $line =~ s/[ \t]$//;
                print $tfh $line;
              }
            }
          }
          # Close all opened files
          foreach my $fh (@fh) {
            close($fh);
          }
          # Delete smaller table files if needed
          if ($FS_DEF->{del_pair_multe_tabf}) { GENERAL::csystem("rm -r $FS_DEF->{dpair_multe}/$s1\_$ri1/"); }
        }
      }
      close($tfh);
      # Document the created energy table
      push(@{$en->{files}}, File::Spec->abs2rel($tablef, $FS_DEF->{dhome}));
    }
  }

  # Save term tracker
  if ($PRANK_DEF == 0) {
    GENERAL::SaveStruct($TERM_TRACKER, $FS_DEF->{term_trackerf});
  }
#  GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF, $NFS_DELAY);
}



=head2

 Title   : writePairwiseDDE_fast
 Usage   :
 Function: Writes and runs a series of charmm scripts for generating pairwise energies
           (of the rotamers defiend in the given rotamer object) using multe.
           Returns an array of script names.
 Returns :
 Args    : 1. ROTAMER object corresponding to the design problem.
           2. Energy scheme.
           3. Resume flag - optional.
           4. Custom rdie to use with DDE - optional.

=cut

sub writePairwiseDDE_fast {

  my $rotamer = shift;
  my $escheme = shift;
  GENERAL::requireArgs($rotamer, $escheme);
  my $resume = shift;
  if (!defined($resume)) { $resume = 0; }
  my $rdie = shift;

  # If resuming a previous run, see if the total table has already been created
  if ($resume && GENERAL::fileCheck($FS_DEF->{pair_multe_tot_tabf}, "e") && (GENERAL::GetFileInfo($FS_DEF->{pair_multe_tot_tabf}, 'lines') == $rotamer->countValidPairTotal())) {
    GENERAL::GLog("Skipping pair DDE calculation - has been completed already...");
    return;
  }

  # -------------------------------------------------------
  # Get the term tracker structure and setup the new entry
  my $en;
  if ($PRANK_DEF == 0) {
   $en = ENERGY::addTrackerTerm($TERM_TRACKER, 'pair', 'pair_dde', 'CHARMM (rot - rot interaction)',
           ("Van Der energy from MultE calculation",
           "Electrostatic energy from MultE with reduced charge",
           "H-bond energy from MulteE calculation"));
    # Do the following for every energy scheme this term belongs to
    ENERGY::addSchemeToTerm($en, '1', (0, 1), ('vdw_sf', 'dde_pair_sf'));
    ENERGY::addSchemeToTerm($en, '1h', (0, 1, 2), ('vdw_sf', 'dde_pair_sf', 'hbond_sf'));
    ENERGY::addSchemeToTerm($en, 'vdw', (0), ('vdw_sf'));
    ENERGY::addSchemeToTerm($en, '2', (0), ('vdw_sf'));
    ENERGY::addSchemeToTerm($en, '2b', (0), ('vdw_sf'));
    ENERGY::addSchemeToTerm($en, '2c', (0), ('vdw_sf'));
    ENERGY::addSchemeToTerm($en, 'c1', (0, 1, 2), ('vdw_sf', 'dde_pair_sf', 'hbond_sf'));
    ENERGY::addSchemeToTerm($en, '3', (0, 1), ('vdw_sf', 'dde_pair_sf'));
    ENERGY::addSchemeToTerm($en, '4', (1), ('dde_pair_sf'));
    ENERGY::addSchemeToTerm($en, '1e', (0), ('vdw_sf'));
    ENERGY::addSchemeToTerm($en, '1m', (0), ('vdw_sf'));
  }
  # -------------------------------------------------------

  # Do we need to run this term under the current escheme?
  if (($escheme eq '1') || ($escheme eq '1h') || ($escheme eq 'vdw') || ($escheme eq '2') ||
      ($escheme eq '2b') || ($escheme eq '2c') || ($escheme eq 'c1') || ($escheme eq '3') || ($escheme eq '4') ||
      ($escheme eq '1e') || ($escheme eq '1m')) {
    my $time = time();
    $rotamer->numify();
    my $sitelist = $rotamer->{dslist};
    # We want to make the number of pair files in a directory to grow linearly with
    # the number of design residues and not quadratically, so we create another layers of directories.
    if ($PRANK_DEF == 0) {
      for (my $s = 0; $s < scalar(@{$rotamer->{dslist}}); $s++) {
        for (my $res = 0; $res < scalar(@{$rotamer->{dslist}->[$s]->{reslist}}); $res++) {
          GENERAL::cmkdir("$FS_DEF->{dpair_multe}/$s\_$res");
        }
      }
    }
    MPI_Barrier($MPI_COMM_WORLD_DEF);
    # Divide the task among processes
    my $a = GENERAL::splitTasks($rotamer->countResPairTotal(), $NPROC_DEF);
    my $N = $a->[$PRANK_DEF]->{end} - $a->[$PRANK_DEF]->{beg} + 1;
#    my $charmminp = "pairmulte$PRANK_DEF.inp";
    my $charmm = CHARMM::new($MULTE_DEF);
    my @terms = ("vdw", "elec");
    if ($escheme eq '1h') {
      $charmm->loadPairParm("scaled hbond");
      push(@terms, "hbond");
    } else { $charmm->loadPairParm("scaled"); }
    if ($N > 0) {
      print "Running pair DDE, VdW calculation (process $PRANK_DEF)...\n";

      my $sitelist = $rotamer->{dslist};
      my $chidef = $rotamer->{chidef};
      my $pairnum = 0;
      # loop over all the pairwise terms
      for (my $s1 = 0; $s1 < scalar(@{$rotamer->{dslist}}); $s1++) {
        my $site1 = $rotamer->{dslist}->[$s1];
        my $chain1 = $site1->{chain};
        my $iresnum1 = $site1->{iresnum};
        for (my $r1 = 0; $r1 < scalar(@{$site1->{reslist}}); $r1++) {
          my $res1 = $site1->{reslist}->[$r1];
          my $resname1 = uc($res1->{resname});

          for (my $s2 = $s1+1; $s2 < scalar(@{$rotamer->{dslist}}); $s2++) {
            my $site2 = $rotamer->{dslist}->[$s2];
            my $chain2 = $site2->{chain};
            my $iresnum2 = $site2->{iresnum};
            for (my $r2 = 0; $r2 < scalar(@{$site2->{reslist}}); $r2++) {
              my $res2 = $site2->{reslist}->[$r2];
              my $resname2 = uc($res2->{resname});

              # every residue pair has a pair.inp
              $pairnum++;
              next if (($pairnum < $a->[$PRANK_DEF]->{beg}) || ($pairnum > $a->[$PRANK_DEF]->{end}));

              # --- Make sure multe does not run on empty files -----------------------------------------------
              if (($res1->{numvalidrot})*($res2->{numvalidrot}) == 0) { next; }
              # -----------------------------------------------------------------------------------------------

              my $energytab = "$FS_DEF->{dpair_multe}/$s1\_$r1/pairmulte$pairnum.tab";
              $energytab =~ s/\%/$pairnum/;
              $energytab = File::Spec->abs2rel($energytab);
              # if resuming a previous run, skip if this pair of residues has been completely treated
              if ($resume && GENERAL::fileCheck($energytab, "e") && (GENERAL::GetFileInfo($energytab, 'lines') == $res1->{numvalidrot}*$res2->{numvalidrot})) {
                print "Skipping pair $pairnum - has already been calculated in a previous run!\n";
                next;
              }
              # Don't run multe for GLY because it freaks out for empty sidechains
              if (($resname1 =~ /GLY/) || ($resname2 =~ /GLY/)) {
                my $fh = GENERAL::GetOutFH($energytab);
                my $np = ($res1->{numvalidrot})*(($res2->{numvalidrot}));
                for (my $i = 1; $i <= $np; $i++) {
                  $fh->printf(("0 "x(scalar(@terms)) . "\n"));
                }
                close($fh);
                next;
              }
              $charmm->setPairSequence($resname1, 'A');
              $charmm->setPairSequence($resname2, 'B');
              $charmm->verbose("ic param all\nhbuild\n\n");
              $charmm->setupMultePairNonBonded($resname1, $resname2, 1, $rdie);
              $charmm->verbose("skipe all exclude " . join(" exclude ", @terms));
              $charmm->selectINTE('sele segid A end', 'sele segid B end', 'keep');
              $charmm->setupMultePair($s1, $r1, $s2, $r2, $energytab, @terms);
              $charmm->verbose("delete atom select all end");
              if (time() - $time > 30*$LOCAL_TIMEOUT) { GENERAL::touchLocalSpace(); $time = time(); }
            }
          }
        }
      }
      $charmm->endofstory();
      # Do we need to delete the script
#      if ($FS_DEF->{del_charmm_inp}) { GENERAL::crm("$charmminp"); }
    }
#    if ($FS_DEF->{del_charmm_inp}) { GENERAL::crmdir("$wdir"); }
    # Wait until all processes are done before putting the data together
    GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF, $NFS_DELAY);

    # Combine all tables into one sorted by rotamer
    if ($PRANK_DEF == 0) {
      my $pairnum=0;
      my $tablef = $FS_DEF->{pair_multe_tot_tabf};
      my $tfh = GENERAL::GetOutFH($tablef, "bzip2");
      for (my $s1 = 0; $s1 <= $#{$sitelist}; $s1++) {
        for (my $ri1 = 0; $ri1 <= $#{$sitelist->[$s1]->{reslist}}; $ri1++) {
          # Collect all files handles for tables involving this residue at this site
          # and all residues with which this residue interacts
          my @fh = ();
          my @f = ();
          my @res2 = ();
          my $res1 = $sitelist->[$s1]->{reslist}->[$ri1];
          for (my $s2 = $s1+1; $s2 <= $#{$sitelist}; $s2++) {
            for (my $ri2 = 0; $ri2 <= $#{$sitelist->[$s2]->{reslist}}; $ri2++) {
              push(@res2, $sitelist->[$s2]->{reslist}->[$ri2]);
              $pairnum++;
              my $energytab = "$FS_DEF->{dpair_multe}/$s1\_$ri1/pairmulte$pairnum.tab";
              $energytab =~ s/\%/$pairnum/;
              push(@f, $energytab);
              push(@fh, GENERAL::GetInFH($energytab));
            }
          }
          # Now merge all tables involving the current residue
          foreach my $rot1 (@{$res1->{rotamerlist}}) {
            next if ($rot1->{valid} == 0);
            for (my $ri2 = 0; $ri2 < scalar(@res2); $ri2++) {
              my $res2 = $res2[$ri2];
              foreach my $rot2 (@{$res2->{rotamerlist}}) {
                next if ($rot2->{valid} == 0);
                my $fh = $fh[$ri2];
                my $line = <$fh>;
                if (!defined($line)) { GENERAL::error("File $f[$ri2] seems to be too short!"); }
                # Get rid of extra blank spaces - waste of space
                $line =~ s/[ \t]+/ /g; $line =~ s/^[ \t]//; $line =~ s/[ \t]$//;
                print $tfh $line;
              }
            }
          }
          # Close all opened files
          for (my $ri2 = 0; $ri2 < scalar(@res2); $ri2++) {
            my $fh = $fh[$ri2];
            if (defined(<$fh>)) { GENERAL::error("File $f[$ri2] seems to be too long!"); }
            close($fh);
          }
#          foreach my $fh (@fh) {
#            close($fh);
#          }
          # Delete smaller table files if needed
          if ($FS_DEF->{del_pair_multe_tabf}) { GENERAL::csystem("rm -r $FS_DEF->{dpair_multe}/$s1\_$ri1/"); }
        }
      }
      close($tfh);
      # Document the created energy table
      push(@{$en->{files}}, File::Spec->abs2rel($tablef, $FS_DEF->{dhome}));
    }
  }

  # Save term tracker
  if ($PRANK_DEF == 0) {
    GENERAL::SaveStruct($TERM_TRACKER, $FS_DEF->{term_trackerf});
  }
}



=head2

 Title   : writePairwiseDDE_min
 Usage   :
 Function: For testing purposes, this function calculates rotamer pairwise interactions (vdW and DDE) after a very short
           minimization step in CHARMM.
 Returns :
 Args    :

=cut

sub writePairwiseDDE_min {

  my $rotamer = shift;
  my $crdf = shift;
  my $escheme = shift;
  GENERAL::requireArgs($rotamer, $escheme);
  my $resume = shift;
  if (!defined($resume)) { $resume = 0; }
  my $rdie = shift;
  
  # If resuming a previous run, see if the total table has already been created
  if ($resume && GENERAL::fileCheck($FS_DEF->{pair_eef_tot_tabf}, "e") && (GENERAL::GetFileInfo($FS_DEF->{pair_eef_tot_tabf}, 'lines') == $rotamer->countValidPairTotal())) {
    GENERAL::GLog("Skipping pair vdW, DDE and EEF calculation - has been completed already...");
    return;
  }

  # -------------------------------------------------------
  # Get the term tracker structure and setup the new entry
  my $en;
  if ($PRANK_DEF == 0) {
    $en = ENERGY::addTrackerTerm($TERM_TRACKER, 'pair', 'pair_dde_min', 'Pairwise CHARMM energy after short minimization',
           ("Van Der Waals energy",
            "DDE Electrostatic energy"));
    # Do the following for every energy scheme this term belongs to
    ENERGY::addSchemeToTerm($en, '1s', (0, 1), ('vdw_sf', 'dde_pair_sf'));
    ENERGY::addSchemeToTerm($en, 'vdws', (0), ('vdw_sf'));
  }
  # -------------------------------------------------------

  # Do we need to run this term under the current escheme?
  if (($escheme eq '1s') || ($escheme eq 'vdws')) {
    my $time = time();
    my @terms = ("vdw", "elec");
    $rotamer->numify();
    $rotamer->extent(); # calculate the maxmimum extent of each side chain and get the coordinates of the base atom
    my $sitelist = $rotamer->{dslist};
    my $cdir = GENERAL::GetDir();
    my $wdir = "$LOCAL_DEF/dde_min_pair_$PRANK_DEF";
    GENERAL::cmkdir($wdir);
    # We want to make the number of pair files in a directory to grow linearly with
    # the number of design residues and not quadratically, so we create another layers of directories.
    if ($PRANK_DEF == 0) {
      for (my $s = 0; $s < scalar(@{$rotamer->{dslist}}); $s++) {
        for (my $res = 0; $res < scalar(@{$rotamer->{dslist}->[$s]->{reslist}}); $res++) {
          GENERAL::cmkdir("$FS_DEF->{dpair_multe}/$s\_$res");
        }
      }
    }
    MPI_Barrier($MPI_COMM_WORLD_DEF);
    my $Nt = $rotamer->countResPairTotal();
    # Divide the task among processes
    my $a = GENERAL::splitTasks($rotamer->countResPairTotal(), $NPROC_DEF);
    my $N = $a->[$PRANK_DEF]->{end} - $a->[$PRANK_DEF]->{beg} + 1;
    if ($N > 0) {
      print "Running pair DDE_min calculation (process $PRANK_DEF)...\n";
      my $sitelist = $rotamer->{dslist};
      my $chidef = $rotamer->{chidef};
      my $pairnum = 0;

#      my $charmminp = "$wdir/pairdde_min$PRANK_DEF.inp";
      my $charmminp = undef;
      my $charmm = CHARMM::new($CHARMM_DEF, $charmminp);
      if (!defined($rdie)) { $rdie = $charmm->{par}->{rdie}; }
      $charmm->loadParm("scaled");
      $charmm->send("skipe all excl bond excl angle excl dihe excl impr " . join(" excld ", @terms) . "\n");
      $charmm->setupFromCRD($crdf);
      $charmm->defineTemplate($rotamer);
      $charmm->send("delete atom select rots end\n");
      # loop over all the pairwise term
      for (my $s1 = 0; $s1 < scalar(@$sitelist); $s1++) {
        my $site1 = $sitelist->[$s1];
        my $chain1 = $sitelist->[$s1]->{chain};
        my $iresnum1 = $sitelist->[$s1]->{iresnum};
        for (my $res1 = 0; $res1 < scalar(@{$sitelist->[$s1]->{reslist}}); $res1++) {
          my $resname1 = $sitelist->[$s1]->{reslist}->[$res1]->{resname};
          my $re1 = $sitelist->[$s1]->{reslist}->[$res1];
          $charmm->patchRes($resname1, $chain1, $iresnum1);

          for (my $s2 = $s1+1; $s2 < scalar(@{$sitelist}); $s2++) {
            my $chain2 = $sitelist->[$s2]->{chain};
            my $iresnum2 = $sitelist->[$s2]->{iresnum};
            for (my $res2 = 0; $res2 < scalar(@{$sitelist->[$s2]->{reslist}}); $res2++) {
              $pairnum++;
              next if (($pairnum < $a->[$PRANK_DEF]->{beg}) || ($pairnum > $a->[$PRANK_DEF]->{end}));
              
              my $resname2 = $sitelist->[$s2]->{reslist}->[$res2]->{resname};
              my $re2 = $sitelist->[$s2]->{reslist}->[$res2];
              $charmm->patchRes($resname2, $chain2, $iresnum2);
              
              my $energytab = "$FS_DEF->{dpair_multe}/$s1\_$res1/pairdde_min$pairnum.tab";
              $energytab =~ s/\%/$pairnum/;
              $energytab = File::Spec->abs2rel($energytab);

              # if resuming a previous run, skip if this pair of residues has been completely treated
              if ($resume && GENERAL::fileCheck($energytab, "e") && (GENERAL::GetFileInfo($energytab, 'lines') == 
                       $sitelist->[$s1]->{reslist}->[$res1]->{numvalidrot}*$sitelist->[$s2]->{reslist}->[$res2]->{numvalidrot})) {
                print "Skipping pair $pairnum - has already been calculated in a previous run!\n";
                next;
              }

              # See if the two residues are close enough to bother calling CHARMM
              my $d = sqrt(($re1->{base_coor}->[0]-$re2->{base_coor}->[0])**2 +
                           ($re1->{base_coor}->[1]-$re2->{base_coor}->[1])**2 +
                           ($re1->{base_coor}->[2]-$re2->{base_coor}->[2])**2);
              if ($d - $re1->{extent} - $re2->{extent} > $charmm->{par}->{cutnb}) {
                #print "Process $PRANK_DEF: pair $pairnum too distant - filling with zeros...\n";
                my $fh = GENERAL::GetOutFH($energytab);
                my $np = ($re1->{numvalidrot})*(($re2->{numvalidrot}));
                for (my $i = 1; $i <= $np; $i++) {
                  $fh->printf(("0 "x(scalar(@terms)) . "\n"));
                }
                close($fh);
                next;
              }

              $charmm->defineTemplate($rotamer);
              $charmm->send("update atom rdie switch eps $rdie e14fac $charmm->{par}->{e14fac} -\n".
                            "cutnb $charmm->{par}->{cutnb} ctofnb  $charmm->{par}->{ctofnb} ctonnb $charmm->{par}->{ctonnb} -\n".
                            "nbxmod $charmm->{par}->{nbxmod} vswitch wmin $charmm->{par}->{wmin} vatom vdistance\n");
              $charmm->openCard($energytab, 18);
              for (my $rot1 = 0; $rot1<=$#{$sitelist->[$s1]->{reslist}->[$res1]->{rotamerlist}}; $rot1++) {
                my $rotamer1 = $sitelist->[$s1]->{reslist}->[$res1]->{rotamerlist}->[$rot1];
                next if ($rotamer1->{valid} == 0);
                $charmm->placeRotamer($chain1, $resname1, $iresnum1, $rotamer->{chidef}, $rotamer1);
                for(my $rot2 = 0; $rot2<=$#{$sitelist->[$s2]->{reslist}->[$res2]->{rotamerlist}}; $rot2++) {
                  my $rotamer2 = $sitelist->[$s2]->{reslist}->[$res2]->{rotamerlist}->[$rot2];
                  next if ($rotamer2->{valid} == 0);
                  $charmm->placeRotamer($chain2, $resname2, $iresnum2, $rotamer->{chidef}, $rotamer2);
                  $charmm->send("cons fix sele template end\n");
                  $charmm->send("mini sd nstep 5 step 0.1 ihbfrq 1 inbfrq 1\n");
                  $charmm->inte("sele .byres. atom $chain1 $iresnum1 * .and. .not. backb end",
                                "sele .byres. atom $chain2 $iresnum2 * .and. .not. backb end", 18, join(" ", @terms));
                }
              }
              $charmm->verbose("dele atom select .byres. atom $chain2 $iresnum2 * .and. .not. (" . CHARMM::backbone() . ")");
              $charmm->closeCard(18);
              # from time to time update the local directory so that it is not lost
              if (time() - $time > 30*$LOCAL_TIMEOUT) { GENERAL::touchLocalSpace(); $time = time(); }
            }
          }
          $charmm->verbose("dele atom select .byres. atom $chain1 $iresnum1 * .and. .not. (" . CHARMM::backbone() . ")");
        }
      }
      $charmm->endofstory();
      # Do we need to delete the script?
      if ($FS_DEF->{del_charmm_inp}) { GENERAL::crm("$charmminp") if (defined($charmminp)); }
    }
    # Wait until all processes are done before putting the data together
    GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF, $NFS_DELAY);
    
    # Combine all tables into one sorted by rotamer
    if ($PRANK_DEF == 0) {
      my $pairnum=0;
      my $tablef = "$FS_DEF->{dpair_multe}/pairdde_min.tab";
      my $tfh = GENERAL::GetOutFH($tablef, "bzip2");
      for (my $s1=0; $s1 <=$#{$sitelist}; $s1++) {
        for (my $ri1=0; $ri1<=$#{$sitelist->[$s1]->{reslist}}; $ri1++) {
          # Collect all files handles for tables involving this residue at this site
          # and all residues with which this residue interacts
          my @fh = ();
          my @f = ();
          my @res2 = ();
          my $res1 = $sitelist->[$s1]->{reslist}->[$ri1];
          for (my $s2=$s1+1; $s2 <=$#{$sitelist}; $s2++) {
            for (my $ri2=0; $ri2<=$#{$sitelist->[$s2]->{reslist}}; $ri2++) {
              push(@res2, $sitelist->[$s2]->{reslist}->[$ri2]);
              $pairnum++;
              my $energytab = "$FS_DEF->{dpair_multe}/$s1\_$ri1/pairdde_min$pairnum.tab";
              $energytab =~ s/\%/$pairnum/;
              push(@f, $energytab);
              push(@fh, GENERAL::GetInFH($energytab));
            }
          }
          # Now merge all tables involving the current residue
          foreach my $rot1 (@{$res1->{rotamerlist}}) {
            next if ($rot1->{valid} == 0);
            for (my $ri2=0; $ri2 < scalar(@res2); $ri2++) {
              my $res2 = $res2[$ri2];
              foreach my $rot2 (@{$res2->{rotamerlist}}) {
                next if ($rot2->{valid} == 0);
                my $fh = $fh[$ri2];
                my $line = <$fh>;
                # Get rid of extra blank spaces - waste of space
                $line =~ s/[ \t]+/ /g; $line =~ s/^[ \t]//; $line =~ s/[ \t]$//;
                print $tfh $line;
              }
            }
          }
          # Close all opened files
          foreach my $fh (@fh) {
            close($fh);
          }
          # Delete smaller table files if needed
          if ($FS_DEF->{del_pair_multe_tabf}) { GENERAL::csystem("rm -r $FS_DEF->{dpair_multe}/$s1\_$ri1/"); }
        }
      }
      close($tfh);
      # Document the created energy table
      push(@{$en->{files}}, File::Spec->abs2rel($tablef, $FS_DEF->{dhome}));
    }
    if ($FS_DEF->{del_charmm_inp}) { GENERAL::crmdir("$wdir"); }
  }

  # Save term tracker
  if ($PRANK_DEF == 0) {
    GENERAL::SaveStruct($TERM_TRACKER, $FS_DEF->{term_trackerf});
  }
}


=head2

 Title   : writePairwiseEPIC
 Usage   : ENERGY::writePairwiseEPIC($rotamer, $escheme)
 Function: Calculates rotamer-rotamer pair interactions by interfacing with EPIC.
           Calculates all terms necessary for the given energy scheme.
 Returns : Nothing
 Args    : 1. rotamer structure
           2. energy scheme
           3. calculation parameter hash table reference
           4. Optional: resume flag

=cut

sub writePairwiseEPIC {
  my $rotamer = shift;
  my $escheme = shift;
  my $parm = shift;
  GENERAL::requireArgs($rotamer, $escheme);
  my $resume = shift;
  if (!defined($resume)) { $resume = 0; }
  my $impdbs = shift; $impdbs = () if (!defined($impdbs));

  # -------------------------------------------------------
  # Get the term tracker structure and setup the new entry
  my %en;
  if ($PRANK_DEF == 0) {
    # all terms EPIC can calculate
    $en{VDW} = ENERGY::addTrackerTerm($TERM_TRACKER, 'pair', 'pair_vdw_epic', 'Pair vdW energies using EPIC', "van der Waals interaction energy");
    ENERGY::addSchemeToTerm($en{VDW}, '1p', (0), ('vdw_sf'));
    ENERGY::addSchemeToTerm($en{VDW}, 'gb-eef', (0), ('vdw_sf')); # test energy function
    ENERGY::addSchemeToTerm($en{VDW}, 'vdwp', (0), ('vdw_sf'));
    $en{EEF1} = ENERGY::addTrackerTerm($TERM_TRACKER, 'pair', 'pair_eef1_epic', 'Pair EEF1 energies using EPIC', "EEF1 pair desolvation energy");
    ENERGY::addSchemeToTerm($en{EEF1}, '1p', (0), ('eef_sf'));
    ENERGY::addSchemeToTerm($en{EEF1}, 'gb-eef', (0), ('eef_sf')); # test energy function
    $en{DDE} = ENERGY::addTrackerTerm($TERM_TRACKER, 'pair', 'pair_dde_epic', 'Pair DDE energies using EPIC', "DDE pair interaction energy");
    ENERGY::addSchemeToTerm($en{DDE}, '1p', (0), ('dde_pair_sf'));
  }
  # -------------------------------------------------------

  # Do we need to run this term under the current escheme?
  my $time = time(); my @terms;
  if (($escheme eq '1p') || ($escheme eq 'gb-eef') || ($escheme eq 'vdwp')) {
    my $tmprdie = $parm->{rdie};
    if ($escheme eq '1p') { @terms = ('VDW', 'EEF1', 'DDE'); }
    if ($escheme eq 'gb-eef') { @terms = ('VDW', 'EEF1'); }
    if ($escheme eq 'vdwp') { @terms = ('VDW'); }
    # We want to make the number of pair files in a directory to grow linearly with
    # the number of design residues and not quadratically, so we create another layer of directories.
    if ($PRANK_DEF == 0) {
      for (my $s = 0; $s < scalar(@{$rotamer->{dslist}}) - 1; $s++) {
        for (my $res = 0; $res < scalar(@{$rotamer->{dslist}->[$s]->{reslist}}); $res++) {
          GENERAL::cmkdir("$FS_DEF->{dpair}/$s\_$res");
        }
      }
    }
    GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF, $NFS_DELAY);
    # Calculate the maxmimum extent of each side chain and get the coordinates of the base atom
    my $epic = undef;
    $rotamer->numify();
    $rotamer->extent();
    my $sitelist = $rotamer->{dslist};
    foreach my $term (@terms) {
      print "Process $PRANK_DEF: Calculating pairwise '$term' using EPIC...\n";
      # Have we done this before?
      my $tablef;
      if ($term eq "VDW") { $tablef = "$FS_DEF->{dpair}/vdw_epic.tab"; }
      elsif ($term eq "DDE") { $tablef = "$FS_DEF->{dpair}/dde_epic.tab"; }
      elsif ($term eq "EEF1") { $tablef = "$FS_DEF->{dpair}/eef1_epic.tab"; }
      else { GENERAL::error("Do not know what to call final table for term '$term'"); }
      # If resuming a previous run, see if the total table has already been created
      if ($resume && GENERAL::fileCheck($tablef, "e") && (GENERAL::GetFileInfo($tablef, 'lines') == $rotamer->countValidPairTotal())) {
        GENERAL::GLog("Skipping pair EPIC '$term' calculation - has been completed already...");
        print "Skipping pair EPIC '$term' calculation - has been completed already...\n";
        next;
      }
      # Read appropriate topology and parameter tables
      my ($P, $T);
      if ($term =~ /EEF/) {
        my $Peef = CHARMM::readEEFParameters(DEFINITIONS::getParmFile("solvpar.inp"));
        $T = CHARMM::readCHARMMTopology(DEFINITIONS::getParmFile("toph19_eef1.inp"));
        $T = CHARMM::readCHARMMTopology(DEFINITIONS::getParmFile("patchtop19_eef1.inp"), $T);
        $P = CHARMM::readCHARMMParameters(DEFINITIONS::getParmFile("param19_eef1.inp"));
        CHARMM::expandCHARMMParameters($P, $T);
        foreach my $at (keys(%$Peef)) {
          foreach my $val (keys(%{$Peef->{$at}})) {
            $P->{eef1}{$at}{$val} = $Peef->{$at}{$val};
          }
        }
      } elsif ($term =~ /(VDW|DDE)/) {
        $T = CHARMM::readCHARMMTopology(DEFINITIONS::getParmFile("toph19.inp"));
        $T = CHARMM::readCHARMMTopology(DEFINITIONS::getParmFile("patchtop19.inp"), $T);
        $P = CHARMM::readCHARMMParameters(DEFINITIONS::getParmFile("param19_90.inp"));
        CHARMM::expandCHARMMParameters($P, $T);
      }

      # Divide the task among processes
      my $Nt = $rotamer->countResPairTotal();
      my $a = GENERAL::splitTasks($rotamer->countResPairTotal(), $NPROC_DEF);
      my $N = $a->[$PRANK_DEF]->{end} - $a->[$PRANK_DEF]->{beg} + 1;
      if ($N > 0) {
        my $pairnum = 0;
        # loop over first site and residue
        for (my $si = 0; $si < scalar(@$sitelist) - 1; $si++) {
          my $sitei = $sitelist->[$si];
          for (my $ri = 0; $ri < scalar(@{$sitei->{reslist}}); $ri++) {
            my $resi = $sitei->{reslist}->[$ri];
            # read one rotamer for this residue
            my $rpdbi = $FS_DEF->{rot_pdbf}; $rpdbi =~ s/%%/$si/; $rpdbi =~ s/%/$ri/;
            my $rotsi = PDB::new($rpdbi, "PDB", "CHARMM", 1);
            my $roti = $rotsi->{chain}->[0]->{res}->[0];
            my @sci = PDB::sidechain($roti);

            # loop over second site and residue
            for (my $sj = $si+1; $sj < scalar(@{$sitelist}); $sj++) {
              my $sitej = $sitelist->[$sj];
              for (my $rj = 0; $rj < scalar(@{$sitej->{reslist}}); $rj++) {
                $pairnum++;
                next if (($pairnum < $a->[$PRANK_DEF]->{beg}) || ($pairnum > $a->[$PRANK_DEF]->{end}));
                my $resj = $sitej->{reslist}->[$rj];

                my $energytab = "$FS_DEF->{dpair}/$si\_$ri/pairepic$pairnum.tab";
                $energytab = File::Spec->abs2rel($energytab);

                # if resuming a previous run, skip if this pair of residues has been completely treated
                if ($resume && (scalar(@$impdbs) == 0) && GENERAL::fileCheck($energytab, "e") && (GENERAL::GetFileInfo($energytab, 'lines') == $resi->{numvalidrot}*$resj->{numvalidrot})) {
                  print "Skipping pair $pairnum - has already been calculated in a previous run!\n";
                  next;
                }

                # see if the two residues are close enough to bother calculating
                my $d = sqrt(($resi->{base_coor}->[0]-$resj->{base_coor}->[0])**2 +
                            ($resi->{base_coor}->[1]-$resj->{base_coor}->[1])**2 +
                            ($resi->{base_coor}->[2]-$resj->{base_coor}->[2])**2);
                if ((scalar(@$impdbs) == 0) && ($d - $resi->{extent} - $resj->{extent} > $parm->{cutnb})) {
                  #print "Process $PRANK_DEF: pair $pairnum too distant - filling with zeros...\n";
                  my $fh = GENERAL::GetOutFH($energytab);
                  my $np = ($resi->{numvalidrot})*(($resj->{numvalidrot}));
                  $fh->printf("0\n" x $np);
                  close($fh);
                  next;
                }

                # read one rotamers of this residue
                my $rpdbj = $FS_DEF->{rot_pdbf}; $rpdbj =~ s/%%/$sj/; $rpdbj =~ s/%/$rj/;
                my $rotsj = PDB::new($rpdbj, "PDB", "CHARMM", 1);
                my $rotj = $rotsj->{chain}->[0]->{res}->[0];
                my @scj = PDB::sidechain($rotj);

                # call EPIC
                if (($term =~ /DDE/) && defined($parm->{rdie_type}) && ($parm->{rdie_type} eq "cc")) {
                  $parm->{rdie} = $parm->{cc_rdie}{$parm->{_npc}{$roti->{resname}}}{$parm->{_npc}{$rotj->{resname}}};
                  if (!defined($parm->{rdie})) { GENERAL::error("rdie not defined for residue pair $roti->{resname} x $rotj->{resname}"); }
                }
                my $inpi = $FS_DEF->{rot_outf}; $inpi =~ s/%%/$si/; $inpi =~ s/%/$ri/;
                my $inpj = $FS_DEF->{rot_outf}; $inpj =~ s/%%/$sj/; $inpj =~ s/%/$rj/;
                $epic = ENERGY::interfaceEPIC($energytab, $term, 'pair', $T, $P, $parm, \@sci, $inpi, 'atomsj', \@scj, 'inpj', $inpj, 'sf', 1.0, 'prog', $epic);
                # Calculate image i-j' and i'-j image terms, if necessary, and add them to the energy table
                my @impairs;
                for (my $im = 0; $im < scalar(@$impdbs); $im++) {
                  my $iminpi = $FS_DEF->{rot_ioutf}; $iminpi =~ s/%%%/$im/; $iminpi =~ s/%%/$si/; $iminpi =~ s/%/$ri/;
                  my $iminpj = $FS_DEF->{rot_ioutf}; $iminpj =~ s/%%%/$im/; $iminpj =~ s/%%/$sj/; $iminpj =~ s/%/$rj/;
                  $epic = ENERGY::interfaceEPIC($energytab, $term, 'pair', $T, $P, $parm, \@sci, $inpi, 'atomsj', \@scj, 'inpj', $iminpj, 'sf', 0.5, 'prog', $epic);
                  $epic = ENERGY::interfaceEPIC($energytab, $term, 'pair', $T, $P, $parm, \@sci, $iminpi, 'atomsj', \@scj, 'inpj', $inpj, 'sf', 0.5, 'prog', $epic);
                }
                if (time() - $time > 30*$LOCAL_TIMEOUT) { GENERAL::touchLocalSpace(); $time = time(); }
              }
            }
          }
        }
      }
      if (defined($epic)) { $epic->close(); $epic = undef; } # close to cause a flush

      # Wait until all processes are done before putting the data together
      GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF, $NFS_DELAY);

#       # Combine all tables into one sorted by rotamer
#       if ($PRANK_DEF == 0) {
#         my $pairnum = 0;
#         GENERAL::crm($tablef) if (-e $tablef);
#         for (my $si = 0; $si < scalar(@{$rotamer->{dslist}}) - 1; $si++) {
#           for (my $ri = 0; $ri < scalar(@{$rotamer->{dslist}->[$si]->{reslist}}); $ri++) {
#             for (my $sj = $si+1; $sj < scalar(@{$rotamer->{dslist}}); $sj++) {
#               for (my $rj = 0; $rj < scalar(@{$rotamer->{dslist}->[$sj]->{reslist}}); $rj++) {
#                 $pairnum++;
#                 my $energytab = "$FS_DEF->{dpair}/$si\_$ri/pairepic$pairnum.tab";
#                 GENERAL::csystem("cat $energytab >> $tablef");
#                 GENERAL::crm($energytab);
#               }
#             }
#           }
#         }
#         GENERAL::csystem("bzip2 $tablef");
#         # Document the created energy table
#         push(@{$en{$term}->{files}}, File::Spec->abs2rel($tablef, $FS_DEF->{dhome}));
#       }
                

      if ($PRANK_DEF == 0) {
        my $pairnum = 0;
        my $tfh = GENERAL::GetOutFH($tablef, "bzip2");
        for (my $s1 = 0; $s1 < scalar(@$sitelist) - 1; $s1++) {
          for (my $ri1 = 0; $ri1 < scalar(@{$sitelist->[$s1]->{reslist}}); $ri1++) {
            # Collect all files handles for tables involving this residue at this site
            # and all residues with which this residue interacts
            my @fh = ();
            my @f = ();
            my @res2 = ();
            my $res1 = $sitelist->[$s1]->{reslist}->[$ri1];
            for (my $s2 = $s1+1; $s2 < scalar(@{$sitelist}); $s2++) {
              for (my $ri2 = 0; $ri2 < scalar(@{$sitelist->[$s2]->{reslist}}); $ri2++) {
                push(@res2, $sitelist->[$s2]->{reslist}->[$ri2]);
                $pairnum++;
                my $energytab = "$FS_DEF->{dpair}/$s1\_$ri1/pairepic$pairnum.tab";
                push(@f, $energytab);
                push(@fh, GENERAL::GetInFH($energytab));
                for (my $im = 0; $im < scalar(@$impdbs); $im++) { }
              }
            }
            # Now merge all tables involving the current residue
            foreach my $rot1 (@{$res1->{rotamerlist}}) {
              next if ($rot1->{valid} == 0);
              for (my $ri2 = 0; $ri2 < scalar(@res2); $ri2++) {
                my $res2 = $res2[$ri2];
                foreach my $rot2 (@{$res2->{rotamerlist}}) {
                  next if ($rot2->{valid} == 0);
                  my $fh = $fh[$ri2];
                  my $line = <$fh>;
                  if (!defined($line)) { GENERAL::error("File $f[$ri2] seems to be too short!"); }
                  print $tfh $line;
                }
              }
            }
            # Close all opened files
            for (my $fi = 0; $fi < scalar(@fh); $fi++) {
              my $fh = $fh[$fi]; my $file = $f[$fi];
              if (defined(<$fh>)) { GENERAL::error("File $file seems to be too long!"); }
              close($fh);
            }
            # Delete smaller table files if needed
            GENERAL::csystem("rm $FS_DEF->{dpair}/$s1\_$ri1/*");
          }
        }
        close($tfh);
        # Document the created energy table
        push(@{$en{$term}->{files}}, File::Spec->abs2rel($tablef, $FS_DEF->{dhome}));
      }

      # Wait until all old files are erased before moving on
      GENERAL::my_MPI_Barrier($MPI_COMM_WORLD_DEF);
    }
    if ($PRANK_DEF == 0) {
      for (my $s = 0; $s < scalar(@{$rotamer->{dslist}}) - 1; $s++) {
        for (my $res = 0; $res < scalar(@{$rotamer->{dslist}->[$s]->{reslist}}); $res++) {
          GENERAL::crmdir("$FS_DEF->{dpair}/$s\_$res");
        }
      }
    }
    $parm->{rdie} = $tmprdie;
  }

  # Save term tracker
  if ($PRANK_DEF == 0) {
    GENERAL::SaveStruct($TERM_TRACKER, $FS_DEF->{term_trackerf});
  }
}


=head2

 Title   :  getSASolvation
 Usage   :  my $sol = ENERGY::getSASolvation(\@atoms, \@areas, $model);
 Function:  Calculates solvation energy based on a solvent exposed surface area model.
 Returns :  Solvation energy.
 Args    :  1. An array of atoms
            2. An array of surface areas corresponding to the array of atoms.
            3. Name of model.
=cut

sub getSASolvation {
  my $atoms = shift;
  my $areas = shift;
  my $model = shift;
  GENERAL::requireArgs($atoms, $areas, $model);
  if (scalar(@$atoms) != scalar(@$areas)) {
    GENERAL::error(sprintf("%d atoms and %d areas were passed", scalar(@$atoms), scalar(@$areas)));
  }
  
  my $ener = 0;
  if (uc($model) =~ /EIS/) {
#    my %ND; # Hash for oxygens on ASP and GLY and nitrogens on ARG
    my $C = 16; # Carbon
    my $NO = -6; # Neutral nitrogen or oxygen
    my $Om = -24; # Charged oxygen (most exposed in ASP and GLU)
    my $Np = -50; # Charged nitrogen (LYS and most exposed in ARG)
    my $S = 21; # Sulfur
    for (my $i = 0; $i < scalar(@$atoms); $i++) {
      my $atom = $atoms->[$i];
      next if ($atom->{atomname} =~ /^H/); # hydrogens don't count
      my $id = "$atom->{atomname} $atom->{residue}->{resname}";
      if ($id =~ /^C/) {
        # Carbon
        $ener += $C*$areas->[$i];
      } elsif ($id =~ /(^N |NE1 TRP|ND2 ASN|NE2 GLN|NE ARG|ND1 HIS|NE2 HIS|HD1 HSD|NE2 HSD|^O |OG SER|OG1 THR|OH TYR|OD1 ASN|OE1 GLN)/) {
        # Neutral oxygen and nitrogen
        $ener += $NO*$areas->[$i];
      } elsif ($id =~ /(OD.* ASP|OE.* GLU)/) {
        # Charged oxygens - split the parameter in 2 since only one is charged
        $ener += 0.5*$Om*$areas->[$i];
      } elsif ($id =~ /(NH.* ARG|ND1 HSC|NE2 HSC)/) {
        # Charged nitrogens - split the parameter in 2 since only one is charged
        $ener += 0.5*$Np*$areas->[$i];
#      } elsif ($id =~ /(OD.* ASP|OE.* GLU|NH.* ARG)/) {
#        # These can be either neutral or charged - will decide later
#        my $str = PDB::resStr($atom->{residue});
#        if (defined($ND{$str}{$atom->{atomname}})) { GENERAL::error("Atom $str not unique"); }
#        $ND{$str}{$atom->{atomname}} = $areas->[$i];
      } elsif ($id =~ /NZ LYS/) {
        # Unambiguously charged nitrogens
        $ener += -$Np*$areas->[$i];
      } elsif ($id =~ /^S/) {
        $ener += $S*$areas->[$i];
      } else {
        GENERAL::warning("Could not find parameter for atom $atom->{atomname} in residue $atom->{residue}->{resname}");
      }
    }
#    for my $as (keys(%ND)) {
#      my ($i, $j);
#      if ($as =~ /ASP/) {
#        $i = "OD1"; $j = "OD2";
#      } elsif ($as =~ /GLU/) {
#        $i = "OE1"; $j = "OE2";
#      } elsif ($as =~ /ARG/) {
#        $i = "NH1"; $j = "NH2";
#      }
#      if (!(defined($ND{$as}{$i}) && defined($ND{$as}{$j}))) {
#        GENERAL::error("No pair ($i, $j) found in residue $as to determine which one is charged!");
#      }
#      if ($as =~ /ASP/) {
#        $ener += $NO*GENERAL::min($ND{$as}{$i}, $ND{$as}{$j}) + $Om*GENERAL::max($ND{$as}{$i}, $ND{$as}{$j});
#      } elsif ($as =~ /GLU/) {
#        $ener += $NO*GENERAL::min($ND{$as}{$i}, $ND{$as}{$j}) + $Om*GENERAL::max($ND{$as}{$i}, $ND{$as}{$j});
#      } elsif ($as =~ /ARG/) {
#        $ener += $NO*GENERAL::min($ND{$as}{$i}, $ND{$as}{$j}) + $Np*GENERAL::max($ND{$as}{$i}, $ND{$as}{$j});
#      }
#    }
    # We need to convert to kcal/mol
    $ener /= 1000;
  } else {
    GENERAL::error("Unknown model '$model'");
  }
  
  return $ener;
}

=head2 

 Title   :   parseSelfEnergy
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut


=head2 

 Title   :  parseSubrotSelfEnergy
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut



=head2 

 Title   : mergeEnergyBin
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut


=head2 

 Title   :  buildHybridPotential
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut





=head2 

 Title   :  calculateASP 
 Usage   : 
 Function:  read in buried surface area  
 Returns : 
 Args    : 
 Note    : the input file should inculde: 
           (A) scaled buried area between template and each rotamer
           (B) scaled buried area between each rotamer pair
           (C) scaling parameter is 0.62

=cut

sub calculateASP {

    my $inp=shift;
    my $out=shift;

    $out='asp.tab' if (! defined $out);

    my $inpf=&GENERAL::GetInFH($inp);
    my $outf=&GENERAL::GetOutFH($out);

    while (<$inpf>) {
    my @line=split(" ", $_);

    # first for C  
    # second for amind N 
    # third for amine N
    # fourth for oxygen carbonyl or hydroxyl
    # fifth for oxygen for carboxyl
    # sixth for thiol
    # note: ASP value methione S is treated as 0. 
    # parameters are taken from eisenberg's ASP value of vacuum to water 
    # the reported energy is the solvation penalty

    my $energy= -0.012*$line[0]+0.116*$line[1]+0.186*$line[2]+0.116*$line[3]+0.175*$line[4]+0.018*$line[5];

    $outf->printf("%8.4f\n", $energy);
    }
}
    


=head2

 Title   :  combineResiduePairTables
 Usage   :  ENERGY::combineResiduePairTables($rotamer, sub {}, $dest_file);
 Function:  Combines individual resi-resj pair tables.
 Returns :  Nothing.
 Args    :  1. Rotamer object
            2. A function reference of a function that returns the name of the res-res pair table
               given site index i, residue index i, site index j, residue index j, residue pair number,
               where i refers to the first residue in a pair, j to the second, and residue pair number
               is index of the residue pair in the natural {site i{res i{site j{res j}}}} enumeration
               order (site and residue indices start with 0, pair numbers start with 1).
            3. Optional: destination file name or file handle. If not specified, the resulting energy
               table will be read into memory and returned.
            4. Optional: flag specifying whether the individual res-res pair tables are to be deleted.
               If 0 (or not specified) they will not be. If 1, they will be.
=cut

sub combineResiduePairTables {
  my $rotamer = shift;
  my $ptn_fp = shift;
  GENERAL::requireArgs($rotamer, $ptn_fp);
  my $df = shift;
  my $ofh = undef;
  if (defined($df)) {
    if (ref($df)) { $ofh = $df; }
    else { $ofh = GENERAL::GetOutFH($df); }
  }
  my $del = shift;
  $del = 0 if (!defined($del));
  my @pt;

  # Combine all tables into one sorted by rotamer
  my $pairnum = 0; my $sitelist = $rotamer->{dslist};
  for (my $s1 = 0; $s1 < scalar(@$sitelist) - 1; $s1++) {
    for (my $ri1 = 0; $ri1 < scalar(@{$sitelist->[$s1]->{reslist}}); $ri1++) {
      # Collect all files handles for tables involving this residue at this site
      # and all residues with which this residue interacts
      my @fh = ();
      my @f = ();
      my @res2 = ();
      my $res1 = $sitelist->[$s1]->{reslist}->[$ri1];
      for (my $s2 = $s1+1; $s2 < scalar(@{$sitelist}); $s2++) {
        for (my $ri2 = 0; $ri2 < scalar(@{$sitelist->[$s2]->{reslist}}); $ri2++) {
          push(@res2, $sitelist->[$s2]->{reslist}->[$ri2]);
          $pairnum++;
#          my $energytab = "$FS_DEF->{dpair}/$s1\_$ri1/pairepic$pairnum.tab";
          my $energytab = $ptn_fp->($s1, $ri1, $s2, $ri2, $pairnum);
          push(@f, $energytab);
          push(@fh, GENERAL::GetInFH($energytab));
        }
      }
      # Now merge all tables involving the current residue
      foreach my $rot1 (@{$res1->{rotamerlist}}) {
        next if ($rot1->{valid} == 0);
        for (my $ri2=0; $ri2 < scalar(@res2); $ri2++) {
          my $res2 = $res2[$ri2];
          foreach my $rot2 (@{$res2->{rotamerlist}}) {
            next if ($rot2->{valid} == 0);
            my $fh = $fh[$ri2];
            my $line = <$fh>;
            if (!defined($line)) { GENERAL::error("File $f[$ri2] seems to be too short!"); }
            if (defined($ofh)) { print $ofh $line; }
            else { push(@pt, $line + 0.0); }
          }
        }
      }
      # Close all opened files
      for (my $ri2 = 0; $ri2 < scalar(@res2); $ri2++) {
        my $fh = $fh[$ri2];
        if (defined(<$fh>)) { GENERAL::error("File $f[$ri2] seems to be too long!"); }
        close($fh);
      }
      # Delete smaller table files if needed
      if ($del) { foreach my $energytab (@f) { GENERAL::csystem("rm $energytab"); } }
    }
  }

  if (defined($df) && !ref($df)) { close($ofh); }
  if (!defined($df)) { return \@pt; }
}


1;
