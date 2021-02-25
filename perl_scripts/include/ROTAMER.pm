# ROTAMER package
# read/write/convert ROTAMER  info
#
#
# Design Module: ROTAMER
#
# POD documentation - main docs before the code

=head1 NAME ROTAMER
    
  ROTAMER object, with features
    
=head1 SYNOPSIS
    
    someone needs to write it
    
=head1 DESCRIPTION

    1. Read/Write Backbone-dependent Rotamer library

    2. Read/Write Backbone-independent Rotamer library
 
    3. Read/Write Rotamer_input

    4. Read/Write Rotamer_input2

    5. Read/Write Rotamer_input3

    6. Subrotamer extension

    7. Rotamer Analysis on Crystal Structure

    8. Compare calculated rotamer with experimental rotamer

    9. GetRotamerState 

    10. Write Rotamer-Self-State file

    11. Generate discretely, randomly and gaussian distributed subrotamer 
        library 


    Data Structure:

    (1) reference to rotamer files:   

    rot_inp1
    
    rot_inp2
    
    rot_inp3

    (2) reference to library file
    
    bbdeplib
    bbindeplib
    
    (3)
    


=head1 EXAMPLES

    someone needs to write it

=head2 Reporting Bugs

    Jiangang Chen

=head1 AUTHOR 

    Jiangang Chen 
    version 0.0.1: Feb 20th, 2002 

=head1 APPENDIX

    The rest of the documentation details each of the object methods. Internal methods are usually preceded with a "_".

=cut


package ROTAMER;
#$DEF_VER_CONTROL->{ver}->{ROTAMER} = "1.0";
$DEF_VER_CONTROL->{date}->{ROTAMER} = "06/17/04";

use strict;
use locale;
use GENERAL;
use PDB;
use DEFINITIONS;



=head2  NAME -- new

 Title   : new
 Usage   : ROTAMER::new('rot_input'); 
 Function: Constructor for an rotamer object
 Returns : none 
 Args    : (1) rot.inp file name 
           (2) rotstate file name: optional

=cut


sub new {

  # first initialize some data
  my $self={};
  $self->{dslist}=[];
  $self->{chidef}=undef;
  $self->{sigma}=undef;
  $self->{subrotamer}={};
  $self->{subrotamer}->{isValid}=0;
  $self->{subrotamer}->{number}=undef;

  bless $self;

  # then read in arguments

  my $rotinp = shift;
  if (defined $rotinp) {
    # Call readRotInput before calling readChiDef in case a custom
    # rotamer library directory is specified in rot.inp, in which case
    # the appropriate chi definition file must be used.
    $self->readRotInput($rotinp);
    $self->{rif} = $rotinp;
    $self->readChiDef();
    my $rotstate=shift;
    $self->readRotState($rotstate) if (defined $rotstate);
  }


  return $self;

}


=head2 

 Title   :   readRotlib
 Usage   :   $lib253ref=$rotamer->readRotlib(253);
             $bbindepref=$rotamer->readRotlib(0);
             $customref=$rotamer->readRotlib('user_specified.lib');
 Function:   Given a backbone index, return the reference to the corresponding rotamer library.
             If backbone index is 0, assumes a backbone independent rotamer library, which means
             that probabilities are not parsed. If the rotamer index is not an integer string,
             it is assumed to be the name of the rotamer library file.
 Returns :   reference to rotamer library
 Args    :   1. backbone index (or name of rotamer library file)

=cut

sub readRotlib {

  my $self=shift;
  my $bbindex=shift;
  if (! defined $bbindex) { GENERAL::error("backbone index not specified!"); }
  my $libfile;
  if (!GENERAL::isInteger($bbindex)) {
    $libfile = $bbindex;
  } else {
    $libfile = "$ROTLIB_DEF/library_$bbindex";
  }

  my $libfh = GENERAL::GetInFH($libfile);
  my $libref={};
  my $numchi;
  my $rotnum;
  my $resname;

  while (<$libfh>) {

    # skip the empty and commenteds lines
    next if (/^\s*$/);
    next if (/^\s*\#/);
    if (/^RESI/ ) {
      my @line=split(" ", $_);
      $resname = $line[1];
      $numchi = $line[3];
      $rotnum = $line[2];
      $libref->{$resname} = {};
      if (($numchi == 0) && ($rotnum == 0)) {
        $libref->{$resname}->{numrot} = 1;
        $libref->{$resname}->{numvalidrot} = 1; # keeps track of the number of rotamers which haven't been eliminated by self calculations
        $libref->{$resname}->{resname} = $resname;
        $libref->{$resname}->{numchi} = 0;
        my $rec={}; $rec->{probability} = 1.0; $rec->{chi} = []; $rec->{valid} = 1;
        push (@{$libref->{$resname}->{rotamerlist}}, $rec);
      } else {
        $libref->{$resname}->{numrot} = $rotnum;
        $libref->{$resname}->{numvalidrot} = $rotnum; # keeps track of the number of rotamers which haven't been eliminated by self calculations
        $libref->{$resname}->{resname} = $resname;
        $libref->{$resname}->{numchi} = $numchi;
      }
    }

    if ( ! /^RESI/) {
      my @line = split (" ", $_);
      my $rec={};
#      if (($bbindex != 0) && (scalar(@line) > $numchi)) {
      if (scalar(@line) > $numchi) {
        $rec->{probability} = $line[$numchi];
        pop @line;
      } else {
        $rec->{probability} = undef;
      }
      $rec->{chi} = \@line;
      $rec->{valid} = 1;
      push (@{$libref->{$resname}->{rotamerlist}}, $rec);
    }
  }

  # Take care of the number of rotamers
  foreach my $res (keys(%$libref)) {
    $libref->{$res}->{numrot} = scalar(@{$libref->{$res}->{rotamerlist}});
    $libref->{$res}->{numvalidrot} = scalar(@{$libref->{$res}->{rotamerlist}});
  }

  return $libref;
}


=head2 

 Title   :  readSubRotlib
 Usage   :  
 Function: 
 Returns : 
 Args    : 

=cut



sub readSubRotlib {

    my $self=shift;
    my $bbindex=shift;
    if (!defined $bbindex) { die "Error in readSubRotlib: backbone index not specified!\n"; }
    #print $bbindex, "\n";
    my $method=shift;
    if (!defined $method) { die "Error in readSubRotlib: method not specified!\n"; }

    my $lib = "sublibrary_$bbindex";
    die "Error in readSubRotlib: you can't use a subrotamer library before calling useSubRotamer!\n" 
	if ($self->{subrotamer}->{isValid} ==0);

    if (! -e $lib) {
	$self->writeSubRotlib($bbindex, $method);
    }

    # Read subrotamer library from the specified file assuming backbone index 0,
    # so that probabilities are set to undef (since we don't have them in subrotamer
    # libraries).
    return $self->readRotlib(0, $lib);
}


=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut


sub useSubRotamer {
    my $self=shift;

    $self->{subrotamer}->{isValid}=1;

    # here i set allowed subrotamers for ARG and LYS as 3
    # which means i only changes chi1 since ARG and LYS 
    # have 81 rotamers already. 

    my %par= ( 'ALA' =>  0,
	       'ARG' =>  9,
	       'ASN' =>  9, 
	       'ASP' =>  9, 
	       'CYS' =>  3,
               'CYD' =>  0,
	       'GLN' =>  9,
	       'GLU' =>  9,
	       'HIS' =>  9,
	       'HSD' =>  9,
	       'ILE' =>  9,
	       'LEU' =>  9,
	       'LYS' =>  9,
	       'MET' =>  9,
	       'PHE' =>  9,
	       'PRO' =>  9,
	       'SER' =>  3,
	       'THR' =>  3,
	       'TRP' =>  9,
	       'TYR' =>  9,
	       'VAL' =>  3,
	       'PRO' =>  0
	       );
    $self->{subrotamer}->{number}=\%par;

}


=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut


sub setSubRotamerNumber {
    my $self=shift;
    my $resname=uc(shift);
    my $number=shift;

    if ( !exists $self->{subrotamer}->{number}->{$resname} ) {
      die "Error in setSubRotamerNumber: couldn't find entry for $resname.\n"; 
    } else {
      $self->{subrotamer}->{number}->{$resname}=$number;
    }
}



=head2

 Title   : getDsite
 Usage   : $rotamer->getDsite($chain_id, $iresnum)
 Function: Returns a reference to the design site corresponding
           to the specified chain id and site number. Gracefully
           quits if the site is not found.
 Returns :
 Args    : 1. Chain id.
           2. Residue number.
           3. Optional logical flag. If true (default) quits when the
              specified site is not found. If false, simply returns -1.

=cut


sub getDsite {
  my $self = shift;
  my $cid = shift;
  my $iresnum = shift;
  my $halt = shift;
  if (!defined($halt)) { $halt = 1; }

  foreach my $ds (@{$self->{dslist}}) {
    if (($ds->{chain} eq $cid) && ($ds->{iresnum} == $iresnum)) {
      return $ds;
    }
  }
  if ($halt) {
    GENERAL::error("Could not find site $cid$iresnum in the ROTAMER object!");
  } else {
    return -1;
  }
}


=head2

 Title   : getResInDsite
 Usage   : ROTAMER::getResInDsite($ds, $resname)
 Function: Returns a reference to the named residue at the given design. Gracefully
           quits if the residue is not found.
 Returns :
 Args    : 1. design site.
           2. residue name.
           3. Optional logical flag. If true (default) quits when the
              specified site is not found. If false, simply returns -1.

=cut


sub getResInDsite {
  my $ds = shift;
  my $resname = shift;
  my $halt = shift;
  if (!defined($halt)) { $halt = 1; }

  foreach my $res (@{$ds->{reslist}}) {
    if ($res->{resname} eq $resname) { return $res; }
  }
  if ($halt) {
    GENERAL::error("Could not find residue '$resname' at design site $ds->{chain}$ds->{iresnum} in the ROTAMER object!");
  } else {
    return -1;
  }
}


=head2 

 Title   : getResInRotlib
 Usage   : $rotamer->getResInRotlib($ds_rec->{rotlibref}, $resname)
 Function: Returns the rotamer library entry corresponding to
           the specified amino acid. Gracefully quits if the residue
           does not exist in the rotamer library.
 Returns : 
 Args    : 1. Reference to the rotamer library for a particular site.
           2. Name of amino acid in question.

=cut


sub getResInRotlib {
    my $self = shift;
    my $libref = shift;
    my $resname = shift;
    GENERAL::error("Can't find $resname in Rotamer Library!") if (!exists $libref->{$resname});
    return $libref->{$resname};
}

=head2 

 Title   : readChiDef
 Usage   :
 Function: Reads dihedral angle definition file. Stores the definitions in
           the hash table $self->{chidef}, where entries are hashed by residue
           name. Each entry is an array of definitions (first one corresponding
           to chi1, the second one to chi2, and so on). Each definition is simply
           an array of atom names corresponding to the atoms forming the dihedral
           atom (example: N CA CB CG1).
 Returns :
 Args    :

=cut


sub readChiDef {
  my $self=shift;
  $self->{chidef} = undef if (defined $self->{chidef});
  my $sanity_check = shift;
  $sanity_check = 1 if (!defined($sanity_check));
  my $fname = shift;

  if (!defined($fname)) {
    $fname = "chi_definition";
    if (-e "$CUSTOM_PARAM_DEF/$fname") {
      $fname = File::Spec->rel2abs("$CUSTOM_PARAM_DEF/$fname");
    } elsif (-e "$ROTLIB_DEF/$fname") {
      $fname = File::Spec->rel2abs("$ROTLIB_DEF/$fname");
    } else {
      GENERAL::error("Could not find $fname file!");
    }
  }
  my $fh = GENERAL::GetInFH($fname);

  my $resname;
  my $chidef={};

  while (<$fh>) {
    s/\!.*$//; # remove comments
    s/\#.*$//; # remove comments
    next if (/^$/);
    if (/^RESI/) {
      my @line=split(" ", $_);
      $resname= $line[1];
      $chidef->{$resname}=[];
    } else {
      my @line=split(" ", $_);
      pop @line;
      push (@{$chidef->{$resname}}, \@line);
    }
  }
  close($fh);
  $self->{chidef}= $chidef;

  # Check to make sure that the chi definitions agree with
  # the entries in the rotamer library
  if ($sanity_check) {
    foreach my $site (@{$self->{dslist}}) {
      foreach my $res (@{$site->{reslist}}) {
        foreach my $rot (@{$res->{rotamerlist}}) {
          if (scalar(@{$rot->{chi}}) != scalar(@{$chidef->{$res->{resname}}})) {
            GENERAL::error("For residue $res->{resname} at $site->{chain}$site->{iresnum} the number of ".
                           "chi angle definitions is " . scalar(@{$chidef->{$res->{resname}}}) . " and the ".
                           "number of angles listed in the rotamer library is " . scalar(@{$rot->{chi}}));
          }

          foreach my $chi (@{$chidef->{$res->{resname}}}) {
            if (scalar(@$chi) != 4) {
              GENERAL::error("Invalid chi definition for residue $res->{resname}: (" . join(", ", @$chi) . "). ".
                             "It takes exactly 4 atoms to defined a dihedral angle!");
            }
          }
        }
      }
    }
  }
}




=head2 

 Title   : writeRotInput
 Usage   : $rotamer->writeRotInp(...)
 Function: Writes a rot.inp file given the ROTAMER object.
 Returns : nothing
 Args    : 1. file name of stride output
           2. flag specifying whether only the WT residue should
              be included in the rot.inp or all residues. REPACK
              means only WT residues. FULL means all residues.
           3. the name of the output file (such as rot.inp).
           4. optional: wild type PDB object. If specified, the
              WT rotamer will be included in the output file.

=cut

sub writeRotInput {
    my $self = shift;
    my $fname1 = shift;
    my $option = shift;
    my $fname2 = shift;
    GENERAL::requireArgs($self, $fname1, $option, $fname2);
    my $npdb = shift;

    my @seq; # optional custom sequence
    if (defined($option) && (uc($option) =~ /SEQUENCE/)) {
      my $seq = $option; $option = "PACKSEQUENCE";
      $seq =~ s/SEQUENCE//; $seq = GENERAL::Trim($seq);
      @seq = split(" ", $seq);
    }
    
    my @alph = ('ALA', 'VAL', 'LEU', 'ILE', 'PHE', 'TRP', 'ASN', 'THR', 'SER', 'TYR',
                'HIS', 'CYS', 'MET', 'GLN', 'ASP', 'GLU', 'ARG', 'LYS', 'PRO', 'GLY');
    my @calph;
    if (uc($option) =~ /CUSTOM_ALPH\s+(.+)$/) {
      @calph = split(";", GENERAL::Trim($1));
      for (my $ii = 0; $ii < scalar(@calph); $ii++) { my @tmp = split(" ", GENERAL::Trim($calph[$ii])); $calph[$ii] = \@tmp; }
    }

    my $initinp = GENERAL::GetInFH($fname1);
    my $newinp = GENERAL::GetOutFH($fname2);

    my $cc = 0;
    while (<$initinp>) {
#      next if ( /GLY/);   # you don't want to repack GLY, no side chain
#      next if ( /PRO/ );  # also proline. this residue is too weird
      chomp($_);
      my $line = $_;
      my @line = split(" ", $line);
      if (defined($option) && (uc($option) =~ /SEQUENCE/)) {
        GENERAL::assert($cc < scalar(@seq), "Custom sequence provided too short!");
        $line =~ s/\s+([^\s]+)$/;  $seq[$cc]/;
      } elsif (defined($option) && (uc($option) =~ /CUSTOM_ALPH/)) {
        $line =~ s/\s+([^\s]+)$//;
      } else {
        $line =~ s/\s+([^\s]+)$/;  $1/;
      }
#      $line = "$line[0]   $line[1]  $line[2];  $line[3]";
      my $resname = $line[3];
      my $cid = $line[0];
      my $iresnum = $line[1];
      # Do we need to include the WT rotamer?
      if (defined($npdb)) {
        my $chain = $npdb->getChain($cid);
        my $res = $chain->{res}->[$iresnum-1];
        # Find the value of each dihedral angle for this residue in the WT structure
        my $chidef = $self->{chidef}->{$resname};
        if (!defined($chidef)) { GENERAL::error("No chi definition found for residue $resname!"); }
        my @chis;
        foreach my $chi (@{$chidef}) {
          if (scalar(@$chi) != 4) { GENERAL::error("Invalid chi angle definition: " . join(", ", @$chi)); }
          push(@chis, PDB::getDihedral($res, $chi->[0], $chi->[1], $chi->[2], $chi->[3]));
        }
        if (scalar(@chis) > 0) {
          $line .= " (";
          for (my $i=0; $i < scalar(@chis); $i++) { $line .= sprintf("%s%.2f", " "x($i != 0), $chis[$i]); }
          $line .= sprintf(")");
        }
      }
      if (uc($option) =~ /PACK/ ) {
        $newinp->print("$line\n");
      } elsif ( uc($option) =~ /FULL/ ) {
        foreach my $i (@alph) {
          if (uc($resname) eq $i) {
            next;
          } else {
            my $k = "; ".$i;
            $line .= $k;
          }
        }
        $newinp->print($line, "\n");
      } elsif (uc($option) =~ /CUSTOM_ALPH/) {
        GENERAL::assert(scalar(@calph) != 0, "Fewer sites in custom alphabet spec than in the site file!");
        foreach my $i (@{$calph[0]}) {
          $line .= ("; ".$i);
        }
        $newinp->print($line, "\n");
        shift @calph;
      }
      $cc++;
    }
    GENERAL::assert($cc == scalar(@seq), "Custom sequence provided too long!") if (defined($option) && (uc($option) =~ /SEQUENCE/));
    GENERAL::assert(scalar(@calph) == 0, "More sites in custom alphabet spec than in the site file!") if (uc($option) =~ /CUSTOM_ALPH/);
}




=head2 

 Title   :  
 Usage   : 
 Function: 
 Returns : 
 Args    : 
 Version : 

=cut



sub readRotInput {
  my $self = shift;
  my $fname = shift;

  if (!defined $fname) { GENERAL::error("Rotamer input file not specified!"); }
  my $rotinp1 = GENERAL::GetInFH($fname);

  while (<$rotinp1>) { last if ($_ !~ /^\s*$/); }
  # is this a definition of rotamer library path?
  if (/rotlib_path/) {
    $_ =~ /^rotlib_path\s*=\s*(\S+)\s*$/;
    $ROTLIB_DEF = $1;
  } else {
    GENERAL::error("First non-empty line in $fname is expected to contain path to global rotamer library!");
  }
  while (<$rotinp1>) {
    chomp();
    next if (/^\s*$/) ;
    my $rec={};
    my @line = split(";", $_);
    my @sinfo = split(" ", $line[0]);
    my @rinfo = @line[1..(scalar(@line)-1)];

    # now we parse the rot_input1
    $rec->{chain} = $sinfo[0];       # get chain id
    $rec->{iresnum} = $sinfo[1];     # get residue number, note: we use iresnum because
                                     # we work with renumbered pdb files
    $rec->{bbindex} = $sinfo[2];     # get backbone index which is determined by the phi, psi value
    if ($self->{subrotamer}->{isValid} == 0) {
      # get proper rotamer library and store it as reference
      $rec->{rotlibref} = $self->readRotlib($rec->{bbindex});
    } elsif ($self->{subrotamer}->{isValid} == 1) {
      $rec->{rotlibref} = $self->readSubRotlib($rec->{bbindex},'discrete');
    }

    if (scalar(@rinfo) == 0) {
      GENERAL::error("Error parsing $fname. Site $rec->{chain}$rec->{iresnum} appears to have no allowed residues!");
    }
    foreach my $r (@rinfo) {
      # parse each residue
      my $resname = GENERAL::Trim($r);
      if ($r =~ m/\((.*)\)/) {
        # if additional rotamers specified, collect information about them
        my @rots = split(",", $1);
        # now, extract residue name
        $r =~ s/\(.*\)//;
        $resname = GENERAL::Trim($r);
        # retrieve the residue's rotamers from rotamer library, stored as reference in array
        push (@{$rec->{reslist}}, $self->getResInRotlib($rec->{rotlibref}, $resname));
        # put in the additional rotamers
        my $res = $rec->{reslist}->[scalar(@{$rec->{reslist}}) - 1];
        foreach my $rot (@rots) {
          my @angles = split(" ", $rot);
          if (scalar(@angles) != $res->{numchi}) {
            GENERAL::error("At least one of additional rotamer specified for residue $res->{resname} in design " .
                           "site $rec->{chain}$rec->{iresnum} has an incorrect number of chi angles!");
          }
          # make sure chi angles are numeric
          foreach my $a (@angles) {
            if (!GENERAL::isNumeric($a)) { GENERAL::error("Invalid chi angle \"$a\" specified for residue $res->{resname} " .
                                                          "in design site $rec->{chain}$rec->{iresnum}!\n"); }
          }
          $res->{numrot}++;
          $res->{numvalidrot}++;
          my $newrot = {};
          $newrot->{probability} = 0.0; # by default, the additional rotamers get a probability of 0
          $newrot->{valid} = 1;
          $newrot->{chi} = \@angles;
          push(@{$res->{rotamerlist}}, $newrot);
        }
      } else {
        # if no additional rotamers specified
        # retrieve each residue's rotamers from rotamer library, stored as reference in array
        push (@{$rec->{reslist}}, $self->getResInRotlib($rec->{rotlibref}, $resname));
      }

      # set parent pointer for each rotamer in each residue of this residue
      my $res = $self->getResInRotlib($rec->{rotlibref}, $resname);
      foreach my $rot (@{$res->{rotamerlist}}) { $rot->{residue} = $res; Scalar::Util::weaken($rot->{residue}); }
    }
    # set parent pointer for each residue in this site
    foreach my $res (@{$rec->{reslist}}) { $res->{dsite} = $rec;  Scalar::Util::weaken($res->{dsite}); }

    # set parent pointer for design site
    $rec->{parent} = $self;
    Scalar::Util::weaken($rec->{parent});

    # finally, we put it into $self->{dslist} data structure
    push (@{$self->{dslist}}, $rec);

  }
}



=head2

 Title   : writeDesignAAInput
 Usage   : $rotamer->writeDesignAAInput();
 Function: Writes the amino acid input file for the design
           code.
 Returns : nothing
 Args    : Name of the file to output to.

=cut


sub writeDesignAAInput {
  my $self=shift;
  my $fname=shift;
  if ( !defined $fname ) { die "Error in writeDesignAAInput - amino acid input file name not specified!\n"; }

  my $aainp = GENERAL::GetOutFH($fname);

  foreach my $site (@{$self->{dslist}}) {
    foreach my $i (@{$site->{reslist}}) {
      $aainp->printf("%3s ", $i->{resname});
    }
    $aainp->print("\n");
  }
  close($aainp);
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut


sub writeDesignRotInput1 {
    
    my $self = shift;
    my $fname = shift;
    if (!defined($fname)) { GENERAL::error("rotamer input file name not specified!"); }

    my $rotinp2 = GENERAL::GetOutFH($fname);

    foreach my $site (@{$self->{dslist}}) {
      $rotinp2->printf("%1s%5d ", $site->{chain}, $site->{iresnum});
      foreach my $i (@{$site->{reslist}}) {
          $rotinp2->printf("%3s ", $i->{resname});
      }
      $rotinp2->print("\n");
    }
    close($rotinp2);
}



sub writeDesignRotInput2 {
    
    my $self = shift;
    my $fname = shift;
    if (!defined($fname)) { GENERAL::error("rotamer input file name not specified!"); }

    my $rotinp2 = GENERAL::GetOutFH($fname); 

    foreach my $site (@{$self->{dslist}}) {
        $rotinp2->printf("%1s%5d", $site->{chain}, $site->{iresnum});
        foreach my $i (@{$site->{reslist}}) {
            $rotinp2->printf(" $i->{resname} $i->{numrot}");
        }
        $rotinp2->print("\n");
    }
    close($rotinp2);
}





=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub writeDesignRotInput3 {
    
    my $self = shift;
    my $fname = shift;
    if (!defined($fname)) { GENERAL::error("rotamer input file name not specified!"); }
    my $rotinp3 = GENERAL::GetOutFH($fname);

    foreach my $site (@{$self->{dslist}}) {
      $rotinp3->printf("%1s%5d", $site->{chain}, $site->{iresnum});
      foreach my $i (@{$site->{reslist}}) {
        $rotinp3->printf(" $i->{resname} $i->{numvalidrot}");
      }
      $rotinp3->print("\n");
    }
    close($rotinp3);
}


=head2 

 Title   : writeRotState
 Usage   : $rot->writeRotState($filename)
 Function: Writes the rotamer states file - each rotamer
           gets either a 0 or 1 depending on whether it was
           eliminated or not.
 Returns : Nothing
 Args    : File name to write to.

=cut

sub writeRotState {
  my $self=shift;
  my $fname=shift;
  if (!defined $fname) { die "Error in writeRotState: rotamer states file name not specified!\n"; }
  my $rotstate=&GENERAL::GetOutFH($fname);

  # loop over the designed sites
  foreach my $site (@{$self->{dslist}}) {
    # loop over residue lists
    my $num=0;
    foreach my $i (@{$site->{reslist}}) {
      # loop over the rotamer states
      #  my $num=0;
      foreach my $j (@{$i->{rotamerlist}}) {
        $rotstate->printf("%2d",$j->{valid});
        $num++;
        if ($num == 20) {
          $num=0;
          $rotstate->print("\n");
        }
      }
      $rotstate->print("\n\n");
    }
    $rotstate->print("\n");
  }
}


=head2

 Title   : writeDesignRotState
 Usage   : $rot->writeDesignRotState($filename)
 Function: Writes the rotamer states file - each rotamer
           gets either a 0 or 1 depending on whether it was
           eliminated or not. The difference between this function
           and writeRotState is that here rotamers for one site are printed
           in each paragrap, while in writeRotState rotamers for one
           residue are printed in each paragraph with more lines separating
           sites from one another.
 Returns : Nothing
 Args    : File name to write to.

=cut

sub writeDesignRotState {
  my $self = shift;
  my $fname = shift;
  if (!defined($fname)) { GENERAL::error("rotamer states file name not specified!"); }
  my $rotstate = GENERAL::GetOutFH($fname);

  # loop over the designed sites
  foreach my $site (@{$self->{dslist}}) {
    # loop over residue lists
    my $num=0;
    foreach my $i (@{$site->{reslist}}) {
      # loop over the rotamer states
      #  my $num=0;
      foreach my $j (@{$i->{rotamerlist}}) {
        if ($num == 20) {
          $num=0;
          $rotstate->print("\n");
        }
        $rotstate->printf("%2d",$j->{valid});
        $num++;
      }
    }
    $rotstate->print("\n\n");
  }
}

=head2

 Title   : readRotState
 Usage   : $rot->readRotState($filename)
 Function: Reads in the rotamer state file produced by writeRotState
           and marks each rotamer in the object as either valid or invalid.
 Returns : Nothing
 Args    : File name to read from.

DEPRECATED DEPRECATED DEPRECATED DEPRECATED DEPRECATED DEPRECATED DEPRECATED
DEPRECATED DEPRECATED DEPRECATED DEPRECATED DEPRECATED DEPRECATED DEPRECATED
DEPRECATED DEPRECATED DEPRECATED DEPRECATED DEPRECATED DEPRECATED DEPRECATED
Use readDesignRotState instead.

=cut

sub readRotState {
    my $self=shift;
    my $fname=shift;

    if (!defined $fname) { die "Error in readRotState: rotamer states file not specified!\n"; }
    my $rotstate = GENERAL::GetInFH($fname);
    my @states=();

    # collect rotamer state info
    my $res = 0; my $new = 1;
    while (<$rotstate>) {
      $_ = GENERAL::Trim($_);
      if ($_ ne "") {
        my @lines=split(" ", $_);
        push (@{$states[$res]}, @lines);
        $new = 0;
      } elsif ($new == 0) {
        $res++;
        $new = 1;
      }
    }
    close($rotstate);

    $res = 0;
    # loop over the designed sites
    foreach my $site (@{$self->{dslist}}) {
      # loop over residue lists
      foreach my $i (@{$site->{reslist}}) {
        if (!defined($states[$res])) {
          die "Error in readRotState: there seem to be only $res residues in $fname
          but residue $site->{chain}$site->{iresnum} $i->{resname} is the " . ($res+1) .
          "-th design residue according to the ROTAMER object (created with $FS_DEF->{rotinpf}).\n";
        }
        if (scalar(@{$states[$res]}) != scalar(@{$i->{rotamerlist}})) {
          die "Error in readRotState: number of rotamers for residue $site->{chain}$site->{iresnum} $i->{resname}
          as defined in $fname (" . scalar(@{$states[$res]}) . ") conflicts with that according to
          the ROTAMER object created with $FS_DEF->{rotinpf} (" . scalar(@{$i->{rotamerlist}}). ")\n";
        }

        my $validrotnum=0;
        my $rot = 0;
        # loop over the rotamer states
        foreach my $j (@{$i->{rotamerlist}}) {
          $j->{valid}=$states[$res]->[$rot];
          $validrotnum++ if ($j->{valid}==1);
          $rot++;
        }
        # now update valid rotamer number
        $i->{numvalidrot}=$validrotnum;
        # next residue (whether it is in the same site or not)
        $res++;
      }
    }
    if (defined($states[$res])) {
      die "Error in readRotState: $fname contains rotamers for more residues
      (". scalar(@states) . ") then there are design residues according to the ROTAMER
      object created with $FS_DEF->{rotinpf} " . ($res) . ".\n";
    }
}



=head2

 Title   : readDesignRotState
 Usage   : $rot->readDesignRotState($filename)
 Function: Reads in the rotamer state file produced by writeDesignRotState
           and marks each rotamer in the object as either valid or invalid.
 Returns : Nothing
 Args    : File name to read from.

=cut

sub readDesignRotState {
    my $self=shift;
    my $fname=shift;
    GENERAL::requireArgs($self, $fname);

    my $rotstate = GENERAL::GetInFH($fname);
    my @states=();

    # collect rotamer state info
    my $s = 0; my $new = 1;
    while (<$rotstate>) {
      $_ = GENERAL::Trim($_);
      if ($_ ne "") {
        my @lines=split(" ", $_);
        push (@{$states[$s]}, @lines);
        $new = 0;
      } elsif ($new == 0) {
        $s++;
        $new = 1;
      }
    }
    close($rotstate);

    $s = 0;
    # loop over the designed sites
    foreach my $site (@{$self->{dslist}}) {
      if (!defined($states[$s])) {
        GENERAL::error("There seem to be only " . ($s) . " sites in $fname but site ".
                       "$site->{chain}$site->{iresnum} is the " . ($s+1) . "-th design site according ".
                       "to the ROTAMER object (created with $FS_DEF->{rotinpf})!");
      }
      # loop over residue lists
      my $rot = 0;
      foreach my $i (@{$site->{reslist}}) {
        my $validrotnum=0;
        # loop over the rotamer states
        foreach my $j (@{$i->{rotamerlist}}) {
          if (!defined($states[$s]->[$rot])) {
            GENERAL::error("There seem to be only $rot rotamers in site number " .
            ($s+1) . " ($site->{chain}$site->{iresnum}) according to $fname, while according to " .
            "the ROTAMER object there are more rotamers in this site!");
          }
          $j->{valid} = $states[$s]->[$rot];
          $validrotnum++ if ($j->{valid}==1);
          $rot++;
        }
        # now update valid rotamer number
        $i->{numvalidrot} = $validrotnum;
        # next residue (whether it is in the same site or not)
      }
      $s++;
    }
    if (defined($states[$s])) {
      GENERAL::error("$fname contains rotamers for more sites (". scalar(@states) . ") then there are design ".
                     "sites according to the ROTAMER object created with $FS_DEF->{rotinpf} ($s)!");
    }
}

=head2

 Title   :
 Usage   :
 Function:
 Returns :
 Args    :

=cut


sub writeDesignInput {
   my $self = shift;
   $self->writeDesignRotState($FS_DEF->{rot_statesf});
   $self->writeDesignRotInput1($FS_DEF->{rotinp1f});
   $self->writeDesignRotInput2($FS_DEF->{rotinp2f});
   $self->writeDesignRotInput3($FS_DEF->{rotinp3f});
   $self->writeDesignAAInput($FS_DEF->{aainputf});
}


=head2 

 Title   : numify
 Usage   : $rotamer->numify()
 Function: Populates hash fields in the provided ROTAMER
           object, which have to do with the numbers of rotamers in
           each residue and site, number of valid rotamers in each,
           and total number of rotamers and valid rotamers.
 Returns : Note, all indices start at 1.
 Args    : 

=cut


sub numify {
  my $self = shift;

  # loop over the designed sites
  my $tnr = 0; # total number of rotamers
  my $nvr = 0; # number of valid rotamers
  my $si = 1; # site number
  foreach my $site (@{$self->{dslist}}) {
    my $stnr = 0; # total number of rotamers in this site
    my $snvr = 0; # number of valid rotamers in this site
    $site->{index} = $si;
    $si++;
    # loop over residue lists
    foreach my $i (@{$site->{reslist}}) {
      my $rtnr = 0; # total number of rotamers in this residue
      my $rnvr = 0; # number of valid rotamers in this residue
      # loop over the rotamer states
      foreach my $j (@{$i->{rotamerlist}}) {
        if ($j->{valid}) {
          $nvr++;
          $snvr++;
          $rnvr++;
          $j->{vrrindx} = $rnvr; # valid rotamer index of this rotamer in this residue
          $j->{vsrindx} = $snvr; # valid rotamer index of this rotamer in this site
          $j->{vtrindx} = $nvr;  # valid rotamer index of this rotamer in this ROTAMER structure
        }
        $tnr++;
        $stnr++;
        $rtnr++;
        $j->{rrindx} = $rtnr; # rotamer index of this rotamer in this residue
        $j->{srindx} = $stnr; # rotamer index of this rotamer in this site
        $j->{trindx} = $tnr; # rotamer index of this rotamer in this ROTAMER structure
      }
      $i->{numvalidrot} = $rnvr;
      $i->{numrot} = $rtnr;
    }
    $site->{numvalidrot} = $snvr;
    $site->{numrot} = $stnr;
  }
  $self->{numvalidrot} = $nvr;
  $self->{numrot} = $tnr;

  # Assign to each site, the number of valid rotamers with which rotamers from this
  # site have interaction energies (in other words, number of valid rotamers, which
  # come "after" the rotamers in this site).
  foreach my $site (@{$self->{dslist}}) {
    $nvr -= $site->{numvalidrot};
    $site->{num_valid_after} = $nvr;
  }
}


=head2

 Title   : getInterIndex
 Usage   : $rotamer->ROTAMER::getInterIndex($rot1, $rot2)
 Function: Returns the interaction index in the pair energy table
           of the given rotamers. So first rotamer of first design
           site and first rotamer of second design site have an
           interactio index of 1.
 Returns : Interaction index of the rotamer pare
 Args    : 1. rotamer 1
           2. rotamer 2
           3. Optional flag: if set, the pair index table
 Notes   : Running time is O(1) w/o calculating the pair start index
           table and O(ns^2) with calculating the pair start index
           table, where ns is the number of sites.

=cut


sub getInterIndex {
  my $self = shift;
  my $roti = shift;
  my $rotj = shift;
  GENERAL::requireArgs($roti, $rotj);
  my $calc = shift;
  if (!defined($calc)) { $calc = 0; }
  
  if ($calc || (!defined($self->{pair_start_index}))) {
    $self->numify();
    # Now, construct a hash table of pair indices. $Pi{i}{j} will contain the
    # index of the interaction of the 1st rotamer at site i with the first
    # rotamer at site j. This is a bit of a memory waste, but gives lots of speedup.
    $self->{pair_start_index} = {};
    my $ind = 0;
    for (my $si = 0; $si < scalar(@{$self->{dslist}}); $si++) {
      my $iind = $ind;
      for (my $sj = $si+1; $sj < scalar(@{$self->{dslist}}); $sj++) {
        $self->{pair_start_index}->{$si}{$sj} = $iind;
        $iind += $self->{dslist}->[$sj]->{numvalidrot};
      }
      $ind += ($self->{dslist}->[$si]->{numvalidrot})*($self->{dslist}->[$si]->{num_valid_after});
    }
  }

  my $si = $roti->{residue}->{dsite}->{index} - 1;
  my $sj = $rotj->{residue}->{dsite}->{index} - 1;
  my $ind = $self->{pair_start_index}->{$si}{$sj} + $rotj->{vsrindx};
  if ($roti->{vsrindx} > 1) { $ind +=  $self->{dslist}->[$si]->{num_valid_after}; }
  if ($roti->{vsrindx} > 2) { $ind += ($roti->{vsrindx} - 2)*($self->{dslist}->[$si]->{num_valid_after}); }
  return $ind;
}


=head2

 Title   : setValidRots
 Usage   : $rotamer->setValidRots(@valid_arr)
 Function: sets the valid flag of all the rotamers based on the binary
           array valid_arr.
 Returns : nothing
 Args    :

=cut
sub setValidRots {
    my $self = shift;
    my $valid_arr = shift;
    GENERAL::requireArgs($valid_arr);

    # loop over the designed sites
    my $c=0;
    foreach my $site (@{$self->{dslist}}) {
      # loop over residue lists
      foreach my $i (@{$site->{reslist}}) {
         my $validrotnum=0;
         # loop over the rotamer states
         foreach my $j (@{$i->{rotamerlist}}) {
           $j->{valid}=$valid_arr->[$c];
           if ($valid_arr->[$c] == 1) {
             $validrotnum++;
           }
           if (($valid_arr->[$c] != 0) && ($valid_arr->[$c] != 1)) {
             GENERAL::error("valid_arr has to be a binary array!");
           }
           $c++;
         }
         $i->{numvalidrot}=$validrotnum;
      }
    }
}


=head2

 Title   : assignValidRotamers
 Usage   : $rotamer->assignValidRotamers(\@selfE, $selfecut)
 Function: Sets each rotamer as either valid or not valid based on the
           received self energy array and the energy cutoff value.
 Returns : nothing
 Args    : 1. Pointer to an array of self energies (has to have as many
              entries as there are rotamers).
           2. Self energy cutoff value.
=cut

sub assignValidRotamers {
  my $self = shift;
  my $selfE = shift;
  my $selfecut = shift;

  # Make sure we have as many self energies as rotamers
  if (scalar(@$selfE) != $self->countSelfTotal()) {
    die "Error in assignValidRotamers: total number of self " .
    "terms is " . scalar(@$selfE) . " while total number of " .
    "rotamers is " . $self->countSelfTotal() . "!\n";
  }

  # loop over the designed sites
  my $c=0; my $fl = 0;
  foreach my $site (@{$self->{dslist}}) {
    # Number of valid rotamers at this site
    my $snvr = 0;
    # loop over residue lists
    foreach my $i (@{$site->{reslist}}) {
      my $validrotnum = 0;
      # loop over the rotamers
      foreach my $j (@{$i->{rotamerlist}}) {
        $j->{valid} = ($selfE->[$c] < $selfecut) ? 1:0;
        if ($j->{valid}) {
          $validrotnum++;
        }
        $c++;
      }
      $i->{numvalidrot} = $validrotnum;
      if ($i->{numvalidrot} == 0) {
        GENERAL::warning("All rotamers of residue $i->{resname} at site $site->{chain}$site->{iresnum} have been eliminated!");
        $fl = 1;
      }
      $snvr += $validrotnum;
    }
    if ($snvr == 0) {
      GENERAL::warning("Fatal problem: all rotamers for site $site->{chain}$site->{iresnum} have been eliminated!");
      $fl = 1;
    }
  }
  if ($fl) {
    GENERAL::error("FIT_PROBLEMS: Some residues were completely eliminated... revise rotamer input file or change elimination cutoff");
  }
}


=head2

 Title   : indexRotamers
 Usage   : $rotamer->indexRotamers($marker)
 Function: Gives each rotamer a new hash member named $marker and assigns
           to it the serial number of the rotamer. That is, all rotamers
           in the object are numbered 1 through N.
 Returns : nothing
 Args    : 1. Marker field name.
           2. Optional: valid flag. If set, will only consider valid rotamers.
              By default considers all.

=cut

sub indexRotamers {
  my $self = shift;
  my $marker = shift;
  GENERAL::requireArgs($self, $marker);
  my $valid = shift;
  if (!defined($valid)) { $valid = 0; }

  # loop over the designed sites
  my $c = 0;
  foreach my $site (@{$self->{dslist}}) {
    # loop over residue lists
    foreach my $res (@{$site->{reslist}}) {
       # loop over the rotamers
       foreach my $rot (@{$res->{rotamerlist}}) {
         next if ($valid && ($rot->{valid} != 1));
         $c++;
         $rot->{$marker} = $c;
       }
    }
  }
}


=head2

 Title   : markRotamers
 Usage   : $rotamer->markRotamers($beg_index, $end_index, $marker)
 Function: Marks the given range or rotamers with the given "marker". This
           means that it creates a field named $marker in each site, residue
           and rotamer covered by the given range.
 Returns : nothing
 Args    : 1. Begining index (starting with 1).
           2. End index.
           3. Marker field name.

=cut
sub markRotamers {
  my $self = shift;
  my $beg = shift;
  my $end = shift;
  my $marker = shift;
  GENERAL::requireArgs($self, $beg, $end, $marker);

  # loop over the designed sites
  my $c = 0;
  foreach my $site (@{$self->{dslist}}) {
    $site->{$marker} = 0;
    # loop over residue lists
    foreach my $res (@{$site->{reslist}}) {
       $res->{$marker} = 0;
       # loop over the rotamers
       foreach my $rot (@{$res->{rotamerlist}}) {
         $c++;
         if (($c >= $beg) && ($c <= $end)) {
           $site->{$marker} = 1;
           $res->{$marker} = 1;
           $rot->{$marker} = 1;
         } else {
           $rot->{$marker} = 0;
         }
       }
    }
  }
}


=head2

 Title   : markValidRotamerPairs
 Usage   : $rotamer->markValidRotamerPairs($beg_index, $end_index, $marker1, $marker2, $marker3, $marker4)
 Function: Given a range, marks the valid rotamer pairs corresponding to that
           range. The function assumes a "normal" nested loop treversal of valid
           rotamer pairs. In order to mark a range of pairs, the function marks the
           first site,residue,rotamer-combination in the outer loop (by setting a new
           hash member named $marker1 to the value 1 for the site, the residue, and
           the rotamer); the first site,residue,rotamer-combination in the inner loop
           (using $marker2); the last site,residue,rotamer-combination in the outer
           loop (using $marker3); and the last site,residue,rotamer-combination
           in the inner loop (using $marker4). In all other sites, residues, and rotamers
           the values corresponding to all the 4 markers is set to 0.
           This makes it easy to tell an outside function exactly
           which pairs of valid rotamers it is to visit. This function can be used
           to break up pairs of rotamers between parallel processes.
 Returns : nothing
 Args    : 1. Begining index (starting with 1).
           2. End index.
           3. Marker1 field name.
           4. Marker2 field name.
           5. Marker3 field name.
           6. Marker4 field name.

=cut
sub markValidRotamerPairs {
  my $self = shift;
  my $beg = shift;
  my $end = shift;
  my $marker1 = shift;
  my $marker2 = shift;
  my $marker3 = shift;
  my $marker4 = shift;
  GENERAL::requireArgs($self, $beg, $end, $marker1, $marker2, $marker3, $marker4);

  # set all visits to zero
  foreach my $site (@{$self->{dslist}}) {
    $site->{$marker1} = 0;
    $site->{$marker2} = 0;
    $site->{$marker3} = 0;
    $site->{$marker4} = 0;
    foreach my $res (@{$site->{reslist}}) {
       $res->{$marker1} = 0;
       $res->{$marker2} = 0;
       $res->{$marker3} = 0;
       $res->{$marker4} = 0;
       foreach my $rot (@{$res->{rotamerlist}}) {
         $rot->{$marker1} = 0;
         $rot->{$marker2} = 0;
         $rot->{$marker3} = 0;
         $rot->{$marker4} = 0;
       }
    }
  }
  
  # loop over the designed sites
  my $c = 0;
  for (my $i1 = 0; $i1 < scalar(@{$self->{dslist}}) - 1; $i1++) {
    my $site1 = $self->{dslist}->[$i1];
    for (my $j1 = 0; $j1 < scalar(@{$site1->{reslist}}); $j1++) {
      my $res1 = $site1->{reslist}->[$j1];
      for (my $k1 = 0; $k1 < scalar(@{$res1->{rotamerlist}}); $k1++) {
        my $rot1 = $res1->{rotamerlist}->[$k1];
        next if ($rot1->{valid} == 0);
        for (my $i2 = $i1+1; $i2 < scalar(@{$self->{dslist}}); $i2++) {
          my $site2 = $self->{dslist}->[$i2];
          my $fre = 0;
          for (my $j2 = 0; $j2 < scalar(@{$site2->{reslist}}); $j2++) {
            my $res2 = $site2->{reslist}->[$j2];
            my $fro = 0;
            for (my $k2 = 0; $k2 < scalar(@{$res2->{rotamerlist}}); $k2++) {
              my $rot2 = $res2->{rotamerlist}->[$k2];
              next if ($rot2->{valid} == 0);
              $c++;
              if ($c == $beg) {
                $site1->{$marker1} = 1;
                $res1->{$marker1} = 1;
                $rot1->{$marker1} = 1;
                $site2->{$marker2} = 1;
                $res2->{$marker2} = 1;
                $rot2->{$marker2} = 1;
              }
              if ($c == $end) {
                $site1->{$marker3} = 1;
                $res1->{$marker3} = 1;
                $rot1->{$marker3} = 1;
                $site2->{$marker4} = 1;
                $res2->{$marker4} = 1;
                $rot2->{$marker4} = 1;
              }
            }
          }
        }
      }
    }
  }
}

=head2

 Title   : unmarkRotamers
 Usage   : $rotamer->markRotamers($marker)
 Function: Unmarks the rotamers from the given marker.
 Returns : nothing
 Args    : 1. Marker field name.

=cut
sub unmarkRotamers {
    my $self = shift;
    my @markers = @_;
    GENERAL::requireArgs($self, @markers);

    foreach my $marker (@markers) {
      # loop over the designed sites
      foreach my $site (@{$self->{dslist}}) {
        if (defined($site->{$marker})) {
          undef($site->{$marker});
        }
        # loop over residue lists
        foreach my $res (@{$site->{reslist}}) {
          if (defined($res->{$marker})) {
            undef($res->{$marker});
          }
          # loop over the rotamers
          foreach my $rot (@{$res->{rotamerlist}}) {
            if (defined($rot->{$marker})) {
              undef($rot->{$marker});
            }
          }
        }
      }
    }
}


=head2

 Title   : writeRotamers
 Usage   : $rotamer->writeRotamers('test.crd', 'dummy pdb');
 Function: Creates and runs a CHARMM script to write out the content of the rotamer
           object in either the PDB or dummy format.
 Returns : nothing
 Args    : 1. The name of the initial CRD file.
           2. String specifying the format. The allowed values
           are 'dummy', 'pdb' or a combination of these (key insensitive). In
           the case of 'dummy', all valid rotamers will be output in dummy format
           (with the first line indicating the number of rotamers and number of
           atoms per rotamer). One file will be created for each residue. In the
           case of 'eef', PDB format will be used. Again, one file per residue will
           be created, which will contain a concatination of the PDB structures of
           all the valid rotamers for that residue. The rotamers will have appropriate
           residue id's in the PDB file to reflect their rotamer numbers.
=cut

sub writeRotamers {
  my $self = shift;
  my $crd = shift;
  my $format = uc(shift);
  my $im = shift;
  my $pdb = 0; my $dummy = 0;

  if ($format =~ /PDB/) { $pdb = 1;}
  if ($format =~ /DUMMY/) { $dummy = 1; }
  if (($format !~ /PDB/) && ($format !~ /DUMMY/)) {
    die "Error in writeRotamers: unknown format $format\n";
  }

  # Divide the task among processes
  my $Nt = $self->countResTotal();
  my $a = GENERAL::splitTasks($Nt, $NPROC_DEF);

  my $chinp;
  if (defined($LOCAL_DEF)) {
    $chinp = "$LOCAL_DEF/write_rots.$PRANK_DEF.inp";
  } else {
    $chinp = "write_rots.$PRANK_DEF.inp";
  }
  my $charmm = CHARMM::new($CHARMM_DEF, $chinp);
  $charmm->loadParm('scaled');
  $charmm->setupFromCRD($crd, "none");
  $charmm->send("hbuild\n");
#  $charmm->defineTemplate($self);
  my $backbone = CHARMM::backbone();

  my $chidef = $self->{chidef};
  my $i = -1; my $num = 0;
  foreach my $site (@{$self->{dslist}}) {

    # for each designed site, get chain id and residue number
    $i++;
    my $chain = $site->{chain};
    my $iresnum = $site->{iresnum};

    my $j=-1; # will keep track of residue number
    foreach my $r (@{$site->{reslist}}) {
      $j++;
      $num++;
      next if (($num < $a->[$PRANK_DEF]->{beg}) || ($num > $a->[$PRANK_DEF]->{end}));
#printf("Process $PRANK_DEF to do $i/$j/$im (out of $Nt)\n");
      my $resname = uc($r->{resname});
      my $rotout = (defined($im)) ? $FS_DEF->{rot_ioutf} : $FS_DEF->{rot_outf};
      $rotout =~ s/%%%/$im/ if (defined($im));
      $rotout =~ s/%%/$i/;
      $rotout =~ s/%/$j/;
      $rotout = File::Spec->abs2rel($rotout);
      my $rotpdb = (defined($im)) ? $FS_DEF->{rot_ipdbf} : $FS_DEF->{rot_pdbf};
      $rotpdb =~ s/%%%/$im/ if (defined($im));
      $rotpdb =~ s/%%/$i/;
      $rotpdb =~ s/%/$j/;
      $rotpdb = File::Spec->abs2rel($rotpdb);
      my $rotoutunit = 14;
      my $rotpdbunit = 15;

      if ($dummy) {
        $charmm->openCard($rotout, $rotoutunit);
      }
      if ($pdb) {
        $charmm->openCard($rotpdb, $rotpdbunit);
      }

      my $k = 1; # will keep track of rotamer number
      my $first = 1; # for the first rotamer, need to print number of atoms and number of rotamers
      
      # Delete design rotamer
#      $charmm->verbose("define rots select .byres. atom $site->{chain} $site->{iresnum} * .and. .not. ($backbone) end\n");
#      $charmm->verbose("delete atom select rots end\n");
      # Patch the new residue
      $charmm->patchRes($resname, $chain, $iresnum);
      foreach my $rot ( @{$r->{rotamerlist}} ) {
        # If we need this rotamer
        if ($rot->{valid}) {
          # Put in the next rotamer
          $charmm->placeRotamer($chain, $resname, $iresnum, $chidef, $rot);

          # For dummy, write out all/only valid (depending on the valid flag)
          # rotamers for a residue into one file in dummy format
          if ($dummy) {
            # Redefine rotamer since patching has cleared the definitions
            $charmm->verbose("define rots select .byres. atom $site->{chain} $site->{iresnum} * .and. .not. ($backbone) end\n");
            if ($first) {
              # need to know how many atoms per rotamer and how many rotamers
              $charmm->writeCard($rotoutunit, "\* ?NSEL $r->{numvalidrot}");
              $first = 0;
            }
            $charmm->writeMulteCRD($rotoutunit, 'sele rots end');
          }
          if ($pdb) {
#            # Rename the residue name (patch does not do so) and renumber to reflect rotamer number
#            my $resid = "$k";
#            $charmm->verbose("rename RESID $resid select RESID $site->{iresnum} .and. SEGID $site->{chain} end\n");
#            $charmm->verbose("write coor PDB select .byres. atom $site->{chain} $resid * end unit $rotpdbunit\n");
#            # Renumber back not to mess up the order of residues
#            $charmm->verbose("rename RESID $site->{iresnum} select RESID $resid .and. SEGID $site->{chain} end\n");
            
            # Renumber the residue to the rotamer number. Since the residue needs to be renumbered back, in order 
            # not to get confused with RESID's, rename the residue to XXX - later to be substituted by the residue
            # name in the output file.
#            my $resid = "$k";
#            $charmm->verbose("rename RESName XXX select RESID $site->{iresnum} .and. SEGID $site->{chain} end\n");
#            $charmm->verbose("rename RESId $resid select .byres. RESName XXX end\n");
#            $charmm->verbose("write coor PDB select .byres. RESName XXX end unit $rotpdbunit\n");
#            $charmm->verbose("rename RESId $site->{iresnum} select .byres. RESName XXX end\n");
#            $charmm->verbose("rename RESName $resname select .byres. RESName XXX end\n");
          
            $charmm->verbose("write coor PDB select RESID $site->{iresnum} .and. SEGID $site->{chain} end unit $rotpdbunit\n");
          }
        }
        $k++;
      }
      if ($dummy) {
        $charmm->closeCard($rotoutunit);
      }
      if ($pdb) {
        $charmm->closeCard($rotpdbunit);
      }
      # If we just patched a proline, the amide hydrogen was removed, so re-initialize the structure 
      if ($resname =~ /PRO/) {
        $charmm->verbose("delete atom select all end");
        $charmm->setupFromCRD($crd, "none");
      }
    }
  }
  $charmm->endofstory();

  # Remove the END's in between rotamers in the output PDB files
  if ($pdb) {
    $num = 0;
    for (my $i=0; $i < scalar(@{$self->{dslist}}); $i++) {
      for (my $j=0; $j < scalar(@{$self->{dslist}->[$i]->{reslist}}); $j++) {
        $num++;
        next if (($num < $a->[$PRANK_DEF]->{beg}) || ($num > $a->[$PRANK_DEF]->{end}));
        my $rotpdb = (defined($im)) ? $FS_DEF->{rot_ipdbf} : $FS_DEF->{rot_pdbf};
        $rotpdb =~ s/%%%/$im/ if (defined($im));
        $rotpdb =~ s/%%/$i/;
        $rotpdb =~ s/%/$j/;
        # 1. remove ENDs
        # 2. renumber residues to correspond to rotamer numbers
        # 3. fix CHARMM's improper alignment of 4-digit residue numbers
        my $arr = GENERAL::file2array($rotpdb);
        my $ofh = GENERAL::GetOutFH($rotpdb);
        my $k = 1; # rotamer number
        foreach (@$arr) {
          next if (/TER/);
          if (/END/) { $k++; }
          elsif (/^\s*$/) {}
          elsif (/ATOM/) {
            $ofh->printf("%s%4d %s\n", substr($_, 0, 22), $k, substr($_, 27));
          } else { $ofh->printf($_."\n"); }
        }
        close($ofh);
#        # remove the END's from the file
#        GENERAL::csystem("echo \"1,\\\$s/^END//g\nw\nq\n\" | ed $rotpdb > /dev/null 2>&1");
#        # rename the residue properly
#        GENERAL::csystem("echo \"1,\\\$s/XXX/$self->{dslist}->[$i]->{reslist}->[$j]->{resname}/g\nw\nq\n\" | ed $rotpdb > /dev/null 2>&1");
      }
    }
  }
  # Do we need to delete the script?
  if ($FS_DEF->{del_charmm_inp}) {
    GENERAL::csystem("rm -f $chinp");
  }
}



=head2

 Title   : toPDB
 Usage   : my $rpdb = $rotamer->toPDB();
 Function: Converst the rotamer object into a PDB structure. There is one chain
           in the PDB structure for every residue at every site site 
           (chain id = "$chain$iresnum$resname" like A1ALA or B12GLU). 
           Residue id of each residue equals to the rotamer number
           of the corresponding rotamer.
 Returns : PDB structure
 Args    : 1. Rotamer object.
=cut

sub toPDB {
  my $self = shift;
  GENERAL::requireArgs($self);

  my $rpdb = PDB::new();
  my $i = -1;  # will keep track of site number
  foreach my $site (@{$self->{dslist}}) {
    $i++;
    my $chain = $site->{chain};
    my $iresnum = $site->{iresnum};
    my $j=-1; # will keep track of residue number
    foreach my $res (@{$site->{reslist}}) {
      $j++;
      my $resname = uc($res->{resname});
      my $rotpdb = $FS_DEF->{rot_pdbf};
      $rotpdb =~ s/%%/$i/;
      $rotpdb =~ s/%/$j/;
      
      if (! -r $rotpdb) { GENERAL::error("Could not convert ROTAMER to PDB because file $rotpdb is not readable/existent!"); }
      my $pdb = PDB::new($rotpdb);
      if (scalar(@{$pdb->{chain}}) != 1) { GENERAL::error(sprintf("File $rotpdb contains %d chains, not 1", scalar(@{$pdb->{chain}}))); }
      $pdb->{chain}->[0]->{id} = "$site->{chain}$site->{iresnum}$resname";
      $rpdb->insertChain($pdb->{chain}->[0]);
    }
  }
  
  return $rpdb;
}


=head2

 Title   : findNeighbors
 Usage   : $rotamer->findNeighbors();
 Function: For each residue in the ROTAMER object, finds all of the other residues
           in the same object within the specified distance cutoff from it. The list
           of these neighboring residues is stored under the 'int_res' key in each
           residue (further hashed by site id and further by aa name). This is
           convenient for limiting pairwise calculations. The distance considered is
           between the surfaces of vdW radii.
 Returns : nothing
 Args    : 0. Rotamer object.
           1. Distance cutoff.
=cut

sub findNeighbRes {
  my $rotamer = shift;
  my $dist = shift;

  # Convert the ROTAMER object into a PDB
  my $rpdb = $rotamer->toPDB();
  my $R = PDB::sizeLookup();
  # For each residue, find all of its neighbors
  my @rres = $rpdb->ConRes();
  my $i = -1;
  foreach my $site (@{$rotamer->{dslist}}) {
    foreach my $res (@{$site->{reslist}}) {
      my @sa;
      foreach my $rot (@{$res->{rotamerlist}}) {
        $i++;
        push(@sa, @{$rres[$i]->{atom}});
      }
      # Clean any existing list of neighbors
      $res->{int_res} = undef;
      my $wresa = $rpdb->resWithin(\@sa, $dist, $R);
      foreach my $wres (@$wresa) {
        if (("$site->{chain}$site->{iresnum}" ne $wres->{chain}->{id}) || ($res->{resname} ne $wres->{resname})) {
          $res->{int_res}->{$wres->{chain}->{id}}->{$wres->{resname}} = 1;
        }
      }
    }
  }
  
  # Explicitly destroy the PDB object (because of circular referencing, it won't be garvage collected)
#  $rpdb->DESTROY();
}


=head2

 Title   : extent
 Usage   : $rotamer->extent();
 Function: For each residue in the ROTAMER object, calculates the longest distance
           from its CB atom (CA for GLY) to the surface of any other atom in the side chain (over all
           rotamers). Stores this number under the hash keys "extent". Also stores the coordinates of
           the CB atom (CA for GLY) as an array hashed by the key "base_coor".
 Returns : nothing
 Args    : 0. Rotamer object.

=cut

sub extent {
  my $rotamer = shift;

  # Convert the ROTAMER object into a PDB
  my $rpdb = $rotamer->toPDB();
  my $R = PDB::sizeLookup();
  # For each residue, find its extent
  my @rres = $rpdb->ConRes();
  my $i = -1;
  foreach my $site (@{$rotamer->{dslist}}) {
    foreach my $res (@{$site->{reslist}}) {
      my @sa;
      foreach my $rot (@{$res->{rotamerlist}}) {
        $i++;
        push(@sa, PDB::sidechain($rres[$i]));
      }
      my $ba; # base atom (either CA for GLY or CB)
      my $f = 0;
      if ($res->{resname} eq "GLY") {
        foreach my $a (@{$rres[$i]->{atom}}) {
          if ($a->{atomname} eq "CA") {
            $ba = $a;
            $f = 1; last;
          }
        }
      } else {
        foreach my $a (@{$rres[$i]->{atom}}) {
          if ($a->{atomname} eq "CB") {
            $ba = $a;
            $f = 1; last;
          }
        }
      }
      if ($f == 0) { GENERAL::error("Failed to find the base atom for residue $res->{resname} at site $site->{chain}$site->{iresnum}."); }
      $res->{extent} = undef;
      foreach my $a (@sa) {
        my $ext = PDB::atomDist($a, $ba) + $R->{$a->{residue}->{resname}}{$a->{atomname}};
        if (!defined($res->{extent})) {
          $res->{extent} = $ext;
          next;
        }
        if ($res->{extent} < $ext) {
          $res->{extent} = $ext;
        }
      }
      if (!defined($res->{extent})) { 
        $res->{extent} = $R->{$ba->{residue}->{resname}}{$ba->{atomname}};
      }
      my @coor = ($ba->{xcoor}, $ba->{ycoor}, $ba->{zcoor});
      $res->{base_coor} = \@coor;
    }
  }
  
  # Explicitly destroy the PDB object (because of circular referencing, it won't be garbage collected)
#  $rpdb->DESTROY();
}


=head2

 Title   : rewriteRotmulte
 Usage   :
 Function:
 Returns :
 Args    :

=cut

sub rewriteRotmulte {

    my $self=shift;
    my $sitelist= $self->{dslist};

    GENERAL::cchdir("rotmulte");
    my $i=-1;
    foreach my $site (@{$sitelist})   {
      $i++;
      my $j=-1;
      foreach my $res ( @{$site->{reslist}})  {
        $j++;
        my $rotout='rotmulte'.$i.'-'.$j.'.out';
        my $rotdat='rotmulte'.$i.'-'.$j.'.dat';
        my $rotoutfh=&GENERAL::GetInFH($rotout);
        my $rotdatfh=&GENERAL::GetOutFH($rotdat);

        my $line=0;
        my $start=999;
        my $end=999;
        my ($atomno, $rotnum);
        my $k=-1;
        my @coor;

        while (<$rotoutfh> ) {
          $line++;
          if ($line==1) {
            my @lines=split(" ", $_);
            $atomno=$lines[0];
            $rotnum=$lines[1];
          } else {
            push (@coor, $_);
          }
        }

        $rotdatfh->printf("%2d%3d\n", $atomno, $res->{numvalidrot});
        foreach my $rot ( @{$res->{rotamerlist}}) {
          $k++;
          if ($rot->{valid} == 1) {
            $start=$k*$atomno;
            $end=($k+1)*$atomno-1;
            for (my $m=$start; $m<=$end; $m++) {
              $rotdatfh->print($coor[$m]);
            }
          }
        }
        close($rotoutfh);
        close($rotdatfh);
      }
    }
    GENERAL::cchdir("../");
}




=head2 

 Title   : readSigma
 Usage   : 
 Function: Reads the deviations from each chi angle in each rotamer of each residue.
           This is for use in creating a subrotamer library.
 Returns : 
 Args    : 

=cut


sub readSigma {

    my $self=shift;
 
    $self->{sigma}=undef if (defined $self->{sigma});
    my $lib = "$ROTLIB_DEF/sigma_0";

    my $libfh = &GENERAL::GetInFH($lib);
    my $sigmaref={};
    my $resname;

    while (<$libfh>) {

      # skip the empty line
      next if (/^$/);
      if (/^RESI/ ) {

                my @line=split(" ", $_);
          # here we parse the line started with RESI
          $resname= $line[1];
          $sigmaref->{$resname}=[];

      }

      if ( ! /^RESI/) {
          my @line = split (" ", $_);
          push (@{$sigmaref->{$resname}}, \@line);
      }
    }

   $self->{sigma}=$sigmaref;
}



=head2 

 Title   :  writeSubRotlib()
 Usage   : 
 Function:  subrotmaer library is created on the fly. 
            
 Returns : 
 Args    : (1) backbone index
           (2) method to expand the rotamer
=cut


sub writeSubRotlib {

    my $self=shift;
    my $bbindex=shift;
    my $method=shift;
    

    # get sigma
    $self->readSigma() if (!defined $self->{sigma});
 
    # get proper rotamer library reference
    my $libref=$self->readRotlib($bbindex);

    # give a name for subrotamer
    my $fname='sublibrary_'.$bbindex;
    if (!defined $method) { die "Error: subrotamer library generation method not specified!\n"; }
    #print " $fname  $method $bbindex  \n";


    if ($method eq 'discrete') {
	my $newlib=$self->_expandDiscreteRotamerlib($libref);
#        my $newlib=$libref;
	$self->_writeNewlib($newlib, $fname);
    } elsif 
     ($method eq 'gaussian') {
	my $newlib=$self->_expandGaussianRotamerlib($libref);
	$self->_writeNewlib($newlib,$fname);
    } elsif ($method eq 'random') {
	my $newlib=$self->_expandRandomRotamerlib($libref);
	$self->_writeNewlib($newlib, $fname);
    } else {
	die "undefined method to expand subrotmaer method!\n";
    }

}




=head2 

 Title   :  _writeNewlib
 Usage   : 
 Function:  write sub rotamer rotamer library, for checking purpose
 Returns : 
 Args    : 

=cut

sub _writeNewlib {

    my $self=shift;
    my $newlib=shift;
    my $fname=shift;

    $newlib->{'HSD'} = $newlib->{'HIS'};

    #my $newlib =$fname;
    my $outf= &GENERAL::GetOutFH($fname);
    
   # $outf->print("RESI ALA  0\n");
    foreach my $key (sort keys %{$newlib}) {
	my $number2=  (defined $newlib->{$key}->{numchi})?$newlib->{$key}->{numchi}:" "; 
	my $number1=  $newlib->{$key}->{numrot};
	$outf->print("RESI $key  $number1 $number2\n");
        #if (defined $newlib->{$key}->{rotamerlist}->{chi}) {
	 foreach my $i (0..$#{$newlib->{$key}->{rotamerlist}}) {
	    foreach my $j (0..$#{$newlib->{$key}->{rotamerlist}->[$i]->{chi}}) {
        #        print "$key  $i   $j   $newlib->{$key}->{rotamerlist}->[$i]->{chi}->[$j]  \n"; 
		 printf $outf "%6.1f ", $newlib->{$key}->{rotamerlist}->[$i]->{chi}->[$j];
	    }
	    $outf->print("\n");
	 }
       # }
    }
    close $outf;
}



=head2 

 Title   : _expandDiscreteRotamerlib
 Usage   : 
 Function: expand rotamer library according to average deviation
           current implementation expand chi1 into 3, chi2 into 3.
           I may change it to expand into arbitrary number of subrotamer
           ( note that: subrotamer approach significantly increasing the computation
           time.) 
 Returns : 
 Args    : 

=cut


sub _expandDiscreteRotamerlib {

    my $self=shift;
    my $libref=shift;
    die "Error in _expandDiscreteRotamerlib: I don't know the deviation!\n" if (!defined $self->{sigma});
    die "Error in _expandDiscreteRotamerlib: you do not seem to have called useSubrotamer!\n" if ($self->{subrotamer}->{isValid} != 1);

    my $newlibref;
						   
    my $sigma= $self->{sigma};
    foreach my $resname (sort keys %{$libref}) {
        #print $resname, "\n";
	$newlibref->{$resname}={};
	$newlibref->{$resname}->{resname}=$resname;
	
	my $rotnum;
	# here we define the total avaible subrotamer
	if ( $self->{subrotamer}->{number}->{$resname} == 0) {
	    # like pro and ala,  they will keep same 
	    $rotnum=$libref->{$resname}->{numrot};
            #print $rotnum, "\n";
	} else {
	    $rotnum=$self->{subrotamer}->{number}->{$resname}*$libref->{$resname}->{numrot};
            #print $rotnum, "\n";
	}
	$newlibref->{$resname}->{numrot}=$rotnum;
	$newlibref->{$resname}->{numvalidrot}=$rotnum;
        $newlibref->{$resname}->{rotamerlist}=[];

	# total entry is reduced by 1 since we don't have probability term
#	$newlibref->{$resname}->{entry}= $libref->{$resname}->{entry}-1 if ( defined $libref->{$resname}->{entry}) ;


	my $ir=-1;
	for (my $i=0; $i<=$#{$libref->{$resname}->{rotamerlist}}; $i++) {
	    # expand rotamerlibrary, use subrotamer par hash to keep track the number of 
            # subrotamers
            #print "rotamer $i \n";

	    if ($self->{subrotamer}->{number}->{$resname} == 9 ) {
		my $var1 = $libref->{$resname}->{rotamerlist}->[$i]->{chi}->[0];
                my $var2 = $sigma->{$resname}->[$i]->[0];
                my $var3 = $libref->{$resname}->{rotamerlist}->[$i]->{chi}->[1];
                my $var4 = $sigma->{$resname}->[$i]->[1];
                my @newarray = @{$libref->{$resname}->{rotamerlist}->[$i]->{chi}};
                shift @newarray;
		shift @newarray;

		for (my $j=-1; $j<=1; $j++ ) {
		    for (my $k=-1; $k<=1; $k++) {
                        $ir++;
                        #print "sub rotamer index:  $ir\n";
                        $newlibref->{$resname}->{rotamerlist}->[$ir]={};
                        $newlibref->{$resname}->{rotamerlist}->[$ir]->{chi}=[];
                        $newlibref->{$resname}->{rotamerlist}->[$ir]->{valid}=undef;
			#print $libref->{$resname}->{rotamerlist}->[$i]->{chi}->[0], "\n";
			#my $var1 = $libref->{$resname}->{rotamerlist}->[$i]->{chi}->[0];
			#my $var2 = $sigma->{$resname}->[$i]->[0];
                        #print "resname $resname  rotamer lib  sigma $i  $var1  $var2 \n";
			my $chi1= $var1 + $j*$var2;
			$chi1 -=360.0 if ( $chi1 >=180.0 );
			$chi1 +=360.0 if ( $chi1 <=-180.0); 
                        #print "chi1 ", $chi1, "\n";
			#print $libref->{$resname}->{rotamerlist}->[$i]->{chi}->[1], "\n";
                        #my $var3 = $libref->{$resname}->{rotamerlist}->[$i]->{chi}->[1];
                        #my $var4 = $sigma->{$resname}->[$i]->[1];
			my $chi2= $var3 + $k* $var4;
			$chi2 -=360.0 if ( $chi2 >=180.0 );
			$chi2 +=360.0 if ( $chi2 <=-180.0); 
                        #print "chi2 ", $chi2, "\n";
			my @chi=($chi1, $chi2, @newarray);
			$newlibref->{$resname}->{rotamerlist}->[$ir]->{chi}=\@chi;
                        #print " my chi array:   @chi\n";
			my $valid=1;
			$newlibref->{$resname}->{rotamerlist}->[$ir]->{valid}=$valid;
			
		    }
		}
		
	    } elsif ($self->{subrotamer}->{number}->{$resname} == 3 ) {
                my $var1 = $libref->{$resname}->{rotamerlist}->[$i]->{chi}->[0];
                my $var2 = $sigma->{$resname}->[$i]->[0];
                my @newarray = @{$libref->{$resname}->{rotamerlist}->[$i]->{chi}};
                shift @newarray;
 
		for (my $j=-1; $j<=1; $j++ ) {

                    $ir++;
                    $newlibref->{$resname}->{rotamerlist}->[$ir]={};
                    $newlibref->{$resname}->{rotamerlist}->[$ir]->{chi}=[];
                    $newlibref->{$resname}->{rotamerlist}->[$ir]->{valid}=undef;

		    # expand chi1
		    my $chi1=$var1+ $j*$var2;
		    # in case chi1, chi2 out of -180.. 180 range 
		    $chi1 -=360.0 if ( $chi1 >=180.0 );
		    $chi1 +=360.0 if ( $chi1 <=-180.0); 
		    my @chi=($chi1,  @newarray);
		    $newlibref->{$resname}->{rotamerlist}->[$ir]->{chi}= \@chi;
		    my $valid=1;
		    $newlibref->{$resname}->{rotamerlist}->[$ir]->{valid}= $valid;
		}
	    } elsif  ($self->{subrotamer}->{number}->{$resname} == 0) {
		$newlibref->{$resname}=$libref->{$resname};
	    }
	}
    }
    return $newlibref;
}


=head2 

 Title   :  _expandGaussianRotamerlib
 Usage   :  
 Function:  generate gausian distributed rotamer library
            sample size is determined by $self->{subrotamer}->{par}
 Returns : 
 Args    : 

=cut


sub _expandGaussianRotamerlib {

    my $self=shift;
    my $libref=shift;
    die "Error in _expandGaussianRotamerlib: I don't know the deviation!\n" if (!defined $self->{sigma});
    die "Error in _expandGaussianRotamerlib: you do not seem to have called useSubrotamer!\n" if ($self->{subrotamer}->{isValid} != 1);

    my $newlibref;

    my $sigma= $self->{sigma};
    foreach my $resname (sort keys %{$libref}) {

        $newlibref->{$resname}={};
        $newlibref->{$resname}->{resname}=$resname;

	my $rotnum;
        # here we define the total avaible subrotamer
        if ( $self->{subrotamer}->{number}->{$resname} == 0) {
            # like pro and ala,  they will keep same
            $rotnum=$libref->{$resname}->{numrot};
        } else {
            $rotnum=$self->{subrotamer}->{number}->{$resname}*$libref->{$resname}->{numrot};
        }
        $newlibref->{$resname}->{numrot}=$rotnum;
        $newlibref->{$resname}->{numvalidrot}=$rotnum;

        # total entry is reduced by 1 since we don't have probability term
#        $newlibref->{$resname}->{entry}=$libref->{$resname}->{entry}-1;

        for (my $i=0; $i<=$#{$libref->{$resname}->{rotamerlist}->{chi}}; $i++) {
            # expand rotamerlibrary, use subrotamer par hash to keep track the number of
            # subrotamers

            if ($self->{subrotamer}->{number}->{$resname} == 9 ) {

                for (my $j=0; $j<=8; $j++ ) {
		    my $random1=&GENERAL::gaussianRand();
		    my $chi1=$libref->{$resname}->{rotamerlist}->{chi}->[$i]->[0]+ $random1*$sigma->{$resname}->[$i]->[0];
		    $chi1 -=360.0 if ( $chi1 >=180.0 );
		    $chi1 +=360.0 if ( $chi1 <=-180.0);
		    my $random2=&GENERAL::gaussianRand();
		    my $chi2=$libref->{$resname}->{rotamerlist}->{chi}->[$i]->[1]+ $random2*$sigma->{$resname}->[$i]->[1];
		    $chi2 -=360.0 if ( $chi2 >=180.0 );
		    $chi2 +=360.0 if ( $chi2 <=-180.0);
		    shift @{$libref->{$resname}->{rotamerlist}->{chi}->[$i]};
		    shift @{$libref->{$resname}->{rotamerlist}->{chi}->[$i]};
		    my @chi=($chi1, $chi2, @{$libref->{$resname}->{rotamerlist}->{chi}->[$i]});
		    push (@{$newlibref->{$resname}->{rotamerlist}->{chi}}, \@chi);
		    my $valid=1;
		    push  (@{$newlibref->{$resname}->{rotamerlist}->{valid}}, $valid);
		}
	    } elsif ($self->{subrotamer}->{number}->{$resname} == 3 ) {
                for (my $j=-1; $j<=1; $j++ ) {

                    my $random=&GENERAL::gaussianRand();
                    # expand chi1
                    my $chi1=$libref->{$resname}->{rotamerlist}->{chi}->[$i]->[0]+ $random*$sigma->{$resname}->[$i]->[0];

                    # in case chi1, chi2 out of -180.. 180 range
                    $chi1 -=360.0 if ( $chi1 >=180.0 );
                    $chi1 +=360.0 if ( $chi1 <=-180.0);
                    shift @{$libref->{$resname}->{rotamerlist}->{chi}->[$i]};

                    my @chi=($chi1,  @{$libref->{$resname}->{rotamerlist}->{chi}->[$i]});
                    push (@{$newlibref->{$resname}->{rotamerlist}->{chi}}, \@chi);
                    my $valid=1;
                    push  (@{$newlibref->{$resname}->{rotamerlist}->{valid}}, $valid);
                }
            } elsif  ($self->{subrotamer}->{number}->{$resname} == 0) {
                $newlibref->{$resname}=$libref->{$resname};
            }
	}
    }
    return $newlibref;
}



=head2 

 Title   :  _expandRandomRotamerlib
 Usage   : 
 Function:  generate random distributed rotamer library
            from -sigma to +sigmae
 Returns : 
 Args    : 

=cut


sub _expandRandomRotamerlib {


    my $self=shift;
    my $libref=shift;
    die "Error in _expandRandomRotamerlib: I don't know the deviation!\n" if (!defined $self->{sigma});
    die "Error in _expandRandomRotamerlib: you do not seem to have called useSubrotamer!\n" if ($self->{subrotamer}->{isValid} != 1);
    
    my $newlibref;

    my $sigma= $self->{sigma};
    foreach my $resname (sort keys %{$libref}) {
    
        $newlibref->{$resname}={};
        $newlibref->{$resname}->{resname}=$resname;

	my $rotnum;
        # here we define the total avaible subrotamer
        if ( $self->{subrotamer}->{number}->{$resname} == 0) {
            # like pro and ala,  they will keep same
            $rotnum=$libref->{$resname}->{numrot};
        } else {
            $rotnum=$self->{subrotamer}->{number}->{$resname}*$libref->{$resname}->{numrot};
        }
        $newlibref->{$resname}->{numrot}=$rotnum;
        $newlibref->{$resname}->{numvalidrot}=$rotnum;

        # total entry is reduced by 1 since we don't have probability term
#        $newlibref->{$resname}->{entry}=$libref->{$resname}->{entry}-1;

        for (my $i=0; $i<=$#{$libref->{$resname}->{rotamerlist}->{chi}}; $i++) {
            # expand rotamerlibrary, use subrotamer par hash to keep track the number of
            # subrotamers

            if ($self->{subrotamer}->{number}->{$resname} == 9 ) {

                for (my $j=0; $j<=8; $j++ ) {
                    my $random1=2*rand()-1.0;;
                    my $chi1=$libref->{$resname}->{rotamerlist}->{chi}->[$i]->[0]+ $random1*$sigma->{$resname}->[$i]->[0];
                    $chi1 -=360.0 if ( $chi1 >=180.0 );
                    $chi1 +=360.0 if ( $chi1 <=-180.0);
                    my $random2=2*rand()-1.0;;
                    my $chi2=$libref->{$resname}->{rotamerlist}->{chi}->[$i]->[1]+ $random2*$sigma->{$resname}->[$i]->[1];
                    $chi2 -=360.0 if ( $chi2 >=180.0 );
                    $chi2 +=360.0 if ( $chi2 <=-180.0);
                    shift @{$libref->{$resname}->{rotamerlist}->{chi}->[$i]};
                    shift @{$libref->{$resname}->{rotamerlist}->{chi}->[$i]};
                    my @chi=($chi1, $chi2, @{$libref->{$resname}->{rotamerlist}->{chi}->[$i]});
                    push (@{$newlibref->{$resname}->{rotamerlist}->{chi}}, \@chi);
                    my $valid=1;
                    push  (@{$newlibref->{$resname}->{rotamerlist}->{valid}}, $valid);
                }
            } elsif ($self->{subrotamer}->{number}->{$resname} == 3 ) {
                for (my $j=-1; $j<=1; $j++ ) {

                    my $random=2*rand()-1.0;;
                    # expand chi1
                    my $chi1=$libref->{$resname}->{rotamerlist}->{chi}->[$i]->[0]+ $random*$sigma->{$resname}->[$i]->[0];

                    # in case chi1, chi2 out of -180.. 180 range
                    $chi1 -=360.0 if ( $chi1 >=180.0 );
                    $chi1 +=360.0 if ( $chi1 <=-180.0);
                    shift @{$libref->{$resname}->{rotamerlist}->{chi}->[$i]};

                    my @chi=($chi1,  @{$libref->{$resname}->{rotamerlist}->{chi}->[$i]});
                    push (@{$newlibref->{$resname}->{rotamerlist}->{chi}}, \@chi);
                    my $valid=1;
                    push  (@{$newlibref->{$resname}->{rotamerlist}->{valid}}, $valid);
                }
            } elsif  ($self->{subrotamer}->{number}->{$resname} == 0) {
                $newlibref->{$resname}=$libref->{$resname};
            }
        }
    }
    return $newlibref;
}


=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub countValidPairTotal {
    my $self=shift;

    $self->numify();

    my $sitelist=$self->{dslist};
    my $counter=0;
    for (my $s1=0; $s1 <=$#{$sitelist}; $s1++) {
      foreach my $res1 (@{$sitelist->[$s1]->{reslist}}) {
        for (my $s2=$s1+1; $s2 <=$#{$sitelist}; $s2++) {
          foreach my $res2 (@{$sitelist->[$s2]->{reslist}}) {
            $counter += ($res1->{numvalidrot})*($res2->{numvalidrot});
          }
        }
      }
    }

   return $counter;
}


=head2 

 Title   : countPairTotal
 Usage   : $rotamer->coutPairTotal
 Function: counts the total number of rotamer pairs of rotamers in the ROTAMER object
 Returns : 
 Args    : 

=cut


sub countPairTotal {
    my $self=shift;
    my $sitelist=$self->{dslist};

    my $counter=0;
      for (my $s1=0; $s1 <=$#{$sitelist}; $s1++) {
        for (my $res1=0; $res1<=$#{$sitelist->[$s1]->{reslist}}; $res1++) {
            for (my $s2=$s1+1; $s2 <=$#{$sitelist}; $s2++) {
                for (my $res2=0; $res2<=$#{$sitelist->[$s2]->{reslist}}; $res2++) {
                    for(my $rot1=0; $rot1<=$#{$sitelist->[$s1]->{reslist}->[$res1]->{rotamerlist}}; $rot1++) {
                        for(my $rot2=0; $rot2<=$#{$sitelist->[$s2]->{reslist}->[$res2]->{rotamerlist}}; $rot2++) {
                           $counter++;
                        }
                    }
                }
            }
        }
    }

   return $counter;
}


=head2 

 Title   : countPairs
 Usage   : $rotameter->countPairs()
 Function: counts the numnber of pairs of residues in the ROTAMER object
 Returns : 
 Args    : 

=cut

sub countPairs {
  my $self=shift;
  my $sitelist=$self->{dslist};

  my $counter=0;
  for (my $s1=0; $s1 <=$#{$sitelist}; $s1++) {
    for (my $res1=0; $res1<=$#{$sitelist->[$s1]->{reslist}}; $res1++) {
      for (my $s2=$s1+1; $s2 <=$#{$sitelist}; $s2++) {
        for (my $res2=0; $res2<=$#{$sitelist->[$s2]->{reslist}}; $res2++) {
          $counter++;
        }
      }
    }
  }
  return $counter;
}

=head2

 Title   : countResTotal
 Usage   : my $num_res = $rotameter->countResTotal()
 Function: Counts the number of residues in all sites together.
 Returns :
 Args    :

=cut

sub countResTotal {
  my $self=shift;
  my $sitelist=$self->{dslist};

  my $counter=0;
  for (my $s1=0; $s1 <=$#{$sitelist}; $s1++) {
    for (my $res1=0; $res1<=$#{$sitelist->[$s1]->{reslist}}; $res1++) {
      $counter++;
    }
  }
  return $counter;
}


=head2

 Title   : countResPairTotal
 Usage   : my $num_res = $rotameter->countResPairTotal()
 Function: Counts the total number of residue pairs in all site pairs.
           Synonymous to ROTAMER::countPairs().
 Returns :
 Args    :

=cut

sub countResPairTotal {
  my $self=shift;
  return $self->countPairs();
}


=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub countValidSelfTotal {
  my $self=shift;
  my $counter=0;
  my $sitelist=$self->{dslist};


  for (my $s1=0; $s1 <=$#{$sitelist}; $s1++) {
    for (my $res1=0; $res1<=$#{$sitelist->[$s1]->{reslist}}; $res1++) {
      for(my $rot1=0; $rot1<=$#{$sitelist->[$s1]->{reslist}->[$res1]->{rotamerlist}}; $rot1++) {
        next if ($sitelist->[$s1]->{reslist}->[$res1]->{rotamerlist}->[$rot1]->{valid} == 0 );
          $counter++;
      }
    }
  }

  return $counter;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut



sub countSelfTotal {
  my $self = shift;
  my $counter = 0;
  my $sitelist = $self->{dslist};


  for (my $s1=0; $s1 <=$#{$sitelist}; $s1++) {
    for (my $res1=0; $res1<=$#{$sitelist->[$s1]->{reslist}}; $res1++) {
      for(my $rot1=0; $rot1<=$#{$sitelist->[$s1]->{reslist}->[$res1]->{rotamerlist}}; $rot1++) {
        $counter++;
      }
    }
  }
  return $counter;
}


=head2 

 Title   : writeSummary
 Usage   : $rotamer->writeSummary()
 Function: 
 Returns : 
 Args    : 

=cut


sub writeSummary {
    my $self=shift;

    my $s1=$self->countValidSelfTotal();
    my $s0=$self->countSelfTotal();
    
    my $p1=$self->countValidPairTotal();
    my $p0=$self->countPairTotal();
    my $t=$s1+$p1;
    my $tf=$self->countPairs();

    print "Before clash elimination, the total number of rotamers is $s0...\n";
    print "After clash elimination, the total number of rotamers is $s1...\n";
    print "Total number of pairs of valid rotamers is $p1...\n";
    print "Thus, total number of energy terms (pair + self) should be $t...\n";
    print "Total number of pair files should be $tf...\n";
}
	

=head2 

 Title   : 
 Usage   : 
 Function: given a sequence. write out a design input file to repack 
           that sequence  
 Returns : 
 Args    : 

=cut

sub writeNextConfig {

    my $self=shift;
    my $seq=shift;
    
    die "no new sequence is received.\n" if (!defined $seq);    
    my $fname=shift;
    $fname='rot_subset' if ( !defined $fname );

    my $dir=shift;
    $dir='Design' if (!defined $dir);
	
    chdir($dir);
    my $rotsub = &GENERAL::GetOutFH($fname);

    my $i=-1;
    foreach my $site (@{$self->{dslist}}) {
	$rotsub->printf("%1s %3d ", $site->{chain}, $site->{iresnum});
        $i++;
	my $res=$seq->[$i];
	$rotsub->printf("%3s %1s ", $res, '*');
	$rotsub->print("\n");
    }

    chdir("../");

}


=head2 

 Title   :  chopIt 
 Usage   :  $rotamer->chopIt(5) 
 Function:  divide the pair job into N pieces 
 Returns :  
 Args    : 

=cut


sub chopIt {
	
     my $self=shift;
     my $n=shift;

     $self->readRotState();
     my $p1=$self->countValidPairTotal();

     my @pos;
     my $length=int($p1/$n);

     for (my $m=1; $m<$n; $m++) {
   	  my $po=$m*$length;
	  push (@pos, $po);
     }

     push (@pos, $p1);
#     foreach (@pos) { print $_, "\n"; }


     my $pend=shift @pos;
     my $co=0;
     my $a={};
     my $sitelist=$self->{dslist};

     for (my $s1=0; $s1 <=$#{$sitelist}; $s1++) {
        for (my $res1=0; $res1<=$#{$sitelist->[$s1]->{reslist}}; $res1++) {
	   for (my $s2=$s1+1; $s2 <=$#{$sitelist}; $s2++) {
	      for (my $res2=0; $res2<=$#{$sitelist->[$s2]->{reslist}}; $res2++) {
		for(my $rot1=0; $rot1<=$#{$sitelist->[$s1]->{reslist}->[$res1]->{rotamerlist}}; $rot1++) {
		    next if ($sitelist->[$s1]->{reslist}->[$res1]->{rotamerlist}->[$rot1]->{valid} ==0 );
		    for(my $rot2=0; $rot2<=$#{$sitelist->[$s2]->{reslist}->[$res2]->{rotamerlist}}; $rot2++) {
			next if ($sitelist->[$s2]->{reslist}->[$res2]->{rotamerlist}->[$rot2]->{valid} ==0 );
			# every residue pair has a pair.inp
			$co++;
			if ( $co< $pend) {
			    next;
			} elsif ($co == $pend) {
			    my @start=($s1, $res1, $co);
			    my @end=($s2, $res2, $co);
			    push (@{$a->{start}}, \@start);
			    push (@{$a->{end}}, \@end);
			    $pend=shift @pos;
			} else {
			    next;
			}
		   }
		}
	     }
	  }
        }
      }

#     foreach my $kk (@{$a->{start}}) {
#	      print $kk->[0], "  ", $kk->[1], "  ", $kk->[2], "\n";
#     }

#     foreach my $kk (@{$a->{end}}) {
#              print $kk->[0], "  ", $kk->[1], "  ", $kk->[2], "\n";
#     }

      my $p30=$a->{start}->[0]->[0];
      my $p31=$a->{start}->[0]->[1];
      my $p40=$a->{end}->[0]->[0];
      my $p41=$a->{end}->[0]->[1];

      my $fnum=0;
      my @fs=(1);
      my $ii=0;

      for (my $s11=0; $s11 <=$#{$sitelist}; $s11++) {
	for (my $res11=0; $res11<=$#{$sitelist->[$s11]->{reslist}}; $res11++) {
	    for (my $s21=$s11+1; $s21 <=$#{$sitelist}; $s21++) {
		for (my $res21=0; $res21<=$#{$sitelist->[$s21]->{reslist}}; $res21++) {
                     $fnum++;			
		     if ($p30 == $s11 && $p31 == $res11 && $p40 == $s21 && $p41 == $res21) {	
			push (@fs, $fnum);
			my $te=$fnum+1;
			push (@fs, $te);
			$ii++;
			$p30 = shift @{$a->{start}->[$ii]};
			$p31 = shift @{$a->{start}->[$ii]};
			$p40 = shift @{$a->{end}->[$ii]};
			$p41 = shift @{$a->{end}->[$ii]};
#                        print $ii, "\n";
#                        print $p30, "  ", $p31, "  ", $p40, "  ", $p41, "\n";
		     } else {
			next;
		     }
                  }
              }
           }
        }

#       print "\n\n", $fnum, "\n\n";

        # get rid of extra one 
        pop @fs;

#        foreach (@fs) {print $_, "\n";}

         my $mm=($#fs+1)/2;

	 if ( $mm != $n) {
		print "try unoptimized way to distribute the pair calculations \n";
        	@fs=();
		my $t=$self->countPairs();
		my $length1=int($t/$n)+1;

        	push (@fs, 1);
        	for (my $i=1; $i<$n; $i++) {
             		my $po=$i*$length1;
             		push (@fs, $po);
	     		$po++;
             		push (@fs, $po);
        	}

        	push (@fs, $t);

	}  else {
#          do nothing
	} 

        

#	foreach (@fs) {print $_, "\n";}



        return \@fs;
}



=head2

 Title   :  siteStr
 Usage   :  ROTAMER::siteStr($site)
 Function:  Given a site returns a unique string representing it.
 Returns :  A string.
 Args    :  Reference to a design site.

=cut

sub siteStr {
  my $site = shift;
  return "$site->{chain}\_$site->{iresnum}";
}

1;
