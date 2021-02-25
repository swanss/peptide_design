# PDB package
# read/write/convert structure info
#

#
# Design Module: PDB 
#
# POD documentation - main docs before the code

=head1 NAME PDB

PDB - Structure object, with features

=head1 SYNOPSIS

someone needs to write it

=head1 DESCRIPTION

Functions: 
(1) read/write pdb file
(2) read/write charmm crd file  
(3) format converting: fix the difference in atom names and histine name 
    between different pdb, crd format
(4) select certain chains, and write out pdb/crd
(5) renumber residue number: for example
    numbering in rcsb pdb file: 54, 54A, 54B will be converted to  
    54, 55, 56. 
(6) translate renumbered pdb back to original style
(7) get interface residues
(8) calculate center of mass
(9) select certain residues, write out selected coord
(10) many more.... 

Data Structures:


(1) molecule information with substructures containing atom and residue
information as well as a residue lookup table and coordinate cache
arrays.

  chain[] -> { id atom[] res[] resinx[]
              xcoor[] ycoor[] zcoor[] }

        atom[] -> { atominx atomname resname resnum iresnum
               chain xcoor ycoor zcoor hyd }
        res[]  -> { resname resnum iresnum start end valid chain }


(2) lookup hash table for multiple chains

    chainlookup -> { chainid ... }

(3) currently selected chainlists

    pickedchain[] ->  


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


package PDB;
#$DEF_VER_CONTROL->{ver}->{PDB} = "1.0";
$DEF_VER_CONTROL->{date}->{PDB} = "06/17/04";

use strict;
use locale;
use IO::File;
use GENERAL;
use DEFINITIONS;
use CHARMM;
use Scalar::Util;
use Devel::Peek;
use PDL;
use PDL::MatrixOps;
use Scalar::Util qw(looks_like_number);
#use SEQUENCE;


use vars qw (  %_seqabbrev %_seqlong );

BEGIN {


  %_seqabbrev = ('GLY', 'G', 'PRO', 'P', 'ALA', 'A', 'VAL', 'V',
                 'LEU', 'L', 'ILE', 'I', 'MET', 'M', 'CYS', 'C',
                 'PHE', 'F', 'TYR', 'Y', 'TRP', 'W', 'HIS', 'H',
                 'HSD', 'H', 'HSE', 'H', 'HSP', 'H', 'HSC', 'H', 'LYS', 'K',
                 'ARG', 'R', 'GLN', 'Q', 'ASN', 'N', 'GLU', 'E', 
                 'ASP', 'D', 'SER', 'S', 'THR', 'T', 'BTN', 'B');

  %_seqlong   = ('G', 'GLY', 'P', 'PRO', 'A', 'ALA', 'V', 'VAL',
                 'L', 'LEU', 'I', 'ILE', 'M', 'MET', 'C', 'CYS',
                 'F', 'PHE', 'Y', 'TYR', 'W', 'TRP', 'H', 'HIS',
                 'K', 'LYS', 'R', 'ARG', 'Q', 'GLN', 'N', 'ASN',
                 'E', 'GLU', 'D', 'ASP', 'S', 'SER', 'T', 'THR',
                 'Z', 'APT', 'X', 'APM', 'B', 'BTN');
}


=head2   FUNCTION -- new

 Title   :  new
 Usage   :  new([PDBfile|CRDfile ])     pdb file or crd file is optional
 Function:  creates a new Molecule object and reads a PDB or CRD structures if a file name is given
 Example :
 Returns :  structure object
 Args    :  1. PDB/CRD file name
            2. Optional: format - "PDB" or "CRD". By default will try to determine from the file extension.
            3. Optional: translate - will be passed to readPDB or readCRD for any translations.
            4. Optional: read only the first n residues. So far only supported for readPDB.

=cut

sub new {

    # read in the file name
    my $farg = shift;
    my $format = shift;
    my $translate = shift;
    if (!defined($translate)) { $translate = ""; }
    my $firstn = shift;

    if (defined($format)) {$format = uc($format);}
    # Determine input file format (if the file is specified)
    # If format specified, it has to be either pdb or crd. If not specified, the file
    # name has to end with either pdb or crd.
    if (defined($farg)) {
      if (defined($format) && ($format !~ /PDB|CRD/)) {
        GENERAL::error("Error initializing PDB object: invalid format name \"$format\".");
      } elsif (!defined($format)) {
        if ($farg =~ /(\.pdb\Z|\.ent\Z)/) { $format = "PDB"; }
        elsif ($farg =~ /\.crd\Z/) { $format = "CRD"; }
        else { GENERAL::error("Error initializing PDB object: could not determine format from filename \"$farg\""); }
      }
    }
    my $self = {};

    # initialize chain list
    $self->{chain} = ();

    # initialize chain lookup table hash
    $self->{chainlookup} = {};

    # picked chains and interface residue lists are set to null
    $self->{pickedchains} = undef;
    $self->{coordfile} = undef;
    $self->{water} = ();

    bless $self;

    # basd on the format, use approriate method to initialized the object
    if (defined($farg)) {
      $self->{coordfile} = $farg;
      $self->readPDB($farg, $translate, $firstn, @_) if ($format =~ /PDB/);
      $self->readCRD($farg, $translate) if ($format =~ /CRD/);
    }
    return $self;
}


=head2   FUNCTION -- DESTROY

 Title   : DESTROY
 Usage   : Called automatically by Perl's Garbage Collection
 Function: Destroys the PDB object. This function is neccesarry because
           PDB objects contain circular references and hence will not be
           properly handeled by the default Garbage Collection system.
 Example :
 Returns : nothing
 Args    : none

=cut

sub DESTROY {
  my $self = shift;
  
  # Go through all chains
  foreach my $c (@{$self->{chain}}) {
    # Go though all residues
    foreach my $res (@{$c->{res}}) {
      # Go through all atoms
      foreach my $atom (@{$res->{atom}}) {
        # destroy the atom
        undef(%$atom);
      }
      # destroy the residue
      undef(%$res);
    }
    # destroy the chain
    undef(%$c);
  }
  # destroy the structure
  undef(%$self);
}


=head2   FUNCTION -- readPDB

 Title   :  readPDB
 Usage   :  readPDB(file,translate)
 Function:  reads a protein structures from a PDB file. translate may be set to
            CHARMM19  for proper recognition of histidine residues.
 Example :  $mol->readPDB('pdb1brs.ent', CHARMM19);
 Returns :  PDB object
 Args    :  1. file --> pdb file name
            2. translate --> optional
            3. n -> read only the first n residues
            4. Various options in a hash:
               'norename'    1 or 0 (default is 0). If set, chains are not going to be renamed even if they continue after TER or multiple chains have the same name.

=cut

sub readPDB {
  my $self = shift;
  my $pdbname = shift;
  my $ifh = GENERAL::GetInFH($pdbname);
  my $translate = shift;
  if (!defined $translate) {
    $translate = "";
  }
  my $firstn = shift;
  my %opts = @_;

  $self->{coordfile} = $pdbname;
  $self->{chain} = ();
  $self->{chainlookup} = {};
  $self->{pickedchains} = undef;


  my $lastresnum = -999999;
  my $lastresname = "";
  my $lasticode = "";
  my $iresnum = 0;
  my $uresnum = 0; # uniform residue number (does not restart for each chain)
  my $atomunx = 0;
  my $lastchain = "";
  my $ter = 0;
  my $chainrec;
  my $resrec;
  my $usesegid = 0;
  my %unnat;
  my $atomflg = "^ATOM";
  $atomflg = "^(ATOM|HETATM)" if ($translate =~ /HETATM/i);
  $usesegid = 1 if ($translate =~ /USESEGID/);

  while (<$ifh>) {
    last if (/^(END|End|\#End)/);
    if (/^TER/) { $ter = 1; next; }
    next if (! /$atomflg/);

    # here we skip the line with alternative atom position
#    next if (substr($_,16,1) =~ /\S/);

    # now read in PDB record
    if (/^$atomflg/) {
        $atomunx++;
        # sometimes PDB lines are too short (if they do not contain some of the
        # last optional columns). We don't want to read past the end of the string!
        $_ .= " " x 100;
        my ($atominx, $atomname, $alt, $resname, $chain, $resnum, $icode, $x, $y, $z, $seg, $B, $occ) = (($translate =~ /CHARMM/) ?
          (substr($_, 6, 5), substr($_, 12, 4), substr($_,16, 1), substr($_,17, 3),
           substr($_, 21, 1), substr($_, 23, 4), " ", substr($_, 30, 8), substr($_, 38,8),
           substr($_,46,8), substr($_, 72, 4), substr($_, 60, 6), substr($_, 54, 6)) :
          (substr($_, 6, 5), substr($_, 12, 4), substr($_,16, 1), substr($_,17, 4), # TIP3, for example, takes 4 characters; according to PDB format, position 21 is not used, so...
           substr($_, 21, 1), substr($_, 22, 4), substr($_, 26, 1), substr($_, 30, 8), substr($_, 38,8),
           substr($_,46,8), substr($_,72, 4), substr($_, 60, 6), substr($_, 54, 6)));
        for ($atominx, $atomname, $alt, $resname, $chain, $resnum, $x, $y, $z, $seg, $B, $occ) { s/^\s+//; s/\s+$//; }
        $unnat{$resname} = 1 if (!exists $_seqabbrev{$resname});
        my $het = (/^HETATM/ ? 1 : 0);

        # use segment ID's instead of chain ID's?
        if ($usesegid) {
          $chain = $seg;
        } elsif (($chain eq "") && ($seg =~ /^[A-Za-z0-9]/)) {
          # Use first character of segment name in two cases:
          # - if this is the first atom (no previous chain name), there is no current chain name and segment name starts with a letter
          # - if previously decided to use segment name and it starts with a letter
          $chain = substr($seg, 0, 1);
        }
        if (($chain eq "") && ($ter == 0)) { $chain = $lastchain; }

        # create a new chain object
        if (($chain eq "") || ($chain ne $lastchain) || $ter) {
          my $ranOut = 0;
          if (($chain eq $lastchain) && $ter) {
            # If the chain name continues after termination, we need to assign a new chain name, but remember the name
            # that was actually read, since this name is what has to be used to determine when the next chain comes
            GENERAL::warning("Chain name '$chain' continues after TER - assigning a new name (in $self->{coordfile})!") if (!defined($opts{norename}) || ($opts{norename} == 0));
            $chainrec = $self->_newChain($chain, $opts{norename}, \$ranOut);
            print "Name assigned: $chainrec->{id}\n" if (!defined($opts{norename}) || ($opts{norename} == 0));
          } else {
            $chainrec = $self->_newChain($chain, $opts{norename}, \$ranOut);
            # Chain needs to be the actual assigned one if $chain is empty
            $chain = $chainrec->{id} if ($chain eq "");
          }
          # start to count residue numbers in this chain
          $iresnum = 0;
          $lastresnum = -999999;
          $lastresname = "";
          $ter = 0;
          $chainrec->{segid} = $seg if (!$ranOut);
          $chainrec->{osegid} = $seg;
        }

        # fixing some generic problems:
#        if (defined $translate && $translate =~ /CHARMM/ ) {
#          $atomname="CD1"
#            if ($resname eq "ILE" && $atomname eq "CD");
#          $atomname="O"
#            if ($atomname eq "OT1" || $atomname eq "O1" || $atomname eq "OCT1");
#          $atomname="OXT"
#            if ($atomname eq "OT2" || $atomname eq "O2" || $atomname eq "OCT2");
#        } else {
#          $atomname="CD"
#            if ($resname eq "ILE" && $atomname eq "CD1");
#        }
        if (defined($translate) && ($translate =~ /CHARMM19/)) {
          $resname =~ s/HSE/HSD/;
          $resname =~ s/HSD/HIS/;
        }
        $atomname = "CD" if ($resname eq "ILE" && $atomname eq "CD1" && !(defined($translate) && ($translate =~ /NONE/)));

        # if necessary, make a new residue
        my $new_f = 1; # is this a truely new atom, as opposed to an alternative position?
        if (($resnum ne $lastresnum) || ($resname ne $lastresname) || (($translate =~ /ICODE/i) && ($icode ne $lasticode)))  {
          $iresnum++; $uresnum++;
          # stop reading if only supposed to read first several residues
          if (defined($firstn) && ($uresnum > $firstn)) { last; }
          my $nresrec = PDB::_newResidue($resname, $resnum, $iresnum, $uresnum, $chainrec);
          $nresrec->{icode} = $icode;
          $resrec = $nresrec;
        } elsif ($alt =~ /\S/) {
          # if this is not a new residue AND the alternative location flag is specified,
          # figure out if another location for this atom has already been given. If not,
          # then treat this as the "primary" location, and whatever other locations
          # are specified will be treated as alternative
          my $a = PDB::getAtomInRes($resrec, $atomname, -1);
          if ($a ne -1) {
            $new_f = 0;
            my $alta = PDB::_newAtom($atominx, $atomunx, $atomname, $x, $y, $z, $B, $occ, 0.0, undef, undef, $seg);
            $alta->{altid} = $alt;
            push(@{$a->{altpos}}, $alta);
          }
        }
        # make a new atom
        PDB::_newAtom($atominx, $atomunx, $atomname, $x, $y, $z, $B, $occ, 0.0, $resrec, $het, $seg) if ($new_f);
        $lastresnum = $resnum;
        $lasticode = $icode;
        $lastresname = $resname;
        $lastchain = $chain;
    }
  }
  close($ifh) if (!ref($pdbname));
#  if (scalar(keys(%unnat)) > 0) {
#    GENERAL::warning("PDB file $pdbname contains some unnatural amino acids:\n" . join("\n", keys(%unnat)) . "\n");
#  }

  $self->selectAllChains();
}


=head2   FUNCTION -- readCRD

 Title   :  readCRD
 Usage   :  readCRD(file,translate)
 Function:  reads a protein structures from a charmm crd file. translate may be set to
            CHARMM19  for proper recognition of histidine residues.
 Example :  $mol->readPDB('pdb1brs.crd', 'CHARMM19');
 Returns :  PDB structure object
 Args    :  file-->crd file name
            translate--> optional

=cut

sub readCRD {

    my $self = shift;
    my $crdname = shift;
    my $fname = GENERAL::GetInFH($crdname);
    my $translate = shift;
    my $translation;
    if (!defined $translate) {
    $translation = "";
    } else {
    $translation = uc $translate;
    }

    $self->{coordfile}=$crdname;
    $self->{chain}=();
    $self->{chainlookup}={};
    $self->{pickedchains}=undef;


    my $lastresnum=-999;
    my $iresnum=0;
    my $uresnum=0; # uniform residue number (does not restart for each chain)
    my $atomunx=0;
    my $lastchain=".";

    my $chainrec;
    my $resrec;
    my $linecounter=0;

  READCRD:
    while(<$fname>) {

      # we skip the empty line and comments
      next if (/^$/);
      next if ( /^\*/);
      next if ( /TER/);
  
      $linecounter++;
      # we also skip total number line
      next if ($linecounter == 1);
  
      # now read in CRD record
      chomp;
  
      $atomunx++;
      my ($atominx, $resname, $atomname, $x, $y, $z, $resnum, $charge)
      = GENERAL::Trim(substr($_, 0,5), substr($_, 11,4),
      substr($_, 16,4), substr($_, 20,10), substr($_, 30,10),substr($_, 40,10), substr($_, 56,6), substr($_, 62,8));
  
      my $chain = GENERAL::Trim(substr($_,51,4));
      my $resid = GENERAL::Trim(substr($_,55,4)) + 0 ;
  
      if ($chain eq "") {
          # in case we don't have segment name, we assign the chain as 'A'
          $chain='A';
      }
  
      # create a new chain object
      if ($chain ne $lastchain) {
          $chainrec= $self->_newChain($chain);
          # start to count residue numbers in this chain
          $iresnum=0;
          $lastresnum = -999999;
      }
  
      # fixing some generic problems:
      $atomname="CD1"
          if ($resname eq "ILE" && $atomname eq "CD");
#      $atomname="O"
#          if ($atomname eq "OT1" || $atomname eq "O1" || $atomname eq "OCT1");
#      $atomname="OXT"
#          if ($atomname eq "OT2" || $atomname eq "O2" || $atomname eq "OCT2");
  
      if ($translation eq "CHARMM19") {
        $resname=~s/HSE/HSD/;
        $resname=~s/HSD/HIS/;
      }
  
      # start to record residue information if we see a new residue
      if ($resnum ne $lastresnum)  {
          $iresnum++; $uresnum++;
      }
  
      # if necessary, make a new residue
      if ($resnum ne $lastresnum) {
        my $nresrec = PDB::_newResidue($resname, $resnum, $iresnum, $uresnum, $chainrec);
  #      my $nresrec = ();
  #      $nresrec->{chain}=$chainrec;
  #      $nresrec->{resname}=$resname;
  #      $nresrec->{resnum}=$resnum;
  #      $nresrec->{iresnum}=$iresnum;
  #      $nresrec->{uresnum}=$uresnum;
  #      $nresrec->{valid}=1;
  #      push(@{$chainrec->{res}}, $nresrec);
        $resrec = $nresrec;
      }
      
      # make new atom information
      my $pdbrec = PDB::_newAtom($atominx, $atomunx, $atomname, $x, $y, $z, 0.0, 1.0, $charge, $resrec);
  #    my $pdbrec={};
  #    $pdbrec->{atominx}=$atominx+0;
  #    $pdbrec->{atomunx}=$atomunx;
  #    $pdbrec->{atomname}=$atomname;
  #    $pdbrec->{xcoor}=$x+0.0;
  #    $pdbrec->{ycoor}=$y+0.0;
  #    $pdbrec->{zcoor}=$z+0.0;
  #    $pdbrec->{B}=0.0;
  #    $pdbrec->{occ}=1.0;
  #    $pdbrec->{hyd}=($atomname=~/^H.*/)?1:0;
  #    $pdbrec->{charge}=$charge+0.0;
  #    $pdbrec->{residue} = $resrec;
  #    push(@{$resrec->{atom}}, $pdbrec);
  
      $lastresnum = $resnum;
      $lastchain = $chain;
    }

    $self->selectAllChains();

    close($fname);
}


=head2   FUNCTION -- _newChain (Internal Function)

 Title   :  _newChain
 Usage   :  _newChain(chain_id)
    Example :  my $c = self->_newChain('A');
 Function:  Initialize a new chain object
 Returns :  a chain object
 Args    :  chain_id

=cut


sub _newChain {
  my $self = shift;
  my $cid = shift;
  my $nonew = shift; $nonew = 0 if (!defined($nonew));
  my $ranoutPtr = shift;

  # if chain name is empty, make a new one
  if (($cid =~ /^\s*$/) || (defined($self->{chainlookup}->{$cid}))) {
    if ((defined($self->{chainlookup}->{$cid}))) {
      if ($nonew == 0) {
        GENERAL::warning("Chain $cid already exists - coming up with new name..."); $cid = "";
     } else {
       return $self->{chainlookup}->{$cid};
     }
    }
    my @names = ('A' .. 'Z', 0 .. 9); my $f = 0;
    foreach my $pc (@names) {
      if (!defined($self->{chainlookup}->{$pc})) { $cid = $pc; $f = 1; last; }
    }
    if ($f == 0) {
      GENERAL::warning("Ran out of chain names for structure" . (defined($self->{coordfile}) ? " $self->{coordfile}" : "") . " - will use more than a character (PDB sctructure will be incorrect upon writing)!\n");
      $$ranoutPtr = 1 if (defined($ranoutPtr));
      foreach my $pc (@names) {
        for (my $i = 0; $i < 999; $i++) {
          if (!defined($self->{chainlookup}->{"$pc$i"})) { $cid = "$pc$i"; $f = 1; last; }
        }
        last if ($f);
      }
      GENERAL::assert($f, "Ran out of even multi-character chain names -- your PDB structure really has more than 36,000 chains???");
    }
  }

  my $chainrec = {};
  $chainrec->{id} = $cid;
  $chainrec->{segid} = $cid;
  my @tmp; $chainrec->{res} = \@tmp;
  $chainrec->{parent} = $self;

  push(@{$self->{chain}}, $chainrec);
  Scalar::Util::weaken($chainrec->{parent}); # need to weaken reference to avoid memory leak

  $self->{chainlookup}->{$chainrec->{id}} = $chainrec;

  return $chainrec;
}


=head2   FUNCTION -- _newResidue (Internal Function)

 Title   :  _newResidue
 Usage   :  _newResidue()
 Example :  my $res = PDB::_newResidue(...);
 Function:  Initializes a new residue object and links it into the given chain (if provided)
 Returns :  a new residue record
 Args    :  

=cut

sub _newResidue {
  my $resname = shift;
  my $resnum = shift;
  my $iresnum = shift;
  my $uresnum = shift;
  GENERAL::requireArgs($resname, $resnum, $iresnum, $uresnum);
  my $chainrec = shift;
  my $valid = shift;
  if (!defined($valid)) { $valid = 1; }
  my $prepend = shift; $prepend = 0 if (!defined($prepend));
          
  my $nresrec = ();
  $nresrec->{chain} = $chainrec if (defined($chainrec));
  $nresrec->{resname} = $resname;
  $nresrec->{resnum} = $resnum;
  $nresrec->{iresnum} = $iresnum;
  $nresrec->{uresnum} = $uresnum;
  $nresrec->{valid} = $valid;
  if (defined($chainrec)) {
    push(@{$chainrec->{res}}, $nresrec) if (!$prepend);
    unshift(@{$chainrec->{res}}, $nresrec) if ($prepend);
  }
  
  return $nresrec;
}

=head2   FUNCTION -- _newAtom (Internal Function)

 Title   :  _newAtom
 Usage   :  _newAtom()
 Example :  my $atom = PDB::_newAtom(...);
 Function:  Initializes a new atom object and links it into the given residue (if provided)
 Returns :  a new atom record
 Args    :  

=cut

sub _newAtom {
  my $atominx = shift;
  my $atomunx = shift;
  my $atomname = shift;
  my $x = shift;
  my $y = shift;
  my $z = shift;
  my $B = shift;
  my $occ = shift;
  my $charge = shift;
  GENERAL::requireArgs($atominx, $atomunx, $atomname, $x, $y, $z, $B, $occ, $charge);
  my $resrec = shift;
  my $het = shift; $het = 0 if (!defined($het));
  my $segID = shift;

  my $pdbrec={};
  $pdbrec->{atominx} = $atominx;
  $pdbrec->{atomunx} = $atomunx;
  $pdbrec->{atomname} = $atomname;
  $pdbrec->{xcoor} = $x+0.0;
  $pdbrec->{ycoor} = $y+0.0;
  $pdbrec->{zcoor} = $z+0.0;
  $pdbrec->{B} = (Scalar::Util::looks_like_number($B) ? $B : 0.0);
  $pdbrec->{occ} = (Scalar::Util::looks_like_number($occ) ? $occ : 0.0);
  $pdbrec->{hyd} = ($atomname=~/^H.*/) ? 1:0;
  $pdbrec->{charge} = $charge+0.0;
  $pdbrec->{hetero} = $het;
  if (defined($resrec)) {
    $pdbrec->{residue} = $resrec;
    push(@{$resrec->{atom}}, $pdbrec);
    $resrec->{altered} = 1;
  }
  $pdbrec->{segID} = $segID if (defined($segID));

  return $pdbrec;
}

=head2

 Title   :  insertChain
 Usage   :  insertChain
 Example :  $pdb->insertChain($chain);
 Function:  Inserts a new chain into the PDB object.
 Returns :  nothing
 Args    :  0. PDB structure.
            1. The new chain object.

=cut

sub insertChain {
  my $self = shift;
  my $chainrec = shift;

  if (defined($self->{chainlookup}->{$chainrec->{id}})) {
    GENERAL::error("Chain with ID $chainrec->{id} already exists in this PDB structure!");
  }
  # If this chain was previously part of some PDB structure, we need to make sure that this structure does not
  # get destroyed after its reference count is decremented due to this chain no longer pointing
  # to it - so we make this chain REALLY point to it in case this reference was weakened. If this reference
  # was not weakened, this will have no effect.
  my $pdb = $chainrec->{parent};
  if (ref($pdb)) {
    $chainrec->{parent} = $pdb;
    # Also, we need to make sure that the old PDB structure doesn't still think it has this chain
    if (defined($pdb->{chainlookup}->{$chainrec->{id}})) {
      delete $pdb->{chainlookup}->{$chainrec->{id}};
    }
    my @chains;
    foreach my $c (@{$pdb->{chain}}) {
      if ($c ne $chainrec) { push(@chains, $c); }
    }
    $pdb->{chain} = \@chains;
  }

  # Insert chain into new PDB
  $chainrec->{parent} = $self;
  push(@{$self->{chain}}, $chainrec);
  Scalar::Util::weaken($chainrec->{parent}); # need to weaken reference to avoid memory leak
  $self->{chainlookup}->{$chainrec->{id}} = $chainrec if ($chainrec->{id} ne "");
}

=head2

 Title   :  renameChain
 Usage   :  renameChain
 Example :  $pdb->renameChain($chain, $new_chain_name);
 Function:  Renames the given chain.
 Returns :  nothing
 Args    :  0. PDB structure.
            1. Chain object to rename or chain name corresponding to it.

=cut

sub renameChain {
  my $self = shift;
  my $chain = shift;
  my $nchname = shift;
  GENERAL::requireArgs($self, $chain, $nchname);

  if (!ref($chain)) {
    $chain = $self->getChain($chain);
  }

  delete($self->{chainlookup}->{$chain->{id}});
  $self->{chainlookup}->{$nchname} = $chain;
  $chain->{id} = $nchname;
}

=head2   FUNCTION -- _coorCache

 Title   :  _coorCache
 Usage   :  $self->_coorCache() or $self->_coorCache('A');
 Function:  cache the xyz coord of each atom (hash) into array
            indexed by atom number
 Returns :
 Args    :

=cut

sub _coorCache {

    my $self=shift;
    my $chainstr=shift;

    if (! defined $chainstr) {
    $self->selectAllChains();
    } else {
    $self->selectChains($chainstr);
    }


    foreach my $c ( @{$self->{pickedchains}} ) {
    my $a=$c->{atom};
    for (my $i=0; $i<=$#{$a}; $i++) {
        $c->{xcoor}->[$i]=$a->[$i]->{xcoor};
        $c->{ycoor}->[$i]=$a->[$i]->{ycoor};
        $c->{zcoor}->[$i]=$a->[$i]->{zcoor};
    }
    }
}


=head2  FUNCTION -- getChain

 Title   :  getChain
 Usage   :  $mol->getChain([$chainid])
 Function:  returns chain for given ID or first chain if
            no argument is given
 Returns :  reference to chain
 Args    :  chain id or null
 Example :  $ChainA = $mol->getChain('A');

=cut


sub getChain {
  my $self = shift;
  my $chainid = shift;
  GENERAL::requireArgs($self, $chainid);
  my $strict = shift;
  if (!defined($strict)) { $strict = 1; }

  if (!defined($self->{chainlookup}->{$chainid})) {
    if ($strict > 0) {
      GENERAL::error("Invalid chain name specified: $chainid!");
    } elsif ($strict == 0) {
      GENERAL::warning("Invalid chain name specified: $chainid!");
    }
    return -1;
  }

  return $self->{chainlookup}->{$chainid};
}



=head2  FUNCTION -- getResByInd

 Title   :  getResByInd
 Usage   :  $pdb->getResByInd($chain, $iresnum)
 Function:  Returns a reference to residue number $iresnum
            in chain $chain.
 Returns :  Reference to residue
 Args    :  1. chain id or chain reference
            2. residue number (starting from 1)
            3. strict flag (error or warning upon failure?)
            4. renumbered flag. If set to 1 (default) assumes the structure is
               renumbered and residue with index $iresnum-1 of the res array
               of the appropriate chain is what you are looking for. If set to 2,
               will go through all residues of the given chain and find the
               first with iresnum equal to what is given. If set to 3, will do the
               same, but look for the resnum entry, not iresnum.

=cut


sub getResByInd {
  my $pdb = shift;
  my $chain = shift;
  my $iresnum = shift;
  GENERAL::requireArgs($chain, $iresnum);
  my $strict = shift;
  if (!defined($strict)) { $strict = 1; }
  my $renum = shift;
  if (!defined($renum)) { $renum = 1; }

  my $cref;
  if (ref($chain)) { $cref = $chain; }
  else { $cref = $pdb->getChain($chain, $strict); }
  if ($cref eq -1) {
    if ($strict > 0) {
      GENERAL::error("chain $cref->{id} not found");
    } elsif ($strict == 0) {
      GENERAL::warning("chain $cref->{id} not found");
    }
    return -1;
  }
  if ($renum == 1) {
    if ((scalar(@{$cref->{res}}) < $iresnum) || ($iresnum <= 0)) {
      if ($strict > 0) {
        GENERAL::error("Residue number $iresnum out of range for chain $cref->{id} (it has " . scalar(@{$cref->{res}}) . " residues)");
      } elsif ($strict == 0) {
        GENERAL::warning("Residue number $iresnum out of range for chain $cref->{id} (it has " . scalar(@{$cref->{res}}) . " residues)");
      }
      return -1;
    }
    return $cref->{res}->[$iresnum-1];
  } elsif ($renum == 2) {
    my $fres = -1;
    foreach my $res (@{$cref->{res}}) {
      if ($res->{iresnum} eq $iresnum) { $fres = $res; last; }
    }
    if ($fres eq -1) {
      if ($strict > 0) {
        GENERAL::error("Could not find residue with iresnum $iresnum in chain $cref->{id}");
      } elsif ($strict == 0) {
        GENERAL::warning("Could not find residue with iresnum $iresnum in chain $cref->{id}");
      }
    }
    return $fres;
  } elsif ($renum == 3) {
    my $fres = -1;
    foreach my $res (@{$cref->{res}}) {
      if ($res->{resnum} eq $iresnum) { $fres = $res; last; }
    }
    if ($fres eq -1) {
      if ($strict > 0) {
        GENERAL::error("Could not find residue with resnum $iresnum in chain $cref->{id} ($pdb->{coordfile})");
      } elsif ($strict == 0) {
        GENERAL::warning("Could not find residue with resnum $iresnum in chain $cref->{id} ($pdb->{coordfile})");
      }
    }
    return $fres;
  } else {
    GENERAL::error("Unrecognized renum setting '$renum'");
  }
}


=head2  FUNCTION -- getAtomInRes

 Title   :  getAtomInRes
 Usage   :  getAtomInRes($res, $atomname)
 Function:  returns a reference to atom named $atom in residue $res
 Returns :  reference to atom
 Args    :  1. reference to residue
            2. atom name string
            3. optional: strict flag. If true, halts on failure. Default is true.

=cut

sub getAtomInRes {
  my $res = shift;
  my $aname = shift;
  GENERAL::requireArgs($res, $aname);
  my $strict = shift;
  if (!defined($strict)) { $strict = 1; }

  # If residue altered since next time, rebuild the name lookup hash
  if (defined($res->{altered})) { delete($res->{atomnamehash}); delete($res->{altered}); }
  if (!defined($res->{atomnamehash})) {
    foreach my $a (@{$res->{atom}}) {
      $res->{atomnamehash}->{$a->{atomname}} = $a;
    }
  }

  if (defined($res->{atomnamehash}->{$aname})) { return $res->{atomnamehash}->{$aname}; }

  if ($strict > 0) {
    GENERAL::error("Could not find atom named \"$aname\" in residue " . PDB::resStr($res));
  } elsif ($strict == 0) {
    GENERAL::warning("Could not find atom named \"$aname\" in residue " . PDB::resStr($res));
  }
  return -1;
}


=head2  FUNCTION -- getAtomI

 Title   :  getAtomI
 Usage   :  my $a = $pdb->getAtomI($chain_index, $residue_index, $atom_index)
 Function:  a simple getter of atoms by indices. E.g., getAtom(0, 0, 0) gets the first atom
            and getAtom(-1, -1, -1) the last one.
 Returns :  reference to atom
 Args    :  1. chain index (0-initiated)
            2. residue index (0-initiated)
            3. atom index (0-initiated)

=cut

sub getAtomI {
  my $self = shift;
  my $ci = shift;
  my $ri = shift;
  my $ai = shift;
  GENERAL::requireArgs($self, $ci, $ri, $ai);

  return $self->{chain}->[$ci]->{res}->[$ri]->{atom}->[$ai];
}


=head2  FUNCTION -- getAtomsMatchingName

 Title   :  getAtomsMatchingName
 Usage   :  getAtomsMatchingName($res, $regexp)
 Function:  returns a list of atoms with names matching the specified regular expression
 Returns :  list of atoms
 Args    :  1. reference to residue
            2. regular expression to be applied to atom names

=cut

sub getAtomsMatchingName {
  my $res = shift;
  my $namere = shift;
  GENERAL::requireArgs($res, $namere);
  my @list;

  foreach my $a (@{$res->{atom}}) {
    push(@list, $a) if ($a->{atomname} =~ /$namere/);
  }

  return @list;
}


=head2  FUNCTION -- getAtom

 Title   :  getAtom
 Usage   :  $pdb->getAtom(chain_name, residue_inumber, atom_name)
 Function:  returns a reference to the specified atom
 Returns :  reference to atom
 Args    :  1. chain name
            2. residue number
            3. atom name

=cut


sub getAtom {
  my $self = shift;
  my $cid = shift;
  my $iresnum = shift;
  my $aname = shift;
  GENERAL::requireArgs($self, $cid, $iresnum, $aname);
  my $strict = shift;
  if (!defined($strict)) { $strict = 1; }
  my $renum = shift;
  if (!defined($renum)) { $renum = 1; }

  my $res = $self->getResByInd($cid, $iresnum, $strict, $renum);
  # let the appropriate string handling be done by PDB::getResByInd. If we get to this point, that means the strict level did not generate an error
  return -1 if ($res eq -1);
  return PDB::getAtomInRes($res, $aname, $strict);
}


=head2  FUNCTION -- isBondedTo

 Title   :  isBondedTo
 Usage   :  PDB::isBondedTo($resi, $resj)
 Function:  determines whether residue $resi is bonded to residue $resj (i.e., $resi-$resj in the N->C direction)
 Returns :  true (1) or false (0), accordingly. if backbone atoms necessary to make this determination are missing, returns -1
 Args    :  1. putative N-terminal residue
            2. putative C-terminal residue

=cut

sub isBondedTo {
  my $resi = shift;
  my $resj = shift;

  my $c = PDB::getAtomInRes($resi, "C", -1);
  my $n = PDB::getAtomInRes($resj, "N", -1);
  return -1 if (($c eq -1) || ($n eq -1));
  return (PDB::atomDist($c, $n) <= 2.0);
}

=head2  FUNCTION -- splitIntoChains

 Title   :  splitIntoChains
 Usage   :  my $pdbNew = $pdb->splitIntoChains()
 Function:  Interprets the given PDB object as a list of residues and splits this list into chains by connectivity between residues.
 Returns :  New PDB object with one or more chains.
 Args    :  None

=cut

sub splitIntoChains {
  my $pdb = shift;
  my $npdb = PDB::new();
  my @residues = $pdb->ConRes();
  return $npdb if (scalar(@residues) == 0);

  my $nc = $npdb->_newChain("");
  my $nr = PDB::_newResidue("", 1, 1, 1, $nc); PDB::copyResidue($residues[0], $nr);
  for (my $i = 1; $i < scalar(@residues); $i++) {
    $nc = $npdb->_newChain("") if (!PDB::isBondedTo($residues[$i-1], $residues[$i]));
    my $nr = PDB::_newResidue("", 1, 1, 1, $nc); PDB::copyResidue($residues[$i], $nr);
  }

  return $npdb;
}


=head2   FUNCTION -- selectChains

 Title   :  selectChains
 Usage   :  $mol->selectChains('ABC');
            activate A, B, C chains
 Function:  select active chains
 Returns :  no return value
 Args    :  a string of chainid

=cut

sub selectChains {
    my $self = shift;
    my $chainstr = shift;
    GENERAL::requireArgs($self, $chainstr);

    chomp $chainstr;
    my @chainlist = split(" ", $chainstr);
    undef(@{$self->{pickedchains}}) if (defined $self->{pickedchains});

    foreach (@chainlist) {
      push (@{$self->{pickedchains}}, $self->getChain($_));
    }
}


=head2   FUNCTION -- selectALLChains

 Title   :  selectAllChains
 Usage   :  $mol->selectAllChains()
            activate all chains
 Function:  select all chains in structure as active chains
 Returns :  set $self->{pickedchains} as all chains
 Args    :  no argement

=cut

sub selectAllChains {
  my $self=shift;

  undef(@{$self->{pickedchains}}) if (defined $self->{pickedchains});

  foreach my $c (@{$self->{chain}}) {
    push (@{$self->{pickedchains}}, $self->getChain($c->{id}));
  }
}



=head2   FUNCTION -- getChainlist

 Title   :  getChainlist
 Usage   :  $mol->getChainlist();
 Function:  return the array of chain names
 Returns :  array of chains
 Args    :  none

=cut

sub getChainList {
    my $self = shift;
    my @chainlist;

    foreach my $c (@{$self->{chain}}) {
        push (@chainlist, $c->{id});
        #print $c->{id},"\n";
    }

    return \@chainlist;
}


=head2   FUNCTION -- atomStr

 Title   :  atomStr
 Usage   :  PDB::atomStr($atom);
 Function:  Returns a string defining the given atom
 Returns :  Returns a string defining the given atom
 Args    :  1. Atom reference
            2. Optinal: if true, the atom string corresponding to the
               renumbered structure will be returned. If false, the
               given structure will be considered. Default is true.
=cut

sub atomStr {
  my $atom = shift;
  my $flag = shift;
  if (!defined($flag)) { $flag = 1; }

  if (defined($atom->{tag})) {
    return (PDB::resStr($atom->{residue}, $flag) . "\_$atom->{atomname}\_$atom->{tag}");
  } else {
    return (PDB::resStr($atom->{residue}, $flag) . "\_$atom->{atomname}");
  }
}

=head2 FUNCTION -- resStr

 Title   : resStr
 Usage   : PDB::resStr($res)
 Function: Returns the standard residue string. The format is:
            "<chain>_<residue name>_<residue number>"
 Returns : a string.
 Args    : 1. Reference to residue
           2. Optinal: if true, the residue string corresponding to the
              renumbered structure will be returned. If false, the
              given structure will be considered. Default is true.

=cut

sub resStr {
  my $res = shift;
  my $flag = shift;
  if (!defined($flag)) { $flag = 1; }

  if ($flag) {
    my $s = "$res->{chain}->{id}\_$res->{resname}\_$res->{iresnum}";
    if (defined($res->{tag})) { $s .= "\_$res->{tag}"; }
    return $s;
  } else {
    my $s = "$res->{chain}->{id}\_$res->{resname}\_$res->{resnum}";
    if (defined($res->{tag})) { $s .= "\_$res->{tag}"; }
    return $s;
  }
}


=head2   FUNCTION -- getResMap

 Title   :  getResMap
 Usage   :  $resmap=$mol->getResMap();
 Example :


 Function:  get the map between old numbering and new numbering
 Returns :  a data structure for looking up residue number
 Args    :  no argument
 Data Structure:

    hashes of hashes of hash

    $resmap ->{chainid}->{'tocrystal'}->{$iresnum};
    or
    $resmap->{chainid}->{'fromcrystal'}->{$resnum};


=cut

sub getResMap {

    my $self=shift;
    my $resmap={};

    $self->selectAllChains();
    foreach my $c (@{$self->{pickedchains}}) {

    # do some initialization
    $resmap->{$c->{id}}={};
    $resmap->{$c->{id}}->{'tocrystal'}={};
    $resmap->{$c->{id}}->{'fromcrystal'}={};

    foreach my $r ( @{$c->{res}} ) {
        # build up two way lookup hash
        my $key=$r->{iresnum};
        my $value=$r->{resnum};
        $resmap->{$c->{id}}->{'tocrystal'}->{$key}=$value;
        $resmap->{$c->{id}}->{'fromcrystal'}->{$value}=$key;
    }
    }

    return $resmap;
}


=head2   FUNCTION -- writePDB

 Title   : writePDB
 Usage   : $mol->writePDB("test.pdb", charmm19noh, "ABC");
           this will write out charmm19 format pdb without hydrogen
           for A, B, C three chains. if no chain argument is given,
           it will write out all chains
 Example :
 Function: writes out the current structure in PDB format
           the output format may be specified through tranlsate.
           Possible formats are CHARMM19, CHARMM22,
           and GENERIC with combination of NOH (No Hydrogen) and Renumber
           (renumber residue number in case we have alternative numbering)
           A chain ID may be given as the third argument
 Returns :
 Args    : output file name, format, chainid strings
           for example -->  "charmm22 noh renumber" can be your format string
           (1) the order of your format string doesn't matter
           (2) don't care if it is uppercase or lowcase
           (3) noh and renumber is optional
           (4) don't care the white space between each word

=cut

sub writePDB {

    my $self = shift;
    my $pdbname = shift;
    my $fname = GENERAL::GetOutFH($pdbname);
    my $translate = shift;
    die "Error in writePDB - must specifiy a pdb format!\n" if (!defined $translate);

    $translate = uc($translate);
    my $chainstr = shift;


    if (!defined $chainstr) {
      $self->selectAllChains();
    } else {
      $self->selectChains($chainstr);
    }

    my $chaincounter=-1;
    foreach my $c (@{$self->{pickedchains}}) {
    if ($#{$c->{res}}>=0) {
        my $lastresnum=$c->{res}->[$#{$c->{res}}]->{iresnum};

        $chaincounter++;
        my $atomcounter=0;
        foreach my $res ( @{$c->{res}} ) {
          # Make a copy for this residue for translation
          my $tres;
          %{$tres}=%{$res};

          for my $a ( @{$res->{atom}} ) {
                # Make a copy for this atom for translation and connect it
        # to the translate residue
        my $ta;
        %{$ta}=%{$a};
        $ta->{residue} = $tres;

        $atomcounter++;
        # here we handle some dirty details of file format converting
        # charmm19 or charmm22
        if ($translate=~/CHARMM/) {
          $ta->{atomname}="CD"
          if ($ta->{residue}->{resname} eq "ILE" && $ta->{atomname} eq "CD1");
          $ta->{atomname}="OT1"
          if ($ta->{atomname} eq "O" && $ta->{residue}->{iresnum} eq $lastresnum);
          $ta->{atomname}="OT2"
          if ($ta->{atomname} eq "OXT" && $ta->{residue}->{iresnum} eq $lastresnum);
          $ta->{residue}->{resname}=~s/HOH/TIP3/;
        }

        if ($translate =~ /CHARMM19/) {
          $ta->{residue}->{resname}=~s/HSD/HIS/;
          $ta->{residue}->{resname}=~s/HSE/HSD/;
        } elsif ($translate =~ /CHARMM22/) {
          GENERAL::warning("Writing in CHARMM22 format. Residue \"HIS\" assumed to be HSD (neutral histidine)!");
          $ta->{residue}->{resname}=~s/HIS/HSD/;
        } elsif ($translate =~ /GENERIC/) {
          $ta->{residue}->{resname}=~s/HSD/HIS/;
          $ta->{residue}->{resname}=~s/HSE/HIS/;
          $ta->{atomname}="CD1"
          if ($ta->{residue}->{resname} eq "ILE" && $ta->{atomname} eq "CD");
        }


        # hydrogen or non-hydrogen
        if ($translate =~ /NOH/ ) {
          printf $fname "%s\n",&PDB::_pdbLine($ta)
          if ( $ta->{hyd} == 0 && ! ($translate  =~ /REN/));
          printf $fname "%s\n",&PDB::_renumberpdbLine($ta)
            if ( $ta->{hyd} == 0 && $translate  =~ /REN/);
        } else {
          printf $fname "%s\n",&PDB::_renumberpdbLine($ta) if ($translate =~ /REN/) ;
          printf $fname "%s\n",&PDB::_pdbLine($ta)  if (!($translate =~ /REN/));
        }

        if ($translate =~ /NOT/ ) {
        } else {
            printf  $fname "TER\n"
            if ($ta->{residue}->{iresnum} eq $lastresnum && $chaincounter <= $#{$self->{pickedchains}}
                &&  $ta->{atomunx} == $ta->{residue}->{atom}->[$#{$ta->{residue}->{atom}}]->{atomunx});
# $atomcounter == $#{$c->{atom}}+1
        }
          }
        }
      }
    }
    close($fname) if (!ref($pdbname));
}


=head2   FUNCTION -- writeCRD

 Title   : writeCRD
 Usage   : $mol->writeCRD("test.crd", charmm19noh, "ABC");
           this will write out charmm19 format crd without hydrogen
           for A, B, C three chains. if no chain argument is given,
           it will write out all chains. By default, all the chain residues
           have been renumbered. a necessary step for charmm input!
 Example :


 Function: writes out the current structure in CRD format
           the output format may be specified through tranlsate.
           Possible formats are CHARMM19, CHARMM22,
           with combination of NOH (No Hydrogen).
           A string of chain IDs may be given as the third argument
 Returns :
 Args    : output file name, format, chainid strings

=cut

sub writeCRD {
  my $self = shift;
  my $fname = GENERAL::GetOutFH(shift);
  my $translate = shift;
  die "Error in writeCRD - must specifiy a crd file format!\n" if (!defined $translate);
  $translate= uc($translate);
  my $chainstr=shift;
  my $total=0;

  if (! defined $chainstr) {
    $self->selectAllChains();
  } else {
    $self->selectChains($chainstr);
  }

  foreach my $c (@{$self->{pickedchains}}) {
    for my $res (@{$c->{res}}) {
      $total += scalar(@{$res->{atom}});
    }
  }

  printf $fname "%5d\n", $total;

  foreach my $c (@{$self->{pickedchains}}) {
    if ($#{$c->{res}}>=0) {
      my $lastresnum=$c->{res}->[$#{$c->{res}}]->{resnum};

      my $atomcounter=0;
      foreach my $res (@{$c->{res}}) {
        # Make a copy for this residue for translation
        my $tres;
        %{$tres}=%{$res};

        for my $a (@{$res->{atom}}) {
          # Make a copy for this atom for translation and connect it
          # to the translate residue
          my $ta;
          %{$ta}=%{$a};
          $ta->{residue} = $tres;

          # here we handle some dirty details of file format converting
          # charmm19 or charmm22
          $ta->{atomname}="CD" if ($ta->{residue}->{resname} eq "ILE" && $ta->{atomname} eq "CD1");

          if ($translate=~/CHARMM/) {
# I am commenting this out. We have to come up with a more consistent way of specifying termini.
# If we do it like here then when the terminal amino acid (having OT's instead of O's) is patched inside CHARMM
# there will be unspecified atoms (namely the O's) because the patch file has O's and not OT's.
# Gevorg 03/05/04.
#            $ta->{atomname}="OT1" if ($ta->{atomname} eq "O" && $ta->{residue}->{resnum} eq $lastresnum);
#            $ta->{atomname}="OT2" if ($ta->{atomname} eq "OXT" && $ta->{residue}->{resnum} eq $lastresnum);
            $ta->{residue}->{resname}=~s/HOH/TIP3/;
          }

          if ($translate =~ /CHARMM19/) {
            # Neutral HIS, proton on ND1 - called HSD in CHARMM22 and HIS in CHARMM19
            $ta->{residue}->{resname}=~s/HSD/HIS/;
            # Neutral HIS, proton on NE2 - called HSE in CHARMM22 and HSD in CHARMM19
            $ta->{residue}->{resname}=~s/HSE/HSD/;
            # Doubly protonated charged HIS - called HSP in CHARMM22 and HSC in CHARMM19
            $ta->{residue}->{resname}=~s/HSP/HSC/;
          } elsif ($translate =~ /CHARMM22/) {
            if ($ta->{residue}->{resname} =~ /HIS/) {
              GENERAL::warning("Writing in CHARMM22 format. Residue \"HIS\" assumed to be HSD (neutral histidine)!");
              $ta->{residue}->{resname}=~s/HIS/HSD/;
            }
          }


          # hydrogen or non-hydrogen
          if ($translate =~ /NOH/ ) {
            printf $fname "%s\n",&PDB::_crdLine($ta)
            if ($ta->{hyd} == 0);
          } else {
            printf $fname "%s\n",&PDB::_crdLine($ta);
          }
        }
      }
    }
  }
  close $fname;
}


=head2   FUNCTION -- writePQR

 Title   : PDB::writePQR
 Usage   : $pdb->writePQR($pqrfilename)
           PDB::writePQR(\@atoms, $pqrfilename)
 Function: Writes out the specified structure in PQR format to the
           specified file. The structure can either be given as a
           PDB object or as a reference to an array of atoms.
 Returns : Nothing.
 Args    : 1. Input structure (either PDB object or reference to an
              array of atoms).
           2. Output file name.
           3. Optional: reference to a double hash table of atomic
              radii. If not specified, it will be read in using
              PDB::sizeLookup().

=cut

sub writePQR {
  my $struct = shift;
  my $pqrfile = shift;
  GENERAL::requireArgs($struct, $pqrfile);
  my $R = shift;
  my $renum = shift;
  $renum = 1 if (!defined($renum));

  # Determine structure type
  my $inmode;
  if (ref($struct) =~ /PDB/) {
    $inmode = 1;
  } elsif (ref($struct) =~ /ARRAY/) {
    $inmode = 2;
  } else {
    GENERAL::error("Unexpected structure type (reference to " . ref($struct) . ", structure $struct)");
  }

  # Build atomic size lookup table
  if (!defined($R)) { $R = PDB::sizeLookup(); }

  # Open pqr file for output
  my $ofh = GENERAL::GetOutFH($pqrfile);

  my ($pdbline, $pqrline, $i);
  my $atoms;
  if ($inmode == 1) {
    $atoms = PDB::conAtoms($struct, undef, 1);
  } elsif ($inmode == 2) {
    $atoms = $struct;
  }

  $i = 0;
  foreach my $atom (@$atoms) {
    if ($renum) { $i++; }
    else { $i = $atom->{atomunx}; }
    # Write modified line to the output pqr file
    $pqrline = sprintf("%-6s%5d %-5s%3s%6d    %8.3f%8.3f%8.3f%10.5f%10.5f", "ATOM", $i, $atom->{atomname},
    $atom->{residue}->{resname}, $atom->{residue}->{iresnum}, $atom->{xcoor}, $atom->{ycoor}, $atom->{zcoor}, $atom->{charge},
    $R->{$atom->{residue}->{resname}}{$atom->{atomname}});
    print $ofh "$pqrline\n";
  }
  close($ofh);
}


=head2   FUNCTION -- writeDUMMY

 Title   : PDB::writeDUMMY
 Usage   : $pdb->writeDUMMY($pqrfilename)
           PDB::writeDUMMY(\@atoms, $pqrfilename)
 Function: Writes out the specified structure in DUMMY format (i.e. only coordinates) to the
           specified file. The structure can either be given as a PDB object or as a
           reference to an array of atoms.
 Returns : Nothing.
 Args    : 1. Input structure (either PDB object or reference to an
              array of atoms).
           2. Output file name.
=cut

sub writeDUMMY {
  my $struct = shift;
  my $outfile = shift;
  GENERAL::requireArgs($struct, $outfile);

  # Determine structure type
  my $atoms;
  if (ref($struct) =~ /PDB/) {
    $atoms = PDB::conAtoms($struct, undef, 1);
  } elsif (ref($struct) =~ /ARRAY/) {
    $atoms = $struct;
  } else {
    GENERAL::error("Unexpected structure type (reference to " . ref($struct) . ", structure $struct)");
  }

  # Write in dummy format
  my $ofh = GENERAL::GetOutFH($outfile);
  $ofh->printf("%d 1\n", scalar(@$atoms));
  foreach my $atom (@$atoms) {
    $ofh->printf("%12.6f%12.6f%12.6f\n", $atom->{xcoor}, $atom->{ycoor}, $atom->{zcoor});
  }
  close($ofh);
}


=head2   FUNCTION -- _pdbline

 Title   :  _pdbline
 Usage   :  _pdbline(reference_to_atom_record)
 Function:  write out pdb line in original PDB format
 Returns :
 Args    :

here is the example of several different pdb format:

with and without alternative site
ATOM   2047  OE1 GLU H 129A     31.454  30.769  23.165  1.00 36.17           O
ATOM    269  CG1 ILE L  30      40.467  20.562  89.890  1.00 31.03           C

with and without hydrogen:
ATOM     40  CG   LEU L  5      46.830  31.553 102.694  1.00  0.00
ATOM    461  HD22 ASN L 47      27.673  39.580  54.062  1.00  0.00


=cut


sub _pdbLine {
    my $pdbrec=shift;

    my $chainid=$pdbrec->{residue}->{chain}->{id};
    $chainid=" " if (!defined $chainid || $chainid eq "");
    my $segid=$pdbrec->{residue}->{chain}->{segid};
    $segid = $chainid if (!defined($segid));

#    my ($resnum, $icode);
#    if ($pdbrec->{residue}->{resnum} =~ /[A-Za-z]/) {
#      ($resnum, $icode) = ($pdbrec->{residue}->{resnum} =~ /([0-9]+)([A-Za-z])/);
#    } else {
#      $resnum=$pdbrec->{residue}->{resnum};
#      $icode=" ";
#    }
    my $icode = defined($pdbrec->{residue}->{icode}) ? $pdbrec->{residue}->{icode} : " ";

    # Atom name placement is different when it is 4 characters long
    my $atomname;
    if (length($pdbrec->{atomname}) < 4) { $atomname = sprintf(" %-3s", $pdbrec->{atomname}); }
    else { $atomname = sprintf("%4s", $pdbrec->{atomname}); }
    
    return sprintf("%6s%5d %-4s%1s%-4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s",
    ($pdbrec->{hetero} ? "HETATM" : "ATOM  "), $pdbrec->{atominx}, $atomname, defined($pdbrec->{altid}) ? $pdbrec->{altid} : " ", $pdbrec->{residue}->{resname},
    substr($chainid, 0, 1), $pdbrec->{residue}->{resnum}, $icode, $pdbrec->{xcoor}, $pdbrec->{ycoor}, $pdbrec->{zcoor}, $pdbrec->{occ}, $pdbrec->{B}, $segid);
    # Note: PDB claims that under the 2.2 format segment ID needs to be left justified, but SPDBV
    # does not understand it that way.
}


=head2   FUNCTION -- _renumberpdbLine

 Title   : _renumberpdbLine
 Usage   : _renumberpdbline(reference_to_atom_record)
 Function: write out pdb line with renumbered resid
 Returns :
 Args    :

=cut
sub _renumberpdbLine {
    my $pdbrec=shift;

    my $chainid=$pdbrec->{residue}->{chain}->{id};
    $chainid=" " if (!defined $chainid || $chainid eq "");
    my $segid=$pdbrec->{residue}->{chain}->{segid};
    $segid = $chainid if (!defined($segid));

    # Atom name placement is different when it is 4 characters long
    my $atomname;
    if (length($pdbrec->{atomname}) < 4) { $atomname = sprintf(" %-3s", $pdbrec->{atomname}); }
    else { $atomname = sprintf("%4s", $pdbrec->{atomname}); }

    return sprintf("%6s%5d %-4s %-4s%1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s",
    ($pdbrec->{hetero} ? "HETATM" : "ATOM  "), $pdbrec->{atomunx}, $atomname, $pdbrec->{residue}->{resname},
    substr($chainid, 0, 1), $pdbrec->{residue}->{iresnum}, $pdbrec->{xcoor}, $pdbrec->{ycoor}, $pdbrec->{zcoor}, $pdbrec->{occ}, $pdbrec->{B}, $segid);
}



=head2   FUNCTION -- _crdLine

 Title   : _crdLine
 Usage   : _crdLine(reference_to_atom_record)
 Function: write out crd line with charmm format
 Returns :
 Args    :

=cut

sub _crdLine {
    my $pdbrec=shift;

    my $charge;
    if (defined($pdbrec->{charge})) {$charge = $pdbrec->{charge};}
    else {$charge = 0;}

    # implementation of crd file is simple
    return sprintf "%5d%5s %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4s%10.5f",
    $pdbrec->{atominx},$pdbrec->{residue}->{iresnum},
    $pdbrec->{residue}->{resname},$pdbrec->{atomname},
    $pdbrec->{xcoor}, $pdbrec->{ycoor}, $pdbrec->{zcoor},
        $pdbrec->{residue}->{chain}->{id}, $pdbrec->{residue}->{iresnum}, $charge;

}


=head2   FUNCTION -- ConRes

 Title   : ConRes
 Usage   : my @residues = $pdb->ConRes()
 Function: Concatinates all the residues of all chains into one array and returns it.
 Returns : An array of all residues.
 Args    :

=cut

sub ConRes {
  my $pdb = shift;
  my $chains = shift;
  if (defined($chains)) {
    $pdb->selectChains($chains);
  } else {
    $pdb->selectAllChains();
  }
  my @residues;

  foreach my $chain (@{$pdb->{pickedchains}}) {
    push(@residues, @{$chain->{res}});
  }
  return @residues;
}


=head2   FUNCTION -- getResidues

 Title   : ConRes
 Usage   : my $residues = $pdb->getResidues("A1 A3 B14")
 Function: Returns a list of specified residues.
 Returns : A reference to an array of selected residues.
 Args    : 1. A string representing a space-separated list of residues to select.

=cut

sub getResidues {
  my $pdb = shift;
  my $resStr = shift;
  GENERAL::requireArgs($pdb, $resStr);
  my $setAsValid = shift;
  my @resStrs= split(" ", GENERAL::Trim($resStr));
  my @res;

  foreach my $res (@resStrs) {
    GENERAL::assert(($res =~ /^(\D+)(\d+)$/) ? 1 : 0, "could not parse position '$res'");
    my $cid = $1;
    my $rid = $2;
    my $res = $pdb->getResByInd($cid, $rid, 1, 3);
    $res->{valid} = 1 if (defined($setAsValid));
    push(@res, $res);
  }

  return \@res;
}

=head2   FUNCTION -- conAtoms

 Title   : conAtom
 Usage   : my @atoms = $pdb->conAtom()
           my @atoms = $pdb->conRes('A B C')
 Function: Concatinates all atoms of all chains (or just the specified
           chains if given) into one array and returns it.
 Returns : An array of atoms.
 Args    : (1) Optional: a space separated list of directives. Directives
               can be either about which chains to include/not to include,
               or which residues to include/not to include, or whether only
               the backbones or only the sidechains of certain residues are
               to be included.
               'i_A' - include chain A only.
               'i_A_14' - include residue A14 only.
               'i_A_14 i_A_16 i_B_20' - include only residue A14, A16, and A20.
               'n_A' - include all chains except A.
               'n_A_14' - include all residues except A14.
               'n_A_14 n_B_5' - include all residues except A14 and B5.
               's_A_14' - for residue A14, only include the sidechain.
               'b_A_20' - for residue A20 only include the backbone.
               All combinations of directives are allowed although some make more
               sense than others.
           (2) Optional: pointer flag - if set, a pointer to an array of atoms will
               be returned.
=cut

sub conAtoms {
  my $pdb = shift;
  my $rules = shift; if (!defined($rules)) { $rules = ""; }
  my $pflag = shift;

  my @rules = split(" ", $rules);
  my %cd;
  my $pcd = 0; # number of positive chain directives
  my $ncd = 0; # number of negative chain directives
  my $cd = 0; # total number of chain directives
  my %rd;
  my $prd = 0; # number of positive residue directives
  my $nrd = 0; # number of negative residue directives
  my $rd = 0; # total number of residue directives

  # Build rules
  foreach my $r (@rules) {
    my @comp = split("_", $r);
    if (scalar(@comp) == 1) {
      GENERAL::error("Invalid rule: $r");
    } elsif (scalar(@comp) == 2) {
      # This is a chain directive, check if the chain exists
      my $c = $pdb->getChain($comp[1]);
      if (!defined($c)) {
        GENERAL::error("Invalid rule: $r. Bad chain name $comp[1]");
      }
      if ($comp[0] eq "i") {
        $cd{$comp[1]} = 1;
        $pcd++;
      } elsif ($comp[0] eq "n") {
        $cd{$comp[1]} = 0;
        $ncd++;
      } else {
        GENERAL::error("Invalid rule: $r. Bad chain modifier $comp[0]");
      }
      $cd++;
    } elsif (scalar(@comp) == 3) {
      # This is a residue rule, check if the chain exists
      my $c = $pdb->getChain($comp[1]);
      if (!defined($c)) {
        GENERAL::error("Invalid rule: $r. Bad chain name $comp[1]");
      }
      # Check if the residue exists
      if ($pdb->getResByInd($c, $comp[2], 0, 2) eq 0) {
#      if (scalar(@{$c->{res}}) < $comp[2]) {
        GENERAL::error("Invalid rule: $r. Bad residue number $comp[2]");
      }
      if ($comp[0] eq "i") {
        $rd{"$comp[1]\_$comp[2]"} = 1;
        $prd++;
      } elsif ($comp[0] eq "n") {
        $rd{"$comp[1]\_$comp[2]"} = 0;
        $nrd++;
      } elsif ($comp[0] eq "s") {
        $rd{"$comp[1]\_$comp[2]"} = 2;
        $nrd++;
      } elsif ($comp[0] eq "b") {
        $rd{"$comp[1]\_$comp[2]"} = 3;
        $nrd++;
      } else {
        GENERAL::error("Invalid rule: $r. Bad residue modifier $comp[0]");
      }
      $rd++;
    }
  }

  my @atoms;
  foreach my $chain (@{$pdb->{chain}}) {
    my $cp = $cd{$chain->{id}}; # chain permission (may be undefined)
    if ((!$cd) || ($pcd && $cp) || (!$pcd && (!defined($cp) || $cp))) {
      foreach my $res (@{$chain->{res}}) {
        my $rp = $rd{"$chain->{id}\_$res->{iresnum}"}; # residue permission (may be undefined)
        if ((!$rd) || ($prd && $rp) || (!$prd && (!defined($rp) || $rp))) {
          if (!defined($rp) || ($rp == 1)) {
            push(@atoms, @{$res->{atom}});
          } elsif ($rp == 2) {
            push(@atoms, PDB::sidechain($res));
          } elsif ($rp == 3) {
            push(@atoms, PDB::backbone($res));
          } else {
            GENERAL::error("Invalid residue permission $rp - internal function problem!");
          }
        }
      }
    }
  }
  if (defined($pflag) && $pflag) {
    return \@atoms;
  } else {
    return @atoms;
  }
}


=head2   FUNCTION -- closest

 Title   : closest
 Usage   : my $dist = PDB::closest($res1, $res2)
 Function: Returns the closest distance between atoms of the two
           given residues.
 Returns : Closest distance.
 Args    : 1. Residue 1
           2. Residue 2
           3. Optional - a double hash table with atomic radii (hashed
              first by residue name and second by atom name). If this
              parameter is specified, the distance between surfaces of
              atoms will be considered instead of the distance between
              atom centers.

=cut

sub closest {
  my $res1 = shift;
  my $res2 = shift;
  GENERAL::requireArgs($res1, $res2);
  my $R = shift;

  my $min = -1;
  foreach my $atom1 (@{$res1->{atom}}) {
    foreach my $atom2 (@{$res2->{atom}}) {
      my $d = (($atom1->{xcoor} - $atom2->{xcoor})**2 + ($atom1->{ycoor} - $atom2->{ycoor})**2 + ($atom1->{zcoor} - $atom2->{zcoor})**2)**0.5;
      if (defined($R)) {
        if (!defined($R->{$res1->{resname}}->{$atom1->{atomname}}) || !defined($R->{$res2->{resname}}->{$atom2->{atomname}})) {
          GENERAL::error("Undefined atom size: $atom1->{atomname} in $res1->{resname} or  $atom2->{atomname} in $res2->{resname} (or both)!");
        }
        $d -= $R->{$res1->{resname}}->{$atom1->{atomname}} + $R->{$res2->{resname}}->{$atom2->{atomname}};
      }
      if (($min < 0) || ($min > $d)) {
        $min = $d;
      }
    }
  }
  return $min;
}


=head2   FUNCTION -- atomsWithin

 Title   : atomsWithin
 Usage   : my @atoms = $pdb->atomsWithin(\@source_atoms, $distance)
 Function: Returns a list of atoms within the specified distance from a set
           of specified source atoms. The returned list naturally at
           least contains the list of source atoms.
 Returns : A list of atoms.
 Args    : 1. A pointer to an array of source atoms.
           2. Distance cutoff
           3. Optional - a double hash table with atomic radii (hashed
              first by residue name and second by atom name). If this
              parameter is not specified, this hash table is read from disk
              using PDB::sizeLookup(). If this parameter is specified as 0,
              distances will be measured between atom centers, not to the vdW surfaces.

=cut

sub atomsWithin {
  my $pdb = shift;
  my $atoms = shift;
  my $dcut = shift;
  GENERAL::requireArgs($pdb, $atoms, $dcut);
  if (scalar(@$atoms) == 0) { my @arr = (); return \@arr; }
  my $R = shift;
  if (!defined($R)) { $R = PDB::sizeLookup(); }
  my @watoms;

  # First, construct a box in which the source atoms are located
  my $fa = $atoms->[0];
  my $far = ($R eq 0) ? 0 : $R->{$fa->{residue}->{resname}}->{$fa->{atomname}};
  my ($xmin,$ymin,$zmin,$xmax,$ymax,$zmax) =
      ($fa->{xcoor}-$far,$fa->{ycoor}-$far,$fa->{zcoor}-$far,$fa->{xcoor}+$far,$fa->{ycoor}+$far,$fa->{zcoor}+$far);
  for (my $i=1; $i < scalar(@$atoms); $i++) {
    my $a = $atoms->[$i];
    my $ar = ($R eq 0) ? 0 : $R->{$a->{residue}->{resname}}->{$a->{atomname}};
    if ($a->{xcoor}-$ar < $xmin) { $xmin = $a->{xcoor}-$ar;}
    if ($a->{xcoor}+$ar > $xmax) { $xmax = $a->{xcoor}+$ar;}
    if ($a->{ycoor}-$ar < $ymin) { $ymin = $a->{ycoor}-$ar;}
    if ($a->{ycoor}+$ar > $ymax) { $ymax = $a->{ycoor}+$ar;}
    if ($a->{zcoor}-$ar < $zmin) { $zmin = $a->{zcoor}-$ar;}
    if ($a->{zcoor}+$ar > $zmax) { $zmax = $a->{zcoor}+$ar;}
  }
  # Don't bother checking closely any atoms outside of that box + distance cutoff
  foreach my $chain (@{$pdb->{chain}}) {
    foreach my $res (@{$chain->{res}}) {
      foreach my $a (@{$res->{atom}}) {
        my $ar = ($R eq 0) ? 0 : $R->{$a->{residue}->{resname}}->{$a->{atomname}};
        # If the atom is within the box, consider it
        if (($a->{xcoor}+$ar+$dcut >= $xmin) && ($a->{xcoor}-$ar-$dcut <= $xmax) &&
            ($a->{ycoor}+$ar+$dcut >= $ymin) && ($a->{ycoor}-$ar-$dcut <= $ymax) &&
            ($a->{zcoor}+$ar+$dcut >= $zmin) && ($a->{zcoor}-$ar-$dcut <= $zmax)) {
          foreach my $sa (@$atoms) {
            my $d = (($a->{xcoor} - $sa->{xcoor})**2 + ($a->{ycoor} - $sa->{ycoor})**2 + ($a->{zcoor} - $sa->{zcoor})**2)**0.5;
            my $dr = ($R eq 0) ? 0 : $R->{$res->{resname}}->{$a->{atomname}} + $R->{$sa->{residue}->{resname}}->{$sa->{atomname}};
            if ($d - $dr <= $dcut) {
              push(@watoms, $a);
              last;
            }
          }
        }
      }
    }
  }
  return \@watoms;
}


=head2   FUNCTION -- resWithin

 Title   : resWithin
 Usage   : my @res = $pdb->resWithin(\@source_atoms, $distance)
 Function: Returns a list of residues within the specified distance from a set
           of specified source atoms. The returned list naturally at
           least contains the list of residues the source atoms are in.
 Returns : A list of residues.
 Args    : 1. A pointer to an array of source atoms.
           2. Distance cutoff
           3. Optional - a double hash table with atomic radii (hashed
              first by residue name and second by atom name). If this
              parameter is specified, this hash table is read from disk
              using PDB::sizeLookup(). If this parameter is specified as 0,
              distances will be measured between atom centers, not to the vdW surfaces.

=cut

sub resWithin {
  my $pdb = shift;
  my $atoms = shift;
  my $dcut = shift;
  GENERAL::requireArgs($pdb, $atoms, $dcut);
  my $R = shift;
  if (!defined($R)) { $R = PDB::sizeLookup(); }
  my @wres;

  # First, construct a box in which the source atoms are located
  if (scalar(@$atoms) == 0) { GENERAL::error("Zero source atoms were specified!"); }
  my $fa = $atoms->[0];
  my $far = ($R eq 0 ? 0 : $R->{$fa->{residue}->{resname}}->{$fa->{atomname}});
  my ($xmin,$ymin,$zmin,$xmax,$ymax,$zmax) =
      ($fa->{xcoor}-$far,$fa->{ycoor}-$far,$fa->{zcoor}-$far,$fa->{xcoor}+$far,$fa->{ycoor}+$far,$fa->{zcoor}+$far);
  for (my $i=1; $i < scalar(@$atoms); $i++) {
    my $a = $atoms->[$i];
    my $ar = ($R eq 0 ? 0 : $R->{$a->{residue}->{resname}}->{$a->{atomname}});
    if ($a->{xcoor}-$ar < $xmin) { $xmin = $a->{xcoor}-$ar;}
    if ($a->{xcoor}+$ar > $xmax) { $xmax = $a->{xcoor}+$ar;}
    if ($a->{ycoor}-$ar < $ymin) { $ymin = $a->{ycoor}-$ar;}
    if ($a->{ycoor}+$ar > $ymax) { $ymax = $a->{ycoor}+$ar;}
    if ($a->{zcoor}-$ar < $zmin) { $zmin = $a->{zcoor}-$ar;}
    if ($a->{zcoor}+$ar > $zmax) { $zmax = $a->{zcoor}+$ar;}
  }
  # Don't bother checking closely any atoms outside of that box + distance cutoff
  foreach my $chain (@{$pdb->{chain}}) {
    foreach my $res (@{$chain->{res}}) {
      my $found = 0;
      foreach my $a (@{$res->{atom}}) {
        my $ar = ($R eq 0 ? 0 : $R->{$a->{residue}->{resname}}->{$a->{atomname}});
        # If the atom is within the box, consider it
        if (($a->{xcoor}+$ar+$dcut > $xmin) && ($a->{xcoor}-$ar-$dcut < $xmax) &&
            ($a->{ycoor}+$ar+$dcut > $ymin) && ($a->{ycoor}-$ar-$dcut < $ymax) &&
            ($a->{zcoor}+$ar+$dcut > $zmin) && ($a->{zcoor}-$ar-$dcut < $zmax)) {
          foreach my $sa (@$atoms) {
            my $d = (($a->{xcoor} - $sa->{xcoor})**2 + ($a->{ycoor} - $sa->{ycoor})**2 + ($a->{zcoor} - $sa->{zcoor})**2)**0.5;
            if ($d - ($R eq 0 ? 0 : $R->{$res->{resname}}->{$a->{atomname}}) - ($R eq 0 ? 0 : $R->{$sa->{residue}->{resname}}->{$sa->{atomname}}) < $dcut) {
              push(@wres, $a->{residue});
              $found = 1;
              last;
            }
          }
        }
        if ($found) { last; }
      }
    }
  }
  return \@wres;
}


=head2   FUNCTION -- box

 Title   : box
 Usage   : my @box = $pdb->box();
           my @box = PDB::box(\@atoms);
 Function: Builds a box around the given atoms and returns it.
 Returns : An array defining the box: (xmin, xmax, ymin, ymax, zmin, zmax).
 Args    : 1. A pointer to an array of source atoms or a PDB structure (all atoms will be considered).
           2. Optional - a double hash table with atomic radii (hashed
              first by residue name and second by atom name). If this
              parameter is not specified, this hash table is read from disk
              using PDB::sizeLookup(). If this parameter is specified as 0,
              the box boundaries will be defined by atom centers, not vdW surfaces.

=cut

sub box {
  my $pdb = shift;
  GENERAL::requireArgs($pdb);
  my $R = shift;
  if (!defined($R)) { $R = PDB::sizeLookup(); }

  my $atoms;
  if (ref($pdb) =~ /ARRAY/) { $atoms = $pdb; }
  elsif (ref($pdb) =~ /PDB/) { $atoms = $pdb->conAtoms(); }
  else { GENERAL::error("Unexpected type of reference (" . ref($pdb)); }

  if (scalar(@$atoms) == 0) { GENERAL::error("Zero source atoms were specified!"); }
  my $fa = $atoms->[0];
  my $far = ($R eq 0) ? 0 : $R->{$fa->{residue}->{resname}}->{$fa->{atomname}};
  my ($xmin,$ymin,$zmin,$xmax,$ymax,$zmax) =
      ($fa->{xcoor}-$far,$fa->{ycoor}-$far,$fa->{zcoor}-$far,$fa->{xcoor}+$far,$fa->{ycoor}+$far,$fa->{zcoor}+$far);
  for (my $i=1; $i < scalar(@$atoms); $i++) {
    my $a = $atoms->[$i];
    my $ar = ($R eq 0) ? 0 : $R->{$a->{residue}->{resname}}->{$a->{atomname}};
    if ($a->{xcoor}-$ar < $xmin) { $xmin = $a->{xcoor}-$ar;}
    if ($a->{xcoor}+$ar > $xmax) { $xmax = $a->{xcoor}+$ar;}
    if ($a->{ycoor}-$ar < $ymin) { $ymin = $a->{ycoor}-$ar;}
    if ($a->{ycoor}+$ar > $ymax) { $ymax = $a->{ycoor}+$ar;}
    if ($a->{zcoor}-$ar < $zmin) { $zmin = $a->{zcoor}-$ar;}
    if ($a->{zcoor}+$ar > $zmax) { $zmax = $a->{zcoor}+$ar;}
  }

  return ($xmin, $xmax, $ymin, $ymax, $zmin, $zmax);
}

=head2   FUNCTION -- atomDist

 Title   : atomDist
 Usage   : my $d = PDB::atomDist()
 Function: Calculates the Euclidian distance between two atoms.
 Returns : Distance (scalar).
 Args    : 1. atom 1 (reference)
           2. atom 2 (reference)
=cut

sub atomDist {
  my $atom1 = shift;
  my $atom2 = shift;

  return (($atom1->{xcoor}-$atom2->{xcoor})**2+($atom1->{ycoor}-$atom2->{ycoor})**2+($atom1->{zcoor}-$atom2->{zcoor})**2)**0.5;
}


=head2   FUNCTION -- sizeLookup

 Title   : sizeLookup
 Usage   : my $R = PDB::sizeLookup()
 Function: Creates a double hash table of VdW radii for atom types
           in residues (uses Charm19 size file). Radii are hashed
           by residue name first and atom type second.
 Returns : Hash table reference.
 Args    : none

=cut

sub sizeLookup {
  my $par = shift;
  if (!defined($par)) { $par = "19"; }
  my ($filename, $fh);

  # Open charge lookup file
  if ($par eq "19") {
    $filename = DEFINITIONS::getParmFile("charmm19.siz");
  } elsif ($par eq "22") {
    $filename = DEFINITIONS::getParmFile("charmm22.siz");
  } else {
    GENERAL::error("Parameter version \"$par\" unknown!");
  }
  $fh = GENERAL::GetInFH($filename);

  # Read the size lookup file into a hash
  my %R; # hash table for radii of atoms by residue
  my %r; # hash table for radii of atoms, for which the radius does not depend on residue
  my ($res, $atom, $radius);
  while (<$fh>) {
    chomp($_);
    s/\!$//; # remove comments
    next if (! /\s+/);
    if (/\s*(\S+)\s+(\S+)\s+(\S+)\s*/) {
      $atom = $1;
      $res = $2;
      $radius = $3;
      $R{$res}{$atom} = $radius;
    } elsif (/\s*(\S+)\s+(\S+)\s*/) {
      $atom = $1;
      $radius = $2;
      $r{$atom} = $radius;
    }
  }
  $fh->close();

  # Each residue should have an entry for every residue independent atom radius
  foreach $res (keys(%R)) {
    foreach $atom (keys(%r)) {
      $R{$res}{$atom} = $r{$atom};
    }
  }

  return \%R;
}

=head2   FUNCTION -- parsePSF

 Title   : parsePSF
 Usage   : $pdb->parsePSF($psffile)
 Function: Parses the a CHARMM PSF file and populates the specified PDB structure with various information.
           Currently, only reads in the charges and masses.
 Returns :
 Args    : 1. PDB structure.
           2. Name of PSF file

=cut
sub parsePSF {
  my $pdb = shift;
  my $psff = shift;
  my $mrk = "_123_mark";

  my $ifh = GENERAL::GetInFH($psff);
  while (<$ifh>) { last if (/\!NATOM/); }
  while (<$ifh>) {
    last if (/\!NBOND/);
    $_ = GENERAL::Trim($_); next if ($_ eq "");
    my @line = split(" ", $_);
    my $res = $pdb->getResByInd($line[1], $line[2]);
    GENERAL::assert($res->{resname} eq $line[3], "Residue in chain $line[1], with iresnum $line[2] is called $res->{resname} in the PDB object and $line[3] in the PSF ($psff)");
    my $a = PDB::getAtomInRes($res, $line[4]);
    $a->{atype} = $line[5];
    $a->{charge} = $line[6];
    $a->{mass} = $line[7];
    $a->{$mrk} = 1;
  }
  close($ifh);

  # check if all atoms were accounted for
  foreach my $chain (@{$pdb->{chain}}) {
    foreach my $res (@{$chain->{res}}) {
      foreach my $at (@{$res->{atom}}) {
        GENERAL::error("Not all atoms were found in PSF (e.g. $at->{atomname} in residue " . PDB::resStr($res) . " was not found)") if (!defined($at->{$mrk}));
        delete($at->{$mrk});
      }
    }
  }
}

=head2   FUNCTION -- charge

 Title   : charge
 Usage   : $pdb->charge()
           PDB::charge(\@atoms)
 Function: Charges the given PDB structure. Assigns charges
           to atoms based on the specified double hash table.
           If no table is passed, creates one with PDB::chargeLookup.
 Returns :
 Args    : 1. Either a PDB structure or a reference to an array of atoms.
           2. Optional: residue specific atomic charge hash
           3. Optional: parameter version (19 or 22)

=cut
sub charge {
  my $struct = shift;
  my $C = shift;
  my $par = shift;
  if (!defined($C)) {
    $C = PDB::chargeLookup($par);
  }

  if (ref($struct) =~ /PDB/) {
    foreach my $c (@{$struct->{chain}}) {
      foreach my $res (@{$c->{res}}) {
        foreach my $atom (@{$res->{atom}}) {
          if (defined($C->{$res->{resname}}->{$atom->{atomname}})) {
            $atom->{charge} = $C->{$res->{resname}}->{$atom->{atomname}};
          } else {
            $atom->{charge} = 0.0;
            GENERAL::warning("Could not find charge for atom " . $atom->{atomname} . " in residue " . $res->{resname} . "!");
          }
        }
      }
    }
  } elsif (ref($struct) =~ /ARRAY/) {
    foreach my $atom (@$struct) {
      my $res = $atom->{residue};
      if (defined($C->{$res->{resname}}->{$atom->{atomname}})) {
        $atom->{charge} = $C->{$res->{resname}}->{$atom->{atomname}};
      } else {
        $atom->{charge} = 0.0;
        GENERAL::warning("Could not find charge for atom " . $atom->{atomname} . " in residue " . $res->{resname} . "!");
      }
    }
  } else {
    GENERAL::error("Unexpected structure type " . ref($struct) . " (structure $struct)!");
  }
}

=head2   FUNCTION -- neutralize

 Title   : charge
 Usage   : $pdb->neutralize()
           PDB::neutralize(\@atoms)
 Function: Zeros down the atomic charges in the given PDB structure.
 Returns :
 Args    : 1. Either a PDB structure or a reference to an array of atoms.

=cut
sub neutralize {
  my $struct = shift;

  if (ref($struct) =~ /PDB/) {
    foreach my $c (@{$struct->{chain}}) {
      foreach my $res (@{$c->{res}}) {
        foreach my $atom (@{$res->{atom}}) {
          $atom->{charge} = 0.0;
        }
      }
    }
  } elsif (ref($struct) =~ /ARRAY/) {
    foreach my $atom (@$struct) {
      $atom->{charge} = 0.0;
    }
  } else {
    GENERAL::error("Unexpected structure type " . ref($struct) . " (structure $struct)!");
  }
}

=head2   FUNCTION -- chargeLookup

 Title   : chargeLookup
 Usage   : my $C = PDB::chargeLookup()
 Function: Creates a double hash table of charges for atom types
           in residues (uses Charm19 charge file). Charges are
           hashed by residue name first and atom type second.
 Returns : Hash table reference.
 Args    : none

=cut

sub chargeLookup {
  my $par = shift;
  if (!defined($par)) { $par = "19"; }
  my ($filename, $fh);

  # Open charge lookup file
  if ($par eq "19") {
    $filename = DEFINITIONS::getParmFile("Charmm19_charge");
  } elsif ($par eq "22") {
    $filename = DEFINITIONS::getParmFile("Charmm22_charge");
  } else {
    GENERAL::error("Parameter version \"$par\" unknown!");
  }
  $fh = GENERAL::GetInFH($filename);

  # Read the charge lookup file into a hash
  my %C; # hash table for charges
  my ($res, $atom, $charge);
  while (<$fh>) {
    if (/RESI\s+(.+)/) {
      chomp($res = $1);
      $C{$res} = {};
      $C{$res}->{'total'} = 0;
      $C{$res}->{'charged'} = 0;
    }
    if (/ATOM\s*(.+)\s+(.+)/) {
      $atom = $1;
      $charge = $2;
      $atom =~ s/\s+$//;
      $charge =~ s/\s+$//;
      if (defined($C{$res}->{$atom})) {
        GENERAL::error("The same atom type ($atom), defined twice in residue ($res)");
      } else {
        $C{$res}->{$atom} = $charge;
        $C{$res}->{'total'} += $charge;
        if ($charge != 0) {$C{$res}->{'charged'} += 1;}
      }
    }
  }
  close($fh);

  return \%C;
}


=head2   FUNCTION -- isRenumbered

 Title   : isRenumbered
 Usage   : my $yes_no = PDB::isRenumbered("test.pdb")
           OR
           my $yes_no = $pdb->isRenumbered();
 Function: Determines if the PDB (either object or file) is renumbered.
 Returns : 0 or 1 accordingly.
 Args    : PDB file name or PDB object.

=cut

sub isRenumbered {
  my $inp = shift;
  my $pdb;
  if (ref($inp) =~ /PDB/) {
    $pdb = $inp;
  } else {
    $pdb = PDB::new($inp);
  }
  
  foreach my $chain (@{$pdb->{chain}}) {
    my $c = 1;
    foreach my $res (@{$chain->{res}}) {
      if ($res->{resnum} ne $c) { return 0; }
      $c++;
    }
  }
  
  return 1;
}

=head2   FUNCTION -- Remap2Crystal

 Title   :  Remap2Crystal
 Usage   :  $mol->Remap2Crystal($resmap, 'ABC')
 Function:  applies the reverse translation from the one
            given in the map argument, chain id strings argument
            is optional
 Example :  during the design, we did every thing on renumbered pdb file
            however, when we get the solution, we want to translate it back to
            original pdb residue naming so we can compare with literature
            discussion.
            here is the procedure how we can do that:

            $crystal = PDB::new('pdbblah.ent');
            $resmap = $crystal->getResMap();  # know how to map between two structure
            $mol= PDB::new('test.crd');     # solution crd file
            $mol->Remap2Crystal($resmap);   # renumber the residue numbering
            $mol->writePDB("blah.pdb", 'generic_noh')  # write out renumbered pdb files

            # can't be more complicated!

 Returns :  none
 Args    :  reference to residue map and chain id strings

=cut


sub Remap2Crystal {
    my $self=shift;
    my $map=shift;
    my $chainstr=shift;

    if (! defined $chainstr) {
    $self->selectAllChains();
    } else {
    $self->selectChains($chainstr);
    }

  foreach my $c (@{$self->{pickedchains}}) {
    my $chainid=$c->{id};
    foreach my $res ($c->{res}) {
# Gevorg: under the new PDB structure, this part becomes uncecessary. Only
# residues need to know that their number has changed, since everyone else (atoms
# or the chain) get this information through the residues.
#    my $res=$c->{res};
#        foreach my $atom (@{$res}) {
#            my $atom=$c->{atom};
#        for (my $i=0; $i<=$#{$atom}; $i++) {
#            my $newnum=$map->{$chainid}->{'tocrystal'}->{$atom->[$i]->{iresnum}};
#            $atom->[$i]->{resnum}=$newnum;
#        }

      my $newnum = $map->{$chainid}->{'tocrystal'}->{$res->{iresnum}};
      $res->{num}=$newnum;
    }
  }
}




=head2   FUNCTION -- writeInterfaceCoord

 Title   :  writeInterfaceCoord
 Usage   :  $mol->writeInterfaceCoord($reference_to_interface_res, "interface.pdb", $resmap)
 Function:  given a reference to interface residue list, write out
            pdb file for interface residues
 Example :  problem defined: we want to translate the interface residue file back to
            crystal residue naming.

            my $surface = NACCESS::new();
            my $refIntRes = $surface->getInterface($Crystal_complex_pdb, $molecule1_pdb, $molecule2_pdb)

            my $mol = PDB::new($SolutionCrd);
            $resmap = $mol->getResMap();

            my $newmol = $mol->clone(1);     # let's work with another copy of molecule object
                                             # you don't want to change the original data
            $newmol->writeInterface($refIntRes,"interface.pdb", $resmap);

            see /usr/programs/ProteinDesign/perllib/test/PDB/Interface/test.pl

 Returns :
 Args    :  (1) reference to interface residues
            (2) output pdb file for the interface.
            (3) if reference to residue naming map is given, then renumber the pdb to the style
                in crystal pdb file.

=cut


sub writeInterfaceCoord {

    my $self=shift;
    my $rotamer=shift;
    my $fname=&GENERAL::GetOutFH(shift);
    my $resmap=shift;

    $self->selectResidues($rotamer);
    my $newmol = $self->clone(1);
    if (! defined $resmap) {
    $newmol->writePDB($fname, 'genericnoh');
    } else {
    $newmol->Remap2Crystal($resmap);
    $newmol->writePDB($fname, 'genericnoh');
    }

}


sub getDesignMol {

    my $self=shift;
    my $rotamer=shift;

    $self->selectResidues($rotamer);
    my $newmol = $self->clone(1);
    return $newmol;
}

sub writeInterfaceFragments {

    my $self=shift;
    my $fraglist=shift;
    my $resmap=shift;

    if (! defined $fraglist) {
    print "the fragment list at the interface is not defined.\n";
    exit;
    }

    &GENERAL::MakeDir("SecStr");
    chdir("SecStr");

    foreach my $fragment ( @$fraglist ) {
       my $fname=$fragment->{'label'}.".pdb";
       my $reslist=[];
       for (my $j=$fragment->{start}; $j <= $fragment->{end}; $j++) {
            my $res={};
            $res->{chain}= $fragment->{chain};
            $res->{iresnum}=$j;
        push (@$reslist, $res);
        }
       $self->selectResidues($reslist);
       my $newmol = $self->clone(1);
       if (! defined $resmap) {
            $newmol->writePDB($fname, 'genericnohnot');
       } else {
            $newmol->Remap2Crystal($resmap);
            $newmol->writePDB($fname, 'genericnohnot');
       }
       system("cat $fname >> total.pdb");

   }
}



=head2   FUNCTION -- saveCoordinates

 Title   : saveCoordinates
 Usage   : $pdb->saveCoordinates() or $pdb->saveCoordinates("saved1")
 Function: Saves the current coordinates contained in the PDB structure into
           temporary fileds within atom entries. This way the coordinates can
           be modified and later restored, if needed, via restorCoordinates().
 Returns : nothing
 Args    : 1. PDB object
           2. optional: name of the save coordinate set. If not set, will use default name.

=cut

sub saveCoordinates {
  my $self = shift;
  my $name = shift;
  $name = "saved_default" if (!defined($name));

  foreach my $ch (@{$self->{chain}}) {
    foreach my $res (@{$ch->{res}}) {
      foreach my $a (@{$res->{atom}}) {
        $a->{savedCoords}->{$name}->{xcoor} = $a->{xcoor};
        $a->{savedCoords}->{$name}->{ycoor} = $a->{ycoor};
        $a->{savedCoords}->{$name}->{zcoor} = $a->{zcoor};
      }
    }
  }
}

=head2   FUNCTION -- restoreCoordinates

 Title   : restoreCoordinates
 Usage   : $pdb->restoreCoordinates() or $pdb->restoreCoordinates("saved1")
 Function: Restrores previously saved coordinates by copying them on top of the current ones.
 Returns : nothing
 Args    : 1. PDB object
           2. optional: name of the saved coordinate set. If not set, will use default name.

=cut

sub restoreCoordinates {
  my $self = shift;
  my $name = shift;
  $name = "saved_default" if (!defined($name));

  foreach my $ch (@{$self->{chain}}) {
    foreach my $res (@{$ch->{res}}) {
      foreach my $a (@{$res->{atom}}) {
        $a->{xcoor} = $a->{savedCoords}->{$name}->{xcoor};
        $a->{ycoor} = $a->{savedCoords}->{$name}->{ycoor};
        $a->{zcoor} = $a->{savedCoords}->{$name}->{zcoor};
      }
    }
  }
}


=head2   FUNCTION -- clone

 Title   : clone
 Usage   : $newmol = $mol->clone(1)  or $newmol = $mol->clone(1, "AB")
 Function: generates a copy of the current molecule
           and returns it. If the valid flag is set, only
           valid residues are copied. a chain string is optional argument
 Returns : new PDB structure object
 Args    : 1. Optional: valid_flag. If set to 0, everything will be copied (default).
              If set to 1, only the valid residues will be copied. If set to
              0.5, only the valid residues will be copied and in those only
              the valid atoms will be copied.
           2. Optional: chain string: 'AB', 'A', 'ABC' etc....

=cut

sub clone {
  my $self = shift;
  GENERAL::requireArgs($self);

  my $valid = shift;
  if (!defined($valid)) { $valid = 0; }
  my $chainstr = shift;

  my $n = {};
  $n->{chain} = ();
  $n->{chainlookup} = {};
  $n->{pickedchains} = undef;
  $n->{coordfile} = $self->{coordfile};
  bless $n;

  if (!defined $chainstr) {
    $self->selectAllChains();
  } else {
    $self->selectChains($chainstr);
  }

  foreach my $c (@{$self->{pickedchains}}) {
    # are there any valid residues in this chain?
    my $valid_in_chain = 0; my $nc;
    foreach my $r (@{$c->{res}}) {
      if (!$valid || (defined($r->{valid}) && $r->{valid})) {
        if ($valid_in_chain == 0) {
          $nc = $n->_newChain($c->{id});
          $valid_in_chain = 1;
        }
        my $nr;
        if ($valid == 0.5) {
          $nr = PDB::cloneResidue($r, 1);
        } else {
          $nr = PDB::cloneResidue($r);
        }
        push(@{$nc->{res}}, $nr);
        $nr->{chain} = $nc;
      }
    }
  }
  return $n;
}


=head2   FUNCTION -- cloneResidue

 Title   : cloneResidue
 Usage   : $newres = $pdb->cloneResidue($res)
 Function: Generates a copy of the specified residue.
 Returns : A copy of the specified residue.
 Args    : 1. Residue reference to copy.
           2. Optional: valid fralg. If set, only valid atoms will be copied.
              The default is to ignore valid flag.

=cut

sub cloneResidue {
  my $res = shift;
  GENERAL::requireArgs($res);
  my $valid = shift;
  if (!defined($valid)) { $valid = 0; }

  my $n = {};
  # copy all non-reference values
#   foreach my $k (keys(%$res)) {
#     if (!ref($res->{$k})) {
#       $n->{$k} = $res->{$k};
#     }
#   }
  $n->{resname} = $res->{resname};
  $n->{resnum} = $res->{resnum};
  $n->{iresnum} = $res->{iresnum};
  $n->{uresnum} = $res->{uresnum};
  $n->{valid} = $res->{valid};
  $n->{chain} = $res->{chain};
  $n->{atom} = ();
  foreach my $atom (@{$res->{atom}}) {
    if (!$valid || (defined($atom->{valid}) && $atom->{valid}) || (!defined($atom->{valid}))) {
      my $na = PDB::cloneAtom($atom);
      push(@{$n->{atom}}, $na);
      $na->{residue} = $n;
    }
  }

  return $n;
}


=head2   FUNCTION -- copyChain

 Title   : copyRChain
 Usage   : PDB::copyChain($source_chain, $desination_chain)
 Function: Copies the contents of the source chain into the destination
           chain. The information in the old destination chain is lost.
 Returns : Nothing
 Args    : 1. Source chain
           2. Destination chain
           3. Optional: specifies whether the atoms of the source chain are
              to be moved (1) or copied (0) to the new chain. Default is copied.
=cut

sub copyChain {
  my $sch = shift;
  my $dch = shift;
  my $move = shift;
  if (!defined($move)) { $move = 0; }

  $dch->{res} = ();
  foreach my $res (@{$sch->{res}}) {
    my $nr = PDB::_newResidue($res->{resname}, $res->{resnum}, $res->{iresnum}, $res->{uresnum}, $dch);
    PDB::copyResidue($res, $nr, $move);
  }
}


=head2   FUNCTION -- copyResidue

 Title   : copyResidue
 Usage   : PDB::copyResidue($source_res, $desination_res)
 Function: Copies the contents of the source residue into the destination
           residue. The information in the old destination residue is lost.
 Returns : Nothing
 Args    : 1. Source residue
           2. Destination residue
           3. Optional: specifies whether the atoms of the source residue are
              to be moved (1) or copied (0) to the new residue. Default is copied.
=cut

sub copyResidue {
  my $sres = shift;
  my $dres = shift;
  my $move = shift;
  if (!defined($move)) { $move = 0; }

  $dres->{resname} = $sres->{resname};
  $dres->{resnum} = $sres->{resnum};
  # iresnum and uresnum do not change
#  $dres->{iresnum} = $sres->{iresnum};
#  $dres->{uresnum} = $sres->{uresnum};
  if (defined($sres->{valid})) {
    $dres->{valid} = $sres->{valid};
  }
  if (defined($sres->{charge})) {
    $dres->{charge} = $sres->{charge};
  }
#  $dres->{chain} = $sres->{chain};
  undef(@{$dres->{atom}});
  $dres->{atom} = ();
  foreach my $atom (@{$sres->{atom}}) {
    my $na;
    if ($move) {
      $na = $atom;
    } else {
      $na = PDB::cloneAtom($atom);
    }
    push(@{$dres->{atom}}, $na);
    $na->{residue} = $dres;
  }

  # since messed with atoms, set the altered flag
  $dres->{altered} = 1;
}


=head2   FUNCTION -- copyAtom

 Title   : copyAtom
 Usage   : PDB::copyAtom($source_atom, $desination_atom)
 Function: Copies the contents of the source atom into the destination
           atom. The information in the old destination atom is lost.
 Returns : Nothing
 Args    : 1. Source residue
           2. Destination residue
=cut

sub copyAtom {
  my $satom = shift;
  my $datom = shift;

  $datom->{atominx} = $satom->{atominx};
  $datom->{atomunx} = $satom->{atomunx};
  $datom->{atomname} = $satom->{atomname};
  $datom->{xcoor} = $satom->{xcoor};
  $datom->{ycoor} = $satom->{ycoor};
  $datom->{zcoor} = $satom->{zcoor};
  $datom->{B} = $satom->{B};
  $datom->{occ} = $satom->{occ};
  $datom->{hyd} = $satom->{hyd};
  $datom->{charge} = $satom->{charge};
  $datom->{valid} = $satom->{valid} if (defined($satom->{valid}));

  # since messed with atoms, set the altered flag
  $datom->{altered} = 1;
}

=head2   FUNCTION -- replaceSidechain

 Title   : replaceSidechain
 Usage   : PDB::replaceSidechain($source_res, $desination_res)
 Function: Changes the identity of the destination residue by copying over only the sidechain
           atoms of the source residue. Properties of the destination residue are also changed
           to those of the source residue (name, charge).
 Returns : Nothing
 Args    : 1. Source residue
           2. Destination residue.
           3. Optional: specifies whether the atoms of the source residue are
              to be moved (1) or copied (0) to the new residue. Default is copied.
=cut

sub replaceSidechain {
  my $sres = shift;
  my $dres = shift;
  my $move = shift;
  if (!defined($move)) { $move = 0; }

  $dres->{resname} = $sres->{resname};
  if (defined($sres->{charge})) {
    $dres->{charge} = $sres->{charge};
  }
  my @ba = PDB::backbone($dres);
  $dres->{atom} = \@ba;
  my @sa = PDB::sidechain($sres);
  foreach my $atom (@sa) {
    my $na;
    if ($move) {
      $na = $atom;
    } else {
      $na = PDB::cloneAtom($atom);
    }
    push(@{$dres->{atom}}, $na);
    $na->{residue} = $dres;
  }

  # since messed with atoms, set the altered flag
  $dres->{altered} = 1;
}


sub cloneAtom {
  my $atom = shift;

  my $n = {};
  $n->{atominx} = $atom->{atominx};
  $n->{atomunx} = $atom->{atomunx};
  $n->{atomname} = $atom->{atomname};
  $n->{residue} = $atom->{residue};
  $n->{xcoor} = $atom->{xcoor};
  $n->{ycoor} = $atom->{ycoor};
  $n->{zcoor} = $atom->{zcoor};
  $n->{B} = $atom->{B};
  $n->{occ} = $atom->{occ};
  $n->{hyd} = $atom->{hyd};
  $n->{charge} = $atom->{charge};
  $n->{hetero} = $atom->{hetero};
  if (defined($atom->{valid})) {
    $n->{valid} = $atom->{valid};
  }

  return $n;
}


# Frees up the memory associated with the given residue (and all
# atoms that belong to it).
sub freeResidue {
  my $res = shift;
  for (my $i=0; $i < scalar(@{$res->{atom}}); $i++) {
    undef(%{$res->{atom}->[$i]});
  }
  undef(@{$res->{atom}});
  undef($res->{atom});
  undef(%$res);
  return;
}


=head2   FUNCTION -- renumberResidues

 Title   : renumberResidues
 Usage   : $pdb->renumberResidues()
 Function: Renumbers residues in all chains such that iresnum and resnum start with 1 in
           each chain and grow sequentially.
 Returns : nothing
 Args    : 1. PDB structure
           2. Flag - renumber internal numbers as well (resnum)?

=cut

sub renumberResidues {
  my $self = shift;
  my $int_num = shift;

  foreach my $chain (@{$self->{chain}}) {
    my $iresnum = 1;
    foreach my $res (@{$chain->{res}}) {
      $res->{iresnum} = $iresnum;
      $res->{resnum} = $iresnum if (defined($int_num) && $int_num);
      $iresnum += 1;
    }
  }
}


=head2   FUNCTION -- renumber

 Title   : renumber
 Usage   : $pdb->renumber()
 Function: Renumbers residues and atoms in all chains such that iresnum and resnum start with 1 in
           each chain and grow sequentially and atom numbers (atominx and atomunx) start with 1 and
           grow sequentially throughout the structure.
 Returns : nothing
 Args    : 1. PDB structure
           2. Flag - renumber internal numbers as well (resnum and atominx)?

=cut

sub renumber {
  my $self = shift;
  my $int_num = shift;

  my $atomi = 1;
  foreach my $chain (@{$self->{chain}}) {
    my $iresnum = 1;
    foreach my $res (@{$chain->{res}}) {
      $res->{iresnum} = $iresnum;
      $res->{resnum} = $iresnum if (defined($int_num) && $int_num);
      $iresnum += 1;
      foreach my $a (@{$res->{atom}}) {
        $a->{atomunx} = $atomi;
        $a->{atominx} = $atomi if (defined($int_num) && $int_num);
        $atomi++;
      }
    }
  }
}


=head2 FUNCTION -- sidechain

 Title   : sidechain
 Usage   : my @sidechain_atoms = PDB::sidechain($res)
 Function: Returns an array of all side chain atoms of the specified residue.
 Returns : An array of atoms.
 Args    : 1. Reference to residue
           2. Flag specifying whether to return a reference or an array.
=cut

sub sidechain {
  my $res = shift;
  my $ptr = shift;
  $ptr = 0 if (!defined($ptr));
  my @sa;
  my $bbn = PDB::backboneA('REGEXP');
  foreach my $atom (@{$res->{atom}}) {
    if ($atom->{atomname} !~ /$bbn/) {
      push(@sa, $atom);
    }
  }
  if ($ptr) { return \@sa; }
  else { return @sa; }
}


=head2 FUNCTION -- backbone

 Title   : backbone
 Usage   : my @backbone_atoms = PDB::backbone($res)
 Function: Returns an array of all backbone atoms of the specified residue.
 Returns : An array of atoms.
 Args    : 1. Reference to residue
           2. Flag specifying whether to return a reference or an array.
=cut

sub backbone {
  my $res = shift;
  my $ptr = shift;
  $ptr = 0 if (!defined($ptr));
  my @ba;
  my $bbn = PDB::backboneA('REGEXP');
  foreach my $atom (@{$res->{atom}}) {
    if ($atom->{atomname} =~ /$bbn/) {
      push(@ba, $atom);
    }
  }
  if ($ptr) { return \@ba; }
  else { return @ba; }
}


=head2 FUNCTION -- selectResidues

 Title   :  selectResidues
 Usage   :  selectResidues(reslist)
 Function:  sets the valid flag to 1
            for residues in the fragment list and to
            0 for residues outside the list
 Example :

 Returns :
 Args    :  An array of reference to residues

=cut

sub selectResidues {

    my $self=shift;
    my $rotamer=shift;
    my $reslist=$rotamer->{dslist};

    $self->resetValidResidues(0);

    foreach my $l ( @{$reslist} ) {
      my $c=$self->getChain($l->{chain});
      my $i=$l->{iresnum}-1;
      my $r=$c->{res}->[$i];
      $r->{valid}=1 if (defined $r);
    }
}


=head2   FUNCTION -- resetValidResidues

 Title   : resetValidResidues
 Usage   : resetValidResidues(1) or resetValidResidues(0)
 Function: sets the valid flag of all residues to the given value (default: 1)
 Returns :
 Args    : 1 or 0

=cut

sub resetValidResidues {
  my $self = shift;
  my $value = shift;

  $value = 1 if (!defined $value);

  $self->selectAllChains();
  foreach my $c ( @{$self->{pickedchains}} ) {
    foreach my $r ( @{$c->{res}} ) {
      $r->{valid} = $value;
      foreach my $a (@{$r->{atom}}) {
        $a->{valid} = $value;
      }
    }
  }
}



=head2  FUNCTION -- centerOfMass

 Title   :  centerOfMass
 Usage   :  my ($cx, $cy, $z) = $mol->centerOfMass() or  $mol->centerOfMass("AB")
 Function:  calculates the center of mass for the current structure. if chain id string
            is given, the center of mass of selected chains will be returned
 Returns :  x,y,z value of center of mass
 Args    :  null or chainid string

=cut

sub centerOfMass {
  my $struct = shift;
  my $chainstr = shift;

  my $N = 0; my $cx = 0; my $cy = 0; my $cz = 0;
  if (ref($struct) =~ /PDB/) {
    if (!defined $chainstr) {
      $struct->selectAllChains();
    } else {
      $struct->selectChains($chainstr);
    }
    foreach my $c ( @{$struct->{pickedchains}} ) {
      foreach my $r (@{$c->{res}}) {
        my $atom = $r->{atom};
        $N += scalar(@$atom);

        for (my $i = 0; $i < scalar(@$atom); $i++) {
          $cx += $atom->[$i]->{xcoor};
          $cy += $atom->[$i]->{ycoor};
          $cz += $atom->[$i]->{zcoor};
        }
      }
    }
  } elsif (ref($struct) =~ /ARRAY/) {
    foreach my $a (@$struct) {
      $N++;
      $cx += $a->{xcoor};
      $cy += $a->{ycoor};
      $cz += $a->{zcoor};
    }
  } else {
    GENERAL::error("Unexpected structure type (reference to " . ref($struct) . ", structure $struct)");
  }

  $cx/=$N;
  $cy/=$N;
  $cz/=$N;

  return ($cx,$cy,$cz);
}


=head2   FUNCTION -- center

 Title   :  center
 Usage   :  $mol->center() or  $mol->center("ABC")
 Function:  centers the current structure by subtracting the center
            of mass from all coordinates
 Returns :
 Args    :  (1) optional: chain selection string (by default all). These chains will be used
                for computing the center of mass

=cut

sub center {
  my $self = shift;
  my $chainstr = shift;
  my $ori = shift;

  if (! defined $chainstr) {
    $self->selectAllChains();
  } else {
    $self->selectChains($chainstr);
  }

  my ($cx,$cy,$cz);
  if (defined($ori)) {
    GENERAL::assert(scalar(@$ori) == 3, "if center specified, must contain 3 values!");
    $cx = $ori->[0]; $cy = $ori->[1]; $cz = $ori->[2];
  } else {
    ($cx,$cy,$cz) = $self->centerOfMass($chainstr);
  }

  foreach my $c (@{$self->{chain}}) {
    foreach my $res (@{$c->{res}}) {
      foreach my $atom (@{$res->{atom}}) {
        $atom->{xcoor} -= $cx;
        $atom->{ycoor} -= $cy;
        $atom->{zcoor} -= $cz;
      }
    }
  }
}


=head2   FUNCTION

 Title   :  rotateAroundZ
 Usage   :  PDB::rotateAroundZ(\@atoms, $angle)
 Function:  rotates the given list of atoms around Z by the given angle (in degrees)
 Returns :
 Args    :  (1) list of atoms to rotate aorund Z
            (2) angle to rotate by

=cut

sub rotateAroundZ {
  my $atoms = shift;
  my $an = shift;
  $an *= atan(1)/45;

  my $M = zeroes(3, 3);
  $M->slice("0,0:2") += pdl([cos($an)], [sin($an)], [0]);
  $M->slice("1,0:2") += pdl([-sin($an)], [cos($an)], [0]);
  $M->slice("2,0:2") += pdl([0], [0], [1]);

  foreach my $a (@$atoms) {
    my $coor = zeroes(1, 3);
    $coor->slice("0,0:2") += pdl([$a->{xcoor}], [$a->{ycoor}], [$a->{zcoor}]);
    $coor = $M x $coor;
    $a->{xcoor} = $coor->at(0, 0); $a->{ycoor} = $coor->at(0, 1); $a->{zcoor} = $coor->at(0, 2);
  }
}


=head2   FUNCTION -- matAlignVectorWithZAxis

 Title   :  matAlignVectorWithZAxis
 Usage   :  my $M = PDB::matAlignVectorWithZAxis($u, $v, $w)
 Function:  Returns the transformation matrix for aligning the Z axis onto the given axis
 Returns :
 Args    :  (1-3) the three components specifying the axis of interest

=cut

sub matAlignVectorWithZAxis {
  my $u = shift;
  my $v = shift;
  my $w = shift;
  my $uv = sqrt($u**2 + $v**2);
  my $uvw = sqrt($u**2 + $v**2 + $w**2);
  if ($uv == 0) {
    return pdl([1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]);
  }

  my $M1 = pdl([$u/$uv,  $v/$uv, 0,         0], [-$v/$uv, $u/$uv, 0, 0], [0,        0, 1,       0], [0, 0, 0, 1]);
  my $M2 = pdl([$w/$uvw, 0,      -$uv/$uvw, 0], [0,       1,      0, 0], [$uv/$uvw, 0, $w/$uvw, 0], [0, 0, 0, 1]);
  return $M2 x $M1;
}


=head2   FUNCTION -- matAlignVectorWithYAxis

 Title   :  matAlignVectorWithYAxis
 Usage   :  my $M = PDB::matAlignVectorWithYAxis($u, $v, $w)
 Function:  Returns the transformation matrix for aligning the Y axis onto the given axis
 Returns :
 Args    :  (1-3) the three components specifying the axis of interest

=cut

sub matAlignVectorWithYAxis {
  my $u = shift;
  my $v = shift;
  my $w = shift;
  my $uw = sqrt($u**2 + $w**2);
  my $uvw = sqrt($u**2 + $v**2 + $w**2);
  if ($uw == 0) {
    return pdl([1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]);
  }

  my $M1 = pdl([$w/$uw, 0, -$u/$uw, 0],  [0, 1, 0, 0], [$u/$uw, 0, $w/$uw, 0], [0, 0, 0, 1]);
  my $M2 = pdl([1, 0, 0, 0], [0, $v/$uvw, $uw/$uvw, 0], [0, -$uw/$uvw, $v/$uvw, 0], [0, 0, 0, 1]);
  return $M2 x $M1;
}

=head2   FUNCTION -- matAlignVectorWithXAxis

 Title   :  matAlignVectorWithXAxis
 Usage   :  my $M = PDB::matAlignVectorWithXAxis($u, $v, $w)
 Function:  Returns the transformation matrix for aligning the X axis onto the given axis
 Returns :
 Args    :  (1-3) the three components specifying the axis of interest

=cut

sub matAlignVectorWithXAxis {
  my $u = shift;
  my $v = shift;
  my $w = shift;
  my $vw = sqrt($v**2 + $w**2);
  my $uvw = sqrt($u**2 + $v**2 + $w**2);
  if ($vw == 0) {
    return pdl([1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]);
  }

# this expression appears to be wrong! Fix later, for now use the hack below.
#  my $M1 = pdl([1, 0, 0, 0], [0, $w/$vw, -$v/$vw, 0],  [0, $v/$vw, $w/$vw, 0], [0, 0, 0, 1]);
#  my $M2 = pdl([$vw/$uvw, 0, -$u/$uvw, 0], [0, 1, 0, 0], [$u/$uvw, 0, $vw/$uvw, 0], [0, 0, 0, 1]);
#  return $M2 x $M1;

  return matRotY(2*atan(1)) x matAlignVectorWithZAxis($u, $v, $w);
}


=head2   FUNCTION -- matRotX

 Title   :  matRotX
 Usage   :  my $M = PDB::matRotX(3.14159)
 Function:
 Returns :  Rotation matrix around X by the given angle
 Args    :  1) angle

=cut

sub matRotX {
  my $a = shift;
  return pdl([1, 0, 0, 0], [0, cos($a), -sin($a), 0], [0, sin($a), cos($a), 0], [0, 0, 0, 1]);
}


=head2   FUNCTION -- matRotY

 Title   :  matRotY
 Usage   :  my $M = PDB::matRotY(3.14159)
 Function:            
 Returns :  Rotation matrix around Y by the given angle
 Args    :  1) angle

=cut

sub matRotY {
  my $a = shift;
  return pdl([cos($a), 0, sin($a), 0], [0, 1, 0, 0], [-sin($a), 0, cos($a), 0], [0, 0, 0, 1]);
}


=head2   FUNCTION -- matRotZ

 Title   :  matRotZ
 Usage   :  my $M = PDB::matRotZ(3.14159)
 Function:            
 Returns :  Rotation matrix around Z by the given angle
 Args    :  1) angle

=cut

sub matRotZ {
  my $a = shift;
  return pdl([cos($a), -sin($a), 0, 0], [sin($a), cos($a), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]);
}


=head2   FUNCTION -- matRotOrigAxis

 Title   :  matRotOrigAxis
 Usage   :  my $M = PDB::matRotOrigAxis($u, $v, $w, $angle)
 Function:  Returns the matrix for rotating around vector (u,v,w) (e.g. axis formed by points (0, 0, 0) and (u, v, w)) by angle a
 Returns :  
 Args    :  1-3) three components defining the axis
            4)   rotation angle

=cut

sub matRotOrigAxis {
  my $u = shift;
  my $v = shift;
  my $w = shift;
  my $a = shift;

  if (sqrt($u**2 + $v**2) < 10**(-8)) {
    return matRotZ($a);
  }

  # matrix for rotating vector (u,v,w) about Z-axis into XY-plane
  my $uv = sqrt($u**2 + $v**2);
  my $Rxz = pdl([$u/$uv, $v/$uv, 0, 0], [-$v/$uv, $u/$uv, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]);

  # matrix for rotating vector (u,v,w) about Y-axis to the Z-axis
  my $uvw = sqrt($u**2 + $v**2 + $w**2);
  my $Rxz2z = pdl([$w/$uvw, 0, -$uv/$uvw, 0], [0, 1, 0, 0], [$uv/$uvw, 0, $w/$uvw, 0], [0, 0, 0, 1]);

  # matrix for rotating around vector (u,v,w) by angle a
  return inv($Rxz) x inv($Rxz2z) x matRotZ($a) x $Rxz2z x $Rxz;
}


=head2   FUNCTION -- matTransX

 Title   :  matTransX
 Usage   :  my $M = PDB::matTransX(3.2)
 Function:            
 Returns :  Translation matrix along X
 Args    :  1) displacement

=cut

sub matTransX {
  my $d = shift;
  return pdl([1, 0, 0, $d], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]);
}


=head2   FUNCTION -- matTransY

 Title   :  matTransY
 Usage   :  my $M = PDB::matTransY(3.2)
 Function:            
 Returns :  Translation matrix along Y
 Args    :  1) displacement

=cut

sub matTransY {
  my $d = shift;
  return pdl([1, 0, 0, 0], [0, 1, 0, $d], [0, 0, 1, 0], [0, 0, 0, 1]);
}


=head2   FUNCTION -- matTransZ

 Title   :  matTransZ
 Usage   :  my $M = PDB::matTransZ(3.2)
 Function:            
 Returns :  Translation matrix along Z
 Args    :  1) displacement

=cut

sub matTransZ {
  my $d = shift;
  return pdl([1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, $d], [0, 0, 0, 1]);
}

=head2   FUNCTION -- matTrans

 Title   :  matTrans
 Usage   :  my $M = PDB::matTrans(1.2, 3.0, 3.2)
 Function:            
 Returns :  Translation matrix by the given vector
 Args    :  1. displacement along X
            2. displacement along Y
            3. displacement along Z

=cut

sub matTrans {
  my $dx = shift;
  my $dy = shift;
  my $dz = shift;
  return matTransZ($dz) x matTransY($dy) x matTransX($dx);
}


=head2   FUNCTION -- applyTransformationToPDB

 Title   :  applyTransformationToPDB
 Usage   :  $pdb->applyTransformationToPDB($M)
 Function:  Given a transformation matrix, apply it to all atoms in the PDB object.
 Returns :
 Args    :  1) transformation matrix

=cut

sub applyTransformationToPDB {
  my $pdb = shift;
  my $M = shift;

  foreach my $ch (@{$pdb->{chain}}) {
    foreach my $res (@{$ch->{res}}) {
      foreach my $a (@{$res->{atom}}) {
        my $coor = $M x pdl([$a->{xcoor}], [$a->{ycoor}], [$a->{zcoor}], [1]);
        $coor = $coor / $coor->at(0, 3); # divide by homo-geneous coordinate
        $a->{xcoor} = $coor->at(0, 0);
        $a->{ycoor} = $coor->at(0, 1);
        $a->{zcoor} = $coor->at(0, 2);
      }
    }
  }
}


=head2   FUNCTION -- applyTransformation

 Title   :  applyTransformation
 Usage   :  $pdb->applyTransformation($M)
            applyTransformation(\@atoms, $M)
 Function:  Same as applyTransformationToPDB, but can be used either for PDB objects or lists of atoms
 Returns :
 Args    :  1) transformation matrix

=cut

sub applyTransformation {
  my $pdb = shift;
  my $M = shift;
  my $atoms = (ref($pdb) =~ /PDB/) ? $pdb->conAtoms(undef, 1) : $pdb;

  foreach my $a (@$atoms) {
    my $coor = $M x pdl([$a->{xcoor}], [$a->{ycoor}], [$a->{zcoor}], [1]);
    $coor = $coor / $coor->at(0, 3); # divide by homo-geneous coordinate
    $a->{xcoor} = $coor->at(0, 0);
    $a->{ycoor} = $coor->at(0, 1);
    $a->{zcoor} = $coor->at(0, 2);
  }
}

=head2   FUNCTION

 Title   :  transformPoint
 Usage   :  PDB::transformPoint($M, 12.3, 45.7, 0.3)
 Function:  Applies the specified transformation to the specified point
 Returns :
 Args    :  1) transformation matrix
            2) x coordinate
            3) y coordinate
            4) z coordinate

=cut

sub transformPoint {
  my $M = shift;
  GENERAL::assert(scalar(@_) % 3 == 0, "expected a number of values divisible by 3 (treated as 3D coordinates)!");
  my @ans;

  for (my $i = 0; $i < @_; $i += 3) {
    my $coor = $M x pdl([$_[$i]], [$_[$i+1]], [$_[$i+2]], [1]);
    $coor = $coor / $coor->at(0, 3); # divide by homo-geneous coordinate
    push(@ans, $coor->at(0, 0));
    push(@ans, $coor->at(0, 1));
    push(@ans, $coor->at(0, 2));
  }
  return @ans;
}

=head2   FUNCTION

 Title   :  appendStructure
 Usage   :  $mol1->appendStructure($mol2)
 Function:  Appends all the chains in $mol2 to the structure in $mol1.
            New chain names are given according to the standard in _newChain()
 Returns :
 Args    :

=cut

sub appendStructure {
  my $dpdb = shift;
  my $spdb = shift;

  foreach my $ch (@{$spdb->{chain}}) {
    my $nch = $dpdb->_newChain("");
    PDB::copyChain($ch, $nch);
  }
}


=head2   FUNCTION

 Title   :  getBackboneDihedrals
 Usage   :  my ($phi, $psi, $omega) = $mol->getBackboneDihedrals()
 Function:  Returns the definitions and the values of backbone dihedral angles
 Returns :
 Args    :

=cut

sub getBackboneDihedrals {
  my $pdb = shift;
  my (@phi, @psi, @omega);
  my $CTre = '^CY$';
  my $Cre = '^C$';
  my $Nre = '^N$';
  my $NTre = '^NT$';
  my $CAre = '^CA$';

  foreach my $chain (@{$pdb->{chain}}) {
    my $L = scalar(@{$chain->{res}});
    for (my $ri = 0; $ri < $L; $ri++) {
      my $pres = undef; my $cres = undef; my $nres = undef;
      $pres = $chain->{res}->[$ri-1] if ($ri > 0);
      $cres = $chain->{res}->[$ri];
      $nres = $chain->{res}->[$ri+1] if ($ri < $L - 1);

      # get phi
      if (($ri > 0) || (scalar(PDB::getAtomsMatchingName($cres, $CTre)) > 0)) {
        my @A = ($ri > 0) ? PDB::getAtomsMatchingName($pres, $Cre) : PDB::getAtomsMatchingName($cres, $CTre); GENERAL::assert(scalar(@A) == 1);
        my @B = PDB::getAtomsMatchingName($cres, $Nre); GENERAL::assert(scalar(@B) == 1);
        my @C = PDB::getAtomsMatchingName($cres, $CAre); GENERAL::assert(scalar(@C) == 1);
        my @D = PDB::getAtomsMatchingName($cres, $Cre); GENERAL::assert(scalar(@D) == 1);
        my @dihe = ($A[0], $B[0], $C[0], $D[0]);
        my %dihe; $dihe{def} = \@dihe; $dihe{val} = GENERAL::Dihedral(@dihe); $dihe{type} = "phi";
        push(@phi, \%dihe);
      }

      # get psi
      if (($ri < $L - 1) || (scalar(PDB::getAtomsMatchingName($cres, $NTre)) > 0)) {
        my @A = PDB::getAtomsMatchingName($cres, $Nre); GENERAL::assert(scalar(@A) == 1);
        my @B = PDB::getAtomsMatchingName($cres, $CAre); GENERAL::assert(scalar(@B) == 1);
        my @C = PDB::getAtomsMatchingName($cres, $Cre); GENERAL::assert(scalar(@C) == 1);
        my @D = ($ri < $L - 1) ? PDB::getAtomsMatchingName($nres, $Nre) : PDB::getAtomsMatchingName($cres, $NTre); GENERAL::assert(scalar(@D) == 1);
        my @dihe = ($A[0], $B[0], $C[0], $D[0]);
        my %dihe; $dihe{def} = \@dihe; $dihe{val} = GENERAL::Dihedral(@dihe); $dihe{type} = "psi";
        push(@psi, \%dihe);
      }

      # get omegga
      if ($ri < $L - 1) {
        my @A = PDB::getAtomsMatchingName($cres, $CAre); GENERAL::assert(scalar(@A) == 1);
        my @B = PDB::getAtomsMatchingName($cres, $Cre); GENERAL::assert(scalar(@B) == 1);
        my @C = PDB::getAtomsMatchingName($nres, $Nre); GENERAL::assert(scalar(@C) == 1);
        my @D = PDB::getAtomsMatchingName($nres, $CAre); GENERAL::assert(scalar(@A) == 1);
        my @dihe = ($A[0], $B[0], $C[0], $D[0]);
        my %dihe; $dihe{def} = \@dihe; $dihe{val} = GENERAL::Dihedral(@dihe); $dihe{type} = "omega";
        push(@omega, \%dihe);
      }
    }
  }

  return (\@phi, \@psi, \@omega)
}


=head2  FUNCTION -- getResChi

 Title   : getChi
 Usage   : $mol->getResChi($rotamer)
 Function: give a ROMATER object, return sets of atoms defining all chi angles
 Returns : sets of atoms defining all chi angles
 Args    : (1) rotamer object

=cut

sub getChi {
  my $self = shift;
  my $rotamer = shift;
  GENERAL::requireArgs($self, $rotamer);

  my @chis;
  my @res = $self->ConRes();
  foreach my $res (@res) {
    # loop over and get chis
    for (my $i = 0; $i < scalar(@{$rotamer->{chidef}->{$res->{resname}}}); $i++) {
      my $c = $rotamer->{chidef}->{$res->{resname}}->[$i];
      my @chi = (PDB::getAtomInRes($res, $c->[0]), PDB::getAtomInRes($res, $c->[1]), PDB::getAtomInRes($res, $c->[2]), PDB::getAtomInRes($res, $c->[3]));
      my %dihe; $dihe{def} = \@chi; $dihe{val} = GENERAL::Dihedral(@chi); $dihe{type} = "chi";
      push (@chis, \%dihe);
    }
  }
  return \@chis;
}


=head2   FUNCTION

 Title   :  getImropersForCB
 Usage   :  my $impr = $mol->getImropersforCB()
 Function:  Returns the definitions and the values of improper angle dihedrals for defining the CB atom, defined by: N-C-CA-CB
 Returns :
 Args    :

=cut

sub getImropersForCB {
  my $pdb = shift;
  my @impr;

  foreach my $chain (@{$pdb->{chain}}) {
    foreach my $res (@{$chain->{res}}) {
      next if ($res->{resname} =~ /GLY/);
      my @dihe = (PDB::getAtomInRes($res, "N"), PDB::getAtomInRes($res, "C"), PDB::getAtomInRes($res, "CA"), PDB::getAtomInRes($res, "CB"));
#      my @dihe = (PDB::getAtomInRes($res, "CA"), PDB::getAtomInRes($res, "N"), PDB::getAtomInRes($res, "C"), PDB::getAtomInRes($res, "CB"));
      my %dihe; $dihe{def} = \@dihe; $dihe{val} = GENERAL::Dihedral(@dihe); $dihe{type} = "impr";
      push(@impr, \%dihe);
    }
  }
  return \@impr;
}



=head2   FUNCTION

 Title   :  buildAtom
 Usage   :  PDB::buildAtom($atom, $distanceA, $angleA, $torsionA, $dist, $angle, $torsion)
 Function:  Builds the given atom using three placed atoms, an distance, angle and a dihedral
 Returns :
 Notes   : Taken from MSL's implementation
 Args    :

=cut

sub buildAtom {
  my $a = shift;
  my $da = shift;
  my $aa = shift;
  my $ta = shift;
  my $di = shift;
  my $an = shift;
  my $th = shift;
  my $pi = 3.141592653589793238462643383279502884;
  $an *= $pi/180;
  $th *= $pi/180;

  # unit vector from _distAtom to _angleAtom (B - C)
  my @uCB = ($da->{xcoor} - $aa->{xcoor}, $da->{ycoor} - $aa->{ycoor}, $da->{zcoor} - $aa->{zcoor});
  @uCB = GENERAL::getUnitVector(\@uCB);

  # vector from _angleAtom to _dihedralAtom (C - D)
  my @dDC = ($aa->{xcoor} - $ta->{xcoor}, $aa->{ycoor} - $ta->{ycoor}, $aa->{zcoor} - $ta->{zcoor});

  my $an2 = $pi - $an;
  my $th2 = $pi + $th;
  my $rsin = $di * sin($an2);
  my $rcos = $di * cos($an2);
  my $rsinsin = $rsin * sin($th2);
  my $rsincos = $rsin * cos($th2);

  my @c1 = GENERAL::cross(\@uCB, \@dDC);
  if (sqrt($c1[0]**2 + $c1[1]**2 + $c1[2]**2) < 10**(-8)) { return 0; }
  @c1 = GENERAL::getUnitVector(\@c1);
  GENERAL::scale(\@c1, $rsinsin);
  my @c2 = @uCB;
  GENERAL::scale(\@c2, -GENERAL::dot(\@dDC, \@uCB));
  @c2 = GENERAL::vecSum(\@dDC, \@c2);
  @c2 = GENERAL::getUnitVector(\@c2);
  GENERAL::scale(\@c2, $rsincos);
  my @c3 = @uCB;
  GENERAL::scale(\@c3, $rcos);
  my @dd = GENERAL::vecSum(\@c1, \@c2);
  @dd = GENERAL::vecSum(\@dd, \@c3);

  # update coordinate of placed atom
  $a->{xcoor} = $da->{xcoor} + $dd[0];
  $a->{ycoor} = $da->{ycoor} + $dd[1];
  $a->{zcoor} = $da->{zcoor} + $dd[2];

  return 1;
}

=head2  FUNCTION -- readPDBWater

 Title   :  readPDBWater
 Usage   :  $mol->readPDBWater( );
 Function:  read waters in PDB files
 Returns :
 Args    :  none

=cut

sub readPDBWater {

    my $self=shift;
    my $pdb=$self->{coordfile};
    my $fname=&GENERAL::GetInFH($pdb);

    my $watercount=0;
    my $iresnum;

    while (<$fname>) {

    next if (! /^HETATM/);

    if (/^HETATM/ && /HOH/) {
        $watercount++ if ( / O /);
        $iresnum=$watercount;

        my ($atominx, $atomname, $resname, $resnum, $x, $y, $z)
            = &GENERAL::Trim(substr($_, 6,5), substr($_, 12,4), substr($_,17,4),
                     substr($_, 22,5), substr($_, 30,8), substr($_, 38,8), substr($_,46,8));


        my $pdbrec={};

        $pdbrec->{atominx}=$atominx+0;
        $pdbrec->{atomname}=$atomname;
        $pdbrec->{resname}=$resname;
        $pdbrec->{resnum}=$resnum;
        $pdbrec->{iresnum}=$iresnum;
        $pdbrec->{seg}='SOLV';
        $pdbrec->{xcoor}=$x+0.0;
        $pdbrec->{ycoor}=$y+0.0;
        $pdbrec->{zcoor}=$z+0.0;
        $pdbrec->{hyd}=($atomname=~/^H.*/)?1:0;
        $pdbrec->{charge}=0.0;

        push (@{$self->{water}->{atom}}, $pdbrec);
        }
    }

    undef $fname;
}



=head2  FUNCTION -- writeCRDWater

 Title   :  writeCRDWater
 Usage   :  $mol->writeCRDWater($crdfilehandle, 'NOH');
 Function:  Append water molecule to the end of crd files.
 Returns :
 Args    :

=cut

sub writeCRDWater {

    my $self=shift;
    my $crd=shift;
    my $fname=&GENERAL::GetAppendOutFH($crd);
    my $format=shift;
    die "must specifiy a crd file format!\n" if (!defined $format);

    my $translate= uc($format);

    foreach my $a ( @{$self->{water}->{atom}} ) {

    my $ta;
    %{$ta}=%{$a};

    $ta->{resname}=~s/HOH/TIP3/;
    $ta->{chain}='SOLV';
    $ta->{atomname}='OH2';

    if ($translate =~ /NOH/ ) {
        printf $fname "%s\n",&_crdLine($ta)
        if ( $ta->{hyd} == 0);
    } else {
        printf $fname "%s\n",&_crdLine($ta);
    }
    }

    &GENERAL::FixCRDTotalNum($crd, 'temp.crd');

    undef $fname;
}


=head2   FUNCTION -- writePDBWater

 Title   :  writePDBWater
 Usage   :  $mol->writePDBWater( );
 Function:  Append water molecule to the end of crd files
 Returns :  none
 Args    :  none
 Note    :  to be implemeted.
            can some one do it?
=cut





=head2  FUNCTION -- sMutate

 Title   :  sMutate (symbolic mutate)
 Usage   :  sMutate($res, $aa_name, $clone_flag)
 Function:  Symbolic mutation. One can only mutate to a smaller residue because
            all this function does is through away unecessary atoms. This function
            was originally made for mutating to GLY or ALA.
 Returns :  a copy of the original residue (if asked for)
 Args    :  1. $res - pointer to the residue to mutate.
            2. $aa_name - the name of the amino acid to mutate to.
            3. $clone_flag - optional parameter. If true, a copy of the original
            residue is returned. Default is false.
            4. Optional: strict flag.  If set to 1 (default) will through an error
               when trying to mutate a GLY or a PRO to anything other than themselves.
               If set to 2, will through a warning but will procede with the mutantion
               ignoring problems. If 3 (or anything else), will procede with the 
               mutation silently.
=cut

sub sMutate {
  my $res = shift;
  my $aa = uc(shift);
  my $clone = shift;
  if (!defined($clone)) {$clone = 0;}
  my $strict = shift;
  if (!defined($strict)) { $strict = 1; }

  my $cres;
  if ($clone) {$cres = PDB::cloneResidue($res);}
  my $bbexp = PDB::backboneA("REGEXP");

  if ($aa =~ /ALA/) {
    # Currently, don't know how to change a GLY into an ALA
    if ((($res->{resname} eq "GLY") || ($res->{resname} eq "PRO")) && !($res->{resname} eq $aa)) {
      GENERAL::error("Don't know how to mutate a $res->{resname} to $aa") if ($strict == 1);
      GENERAL::warning("Don't know how to mutate a $res->{resname} to $aa") if ($strict == 2);
    }
    $res->{resname} = "ALA";
    my $k = 0;
    my $atom_arr;
    foreach my $atom (@{$res->{atom}}) {
      # Do not save unwanted atoms
      if ($atom->{atomname} =~ /($bbexp|CB)/) {
        $atom_arr->[$k] = $atom;
        $k++;
      }
    }
    undef(@{$res->{atom}});
    $res->{atom} = $atom_arr;
  } elsif ($aa =~ /GLY/) {
    $res->{resname} = "GLY";
    my $k = 0;
    my $atom_arr;
    foreach my $atom (@{$res->{atom}}) {
      # Do not save unwanted atoms
      if ($atom->{atomname} =~ /$bbexp/) {
        $atom_arr->[$k] = $atom;
        $k++;
      }
    }
    undef(@{$res->{atom}});
    $res->{atom} = $atom_arr;
  } elsif ($aa =~ /NONE/) {
    # This means remove the side chain (the same as mutating to GLY, but do not rename)
    my $k = 0;
    my $atom_arr;
    foreach my $atom (@{$res->{atom}}) {
      # Do not save unwanted atoms
      if ($atom->{atomname} =~ /$bbexp/) {
        $atom_arr->[$k] = $atom;
        $k++;
      }
    }
    undef(@{$res->{atom}});
    $res->{atom} = $atom_arr;
  } elsif ($aa =~ /GAP/) {
    @{$res->{atom}} = ();
  } else {
    GENERAL::error("Don't know how to mutate to $aa");
  }

  if ($clone) {return $cres;}
}


=head2  FUNCTION -- createUnfoldmer

 Title   :  createUnfoldmer
 Usage   :  PDB::createUnfoldmer($res, $unfoldmer_length)
 Function:  Creates a PDB structure corresponding to an unfoldmer of a given length
            with the given residue being the design residue.
 Returns :  none
 Args    :  (1) reference to residue
            (2) length of unfoldmer

=cut

sub createUnfoldmer {
  my $res = shift;
  my $unfoldmer = shift;
  
  if ($unfoldmer%2 == 0) { GENERAL::error("Unfoldmer length must be odd - is $unfoldmer!"); }
  
  my $iresnum = $res->{iresnum};
  my $chain = $res->{chain};
  my $pdb = $res->{chain}->{parent};
  
  my $off = int($unfoldmer/2);
  $pdb->resetValidResidues(0);
  my $cr = 0; my $len = 0;
  for (my $i = $iresnum - $off; $i <= $iresnum + $off; $i++) {
    unless (($i <= 0) || ($i > scalar(@{$chain->{res}}))) {
      $chain->{res}->[$i-1]->{valid} = 1;
      $len += 1;
    }
    if ($i == $iresnum) { $cr = $len; }
  }
  my $updb = $pdb->clone(1);
  for (my $i = 1; $i <= $len; $i++) {
    $updb->{chain}->[0]->{res}->[$i-1]->{iresnum} = $i;
    unless ($i == $cr) { PDB::sMutate($updb->{chain}->[0]->{res}->[$i-1], 'GLY'); }
  }
  return ($updb, $cr, $len);
}

=head2  FUNCTION -- writeUnfoldFragment

 Title   :  writeUnfoldFragment
 Usage   :  $mol->writeUnfoldFragment($interfacereslist)
 Function:  write the pentapetide fragment pdb files
 Returns :  none
 Args    :  (1) rotamer object from ROTAMER.pm
            (2) $resmap (for converting native naming in RCSB files
             and renumbered PDB files); optional

=cut


sub writeUnfoldFragment {
  my $self = shift;
  my $rotamer = shift;
  my $unfoldmer = shift;
  my $resmap = shift;

  if (!defined $unfoldmer) { die "Error in writeUnfoldFragment: unfoldmer size not specified!\n"; }

  # now loop over rotamer object, first residue list
  foreach my $s (@{$rotamer->{dslist}}) {
    my $chain = $s->{chain};
    my $iresnum = $s->{iresnum};
    my $pdb = $chain.$iresnum.'.pdb';
    my $crd = $chain.$iresnum.'.crd';
    my $dum = $chain.$iresnum.'.out';

    # first reset all the residue in the whole molecule as invalid
    my $c = $self->getChain($chain);
    $self->resetValidResidues(0);

    # now only set fragment: pentapeptide  as valid
    my $rescounter = 0;
    my $ii = ($unfoldmer+1)*0.5;
    my $jj = $unfoldmer-$ii-1;
    for (my $i = $iresnum-$ii; $i <= $iresnum+$jj; $i++) {
      # now check the integrity of $i in case it is close to terminal
      if ( $i < 0 || $i > $#{$c->{res}}) {
        # do nothing!!!!
      } else {
        my $r = $c->{res}->[$i];
        $r->{valid}=1 if (defined $r);
        $rescounter++;
      }
    }

    # only copy those valid pentapeptide
    my $newmol = $self->clone(1);

    if (!defined $resmap) {
      # replace all residues in the fragment with GLY's except the central one
#1      my $cres;
      for (my $i = 0; $i < scalar(@{$newmol->{chain}->[0]->{res}}); $i++) {
        my $res = $newmol->{chain}->[0]->{res}->[$i];
#1        # Mutate all but the central one
#1        if ($iresnum != $res->{iresnum}) { PDB::sMutate($res, 'GLY', 0); }
#1        else { $cres = $res; }
        PDB::sMutate($res, 'GLY', 0);
        $res->{iresnum} = $i+1;
        $res->{resnum} = $i+1;
      }
      my $fname = GENERAL::GetOutFH($pdb);
      # we add a comment line in case we want to check it later on
      $fname->print("REMARK  $s->{bbindex} $chain$iresnum\n");
      # write it out as generic format without any hydrogen and renumber it
      $newmol->writePDB($fname, 'genericrenum');
      close ($fname);
      $newmol->writeCRD($crd, '');
#1      PDB::sMutate($cres, 'GLY', 0); # make all GLY for the dummy files
      $newmol->writeDUMMY($dum);
    } elsif (defined $resmap) {
# THIS PART DEPRICATED!!!!!!!!!!!!!!!!!!!!!!!
GENERAL::error("This part of the code is not supposed to be visited!");
      my $fname = GENERAL::GetOutFH($pdb);
      # if we specify the residue map. we want to write out a file with
      # numbering from crystal structure
      $fname->print("REMARK  $s->{bbindex} $chain$iresnum\n");
      $newmol->Remap2Crystal($resmap);
      # replace all residues in the fragment with GLY's except the central one
      for (my $i=0; $i < scalar(@{$newmol->{chain}->[0]->{res}}); $i++) {
        next if (2*$i == $unfoldmer-1);
        PDB::sMutate($newmol->{chain}->[0]->{res}->[$i], 'GLY', 0);
      }
      $newmol->writePDB($fname, 'genericnoh');
      close ($fname);
      $newmol->writeCRD($crd, '');
    }
  }
}



=head2  backboneA

 Title   :  backboneA()
 Usage   :  PDB::backboneA()
 Function:  String definition of backbone atom name. Either an array
            of atom names (with asteriscs where necessary)
            or a regular expression is returned.
 Returns :  String definition of backbone atom name
 Args    :  1. string swtich
                - if contains "WATER", water molecules will also be
                  included in the backbone.
                - if contains "REGEXP", the definition of a backbone
                  atom will be returned as a regular expression.

=cut
sub backboneA {
  my $switch = shift;
  if (!defined($switch)) { $switch = ""; }
  my @bba = ('N', 'NT', 'CA', 'C', 'CY', 'OY', 'CAY', 'O', 'OCT*', 'H', 'HY*', 'HA1', 'HN', 'HT*', 'OXT', 'OT1', 'OT2');
  my @wbba = ('OH2', 'H1', 'H2');

  my $bckbn;
  if (uc($switch) =~ /REGEXP/) {
    $bckbn = join("|", @bba);
    if (uc($switch) =~ /WAT/) {
      $bckbn = "$bckbn|" . join("|", @wbba);
    }
    $bckbn = "\^($bckbn)\$";
    $bckbn =~ s/\*/\.\*/g;
    return $bckbn;
  } else {
    if (uc($switch) =~ /WAT/) {
      push(@bba, @wbba);
    }
    return @bba;
  }
}


sub writeUnfoldFragment_old {

    my $self=shift;
    my $rotamer=shift;
    my $unfoldmer=shift;
    my $resmap=shift;

    if (!defined $unfoldmer) { die "Error in writeUnfoldFragment: unfoldmer size not specified!\n"; }
    # make a special directory for unfolding reference state calculations
    &GENERAL::cmkdir("unfold");
    GENERAL::cchdir('unfold');

    # now loop over rotamer object, first residue list
    foreach my $s (@{$rotamer->{dslist}}) {
    my $chain=$s->{chain};
    my $iresnum=$s->{iresnum};
    my $j=0;
    foreach my $r (@{$s->{reslist}}) {

        # give the fragment a name: H45_0.pdb, H45_1.pdb etc
        # H: chain id   45: residue ID (renumbered)
        # 0, 1, 2, 3    list of posible residues at design site H45

        my $resname=$r->{resname};
        my $pdb=$chain.$iresnum.'.pdb';
        my $crd=$chain.$iresnum.'.crd';
        my $fname = &GENERAL::GetOutFH($pdb);

        # first reset all the residue in the whole molecule as invalid
        my $c=$self->getChain($chain);
        $self->resetValidResidues(0);


        my $rescounter=0;

        # now only set fragment: pentapeptide  as valid
        my $ii=($unfoldmer+1)/2;
            my $jj=$unfoldmer-$ii-1;

        for (my $i=$iresnum-$ii; $i <=$iresnum+$jj; $i++) {
        # now check the integrity of $i in case it is close to terminal

        if ( $i < 0 || $i > $#{$c->{res}}) {
            # do nothing!!!!
        } else {
            my $r=$c->{res}->[$i];
            $r->{valid}=1 if (defined $r);
            $rescounter++;
        }
        }

        # only copy those valid pentapeptide
        my $newmol = $self->clone(1);

        if (! defined $resmap) {

        # we add a comment line in case we want to check it later on
        $fname->print("REMARK   $s->{bbindex}  $iresnum  $resname \n");
        # write it out as generic format without any hydrogen and renumber it
        $newmol->writePDB($fname, 'genericrenum');
        close ($fname);

        # we have to do some post-process to clean up the pdb files.
                # (1) residue numbering
                # (2) residue naming
                # (3) chop off all the unnecessary things

        PDB::rewriteFragmentPDB($pdb, $iresnum);

        # now convert it into crd file
                # I need to add a comment line.-----to do list
        my $newmol2=PDB::new($pdb);
        $newmol2->writeCRD($crd, '');
                $newmol2->writePDB('tmp.pdb', 'GENERIC_RENUM');
        GENERAL::csystem("rm $pdb");
                GENERAL::csystem("mv tmp.pdb $pdb");

        } elsif ( defined $resmap) {

        # if we specify the residue map. we want to write out a file with
        # numbering from crystal structure

        $fname->print("REMARK   $r->{bbindex}  $iresnum  $resname \n");
        $newmol->Remap2Crystal($resmap);
        $newmol->writePDB($fname, 'genericnoh');
        close ($fname);
        PDB::rewriteFragmentPDB($pdb, $iresnum);

        # now convert it into crd file
        my $newmol2=PDB::new($pdb);
        $newmol2->writeCRD($crd, "");
#        &GENERAL::Remove($pdb);

        }
        $j++;
    }
    }

    # go home. don't get lost in the maze of directory tree.
    GENERAL::cchdir("../");
}


sub rewriteFragmentPDB {

    my $pdb=shift;
    my $iresnum=shift;
    my $res=shift;
    $res='GLY' if (!defined $res);
    my $fname=&GENERAL::GetInFH("$pdb");
    my $temp='temp'.$iresnum.'.pdb';
    my $temppdb=&GENERAL::GetOutFH($temp);
    #print $iresnum, "\n";

    while ( <$fname> ) {

        # print out the comment line
        $temppdb->print($_)   if ( /^REMARK/);

        if ( /^ATOM/ ) {

            # retrienve the information of residue name and indexed residue number
            my ($resname, $resnum)= &GENERAL::Trim(substr($_,17,4), substr($_, 22,4));
            #print $resnum,"\n";

            if ($res eq 'GLY' ) {
               # we want to make everything is GLY except the middle one.
               if ( $resnum != $iresnum && (/CA | C  | O | N | H /) ) {
                  substr($_, 17,4) =sprintf("%4s", 'GLY ');
                  my $resno= $resnum-$iresnum+4;
                  substr($_,22,4) =sprintf("%4d", $resno);
                  $temppdb->print($_);
                } elsif ( $resnum == $iresnum && ( /CA | C | O | N | H / ) ) {
                  my $resno= $resnum-$iresnum+4;
                  substr($_,22,4) =sprintf("%4d", $resno);
                  $temppdb->print($_);
                }
            } elsif ($res eq 'ALA') {
               if ( $resnum != $iresnum && (/CA | C  | O | N | H | CB /) ) {
                  substr($_, 17,4) =sprintf("%4s", 'ALA');
                  my $resno= $resnum-$iresnum+4;
                  substr($_,22,4) =sprintf("%4d", $resno);
                  $temppdb->print($_);
                } elsif ( $resnum == $iresnum && ( /CA | C | O | N | H| CB /) ) {
                  my $resno= $resnum-$iresnum+4;
                  substr($_,22,4) =sprintf("%4d", $resno);
                  $temppdb->print($_);
                }
            } else {
            }
        }
    }

    GENERAL::csystem ("mv $temp $pdb");
}


=head2   FUNCTION -- writeInpSeq

 Title   :  writeInpSeq
 Usage   :  $pdb->writeInpSeq('blah.seq', 'AB');
 Function:  write a outfile close to rot_input for all the residues
 Returns :  none
 Args    :  (1) file name to containing the sequence
            (2) chain ids (optional)

=cut



sub writeInpSeq {

    my $self=shift;
    my $fname=shift;
    my $chainstr=shift;

    my $outf=&GENERAL::GetOutFH($fname);


    if (! defined $chainstr) {
        $self->selectAllChains();
    } else {
        $self->selectChains($chainstr);
    }

    foreach my $c ( @{$self->{pickedchains}}) {
        my $chain=$c->{id};
    foreach my $r (@{$c->{res}}) {
       my $res=$r->{resname};
           my $index=$r->{iresnum};
           $outf->printf("%1s%4d%5s\n", $chain, $index, $res);
        }
    }

}


=head2   FUNCTION -- writeFasta

 Title   :  writeFasta
 Usage   :  $mol->writeFasta('blah.fasta')
 Function:  finally, i add this functions.
            read pdb file, write out a fasta files
 Returns :  none
 Args    :  fasta file to be write out

=cut

sub writeFasta {

    my $self=shift;
    my $fname=shift;
    my $name = shift;
    my $chainstr=shift;

    my $outf=&GENERAL::GetOutFH($fname);


    if (! defined $chainstr) {
        $self->selectAllChains();
    } else {
        $self->selectChains($chainstr);
    }

    foreach my $c ( @{$self->{pickedchains}}) {
        my $chain=$c->{id};
    my $total=$#{$c->{res}}+1;
    $outf->print("> $name chain $chain\n");
    my $counter=0;
    foreach my $r (@{$c->{res}}) {
        $counter++;
        my $res=$r->{resname};
        my $index=$r->{iresnum};
        my $resabbrev=$_seqabbrev{$res};
        $outf->printf("%1s", $resabbrev);
        # every line contains 60 residues
        $outf->print("\n") if ($counter%60 == 0);
        }
        $outf->print("\n");
    }

}


=head2   FUNCTION -- getAtomCoord

 Title   : getAtomCoord
 Usage   : $mol->getAtomCoord(H, 15, CB)
 Function: given a chain id and residue name and atom name, return coordinates
 Returns : reference to atom
 Args    : none

=cut



sub getAtomCoord {

    my $self=shift;
    my $chain=shift;
    my $iresnum=shift;
    my $atomname=shift;
   # my @atomcoord=[];

    my $c=$self->getChain($chain);
    my $counter=0;

    # loop oever all the residues
    foreach my $r (@{$c->{res}}) {

    # find the one with right index residue number
    if ($r->{iresnum} == $iresnum) {
            #print "the residue number is $r->{iresnum}.\n";
        # loop over residue atoms
        for (my $i=$r->{start}; $i<=$r->{end}; $i++) {

        if ($c->{atom}->[$i]->{atomname} eq $atomname) {

                    #print "the atom name is $c->{atom}->[$i]->{atomname}.\n";
            $counter++;
            # print return an atom object
            return $c->{atom}->[$i];
        }
        }
    }
    }


    if ($counter==0) {
    warn "$chain  $iresnum  $atomname couldn't be found.\n";
        return undef;
    }
}


=head2   FUNCTION -- getDihedral

 Title   :  getDihedral
 Usage   :  $mol->getDihedral(CA, CB, CG, CD)
 Function:  given 4 atoms, spit the dihedral value
 Returns :  chi, phi, psi value
 Args    :  4 atom names specificy the dihedral

=cut


sub getDihedral {
  my $res = shift;
  my $chiatom1 = shift;
  my $chiatom2 = shift;
  my $chiatom3 = shift;
  my $chiatom4 = shift;
  GENERAL::requireArgs($res, $chiatom1, $chiatom2, $chiatom3, $chiatom4);

  # Find the actual atoms
  my ($atom1, $atom2, $atom3, $atom4) = (undef, undef, undef, undef);
  foreach my $atom (@{$res->{atom}}) {
    if ($atom->{atomname} eq $chiatom1) {
      if (defined($atom1)) { GENERAL::error("Found two occurances of atom named $chiatom1 in residue " . PDB::resStr($res)); }
      $atom1 = $atom;
    } elsif ($atom->{atomname} eq $chiatom2) {
      if (defined($atom2)) { GENERAL::error("Found two occurances of atom named $chiatom2 in residue " . PDB::resStr($res)); }
      $atom2 = $atom;
    } elsif ($atom->{atomname} eq $chiatom3) {
      if (defined($atom3)) { GENERAL::error("Found two occurances of atom named $chiatom3 in residue " . PDB::resStr($res)); }
      $atom3 = $atom;
    } elsif ($atom->{atomname} eq $chiatom4) {
      if (defined($atom4)) { GENERAL::error("Found two occurances of atom named $chiatom4 in residue " . PDB::resStr($res)); }
      $atom4 = $atom;
    }
    if (defined($atom1) && defined($atom2) && defined($atom3) && defined($atom4)) {last;}
  }
  if (!(defined($atom1) && defined($atom2) && defined($atom3) && defined($atom4))) {
    GENERAL::error("Could not find one of $chiatom1, $chiatom2, $chiatom3 or $chiatom4 in residue " . PDB::resStr($res));
  }
  my $dihe = GENERAL::Dihedral($atom1, $atom2, $atom3, $atom4);
  return $dihe;
}


=head2  FUNCTION -- getResChi

 Title   : getResChi
 Usage   : $mol->getResChi(H, 15, ASN, $rotamer)
 Function: give a residue, spit out chi values
 Returns : chi values for the residue
 Args    : (1) reference to chain
           (2) indexed residue number
           (3) residue name
           (4) reference to rotamer object

=cut

sub getResChi {
  my $res = shift;
  my $reschidef = shift;
  GENERAL::requireArgs($res, $reschidef);

  my @chis;
  # loop over and get chis
  for (my $i=0; $i<=$#{$reschidef}; $i++) {
    my $c = $reschidef->[$i];
    my $di = PDB::getDihedral($res, $c->[0], $c->[1], $c->[2], $c->[3]);
    if (defined $di) {
      push (@chis, $di);
    } else {
      GENERAL::error("Dihedral angle \"$c->[0] $c->[1] $c->[2] $c->[3]\" for residue " . PDB::resStr($res) . " could not be found!");
    }
  }
  return \@chis;
}


=head2   FUNCTION -- getRotChi

 Title   : getRotchi
 Usage   : $mol->getRotChi($rotamer)
 Function: given a rotamer object, spit out all the chi values
 Returns : write out all the chi values in a file called CrystalChi
 Args    : reference to rotamer

=cut

sub getRotChi {
  my $self=shift;
  my $rotamer=shift;
  GENERAL::requireArgs($self, $rotamer);

  my $sitelist = $rotamer->{dslist};
  my @chis;

  foreach my $s (@{$sitelist}) {
    my $cid = $s->{chain};
    my $iresnum = $s->{iresnum};
    my $chain = $self->getChain($cid);
    my $res = $chain->{res}->[$iresnum-1];
    my $resname = $res->{resname};
    my $chidef = $rotamer->{chidef};
    my $reschidef = $chidef->{$resname};
    my $bbindex = $s->{bbindex};
    my $chi = PDB::getResChi($res, $reschidef);
    push (@chis, $chi);
  }
  
  return \@chis;
}


=head2   FUNCTION -- chiDiff

 Title   :  chiDiff
 Usage   :  ROTAMER::PDB($angle1, $angle2, $chi_angle_number, $resname, $neg_deviation, $pos_deviation)
 Function:  Determines whether the two specified chi angle values are the same within the specified deviations.
            Checks for flips where appropriate given the number of the chi angle (chi2 or chi3) and the residue name.
 Returns :  An array consisting of two values - the final difference in angles (between -180 and 180) taking flips into
            account and a flag designating the result of the comparison. 1 means the angles matched, 2 means the angles
            matched after a simmetrical flip, 3 means that angles matched after a flip of chi2 for ASN and chi3 for GLN
            (not really simmetrical, but in cristal structures often the difference between the nitrogen and the oxygen
            can not be seen), 0 means the angles did not match.
 Args    :  (1) first chi angle value
            (2) second chi angle value
            (3) chi angle number (1, 2, 3, ...)
            (4) residue name
            (5) positive deviation
            (6) negative deviation
=cut

sub chiDiff {
  my $chi1 = shift;
  my $chi2 = shift;
  my $chin = shift;
  my $resn = uc(shift);
  my $pdev = shift;
  my $ndev = shift;
  my $cor = 0;

  my $da = $chi1 - $chi2;
  $da += 360.0 if ($da < -180.0);
  $da -= 360.0 if ($da > 180.0);
  
  if ($da >= $ndev && $da <= $pdev) {
      $cor = 1;
  
  # try flipping for appropriate residues - JG used to have flipping of chi2 for HIS, but I don't think that's correct
  } elsif ((($chin == 2) && ($resn =~ /(ASP|PHE|TYR)/)) || (($chin == 3) && ($resn =~ /GLU/))) {
    my $da1 = $da + 180.0;
    $da1 += 360.0 if ($da1 < -180.0);
    $da1 -= 360.0 if ($da1 > 180.0);
    if ($da1 >= $ndev && $da1 <= $pdev) {
      $cor = 2;
      $da = $da1;
    }

  } elsif ((($chin == 2) && ($resn =~ /(ASN)/)) || (($chin == 3) && ($resn =~ /GLN/))){
    my $da1 = $da + 180.0;
    $da1 += 360.0 if ($da1 < -180.0);
    $da1 -= 360.0 if ($da1 > 180.0);
    $cor = 3 if ($da1 >= $ndev && $da1 <= $pdev);
  
  } else {
    $cor = 0;
  }
  return ($da, $cor);
}

=head2   FUNCTION -- getClosestLibrot

 Title   :  getClosestLibrot
 Usage   :  $mol->getClosestLibrot($rotamer)
 Function:  Return a structure containing the rotamers out of the rotamer
            library, which are closest to the given PDB structure rotamers.
 Returns :  none
 Args    :  (1) reference to rotamer object
            (2) flag: should the function check for flipped rotamers? Default is 0.
            (3) positive deviation (optional: default 40)
            (4) negative deviation (optional: defailt -40)

=cut

sub getClosestLibrot {

  my $self = shift;
  my $rotamer = shift;
  GENERAL::requireArgs($self, $rotamer);
  my $checkflip = shift;
  if (!defined($checkflip)) { $checkflip = 0; }
  my $pdeviation = shift;
  my $ndeviation = shift;
  $pdeviation = 40.0 if (!defined $pdeviation);
  $ndeviation = -40.0 if (!defined $ndeviation);

  my $chis = $self->getRotChi($rotamer);
  my $sitelist = $rotamer->{dslist};
  my $chidef = $rotamer->{chidef};
  my $crots = [];

  my $index = -1;
  foreach my $site (@{$sitelist}) {
    $index++;
    my $cid = $site->{chain};
    my $chain = $self->getChain($cid);
    my $iresnum = $site->{iresnum};
    my $bbindex = $site->{bbindex};
    my $res = $chain->{res}->[$iresnum-1];
    my $resname = $res->{resname};
    my $S;
    # will hold the match (in the general case, the array of matches)
    # because more than one rotamer out of the library may match
    $S = {}; $S->{librots} = []; $S->{chi} = [];
    for (my $j=0; $j < scalar(@{$chis->[$index]}); $j++) {
      push(@{$S->{chi}}, $chis->[$index]->[$j]);
    }
    push(@$crots, $S);

    if ($resname eq "ALA" || $resname eq "GLY") {
      my $s = {};
      $s->{rotamer} = 0;
      $s->{libchi} = [];
      # was a flip neccessary to match this rotamer?: 
      # 1 - no flip was neccesarry, 2 - fair flip, 3 - unfair flip (chi2 for ASN and chi3 for GLN)
      $s->{flip} = 1;
      push(@{$S->{librots}}, $s);
      $S->{match} = 1;
      next;
    }

    # retrieve right back-dependent rotamer library from rotamer object
    my $libref = $site->{rotlibref};
    # retrieve the resdiues in rotamer library from reference to rotamer library
    my $resref = $rotamer->getResInRotlib($libref, $resname);
    my $found = 0;
    my $flip = 1;

    # now loop over library rotamers and compare them to the conformation at the current site in the cristal structure
    for (my $i=0; $i < scalar(@{$resref->{rotamerlist}}); $i++) {

      # retrieve the chi values
      my $libchi = $resref->{rotamerlist}->[$i]->{chi};
      my $correct = 0;
      for (my $j = 0; $j < scalar(@$libchi); $j++) {
        my ($da, $f) = PDB::chiDiff($libchi->[$j], $chis->[$index]->[$j], $j+1, $resname, $pdeviation, $ndeviation);
        if (($f == 1) || ($checkflip && ($f > 1))) { $correct++; $flip = GENERAL::max($flip, $f); }
        else { last; }
      }
        
      # if every chi value agrees with library one, we say aha, that's it.
      if ($correct == scalar(@$libchi)) {
        my $s = {};
        $s->{rotamer} = $i;
        $s->{flip} = $flip;
        for (my $j=0; $j < scalar(@{$libchi}); $j++) {
          push(@{$s->{libchi}}, $libchi->[$j]);
        }
        push(@{$S->{librots}}, $s);
        $S->{match} = 1;
        $found = 1;
      }
    }
    # if no rotamer matched the crystal structure, mark it as "x"
    if ($found == 0) {
      my $s = {};
      $s->{rotamer} = "x";
      $s->{libchi} = [];
      $s->{flip} = 1;
      push(@{$S->{librots}}, $s);
      $S->{match} = 0;
    }
  }
  
  return $crots;
}


=head2   FUNCTION -- createPatchedResidue()

 Title   :  createPatchedResidue
 Usage   :  PDB::createPatchedResidue($Topology, $resname, $patchname, $newresiduename)
 Function:  Creates a new version of the given residue, patched with the given patch.
 Returns :  Nothing, the new residue is automatically inserted into the topology.
 Args    :  (1) Topology structure.
            (2) Residue name to patch.
            (3) Patch name.
            (4) Optional: name to call the newly created residue. Default is 'resname_patch_patchname'

=cut

sub createPatchedResidue {
  my $T = shift;
  my $rn = shift;
  my $pn = shift;
  GENERAL::requireArgs($T, $rn, $pn);
  my $nrn = shift;
  my $strict = shift; $strict = 1 if (!defined($strict));
  $nrn = "$rn\_patch\_$pn" if (!defined($nrn));

  GENERAL::assert(defined($T->{residues}->{$rn}), "Residue named '$rn' not found in topology!");
  GENERAL::assert(defined($T->{residues}->{$pn}), "No patch named '$pn' found in topology!") if (uc($pn) ne "NONE");
  GENERAL::assert($T->{residues}->{$pn}->{patchf} == 1, "Residue '$pn' is not marked as a patch in topology!") if (uc($pn) ne "NONE");
  if (defined($T->{residues}->{$nrn})) {
    GENERAL::error("Patched residue $nrn already exists!") if ($strict);
    return;
  }

  # first, copy the target residue
  my $pres = $T->{residues}->{$pn};
  my $ores = $T->{residues}->{$rn};
  my $res = {};
  $res->{name} = $nrn;
  $res->{patchf} = $ores->{patchf};
  $res->{atoms} = {}; $res->{groups} = {};
  foreach my $a (@{$ores->{atoma}}) {
    my $na = {}; $na->{name} = $a->{name};
    $na->{group} = $a->{group};
    $na->{type} = $a->{type};
    $na->{charge} = $a->{charge};
    $res->{atoms}->{$na->{name}} = $na;
    push(@{$res->{groups}->{$na->{group}}}, $na);
    push(@{$res->{atoma}}, $na);
  }
  foreach my $ai (keys(%{$ores->{bonds}})) {
    foreach my $aj (keys(%{$ores->{bonds}->{$ai}})) {
      $res->{bonds}->{$ai}->{$aj} = 1;
      $res->{bonds}->{$aj}->{$ai} = 1;
    }
  }
  $T->{residues}->{$nrn} = $res;
  return if (uc($pn) eq "NONE");

  # delete necessary atoms
  foreach my $an (@{$pres->{delete}}) {
    GENERAL::assert(defined($res->{atoms}->{$an}), "Did not find atom named '$an' in topology of $rn");
    my $a = $res->{atoms}->{$an};
    delete $res->{atoms}->{$an}; # delete atom from atom hash
    my $i = GENERAL::findin($a, @{$res->{atoma}});
    splice(@{$res->{atoma}}, $i, 1); # delete atom from atom array
    my $j = GENERAL::findin($a, @{$res->{groups}->{$a->{group}}});
    splice(@{$res->{groups}->{$a->{group}}}, $j, 1);
    delete $res->{bonds}->{$an} if (defined($res->{bonds}->{$an})); # delete any bonds with this atom
    foreach my $ani (keys(%{$res->{bonds}})) {
      delete $res->{bonds}->{$ani}->{$an} if (defined($res->{bonds}->{$ani}->{$an}));
    }
  }

  # insert new atoms, one group at a time
  foreach my $gi (keys(%{$pres->{groups}})) {
    # find an atom in this group of the patch that already exists in the residue
    my $ca = undef; # common atom
    foreach my $ga (@{$pres->{groups}->{$gi}}) {
      if (defined($res->{atoms}->{$ga->{name}})) { $ca = $ga; last; }
    }
    # put atoms of group $gi of patch into group $gj of residue
    my $gj;
    if (defined($ca)) { $gj = $res->{atoms}->{$ca->{name}}->{group}; }
    else {
      $gj = GENERAL::max(keys(%{$res->{groups}})) + 1;
    }
    foreach my $na (@{$pres->{groups}->{$gi}}) {
      # for atom that are already in place (the common atom and possibly more), only need to adjust the charge
      if (defined($res->{atoms}->{$na->{name}})) {
        $res->{atoms}->{$na->{name}}->{charge} = $na->{charge}; next;
      }
      push(@{$res->{atoma}}, $na); # add new atom to atom array
      $res->{atoms}->{$na->{name}} = $na; # add new atom to atom hash
      push(@{$res->{groups}->{$gj}}, $na); # add new atom into the group
    }
  }

  # insert any new bonds
  foreach my $ai (keys(%{$pres->{bonds}})) {
    foreach my $aj (keys(%{$pres->{bonds}->{$ai}})) {
      $res->{bonds}->{$ai}->{$aj} = 1;
      $res->{bonds}->{$aj}->{$ai} = 1;
    }
  }
}

=head2   FUNCTION -- exclusionList()

 Title   :  exclusionList
 Usage   :  PDB::exclusionList(\@atoms, "1-2 1-3 1-4", $T)
 Function:  Returns the exclusion list as an array of atom pairs (i, j). Atom i belongs
            to the list passed to the function and atom j is part of the exclusion list
            for atom i. Atom j may also belong to the list passed to the function. To
            ignore such pairs, the last optional parameter can be used. Which particular
            list to generate depends on the second input. Assumes that the atoms are part
            of one molecule and that the {iresnum} field of every residue correctly
            represents the order of the residue in the respective chain.
 Returns :  A list of exclusions (atom pairs). Either breaks down different exclusion types
            or returns all in one list (depends on wantarray).
 Args    :  (1) Array of atoms to generate an exclusion list for.
            (2) String specifying the exclusion list (can be a concatenation).
            (3) Topology hash table.
            (4) Optional: this argument can be used to specify which pairs are to be selectively
                ignored in forming the exclusion list. If numeric, the meaning is the following:
                0 - include self pairs (i.e. pairs i-j, where both i and j are within the
                    specified array of atoms) - default.
                1 - ignore all self pairs.
                2 - consider ONLY self pairs.
                3 - separate out the self and non-self pairs into different lists.
                The argument can also be a string (non-numeric). If so, this string will be
                interpreted as the name of the field of the atom hash, which if exists for any
                atom, signifies that pairs where this atom is j are not to be considered.
            (5) Optional: if this flag is set (not set by default), each entry in exclusion lists
                will contain all atoms on the path between the two terminal atoms. Each unique
                atom sequence will then be considered as a unique list entry.

=cut

sub exclusionList {
  my $atoms = shift;
  my $ltype = shift;
  my $T = shift;
  GENERAL::requireArgs($atoms, $ltype, $T);
  my $sxf = shift; # self exclusion flag
  my $long = shift;
  $long = 0 if (!defined($long));
  my $tag; my @tatoms; # atoms that were tagged in this function (so that they can be untagged at the end)
  if (defined($sxf)) {
    if (GENERAL::isNumeric($sxf)) {
      if (($sxf == 1) || ($sxf == 2) || ($sxf == 3)) {
        $tag = "_exl_flg";
        foreach my $a (@$atoms) { $a->{$tag} = 1; push(@tatoms, $a); }
      } else {
        GENERAL::error("Argument sxf must be either 1, 2 or 3!");
      }
    } else { $tag = $sxf; $sxf = 1; }
  } else { $sxf = 0; }

  if (($sxf == 3) && !wantarray) { GENERAL::error("Function cannot be used in scalar environment with sxf = $sxf"); }

  my (@list12, @list12s);
  # Get the 1-2 list
  if ($ltype =~ /(1-2|1-3|1-4)/) {
    foreach my $a (@$atoms) {
      # guard against bonds being defined more than once (i.e. "C +N" and "-C N") - should not happen, but who knows
      my %bh;
      # find all atom names bonded to this one
      my $rn = $a->{residue}->{resname};
      if (defined($a->{residue}->{patch})) { # see if residue is patched
        $rn = $a->{residue}->{resname} . "_patch_" . $a->{residue}->{patch};
        # if not yet done, create topology for patched residue
        PDB::createPatchedResidue($T, $a->{residue}->{resname}, $a->{residue}->{patch}, $rn) if (!defined($T->{residues}->{$rn}));
      }
      my $an = $a->{atomname};
      if (!defined($T->{residues}->{$rn}->{bonds}->{$an})) {
        GENERAL::warning("No connectivity information found for atom $an in residue $rn");
        next;
      }
      # 1. check for all bonds {an} - {x} defined within this residue and decide where x is
      if (!defined($T->{residues}->{$rn}->{bonds}->{$an})) { GENERAL::error("Could not find connctivity info for atom $an in residue $rn"); }
      my @bns = keys(%{$T->{residues}->{$rn}->{bonds}->{$an}});
      if (scalar(@bns) == 0) { GENERAL::error("No atoms bond with atom $an in residue $rn"); }
      foreach my $bn (@bns) {
        my $res; # residue that the bonded atom belongs to
        if ($bn =~ /\+/) {
          $bn =~ s/\+//;
          if ($a->{residue}->{iresnum} == scalar(@{$a->{residue}->{chain}->{res}})) { next; }
          $res = $a->{residue}->{chain}->{res}->[$a->{residue}->{iresnum}];
        } elsif ($bn =~ /\-/) {
          $bn =~ s/\-//;
          if ($a->{residue}->{iresnum} == 1) { next; }
          $res = $a->{residue}->{chain}->{res}->[$a->{residue}->{iresnum}-2];
        } else {
          $res = $a->{residue};
        }
        my $b = PDB::getAtomInRes($res, $bn, 0);
        if ($b == 0) { GENERAL::error("Could not find second atom of the bond $an with $bn in residue $rn"); }
        next if (defined($bh{PDB::atomStr($a)}{PDB::atomStr($b)}));
        my @arr = ($a, $b); push(@list12, \@arr);
        $bh{PDB::atomStr($a)}{PDB::atomStr($b)} = 1; $bh{PDB::atomStr($b)}{PDB::atomStr($a)} = 1;
      }
      # 2. check for all bonds {-an} - {x} defined within the next residue and decide where x is
      if ($a->{residue}->{iresnum} < scalar(@{$a->{residue}->{chain}->{res}})) {
        my $nres = $a->{residue}->{chain}->{res}->[$a->{residue}->{iresnum}];
        my $nrn = $nres->{resname};
        if (defined($T->{residues}->{$nrn}->{bonds}->{"-$an"})) {
          my @bns = keys(%{$T->{residues}->{$nrn}->{bonds}->{"-$an"}});
          foreach my $bn (@bns) {
            if ($bn =~ /[\+\-]/) {
              GENERAL::error("Invalid bond between -$an and $bn in residue $nrn");
            }
            my $b = PDB::getAtomInRes($nres, $bn, 0);
            if ($b == 0) { GENERAL::error("Could not find second atom of the bond $an with $bn in residue $nrn"); }
            next if (defined($bh{PDB::atomStr($a)}{PDB::atomStr($b)}));
            my @arr = ($a, $b); push(@list12, \@arr);
            $bh{PDB::atomStr($a)}{PDB::atomStr($b)} = 1; $bh{PDB::atomStr($b)}{PDB::atomStr($a)} = 1;
          }
        }
      }

      # 3. check for all bonds {+an} - {x} defined within the previous residue and decide where x is
      if ($a->{residue}->{iresnum} > 1) {
        my $nres = $a->{residue}->{chain}->{res}->[$a->{residue}->{iresnum}-2];
        my $nrn = $nres->{resname};
        if (defined($T->{residues}->{$nrn}->{bonds}->{"+$an"})) {
          my @bns = keys(%{$T->{residues}->{$nrn}->{bonds}->{"+$an"}});
          foreach my $bn (@bns) {
            if ($bn =~ /[\+\-]/) {
              GENERAL::error("Invalid bond between -$an and $bn in residue $nrn");
            }
            my $b = PDB::getAtomInRes($nres, $bn, 0);
            if ($b == 0) { GENERAL::error("Could not find second atom of the bond $an with $bn in residue $nrn"); }
            next if (defined($bh{PDB::atomStr($a)}{PDB::atomStr($b)}));
            my @arr = ($a, $b); push(@list12, \@arr);
            $bh{PDB::atomStr($a)}{PDB::atomStr($b)} = 1; $bh{PDB::atomStr($b)}{PDB::atomStr($a)} = 1;
          }
        }
      }
    }
  }

  # To get the 1-3, take all 1-2 pairs (i-j) and find all 1-2 pairs j-k, where k is not in {i}
  my (@list13, @list13s);
  if ($ltype =~ /(1-3|1-4)/) {
    foreach my $exij (@list12) {
      $exij->[0]->{"_ignore_13"} = 1;
      my @na = ($exij->[1]);
      my $add = PDB::exclusionList(\@na, "1-2", $T, "_ignore_13");
      foreach my $exjk (@$add) {
        my @arr = ($exij->[0], $exij->[1], $exjk->[1]);
        push(@list13, \@arr);
      }
      delete($exij->[0]->{"_ignore_13"});
    }
  }

  # To get the 1-4, take all 1-3 pairs (i-k) and find all 1-2 pairs k-l, where l is not in {i} or {j}
  my (@list14, @list14s);
  if ($ltype =~ /1-4/) {
    foreach my $exik (@list13) {
      $exik->[0]->{"_ignore_14"} = 1;
      $exik->[1]->{"_ignore_14"} = 1;
      my @na = ($exik->[2]);
      my $add = PDB::exclusionList(\@na, "1-2", $T, "_ignore_14");
      foreach my $exkl (@$add) {
        my @arr = ($exik->[0], $exik->[1], $exik->[2], $exkl->[1]);
        push(@list14, \@arr);
      }
      delete($exik->[0]->{"_ignore_14"});
      delete($exik->[1]->{"_ignore_14"});
    }
  }

  # Ignore undesired exclusions, remove duplicates
  my %exl;
  if ($ltype =~ /1-2/) {
    my (@tmp, @tmps);
    foreach my $ex (@list12) {
      next if ( (($sxf == 1) && (defined($ex->[1]->{$tag}))) || (($sxf == 2) && (!defined($ex->[1]->{$tag}))) );
      my $stri = PDB::atomStr($ex->[0]) . PDB::atomStr($ex->[1]);
      my $strj = PDB::atomStr($ex->[1]) . PDB::atomStr($ex->[0]);
      next if (defined($exl{$stri}) || defined($exl{$strj}));
      my @arr = ($ex->[0], $ex->[1]);
      if (($sxf == 3) && (defined($ex->[1]->{$tag}))) { push(@tmps, \@arr); }
      else { push(@tmp, \@arr); }
      $exl{$stri} = 1;
    }
    @list12 = @tmp;
    @list12s = @tmps;
  }
  if ($ltype =~ /1-3/) {
    my (@tmp, @tmps);
    foreach my $ex (@list13) {
      next if ( (($sxf == 1) && (defined($ex->[2]->{$tag}))) || (($sxf == 2) && (!defined($ex->[2]->{$tag}))) );
      my $stri; my $strj;
      if ($long) {
        $stri = PDB::atomStr($ex->[0]) . PDB::atomStr($ex->[1]) . PDB::atomStr($ex->[2]);
        $strj = PDB::atomStr($ex->[2]) . PDB::atomStr($ex->[1]) . PDB::atomStr($ex->[0]);
      } else {
        $stri = PDB::atomStr($ex->[0]) . PDB::atomStr($ex->[2]);
        $strj = PDB::atomStr($ex->[2]) . PDB::atomStr($ex->[0]);
      }
      next if (defined($exl{$stri}) || defined($exl{$strj}));
      my @arr;
      if ($long) { @arr = ($ex->[0], $ex->[1], $ex->[2]); }
      else { @arr = ($ex->[0], $ex->[2]); }
      if (($sxf == 3) && (defined($ex->[2]->{$tag}))) { push(@tmps, \@arr); }
      else { push(@tmp, \@arr); }
      $exl{$stri} = 1;
    }
    @list13 = @tmp;
    @list13s = @tmps;
  }
  if ($ltype =~ /1-4/) {
    my (@tmp, @tmps);
    foreach my $ex (@list14) {
      next if ( (($sxf == 1) && (defined($ex->[3]->{$tag}))) || (($sxf == 2) && (!defined($ex->[3]->{$tag}))) );
      my $stri; my $strj;
      if ($long) {
        $stri = PDB::atomStr($ex->[0]) . PDB::atomStr($ex->[1]) . PDB::atomStr($ex->[2]) . PDB::atomStr($ex->[3]);
        $strj = PDB::atomStr($ex->[3]) . PDB::atomStr($ex->[2]) . PDB::atomStr($ex->[1]) . PDB::atomStr($ex->[0]);
      } else {
        $stri = PDB::atomStr($ex->[0]) . PDB::atomStr($ex->[3]);
        $strj = PDB::atomStr($ex->[3]) . PDB::atomStr($ex->[0]);
      }
      next if (defined($exl{$stri}) || defined($exl{$strj}));
      my @arr;
      if ($long) { @arr = ($ex->[0], $ex->[1], $ex->[2], $ex->[3]); }
      else { @arr = ($ex->[0], $ex->[3]); }
      if (($sxf == 3) && (defined($ex->[3]->{$tag}))) { push(@tmps, \@arr); }
      else { push(@tmp, \@arr); }
      $exl{$stri} = 1;
    }
    @list14 = @tmp;
    @list14s = @tmps;
  }

  # Remove tags
  foreach my $a (@tatoms) {
    delete($a->{$tag});
  }

  # Concatenate all arrays
  my @ans;
  if (wantarray) {
    push(@ans, \@list12) if ($ltype =~ /1-2/);
    push(@ans, \@list13) if ($ltype =~ /1-3/);
    push(@ans, \@list14) if ($ltype =~ /1-4/);
    if ($sxf == 3) {
      push(@ans, \@list12s) if ($ltype =~ /1-2/);
      push(@ans, \@list13s) if ($ltype =~ /1-3/);
      push(@ans, \@list14s) if ($ltype =~ /1-4/);
    }
    return @ans;
  } else {
    push(@ans, @list12) if ($ltype =~ /1-2/);
    push(@ans, @list13) if ($ltype =~ /1-3/);
    push(@ans, @list14) if ($ltype =~ /1-4/);
    return \@ans;
  }
}

sub printExclusionList {
  my $list = shift;
  
  foreach my $ex (@$list) {
    print(PDB::atomStr($ex->[0]));
    for (my $i = 1; $i < scalar(@$ex); $i++) { print(" : " . PDB::atomStr($ex->[$i])); }
    if (defined($ex->[0]->{ai})) { 
      print(" ($ex->[0]->{ai}");
      for (my $i = 1; $i < scalar(@$ex); $i++) { print(" : " . $ex->[$i]->{ai}); }
      print(")");
    }
    print "\n";
  }
}

=head2   FUNCTION -- dummyPDBfromSequence()

 Title   :  dummyPDBfromSequence
 Usage   :  PDB::dummyPDBfromSequence($pdbf, @chains)
 Function:  Writes a dummy CA trace given a sequence of amino acids.
 Returns :  nothing.
 Args    :  (1) Name of output PDB file.
            (2) Array of arrays, each representing a sequence of amino acids.

=cut

sub dummyPDBfromSequence {
  my $pdbf = shift;
  
  my $pdb = PDB::new();
  my $uresnum = 1;
  my $atomunx = 1;
  foreach my $cseq (@_) {
    my $chain = $pdb->_newChain("");
    my $iresnum = 1;
    my $atominx = 1;
    foreach my $aa (@$cseq) {
      my $res = PDB::_newResidue($aa, $iresnum, $iresnum, $uresnum, $chain, 1);
      my $a = PDB::_newAtom($atominx, $atomunx, "CA", 0, 0, 0, 0, 1, 0, $res);
      $atominx++;
      $atomunx++;
      $uresnum++;
      $iresnum++;
    }
  }
  $pdb->writePDB($pdbf, "");
}


=head2

 Title   :  getSolenoidParameters
 Usage   :  $pdb->getSolenoidParameters
 Function:  Assumes that the chains in the structure are forming an ideal solenoid
            (i.e., a helix) and computes its parameters.
 Returns :  Solenoid parameters and axis direction and location
 Args    :  (1) PDB object or list of atoms

=cut

sub getSolenoidParameters {
  my $pdb = shift;
  GENERAL::assert(scalar(@{$pdb->{chain}}) >= 4, "need at least four sub-units (chains) to unambiguously determine the solenoid axis!");
  my @cA = $pdb->centerOfMass($pdb->{chain}->[0]->{id});
  my @cB = $pdb->centerOfMass($pdb->{chain}->[1]->{id});
  my @cC = $pdb->centerOfMass($pdb->{chain}->[2]->{id});
  my @cD = $pdb->centerOfMass($pdb->{chain}->[3]->{id});

  # find centroid-to-centroid distances
  my @vecAB = GENERAL::vecDiff(\@cB, \@cA);
  my @vecBC = GENERAL::vecDiff(\@cC, \@cB);
  my @vecCD = GENERAL::vecDiff(\@cD, \@cC);

  # differences between adjacent centroid-to-centroid distances should
  # point directly into the helical axis
  my @u = GENERAL::vecDiff(\@vecBC, \@vecAB); @u = GENERAL::getUnitVector(\@u);
  my @v = GENERAL::vecDiff(\@vecCD, \@vecBC); @v = GENERAL::getUnitVector(\@v);

  # find helical radius by solving for how far along the two obtained
  # vectors, starting from their corresponding centroids, one must travel
  # to hit an axis that is perpendicular to both (the helical axis):
  # * the central axis is between points P + r*u and Q + r*v
  #   (where r is the radius and P and Q are the centroids corresponding to u and v)
  # * therefore, the central axis has the vector (P - Q) + r*(u - v)
  # * thus, we must simply solve dot((P - Q) + r*(u - v), v) = 0
  # * or, alternatively, dot((P - Q) + r*(u - v), u) = 0 (equivalent)
  my @PQ = GENERAL::vecDiff(\@cB, \@cC);
  my $R = GENERAL::dot(\@PQ, \@u)/(GENERAL::dot(\@u, \@v) - 1); # this is the radius of the helix

  # the direction of the axis itself can be simply found as the cross
  # product of the two vectors pointing into the helical axis
  my @H = GENERAL::cross(\@u, \@v);

  # also find d (rise per unit) and w (helical frequency)
  my $sgn = (GENERAL::dot(\@H, \@vecBC) > 0) ? 1 : -1; # is it right-handed (positive) or left-handed (negative)
  my $w = $sgn * abs(acos(GENERAL::dot(\@u, \@v)));
  my $d = GENERAL::vecNorm(GENERAL::arrayRef(GENERAL::vecSum(\@PQ, GENERAL::arrayRef(GENERAL::vecScale(GENERAL::arrayRef(GENERAL::vecDiff(\@u, \@v)), $R)))));

  # first, center solenoid on the the axis
  my @o = GENERAL::vecSum(\@cB, GENERAL::arrayRef(GENERAL::vecScale(\@u, $R)));

  return ($R, $w, $d, \@H, \@o);
}

=head2

 Title   :  alignWithSolenoidAxis
 Usage   :  $pdb->alignWithSolenoidAxis()
 Function:  Reorients the structure of a solenoid to place the solenoid axis along Z
 Returns :  Nothing
 Args    :  (1) PDB object. Assumes individual chains are arranged in a solenoid (helix).

=cut

sub alignWithSolenoidAxis {
  my $pdb = shift;
  my ($R, $w, $d, $H, $o) = $pdb->getSolenoidParameters();

  # first, center solenoid on the the axis
  $pdb->center(undef, $o);

  # and then, align the solenoid along the helical axis
  $pdb->applyTransformation(PDB::matAlignVectorWithZAxis(@$H));

  return ($R, $w, $d);
}

=head2

 Title   :  alignAlongPrincipalComponents
 Usage   :  $pdb->alignAlongPrincipalComponents() or $pdb->alignAlongPrincipalComponents(1, 2)
 Function:  Reorients the structure to place principal components along laboratory axes.
 Returns :  Nothing
 Args    :  (1) PDB object or list of atoms
            {2} Backbone flag, the same meaning as in PDB::getPrincipalComponents
            {3} Which axis to align the primary principal component along:
                1 - X-axis
                2 - Y-axis
                3 - Z-axis

=cut

sub alignAlongPrincipalComponents {
  my $pdb = shift;
  my $bbf = shift; $bbf = 1 if (!defined($bbf));
  my $ord = shift; $ord = 3 if (!defined($ord));
  GENERAL::assert(($ord == 1) || ($ord == 2) || ($ord == 3), "order parameter must be 1, 2, or 3");

  # compute principal components
  my @P = $pdb->getPrincipalComponents($bbf);

  # convert to new coordinate system
  my $M;
  if ($ord == 1) {
    $M = PDB::matAlignVectorWithXAxis($P[0]->[0], $P[0]->[1], $P[0]->[2]);
    my @pri = (1, 0, 0);
    my @sec = (0, 1, 0);
    my @ax2 = PDB::transformPoint($M, $P[1]->[0], $P[1]->[1], $P[1]->[2]);
    my $sgn = (GENERAL::dot(GENERAL::cross(\@ax2, GENERAL::arrayRef(@sec), 1), GENERAL::arrayRef(@pri)) > 0) ? 1 : -1;
    $M = PDB::matRotX($sgn * acos($ax2[1])) x $M;
  } elsif ($ord == 2) {
    $M = PDB::matAlignVectorWithYAxis($P[0]->[0], $P[0]->[1], $P[0]->[2]);
    my @pri = (0, 1, 0);
    my @sec = (0, 0, 1);
    my @ax2 = PDB::transformPoint($M, $P[1]->[0], $P[1]->[1], $P[1]->[2]);
    my $sgn = (GENERAL::dot(GENERAL::cross(\@ax2, GENERAL::arrayRef(@sec), 1), GENERAL::arrayRef(@pri)) > 0) ? 1 : -1;
    $M = PDB::matRotY($sgn * acos($ax2[2])) x $M;
  } elsif ($ord == 3) {
    $M = PDB::matAlignVectorWithZAxis($P[0]->[0], $P[0]->[1], $P[0]->[2]);
    my @pri = (0, 0, 1);
    my @sec = (1, 0, 0);
    my @ax2 = PDB::transformPoint($M, $P[1]->[0], $P[1]->[1], $P[1]->[2]);
    my $sgn = (GENERAL::dot(GENERAL::cross(\@ax2, GENERAL::arrayRef(@sec), 1), GENERAL::arrayRef(@pri)) > 0) ? 1 : -1;
    $M = PDB::matRotZ($sgn * acos($ax2[0])) x $M;
  }
  $pdb->applyTransformation($M);
}


=head2   FUNCTION PCA

 Title   :  getPrincipalComponents
 Usage   :  $pdb->getPrincipalComponents() or PDB::getPrincipalComponents(\@atoms)
 Function:  Performs Principal Component Analysis on the coordinates of the given PDB or list of atoms.
            Returns the three orthogonal components in the order of importance. The first component can
            often be used as the "axis" of the molecule for the purposes of aligning it.
 Returns :  Array of 3 array references, each with 3 elements.
 Args    :  (1) PDB object or list of atoms
            {2} Backbone flag:
                1 - consider backbone atoms only
                2 - consider CA atoms only
                3 - consider all non-hydrogen atoms

=cut

sub getPrincipalComponents {
  my $pdb = shift;
  my $bbf = shift;
  $bbf = 0 if (!defined($bbf));
  my @atoms;
  if (ref($pdb) =~ /ARRAY/) { @atoms = @$pdb; }
  elsif (ref($pdb) =~ /PDB/) { @atoms = $pdb->conAtoms(); }
  else { GENERAL::error("Unknown input parameter type: " . ref($pdb)); }
  my $bbn = ".";
  if ($bbf == 1) { $bbn = PDB::backboneA('REGEXP'); }
  elsif ($bbf == 2) { $bbn = "CA"; }
  elsif ($bbf == 3) { $bbn = "^[^H]"; }

  my @xyz;
  foreach my $a (@atoms) {
    next if ($bbf && ($a->{atomname} !~ /$bbn/));
    push(@{$xyz[0]}, $a->{xcoor}); push(@{$xyz[1]}, $a->{ycoor}); push(@{$xyz[2]}, $a->{zcoor});
  }

  # create covariance matrix
  my $M = zeros(3, 3);
  for (my $i = 0; $i < 3; $i++) {
    for (my $j = 0; $j < 3; $j++) {
      $M->slice("$i,$j") .= GENERAL::cov($xyz[$i], $xyz[$j]);
    }
  }

  # calculate the principal component of the covariance matrix
  my ($ve, $va) = eigens_sym $M;
  my @ord = sort { $va->at($b) <=> $va->at($a) } 0..2;
  my @p;
  for (my $i = 0; $i < scalar(@ord); $i++) {
    my $k = $ord[$i];
    my $mag = sqrt(($ve->at($k, 0))**2 + ($ve->at($k, 1))**2 + ($ve->at($k, 2))**2);
    my @tmp = ($ve->at($k, 0)/$mag, $ve->at($k, 1)/$mag, $ve->at($k, 2)/$mag);
    push(@p, \@tmp);
  }

  # make sure system is right handed
  if (GENERAL::dot(GENERAL::cross($p[0], $p[1], 1), $p[2]) < 0) {
      $p[0]->[0] = -$p[0]->[0];
      $p[0]->[1] = -$p[0]->[1];
      $p[0]->[2] = -$p[0]->[2];
  }
  return @p;
}


=head2   FUNCTION findHelixAxis

 Title   :  findHelixAxis
 Usage   :  PDB::findHelixAxis(\@res, $pdb)
 Function:  Finds the helix axis.
 Returns :  Returns the helix axis as six values - three defining an origin and the next three defining the direction.
 Args    :  (1) reference to an array of residue or array of CA atoms
            (2) optional: PDB structure. If defied, then the axis points used in building the axis will be added as
                an new chain and a new residue to the PDB structure.
            (3) optional: reference to a hash, where additional information will be returned.

=cut
sub findHelixAxis {
  my $resArr = shift;
  GENERAL::requireArgs($resArr);
  my $pdb = shift;
  my $ret = shift;
  my $R = 2.26;

  my (@ca, @x, @y, @z, $nch, $nres);
  foreach my $res (@$resArr) {
    if (defined($res->{xcoor}) && defined($res->{ycoor}) && defined($res->{zcoor})) {
      push(@ca, $res);
    } else {
      push(@ca, PDB::getAtomInRes($res, "CA"));
    }
  }
  GENERAL::assert(scalar(@ca) >= 3, "helix has to be at least three residues long");
  if (defined($pdb)) {
    $nch = $pdb->PDB::_newChain("");
    $nres = PDB::_newResidue("AXS", 1, 1, 1, $nch);
  }
  for (my $i = 1; $i < scalar(@ca)-1; $i++) {
    my @r1 = ($ca[$i-1]->{xcoor} - $ca[$i]->{xcoor}, $ca[$i-1]->{ycoor} - $ca[$i]->{ycoor}, $ca[$i-1]->{zcoor} - $ca[$i]->{zcoor});
    my @r2 = ($ca[$i+1]->{xcoor} - $ca[$i]->{xcoor}, $ca[$i+1]->{ycoor} - $ca[$i]->{ycoor}, $ca[$i+1]->{zcoor} - $ca[$i]->{zcoor});
    my @r = ($r1[0] + $r2[0], $r1[1] + $r2[1], $r1[2] + $r2[2]);
    my $rl = sqrt($r[0]**2 + $r[1]**2 + $r[2]**2);
    $r[0] /= $rl; $r[1] /= $rl; $r[2] /= $rl;

    my $nx = $ca[$i]->{xcoor} + $r[0]*$R;
    my $ny = $ca[$i]->{ycoor} + $r[1]*$R;
    my $nz = $ca[$i]->{zcoor} + $r[2]*$R;
    push(@x, $nx); push(@y, $ny); push(@z, $nz);

    my $na = PDB::_newAtom($i, $i, (($i == 1) ? "AXN" : (($i == scalar(@ca)-2) ? "AXC" : "A$i")), $nx, $ny, $nz, 0.0, 1.0, 0.0, $nres) if (defined($pdb));
  }

  if (defined($ret)) { $ret->{xcoor} = \@x; $ret->{ycoor} = \@y; $ret->{zcoor} = \@z; }

  return GENERAL::leastSquares3D(\@x, \@y, \@z);
}


=head2   FUNCTION superimpose

 Title   :  superimpose
 Usage   :  PDB::superimpose(\@A, \@B)
 Function:  Optimally super-imposes first atom array onto the second one and reports the RMSD.
            Optionally, also updates the coordinates of the first array to superimpose onto the second.
            Also, if the move flag is an array of atoms, that array of atoms is moved, rather than the first array.
 Returns :  RMSD of the best-fit superposition
 Args    :  (1) first array of atoms
            (2) second array of atoms
            {3} move flag - 0, 1 or an array of atoms, or a PDB structure (default 0)

=cut
sub superimpose {
  my $a = shift;
  my $b = shift;
  my $move = shift;
  my $N = scalar(@$a);
  GENERAL::assert(scalar(@$b) == $N, "the two atom vectors must be of the same length!");

  if ($N == 1) { # trivial case
    if (defined($move)) {
      GENERAL::assert((ref($move) ne "ARRAY") && (ref($move) ne "PDB"), "The super-position of one atom onto another one is too ambiguous to extract a transformation from it for applying to other atoms.");
      for (my $i = 0; $i < $N; $i++) {
        $a->[$i]->{xcoor} = $b->[$i]->{xcoor};
        $a->[$i]->{ycoor} = $b->[$i]->{ycoor};
        $a->[$i]->{zcoor} = $b->[$i]->{zcoor};
      }
    }
    return 0;
  }

  # copy atom arrays to 3xN matrices
  my $A = zeroes($N, 3); my $B = zeroes($N, 3);
  for (my $i = 0; $i < $N; $i++) {
    $A->slice("$i,0:2") += pdl([$a->[$i]->{xcoor}], [$a->[$i]->{ycoor}], [$a->[$i]->{zcoor}]);
    $B->slice("$i,0:2") += pdl([$b->[$i]->{xcoor}], [$b->[$i]->{ycoor}], [$b->[$i]->{zcoor}]);
  }

  # center
  my $cA = transpose(sumover($A)/$N);
  for (my $i = 0; $i < $N; $i++) {
    $A->slice("$i,0:2") -= $cA;
  }
  my $cB = transpose(sumover($B)/$N);
  for (my $i = 0; $i < $N; $i++) {
    $B->slice("$i,0:2") -= $cB;
  }

  # correlation matrix, SVD, etc.
  my $R = $A x transpose($B);
  my ($V, $S, $W) = svd($R);
  my $I = identity(3);
  my $sgn = (det($R) > 0) ? 1 : ((det($R) < 0) ? -1 : 0);
  $I->slice("2,2") .= $sgn;
  my $M = $W x $I x transpose($V); # optimal alignment matrix
  my $rmsd = sqrt(sum(($B - $M x $A)**2)/$N);

  if ($move) {
    if (ref($move) eq "PDB") {
      my @arr = $move->conAtoms();
      $move = \@arr;
    }
    # If move is set, either move A, or, if move is an array of atoms, replace A with the contents of this array and move these atoms instead.
    # In the latter case, first need to apply the same translation that was applied to A
    if (ref($move) eq "ARRAY") {
      $N = scalar(@$move);
      # copy (replace the A matrix with whatever needs to be moved)
      my $m = zeroes($N, 3);
      for (my $i = 0; $i < $N; $i++) {
        $m->slice("$i,0:2") += pdl([$move->[$i]->{xcoor}], [$move->[$i]->{ycoor}], [$move->[$i]->{zcoor}]);
      }

      # apply the same move that was applied to original A
      for (my $i = 0; $i < $N; $i++) {
        $m->slice("$i,0:2") -= $cA;
      }

      # rotate A to optimally align with origin-centered B, then move its center into the center of original B
      $m = $M x $m;
      for (my $i = 0; $i < $N; $i++) {
        $move->[$i]->{xcoor} = $m->at($i, 0) + $cB->at(0, 0);
        $move->[$i]->{ycoor} = $m->at($i, 1) + $cB->at(0, 1);
        $move->[$i]->{zcoor} = $m->at($i, 2) + $cB->at(0, 2);
      }
    } else {
      # rotate A to optimally align with origin-centered B, then move its center into the center of original B
      $A = $M x $A;
      for (my $i = 0; $i < $N; $i++) {
        $a->[$i]->{xcoor} = $A->at($i, 0) + $cB->at(0, 0);
        $a->[$i]->{ycoor} = $A->at($i, 1) + $cB->at(0, 1);
        $a->[$i]->{zcoor} = $A->at($i, 2) + $cB->at(0, 2);
      }
    }
  }
  return $rmsd;
}

=head2   FUNCTION TMscore

 Title   :  TMscore
 Usage   :  my $score = PDB::TMscore($pdbFile1, $pdbFile2)
            my $score = PDB::TMscore($pdbFile1, $pdbFile2, $pdb1)
            my $score = PDB::TMscore($pdbFile1, $pdbFile2, \@atoms)
 Function:  Computes the TM score between two structures of the same protein. For example, a model versus
            the native structure (if this is the case, the model should come first and then the native --
            that's how it is specified in the help message of the TMscore program).
            Optionally, also applies the alignment (that takes superimposes the first structure onto the
            second one) to the specified set of atoms or structure.
 Returns :  The TM score.
 Args    :  (1) first structure (e.g. the model): name of a PDB file.
            (2) second structure (e.g. the native): name of a PDB file.
            {3} PDB structure or a reference to an array of atoms to apply the transformation to (the transformation
            for going from the first structure to the second).

=cut

sub TMscore {
  my $pdbf1 = shift;
  my $pdbf2 = shift;
  my $applyTo = shift;
  my $line;

  # run TMscore and parse output
  my $out = `$BIN_DEF/TMscore $pdbf1 $pdbf2`;
  my @out = split("\n", GENERAL::Trim($out));
  do { $line = shift @out; } while (defined($line) && ($line !~ /^TM-score\s+=\s+(\S+)\s/));
  GENERAL::assert(defined($line), "could not parse TMscore output (1)");
  my $TMscore = $1;
  GENERAL::assert(GENERAL::isNumeric($TMscore), "could not parse TMscore output (2)");
  do { $line = shift @out; } while (defined($line) && ($line !~ /rotation matrix to rotate Chain-1/));
  GENERAL::assert(defined($line), "could not parse TMscore output (3)");
  shift @out;
  my $m1 = shift @out;
  my $m2 = shift @out;
  my $m3 = shift @out;
  GENERAL::assert(defined($m1) && defined($m2) & defined($m3), "could not parse TMscore output (4)");
  my @m1 = split(" ", $m1);
  my @m2 = split(" ", $m2);
  my @m3 = split(" ", $m3);
  GENERAL::assert((scalar(@m1) == 5) && (scalar(@m2) == 5) && (scalar(@m3) == 5), "could not parse TMscore output (5)");
  GENERAL::assert(($m1[0] == 1) && ($m2[0] == 2) && ($m3[0] == 3), "could not parse TMscore output (5)");
  my $M = PDB::matTransX($m1[1]) x PDB::matTransY($m2[1]) x PDB::matTransZ($m3[1]) x pdl([$m1[2], $m1[3], $m1[4], 0], [$m2[2], $m2[3], $m2[4], 0], [$m3[2], $m3[3], $m3[4], 0], [0, 0, 0, 1]);

  # apply transformation, if asked
  PDB::applyTransformation($applyTo, $M) if (defined($applyTo));

  return ($TMscore, $M);
}

=head2   FUNCTION waterTIP3

 Title   :  waterTIP3
 Usage   :  my $pdb = PDB::waterTIP3()
 Function:  Creates a PDB structure corresponding to TIP3 water
 Returns :  PDB structure of TIP3 water
 Args    :

=cut

sub waterTIP3 {
  my $pdb = PDB::new();

  my $c = $pdb->_newChain("");
  my $res = PDB::_newResidue("TIP3", 1, 1, 1, $c, 1);
  PDB::_newAtom(1, 1, "OH2", 0.340, -0.222,  0.027, 0.00, 1.00, 0.00, $res);
  PDB::_newAtom(2, 2, "H1", -0.104,  0.147, -0.776, 0.00, 1.00, 0.00, $res);
  PDB::_newAtom(3, 3, "H2", -0.236,  0.075,  0.750, 0.00, 1.00, 0.00, $res);

  return $pdb;
}

1;
