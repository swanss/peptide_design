# Sequence package
# read/write/convert sequence information
#

package SEQUENCE;

use strict;

#use FileHandle;
use GENERAL;
use PDB;

## data: sequence[] -> { residue secondary index valid }
## sequence data

## data: resinx[]
## lookup table for residue index


use vars qw ( %_sectrans %_seqabbrev %_seqlong %_secnum %_secletter %_seqindx %_indxseq);


BEGIN {
  %_sectrans  = ( ALPHA   => 'H', A => 'H',
                  BETA    => 'E', B => 'E', SHEET => 'E', EXTENDED => 'E',
                  COIL    => 'C',
                  OTHER   => 'U', UNKNOWN => 'U' );

  %_secnum    = ( H => 2, E => 4, C => 1, U => 1 ); 

  %_secletter = ( 1 => 'U', 2 => 'H', 4 => 'E' );

  %_seqabbrev = ('GLY', 'G', 'PRO', 'P', 'ALA', 'A', 'VAL', 'V',
                 'LEU', 'L', 'ILE', 'I', 'MET', 'M', 'CYS', 'C',
                 'PHE', 'F', 'TYR', 'Y', 'TRP', 'W', 'HIS', 'H',
                 'HSD', 'H', 'HSE', 'H', 'HSC', 'H', 'HSP', 'H',
                 'LYS', 'K', 'ARG', 'R', 'GLN', 'Q', 'ASN', 'N', 
                 'GLU', 'E', 'ASP', 'D', 'SER', 'S', 'THR', 'T', 'MSE', 'M');

  %_seqlong   = ('G', 'GLY', 'P', 'PRO', 'A', 'ALA', 'V', 'VAL', 
                 'L', 'LEU', 'I', 'ILE', 'M', 'MET', 'C', 'CYS',
                 'F', 'PHE', 'Y', 'TYR', 'W', 'TRP', 'H', 'HIS',
                 'K', 'LYS', 'R', 'ARG', 'Q', 'GLN', 'N', 'ASN',
                 'E', 'GLU', 'D', 'ASP', 'S', 'SER', 'T', 'THR');
  
  %_seqindx   = ('A', 0, 'C', 1, 'D', 2, 'E', 3, 'F', 4, 'G', 5,
                 'H', 6, 'I', 7, 'K', 8, 'L', 9, 'M',10, 'N',11,
                 'P',12, 'Q',13, 'R',14, 'S',15, 'T',16, 'V',17,
                 'W',18, 'Y',19);
  
  %_indxseq   = (0, 'A', 1, 'C', 2, 'D', 3, 'E', 4, 'F', 5, 'G',
                 6, 'H', 7, 'I', 8, 'K', 9, 'L',10, 'M',11, 'N',
                 12,'P',13, 'Q',14, 'R',15, 'S',16, 'T',17, 'V',
                 18,'W',19, 'Y');

}


## constructor: new([seqstring|molecule][,fraglist])
## creates a new Sequence object. An abbreviated
## sequence string or a Molecule object can be
## given as arguments to extract the sequence.
## A fragment list can be given as a second argument
## to complete fragments missing in a molecule structure
## from abbreviated sequence strings.

sub new {

    my $arg=shift;
    my $slist=shift;

    my $self={};
    
    if (!defined $arg) {
	$self->{sequence}=();
    } elsif (!ref($arg)) {
	$arg=~s/[ \n\t]+//;
	for (my $i=0; $i<length($arg); $i++) {
	    my $seqrec={};
	    my $char=uc substr($arg,$i,1);
	    die "unknown residue name abbreviation $char found" 
		if (!exists $_seqlong{$char});
	    $seqrec->{index}=$i+1;
	    $seqrec->{residue}=$_seqlong{$char};
	    $seqrec->{secondary}='U';
	    $seqrec->{valid}=1;
	    push(@{$self->{sequence}},$seqrec);
	}
    } else {
	
	my $c=$arg->activeChains()->[0];
	
	my $ts=();
	foreach my $r ( @{$c->{res}} ) {
	    my $tsrec={};
	    $tsrec->{index}=$r->{num};
	    $tsrec->{residue}=$r->{name};
	    
	    push(@{$ts},$tsrec);
	}
	
	foreach my $s ( @{$slist} ) {
	    for (my $i=0; $i<length($s->{seq}); $i++) {
		my $char=uc substr($s->{seq},$i,1);
		die "unknown residue name abbreviation $char found" 
		    if (!exists $_seqlong{$char});
		my $tsrec={};
		$tsrec->{index}=$i+$s->{inx};
		$tsrec->{residue}=$_seqlong{$char};
		push(@{$ts},$tsrec);
	    }
	}
	
	my $lastinx;
	foreach my $t ( sort {$a->{index}<=>$b->{index}} @{$ts} ) {
	    my $seqrec={};
	    $seqrec->{index}=$t->{index};
	    $seqrec->{residue}=$t->{residue};
	    $seqrec->{secondary}='U';
	    $seqrec->{valid}=1;
      
	    printf STDERR "warning: non-continuous sequence\n"
		if (defined $lastinx && $seqrec->{index}!=$lastinx+1);
	    
	    push(@{$self->{sequence}},$seqrec);
	    $lastinx=$seqrec->{index};
	}
    }
    
    bless($self);
    return $self;
}


## method: $string = abbrevSeq()
## generates an abbreviated sequence string

sub abbrevSeq {
    my $self=shift;
    my $tstr="";
    
    for (my $i=0; $i<=$#{$self->{sequence}}; $i++) {
	my $resname=uc $self->{sequence}->[$i]->{residue};
	die "unknown residue name" 
	    if (!exists $_seqabbrev{$resname});
	$tstr.=$_seqabbrev{$resname};
    }
    return $tstr;
}

## method: $string = abbrevSec()
## generates an abbreviated string for the secondary
## structure information

sub abbrevSec {
    my $self=shift;
    my $tstr="";
    
    for (my $i=0; $i<=$#{$self->{sequence}}; $i++) {
	$tstr.=$self->{sequence}->[$i]->{secondary};
    }
    return $tstr;
}

## method: setSecondary(index,type)
## sets the secondary structure information
## for a specific residue

sub setSecondary {
    my $self=shift;
    my $num=shift;
    my $type=uc shift;

    my $sec = (exists $_secnum{$type}) ? $type : $_sectrans{$type};
    
    die "invalid secondary structure specification" 
	if (!defined $sec || $sec eq "");
    
  die "index out of range" 
      if (!defined $self->{sequence}->[$num]);
    
    $self->{sequence}->[$num]->{secondary}=$sec;
}

## method: secFromPredict(filelist)
## extracts secondary structure information from
## a list of file names containing output from
## common prediction servers

sub secFromPredict {
    my $self=shift;
    my @fileArr=@_;
    
    my @allsec=();
    foreach my $filename (@fileArr) {
	my $sec=&_readPredictOutput($filename,$self->{sequence});
	push (@allsec,$sec);
    }
    
    for (my $i=0; $i<=$#{$self->{sequence}}; $i++) {
	my $cntE=0;
	my $cntH=0;
	my $cntU=0;
	
    for (my $j=0; $j<=$#allsec; $j++) {
	my $c=substr($allsec[$j],$i,1);
	if ($c eq "E") {
	    $cntE++;
	} elsif ($c eq "H") {
	    $cntH++;
	} else {
	    $cntU++;
	}
    }
	
	if ($cntE>$cntH && $cntE>=$cntU) {
	    $self->{sequence}->[$i]->{secondary}="E";
	} elsif ($cntH>$cntE && $cntH>=$cntU) {
	    $self->{sequence}->[$i]->{secondary}="H";
	} else {
	    $self->{sequence}->[$i]->{secondary}="U";
	}
    }
}

## method: secFromDSSP(molecule)
## sets secondary structure information from the
## output of the external DSSP program for the
## protein structure given as Molecule object argument 

sub secFromStride {

    my $self=shift;
    my $mol=shift;

    use Sys::Hostname;

    &GENERAL::SetLogFile("stride.log");
    &GENERAL::Log("Sequence::secFromStride");
    
    my $stridebin=&GENERAL::FindExecutable("stride");
    die "Cannot find stride executable"
	if (!defined $stridebin);

    my $tempfile = hostname."temp".".pdb";
    my $pdbfile = &GENERAL::GetOutFH($tempfile);
    $mol->writePDB($pdbfile,"GENERIC");
   
    system "$stridebin $pdbfile > stride.out";
 
    my $strideout = GENERAL::GetInFH("stride.out");

    my $readdata=0;
    while(<$strideout>) {
	if ($readdata) {
	    my $num=substr($_,0,5)+0;
	    my $ab=uc substr($_,13,1);
	    my $sec=uc substr($_,16,1);
      
	    my $eab=$_seqabbrev{$self->{sequence}->[$num-1]->{residue}};
	    die "Residues do not match (dssp: $ab, expected: $eab)"
		if ( $eab ne $ab);
	    
	    $self->{sequence}->[$num-1]->{secondary}=($sec eq "H")?"H":(($sec eq "E")?"E":"U");
	    
	} elsif (/^ +\#  RESIDUE AA.*/) {
	    $readdata=1;
	}
    }
   close $strideout;

}

## method: number = firstResNum()
## number of first residue

sub firstResNum {
    my $self=shift;
    return $self->{sequence}->[0]->{index};
}

## method: number = lastResNum()
## number of last residue

sub lastResNum {
    my $self=shift;
    return $self->{sequence}->[$#{$self->{sequence}}]->{index};
}

## method: setValidResidues(fraglist)
## sets the <mark>valid</mark> flag to 1
## for residues in the fragment list and to
## 0 for residues outside the list

sub setValidResidues {
    my $self=shift;
    my $fraglist=shift;
    my $exclmode=shift;

    $exclmode=0 if (!defined $exclmode);
    
    $self->resetValidResidues($exclmode?1:0);
    
    foreach my $l ( @{$fraglist} ) {
    if (!defined $l->{from}) {
	$l->{from}=$self->firstResNum();
	$l->{to}=$self->lastResNum();
    } elsif (!defined $l->{to}) {
	$l->{to}=$l->{from};
    }
    for (my $i=$l->{from}; $i<=$l->{to}; $i++) {
      my $r=$self->getResidue($i);
      $r->{valid}=($exclmode?0:1) if (defined $r);
  }    
  }
}

## method: resetValidResidues([value])
## sets the <mark>valid</mark> flag of all
## residues to the given value (default: 1)

sub resetValidResidues {
    my $self=shift;
    my $value=shift;
    
    $value=1 if (!defined $value);
    
    foreach my $r ( @{$self->{sequence}} ) {
    $r->{valid}=$value;
}
}

## method: $list = listFromValid(force)
## returns a list of residues from the residues 
## previously set with <mark>setValidResidues</mark>.

sub listFromValid {
    my $self=shift;
    
    my $retlist=();
    my @arr;
    foreach my $r ( @{$self->{sequence}} ) {
	push(@arr,$r->{index}) if ($r->{valid});
    }
    if ($#arr>=0) {
	foreach my $tlist ( @{&GenUtil::fragListFromArray(\@arr)} ) {
	    push(@{$retlist},$tlist);
	}
    } else {
    my $rec={};
    $rec->{from}=$self->firstResNum();
    $rec->{to}=$self->lastResNum();
    push(@{$retlist},$rec);
}
    
    return $retlist;
}

## method: $index = getResidue(resnum)
## returns the residue index for a residue number. 
## A reference to a chain structure may be given as 
## additional argument for multi-domain structures

sub getResidue {
    my $self=shift;
    my $inx=shift;
    
    if (!defined $self->{resinx}) {
    $self->{resinx}={};
    for (my $ir=0; $ir<=$#{$self->{sequence}}; $ir++) {
	my $key=$self->{sequence}->[$ir]->{index};
	$self->{resinx}->{$key}=$ir;
    }
}
    
    return $self->{sequence}->[$self->{resinx}->{$inx}];
}  


sub _readPredictOutput {
    my $filename=shift;
    my $seq=shift;
    my $secstr="";
    
    die "cannot read prediction output filename $filename" 
	if (!-r $filename);
    
    my $current=0;
    my $sspro=0;
    my $phd=0;
    
    open PREDINP,"$filename";
  PREDIO:
    while (<PREDINP>) {
	chomp;
	
	my $sa="";
	my $sec="";
	
	if (/SSpro prediction/) {
	    $sspro=1;
	} elsif ($sspro==1 && /^Prediction:/) {
	    $sspro=2;
	} elsif ($sspro==2 && /^[A-Z]+/) {                         # sspro
	    ($sa=$_)=~s/ +//;
	    $sec=<PREDINP>;
	    chomp $sec;
	    $sec=~s/ +//;
	} elsif (/PHD output/) {
	    $phd=1;
	} elsif ($phd==1 && /^[ \t]+protein:[ \t]+predict[ \t]+length/) {
	    $phd=2;
	} elsif ($phd==2 && /^[ \t]+AA[ \t]+\|/) {                 # phd
	    ($sa=$_)=~s/^[ \t]+AA[ \t]+\|//;
	    $sa=~s/\|.*$//;
	    $sec=<PREDINP>;
	    chomp $sec;
	    $sec=~s/^[ \t]+PHD[ \ta-z]+\|//;
	    $sec=~s/\|.*$//;
	} elsif (/^Pred: [CEH]+/) {                                # psipred
	    ($sec=$_)=~s/^Pred: //;
	    $sec=~s/ +//;
	    $sa=<PREDINP>;
	    chomp $sa;
	    $sa=~s/^  AA: //;
	    $sa=~s/ +//;
	} elsif (/^[ \t]*[A-Z]+[ \t]+[CEHU][ \t]+[01]\.[0-9]/) {   # jpred2, pssp, prof2, pred2ary
	    s/^[ \t]+//;
	    my @f=split(/[ \t]+/);
	    
	    $sa=(exists $_seqabbrev{$f[0]}) ? $_seqabbrev{$f[0]} : $f[0];
	    $sec=$f[1];
	} elsif (/^ *[0-9]+ [A-Z] . [HE-] [0-9]/) {                # pred2ary
	    s/^[ \t]+//;
	    my @f=split(/[ \t]+/);
	    
	    $sa=(exists $_seqabbrev{$f[1]}) ? $_seqabbrev{$f[1]} : $f[1];
	    $sec=$f[3];
	}
	
	for (my $i=0; $i<length($sa); $i++) {
	    my $sai=substr($sa,$i,1);
	    my $seci=substr($sec,$i,1);
	    if ($sai eq $_seqabbrev{$seq->[$current]->{residue}}) {
		$secstr.=$seci;
		last PREDIO if (++$current>$#{$seq}); 
	    }
	}
    }
    
    die "did not read secondary prediction for complete sequence in $filename"
	if ($current<=$#{$seq});
    
    close PREDINP;
    
    return $secstr;
}



# New functions added by Gevorg starting 11/09/2004

=head2 

 Title   : allSingles
 Usage   : SEQUENCE::allSingles()
 Function: Returns an array of 20 natural amino acids.
 Returns : 
 Args    : 

=cut

sub allSingles {
  return keys(%_seqlong);
}

=head2 

 Title   : s2t
 Usage   : SEQUENCE::s2t($sequence)
 Function: Translate from single letter to triple letter amino acid code
 Returns : Translated triple letter amino acid
 Args    : 1. A character representing an single letter amino acid.

=cut

sub s2t {
  my $s = shift;
  my $strict = shift; $strict = 1 if (!defined($strict));

  if (!defined($_seqlong{$s})) {
    if ($strict eq 1) {
      GENERAL::error("Could not find the three letter code for amino acid \"$s\"");
    } elsif ($strict eq 0) {
      GENERAL::warning("Could not find the three letter code for amino acid \"$s\"");
      return -1;
    } else {
      return $strict;
    }
  }
  return $_seqlong{$s};
}

=head2 

 Title   : t2s
 Usage   : SEQUENCE::t2s($sequence)
 Function: Translate from triple letter to single letter amino acid code
 Returns : Translated single letter amino acid
 Args    : 1. A string representing a triple letter amino acid.

=cut

sub t2s {
  my $t = shift;
  my $strict = shift; $strict = 1 if (!defined($strict));

  if (!defined($_seqabbrev{$t})) {
    if ($strict eq 1) {
      GENERAL::error("Could not find the single letter code for amino acid \"$t\"");
    } elsif ($strict eq 0) {
      GENERAL::warning("Could not find the single letter code for amino acid \"$t\"");
      return -1;
    } else {
      return $strict;
    }
  }
  return $_seqabbrev{$t};
}

=head2 

 Title   : aa2index
 Usage   : SEQUENCE::aa2index($aa)
 Function: Translate from single or tripple letter code amino acid into a unique index.
 Returns : 
 Args    : 1. Amino acid name (1 or 3 characters)
           2. Strictness flag. 1 (default) - error on unknown aa; 0 - warning on unknown aa and return -1;
              -1 or other negative - return the flag on unknown aa.

=cut

sub aa2index {
  my $aa = shift;
  my $strict = shift;
  $strict = 1 if (!defined($strict));

  if (length($aa) == 3) {
    if (!defined($_seqabbrev{$aa}) || !defined($_seqindx{$_seqabbrev{$aa}})) {
      GENERAL::error("Unknown amino acid $aa!") if ($strict > 0);
      if ($strict == 0) { GENERAL::warning("Unknown amino acid $aa!"); return -1; }
      return $strict;
    }
    return $_seqindx{$_seqabbrev{$aa}};
  } elsif (length($aa) == 1) {
    if (!defined($_seqindx{$aa})) {
      GENERAL::error("Unknown amino acid $aa!") if ($strict > 0);
      if ($strict == 0) { GENERAL::warning("Unknown amino acid $aa!"); return -1; }
      return $strict;
    }
    return $_seqindx{$aa};
  } else {
    GENERAL::error("Amino acid name '$aa' is neither in single nor triple code!") if ($strict > 0);
    if ($strict == 0) { GENERAL::error("Amino acid name '$aa' is neither in single nor triple code!"); return -1; }
    return $strict;
  }
}

=head2 

 Title   : index2s
 Usage   : SEQUENCE::index2s($i)
 Function: Converts from amino acid index to single code.
 Returns : 
 Args    : 1. Amino acid index.

=cut

sub index2s {
  my $i = shift;
  if (!defined($_indxseq{$i})) {
    GENERAL::error("Unknown amino acid index $i!");
  }
  return $_indxseq{$i};
}

=head2 

 Title   : index2t
 Usage   : SEQUENCE::index2t($i)
 Function: Converts from amino acid index to triplet code.
 Returns : 
 Args    : 1. Amino acid index.

=cut

sub index2t {
  my $i = shift;
  return SEQUENCE::s2t(SEQUENCE::index2s($i));
}


1;



