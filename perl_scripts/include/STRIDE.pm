#
# Design module for Stride
#
#
# POD documentation - main docs before the code
=head1 NAME

STRIDE.pm - Module for running/parsing/accessing stride output


=head1 SYNOPSIS

 my $StrideObj = STRIDE::new("blah.pdb" );
 my $StrideObj = STRIDE::new('blah.stride' );


=head1 DESCRIPTION

STRIDE is a module for objectifying STRIDE output.  STRIDE is a
program (similar to DSSP) for assigning secondary structure to
individual residues of a pdb structure file.

    ( Knowledge-Based Protein Secondary Structure Assignment,
    PROTEINS: Structure, Function, and Genetics 23:566-579 (1995) )

STRIDE is available here:
http://www.embl-heidelberg.de/argos/stride/down_stride.html

Methods are then available for extracting all of the infomation
present within the output or convenient subsets of it.

Although they are very similar in function, DSSP and STRIDE differ
somewhat in output format.  Thes differences are reflected in the
return value of some methods of these modules.  For example, both
the STRIDE and DSSP parsers have resSecStr() methods for returning
the secondary structure of a given residue.  However, the range of
return values for DSSP is ( H, B, E, G, I, T, and S ) whereas the
range of values for STRIDE is ( H, G, I, E, B, b, T, and C ).  See
individual methods for details.

=head1 AUTHOR - Jiangang Chen

Email chen2001@mit.edu

=head1 APPENDIX

The Rest of the documentation details each method.
Internal methods are preceded with a _.

=cut




package STRIDE;
#$DEF_VER_CONTROL->{ver}->{ROTAMER} = "1.0";
$DEF_VER_CONTROL->{date}->{STRIDE} = "06/17/04";

use GENERAL;
use DEFINITIONS;


=head2 

 Title   :  new
 Usage   :  new(PDBfile)     pdb file
 Function:  creates a new Stride object and reads a PDB if a file name is given
 Example :    
 Returns :  stride object
 Args    :  a pdb file name or stride output is optional 
            however, pdb file must end with pdb or ent
            stride output must end up with stride
=cut




sub new {
    
  my $self={};
  my $fname=shift;

  $self->{'interface'}=undef;
  $self->{'totalarea'}= 0.0;
  $self->{'stride'}=undef;
  $self->{'pdb'}=undef;
  $self->{'ASG'} = {};
  $self->{'LOC'} = [];

  bless $self;

  # here we create a object by
  #      either reading in a pdb file
  #          or reading in a stride output
  if (defined $fname &&  $fname =~ /\.pdb|\.ent/ ) {
    $self->{'pdb'}=$fname;
    # run stride calculations
    $self->runStride();

    # get stride ouput
    my $prefix = GENERAL::GetBase($fname);
    my $strideout = $prefix.'.stride';
    $self->readStride($strideout);

  } elsif  (defined $fname &&  $fname =~ /\.stride/ ) {
    # if stride output is provided, parse the output
    $self->readStride($fname);
  }

  return $self;
}


=head2 

 Title   :  readStride
 Usage   :  $strideobj->readStride('blah.stride') 
 Function:  parse the stride output file into a data structure
 Returns : 
 Args    :  a stride output file name 

=cut


sub readStride {
    
    my $self=shift;
    my $inp=shift;
    my $fname = GENERAL::GetInFH($inp);

    # initialize the self data structure

    $self->{'ASG'}={};
    $self->{'stride'} = $inp;
    $self->{'LOC'}=[];

 
    while (<$fname>) {

	# here we speed up the line reading 
	# by using next if techniques

	next if ( /^REM/ );
	next if ( /^HDR/);
	next if ( /^SEQ/);
	next if ( /^STR/);
	
	# Here we parsed LOC record

	if ( /^LOC/ ) {
	    
	    substr($_,38, 1) = 'A' if  (substr($_,38, 1) eq " ");

	    my @line = split(" ",$_);
	    my $rec;

	    # we record the chain id, secondary structure and start, end position

	    $rec->{secdetail}= $line[1];

	    # we don't want to collect the disulfide bond information

      next if ( $rec->{secdetail} =~ /Disulfide/ ); 
	    $rec->{secstr}='T' if ( $rec->{secdetail} =~ /Turn/);
	    $rec->{secstr}='H' if ( $rec->{secdetail} =~ /AlphaHelix/);
	    $rec->{secstr}='G' if ( $rec->{secdetail} =~ /310Helix/);
	    $rec->{secstr}='E' if ( $rec->{secdetail} =~ /Strand/);
	    $rec->{secstr}='T' if ( $rec->{secdetail} =~ /Gamma/);
	    $rec->{chain}=$line[4];
	    $rec->{start}=$line[3];
	    $rec->{end}=$line[6];
	    my $label=$rec->{secstr}."_".$rec->{chain}."_".$rec->{start}."_". $rec->{end};
	    $rec->{label}=$label;

	    push (@{$self->{'LOC'}}, $rec);

	}
       
	if ( /^ASG/ ) {

	    # incase we have empty field, a single chain domain... 
	    # we assign 'A' to that domain.

	    substr($_, 9, 1)='A' if ( substr($_, 9, 1) eq "-");
	    my @line=split(" ", $_);
	    my $chain=$line[2];
	    
	    my $asgrec={};
    
      # now parse the input into hash
      $asgrec->{resname}=$line[1];
      $asgrec->{chain} =$line[2];
      $asgrec->{resnum}=$line[3];
      $asgrec->{iresnum}=$line[4]+0;
      $asgrec->{secstr}=$line[5];
      $asgrec->{secdetail} =$line[6];
      $asgrec->{phi}=$line[7]+0.0;
      $asgrec->{phi}=0.0  if ($asgrec->{phi} == 360.00);
      $asgrec->{psi}=$line[8]+0.0;
      $asgrec->{psi}=0.0  if ($asgrec->{psi} == 360.00);
      $asgrec->{area}=$line[9]+0.0;

      my $i = GENERAL::RoundIt($asgrec->{phi});
      my $j = GENERAL::RoundIt($asgrec->{psi});
      my $index = ($i/10+18)*37 + ($j/10+19);
      $asgrec->{bbindex} = $index; 

	    # store the record into array
	    push (@{$self->{'ASG'}->{$chain}}, $asgrec);
	}	
    }


}



=head2

 Title   :  getResname
 Usage   :  $self->getResname($chain, $iresenum)
 Function:  obtain residue's name
 Returns :  residue name
 Args    :  $chain: chain ID
            $iresnum: indexed residue number

=cut

sub getResname {

    my $self=shift;
    my $chain=shift;
    my $iresnum=shift;
    my $index=$iresnum-1;
    my $resname=$self->{'ASG'}->{$chain}->[$index]->{resname};

    return $resname;
}


=head2

 Title   :  getResSecStr
 Usage   :  $self->getResSecStr($chain, $iresenum)
 Function:  obtain residue's secondary structure
 Returns :  residue secondary structure
 Args    :  $chain: chain ID
            $iresnum: indexed residue number

=cut

sub getResSecStr {

    my $self=shift;
    my $chain=shift;
    my $iresnum=shift;
    my $index=$iresnum-1;
    my $secstr = $self->{'ASG'}->{$chain}->[$index]->{secstr};

    return $secstr;
}



=head2

 Title   :  getCustomizedInp
 Usage   :
 Function:
 Returns :
 Args    :

=cut

sub getCustomizedInp  {

    my $self=shift;
    my $fname=shift;
    my $outf=shift;
    my $option=shift;

    if (!defined $outf) { die "Error in getCustomizedInp: output file not defined!\n"; }
    if (!defined $option) { die "Error in getCustomizedInp: bbdep/bbindep not defined!\n"; }

    my $inpfh = GENERAL::GetInFH($fname);
    my $outfh = GENERAL::GetOutFH($outf);

    while (<$inpfh>) {
        chomp ($_);
        next if ($_ eq "");
        my $chain = substr($_, 0,1);
        my $iresnum = GENERAL::Trim(substr($_,1));
        GENERAL::assert(GENERAL::isInteger($iresnum), "string '$_' from site file could not be parsed as a chain id/residue id pair");

        if ($iresnum > $#{$self->{'ASG'}->{$chain}}+1) {
          GENERAL::error("Site $chain$iresnum from the site file has no STRIDE output (is the site present in the PDB file?).");
        } else {
           my $resname = $self->getResname($chain, $iresnum);
           my $bbindex = $self->getResPhiPsiIndex($chain, $iresnum);
           my $secstr = $self->getResSecStr($chain, $iresnum);
           my $secstrind;
           if ($secstr =~ /H/) {
             $secstrind = 1;
           } else {
             $secstrind = 3;
           }
           if (lc($option) =~ /bbdep/) {
             $outfh->printf("%1s%4d%5d%5s\n",  $chain, $iresnum, $bbindex, $resname);
           } elsif (lc($option) =~ /bbindep/) {
             $outfh->printf("%1s%4d%5s\n",  $chain, $iresnum, $resname)
           } elsif (lc($option) =~ /secstrdep/) {
             $outfh->printf("%1s%4d%5d%5s\n",  $chain, $iresnum, $secstrind, $resname);
           } else {
             die "Error in getCustomizedInp: bbdep/bbindep option ($option) not valid!\n";
           }
        }
    }

}




=head2

 Title   :  getCustomizedInp2
 Usage   :  
 Function:
 Returns :
 Args    :
 
=cut

sub getCustomizedInp2  {

    my $self=shift;
    my $fname=shift;
    my $outf=shift;
    my $option=shift;
    
    $outf='rot_input' if (!defined $outf);
    $option='bbdep' if (!defined $option);

    my $inpfh=&GENERAL::GetInFH($fname);
    my $outfh=&GENERAL::GetOutFH($outf);

    while (<$inpfh>) {

        chomp ($_);
        my @lines=split(" ", $_);
        my $chain=$lines[0];
        my $iresnum=$lines[1];

#       print $chain, $iresnum, "\n";
        if ( $iresnum > $#{$self->{'ASG'}->{$chain}}+1 ) {
            return;
        } else {
           my $resname=$self->getResname($chain, $iresnum);
           my $bbindex=$self->getResPhiPsiIndex($chain, $iresnum);
           $outfh->printf("%1s%4d%5d%5s\n",  $chain, $iresnum, $bbindex, $resname)
                    if (lc($option) =~ /bbdep/)  ;
           $outfh->printf("%1s%4d%5s\n",  $chain, $iresnum, $resname)
                    if (!$option);
        }
    }

}



=head2 

 Title   :  getResPhiPsi 
 Usage   :  $self->getResPhiPsi($chain, $iresenum)
 Function:  obtain the phi, psi values for given chain id and residue name
 Returns :  phi, psi values
 Args    :  $chain: chain ID
            $iresnum: indexed residue number

=cut



sub getResPhiPsi {

    my $self=shift;
    my $chain=shift;
    my $iresnum=shift;

    my $phi = $self->{'ASG'}->{$chain}->[$iresnum-1]->{phi};
    my $psi = $self->{'ASG'}->{$chain}->[$iresnum-1]->{psi};
    
    return ($phi, $psi);
}



=head2 

 Title   :    getResPhiPsiIndex
 Usage   :    $self->getResPhiPsiIndex($chain, $iresnum);
 Function:    obtain the backbone index for given chain id and residue name
 Returns :    backbone index value
 Args    :    $chain: chain ID
              $iresnum: indexed residue number
=cut


sub getResPhiPsiIndex {

    my $self=shift;
    my $chain=shift;
    my $iresnum=shift;

    return $self->{'ASG'}->{$chain}->[$iresnum-1]->{bbindex};
}



=head2 

 Title   :    getResArea
 Usage   :    $self->getResArea($chain, $iresnum);
 Function:    obtain the accessible area for given chain id and residue name
 Returns :    res's accessible area
 Args    :    $chain: chain ID
              $iresnum: indexed residue number
=cut


sub getResArea {

    my $self=shift;
    my $chain=shift;
    my $iresnum=shift;

    return $self->{'ASG'}->{$chain}->[$iresnum-1]->{area};
}


=head2 

THIS VERSION OF getResSecStr HAS BEEN REMOVED (IT WAS NOT USED ANYWHERE
AT THE TIME OF REMOVAL).
Gevorg 04/30/04

 Title   :    getResSecStr
 Usage   :    $self->getResSecStr($chain, $iresnum);
 Function:    obtain the secondary structure for given chain id and residue name
 Returns :    res's secondary structure or 
 Args    :    $chain: chain ID
              $iresnum: indexed residue number
              $detail flag: optional. if this argument is true, it returns the 
              detailed secondary structure description.


sub getResSecStr {

    my $self=shift;
    my $chain=shift;
    my $iresnum=shift;
    my $detail_or_not = shift;

    if ( !defined $detail_or_not) { 
      $ss_char=$self->{'ASG'}->{$chain}->[$iresnum-1]->{secstr};
      if ( $ss_char eq 'H' || $ss_char eq 'G' || $ss_char eq 'I' ) {
          return 'H';
      }
      if ( $ss_char eq 'E' || $ss_char eq 'B' || $ss_char eq 'b' ) {
          return 'B';
      }
      if ( $ss_char eq 'T' ) {
          return 'T';
      }
      else {
          return 'C';
      }
        } else {
      return $self->{'ASG'}->{$chain}->[$iresnum-1]->{secdetail};
    }
}

=cut


=head2 

 Title   :    getResSecFragment
 Usage   :    $self->getResSecFragment($chain, $iresnum);
 Function:    obtain the accessible area for given chain id and residue name
 Returns :    return the secondary structure fragment that residue belong to 
              it only gives helix, sheet and trun. no coil
 Args    :    $chain: chain ID
              $iresnum: indexed residue number
=cut


sub getResSecFragment {

     my $self=shift;
     my $chain=shift;
     my $iresnum=shift;

     foreach my $fragment (@{$self->{'LOC'}}) {
	 next if ($fragment->{chain} ne $chain);
	 if ( $iresnum >= $fragment->{start} && $iresnum <= $fragment->{end} )  {
	     return $fragment;
	 }
     }
 
     return undef;

 }


=head2 

 Title   :  getInterfaceSecFragment 
 Usage   :  
 Function: 
 Returns : 
 Args    : 

=cut

sub getInterfaceSecFragments {

    my $self=shift;
   
    die " you must first calculate interface residues. \n" 
	if (!defined $self->{'interface'});

    my %seen;
    my @fragments;

    foreach my $r (@{$self->{'interface'}}) {

	my $chain=$r->{chain};
	my $iresnum=$r->{iresnum};
	my $fragment = $self->getResSecFragment($chain, $iresnum);

	if ( defined $fragment ) {
            if (! exists $seen{$fragment->{label}}) {
	    push (@fragments, $fragment);
	    $seen{$fragment->{label}} = 1;
	    } 

        }
    }


return \@fragments;
}		   
		   	    

=head2 

 Title   :  runStride 
 Usage   :  $self->runStride()
 Function:  run stride calculations
 Returns :  None
 Args    :  None

=cut

sub runStride {

  my $self=shift;
  my $fname=$self->{pdb};
  my $prefix = GENERAL::GetBase($fname);
  my $outfnameA= $prefix.".stride";
  my $outfnameB= $prefix.".fasta";

  my $option=shift;
  $option='S' if (!defined $option);


  die "Error in runStride: you must specify the input file name for stride.\n"
  if (!defined $fname);

  if (!defined($self->{output})) {
    $self->{output} = ();
  }
  if ( $option =~ /S|s/ ) {
    GENERAL::csystem("$STRIDE_DEF $fname -f$outfnameA 1>& /dev/null", undef, 256);
    push(@{$self->{output}}, $outfnameA);
  }                                                                            
  if ( $option =~ /F|f/ ) {
    GENERAL::csystem("$STRIDE_DEF -q$outfnameB 1>& /dev/null", undef, 256);
    push(@{$self->{output}}, $outfnameB);
  } 
}

=head2

 Title   :  cleanup
 Usage   :  $self->cleanup()
 Function:  Cleans up the output of runStride
 Returns :  None
 Args    :  None

=cut

sub cleanup {
  my $self = shift;

  foreach my $ofile (@{$self->{output}}) {
    GENERAL::csystem("rm -f $ofile");
  }
}

=head2 

 Title   :  getFasta 
 Usage   :  $self->getFasta();
 Function:  write out fasta files
 Returns :  non
 Args    :  optional: chain strings for example: 'ABC' for ABC three chains or 'A'
            for 'A' chain. default is for all. 
=cut

 sub getFasta {

     my $self=shift;
     my $fname=$self->{pdb};
     die "you need to specify the pdb file through
	  setPDB method. for example: $stride->setPDB('pdb1ets.ent').\n"
	      if (!defined $fname);

     my $chainstr=shift;

     my $prefix=&GENERAL::GetBase($fname);

     if (defined $chainstr) {
	 my $chains=uc $chainstr;
	 my $outfname= $prefix.'_'.$chains.'.fasta';
	 unlink($outfname) if ( -e $outfname);
	 system("stride $fname -q$outfname 1>& /dev/null");
     } else {
	 my $chains='ALL Chains';
	 my $outfname= $prefix.'_all.fasta';
	 unlink($outfname) if ( -e $outfname);
	 system("stride $fname -q$outfname 1>& /dev/null");
     }

     print "The fasta file for $chains is $outfname.\n";
 }


=head2 

 Title   :  setPDB
 Usage   :  $self->setPDB('pdb1ent.ent');
 Function:  set the value of $self->{'pdb'}
 Returns :  
 Args    :  a string: file name

=cut


sub setPDB {

    my $self=shift;
    my $fname=shift;
    
    die "you need to give a name!\n" if (!defined $fname);
    $self->{pdb}=$fname;
    
}


=head2 

 Title   : calculateInterfaceResByNum
 Usage   : $self->calculateInterfaceResByNum('A', 'B', 5.0) 
 Function: get the reference to the interface residues between given chains 
 Returns : reference to interface residues 
 Args    : two chain string: 'A', and 'B'
           it could aso be: 'AB' and 'CDEF' which gives 
 
           one criterion: 1.0

=cut


sub calculateInterfaceResByNum {

    my $self=shift;
    my $fname=$self->{pdb};
    die "Error in calculateInterfaceResByNum: you need to specify the pdb file through setPDB method.\n"
    if (!defined $fname);

    # get chains
    my $chainA =uc(shift);
    my $chainB =uc(shift);
    my $twochains=$chainA.$chainB;
    my $criterion =shift; 
    $criterion=1.0 if (! defined $criterion);

    system("stride $fname -r$twochains -fTempAB.stride 1>& /dev/null");
    system("stride $fname -r$chainA -fTempA.stride 1>& /dev/null");
    system("stride $fname -r$chainB -fTempB.stride 1>& /dev/null");

    my $TempAB= &STRIDE::new('TempAB.stride');
    my $TempA = &STRIDE::new('TempA.stride');
    my $TempB = &STRIDE::new('TempB.stride');

    my $interfaceRes = $TempAB->_compareResArea($TempA, $TempB, $criterion);

    &GENERAL::Remove('TempAB.stride');
    &GENERAL::Remove('TempA.stride');
    &GENERAL::Remove('TempB.stride');
    
    $self->{'interface'} = $interfaceRes;
}



=head2 

 Title   :  getInterfaceRes
 Usage   :  $interfacereslist = $stride->getInterfaceRes()
 Function: 
 Returns : 
 Args    : 

=cut


sub getInterfaceRes {

    my $self=shift;
    if (defined $self->{'interface'} ) {
	return $self->{'interface'};
    } else {
	die "You need to run calculateInterfaceResByNum method first.\n";
    }

}


=head2 

 Title   :  _compareResArea
 Usage   : 
 Function:  get the interface residue
 Returns : 
 Args    : 

=cut


sub _compareResArea {

    my $self=shift;
    my $chainA=shift;
    my $chainB=shift;
    my $criterion=shift;

    my @reslist;

    my %newASG = (%{$chainA->{'ASG'}}, %{$chainB->{'ASG'}}); 

    foreach $key ( sort keys %newASG) { 
        for (my $i=0; $i <=$#{$newASG{$key}}; $i++) {
	
            my $diff= $newASG{$key}->[$i]->{area} - $self->{'ASG'}->{$key}->[$i]->{area};
	    if ($diff >= $criterion ) { 
                   $self->{'ASG'}->{$key}->[$i]->{buried}= $diff;
	          push (@reslist, $self->{'ASG'}->{$key}->[$i]);
	    }
        }
    }
    return \@reslist;
}


=head2 

 Title   :  selectInterfaceResidues 
 Usage   :  
 Function: 
 Returns : 
 Args    : 

=cut


sub selectInterfaceResidues {

    my $self=shift;
    my $chainstr=shift;

    die "you must specify the chain ids you wnat to design\n
        for example: \'A\' or \'AB\' \n" if (! defined $chainstr);

    my $interface= $self->{'interface'};
    my @selectedres;
    my $i=shift;
    my ($chain, $resnum);
    my %seen;
    my $newinterface;

    if (defined $i) {
	@selectedres=@_;
	foreach my $j (@selectedres) {
	    ($chain, $resnum) = split (/:/, $j);
	    my $label= $chain." ".$resnum;
	    $seen{$label}=1;
	}

	foreach my $j (@{$interface}) {
	    my $chain=$j->{chain};
	    my $resnum=$j->{resnum};
	    my $label= $chain." ".$resnum;
	    push (@{$newinterface}, $j) if ( $seen{$label});
	}

	$self->{'interface'}=$newinterface;
    } else {
	my @chains = split ("", $chainstr);
	foreach my $j (@chains) {
	    $seen{$j}=1;
	}
	foreach my $j (@{$interface}) {
	    my $chain=$j->{chain};
	    push (@{$newinterface}, $j) if ( $seen{$chain});
	}
    }
}




=head2

 Title   :  writeInterfaceInput
 Usage   :  writeInterfaceInput('A', 'B', 5.0);
 Example : 
         
         Problem: you want to write input for design and self, pairwise interaction calculation  
                  however, there are many choices: 
                  (A) you want to use backbone dependent library
                  (B) you want to use backbone independent library
                  (C) you have to eliminate backbone index in order to provide input file for design
                  every file has slightly different file format. 
                  (D) you also want to be able to choose the interface residues between different chains
                      in case you have a complex between multidomain proteins. for example: antibody (has 
                      heavy chain and light chain); VEGF has 4 chains.  
                  (E) you want to be able to select interface residues based on different criterions
                      sometimes, you want a tight one. sometime you want a loose one. 
        Solution: STRIDE module


        Here is an example:
 
         use vars qw ( $perllibdir );

	 BEGIN {
	   $perllibdir='/home/jchen/Module';
	 }

	 use lib $perllibdir;
	 use GENERAL;
	 use STRIDE;

         # read in pdb file name
         my $inp=shift;
         # create a new stride object
         $stride=STRIDE::new($inp);

         # Calculate the interface residues between A and B chain, 5.0 is the 
         # is the criterion: any residue being buried by 5 square angstrom 
         # will be considered as interface residues
         # Note: you can also use 'AB', 'VJKW', which will calculate the interface residue 
         # between AB chains and VJWK chains

         $stride->calculateInterfaceResByNum('A', 'B', 5.0);

         # Write the input file for self and pairwise energy calculations
         # Note: this input contains backbone index, which is suitable for backbone dependenet 
         # rotamer library. Alternatively, you can write 
         # (A) backbone independent interface input, $stride->writeInterfaceInput('test.out', 'bbindep');
         # (B) interface input file compatible with old scripts and danny's design code: 
         #      $stride->writeInterfaceInput('test.out', 'design')

         $stride->writeInterfaceInput('test.out');
         $stride->writeInterfaceReport('test2.out');

 Returns :
 Args    :

=cut


sub writeInterfaceInput {

    my $self=shift;
    my $fname=&GENERAL::GetOutFH(shift);
    my $option=shift;

    if (! defined $option) {
      foreach my $site (@{$self->{'interface'}}) {
        $fname->printf("%1s%4d%5d%5s \n", $site->{chain}, $site->{iresnum}, $site->{bbindex}, $site->{resname});
      }
    } elsif ( $option eq 'bbindep') {
      foreach my $site (@{$self->{'interface'}}) {
        $fname->printf("%1s%4d%5d%5s \n", $site->{chain}, $site->{iresnum}, 0, $site->{resname});
      }
    } elsif ( $option eq 'design') {
      foreach my $site (@{$self->{'interface'}}) {
        $fname->printf("%1s%4d%5s \n", $site->{chain}, $site->{iresnum},$site->{resname});
      }
    } else {
      print "Error in writeInterfaceInput: $option is not a valid option!\n";
    }

}




sub writeInterfaceReport {

    my $self=shift;
    my $fname=&GENERAL::GetOutFH(shift);

    $fname->print("\n");
    $fname->print("#Column 1:  Chain ID\n");
    $fname->print("#Column 2:  Residue  Number (after renumbering)\n");
    $fname->print("#Column 3:  Residue Name \n"); 
    $fname->print("#Column 4:  Backbone Index \n");
    $fname->print("#Column 5:  Secondary Structure \n");
    $fname->print("#Column 6:  Phi value \n");
    $fname->print("#Column 7:  Psi value \n");
    $fname->print("#Column 8:  Accessible Area in complex \n");
    $fname->print("#Column 9:  Buried Area During Complex \n");
    $fname->print("\n");

    foreach my $site (@{$self->{'interface'}}) {
        $fname->printf("%1s%4d%5s%5d %1s %8.3f%8.3f%8.3f%8.3f\n", $site->{chain}, $site->{iresnum},
                        $site->{resname}, $site->{bbindex}, $site->{secstr}, $site->{phi}, $site->{psi},
                        $site->{area}, $site->{buried});
    }
}






1;



