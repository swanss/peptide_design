#
# module 
#
# POD documentation - main docs before the code

=head1 NAME

DESIGN package

=head1 SYNOPSIS

=head1 DESCRIPTION

a tiny module to interface with danny's design code
(1) write optimal design input 
(2) run design code
(3) parse output

=head1 EXAMPLES

=head2 Reporting Bugs

=head1 AUTHOR 

Jiangang Chen

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a "_".

=cut

=head1 Methods
 
Here are the collections of all methods

=cut


package DESIGN;
#$DEF_VER_CONTROL->{ver}->{DESIGN} = "1.0";
$DEF_VER_CONTROL->{date}->{DESIGN} = "06/17/04";

use GENERAL;
use ROTAMER;
use PDB;
use CHARMM;



=head2 Constructor 

 Title   :  constructor for the object DESIGN 
 Usage   :  DESIGN::new()
 Function:  
 Returns : 
 Args    :  (1) method design program is taking:    dee-astar, dee-scmf, mc default is dee-astar
            (2) file name for the control input file:       
            (3) binary file name of energy data:    i.e.    energy.bin

=cut

sub new {

    my $self={};
    my $method=shift;
    my $inp=shift;
    my $ebin=shift;
    $self->{files}={};
    $self->{handles}={};
    $self->{files}->{inp}=undef;
    $self->{handles}->{'control'}=undef;
    $self->{files}->{'output'}=undef;
    $self->{handles}->{'output'}=undef;
    $self->{'solution'}={};
    $self->{solution}->{results}=[];
    $self->{solution}->{probability}=[];
    $self->{'method'}=undef;
    
    
    # define some input parameters
    # 'initial' is for mean field method.
    # when initial =1, is flat distribution. 
    # when initial=0, it is a zero distribution
    # ecut depends on the energy landscape, if the energy difference 
    # is very small, which means a very flat energy landscape. 
    # you should use small ecut. 
    # here is the suggested value for ecut in the different types of 
    # calculations
    #         (A) vdw packing: ecut=0.5 
    #         (B) vdw + solvation + electrostatics: ecut=0.05
    #         (C) vdw + solvation + electrostatics + hybrid potential: ecut=0.1
    
    my %par=( 'ecut' => 0.5 ,
              'ebin'  => undef,
              'initial' => 0,
              'cycles' =>100,
              'steps' => 10000,
              'initialT' => 4000,
              'finalT' => 150
             );
    $self->{par}=\%par;

    bless $self;
    
    my $newmethod = lc($method);
    $self->{'method'} = $newmethod;
    $self->{files}->{inp}=$inp;
    $self->{handles}->{'control'} = GENERAL::GetOutFH($inp);
    $self->{par}->{ebin} = $ebin;

    return $self;
}



=head2 

 Title   :  writeDesignInput 
 Usage   :  $design->writeDesignInput(0.02, 1)
 Function:  write design input file  
 Returns : 
 Args    :  (1) ecut for the dee-astar
            (2) initialization method for scmf

=cut

sub writeDesignInput {

    my $self=shift;
    my $secondrun=shift;
    my $ecut=shift;
    my $initial=shift;
    my $rot_sub=shift;

    $secondrun=0 if (!defined $secondrun);
    $self->{par}->{ecut}=$ecut if (defined $ecut);
    $self->{par}->{initial}=$initial if (defined $initial);
     
    GENERAL::error("You must define the control file name and method for design code.") 
    if (!defined $self->{files}->{inp} || !defined $self->{handles}->{'control'} || !defined $self->{'method'});
    
    my $method =$self->{'method'};
    if ($method eq 'dee-astar' ) {
      if (!defined $rot_sub) {
        $self->_writeDeeAstar($secondrun);
      } else {
        $self->_writeDeeAstar($secondrun, $rot_sub);
      }
    } elsif ($method eq 'dee-scmf' ) {
        if (!defined $rot_sub) { 
            $self->_writeDeeSCMF($secondrun);
        } else {
            $self->_writeDeeSCMF($secondrun, $rot_sub);
        }
    } elsif ($method =~ /mc/ ) {
      $self->_writeDeeMC();
    } else {
      GENERAL::error("Invalid method $method!");
    }
    $self->{handles}->{control}->close();
}

=head2 

 Title   : setMC
 Usage   : $design->setMC(20, 100000, 4000, 150)
 Function: set Monte Carlo simmulated annealing parameters
           (1) number of cycles of heating and quenching procedure. ( quick calculation 20 cycles)
               serious calculation 1000 cycles
           (2) number of monte carlo running at each cycle. typically
               set as 1000000 (one million)
           (3) inital temperature. it should be high. the main purpose 
               is to help the rotamers to overcome some energy barrier. 
               4000K is good value. 
           (4) final temperature 150K               
 Returns : no return value
 Args    :  (1) cycles (default 100)
            (2) steps  (default 10000)
            (3) initial Temperature (default 4000)
            (4) final temperature (default 150)

=cut


sub setMC {

    my $self=shift;
    my $cycles=shift;
    $cycles=100 if (!defined $cycles);
    my $steps=shift;
    $steps=10000 if(!defined $steps);
    my $initialT=shift;
    $initialT=4000 if (!defined $initialT);
    my $finalT=shift;
    $finalT=150 if (!defined $finalT);

    $self->{par}->{'cycles'}= $cycles;
    $self->{par}->{'steps'}=$steps;
    $self->{par}->{'initialT'}=$initialT;
    $self->{par}->{'finalT'}=$finalT;
}


=head2 

 Title   :  _writeDeeAstar
 Usage   :  $self->_writeDeeAstar()
 Function:  Internal method: write dee astar method
 Returns :  none
 Args    :  (1) second run bit
	    (2) rot_subset bit

=cut


sub _writeDeeAstar {

    my $self=shift;
    my $secondrun=shift;
    my $rot_subset=shift;
    $rot_subset='rot_subset' if (!defined $rot_subset);
    my $fhandle = $self->{handles}->{control};
    
    # Need to be aware of file space
    GENERAL::setFileSpace($FS_DEF, 1, 0);

    my $ebin= $self->{par}->{ebin};
    my $ecut= $self->{par}->{ecut};

    $fhandle->print("\# tweaked to avoid 1st order pairs \n");
    $fhandle->print("\#",  "\n\n");
    $fhandle->print("rot_input2     $FS_DEF->{rotinp1f}",  "\n");
    $fhandle->print("rot_input3     $FS_DEF->{rotinp1f}",  "\n");
    $fhandle->print('rot_states     $FS_DEF->{rot_statesf}',  "\n");
    $fhandle->print("energies       $ebin",  "\n");
    $fhandle->print("\#rot_subset     rot_subset \n\n") if ($secondrun ==0);
    $fhandle->print("rot_subset     $rot_subset \n\n") if ($secondrun ==1);
    $fhandle->print("\#------------------sequence below----------------" ,  "\n");
    $fhandle->print(  "\n");
    $fhandle->print('begin-sequence',  "\n\n");
    $fhandle->print("  set-ecut     $ecut \n\n");
    $fhandle->print("  \# skim a few cheap ones off the top, each step takes half a min",  "\n");
    $fhandle->print('  dee_singles_zeroth  repeat',  "\n");
    $fhandle->print('  dee_singles_mb      once',  "\n");
    $fhandle->print('  dee_singles_mb      once',  "\n\n");
    $fhandle->print("  \# heavy singles, light pairs, each step takes half an hr \n");
    $fhandle->print("  begin-loop LOOP-1 1e12          \# LOOP-1 begin \n\n");
    $fhandle->print("    dee_singles_first    repeat \n");
    $fhandle->print("    dee_singles_split_s1 repeat tryall \n\n");
    $fhandle->print('    begin-alternate LOOP-1',  "\n");
    $fhandle->print('      dee_singles_split_s2 once   fixedbylooger',  "\n");
    $fhandle->print('      dee_singles_split_s2 once   fixedbymayo',  "\n");
    $fhandle->print('      dee_singles_split_s2 repeat fixedbylooger',  "\n");
    $fhandle->print('    end-alternate',  "\n\n");
    $fhandle->print('    mem-condense',  "\n");
    $fhandle->print(  "\n");
    $fhandle->print('    dee_pairs_mb    repeat',  "\n");
    $fhandle->print(  "\n");
    $fhandle->print("  end-loop					\# LOOP-1 end",  "\n");
    $fhandle->print(  "\n");
    $fhandle->print("   \# heavier singles, heavy pairs, each step takes a day or so",  "\n");
    $fhandle->print("  begin-loop LOOP-2	1e12          \# LOOP-2 begin",  "\n\n");
    $fhandle->print('    begin-alternate LOOP-2',  "\n");
    $fhandle->print('      dee_pairs_first once',  "\n");
    $fhandle->print('      dee_pairs_mb    repeat',  "\n");
    $fhandle->print('    end-alternate',  "\n\n");
    $fhandle->print('    dee_singles_first    repeat',  "\n");
    $fhandle->print('    dee_singles_split_s1 repeat tryall',  "\n");
    $fhandle->print(  "\n");
    $fhandle->print('    begin-alternate LOOP-2',  "\n");
    $fhandle->print('      dee_singles_split_s2 once   fixedbylooger',  "\n");
    $fhandle->print('      dee_singles_split_s2 once   fixedbymayo',  "\n");
    $fhandle->print('      dee_singles_split_s2 repeat fixedbylooger',  "\n");
    $fhandle->print('    end-alternate',  "\n");
    $fhandle->print(  "\n");
    $fhandle->print('    mem-condense',  "\n");
    $fhandle->print(  "\n");
    $fhandle->print("  end-loop					\# LOOP-2 end",  "\n");
    $fhandle->print(  "\n");
    $fhandle->print('  bnb_astar	1024	leach_pos leach_h*',  "\n");
    $fhandle->print(  "\n");
    $fhandle->print('end-sequence',  "\n");
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub _writeDeeSCMF {

    my $self=shift;

    my $secondrun=shift;
    $secondrun=0 if (!defined $secondrun);
    my $rot_subset=shift;
    $rot_subset='rot_subset' if (!defined $rot_subset);

    my $initial=$self->{par}->{'initial'};
    my $scmf;
    my $fhandle=$self->{handles}->{control};
    my $ebin= $self->{par}->{ebin};

    if ($initial ==0) {
      $scmf="     rand_scmf     init-zero \n";
    } elsif ($initial==1) {
      $scmf="     rand_scmf     init-flat \n";
    } elsif ($initial ==2) {
      $scmf="     rand_scmf     init-zero \n     rand_scmf     init-flat \n";
    } else {
      die " invalid choice for mean field initialization method.\n";
    }
    $fhandle->print("\#",  "\n");
    $fhandle->print("\# standard battery of routines",  "\n");
    $fhandle->print("\#",  "\n");
    $fhandle->print(  "\n");
    $fhandle->print("rot_input2     rot_input2            \n");
    $fhandle->print("rot_input3     rot_input3            \n");
    $fhandle->print('rot_states     rot_states.self',  "\n");
    $fhandle->print("energies       $ebin \n");
    $fhandle->print("\#rot_subset       rot_subset \n\n") if ($secondrun ==0);
    $fhandle->print("rot_subset $rot_subset \n\n") if ($secondrun ==1);
    $fhandle->print(  "\n");
    $fhandle->print("\#-------------------sequence below ------------------\n");
    $fhandle->print(  "\n");
    $fhandle->print('begin-sequence',  "\n");
    $fhandle->print(  "\n");
    $fhandle->print('  dee_singles_zeroth   repeat',  "\n");
    $fhandle->print('  dee_singles_mb      once', "\n");
    $fhandle->print('  dee_singles_mb      once', "\n");
    $fhandle->print('  mem-condense',  "\n");
    $fhandle->print(  "\n");
    $fhandle->print("$scmf \n");
    $fhandle->print(  "\n");
    $fhandle->print('end-sequence',  "\n");
}


=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut


sub _writeDeeMC {

    my $self=shift;
    my $cycles = $self->{par}->{'cycles'};
    my $steps =$self->{par}->{'steps'};
    my $initialT = $self->{par}->{'initialT'};
    my $finalT = $self->{par}->{'finalT'};
    my $fhandle=$self->{handles}->{control};
    my $ebin= $self->{par}->{ebin};

    $fhandle->print("\#",  "\n");
    $fhandle->print("\# standard battery of routines",  "\n");
    $fhandle->print("\#",  "\n");
    $fhandle->print(  "\n");
    $fhandle->print("rot_input2     rot_input2            \n");
    $fhandle->print("rot_input3     rot_input3            \n");
    $fhandle->print('rot_states     rot_states.self',  "\n");
    $fhandle->print("energies       $ebin\n");
    $fhandle->print(  "\n");
    $fhandle->print("\#-------------------sequence below ------------------\n");
    $fhandle->print(  "\n");
    $fhandle->print('begin-sequence',  "\n");
    $fhandle->print(  "\n");
    $fhandle->print('  dee_singles_zeroth   repeat',  "\n");
    $fhandle->print('  dee_singles_mb      once', "\n");
    $fhandle->print('  dee_singles_mb      once', "\n");
    $fhandle->print('  mem-condense',  "\n");
    $fhandle->print(  "\n"); 

    $fhandle->print("  rand_mc      $cycles $steps      $initialT  $finalT        \n");
    $fhandle->print(  "\n");
    $fhandle->print('end-sequence',  "\n");
}



=head2 

 Title   :   runDesign
 Usage   :   $design->runDesign()
 Function:   run the design code
 Returns :   none
 Args    :   (1) the directory containing design input file and energy bin 
                files

=cut


sub runDesign {
    
    my $self=shift;
    my $design=shift;
    $design='design' if (!defined $design);
    my $cur=`pwd`;
    chomp $cur;
    my $dir="$cur/Design" if (!defined $dir);
    my $ebin = $self->{par}->{ebin};

    chdir($dir);
    system("cp ../EnergyTable/$ebin .") if (! -e $ebin); 
    my $inp=$self->{files}->{inp};
    my $prefix=&GENERAL::GetBase($inp);
    my $method=$self->{'method'};
    my $out=$prefix.'.out';
    $self->{files}->{output}=$out;

    system("mv ../$inp .");
    system("$design $inp > $out");

    chdir($cur);
}


=head2 

 Title   :  parseOutput
 Usage   :  $design->parseOutput()
 Function:  parse design output
 Returns :  none
 Args    :  directory containing design output data

=cut

    
sub parseOutput {

    my $self=shift;
    my $method=$self->{'method'};
    my $dir=shift;
    my $cur=`pwd`;
    chomp $cur;
    $dir="$cur/Design" if (!defined $dir);

    my $inp=$self->{files}->{inp};
    my $prefix=&GENERAL::GetBase($inp);
    my $out=$prefix.'.out';
    $self->{files}->{output}=$out if (!defined $self->{files}->{output});

    chdir($dir);

    if ($method eq 'dee-astar' ) {
      $self->_parseDeeAstar();
    } elsif ($method eq 'dee-scmf' ) {
      $self->_parseDeeSCMF();
    } elsif ($method eq 'dee-mc') {
      $self->_parseDeeMC();
    }
    chdir($cur);
}


=head2 

 Title   :  _parseDeeSCMF 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub _parseDeeSCMF {

    my $self=shift;
    my $inital= $self->{par}->{'inital'};

    my $outfh=&GENERAL::GetInFH($self->{files}->{output});
    $self->{handles}->{output}=$outfh;

    my $havesolution=0;
    while (<$outfh>) {
        chomp $_;
	
	if (/stoc_scmf best E/) {
	    $self->{solution}={};
	    $self->{solution}->{probability}=[];
	    $self->{solution}->{results}=[];
	    $havesolution=1;
	}

	# now parse the output in the scmf probability matrix

	if ($havesolution==1 && /^\s+/ ) {
	    my $probrec={};
	    my @lines=split(" ", $_);
	    $probrec->{resid}=shift @lines;
	    $probrec->{probmatrix}=[];
	    while (@lines) {
		my $temprec={};
		my $resrotamer=shift @lines;
		my ($res, $rot)= split(/\./, $resrotamer);

		# we have to do a little bit process to get rid of 
		# leading 0 for each rotamer number
		my @temp=split("", $rot);
		if ($temp[0]==0) {
		    # get rid of leading 0
		    shift @temp;
		}

		# now get a new rotnumber without leading 0
		my $rotnumber=join("", @temp);

		my $prob=shift @lines;
		$temprec->{resname}=$res;
		$temprec->{rotnumber}=$rotnumber;
		$temprec->{prob}=$prob;
		push (@{$probrec->{probmatrix}}, $temprec);
	    }
	    push (@{$self->{solution}->{probability}}, $probrec);
	}

	# now parse the scmf solution
	if ( $havesolution==1 && !/^stoc_scmf/ && !/^\s+/) {
	    my @line2=split(" ", $_);
	    my $srec={};
	    $srec->{energy}=shift @line2;
	    $srec->{config}=[];
	    while (@line2) {
		my ($resname2, $rotnum2)= split(/\./, shift @line2);
		my @temp2=split("", $rotnum2);
		if ($temp2[0]==0) {
		    # get rid of leading 0
		    shift @temp2;
		}
		my $rotnumber2=join("", @temp2);
		my $temp2rec={};
		$temp2rec->{resname}=$resname2;
		$temp2rec->{rotnumber}=$rotnumber2;
		push (@{$srec->{config}}, $temp2rec);
	    }
	    push (@{$self->{solution}->{results}}, $srec);
	}

	if (/stoc_scmf end/ ) {
	    $havesolution=0;
#            print "solution!\n";
	}

    }
}




=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut


sub _parseDeeAstar {

    my $self=shift;
    my $outfh=&GENERAL::GetInFH($self->{files}->{output});
    $self->{handles}->{output}=$outfh;

    my $havesolution=0;
   
    while (<$outfh>) {
	
	if (/branch ordering/) {
	    $self->{solution}={};
	    $self->{solution}->{results}=[];
	    $havesolution=1;
	}

	if ( $havesolution==1 && !/^bnb_astar/ && !/^\s+/) {
	    my @line2=split(" ", $_);
	    my $srec={};
	    $srec->{energy}=shift @line2;
	    $srec->{config}=[];
	    while (@line2) {
                my $temp = shift @line2;
		my ($resname2, $rotnum2)= split(/\./, $temp);
		my @temp2=split("", $rotnum2);
		if ($temp2[0]==0) {
		    # get rid of leading 0
		    shift @temp2;
		}
		my $rotnumber2=join("", @temp2);
		my $temp2rec={};
		$temp2rec->{resname}=$resname2;
		$temp2rec->{rotnumber}=$rotnumber2;
		push (@{$srec->{config}}, $temp2rec);
	    }
	    push (@{$self->{solution}->{results}}, $srec);
	}

	if (/bnb_astar end/ ) {
	    $havesolution=0;
	}

    }
}

=head2 

 Title   : _parseMC
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub _parseDeeMC {
    
    my $self=shift;
    my $outfh=&GENERAL::GetInFH($self->{files}->{output});
    $self->{handles}->{output}=$outfh;

    my $havesolution=0;
   
    while (<$outfh>) {
	
	if (/^rand_mc/ && /cycles/ ) {
	    $self->{solution}={};
	    $self->{solution}->{results}=[];
	    $havesolution=1;
	}

	if ( $havesolution==1 && !/^rand_mc/) {
	    my @line2=split(" ", $_);
	    my $srec={};
	    $srec->{energy}=shift @line2;
	    $srec->{config}=[];
	    while (@line2) {
		my ($resname2, $rotnum2)= split(/\./, shift @line2);
		my @temp2=split("", $rotnum2);
		if ($temp2[0]==0) {
		    # get rid of leading 0
		    shift @temp2;
		}
		my $rotnumber2=join("", @temp2);
		my $temp2rec={};
		$temp2rec->{resname}=$resname2;
		$temp2rec->{rotnumber}=$rotnumber2;
		push (@{$srec->{config}}, $temp2rec);
	    }
	    push (@{$self->{solution}->{results}}, $srec);
	}

	if (/^rand_mc end/ ) {
	    $havesolution=0;
	}

    }
}


=head2 

 Title   :   writeSolution
 Usage   :   $design->writeSolution()
 Function:   write out solutions
 Returns : 
 Args    :   (A) directory of design output. default is top of working dir. 

=cut



sub writeSolution {


    my $self=shift;
    my $method=$self->{'method'};
    my $suffix=shift;
    $suffix="default" if (!defined $suffix);
    $suffix ='-'."$suffix";
    $method .= $suffix; 
    my $round=shift;
    $round=1 if (!defined $round);
    my $dir="Design";

    # read in site information, we need the info
    # on the chain and residue id

    my $site=&GENERAL::GetInFH('rot.inp');
    my @sites=<$site>;
    chomp @sites;
    my $sites=[];
    for (@sites) {
	chomp $_;
	my @lines=split(" ", $_);
	pop @lines;
	pop @lines;
	push (@{$sites}, \@lines);
    }


    # we need create a unique directory for the design output

    &GENERAL::MakeDir("$dir/$method");
    chdir("$dir/$method");

    my $designout=$self->{files}->{output};
    my $designinp=$self->{files}->{inp};
    # move the output input into that directory
    system("mv ../$designout ../$designinp .");

    my $results=$self->{solution}->{results};
    my $i=-1;
    foreach my $r (@{$results}) {
	$i++;
        print $i, "\n";
	my $config=$r->{config};
	my $energy=$r->{energy};
	
	my $out='solution-'."$method".'-'.$i;
        my $pmat;
        my $pmatfh;
	my $rotsubset;
	my $rotsubsetfh;

        # if we use dee-scmf, create a prob matrix files

        if ($self->{'method'} eq 'dee-scmf') {
            $pmat='probmat-'.$i;
            $pmatfh=&GENERAL::GetOutFH($pmat);
	    $rotsubset='rot_subset';;
	    $rotsubsetfh=&GENERAL::GetOutFH($rotsubset);
       }

        my $outfh;
        if ($self->{'method'} ne 'dee-mc') {
	   $outfh=&GENERAL::GetOutFH($out);
	   $outfh->print("# The energy is $energy \n");
        }

        my $outmc;
        my $outmcfh;
        if ($self->{'method'} eq 'dee-mc' &&  $i == $#{$results}) {
           $outmc='solution-'."$method".'-'.$round;
           $outmcfh=&GENERAL::GetOutFH($outmc);
           $outmcfh->print("# The energy is $energy \n");
        } elsif ($self->{'method'} eq 'dee-mc' && $i < $#{$results})  {
           next;
        }
	my $j=-1;
	foreach my $conf ( @{$config}) {
	    $j++;
	    $outfh->printf("%1s%5d%5s%3d\n", $sites->[$j]->[0], $sites->[$j]->[1], $conf->{resname}, $conf->{rotnumber})
                    if ($self->{'method'} eq 'dee-astar');
            if ($self->{'method'} eq 'dee-mc') {
                   $outmcfh->printf("%1s%5d%5s%3d\n", $sites->[$j]->[0], $sites->[$j]->[1], $conf->{resname}, $conf->{rotnumber});
            }
            $outfh->printf("%1s%5d%5s%3d  %4.2f\n", $sites->[$j]->[0], $sites->[$j]->[1], $conf->{resname}, $conf->{rotnumber}, 
                   $self->{solution}->{probability}->[$j]->{probmatrix}->[0]->{prob}) if ($self->{'method'} eq 'dee-scmf');
            if ($self->{'method'} eq 'dee-scmf') {
                 $pmatfh->printf("%1s%5d", $sites->[$j]->[0], $sites->[$j]->[1]);
		 $rotsubsetfh->printf("%1s%5d", $sites->[$j]->[0], $sites->[$j]->[1]);
		 my $probmatrix=$self->{solution}->{probability}->[$j]->{probmatrix};
                 for (my $k=0; $k<=$#{$probmatrix}; $k++) {
                      $pmatfh->printf("%5s%3d%6.2f", $probmatrix->[$k]->{resname}, $probmatrix->[$k]->{rotnumber}, $probmatrix->[$k]->{prob});
                      if ( $k == $#{$probmatrix} || $k == 9) {
			  $pmatfh->print("\n");
			  last;
                      }
                  }  

		for (my $l=0; $l <=$#{$probmatrix}; $l++) {
		      $rotsubsetfh->printf("%5s%3d", $probmatrix->[$l]->{resname}, $probmatrix->[$l]->{rotnumber});
                      if ( $l == $#{$probmatrix} || $l == 9) {
                          $rotsubsetfh->printf("\n");
                          last;
                      }
		}	
	     }
	}
    }

    chdir("../../");
}


=head2 

 Title   : convert2PDB
 Usage   : $design->convert2PDB();
 Function: convert the design solution into PDB file
           convert to pdb without any hydrogen
           (A) based on the design solution
           (B) according to backbone dependent rotamer library 
 Returns : 
 Args    : (A) the reference to rotamer
           (B) crd file name for crystal structure
           (C) suffix for the solution directory

=cut


sub convert2PDB {

    my $self=shift;
    my $rotamer=shift;

    # we need crd file in order to run charmm 
    my $crd=shift;

    # we need suffix inorder to determine which directory to enter
    my $suffix=shift;
    $suffix='default' if (!defined $suffix);
    $suffix ='-'."$suffix";      
    my $dir=shift;
    my $cur=`pwd`;
    chomp $cur;
    $dir="$cur/Design" if (!defined $dir);
    
    # first contruct proper directory. 
    # we don't want to mess up whole diretcory

    # first get all the 
    my $sitelist=$rotamer->{dslist};
    my $chidef=$rotamer->{chidef};

    my $method=$self->{'method'};

    # add suffix to method. This is where all the solutions sit
    # note: suffix in this subroutine must be same as 
    # the one in above! 

    $method .=$suffix;
    &GENERAL::MakeDir("$dir/$method");
    system "cp $crd $dir/$method";
    chdir("$dir/$method");

    # build charmm input files
    # (1) loop over list of solutions 

    my $results=$self->{solution}->{results};
    my $i=-1;
    foreach my $r (@{$results}) {
	$i++;
	my $config=$r->{config};
	my $energy=$r->{energy};

        if ($self->{'method'} eq 'dee-mc' &&  $i == $#{$results} || $self->{'method'} ne 'dee-mc') {

	    # first step: delete all the original sites in crystal structure
	    my $charmm=CHARMM::new('charmm', 'temp.inp');
	    $charmm->loadParm('param', 19);
	    # here we read crd file, set hbond flag as 0 and set terminal 
	    $charmm->setupFromCRD($crd,0, "ter");
	    $charmm->verbose("ic build\n\n");
	    $charmm->defineTemplate($rotamer);
	    $charmm->deleteDesignSite();
	    
	    # for each designed site, get chain id and residue number
	    
	    $charmm->verbose("bomlev -2");
	    
	    my $index0=-1;	
	    foreach my $res ( @{$config}) {
		$index0++;
		my $chain=$sitelist->[$index0]->{chain};
		my $iresnum=$sitelist->[$index0]->{iresnum}; 
		my $resname=$res->{resname};
		$charmm->patchRes($resname, $chain, $iresnum);
	    }
	    
	    $charmm->defineTemplate($rotamer);	
	    $charmm->verbose("ic param\n\ncoor init select rots end\nic param\nic edit\n\n"); 
 
	    
	    # loop over list of residues
	    my $index=-1;
	    foreach my $res (@{$config}) {
	    $index++;
	    my $resname=$res->{resname};
	    
	    # obtain right chain id and indexed residue number
	    my $chain=$sitelist->[$index]->{chain};
	    my $iresnum=$sitelist->[$index]->{iresnum};
	    
	    # retrieve right back-dependent rotamer library from rotamer object
	    my $libref=$sitelist->[$index]->{rotlibref};

	    # retrieve the resdiues in rotamer library from reference to rotamer library
	    my $resref=$rotamer->getResInRotlib($libref,$resname);
	    
	    # retrieve the chi values
	    my $chi = $resref->{rotamerlist}->[$res->{rotnumber}]->{chi};
	    # we use manyflag string to charmm object don't build this rotamer
	    $charmm->rebuildRotamer($chidef, $chi, $resname, $chain, $iresnum, "manyflag");
	}
	    
	    # now tell charmm to build angles
	    $charmm->verbose("end\nic build\n");	    
	    
	    # tell charmm to write pdb file
	    $charmm->writeCRD('temp.crd');
	    
	    # run charmm. we use backtick to avoid annoying fortran stop stderr
	    `charmm < temp.inp 2>&1 1>& temp.out`;
	    
	    # now clean up a little bit
 #	    &GENERAL::Remove('temp.inp');
#	    &GENERAL::Remove('temp.out');
	    
	    # postprocess pdb files. remove the hydrogens
	    my $pdb=PDB::new('temp.crd');
	    my $newpdb='solution-'."$method".'_'."$i".'.pdb';
	    $pdb->writePDB($newpdb, "GENERIC_NOH_RENUMBERED");
            &GENERAL::Remove('temp.crd');
	    
        } elsif ($self->{'method'} eq 'dee-mc' && $i < $#{$results})  {
	    next;
        }
    }
    chdir("../../");
}


=head2 

 Title   :  getStatistics
 Usage   : 
 Function:  given crystal chi values, spit out the statistcs for the rotamer 
            prediction.  
 Returns : 
 Args    : 

=cut


sub getStatistics {

    my $self=shift;

    # get rotamer object to retrieve the chi values
    my $rotamer=shift;
    my $sitelist=$rotamer->{dslist};

    # get crsytal chi values
    my $crystalchis=shift;


    # we need suffix inorder to determine which directory to enter
    my $suffix=shift;
    $suffix='default' if (!defined $suffix);
    $suffix ='-'."$suffix"; 

    my $dir=shift;
    my $cur=`pwd`;
    chomp $cur;
    $dir="$cur/Design" if (!defined $dir);

    my $method=$self->{'method'};

    # add suffix to method. This is where all the solutions sit
    # note: suffix in this subroutine must be same as 
    # the one in above! 

    $method .=$suffix;
    chdir("$dir/$method");

    my $results=$self->{solution}->{results};
    my $i=-1;
    foreach my $r (@{$results}) {
	$i++;
	my $config=$r->{config};

	if ($self->{'method'} eq 'dee-mc' &&  $i == $#{$results} || $self->{'method'} ne 'dee-mc') {
	    
	    # loop over list of residues

	    # here are some variables to hold the statistics

	    my $index=-1;
	    my $chi1correct=0;
	    my $chi12correct=0;
	    my $allcorrect=0;
	    my $total=$#{$config}+1.0;

	    # here is the reference to list of list to hold the comparison result
	    my $cmatrix=[];

	    # some file handles to hold output
	    my $outf='solution-'.$i.'.summary';
	    my $outfh=&GENERAL::GetOutFH($outf);

	    foreach my $res (@{$config}) {


                #now compare two sets of chi vaules.
                #varible to hold the comparison result and flipping process
                my $cflag;
                my $flip;

		$index++;
		my $resname=$res->{resname};
		
		# obtain right chain id and indexed residue number
		my $chain=$sitelist->[$index]->{chain};
		my $iresnum=$sitelist->[$index]->{iresnum};
                # print $resname, $chain, $iresnum, "\n";

                if ($resname eq "ALA" || $resname eq "GLY") {
                        $cflag=1;
                        push (@{$cmatrix->[$index]}, $cflag);
                }
		
		# retrieve right back-dependent rotamer library from rotamer object
		my $libref=$sitelist->[$index]->{rotlibref};

		# retrieve the resdiues in rotamer library from reference to rotamer library
		my $resref=$rotamer->getResInRotlib($libref,$resname);
	    
		# retrieve the chi values
		my $chi = $resref->{rotamerlist}->[$res->{rotnumber}]->{chi};

		# get crystal chi value for that residue
		my $crystalchi=$crystalchis->[$index];


		# loop over chi values
		for (my $j=0; $j<=$#{$chi}; $j++) {
                
		    my $da = $chi->[$j] - $crystalchi->[$j];
		    $da+=360.0 if ($da<-180.0);
		    $da-=360.0 if ($da>180.0);

		    if ( $da>=-40.0 && $da<=40.0 ) {
			$cflag=1;
			push (@{$cmatrix->[$index]}, $cflag);

			# now i have to correct chi2 stuff...
			# note: no need for TRP. libray already included one!
		    } elsif ( ($resname =~ /AS/ || $resname eq "HIS" || $resname eq "PHE"
			       || $resname eq "TYR") && ($j ==1) && ($da<-40.0 || $da>40.0) ) {

			# add 180 degree to calculated chi2. do comparision with dunbrack library
			$da +=180.0;
			$da +=360.0 if ($da<-180.0);
			$da -=360.0 if ($da>180.0);
			$cflag=1 if ($da>-40.0 && $da<40.0);
			$flip=1 if ($da>-40.0 && $da<40.0);
			push (@{$cmatrix->[$index]}, $cflag) if ($da>-40.0 && $da<40.0)  ;

			# now flipping chi3 of GLU and GLN
		    } elsif ( ($resname eq "GLN" || $resname eq "GLU" ) && ($j ==2) && ($da<-40.0 || $da>40.0) ) {
			$da +=180.0;
			$da +=360.0 if ($da<-180.0);
			$da -=360.0 if ($da>180.0);
			$cflag=1 if ($da>-40.0 && $da<40.0);
			$flip=1 if ($da>-40.0 && $da<40.0);
			push (@{$cmatrix->[$index]}, $cflag) if ($da>-40.0 && $da<40.0)  ;

		    } else {

			# if we try every thing, still far away from crystal chi value
			
			$cflag=0;
			push (@{$cmatrix->[$index]}, $cflag);
		    }
		}
	    }
	
	    # after all these steps, we have a cmatrix contains results for all the residues.
	    # now let's get the statistics

	    my $k=-1;
	    my $chi1flag=0;
	    my $chi12flag=0;
	    my $allflag=0;
	    
	    foreach my $residue (@{$config}) {
		$k++;
		$chi1correct++ if ($cmatrix->[$k]->[0] == 1);
		$chi1flag=1  if ($cmatrix->[$k]->[0] == 1);
                $chi1flag=0 if ($cmatrix->[$k]->[0] != 1);


		if (defined $cmatrix->[$k]->[1]) {
		         if ($cmatrix->[$k]->[0]==1 && $cmatrix->[$k]->[1]==1)  {
                                       $chi12correct++;
		                       $chi12flag=1;
                         } else {
				       $chi12flag=0;
                         }
                }

                if (!defined $cmatrix->[$k]->[1]) {
                         if ($cmatrix->[$k]->[0]==1) {
                                        $chi12correct++;
					$chi12flag=1;
                         } else {
                                        $chi12flag=0;
                         }
                }


		# this step will print the summary of comparision. 
		# 1 is for correct. 0 is incorrect

		my $sum=-1;
		foreach my $item (@{$cmatrix->[$k]}) {
		    $sum +=$item;
		}

		if ($sum == $#{$cmatrix->[$k]}) {
		    $allcorrect++;
		    $allflag=1;
		} else {
	            $allflag=0;
                }

		my $resname=$residue->{resname};
		
		# obtain right chain id and indexed residue number
		my $chain=$sitelist->[$k]->{chain};
		my $iresnum=$sitelist->[$k]->{iresnum};

		$outfh->printf("%1s  %4d   %3s  %1d  %1d  %1d\n", $chain, $iresnum, $resname, $chi1flag, $chi12flag, $allflag);
	    }

	    # print out the summary of final output

	    $outfh->printf("Total Correctness for Chi1 is %5.2f:  %3d out of %3d.\n", $chi1correct/$total, $chi1correct, $total);
	    $outfh->printf("Total Correctness for Chi1 and Chi2 is %5.2f:  %3d out of %3d.\n", $chi12correct/$total, $chi12correct, $total);
	    $outfh->printf("Total Correctness for All Chi is %5.2f:  %3d out of %3d.\n", $allcorrect/$total, $allcorrect, $total);
		
	}
    }

    chdir("../../");
	    
}


1;













