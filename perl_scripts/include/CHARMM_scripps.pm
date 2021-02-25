# CHARMM package
# interact with CHARMM program
#
# http://mmtsb.scripps.edu/doc/CHARMM.pm.html
# 2001, Michael Feig, Brooks group, TSRI
# 2001, John Karanicolas, Brooks group, TSRI
# 
# Downloaded on 07/19/05 by Gevorg Grigoryan.

package CHARMM;

require 5.004;

use strict;

use IO::Handle;
use IO::File;
use IPC::Open2;
use Sys::Hostname;
use Fcntl;

use GenUtil;
use SICHO;
use Molecule;

## data: handle -> {outlog, cmdlog, enerlog, pertlog, fromcharmm, tocharmm }
## file handles used for logging and communication with CHARMM

## data: par -> { ... } 
## CHARMM parameters used in various commands

## data: mdener[] -> { step, time, temp, total, kine, 
## data:               pot, vdwaals, elec, gb, asp, constr }
## energy data extracted from the last molecular dynamics run

## data: molecule
## Molecule object for current molecule structure

use vars qw ( $exec $datadir );

BEGIN {
  if ($ENV{'CHARMMEXEC'} ne "") {
    $exec=$ENV{'CHARMMEXEC'};
  } else {
    $exec=&GenUtil::findExecutable("charmm");
  }
  
  if ($ENV{'CHARMMDATA'} ne "") {
    $datadir=$ENV{'CHARMMDATA'};
  } elsif ($datadir eq "" && $exec=~/exec/) {
    ($datadir=$exec)=~s/exec.*//;
    $datadir.="toppar";
  }
}

## constructor: new([logfile[,commandlogfile]])
## creates a new CHARMM object. A log file for 
## CHARMM output and a command log file for commands
## sent to CHARMM may be given as an option.
## The constructor opens a connection to CHARMM and
## sets all parameters to default values.

sub new {
  my $self={};
  my $logoutfile=shift;
  my $logcmdfile=shift;
  my $prunnodes=shift;
  my $prunbase=shift;

  srand(time ^ ($$ + ($$ << 15)));

  bless $self;

  $self->{handle}={};

  if (defined $logoutfile) {
    $self->{handle}->{outlog}=&GenUtil::getOutputFile($logoutfile);
    $self->{handle}->{outlog}->autoflush(1);
  }

  if (defined $logcmdfile) {
    $self->{handle}->{cmdlog}=&GenUtil::getOutputFile($logcmdfile);
    $self->{handle}->{cmdlog}->autoflush(1);
  }

  die "cannot execute binary. $exec" if ((!-x $exec) && ($exec !~ "mpi"));

  if (defined $prunnodes && defined $prunbase) {
    $exec="prun -n $prunnodes -B $prunbase $exec";
  }

#  printf STDERR "executing >$exec<\n";
  my $o2read=new IO::Handle;
  my $o2write=new IO::Handle;

  $self->{_charmmpid}=open2($o2read,$o2write,$exec) || 
    die "open2 for $exec failed";

  printf $o2write "* title\n*\n\n";
  printf $o2write "UNBUFIO\n";

  $self->{handle}->{fromcharmm}=$o2read;
  $self->{handle}->{tocharmm}=$o2write;

  $self->{_lastOutput}="";
  $self->_getCHARMMOutput();

  my %parhash = ( 
   bomlev    => undef,          # set bomb level
   prnlev    => undef,          # set print level

   param     => 22,             # parameter set: 19 (CHARMM19), 22 (CHARMM22), opls, a94 (amber 94), eef1

   xtop      => undef,          # extra topology files
   xpar      => undef,          # extra parameter files

   blocked   => undef,          # blocked terminii
   nter      => undef,          # N-terminus, requires blocked
   cter      => undef,		# C-terminus, requires blocked
   nodeoxy   => undef,          # no deoxy patch for nucleic acids

   hsd       => "",             # list of HSD residues
   hse       => "",             # list of HSE residues
   hsp       => "",             # list of HSP residues

   resmod    => "",             # list of residue name modifications		 

   patch     => "",             # list of patches

   buildall  => 1,              # rebuild all missing atoms

   dielec    => 'CDIE',         # dielectric response CDIE (const.) or RDIE (distant dep.)
   epsilon   => 1.0,            # dielectric constant for non-bonded interactions
   cutnb     => undef,          # cutoff for non-bonded list generation
   cutoff    => undef,          # cutoff for non-bonded interactions
   cuton     => undef,          # onset of switching function for non-bonded interaction
   cut       => 1,              # use cutoff
   trunc     => undef,          # shift/switch

   user      => 1,              # USER energy term
   cmap      => 0,              # CMAP energy term

   echeck    => 20.0,           # energy tolerance check for dynamics

   ace       => 0,              # ACE implicit solvation
   aceieps   => 1.0,            # ACE internal dielectric
   aceseps   => 80.0,           # ACE external dielectric
   acealpha  => 1.3,            # ACE alpha
   acesigma  => undef,          # ACE sigma
   acevscale => 1.0,            # ACE volume scale

   sasagamma => 0.00542,        # gamma for SASA
   sasadelta => 0.920,          # delta for SASA

   gb        => undef,          # Generalized Born implicit solvent: bgb, gbmf, gbmf2, gbmva, gbmvg
   gbeps     => 80.0,           # epsilon in GB prefactor

   gbcrad    => 0.91,           # BGB radius scaling of main chain carbon w/ param19
   gbhrad    => 0.8,            # BGB radius scaling of hydrogen atoms w/ param22
   gblambda  => undef,          # BGB lambda factor

   gbmvad     => undef,         # GBMVA shift
   gbmvade    => undef,         # GBMVA eshift
   gbmvas     => undef,         # GBMVA slope
   gbmval1    => 0.5,           # GBMVA lambda 1
   gbmvap1    => 0.45,          # GBMVA P1
   gbmvap2    => 1.25,          # GBMVA P2
   gbmvap3    => 0.70,          # GBMVA P3
   gbmvap4    => undef,         # GBMVA P4
   gbmvap5    => undef,         # GBMVA P5
   gbmvap6    => 8.0,           # GBMVA P6, Still's factor
   gbmvap7    => undef,         # GBMVA P7, VSA adjustment
   gbmvaig    => undef,         # GBMVA integration grid size
   gbmvafrq   => undef,         # GBMVA alpha update frequency
   gbmvdrfrq  => undef,         # GBMVA surface integration point update frequency
   gbmvaimp   => undef,         # GBMVA impulse derivatives
   gbmvaemp   => undef,         # GBMVA exponentially damped impulse
   gbmvarith  => undef,         # GBMVA arithmetic mean: alpha_ij = alpha_i+alpha_j
   gbmvabuf   => 0.1,           # GBMVA buffer for dynamics
   gbmvabeta  => -20,           # GBMVA beta parameter
   gbmvaonx   => 1.9,           # GBMVA OnX parameter
   gbmvaoffx  => 2.1,           # GBMVA OffX parameter
   gbmvacubic => undef,         # GBMVA cubic VSA
   gbmvacorr  => 1,             # GBMVA r7 correction (1) or r5 correction (0)
   gbmvaextra => undef,         # GBMVA extra commands

   gbmva1     => undef,         # GBMV A1 parameter
   gbmva2     => undef,         # GBMV A2 parameter
   gbmva3     => undef,         # GBMV A3 parameter
   gbmvdf     => undef,         # GBMV adjustment factor for shift
   gbmvdoff   => undef,         # GBMV epsilon offset for shift adjustment
   gbmvmstl   => undef,         # GBMV modified Still
   gbmvms0    => 1.0,           # GBMV MS0 parameter
   gbmvms1    => 0.0,           # GBMV MS1 parameter
   gbmvms2    => 0.0,           # GBMV MS2 parameter
   gbmvms3    => 0.0,           # GBMV MS3 parameter

   gbmveps    => undef,         # GBMV solvent dielectric
   gbmvepsu   => undef,         # GBMV solute ref dielectric

   gbmvdu     => undef,         # GBMV adj. delta for solute dielectric

   gbmvsa     => 0.005,         # GBMV hydrophobic SASA term factor
   gbmvsb     => 0.900,         # GBMV hydrophobic SASA term offset 
   gbmvsagb   => undef,         # GBMV hydrophobic SASA term on a selection basis
                                # example: "all:0.00542;type n:0.01;type o:0.02"
   gbmvason   => 1.2,           # GBMV SASA internal parameter, old: 0.5
   gbmvasoff  => 1.5,           # GBMV SASA internal parameter, old: 1.75

   gbswsw     => 0.3,           # GBSW switching length
   gbswrmax   => 20.0,          # GBSW maximum integration radius
   gbswnang   => 38,            # GBSW number of angular integration points
   gbswnrad   => 0,             # GBSW number of radial integration points
   gbswsgamma => 0.00,          # GBSW surface tension coefficient
   gbswdgp    => 1.5,           # GBSW grid spacing for lookup table
   gbswrbuffer => 0.0,          # GBSW buffer region for lookup table
   gbswtmemb  => 0.0,           # GBSW thickness of membrane
   gbswmsw    => 0.3,           # GBSW membrane switching length
   gbswms     => 0,             # GBSW molecular surface approx.

   scalerad  => 0,              # scaled radii for GB/PB

   dcel      => 0.5,            # PB grid spacing
   epsw      => 80.0,           # epsilon for solvent
   epsp      => 1.0,            # epsilon for solute
   epspr     => undef,          # epsilon for solute reference if set
   epsr      => 1.0,            # epsilon for reference environment
   smooth    => 0,              # smooth boundary
   pbdelta   => 4.5,            # distance to boundary
   proberad  => 1.4,            # water probe radius

   sradfac   => 1.0,            # overall scaling factor for GB/PB radii
   sbackn    => 1.0, # 0.90,    # peptide backbone N scaling factor
   sbacko    => 1.0, # 0.94,    # peptide backbone O scaling factor
   sbackca   => 1.0,            # peptide backbone CA scaling factor

   mfnh     => undef,          # radius for H*
   mfn3     => undef,          # radius for ch3
   mfn2     => undef,          # radius for cb
   mfnb10   => undef,          # radius for backbone CA (GLY)
   mfnb9    => undef,          # radius for backbone CA
   mfnb6    => undef,          # radius for backbone C
   mfnb7    => undef,          # radius for backbone O
   mfnb8    => undef,          # radius for backbone N
   mfn5     => undef,          # radius for ch
   mfn13    => undef,          # radius for roh   
   mfn12    => undef,	       # radius for oh
   mfn4     => undef,          # radius for cet   
   mfn19    => undef,          # radius for crh
   mfn18    => undef,          # radius for cpro
   mfn7     => undef,          # radius for ctrp
   mfn16    => undef,          # radius for ntrp
   mfn8     => undef,          # radius for css      
   mfn17    => undef,          # radius for sulfur
   mfn25    => undef,          # radius for s2
   mfn6     => undef,          # radius for cam
   mfn15    => undef,          # radius for nam
   mfn11    => undef,          # radius for oam   
   mfn10    => undef,          # radius for ocar
   mfn9     => undef,          # radius for cze
   mfn14    => undef,          # radius for nal
   mfn21    => undef,          # radius for nhis

   nina9     => undef,          # radius for ARG CZ, LYS CE
   nina10    => undef,          # radius for GLU/ASP OE*/OD*
   nina11    => undef,          # radius for ASN/GLN OE*/OD*
   nina14    => undef,          # radius for ARG NH*/NE, LYS NZ 
   nina7     => undef,          # radius for TRP CE*/CD*/CZ*/CH2
   nina16    => undef,          # radius for TRP NE1
   nina12    => undef,		# radius for OG*
   nina3     => undef,          # radius for CG*
   ninah     => undef,          # radius for H*
   ninab1    => undef,          # radius for backbone CAY/CAT
   ninab2    => undef,          # radius for backbone CY
   ninab3    => undef,          # radius for backbone OY
   ninab4    => undef,          # radius for backbone NT
   ninab5    => undef,          # radius for backbone OT*
   ninab6    => undef,          # radius for backbone C
   ninab7    => undef,          # radius for backbone O
   ninab8    => undef,          # radius for backbone N
   ninab9    => undef,          # radius for backbone CA
   ninab10   => undef,          # radius for backbone CA (GLY)
   nina2     => undef,          # radius for CB
   nina4     => undef,          # radius for GLU CG
   nina5     => undef,          # radius for CD*
   nina6     => undef,          # radius for GLN CD/ASN CG/GLU CD/ASP CG
   nina18    => undef,          # radius for PRO GB/CG/CD
   nina19    => undef,          # radius for TYR/PHE CE*/CD*/CZ
   nina8     => undef,          # radius for MET CE
   nina20    => undef,          # radius for HSD CE1/HSD CD2
   nina13    => undef,          # radius for TYR OH
   nina21    => undef,          # radius for HSD NE2/ND1
   nina22    => undef,          # radius for HSP NE2/ND1
   nina15    => undef,          # radius for GLN NE2/ASN ND2
   nina17    => undef,          # radius for S*

   asp       => 0,              # flag to switch on ASP energy terms
   aspfile   => undef,          # external input file for ASP energy term
   aspval    => 5.42,           # ASP solvation parameter
   aspvalc   => undef,          # ASP solvation parameter for carbon atoms
   aspvaln   => undef,          # ASP solvation parameter for nitrogen atoms
   aspvalo   => undef,          # ASP solvation parameter for oxygen atoms
   aspvalh   => undef,          # ASP solvation parameter for hydrogen atoms
   aspvals   => undef,          # ASP solvation parameter for sulphur atoms
   aspradh   => undef,          # ASP hydrogen radius 
   aspref    => 0.0,            # ASP reference area

   og        => 0,              # flag to switch on Olgun's energy score

   eef1      => 0,              # flag to switch on EEF1 energy terms  
   eef1file  => "solvpar.inp",  # EEF1 parameter file

   sasa      => 0,              # flag to switch on SASA energy terms
   sasaeps   => 2.0,            # SASA dielectric constant
   sasasig1  => -0.060,         # SASA sig 1 value
   sasasig2  => 0.012,          # SASA sig 2 value
   sasasig3  => 0.00,           # SASA sig 3 value

   sdsteps   => 50,             # number of minimization steps for steepest descent min.
   sdstepsz  => 0.005,          # initial step size for steepest descent min.
   updnbsd   => -1,             # update freq. for non-bonded list in steepest descent min.

   minmode   => "abnr",         # minimization protocol
   minsteps  => 500,            # number of minimization steps for min.
   minstepsz => 0.005,          # initial step size for min.
   minetol   => 1.0E-5,         # energy tolerance for min.
   minupdnb  => -1,             # update freq. for non-bonded list in min.
   minoutfrq => 10,             # output freq. for minimization runs

   shake     => 0,              # use shake
   shaketol  => 1.0E-8,         # shake tolerance
   shakemode => 'hyd',          # 'hyd'rogen or 'all' atoms are restrained
   shakefast => 0,              # use fast vector/parallel version of shake

   periodic  => 1,              # use periodic boundaries for explicit solvent
   explicit  => 1,              # do not use implicit solvent for explicit solvent
   ewald     => 1,              # use Ewald for explicit solvent
   pmekappa  => 0.32,           # PME distribution width
   npme      => 64,             # number of grid points for PME grid

   boxshape  => "cubic",        # box shape for periodic boundaries
   boxsize   => undef,          # box size

   dynens    => 'NVT',          # dynamics ensemble 'NVT' or 'NPT'
   dynber    => 0,              # use Berendsen thermostat
   dynbertc  => 5.0,            # coupling constant for Berendsen
   dynnose   => undef,          # use Nose-Hoover thermostat
   dynnoseq  => 50.0,           # Nose-Hoover qref
   dynnosen  => 5,              # Nose-Hoover ncyc
   dyntstep  => 0.002,          # dynamics time step in picoseconds
   dynsteps  => 100,            # dynamics steps
   dyntemp   => 298,            # dynamics temperature
   dynitemp  => undef,          # initial temperature for heating/equilibration 
   dynseed   => undef,          # dynamics random seed
   dyndeltat => 10.0,           # temperature change during heating/cooling
   dynhtfrq  => -1,             # temperature change frequency
   dyntrfrq  => undef,          # trans/rot removal frequency
   dyneqfrq  => 200,            # equilibration frequency
   dyntwin   => 5.0,            # temperature window during equilibrition
   dynscale  => undef,          # temperature scaling factor for restarts
   dynoutfrq => undef,          # frequency of energy/trajectory output
   dynupdnb  => -1,             # update of non-bonded list
   dynupdimg => -1,             # update of images
   dynpress  => 1.00,           # the pressure (in atm) for NPT simulations
   nblisttype=> undef,         # non-bonded list generation
   dyncons   => undef,          # KDYN dynamics
   dynkmas   => 0.1,            # KDYN dynamics mass
   dynkfkc   => undef,          # KDYN FKC
   dynkdkc   => undef,          # KDYN DKC

   sgmd      => 0,              # self-guided MD
   sgmddrag  => 0.2,            # SGMD drag force
   sgmdavg   => 1.0,            # SGMD average?
   sgmdvfb   => -1.0,           # SGMD ?
   sgmdfrq   => 1,              # SGMD update frequency

   noerest   => undef,          # file name with NOE restraint commands
   xnoerest  => undef,          # X-plor file name with NOE restraint commands
   noerexp   => 3,              # exponent for distance averaging
   noekmin   => 1.0,            # force constant for small distance harmonic potential
   noekmax   => 1.0,            # force constant for large distance harmonic potential
   noefmax   => 9999.0,         # maximum force
   noescale  => 1.0,            # overall scaling factor

   nlambda   => 10,             # number of windows for PERT sampling
   pequi     => 1000,           # number of equilibration steps
   pprod     => 4000,           # number of production steps
   qcor      => undef,          # Ewald charge correction
   whamtol   => 0.05,           # WHAM tolerance
   whamsteps => 200,            # WHAM iteration steps
   
   lang      => 0,              # flag for running Langevin dynamics
   langrbuf  => 0.0,            # radius of Langevin inner buffer sphere
   langfbeta => 1.0,            # friction coefficient for Langevin dynamics
   langsel   => undef,          # CHARMM atom selection to apply Langevin friction
   langupd   => 50              # update of Langevin region
	       );

  $self->{par}=\%parhash;

  $self->{_lastCommand}="";
  $self->{_havehmcm}=0;
  $self->{_havermsd}=0;

  $self->{_setup}=0;
  $self->{_setupPeriodic}=0;

  $self->{_keepRestraint}=0;

  $self->{_biasstatus}={};

  $self->{explicitWater}=0;

  $self->{_haveGB}=0;
  $self->{_haveASP}=0;
  $self->{_haveOG}=0;
  $self->{_haveACE}=0;
  $self->{_haveSHAKE}=0;
  
  $self->{_loadedref}=0;

  return $self;
}

DESTROY {
  my $self=shift;
  
  $self->finish() 
    if (!defined $self->{_finished} || !$self->{_finished});
}


## method: setEnergyLogFile(file)
## sets an energy log file and switches on energy logging

sub setEnergyLogFile {
  my $self  = shift;

  undef $self->{handle}->{enerlog}
    if (defined $self->{handle}->{enerlog});
  
  $self->{handle}->{enerlog}=&GenUtil::getOutputFile(shift);
  $self->{handle}->{enerlog}->autoflush(1);
} 

## method: setPertLogFile(file)
## sets an pert log file

sub setPertLogFile {
  my $self  = shift;

  undef $self->{handle}->{pertlog}
    if (defined $self->{handle}->{pertlog});
  
  $self->{handle}->{pertlog}=&GenUtil::getOutputFile(shift);
  $self->{handle}->{pertlog}->autoflush(1);
} 


## method: logEnergy(tag)
## writes the energy output from the last 
## energy evaluation or minimization run to the
## energy log file if one has been set before
## with <mark>setEnergyLogFile</mark>. A tag is 
## required to provide
## identification of the output in the log file.

sub logEnergy {
  my $self=shift;
  my $tag=shift;

  if (defined $self->{handle}->{enerlog} ) {
    my $earr=$self->_processEneOutput();
    foreach my $e (@{$earr}) {
      $self->{handle}->{enerlog}->printf(
	"%-15s %4d %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",
        "$tag:",$e->{step},$e->{total},$e->{vdwaals},
	$e->{elec},$e->{gb},$e->{asp},$e->{constr});
    }
  }
}

## method: logMDEnergy(tag)
## writes energy output from the last dynamics run
## to the energy log file if one has been set before.
## A tag is required as in <mark>logEnergy</mark>.

sub logMDEnergy {
  my $self=shift;
  my $tag=shift;

  if (defined $self->{handle}->{enerlog} &&
      defined $self->{mdener}) {
    foreach my $e (@{$self->{mdener}}) {
      $self->{handle}->{enerlog}->printf(
 	"%-15s %4d %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.6f %10.6f\n",
        "$tag:",$e->{step},$e->{time},$e->{temp},$e->{total},$e->{kine},$e->{pot},
        $e->{vdwaals},$e->{elec},$e->{gb},$e->{asp},$e->{constr},$e->{volume},$e->{boxsize});      
    }
  }
}

## method: logPertOutput()
## writes pert output to the pert log file if one has been set before.

sub logPertOutput {
  my $self=shift;

  if (defined $self->{handle}->{pertlog} &&
      defined $self->{pertdata}) {
    foreach my $e (@{$self->{pertdata}}) {
      $self->{handle}->{pertlog}->printf(
 	"%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n",
        $e->{llast},$e->{lstart},$e->{lstop},$e->{tptot},$e->{titot},$e->{difave});
    }
  }
}


## method: finish()
## stops CHARMM and closes the connection 
## and all log files. After this function
## has been called CHARMM commands can no
## longer be used. It is called automatically
## by the package destructor if it has not
## called explicitly.

sub finish {
  my $self=shift;

  $self->_sendCommand("stop",1);
  sleep 1;
  $self->_getCHARMMOutput();

  $self->_closeAll();

  $self->{_finished}=1;
}

## method: closeEnergy()
## closes the energy log file

sub closeEnergy {
  my $self=shift;
  
  undef $self->{handle}->{enerlog} 
    if (defined $self->{handle}->{enerlog});
}

## method: loadParameters([parameters])
## CHARMM load topology and parameter files into CHARMM.
## This should be called as one of the first CHARMM
## commands before any structures are loaded.
## Depending on the parameter <mark>param</mark> CHARMM19
## or CHARMM22 parameters are loaded. The parameter may
## be set through a hash-type key=>value pair argument.

sub loadParameters {
  my $self=shift;

  $self->_getpar(@_);

  my $havertf=0;
  my $havepara=0;

  my ($topfile,$parfile);

  $self->_sendCommand("bomlev $self->{par}->{bomlev}") 
    if (defined $self->{par}->{bomlev});

  $self->_sendCommand("prnlev $self->{par}->{prnlev}") 
    if (defined $self->{par}->{prnlev});


  if ($self->{par}->{param} eq "19") {
    $topfile="$datadir/toph19.inp";
    $parfile="$datadir/param19.inp";
  } elsif ($self->{par}->{param} eq "eef1") {
    $topfile="$datadir/toph19_eef1.inp";
    $parfile="$datadir/param19_eef1.inp";
  } elsif ($self->{par}->{param} eq "22") {
    if ($self->{par}->{cmap} && 
	-r "$datadir/top_all22_prot_cmap.inp" &&
        -r "$datadir/par_all22_prot_cmap.inp") {
      $topfile="$datadir/top_all22_prot_cmap.inp";
      $parfile="$datadir/par_all22_prot_cmap.inp";
    } else {
      $topfile="$datadir/top_all22_prot.inp";
      $parfile="$datadir/par_all22_prot.inp";
    }
  } elsif ($self->{par}->{param} eq "opls") {
    $topfile="$datadir/opls.rtf";
    $parfile="$datadir/opls.prm";
  } elsif ($self->{par}->{param} eq "a94") {
    $topfile="$datadir/cornell_all.rtf";
    $parfile="$datadir/cornell_all.prm";
  } elsif ($self->{par}->{param} eq "lpdb") {
    $topfile="$datadir/aminoh.rtf-lpdb";
    $parfile="$datadir/parm.prm-lpdb";
  } elsif ($self->{par}->{param} eq "simdb") {
    $topfile="$datadir/simdb.rtf";
    $parfile="$datadir/simdb.prm";
  } elsif ($self->{par}->{param} eq "27") {
    if ($self->{par}->{cmap} && 
	-r "$datadir/top_all27_prot_na_cmap.rtf" &&
        -r "$datadir/par_all27_prot_na_cmap.prm") {
      $topfile="$datadir/top_all27_prot_na_cmap.rtf";
      $parfile="$datadir/par_all27_prot_na_cmap.prm";
    } else {
      $topfile="$datadir/top_all27_prot_na.rtf";
      $parfile="$datadir/par_all27_prot_na.prm";
    }
  }    

  if (defined $topfile && -r $topfile) {
    $self->_sendCommand("open unit 10 read form name \"$topfile\"");
    $self->_sendCommand("read rtf card unit 10");
    $self->_sendCommand("close unit 10");
    $havertf=1;
  }

  if (defined $self->{par}->{xtop}) {
    foreach my $n ( split(/:/,$self->{par}->{xtop}) ) {
      $n="$datadir/$n" if (!-r $n);
      if (-r $n) {
	$self->_sendCommand("open unit 10 read form name \"$n\"");
	if ($havertf) {
	  $self->_sendCommand("read rtf card unit 10 append");
	} else {
	  $self->_sendCommand("read rtf card unit 10");
	  $havertf=1;
	}
	$self->_sendCommand("close unit 10");
      }
    }
  }

  die "No topology information available" if (!$havertf);


  if (defined $parfile && -r $parfile) { 
    $self->_sendCommand("open unit 10 read form name \"$parfile\"");
    $self->_sendCommand("read para card unit 10");
    $self->_sendCommand("close unit 10");
    $self->_sendCommand("faster on");
    $havepara=1;
  }

  if (defined $self->{par}->{xpar}) {
    foreach my $n ( split(/:/,$self->{par}->{xpar}) ) {
      $n="$datadir/$n" if (!-r $n);
      if (-r $n) {
	$self->_sendCommand("open unit 10 read form name \"$n\"");
	if ($havepara) {
	  $self->_sendCommand("read para card unit 10 append");
	} else {
	  $self->_sendCommand("bomlev -1");
	  $self->_sendCommand("read para card unit 10");
	  $havepara=1;
	}
	$self->_sendCommand("close unit 10");
      }
    }
  }

  if (defined $self->{par}->{ace} && $self->{par}->{ace} ne "0") {
    my $acefile=($self->{par}->{ace} eq "1")?"ace2parpx.inp":$self->{par}->{ace};
    $acefile="$datadir/$acefile" if (!-r $acefile);
    die "cannot read ACE parameter file $acefile" if (!-r $acefile);

    $self->_sendCommand("open read unit 10 card name \"$acefile\"");
    $self->_sendCommand("read ACEParameters card unit 10");
    $self->_sendCommand("close unit 10");
  }

  die "No parameter information available" if (!$havepara);
}

## method: clearEnergy([parameters])
## clears energy components

sub clearEnergy {
  my $self=shift;

  $self->_getpar(@_);

  $self->clearGB() if ($self->{_haveGB});
  $self->clearASP() if ($self->{_haveASP});
}    

## method: setupEnergy([parameters]) 
## setups all energy components

sub setupEnergy {
  my $self=shift;

  $self->_getpar(@_);

  $self->{par}->{cuton}=$self->{par}->{cutoff}-2
    if (defined $self->{par}->{cutoff} && !defined $self->{par}->{cuton});
  $self->{par}->{cutnb}=$self->{par}->{cutoff}+3
    if (defined $self->{par}->{cutoff} && !defined $self->{par}->{cutnb});

  if ($self->{explicitWater} && $self->{par}->{explicit}) {
    if (!$self->{_setupPeriodic} && $self->{par}->{periodic}) {
      if (!defined $self->{par}->{cuton} || !defined $self->{par}->{cutoff} || 
	  !defined $self->{par}->{cutnb}) {
	$self->{par}->{cutnb}=10;
	$self->{par}->{cutoff}=9;
	$self->{par}->{cuton}=8;
      }
      $self->setupNonBonded();
      $self->periodicBoundaries();
      $self->ewald() if ($self->{par}->{ewald});
      $self->{par}->{dyntrfrq}=200 if ($self->{par}->{ewald});
      if (!defined $self->{par}->{dynnose}) {
	$self->{par}->{dynnose}=1;
	$self->{par}->{dynnoseq}=1000.0;
      }
      $self->{_setupPeriodic}=1;
    } else {
      if (!defined $self->{par}->{cuton} || !defined $self->{par}->{cutoff} || 
	  !defined $self->{par}->{cutnb}) {
	$self->{par}->{cutnb}=20;
	$self->{par}->{cutoff}=18;
	$self->{par}->{cuton}=16;
      }
      $self->setupNonBonded();
    }
  } else {
    if (!defined $self->{par}->{cuton} || !defined $self->{par}->{cutoff} || 
	!defined $self->{par}->{cutnb}) {
      $self->{par}->{cutnb}=20;
      $self->{par}->{cutoff}=18;
      $self->{par}->{cuton}=16;
    }
    if (defined $self->{par}->{gb} && $self->{par}->{gb} ne "0" && !$self->{_haveGB}) {
      $self->setupGB();
    } elsif ($self->{_haveASP}) {
      $self->_sendCommand("skip excl asp");
    } elsif ($self->{par}->{asp}) {
      $self->setupASP();
    } elsif ($self->{par}->{eef1} && !$self->{_haveASP}) {
      $self->setupEEF1();
    } elsif ($self->{par}->{sasa} && !$self->{_haveASP}) {
      $self->setupSASA();
    } elsif ($self->{par}->{og} && !$self->{_haveOG}) {
      $self->setupOG();
    } elsif (defined $self->{par}->{ace} && $self->{par}->{ace} ne "0" && !$self->{_haveACE}) {
      $self->setupACE();
    } else {
      $self->setupNonBonded();
    }
  }

  if (!$self->{par}->{user}) {
    $self->_sendCommand("skip user");
  }
#  if ($self->{par}->{user}) {
#    system "cp $datadir/ala.map ." if (!-r "ala.map");
#    system "cp $datadir/gly.map ." if (!-r "gly.map");
#  } else {
#    $self->_sendCommand("skip user");
#  }

#  $self->_sendCommand("skip cmap") unless ($self->{par}->{cmap});
}

## method: setupFromPSF(psffile,crdfile)
## reads in a PSF structure and reads coordinates from
## a CHARMM coordinate file

sub setupFromPSF {
  my $self=shift;
  my $psf=shift;
  my $crd=shift;

  $self->_sendCommand("open unit 10 read form name \"$psf\"");
  $self->_sendCommand("read psf card unit 10");

  my $mol=Molecule::new();

  if (defined $crd) {
    $self->_sendCommand("open unit 10 read form name \"$crd\"");
    $self->_sendCommand("read coor card unit 10");
    $mol->readCRD($crd);
  } else {
    $mol->readPSF($psf); 
  }

  $self->{molecule}=$mol;

  foreach my $c ( @{$mol->{chain}} ) {
    for my $r ( @{$c->{res}} ) {
      if ($r->{name} eq "TIP3" || $r->{name} eq "HOH"){ 
	$self->{explicitWater}=1;
      }
    }
  }

  $self->{_setup}=1;
}

## method: setupFromPDB(pdbfile[,terminalgroup])
## generates the CHARMM PSF structure from a protein structure
## in a PDB file. It also rebuilds missing atoms and adds
## hydrogens, if necessary. The type of terminal groups
## may be selected through the second argument. Possible
## options are <mark>none</mark>, <mark>ends</mark>, and
## <mark>all</mark>.

sub setupFromPDB {
  my $self=shift;
  my $pdb=shift;
  my $terminal=shift;
  my $noSS=shift;
  my $findChainBreaks=shift;

  $findChainBreaks=0 if (! defined $findChainBreaks);

  my $mol=Molecule::new();
  $mol->readPDB($pdb);
  $mol->findSSBonds() if ((! defined $noSS) || (! $noSS));
  $mol->fixHistidine($self->{par}->{hsd},$self->{par}->{hse},$self->{par}->{hsp});
  $mol->changeResName($self->{par}->{resmod});

  if ($findChainBreaks) {
      $mol->generateSplitSegNames();
  } else {
      $mol->generateSegNames();
  }

  $self->setupFromMolecule($mol,$terminal);
}

## method: setupFromMolecule(mol[,terminalgroup])
## generates the CHARMM PSF structure from a protein structure
## in a Molecule object. 

sub setupFromMolecule {
  my $self=shift;
  my $mol=shift;
  my $terminal=shift;
  my $readonly=shift;

  $readonly=0 if (!defined $readonly);

  $self->{molecule}=$mol;

  my $slist=$self->{molecule}->getSegNames();

  die "no segments found in PDB"
    if (!defined $slist || $#{$slist}<0);

  if (!defined $terminal) {
       if (defined $self->{par}->{blocked} && $self->{par}->{blocked}) {
         $terminal="blocked";
       } else {
         $terminal="all";
       }
  } elsif ($terminal eq "auto") {
       $terminal="all";
       my $c=$mol->getChain($slist->[0]->{chain});
       my $a=$c->{atom};
       my $tres=$mol->getResidueInChain($slist->[0]->{from},$c);
       
       for (my $ti=$tres->{start}; $ti<=$tres->{end}; $ti++) {
         $terminal="blocked" 
           if ($a->[$ti]->{atomname} =~ /CAY|CAT|CY|NT|CLM|CAM|OAM/);   
       }
  } 

  if ($terminal eq "blocked" && !$readonly) {
    if ($self->{par}->{param} eq "19") {
      my $patch="read rtf card append\n";
      $patch.="* patches\n*\n";
      $patch.="   20    1\n";
      $patch.="PRES AMNP         0.00 \n";
      $patch.="GROUP                  \n";
      $patch.="ATOM CLM  CH3E    0.0  \n";
      $patch.="GROUP                  \n";
      $patch.="ATOM CAM  C       0.55 \n";
      $patch.="ATOM OAM  O      -0.55 \n";
      $patch.="BOND CLM  CAM       CAM  N         CAM  OAM\n";
      $patch.="DIHE C    CA   N    CAM  CA   N    CAM  CLM\n";
      $patch.="IMPH CAM  CLM  N    OAM\n";
      $patch.="IMPH N    CAM  CA   H  \n";
      $patch.="ACCE OAM CAM\n";
      $patch.="IC   C    CA   N    CAM   0.0000    0.00  180.00    0.00   0.000 \n";
      $patch.="IC   CA   N    CAM  CLM   0.0000    0.00  180.00    0.00   0.000 \n";
      $patch.="IC   N    CLM  *CAM OAM   0.0000    0.00  180.00    0.00   0.000 \n";
      $patch.="IC   CA   CAM  *N   H     0.0000    0.00  180.00    0.00   0.000 \n";
      $patch.="\n";
      $patch.="PRES CBXP         0.00\n";
      $patch.="GROUP                 \n";
      $patch.="ATOM NCB  NH1    -0.35\n";
      $patch.="ATOM HCB  H       0.25\n";
      $patch.="ATOM CAC  CH3E    0.10\n";
      $patch.="BOND NCB  CAC       NCB  HCB       C    NCB\n";
      $patch.="DIHE N    CA   C    NCB\n";
      $patch.="DIHE CA   C    NCB  CAC\n";
      $patch.="IMPH NCB  C    CAC  HCB\n";
      $patch.="IMPH C    CA   NCB  O\n";
      $patch.="DONO HCB  NCB\n";
      $patch.="IC   C    CAC  *NCB HCB    0.0000    0.00  180.00    0.00   0.0000\n";
      $patch.="IC   N    CA   C    NCB    0.0000    0.00  180.00    0.00   0.0000\n";
      $patch.="IC   CA   C    NCB  CAC    0.0000    0.00  180.00    0.00   0.0000\n";
      $patch.="IC   NCB  CA   *C   O      0.0000    0.00  180.00    0.00   0.0000\n";
      $patch.="END\n";
      $self->_sendCommand($patch);
    }
  }

  my @ppredo=();
  for (my $if=0; $if<=$#{$slist}; $if++) {
    my $first="";
    my $last="";
    my $f=$slist->[$if];

    if ($f->{name} =~ /^N/) {
      $first=(defined $self->{par}->{nter})?"first $self->{par}->{nter}":"first 5TER" ;
      $last=(defined $self->{par}->{cter})?"last $self->{par}->{cter}":"last 3TER" ;
    } elsif ($terminal eq "none") {
      $first="firs none";
      $last="last none";
    } elsif ($terminal eq "ends") {
      $first="firs none" if (!$f->{first});
      $last="last none" if (!$f->{last});
    } elsif ($terminal eq "blocked") {
      if (defined $self->{par}->{nter} && defined $self->{par}->{cter}) {
	$first="firs $self->{par}->{nter}";
	$last="last $self->{par}->{cter}";
      } elsif ($self->{par}->{param} eq "19") {
	$first="firs amnp";
	$last="last cbxp";
      } elsif ($self->{par}->{param} =~ /22/) {
	$first="firs ace";
	$last="last ct3";
      }	
    }

    my $fname=hostname."-pdb$$-".$f->{name};

    $self->{molecule}->setValidSegment($f->{name});
    my $segmol=$self->{molecule}->clone(1);

    my $hetero=($segmol->{chain}->[0]->{id} eq "+")?1:0;

    $segmol->setChain(" ") if ($hetero);

    $segmol->fixCOO($self->{par}->{blocked}) if (!$hetero);
    $segmol->writePDB($fname,translate=>getConvType($self->{par}->{param}),ssbond=>0);

    if (!$readonly) {
      $self->_sendCommand("open unit 10 read form name \"$fname\"");
      $self->_sendCommand("read sequ pdb unit 10");

      if ($f->{name} eq "TIP3") {
	$self->_sendCommand("generate tip3 setup noangl nodihe");
	$self->{explicitWater}=1;
      } elsif ($hetero) {
	$self->_sendCommand("generate $f->{name} setup warn");
      } else {
	$self->_sendCommand("generate $first $last $f->{name} setup warn");
      }
      $self->_sendCommand("close unit 10");

      if (defined $self->{par}->{patch} && $self->{par}->{patch} ne "") {
	foreach my $pp ( split(/_/,$self->{par}->{patch}) ) {
	  my @ppl=split(/:/,$pp);
	  my $pcmd="patch $ppl[0] ";
	  my $notok=0;
	  for (my $ippl=1; $ippl<=$#ppl; $ippl++) {
	    my $flist=&GenUtil::fragListFromOption($ppl[$ippl]);
	    my $fres=$segmol->getResidue($flist->[0]->{from},$flist->[0]->{chain});
	    if (defined $fres) {
	      $pcmd.=" " if ($ippl>1);
	      $pcmd.="$fres->{seg} $fres->{num}";
	    } else {
	      $notok=1;
	    }
	  }
	  if ($notok) {
	    push(@ppredo,$pp);
	  } else {
	    $self->_sendCommand($pcmd);
	  }
	}
      } elsif ($f->{name} =~ /^N/ && (!defined $self->{par}->{nodeoxy} || !$self->{par}->{nodeoxy}) && $f->{name}!~/R/) {
	foreach my $r ( @{$segmol->{chain}->[0]->{res}} ) {
	  $self->_sendCommand(sprintf("patch deo%d %s %d",
				      ($r->{name}=~/GUA|ADE/)?2:1,
				      $f->{name},$r->{num}));
	}
      }
    }
  }

  foreach my $pp ( @ppredo ) {
    my @ppl=split(/:/,$pp);
    my $pcmd="patch $ppl[0] ";
    my $notok=0;
    for (my $ippl=1; $ippl<=$#ppl; $ippl++) {
      my $flist=&GenUtil::fragListFromOption($ppl[$ippl]);
      my $c=$mol->getChain($flist->[0]->{chain});
      my $fres=$mol->getResidueInChain($flist->[0]->{from},$c);
      if (defined $fres) {
	$pcmd.=" " if ($ippl>1);
	$pcmd.="$fres->{seg} $fres->{num}";
      } else {
	$notok=1;
      }
    }
    $self->_sendCommand($pcmd) unless ($notok);
  }

  foreach my $ss ( @{$mol->{ssbond}} ) {
    my $c1=$mol->getChain($ss->{chain1});
    my $c2=$mol->getChain($ss->{chain2});
    if (defined $c1 && defined $c2) {
      my $r1=$mol->getResidueInChain($ss->{resnum1},$c1);
      my $r2=$mol->getResidueInChain($ss->{resnum2},$c2);
      if (defined $r1 && defined $r2) {
	my $pcmd=sprintf("patch DISU %s %d %s %d\n",$r1->{seg},$r1->{num},$r2->{seg},$r2->{num});
	$self->_sendCommand($pcmd);
      }
    }
  }

  foreach my $f ( @{$slist} ) {
    my $fname=hostname."-pdb$$-".$f->{name};

    $self->_sendCommand("open unit 10 read form name \"$fname\"");
    $self->_sendCommand("read coor pdb unit 10 resi");
    $self->_sendCommand("close unit 10");

    &GenUtil::remove($fname) unless (defined $self->{handle}->{cmdlog});
  }

  $self->_sendCommand("bomlev -2");

  if (!$self->{explicitWater} && $self->{par}->{buildall}) {
#  $self->_sendCommand("ic fill ");
    $self->_sendCommand("ic param");
    if (($self->{par}->{param} eq "19" || $self->{par}->{param} =~ /22/ || 
	 $self->{par}->{param} eq "eef1")) {
      $self->_sendCommand("coor copy comp");
      $self->_sendCommand("ic build comp");
      $self->_sendCommand("coor copy select .not. hydrogen end");
      $self->_sendCommand("hbuild atom cdie eps 80.0 cutnb 10.0 ctofnb 7.5 ctonnb 6.5 shift vshift bygr");
    } else {
      $self->_sendCommand("ic build");
    }
  } else {
    $self->_sendCommand("hbuild atom cdie eps 80.0 cutnb 10.0 ctofnb 7.5 ctonnb 6.5 shift vshift bygr");
  }

  $self->{_setup}=1;
}

## method: readFromPDB(pdbfile)
## reads coordinates for a protein structure from a PDB file
## if the PSF was setup before with <mark>setupFromMolecule</mark>.
## This is useful for repeating a modeling procedure
## for different conformations of the same protein.
## As <mark>setupFromMolecule</mark> this method also rebuilds
## missing atoms and adds hydrogens if necessary.

sub OldreadFromPDB {
  my $self=shift;
  my $pdb=shift;

  die "need to setup protein structure first"
    if (!defined $self->{molecule});

  my $mol=Molecule::new();
  $mol->readPDB($pdb);
  $mol->fixHistidine($self->{par}->{hsd},$self->{par}->{hse},$self->{par}->{hsp});
  $mol->changeResName($self->{par}->{resmod});
  $mol->copySegNames($self->{molecule});

  my $fname=hostname."-pdb$$";

  $mol->fixCOO($self->{par}->{blocked});
  $mol->writePDB($fname,translate=>getConvType($self->{par}->{param}));

  $self->_sendCommand("open unit 10 read form name \"$fname\"");
  $self->_sendCommand("read coor pdb unit 10 resi");
  $self->_sendCommand("close unit 10");
#  $self->_sendCommand("ic fill");
  $self->_sendCommand("ic param");
  $self->_sendCommand("coor copy comp");
  $self->_sendCommand("ic build comp");
  $self->_sendCommand("coor copy select .not. hydrogen end");
  $self->_sendCommand("hbuild atom cdie eps 80.0 cutnb 10.0 ctofnb 7.5 ctonnb 6.5 shift vshift bygr");

  &GenUtil::remove($fname) unless (defined $self->{handle}->{cmdlog});
}

## method: readFromPDB(pdbfile)
## reads coordinates for a protein structure from a PDB file
## if the PSF was setup before with <mark>setupFromMolecule</mark>.
## This is useful for repeating a modeling procedure
## for different conformations of the same protein.
## As <mark>setupFromMolecule</mark> this method also rebuilds
## missing atoms and adds hydrogens if necessary.

sub readFromPDB {
  my $self=shift;
  my $pdb=shift;
  my $terminal=shift;

  die "need to setup protein structure first"
    if (!defined $self->{molecule});

  my $mol=Molecule::new();
  $mol->readPDB($pdb);
  $mol->fixHistidine($self->{par}->{hsd},$self->{par}->{hse},$self->{par}->{hsp});
  $mol->changeResName($self->{par}->{resmod});
  $mol->copySegNames($self->{molecule});

  $self->setupFromMolecule($mol,$terminal,1);
}

## method: loadReference(molecule)
## reads coordinates for a protein structure from
## a PDB file into the alternate COMP coordinate set

sub loadReference {
  my $self=shift;
  my $pdb =shift;

  die "need to setup protein structure from PDB first"
    if (!defined $self->{molecule});

  my $mol=Molecule::new();
  $mol->readPDB($pdb);
  $mol->fixHistidine($self->{par}->{hsd},$self->{par}->{hse},$self->{par}->{hsp});
  $mol->changeResName($self->{par}->{resmod});
  $mol->copySegNames($self->{molecule});

  $self->_sendCommand("coor set xdir 1.0 dist 0.0 comp select all end");

  my $slist=$mol->getSegNames();

  for (my $if=0; $if<=$#{$slist}; $if++) {
    my $first="";
    my $last="";
    my $f=$slist->[$if];

    my $fname=hostname."-pdb$$-".$f->{name};

    $mol->setValidSegment($f->{name});
    my $segmol=$mol->clone(1);

    my $hetero=($segmol->{chain}->[0]->{id} eq "+")?1:0;

    $segmol->setChain(" ") if ($hetero);

    $segmol->fixCOO($self->{par}->{blocked}) if (!$hetero);
    $segmol->writePDB($fname,translate=>getConvType($self->{par}->{param}),ssbond=>0);

    $self->_sendCommand("open unit 10 read form name \"$fname\"");
    $self->_sendCommand("read coor pdb unit 10 comp resi");
    $self->_sendCommand("close unit 10");

    &GenUtil::remove($fname) unless (defined $self->{handle}->{cmdlog});
  }
}

## method: initCoordinates()
## initializes all coordinates

sub initCoordinates {
  my $self=shift;
  $self->_sendCommand("coor init select all end");
}

## method: setupNonBonded([parameters])
## calls the CHARMM command <mark>update</mark> 
## to set non-bonded interaction parameters. 
## Options are constant (dielec=>"CDIE") or 
## distance-dependent dielectric (dielec=>"RDIE"),
## the dieletric constant epsilon and interaction
## and list cutoffs (cutnb, cuton, cutoff).
## Force shifting is set for electrostatics as well
## as van der Waals interactions as the default
## for vacuum minimizations.

sub setupNonBonded {
  my $self=shift;

  $self->_getpar(@_);

  if ($self->{par}->{cut}) {
    my $trunc="shift vshift";
    if (defined $self->{par}->{trunc} && $self->{par}->{trunc} eq "switch") {
      $trunc="switch vswitch";
    }
    $self->_sendCommand("update atom $self->{par}->{dielec} eps $self->{par}->{epsilon} cutnb $self->{par}->{cutnb} ctofnb $self->{par}->{cutoff} ctonnb $self->{par}->{cuton} $trunc bygr");
  } else {
    my $trunc="switch vswitch";
    if (defined $self->{par}->{trunc} && $self->{par}->{trunc} eq "shift") {
      $trunc="shift vshift";
    }
    $self->_sendCommand("update atom $self->{par}->{dielec} eps $self->{par}->{epsilon} cutnb 9999.0 ctofnb 9000.0 ctonnb 9000.0 $trunc bygr");
  }
}

## method: setupEEF1([parameters]) 
## sets up EEF1 

sub setupEEF1 {
  my $self=shift;
  
  $self->_getpar(@_);

  my $fname;
  if (-r $self->{par}->{eef1file}) {
    $fname=$self->{par}->{eef1file};
  } elsif (-r "$datadir/$self->{par}->{eef1file}") {
    $fname="$datadir/$self->{par}->{eef1file}";
  } elsif (-r "$ENV{MMTSBDIR}/data/$self->{par}->{eef1file}") {
    $fname="$ENV{MMTSBDIR}/data/$self->{par}->{eef1file}";
  } else {
    die "cannot find EEF1 input file $self->{par}->{eef1file}";
  }

  $self->_sendCommand("eef1 setup temp $self->{par}->{dyntemp} unit 93 name \"$fname\"");
  $self->_sendCommand("update ctonnb 7.0 ctofnb 9.0 cutnb 10.0 group rdie switch vswitch");

  $self->{_haveASP}=1;
}

## method:: setupSASA([parameters])
## sets up surface-area based fast implict solvation model

sub setupSASA {
  my $self=shift;
  $self->_getpar(@_);

  $self->_sendCommand("nbond nbxmod 5 atom rdiel shift vatom vdistance vshift -\ncutnb 8.0 ctofnb 7.5 ctonnb 6.5 eps $self->{par}->{sasaeps} e14fac 0.4 wmin 1.5\n");
#  if ($self->{par}->{sasasig3}>0) {
#    $self->_sendCommand("sasa selection all end SIG1 $self->{par}->{sasasig1} SIG2 $self->{par}->{sasasig2} SIG3 $self->{par}->{sasasig3}\n");
#  } else {
    $self->_sendCommand("sasa selection .not. hydrogen end SIG1 $self->{par}->{sasasig1} SIG2 $self->{par}->{sasasig2}\n");
#  }
}

## method: setupHBond([parameters]) 
## sets up hydrogen bonding list

sub setupHBond {
  my $self=shift;
  $self->_getpar(@_);
  $self->_sendCommand("bomlev -2");
  $self->_sendCommand("hbond hbnoexclusions all -\ncuthb 6.0 ctofhb 5.0 ctonhb 4.0 cutha 100.0 ctofha 90.0 ctonha 90.0");
  $self->_sendCommand("bomlev 0");
}


## method: setupOG([parameters]) 
## sets up Olgun's scoring function

sub setupOG {
  my $self=shift;

  $self->_getpar(@_);
  $self->setupNonBonded();
  $self->_sendCommand("skip all excl dihe user bond angle urey impr");

  $self->{_haveOG}=1;
}  

## method: setupACE([parameters]) 
## sets up ACE

sub setupACE {
  my $self=shift;

  $self->_getpar(@_);

  my $sig=(defined $self->{par}->{acesigma})?sprintf("sigma %f",$self->{par}->{acesigma}):"";  
  if ($self->{par}->{cut}) {
    $self->_sendCommand("update atom ace2 ieps $self->{par}->{aceieps} seps $self->{par}->{aceseps} alpha $self->{par}->{acealpha} $sig switch vdis vswi cutnb $self->{par}->{cutnb} ctofnb $self->{par}->{cutoff} ctonnb $self->{par}->{cuton} FVSCale $self->{par}->{acevscale}");
  } else {
    $self->_sendCommand("update atom ace2 ieps $self->{par}->{aceieps} seps $self->{par}->{aceseps} alpha $self->{par}->{acealpha} $sig switch vdis vswi cutnb 9999.0 ctofnb 9000.0 ctonnb 8000.0 FVSCale $self->{par}->{acevscale}");
  }

  $self->{_haveACE}=1;
}  


## method: setupASP([parameters])
## sets up atomic solvation parameter based energy function

sub setupASP {
  my $self=shift;
  
  $self->_getpar(@_);

  $self->setupNonBonded();

  my $fname;
  if (defined $self->{par}->{aspfile}) {
    if (-r $self->{par}->{aspfile}) {
      $fname=$self->{par}->{aspfile};
    } elsif (-r "$datadir/$self->{par}->{aspfile}") {
      $fname="$datadir/$self->{par}->{aspfile}";
    } elsif (-r "$ENV{MMTSBDIR}/data/$self->{par}->{aspfile}") {
      $fname="$ENV{MMTSBDIR}/data/$self->{par}->{aspfile}";
    }
  }
  
  my $rmfile=0;
  if (!defined $fname) {
    $fname=hostname."-asp$$";
    $rmfile=1;
    
    my $out=&GenUtil::getOutputFile($fname);
    
    $self->_sendCommand("scalar radius show");
    $self->_getCHARMMOutput("PERLDONE\n  \n");
    
    my $data=();
    foreach my $l ( split(/\n/,$self->{_lastOutput}) ) {
      if ($l=~/\( .* \) /) {
	$l=~s/^ +//;
	my @fl=split(/ +/,$l);
	my $rec={};
	$rec->{seg}=$fl[1];
	$rec->{resname}=$fl[2];
	$rec->{resnum}=$fl[3];
	$rec->{atomname}=$fl[4];
	$rec->{radius}=$fl[6];
	push(@{$data},$rec);
      }
    }
    
    my %have;
    
    printf $out "* ASP input\n*\n%f\n",1.4;
    foreach my $d ( @{$data} ) {
      my $aspval=$self->{par}->{aspval};
      
      $aspval=$self->{par}->{aspvalc} 
	if (defined $self->{par}->{aspvalc} && $d->{atomname}=~/^[0-9]*C/);
      $aspval=$self->{par}->{aspvaln} 
	if (defined $self->{par}->{aspvaln} && $d->{atomname}=~/^[0-9]*N/);
      $aspval=$self->{par}->{aspvalo} 
	if (defined $self->{par}->{aspvalo} && $d->{atomname}=~/^[0-9]*O/);
      $aspval=$self->{par}->{aspvalh} 
	if (defined $self->{par}->{aspvalh} && $d->{atomname}=~/^[0-9]*H/);
      $aspval=$self->{par}->{aspvals} 
	if (defined $self->{par}->{aspvals} && $d->{atomname}=~/^[0-9]*S/);
      
      my $asprad=$d->{radius};
      
      $asprad=$self->{par}->{aspradh}
	if (defined $self->{par}->{aspradh} && $d->{atomname}=~/^[0-9]*H/);
      
      printf $out "%s %s %f %f %f\n",
	$d->{resname},$d->{atomname},$aspval,$asprad,$self->{par}->{aspref}
	  if (!defined $have{"$d->{resname}:$d->{atomname}"});
      $have{"$d->{resname}:$d->{atomname}"}=1;
    }
    printf $out "END\n";
    undef $out;
  }

  $self->_sendCommand("open unit 29 read form name \"$fname\"");
  $self->_sendCommand("read surf unit 29");
#    $self->_sendCommand("close unit 29");
    
  $self->_getCHARMMOutput("PERLDONE\n  \n");
  &GenUtil::remove($fname) if ($rmfile && !defined $self->{handle}->{cmdlog});
  
  $self->{_haveASP}=1;
}

## method: clearASP()
## switch off ASP energy terms

sub clearASP {
  my $self=shift;

  if (defined $self->{_haveASP}) {
    $self->_sendCommand("skip asp");
  }
}


## method: _scaleRadii(name)
## scales radii for GB/PB

sub _scaleRadii {
  my $self=shift;
  my $tag=shift;

  return if (!defined $self->{par}->{scalerad} || !$self->{par}->{scalerad});

  if (-r $self->{par}->{scalerad}) {
    open INP,"$self->{par}->{scalerad}";
    while (<INP>) {
      chomp;
      $self->_sendCommand($_);
    }
    close INP;
  } elsif ($self->{par}->{scalerad} =~ /scal/) { 
    foreach my $n ( split(/:/,$self->{par}->{scalerad}) ) {
      $self->_sendCommand($n);
    }
  } else {
    if ($self->{par}->{scalerad} eq "bondi") {
      $self->_sendCommand(sprintf("scalar wmain set 1.70 sele type C* end"));
      $self->_sendCommand(sprintf("scalar wmain set 1.55 sele type N* end"));
      $self->_sendCommand(sprintf("scalar wmain set 1.20 sele type H* end"));
      $self->_sendCommand(sprintf("scalar wmain set 1.80 sele type S* end"));
      $self->_sendCommand(sprintf("scalar wmain set 1.50 sele type O* end"));
    } elsif ($self->{par}->{scalerad} eq "mbondi") {
      $self->_sendCommand(sprintf("scalar wmain set 1.70 sele type C* end"));
      $self->_sendCommand(sprintf("scalar wmain set 1.55 sele type N* end"));
      $self->_sendCommand(sprintf("scalar wmain set 1.80 sele type S* end"));
      $self->_sendCommand(sprintf("scalar wmain set 1.50 sele type O* end"));
      $self->_sendCommand(sprintf("scalar wmain set 1.30 sele type H* end"));
      $self->_sendCommand(sprintf("scalar wmain set 0.80 sele type HG1 .or. type HH end"));
    } elsif ($self->{par}->{scalerad} eq "rebel") {
      $self->_sendCommand(sprintf("scalar wmain set 1.0 sele prop wmain .lt. 1.0 end"));
    } elsif ($self->{par}->{scalerad} eq "amber") {
      $self->_sendCommand(sprintf("scalar wmain mult 0.890899 sele all end"));
    } elsif ($self->{par}->{scalerad} eq "impact") {
      $self->_sendCommand(sprintf("scalar wmain mult 0.890899 sele all end"));
      $self->_sendCommand(sprintf("scalar wmain set 1.0 sele prop wmain .lt. 0.01 end"));
    } else { 
      if ($self->{par}->{param} =~ /19/) {
      } elsif ($self->{par}->{param} =~ /22/) {
	if ($self->{par}->{scalerad} =~/mfn/) {
	  my %amfn;
	  my %mfndef = (    mfnh     => 0.0,
			     mfn3     => 2.10,
			     mfn2     => 2.50,
			     mfnb10   => 2.50,
			     mfnb9    => 2.75,
			     mfnb6    => 2.00,
			     mfnb7    => 1.51,
			     mfnb8    => 2.20,
			     mfn5     => 2.40,
			     mfn13    => 1.83,
			     mfn12    => 1.65,
			     mfn4     => 2.65,
			     mfn19    => 1.90,
			     mfn18    => 1.94,
			     mfn7     => 1.70,
			     mfn16    => 2.05,
			     mfn8     => 1.75,
			     mfn17    => 1.85,
			     mfn25    => 1.00,
			     mfn6     => 1.95,
			     mfn15    => 2.10,
			     mfn11    => 1.48,
			     mfn10    => 1.53,
			     mfn9     => 2.60,
			     mfn14    => 2.05,
			     mfn21    => 1.85 );

	  foreach my $k ( keys %mfndef ) {
	    $amfn{$k}=(defined $self->{par}->{$k})?$self->{par}->{$k}:$mfndef{$k};
	  }

	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn3} sele (type CAY .or. type CAT) end"));
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfnb6} sele type CY end"));
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfnb7} sele type OY end"));
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfnb8} sele type NT end"));
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn11} sele type OT* end"));
	
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfnb6} sele type C  end"));
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfnb7} sele type O  end"));
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfnb8} sele type N  end"));
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfnb9} sele type CA  end"));
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfnb10} sele (resnam GLY .and. type CA) end"));

	  if ($self->{par}->{scalerad}=~/mfns/) {
	    my $snback=$self->{par}->{sbackn}; 
	    my $soback=$self->{par}->{sbacko}; 
	    my $scaback=$self->{par}->{sbackca}; 
	    $self->_sendCommand(sprintf("scalar wmain mult $soback sele type O  end"));
	    $self->_sendCommand(sprintf("scalar wmain mult $snback sele type N  end"));
	    $self->_sendCommand(sprintf("scalar wmain mult $scaback sele type CA  end"));
	  }

	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfnh} sele type H* end"));

	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn2} sele type CB end"));
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn3} sele type CB .and. resname ALA end"));
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn4} sele type CB .and. resname SER end"));

	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn3} sele type CG* .and. resname VAL end "));
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn3} sele type CD* .and. resname LEU end "));
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn3} sele type CG2 .and. resname ILE end "));
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn3} sele type CD .and. resname ILE end "));
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn3} sele type CG2 .and. resname THR end "));
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn3} sele type CE .and. resname MET end "));

	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn5} sele type CG .and. resname LEU end "));
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn5} sele type CG .and. resname GLN end "));
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn5} sele type CG .and. resname GLU end "));
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn5} sele type CG1 .and. resname ILE end "));
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn5} sele (type CG .or. type CD) .and. resname LYS end "));
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn5} sele (type CG .or. type CD) .and. resname ARG end "));


	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn19} sele (resnam TYR .and. (type CG .or. type CE* .or. type CD* .or. -\n  type CZ)) .or. (resnam PHE .and. (type CG .or. type CE* .or. - \n  type CD* .or. type CZ))  end"));
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn13} sele (resnam TYR .and. type OH) end"));
	  
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn18} sele (resnam PRO .and. (type CB .or. type CG .or. type CD)) end"));
	  
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn7} sele (resnam TRP .and. (type CG .or. type CE* .or. type CD* .or. - \n  type CZ* .or. type CH2)) end"));
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn16} sele resnam TRP .and. type NE1 end"));

	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn12} sele type OG* end"));	  

	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn25} sele type SD .and. resname MET end "));
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn8} sele type CG .and. resname MET end "));

	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn17} sele type SG .and. resname CYS end"));

	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn6} sele (resnam GLN .and. type CD) .or. (resnam ASN .and. type CG) .or. - \n  (resnam GLU .and. type CD) .or. (resnam ASP .and. type CG) end  "));
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn11} sele (resnam ASN .or. resnam GLN) .and. (type OE* .or. type OD*) end"));
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn15} sele (resnam GLN .and. type NE2) .or. (resnam ASN .and. type ND2) end "));
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn10} sele (resnam GLU .or. resnam ASP) .and. (type OE* .or. type OD*) end"));

	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn9} sele (resnam ARG .and. type CZ) .or. (resnam LYS .and. type CE) end"));
	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn14} sele resnam ARG .and. (type NH* .or. type NE) .or. -\n (resnam LYS .and. type NZ) end"));

	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn19} sele (resnam HSD .or. resnam HSE .or. resnam HSP) .and. (type CE1 .or. type CD2 .or. type CG) end"));

	  $self->_sendCommand(sprintf("scalar wmain set $amfn{mfn21} sele (resnam HSD .or. resnam HSP .or. resnam HSE) .and. (type NE2 .or. type ND1) end"));

	} elsif ($self->{par}->{scalerad} =~/nina/) {
	  my %anina;
	  my %mninadef = ( nina14  => 2.09,
			   nina10  => 1.37,
			   nina11  => 1.39,
			   nina9   => 2.75,
			   nina7   => 1.75,
			   nina16  => 2.35,
			   nina12  => 1.61,
			   nina3   => 2.41,
			   ninah   => 0.0,
			   ninab1  => 2.02,
			   ninab2  => 2.00,
			   ninab3  => 1.49,
			   ninab4  => 2.19,
			   ninab5  => 1.37,
			   ninab6  => 2.00,
			   ninab7  => 1.49,
			   ninab8  => 2.19,
			   ninab9  => 2.80,
			   ninab10 => 2.33,
			   nina2   => 2.62,
			   nina4   => 2.72,
			   nina5   => 2.39,
			   nina6   => 1.94,
			   nina18  => 1.94,
			   nina19  => 1.96,
			   nina8   => 2.06,
			   nina20  => 1.94,
			   nina13  => 1.81,
			   nina21  => 1.76,
			   nina22  => 2.25,
			   nina15  => 2.11,
			   nina17  => 1.96 );

          my %ninadef  = ( nina14  => 2.13,
			   nina10  => 1.40,
			   nina11  => 1.42,
			   nina9   => 2.80,
			   nina7   => 1.78,
			   nina16  => 2.40,
			   nina12  => 1.64,
			   nina3   => 2.46,
			   ninah   => 0.0,
			   ninab1  => 2.06,
			   ninab2  => 2.04,
			   ninab3  => 1.52,
			   ninab4  => 2.23,
			   ninab5  => 1.40,
			   ninab6  => 2.04,
			   ninab7  => 1.52,
			   ninab8  => 2.23,
			   ninab9  => 2.86,
			   ninab10 => 2.38,
			   nina2   => 2.67,
			   nina4   => 2.77,
			   nina5   => 2.44,
			   nina6   => 1.98,
			   nina18  => 1.98,
			   nina19  => 2.00,
			   nina8   => 2.10,
			   nina20  => 1.98,
			   nina13  => 1.85,
			   nina21  => 1.80,
			   nina22  => 2.30,
			   nina15  => 2.15,
			   nina17  => 2.00 );

	  if ($self->{par}->{scalerad} =~/mnina/) {
	    foreach my $k ( keys %mninadef ) {
	      $anina{$k}=(defined $self->{par}->{$k})?$self->{par}->{$k}:$mninadef{$k};
	    }
	  } else {
	    foreach my $k ( keys %ninadef ) {
	      $anina{$k}=(defined $self->{par}->{$k})?$self->{par}->{$k}:$ninadef{$k};
	    }
	  }

	  $self->_sendCommand(sprintf("scalar wmain set $anina{ninab1} sele (type CAY .or. type CAT) end"));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{ninab2} sele type CY end"));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{ninab3} sele type OY end"));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{ninab4} sele type NT end"));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{ninab5} sele type OT* end"));
	
	  $self->_sendCommand(sprintf("scalar wmain set $anina{ninab6} sele type C  end"));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{ninab7} sele type O  end"));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{ninab8} sele type N  end"));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{ninab9} sele type CA  end"));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{ninab10} sele (resnam GLY .and. type CA) end"));

	  if ($self->{par}->{scalerad}=~/ninas/) {
	    my $snback=$self->{par}->{sbackn}; 
	    my $soback=$self->{par}->{sbacko}; 
	    my $scaback=$self->{par}->{sbackca}; 
	    $self->_sendCommand(sprintf("scalar wmain mult $soback sele type O  end"));
	    $self->_sendCommand(sprintf("scalar wmain mult $snback sele type N  end"));
	    $self->_sendCommand(sprintf("scalar wmain mult $scaback sele type CA  end"));
	  }

	  $self->_sendCommand(sprintf("scalar wmain set $anina{ninah} sele type H* end"));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{nina2} sele type CB end"));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{nina3} sele type CG* end "));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{nina4} sele type CG .and. resnam GLU end"));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{nina5} sele type CD* end "));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{nina6} sele (resnam GLN .and. type CD) .or. (resnam ASN .and. type CG) .or. - \n  (resnam GLU .and. type CD) .or. (resnam ASP .and. type CG) end  "));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{nina18} sele (resnam PRO .and. (type CB .or. type CG .or. type CD)) end"));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{nina19} sele (resnam TYR .and. (type CE* .or. type CD* .or. -\n  type CZ)) .or. (resnam PHE .and. (type CE* .or. - \n  type CD* .or. type CZ))  end"));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{nina7} sele (resnam TRP .and. (type CE* .or. type CD* .or. - \n  type CZ* .or. type CH2)) end"));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{nina8} sele (resname MET .and. type CE) end"));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{nina9} sele (resnam ARG .and. type CZ) .or. (resnam LYS .and. type CE) end"));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{nina20} sele (resnam HSD .and. type CE1) .or. (resnam HSD .and. type CD2) end"));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{nina10} sele (resnam GLU .or. resnam ASP) .and. (type OE* .or. type OD*) end"));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{nina11} sele (resnam ASN .or. resnam GLN) .and. (type OE* .or. type OD*) end"));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{nina12} sele type OG* end"));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{nina13} sele (resnam TYR .and. type OH) end"));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{nina21} sele resnam HSD .and. (type NE2 .or. type ND1) end"));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{nina22} sele resnam HSP .and. (type NE2 .or. type ND1) end"));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{nina14} sele resnam ARG .and. (type NH* .or. type NE) .or. -\n (resnam LYS .and. type NZ) end"));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{nina15} sele (resnam GLN .and. type NE2) .or. (resnam ASN .and. type ND2) end "));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{nina16} sele resnam TRP .and. type NE1 end"));
	  $self->_sendCommand(sprintf("scalar wmain set $anina{nina17} sele type S* end"));
	} else {
	  # backbone
	  my $snback=$self->{par}->{sbackn};  
	  my $soback=$self->{par}->{sbacko}; 
	  my $sbackca=$self->{par}->{sbackca};
	  $self->_sendCommand(sprintf("scalar $tag mult %1.4f sele type ca .and. .not. resnam GLY end",$sbackca));
	  $self->_sendCommand(sprintf("scalar $tag mult %1.4f sele type o end",$soback));	
	  $self->_sendCommand(sprintf("scalar $tag mult %1.4f sele type oy end",$soback));	
	  $self->_sendCommand(sprintf("scalar $tag mult %1.4f sele type n end",$snback));
	  $self->_sendCommand(sprintf("scalar $tag mult %1.4f sele type c end",1.0));
	  $self->_sendCommand(sprintf("scalar $tag mult %1.4f sele type ca .and. resnam GLY end",1.0));
	  $self->_sendCommand(sprintf("scalar $tag mult %1.4f sele type nt end",1.0));
	  $self->_sendCommand(sprintf("scalar $tag mult %1.4f sele type cay end",1.0));
	  $self->_sendCommand(sprintf("scalar $tag mult %1.4f sele type cat end",1.0));
	  
	  my ($soon,$snind,$snh2,$snh2p,$snh3p,$sroh,$soh,$ssulfur,$snimi,$ss2);
	  
	  $soon=0.96;   ### 0.999, 1.002;  #0.95;
	  $snind=0.84;  ### 0.84, 0.86;  #0.90
	  $snh2=0.925;  ### 0.915, 0.927;  #0.90
	  $snh2p=1.010; ### 1.014, 1.022; #1.01
	  $snh3p=1.01;  ### 1.08, 1.087; #1.08
	  $sroh=0.97;   ### 0.98, 1.012;
	  $soh=0.95;    ### 0.955, 0.965;
	  $ssulfur=0.985;
	  $snimi=0.99;  
	  $ss2=1.0;     ### 0.8;
	  
	  #ASP/GLU
	  $self->_sendCommand(sprintf("scalar $tag mult %1.4f sele (type od1 .or. type od2) .and. resnam ASP end",$soon));
	  $self->_sendCommand(sprintf("scalar $tag mult %1.4f sele (type oe1 .or. type oe2) .and. resnam GLU end",$soon));
	  
	  #TRP
	  $self->_sendCommand(sprintf("scalar $tag mult %1.4f sele type ne1 .and. resnam TRP end",$snind));
	
	  #ASN/GLN
	  $self->_sendCommand(sprintf("scalar $tag mult %1.4f sele type nd2 .and. resnam ASN end",$snh2));
	  $self->_sendCommand(sprintf("scalar $tag mult %1.4f sele type ne2 .and. resnam GLN end",$snh2));
	  
	  #ARG
	  $self->_sendCommand(sprintf("scalar $tag mult %1.4f sele (type nh1 .or. type nh2) .and. resnam ARG end",$snh2p));
	  
	  #LYS
	  $self->_sendCommand(sprintf("scalar $tag mult %1.4f sele type nz .and. resnam LYS end",$snh3p));
	  
	  #SER/THR
	  $self->_sendCommand(sprintf("scalar $tag mult %1.4f sele type og .and. resnam SER end",$sroh));
	  $self->_sendCommand(sprintf("scalar $tag mult %1.4f sele type og1 .and. resnam THR end",$sroh));
	  
	  #TYR
	  $self->_sendCommand(sprintf("scalar $tag mult %1.4f sele type oh .and. resnam TYR end",$soh));
	  
	  #CYS
	  $self->_sendCommand(sprintf("scalar $tag mult %1.4f sele type sg .and. resnam CYS end",$ssulfur));
	  
	  #MET
	  $self->_sendCommand(sprintf("scalar $tag mult %1.4f sele type sd .and. resnam MET end",$ss2));
	  
	  #HIS
	  $self->_sendCommand(sprintf("scalar $tag mult %1.4f sele (type nd1 .or. type ne2) .and. (resnam HSD .or. resnam HSE .or. resname HSP) end",$snimi));
	  
	  # $self->_sendCommand("scalar $tag set 0.3 sele prop $tag .lt. 0.3 end");
	}
      }
    }
  }
  $self->_sendCommand("scalar $tag mult $self->{par}->{sradfac}") 
     if (defined $self->{par}->{sradfac} && $self->{par}->{sradfac}>0.0);
#  $self->_sendCommand("scalar $tag show");
}


## method: setupGB([parameters])
## sets up the Generalized Born solvent approximation
## The non-bonded interactions are set to constant dielectric,
## epsilon = 1.0, and force switching.
## GB is then initialized and will be used in all
## subsequent energy evaluations. 

sub setupGB {
  my $self=shift;

  $self->_getpar(@_);

  if ($self->{par}->{cut}) {
    $self->_sendCommand("update atom CDIE eps $self->{par}->{epsilon} cutnb $self->{par}->{cutnb} ctofnb $self->{par}->{cutoff} ctonnb $self->{par}->{cuton} switch vswitch");
  } else {
    $self->_sendCommand("update atom CDIE eps $self->{par}->{epsilon} cutnb 999.0 ctofnb 990.0 ctonnb 980.0 switch vswitch");
  }    

  if ($self->{par}->{gb} eq "gbmvg" || $self->{par}->{gb} eq "gbmvgrid") {
    $self->_setupGBMVG();
  } elsif ($self->{par}->{gb} eq "gbmva") { 
    $self->_setupGBMVA();
  } elsif ($self->{par}->{gb} eq "bgb") {
    $self->_setupOrigGB();
  } elsif ($self->{par}->{gb} eq "gbsw") {
    $self->_setupGBSW();
  } elsif ($self->{par}->{gb}) {
    $self->{par}->{gb}="gbmva";
    $self->_setupGBMVA();
  }

  $self->{par}->{dynnose}=1 unless (defined $self->{par}->{dynnose});
  $self->{par}->{dyntrfrq}=200 if (!defined $self->{par}->{dyntrfrq});

  $self->{_haveGB}=1;
}

### setup Wonpil's GB ######

sub _setupGBSW {
  my $self=shift;

  $self->_sendCommand("scalar wmain = radius");
  if (defined $self->{par}->{scalerad} && $self->{par}->{scalerad}) {
    $self->_scaleRadii("wmain");
  }

  if ($self->{par}->{gbswms}) {
    $self->_sendCommand("GBSW sgamma $self->{par}->{gbswsgamma} dgp $self->{par}->{gbswdgp} molsurf");
  } else {
    $self->_sendCommand("GBSW sw $self->{par}->{gbswsw} nrad $self->{par}->{gbswnrad} rmax $self->{par}->{gbswrmax} nang $self->{par}->{gbswnang} sgamma $self->{par}->{gbswsgamma} dgp $self->{par}->{gbswdgp} rbuffer $self->{par}->{gbswrbuffer} tmemb $self->{par}->{gbswtmemb} msw $self->{par}->{gbswmsw}");

  }
  }
  

### setup original GB ######

sub _setupOrigGB {
  my $self=shift;

  $self->_sendCommand("scalar wmain = radius");

  if (defined $self->{par}->{scalerad} && $self->{par}->{scalerad}) {
    $self->_scaleRadii("wmain");
  } else {
    if ($self->{par}->{param} == 19) {
      $self->_sendCommand("scalar wmain mult $self->{par}->{gbcrad} select type c end");
    } else {
      $self->_sendCommand("scalar wmain set $self->{par}->{gbhrad} select ( prop radius .lt. 0.8 )  end");
    }
  }      

  my $lambda;
  if ($self->{par}->{param} == 19) {
    $lambda=(defined $self->{par}->{gblambda})?$self->{par}->{gblambda}:0.7295;
    $self->_sendCommand("gborn p1 0.4152 p2 0.2386 p3 1.7558 p4 10.5082 p5 1.1 lambda $lambda Epsilon $self->{par}->{gbeps} weight");
  } else {
# ($self->{par}->{param} =~/22|simdb/) {
    $lambda=(defined $self->{par}->{gblambda})?$self->{par}->{gblambda}:0.705;
    $self->_sendCommand("gborn p1 0.448 p2 0.173 p3 0.013 p4 9.015 p5 0.9 lambda $lambda Epsilon $self->{par}->{gbeps} weight");
#  } else {
#    die "GB not available for $self->{par}->{param} paramters";
  }
}

### setup GBMV grid based method ######

sub _setupGBMVG {
  my $self=shift;

  $self->{par}->{gbmvad}=-0.007998 if (!defined $self->{par}->{gbmvad});
  $self->{par}->{gbmvas}=0.9026 if (!defined $self->{par}->{gbmvas});
  $self->{par}->{gbmvaig}=99 if (!defined $self->{par}->{gbmvaig});

  my $cmd="GBMV GRID EPSILON $self->{par}->{gbeps} DN 0.2 watr 1.4 -\n";
  $cmd.="P6 $self->{par}->{gbmvap6} KAPPA 0.0 CONV -\n";
  $cmd.="SLOPE $self->{par}->{gbmvas} -\n" if (defined $self->{par}->{gbmvas});
  $cmd.="SHIFT $self->{par}->{gbmvad} -\n" if (defined $self->{par}->{gbmvad});
  $cmd.="ESHIFT $self->{par}->{gbmvade} -\n" if (defined $self->{par}->{gbmvade});
  $cmd.="CORR $self->{par}->{gbmvacorr} -\n" if (defined $self->{par}->{gbmvacorr} && $self->{par}->{gbmvacorr});
  $cmd.="A1 $self->{par}->{gbmva1} -\n" if (defined $self->{par}->{gbmva1});
  $cmd.="A2 $self->{par}->{gbmva2} -\n" if (defined $self->{par}->{gbmva2});
  $cmd.="A3 $self->{par}->{gbmva3} -\n" if (defined $self->{par}->{gbmva3});


  if ($self->{par}->{gbmvaig}==12) {
    $cmd.="WTYP 2 NPHI 12";
  } elsif ($self->{par}->{gbmvaig}==20) {
    $cmd.="WTYP 0";
  } elsif ($self->{par}->{gbmvaig}==26) {
    $cmd.="WTYP 2 NPHI 26";
  } elsif ($self->{par}->{gbmvaig}==32) {
    $cmd.="WTYP 1 NPHI 5";
  } elsif ($self->{par}->{gbmvaig}==38) {
    $cmd.="WTYP 2 NPHI 38";
  } elsif ($self->{par}->{gbmvaig}==99) {
    $cmd.="WTYP 0 NPHI 10";
  } elsif ($self->{par}->{gbmvaig}==95) {
    $cmd.="WTYP 1 NPHI 8";
  } elsif ($self->{par}->{gbmvaig}==96) {
    $cmd.="WTYP 1 NPHI 10";
  } elsif ($self->{par}->{gbmvaig}==97) {
    $cmd.="WTYP 1 NPHI 12";
  }


  if (defined $self->{par}->{scalerad} && $self->{par}->{scalerad}) {
    $self->_sendCommand("scalar wmain = radius");
    $self->_scaleRadii("wmain");
    $cmd.=" WEIGHT";
  }

  $self->_sendCommand($cmd);
}

### setup GBMV analytical method ######

sub _setupGBMVA {
  my $self=shift;

  if (defined $self->{par}->{gbmveps}) {
    if ($self->{par}->{gbmveps}>$self->{par}->{epsilon}) {
      $self->{par}->{gbmvacorr}=2;
      $self->{par}->{gbmvad}=-0.14 unless (defined $self->{par}->{gbmvad});
      $self->{par}->{gbmvdoff}=1 unless (defined $self->{par}->{gbmvdoff});
      $self->{par}->{gbmvdf}=-0.15 unless (defined $self->{par}->{gbmvdf});
      $self->{par}->{gbmvad}+=$self->{par}->{gbmvdf}/($self->{par}->{gbmveps}+$self->{par}->{gbmvdoff});
      $self->{par}->{gbmvade}=0.0 unless (defined $self->{par}->{gbmvade});
      $self->{par}->{gbmvas}=1.0 unless (defined $self->{par}->{gbmvas});
      $self->{par}->{gbmva1}=0.3255 unless (defined $self->{par}->{gbmva1});
      $self->{par}->{gbmva2}=0.0 unless (defined $self->{par}->{gbmva2});
      $self->{par}->{gbmva3}=1.085 unless (defined $self->{par}->{gbmva3});
      $self->{par}->{gbmva3}*=(3.0*$self->{par}->{gbmveps}/(3.0*$self->{par}->{gbmveps}+2.0*$self->{par}->{epsilon}));
      $self->{par}->{gbeps}=$self->{par}->{gbmveps};

#      $self->{par}->{gbmvad}=-0.145 unless (defined $self->{par}->{gbmvad});
#      $self->{par}->{gbmvdf}=-0.40 unless (defined $self->{par}->{gbmvdf});
#      $self->{par}->{gbmvdoff}=0 unless (defined $self->{par}->{gbmvdoff});
#      $self->{par}->{gbmvad}+=$self->{par}->{gbmvdf}/($self->{par}->{gbmveps}+$self->{par}->{gbmvdoff});
#      $self->{par}->{gbmvade}=0.0 unless (defined $self->{par}->{gbmvade});
#      $self->{par}->{gbmvas}=1.0 unless (defined $self->{par}->{gbmvas});
#      my $efac=(1.0-1.0/$self->{par}->{gbmveps});
#      $self->{par}->{gbmva1}=(1.0-1.0/sqrt(2.0)) unless (defined $self->{par}->{gbmva1});
#      $self->{par}->{gbmva2}=0.6169 unless (defined $self->{par}->{gbmva2});
#      $self->{par}->{gbmva3}=0.5754 unless (defined $self->{par}->{gbmva3});
#      $self->{par}->{gbmva1}/=$efac;
#      $self->{par}->{gbmva2}*=$efac;
#      $self->{par}->{gbmva3}*=$efac;
#      $self->{par}->{gbeps}=$self->{par}->{gbmveps};

#      if (defined $self->{par}->{gbmvepsu}) {
#	$self->{par}->{gbmvade}+=0.0;
#	$self->{par}->{gbmvad}+=(1/$self->{par}->{gbmvepsu}-1)*$self->{par}->{gbmvdu};
#      }
    } elsif ($self->{par}->{gbmveps}==1 && $self->{par}->{epsilon}>1) {
#      my $efac=(1.0/$self->{par}->{epsilon}-1.0);
#      $self->{par}->{gbmvacorr}=2;
#      $self->{par}->{gbmvad}=0.647993 unless (defined $self->{par}->{gbmvad});
#      $self->{par}->{gbmvad}=-0.4860 unless (defined $self->{par}->{gbmvad});
#      $self->{par}->{gbmvade}=-0.030415 unless (defined $self->{par}->{gbmvade});
#      $self->{par}->{gbmvade}*=$efac*$efac;
#      $self->{par}->{gbmvade}*=$efac;
#      $self->{par}->{gbmvas}=-1.0 unless (defined $self->{par}->{gbmvas});
#      $self->{par}->{gbmvas}*=$efac;
#      $self->{par}->{gbmva1}=1.1814 unless (defined $self->{par}->{gbmva1});
#      $self->{par}->{gbmva1}*=$efac;
#      $self->{par}->{gbmva2}=0.95538 unless (defined $self->{par}->{gbmva2});
#      $self->{par}->{gbmva2}*=$efac*$efac;
#      $self->{par}->{gbmva2}*=$efac;
#      $self->{par}->{gbmva3}=0.0 unless (defined $self->{par}->{gbmva3});
#      $self->{par}->{gbmva3}*=$efac*$efac;
#      $self->{par}->{gbmva3}*=$efac;

      $self->{par}->{gbmvacorr}=2;
      $self->{par}->{gbmvad}=-0.55 unless (defined $self->{par}->{gbmvad});
      $self->{par}->{gbmvdf}=0.20 unless (defined $self->{par}->{gbmvdf});
      $self->{par}->{gbmvad}+=$self->{par}->{gbmvdf}*$self->{par}->{epsilon};
      $self->{par}->{gbmvade}=0.0 unless (defined $self->{par}->{gbmvade});
      $self->{par}->{gbmvas}=1.0 unless (defined $self->{par}->{gbmvas});
      $self->{par}->{gbmva1}=1.2390 unless (defined $self->{par}->{gbmva1});
      $self->{par}->{gbmva2}=0.0 unless (defined $self->{par}->{gbmva2});
      $self->{par}->{gbmva3}=0.4592 unless (defined $self->{par}->{gbmva3});
      $self->{par}->{gbmva3}*=((1.0/(1.0+$self->{par}->{epsilon}))-1);
      $self->{par}->{gbeps}=1.0;
    } else {
      die "invalid GBMV dielectric solvent: $self->{par}->{gbmveps}, solute: $self->{par}->{epsilon}";
    }
  } else {


#    if ($self->{par}->{param} =~ /19/) {
#      $self->{par}->{gbmvap6}=11.58 if (!defined $self->{par}->{gbmvap6});
#      $self->{par}->{gbmvad}=-0.1765 if (!defined $self->{par}->{gbmvad});
#      $self->{par}->{gbmvade}=-1.2584 if (!defined $self->{par}->{gbmvade});
#    } elsif ($self->{par}->{param} =~ /22/ || $self->{par}->{param} eq "opls") {
#      $self->{par}->{gbmvas}=0.9085 if (!defined $self->{par}->{gbmvas});
#      $self->{par}->{gbmvad}=-0.052 if (!defined $self->{par}->{gbmvad});   
#      -0.052 for 0.0 PB, -0.072 for 0.1 PB, -0.082 for 0.15 PB, -0.102 for 0.25 PB
#    } else {
#      $self->{par}->{gbmvas}=0.9085 if (!defined $self->{par}->{gbmvas});
#      $self->{par}->{gbmvad}=-0.102 if (!defined $self->{par}->{gbmvad});   
#    }    

    $self->{par}->{gbmvas}=0.900 if (!defined $self->{par}->{gbmvas});
    $self->{par}->{gbmvad}=-0.100 if (!defined $self->{par}->{gbmvad});
  }


  $self->{par}->{gbmvaig}=32 if (!defined $self->{par}->{gbmvaig});


  my $cmd="GBMV TOL 1E-8 MEM 20 CUTA 20 DN 1.0 -\n";
  $cmd.="BUFR $self->{par}->{gbmvabuf} -\n" if (defined $self->{par}->{gbmvabuf});
  $cmd.="EPSILON $self->{par}->{gbeps} -\n" if (defined $self->{par}->{gbeps});
  $cmd.="BETA $self->{par}->{gbmvabeta} -\n" if (defined $self->{par}->{gbmvabeta});
  $cmd.="SHIFT $self->{par}->{gbmvad} -\n" if (defined $self->{par}->{gbmvad});
  $cmd.="ESHIFT $self->{par}->{gbmvade} -\n" if (defined $self->{par}->{gbmvade});
  $cmd.="SLOPE $self->{par}->{gbmvas} -\n" if (defined $self->{par}->{gbmvas});
  $cmd.="LAMBDA1 $self->{par}->{gbmval1} -\n" if (defined $self->{par}->{gbmval1});
  $cmd.="P1 $self->{par}->{gbmvap1} -\n" if (defined $self->{par}->{gbmvap1});
  $cmd.="P2 $self->{par}->{gbmvap2} -\n" if (defined $self->{par}->{gbmvap2});
  $cmd.="P3 $self->{par}->{gbmvap3} -\n" if (defined $self->{par}->{gbmvap3});
  $cmd.="P4 $self->{par}->{gbmvap4} -\n" if (defined $self->{par}->{gbmvap4});
  $cmd.="P5 $self->{par}->{gbmvap5} -\n" if (defined $self->{par}->{gbmvap5});
  $cmd.="P6 $self->{par}->{gbmvap6} -\n" if (defined $self->{par}->{gbmvap6});
  $cmd.="P7 $self->{par}->{gbmvap7} -\n" if (defined $self->{par}->{gbmvap7});
  $cmd.="ONX $self->{par}->{gbmvaonx} -\n" if (defined $self->{par}->{gbmvaonx});
  $cmd.="OFFX $self->{par}->{gbmvaoffx} -\n" if (defined $self->{par}->{gbmvaoffx});
  $cmd.="CUBIC $self->{par}->{gbmvacubic} -\n" if (defined $self->{par}->{gbmvacubic});
  $cmd.="CORR $self->{par}->{gbmvacorr} -\n" if (defined $self->{par}->{gbmvacorr} && $self->{par}->{gbmvacorr});
  $cmd.="A1 $self->{par}->{gbmva1} -\n" if (defined $self->{par}->{gbmva1});
  $cmd.="A2 $self->{par}->{gbmva2} -\n" if (defined $self->{par}->{gbmva2});
  $cmd.="A3 $self->{par}->{gbmva3} -\n" if (defined $self->{par}->{gbmva3});
  $cmd.="ARITH -\n" if (defined $self->{par}->{gbmvarith} && $self->{par}->{gbmvarith});
  $cmd.="ALFRQ $self->{par}->{gbmvafrq} -\n" if (defined $self->{par}->{gbmvafrq});
  $cmd.="DRFRQ $self->{par}->{gbmvdrfrq} -\n" if (defined $self->{par}->{gbmvdrfrq});
  $cmd.="LIMP -\n" if (defined $self->{par}->{gbmvafrq} && defined $self->{par}->{gbmvaimp} && $self->{par}->{gbmvaimp});
  $cmd.="EMP $self->{par}->{gbmvaemp} -\n" if (defined $self->{par}->{gbmvafrq} && defined $self->{par}->{gbmvaemp});
  $cmd.="SON $self->{par}->{gbmvason} SOFF $self->{par}->{gbmvasoff} -\n" if (defined $self->{par}->{gbmvason} && defined $self->{par}->{gbmvasoff});
  if (defined $self->{par}->{gbmvaextra}) {
    my $t=$self->{par}->{gbmvaextra};
    $cmd.="$t -\n";
  }

  if (defined $self->{par}->{gbmvmstl} && $self->{par}->{gbmvmstl}) {
    $cmd.="MODSTILL 1 MS0 $self->{par}->{gbmvms0} MS1 $self->{par}->{gbmvms1} MS2 $self->{par}->{gbmvms2} MS3 $self->{par}->{gbmvms3} -\n";
  }

  if ($self->{par}->{gbmvaig}==12) {
    $cmd.="WTYP 2 NPHI 12";
  } elsif ($self->{par}->{gbmvaig}==20) {
    $cmd.="WTYP 0";
  } elsif ($self->{par}->{gbmvaig}==26) {
    $cmd.="WTYP 2 NPHI 26";
  } elsif ($self->{par}->{gbmvaig}==32) {
    $cmd.="WTYP 1 NPHI 5";
  } elsif ($self->{par}->{gbmvaig}==38) {
    $cmd.="WTYP 2 NPHI 38";
  } elsif ($self->{par}->{gbmvaig}==93) {
    $cmd.="WTYP 1 NPHI 6";
  } elsif ($self->{par}->{gbmvaig}==94) {
    $cmd.="WTYP 1 NPHI 7";
  } elsif ($self->{par}->{gbmvaig}==95) {
    $cmd.="WTYP 1 NPHI 8";
  } elsif ($self->{par}->{gbmvaig}==96) {
    $cmd.="WTYP 1 NPHI 10";
  } elsif ($self->{par}->{gbmvaig}==97) {
    $cmd.="WTYP 1 NPHI 12";
  } elsif ($self->{par}->{gbmvaig}==99) {
    $cmd.="WTYP 0 NPHI 10";
  }

  $cmd.=" SA $self->{par}->{gbmvsa}" if (defined $self->{par}->{gbmvsa});
  $cmd.=" SB $self->{par}->{gbmvsb}" if (defined $self->{par}->{gbmvsb});

  if (defined $self->{par}->{scalerad} && $self->{par}->{scalerad}) {
    $self->_sendCommand("scalar wmain = radius");
    $self->_scaleRadii("wmain");
    $cmd.=" WEIGHT";
  }

  if (defined $self->{par}->{gbmvsagb}) {
    $self->_sendCommand("scalar sagb set 0.0 select all end");
    foreach my $saf ( split(/;/,$self->{par}->{gbmvsagb}) ) {
      my @tf=split(/:/,$saf);
      $self->_sendCommand("scalar sagb set $tf[1] select $tf[0] end");
    }
    $self->_sendCommand("scalar sagb show");
  }

  $self->_sendCommand($cmd);
}

## method: clearGB()
## clear and switch off Generalized Born solvent
## approximation 

sub clearGB {
  my $self=shift;
  if (defined $self->{par}->{gb}) {
    if ($self->{par}->{gb} eq "bgb") {
      $self->_sendCommand("gborn clear");
    } elsif ($self->{par}->{gb} eq "gbmv") {
      $self->_sendCommand("gbmv clear");
    }
  }
  $self->{_haveGB}=0;
}

## method: noeRestraints([parameters])
## sets up NOE restraints with CHARMM input from file
## in noe parameter

sub noeRestraints {
  my $self=shift;
  $self->_getpar(@_);

  my $noelist=();
  if (defined $self->{par}->{noerest} && -r $self->{par}->{noerest}) {
    my $inp=&GenUtil::getInputFile($self->{par}->{noerest});
    while(<$inp>) {
      push(@{$noelist},$_);
    }
    undef $inp;
  } elsif (defined $self->{par}->{xnoerest} && -r $self->{par}->{xnoerest}) {
    my $inp=&GenUtil::getInputFile($self->{par}->{xnoerest});
    while(<$inp>) {
      if (/assi[a-z]* *\((.*)\) *\((.*)\) *([0-9\.]+) +([0-9\.]+) +([0-9\.]+)/) {
	my $sel1=$1;
	my $sel2=$2;
	my $r=$3;
	my $dmin=$4;
	my $dmax=$5;
	
	$sel1=~s/and/.and./g;
	$sel1=~s/or/.or./g;
	$sel1=~s/name/type/g;
	$sel1=~s/\#/\*/g;

	$sel2=~s/and/.and./g;
	$sel2=~s/or/.or./g;
	$sel2=~s/name/type/g;
	$sel2=~s/\#/\*/g;

	push(@{$noelist},
	     sprintf("assign rexp %d kmin %f kmax %f fmax %f rmin %f rmax %f select %s end select %s end\n",
		     $self->{par}->{noerexp},$self->{par}->{noekmin},$self->{par}->{noekmax},
		     $self->{par}->{noefmax},$r-$dmin,$r+$dmax,$sel1,$sel2));
      }
    }
  } else {
    return;
  }

  $self->_sendCommand("bomlev -2");
  $self->_sendCommand("NOE\nRESET\nSCALE $self->{par}->{noescale}\nEND");
  for (my $in=0; $in<=$#{$noelist}; $in+=20) {
    my $cmd="";
    for (my $j=$in; $j<=$#{$noelist} && $j<$in+20; $j++) {
      $cmd.=$noelist->[$j];
    }

    $self->_sendCommand("NOE\n".$cmd."\nEND");
  }
  $self->_sendCommand("bomlev 0");
}


## method: getNOEAnalysis()
## obtain NOE distance restraint violations from CHARMM

sub getNOEAnalysis {
  my $self=shift;

  $self->_sendCommand("NOE\nPRINT ANAL\nEND");
  $self->_getCHARMMOutput("PERLDONE\n  \n");

  my $nrec;
  my $flag=0;
  my $nset=();
  foreach my $l ( split(/\n/,$self->{_lastOutput}) ) {
    if ($l=~/RESTRAINT: *([0-9]+) *([0-9]+) *([0-9]+)/) {
      $nrec={};
      $nrec->{inx}=$1;
      $nrec->{natom1}=$2;
      $nrec->{natom2}=$3;
    } elsif ($l=~/FIRST SET ATOM: +([0-9]+) +([A-Z].*) +([0-9]+) /) {
      $nrec->{res1}=$3;
    } elsif ($l=~/SECOND SET ATOM: +([0-9]+) +([A-Z].*) +([0-9]+) /) {
      $nrec->{res2}=$3;
      $flag=1;
    } elsif ($flag==1) {
      $l=~s/^[ \t]+//;
      my @f=split(/[ \t]+/,$l);
      $nrec->{rmin}=$f[0];
      $nrec->{rmax}=$f[2];
      $flag=2;
    } elsif ($flag==2) {
      $l=~s/^[ \t]+//;
      my @f=split(/[ \t]+/,$l);
      $nrec->{ractual}=$f[0];
      $nrec->{violation}=$f[1];
      $flag=0;
      push(@{$nset},$nrec);
    }
  }
  return $nset;
}

## method: minimizeSD([parameters])
## runs a steepest descent minimization.
## Parameters are the non-bonded list update 
## frequency <mark>updnbsd</mark>, the
## number of minimization steps <mark>itersd</mark>
## and the initial step size <mark>stepsd</mark>.

sub minimizeSD { 
  my $self=shift;
  $self->_getpar(@_);
  $self->_sendCommand("mini sd nprint $self->{par}->{minoutfrq} inbfreq $self->{par}->{updnbsd} step $self->{par}->{sdstepsz} nstep $self->{par}->{sdsteps}");
}

## method: minimize([parameters])
## runs minimization protocol
## Parameters are the non-bonded list update 
## frequency <mark>minupdnb</mark>, the
## number of minimization steps <mark>miniter</mark>,
## the initial step size <mark>minstep</mark> and
## an energy tolerance <mark>minetol</mark> to limit
## the minimization run if the tolerance is reached
## in less than the requested number of runs.
## The minimization algorithm can be changed by
## setting <mark>minmode</mark>.

sub minimize {
  my $self=shift;
  $self->_getpar(@_);
  $self->_sendCommand("mini $self->{par}->{minmode} nprint $self->{par}->{minoutfrq} inbfreq $self->{par}->{minupdnb} step $self->{par}->{minstepsz} nstep $self->{par}->{minsteps} tolenr $self->{par}->{minetol}");
}

sub clearShake {
  my $self=shift;
  
  $self->_getpar(@_);
  $self->_sendCommand("shake off");
  $self->{_haveSHAKE}=0;
}

## method: shake([parameters])
## turns on SHAKE to fix bond distances. 
## A fast mode is available if the parameter <mark>shakefast</mark>
## is set. All bonds, including hydrogens, are fixed if <mark>shakemode</mark>
## is set to <mark>hyd</mark>, otherwise bonds involving hydrogens are
## excluded. The SHAKE tolerance is determined by <mark>shaketol</mark>.

sub shake {
  my $self=shift;

  $self->_getpar(@_);

  return if (!$self->{par}->{shake});
  
  my $cmd="shake";
  $cmd.=" fast" if ($self->{par}->{shakefast});
  $cmd.=(lc $self->{par}->{shakemode} eq "hyd")?" bonh":" bond";
  $cmd.=" tol $self->{par}->{shaketol}";
  $cmd.=" param";
  $self->_sendCommand($cmd);

  $self->{_haveSHAKE}=1;
}


## method: periodicBoundaries([parameters]) 
## turns on periodic boundaries.

sub periodicBoundaries {
  my $self=shift;
  $self->_getpar(@_);

  if (!defined $self->{par}->{boxsize}) {
    my $cstat=$self->coorStats();
    my $maxx=($cstat->{xmax}-$cstat->{xmin});
    my $maxy=($cstat->{ymax}-$cstat->{ymin});
    my $maxz=($cstat->{zmax}-$cstat->{zmin});
    if ($maxx>$maxy && $maxx>$maxz) {
      $self->{par}->{boxsize}=$maxx-0.5;
    } elsif ($maxy>$maxx && $maxy>$maxz) {
      $self->{par}->{boxsize}=$maxy-0.5;
    } else {
      $self->{par}->{boxsize}=$maxz-0.5;
    }
    $self->{par}->{boxsize}/=sqrt(4.0/3.0) 
      if ($self->{par}->{boxshape} =~ /octa/);
  }

  if ($self->{par}->{boxshape} =~ /octa/) {
    $self->_sendCommand("crystal defined octa $self->{par}->{boxsize} $self->{par}->{boxsize} $self->{par}->{boxsize} 109.4712206344907 109.4712206344907 109.4712206344907");
  } elsif ($self->{par}->{boxshape} eq "cubic") {
    $self->_sendCommand("crystal defined cubic $self->{par}->{boxsize} $self->{par}->{boxsize} $self->{par}->{boxsize} 90.0 90.0 90.0");
  } else {
    die "unknown box shape $self->{par}->{boxshape}";
  }

  $self->_sendCommand("crystal build cutoff $self->{par}->{cutnb} noper 0");
  $self->_sendCommand("image byseg xcen 0.0 ycen 0.0 zcen 0.0 select .not. resname tip3 .and. .not. resname cla .and. .not. resname sod end");
  $self->_sendCommand("image byres xcen 0.0 ycen 0.0 zcen 0.0 select resname tip3 .or. resname cla .or. resname sod end");
  my $cutim=$self->{par}->{cutnb}+2.0;
  $self->_sendCommand("update $self->{par}->{nblisttype} cutim $cutim imgfrq -1");
}


## method: ewald([parameters])
## turns on Ewald summation (PME) for non-bonded electrostatics. 
## The parameters <mark>pmekappa</mark> is used to set
## the kappa value.

sub ewald {
  my $self=shift;
  $self->_getpar(@_);

  my $qcor=(defined $self->{par}->{qcor})?"qcor $self->{par}->{qcor}":"";
  $self->_sendCommand("update ewald kappa $self->{par}->{pmekappa} pmewald order 4 fftx $self->{par}->{npme} ffty $self->{par}->{npme} fftz $self->{par}->{npme} $qcor");
  $self->{par}->{dyntrfrq}=200 if ((! defined $self->{par}->{dyntrfrq}) || ($self->{par}->{dyntrfrq} <= 0));
}

## method: runDynamics(restin,restout,trajout,enerout[,parameters])
## runs a molecular dynamics simulation. The first four arguments
## are file names for restart input and output, trajectory output and
## CHARMM energy output files. <mark>undef</mark> can be used instead
## of a file name if the corresponding file should not be created.
## If <mark>undef</mark> is given for the restart input file, the
## simulation is started with a random velocity assignment.
## A large number of parameters are available to control the simulation
## run. The most important values are <mark>dynsteps</mark> (number of
## dynamics steps), <mark>dyntemp</mark> (temperature), <mark>dynens</mark>
## (ensemble: "NVT" or "NPT"), <mark>dynoutfrq</mark> (frequency of
## energy output), <mark>dyneqfrq</mark> (frequency of velocity reassignment
## to maintain constant temperature in NVT), and <mark>dynpress</mark> (pressure
## for NPT ensemble).

sub runDynamics {
  my $self=shift;
  my $restartin=shift;
  my $restartout=shift;
  my $trajout=shift;
  my $enerout=shift;

  $self->_getpar(@_);
  
  my $opt=(defined $restartin && -r $restartin)?
    $self->_getRestartDynOptions():$self->_getStartDynOptions();

  my $cmd;
  if (!$self->{par}->{lang}) {
    if ($self->{par}->{sgmd}) {
      $cmd="dyna leap sgmd ydrag $self->{par}->{sgmddrag} nvufrq $self->{par}->{sgmdfrq} -\n";
      $cmd.="yavg $self->{par}->{sgmdavg} vfb0 $self->{par}->{sgmdvfb} -\n";
      $cmd.="tconst tcoup -1.0 tref $self->{par}->{dyntemp} echeck 100 -\n";
      $cmd.=$opt;
    } else {
      if ($self->{par}->{dynber}) {
	$cmd="dyna leap cpt tcons tcoupling $self->{par}->{dynbertc} tref $self->{par}->{dyntemp} ";
	$cmd.=$opt;
      } elsif (defined $self->{par}->{dynnose} && $self->{par}->{dynnose}) {
	if (uc $self->{par}->{dynens} eq "NPT") {
	  $cmd="dyna leap cpt ";
	  $cmd.="hoover reft $self->{par}->{dyntemp} tmass $self->{par}->{dynnoseq} ";
	  $cmd.=$opt;
	  $cmd.=$self->_getNPTDynOptions();
	} else {
	  $cmd="dyna vver nose rstn qref $self->{par}->{dynnoseq} tref $self->{par}->{dyntemp} ncyc $self->{par}->{dynnosen} ";
	  $cmd.=$opt;
	}
      } else { 
	$cmd="dyna leap ";
	$cmd.="cpt " 
	  if (uc $self->{par}->{dynens} eq "NPT");
	$cmd.=$opt;
	$cmd.=$self->_getNPTDynOptions() 
	  if (uc $self->{par}->{dynens} eq "NPT");
      }
    }
  } else {
    if (defined $self->{par}->{langsel}) {
      $self->_sendCommand("scalar fbeta set $self->{par}->{langfbeta} sele $self->{par}->{langsel} end");
    } else {
      $self->_sendCommand("scalar fbeta set $self->{par}->{langfbeta} sele .not. type H* end");
    }
    $cmd="dyna leap lang $opt";
    $cmd.="tbath $self->{par}->{dyntemp} rbuf $self->{par}->{langrbuf} ilbfreq $self->{par}->{langupd} -\n";
  }

  $cmd.="echeck $self->{par}->{echeck} -\n";

  if (defined $self->{par}->{dyncons} && $self->{par}->{dyncons}) {
    $cmd.="kdyn ";
    $cmd.=sprintf("kmas %f ",$self->{par}->{dynkmas}) 
       if (defined $self->{par}->{dynkmas});
    $cmd.=sprintf("fkc %f ",$self->{par}->{dynkfkc}) 
       if (defined $self->{par}->{dynkfkc});
    $cmd.=sprintf("dkc %f ",$self->{par}->{dynkdkc}) 
       if (defined $self->{par}->{dynkdkc});
    $cmd.=" -\n";
  }	

  $cmd.=$self->_getDynIO($restartin,$restartout,$trajout,$enerout);
  $self->_sendCommand($cmd);
  $self->{mdener}=$self->_processMDEneOutput();
  $self->_closeDynIO($restartin,$restartout,$trajout,$enerout);
  if (defined $self->{par}->{dyncons} && $self->{par}->{dyncons}) {
    $self->_sendCommand("scalar cons show");
  }
  $self->{_loadedref}=0;
}

## method: setupPert(select)
## setup PERT

sub setupPert {
  my $self=shift;
  my $sel=shift;
  my $pertfile=shift;

  $self->_sendCommand("cons fix select $sel end");
  $self->_sendCommand("PERT select $sel end");
  $self->_sendCommand("scalar charge set 0.0 select $sel end");

  my $out=&GenUtil::getOutputFile($pertfile);
  
  printf $out "* PERT input nlambda: %d, equi: %d, prod: %d\n",
    $self->{par}->{nlambda},$self->{par}->{pequi},$self->{par}->{pprod};
  printf $out "*\n";
  
  my $delta=1.0/$self->{par}->{nlambda};
  
  my $sum=0;
  for (my $i=0; $i<=$self->{par}->{nlambda}; $i++) {
    printf $out " LSTART %5.2f LAMBDA %5.2f LSTOP %5.2f PSTART %8d PSTOP %8d PWIND\n",
      ($i>0)?$delta*($i-0.5):0.0,$delta*$i,($i<$self->{par}->{nlambda})?$delta*($i+0.5):1.0,
      $sum+$self->{par}->{pequi},$sum+$self->{par}->{pequi}+$self->{par}->{pprod};
    $sum+=$self->{par}->{pequi}+$self->{par}->{pprod};
  }

  $self->{par}->{dynsteps}=$sum;

  printf $out " END\n";
  undef $out;
}

## method: runPertWHAM(filename) 

sub runPertWHAM {
  my $self=shift;
  my $fname=shift;
  
  $self->_sendCommand("open read card unit 54 name $fname");
  $self->_sendCommand("wham maxwindow 100 maxtime 100000 unit 54 tol $self->{par}->{whamtol} nstep $self->{par}->{whamsteps}");
  $self->_getCHARMMOutput("PERLDONE\n  \n")
    if ($self->{_lastOutput} eq "");

  my $final=0.0;
  foreach my $l ( split(/\n/,$self->{_lastOutput}) ) {
    if ($l=~/^.+Number_of_data.+F\(\) +([0-9\.]+)/) {
      $final=$1;
    }
  }

  $self->_sendCommand("close unit 54");
  return $final;
}

## method: runChargeFEP(trajout,wham,pertfile[,parameters])
## runs PERT to obtain charging free energy

sub runChargeFEP {
  my $self=shift;
  my $trajout=shift;
  my $wham=shift;
  my $pertfile=shift;

  $self->_getpar(@_);
  
  my $opt=$self->_getStartDynOptions();

  my $cmd;
  if (!$self->{par}->{lang}) {
    if ($self->{par}->{dynber}) {
      $cmd="dyna leap cpt tcons tcoupling $self->{par}->{dynbertc} tref $self->{par}->{dyntemp} ";
      $cmd.=$opt;
    } elsif (defined $self->{par}->{dynnose} && $self->{par}->{dynnose}) {
      if (uc $self->{par}->{dynens} eq "NPT") {
	$cmd="dyna leap cpt ";
	$cmd.="hoover reft $self->{par}->{dyntemp} tmass $self->{par}->{dynnoseq} ";
	$cmd.=$opt;
	$cmd.=$self->_getNPTDynOptions();
      } else {
	$cmd="dyna vver nose rstn qref $self->{par}->{dynnoseq} tref $self->{par}->{dyntemp} ncyc $self->{par}->{dynnosen} ";
	$cmd.=$opt;
      }
    } else { 
      $cmd="dyna leap ";
      $cmd.="cpt " 
	if (uc $self->{par}->{dynens} eq "NPT");
      $cmd.=$opt;
      $cmd.=$self->_getNPTDynOptions() 
	if (uc $self->{par}->{dynens} eq "NPT");
    }
  } else {
    if (defined $self->{par}->{langsel}) {
      $self->_sendCommand("scalar fbeta set $self->{par}->{langfbeta} sele $self->{par}->{langsel} end");
    } else {
      $self->_sendCommand("scalar fbeta set $self->{par}->{langfbeta} sele .not. type H* end");
    }
    $cmd="dyna leap lang $opt";
    $cmd.="tbath $self->{par}->{dyntemp} rbuf $self->{par}->{langrbuf} ilbfreq $self->{par}->{langupd} -\n";
  }

  $cmd.="echeck $self->{par}->{echeck} -\n";

  $cmd.=$self->_getDynIO(undef,undef,$trajout,undef);
  $self->_sendCommand("open read card unit 88 name $pertfile");
  $cmd.=" punit 88";
  if (defined $wham) {
    $self->_sendCommand("open write card unit 89 name $wham");
    $cmd.=" wham 89";
  }
  $self->_sendCommand($cmd);
  $self->{pertdata}=$self->_processPertOutput();
  $self->_closeDynIO(undef,undef,$trajout,undef);
  $self->_sendCommand("close unit 88");
  $self->_sendCommand("close unit 89") if (defined $wham);
}

### _closeDynIO ######

sub _closeDynIO {
  my $self=shift;
  my $restin=shift;
  my $restout=shift;
  my $trajout=shift;
  my $enerout=shift;

  $self->_sendCommand("close unit 11") 
    if (defined $restin);
  $self->_sendCommand("close unit 14") 
    if (defined $enerout);
  $self->_getCHARMMOutput("PERLDONE\n  \n")
    if (defined $restin || defined $enerout);

}

### _getDynIO ######

sub _getDynIO {
  my $self=shift;
  my $restin=shift;
  my $restout=shift;
  my $trajout=shift;
  my $enerout=shift;

  my $opt="";

  if (defined $restin && $restin ne "" && -r $restin) {
    $self->_sendCommand("open unit 11 read form name \"$restin\"");
    $opt.="iunrea 11";
  } else {
    $opt.="iunrea -1";
  }

  if (defined $restout && $restout ne "") {
    $self->_sendCommand("open unit 12 write form name \"$restout\"");
    $opt.=" iunwri 12";
  } else {
    $opt.=" iunwri -1";
  }
  
  if (defined $trajout && $trajout ne "") {
    $self->_sendCommand("open unit 13 write unform name \"$trajout\"");
    $opt.=" iuncrd 13";
  } else {
    $opt.=" iuncrd -1";
  }
   
  if (defined $enerout && $enerout ne "") {
    $self->_sendCommand("open unit 14 write form name \"$enerout\"");
    $opt.=" kunit 14";
  } else {
    $opt.=" kunit -1";
  }
  return $opt;
}

### _getStartDynOptions ######

sub _getStartDynOptions {
  my $self=shift;

  my $opt="start ".$self->_getGenDynOptions();
 
  $self->{par}->{dynseed}=2*int(rand 10000000)+1
    if (!defined $self->{par}->{dynseed});

  if (defined $self->{par}->{dynitemp} && 
      $self->{par}->{dynitemp} != $self->{par}->{dyntemp}) {
    if ($self->{par}->{dynhtfrq}>0) {
      $opt.="tstru $self->{par}->{dynitemp} firstt $self->{par}->{dynitemp} teminc $self->{par}->{dyndeltat} ihtfrq $self->{par}->{dynhtfrq}";
    } else {
      $opt.="tstru $self->{par}->{dynitemp} firstt $self->{par}->{dyntemp}";
    }
  } else {
    $opt.="tstru $self->{par}->{dyntemp} firstt $self->{par}->{dyntemp}";
  }
  $opt.=" iseed $self->{par}->{dynseed} -\n";

  $self->{par}->{dynseed}=undef;

  return $opt;
}

### _getRestartDynOptions ######

sub _getRestartDynOptions {
  my $self=shift;

  my $opt="rest ".$self->_getGenDynOptions();

  if (defined $self->{par}->{dynscale}) {
    $opt.="iscale 1 scale $self->{par}->{dynscale} -\n";
  }
 
  return $opt;
}

### _getGenDynOptions ######

sub _getGenDynOptions {
  my $self=shift;
  
  my $opt="timestep $self->{par}->{dyntstep} nstep $self->{par}->{dynsteps} finalt $self->{par}->{dyntemp} -\n";
  $opt.="inbfreq $self->{par}->{dynupdnb} imgfreq $self->{par}->{dynupdimg} $self->{par}->{nblisttype} -\n";
  $opt.="iasors 1 iasvel 1 -\n"; 
  $opt.="ntrfrq $self->{par}->{dyntrfrq} -\n" if (defined $self->{par}->{dyntrfrq});

  if ($self->{par}->{dyneqfrq}>0) {
    $opt.="ichecw 1 ieqfrq $self->{par}->{dyneqfrq} twindl -$self->{par}->{dyntwin} twindh +$self->{par}->{dyntwin} -\n";
  } else {
    $opt.="ichecw 0 -\n";
  }

  if (defined $self->{par}->{dynoutfrq}) {
    $opt.="nprint $self->{par}->{dynoutfrq} iprfrq $self->{par}->{dynsteps} nsavc $self->{par}->{dynoutfrq} nsavv 0 iunvel -1 -\n";
  } else {
    $opt.="nprint $self->{par}->{dynsteps} iprfrq $self->{par}->{dynsteps} nsavc $self->{par}->{dynsteps} nsavv 0 iunvel -1 -\n";
  }

  return $opt;
}

### _getNPTDynOptions ######

sub _getNPTDynOptions {
  my $self=shift;
  my $opt="pconstant pmass 100 pref $self->{par}->{dynpress} pgamma 10 -\n";
  return $opt;
}

## method: setupRestraints(fscale,conslist)
## sets up harmonic restraints from a list of 
## <mark>cons</mark> data structures (see <mark>harmonicRestraint</mark>).
## The first argument determines the scaling factor for
## all force constants.

sub setupRestraints {
  my $self=shift;
  my $fscale=shift;
  my $cons=shift;

  die "need to setup protein structure first"
    if (!defined $self->{molecule});

  my $refloaded;
  foreach my $n (@{$cons}) {
    
    if (!$self->{_keepRestraint}) {
      if (defined $n->{exclmode} && $n->{exclmode}) {
	$self->{molecule}->setValidResidues($n->{list},1);
	$n->{list}=&GenUtil::gradForceList($self->{molecule}->listFromValid(),$n->{force});
	$n->{exclmode}=0;
      }
    }

    if ($n->{type} eq "self") {
      $self->harmonicRestraint($n,$fscale,$self->{_keepRestraint});
    } elsif ($n->{type} eq "ref" || $n->{type} eq "rmsd") {
      if (!$self->{_keepRestraint}) {
	die "need reference file for restraint"
	  if (!defined $n->{reffile});

	if (!defined $refloaded) {
	  $self->loadReference($n->{reffile});
	  $refloaded=$n->{reffile};
	} else {
	  die "can have only one reference file"
	    if ($n->{reffile} ne $refloaded);
	}
      }
      if ($n->{type} eq "ref") {
        $self->harmonicRestraint($n,$fscale,$self->{_keepRestraint});
      } elsif ($n->{type} eq "rmsd") {
        $self->rmsdRestraint($n,$fscale);
      }
    } elsif ($n->{type} eq "hmcm") {
      die "need reference file for HMCM restraint"
	if (!defined $n->{reffile});
      $self->hmcmRestraint($n,$fscale);
    }
  }

  $self->{_keepRestraint}=1;
}

## method: harmonicRestraint(cons[,fscale[,keep]])
## sets up harmonic restraints according to the <mark>cons</mark>
## argument. It is expected to have the following data structure:
## <mark>cons -> {type, list[] -> {from, to, force}, sel}</mark>.
## The <mark>type</mark> field is either <mark>self</mark> 
## (current coordinates) 
## or <mark>ref</mark> (alternate coordinates in COMP). 
## <mark>list</mark> contains
## a list of protein segments that are supposed to be restrained.
## <mark>sel</mark> identifies which atoms are restrained for each
## residue in the selection list (<mark>ca</mark>, <mark>cb</mark>,
## <mark>cab</mark>, <mark>heavy</mark>, or <mark>all</mark>).
## An additional argument may be given to scale all restraint forces by
## a factor. If the last argument is set to 1 previously set 
## reference coordinates will be used and only force constants are updated.

sub harmonicRestraint {
  my $self=shift;
  my $cons=shift;
  my $fscale=shift;
  my $keep=shift;

  $fscale=1.0 unless (defined $fscale);

  my $compstr=($cons->{type} eq "ref" && !$keep)?"COMP":"";
  my $keepstr=(defined $keep && $keep)?"KEEP":"";
  
  my $atomsel=&_getSel($cons->{sel});

  my $fh={};
  foreach my $c ( @{$cons->{list}} ) {
    my $force=$c->{force}*$fscale;

    my $ch=$self->{molecule}->getChain($c->{chain});
    die "invalid restraint list: cannot find chain $c->{chain}"
      if (!defined $ch);

    my $rfirst=$self->{molecule}->firstResNum($ch);
    my $rlast=$self->{molecule}->lastResNum($ch);

    if ($c->{to}>=$rfirst && $c->{from}<=$rlast) {
      $c->{from}=$rfirst if ($c->{from}<$rfirst);
      $c->{to}=$rlast    if ($c->{to}>$rlast);

      my $rfrom=$self->{molecule}->getResidueInChain($c->{from},$ch);
      die "invalid restraint list: cannot find first residue $c->{chain}$c->{from} in molecule"
	if (!defined $rfrom);
    
      my $r=$ch->{res};

      my $rfrominx=$ch->{resinx}->{$c->{from}};
      my $rtoinx=$ch->{resinx}->{$c->{to}};

      die "invalid restraint list: cannot find last residue $c->{chain}$c->{to} in molecule"
	if (!defined $rtoinx);
    
      my $sel;
      my $first=$rfrominx;
      my $segid=$rfrom->{seg};
      for (my $ir=$rfrominx+1; $ir<=$rtoinx; $ir++) {
	if ($r->[$ir]->{seg} ne $segid) {
	  my $n1=$r->[$first]->{num};
	  my $n2=$r->[$ir-1]->{num};
	  $sel="( resid ".(($first!=$ir-1)?"$n1:$n2":"$n1")." .and. segid ".$segid." )";
	  push(@{$fh->{$force}},$sel);
	  $first=$ir;
	  $segid=$r->[$ir]->{seg};
	}
      }

      my $n1=$r->[$first]->{num};
      $sel="( resid ".(($first!=$rtoinx)?"$n1:$c->{to}":"$n1")." .and. segid ".$segid." )";
      push(@{$fh->{$force}},$sel);
    }
  }

  foreach my $f ( keys %{$fh} ) {
    foreach my $ff ( @{$fh->{$f}}) {
      my $cmd="cons harm force $f mass select ( $ff ) .and. ( $atomsel ) end $compstr $keepstr";
#    my $cmd="cons harm force $f mass select ( -\n";
#    $cmd.=join(" .or. -\n",@{$fh->{$f}});
#    $cmd.=" ) .and. ( $atomsel ) end $compstr $keepstr";
    $self->_sendCommand($cmd);
  }
  }
}  

sub rmsdRestraint {
  my $self=shift;
  my $cons=shift;
  my $fscale=shift;

  $fscale=1.0 unless (defined $fscale);

  $self->_sendCommand("rgyr reset")
    if (defined $self->{_havermsd} && $self->{_havermsd});

  my $oriestr=(defined $cons->{orient} && $cons->{orient})?"ORIE":"";
  my $refe=(defined $cons->{refe} && $cons->{refe}>0)?$cons->{refe}:0.0;

  my $atomsel=&_getSel($cons->{sel});

  my $fh={};
  foreach my $c ( @{$cons->{list}} ) {
    my $force=$c->{force}*$fscale;

    my $ch=$self->{molecule}->getChain($c->{chain});
    die "invalid restraint list: cannot find chain $c->{chain}"
      if (!defined $ch);

    my $rfirst=$self->{molecule}->firstResNum($ch);
    my $rlast=$self->{molecule}->lastResNum($ch);

    if ($c->{to}>=$rfirst && $c->{from}<=$rlast) {
      $c->{from}=$rfirst if ($c->{from}<$rfirst);
      $c->{to}=$rlast    if ($c->{to}>$rlast);

      my $rfrom=$self->{molecule}->getResidueInChain($c->{from},$ch);
      die "invalid restraint list: cannot find first residue $c->{chain}$c->{from} in molecule"
	if (!defined $rfrom);
    
      my $r=$ch->{res};

      my $rfrominx=$ch->{resinx}->{$c->{from}};
      my $rtoinx=$ch->{resinx}->{$c->{to}};

      die "invalid restraint list: cannot find last residue $c->{chain}$c->{to} in molecule"
	if (!defined $rtoinx);
    
      my $sel;
      my $first=$rfrominx;
      my $segid=$rfrom->{seg};
      for (my $ir=$rfrominx+1; $ir<=$rtoinx; $ir++) {
	if ($r->[$ir]->{seg} ne $segid) {
	  my $n1=$r->[$first]->{num};
	  my $n2=$r->[$ir-1]->{num};
	  $sel="( resid ".(($first!=$ir-1)?"$n1:$n2":"$n1")." .and. segid ".$segid." )";
	  push(@{$fh->{$force}},$sel);
	  $first=$ir;
	  $segid=$r->[$ir]->{seg};
	}
      }

      my $n1=$r->[$first]->{num};
      $sel="( resid ".(($first!=$rtoinx)?"$n1:$c->{to}":"$n1")." .and. segid ".$segid." )";
      push(@{$fh->{$force}},$sel);
    }
  }

  my $pres="(prop xcomp .gt. 0.0001 .or. prop xcomp .lt. -0.0001) .and. -\n(prop ycomp .gt. 0.0001 .or. prop ycomp .lt. -0.0001) .and. -\n(prop zcomp .gt. 0.0001 .or. prop zcomp .lt. -0.0001)";

  foreach my $f ( keys %{$fh} ) {
    foreach my $ff ( @{$fh->{$f}}) {
      my $cmd="rgyr force $f refe $refe rmsd comp $oriestr select ( $ff ) .and. ( $atomsel ) -\n .and. ( $pres ) end";
      $self->_sendCommand($cmd);
    }
  }
  $self->{_havermsd}=1;
}  

## method: simpleRestraint(force,select,keep,comp)
## requests harmonic restraints with a simpler interface as
## <mark>harmonicRestraint</mark>. The arguments are a
## force constant, a CHARMM selection expression and
## flags whether to keep previous coordinates and whether
## to use coordinates from the alternate coordinate set

sub simpleRestraint {
  my $self=shift;
  my $force=shift;
  my $select=shift;
  my $keep=shift;
  my $comp=shift;

  my $compstr=(defined $comp && $comp && !$keep)?"COMP":"";
  my $keepstr=(defined $keep && $keep)?"KEEP":"";
  
  my $cmd="cons harm force $force mass select $select end $compstr $keepstr\n";
  $self->_sendCommand($cmd);
}  

## method: hmcmRestraint(cons[,fscale])
## sets up harmonic center of mass restraints according to
## side chain centers from a SICHO chain file given in the
## <mark>reffile</mark> field of the <mark>cons</mark>
## argument. 

sub hmcmRestraint {
  my $self=shift;
  my $cons=shift;
  my $fscale=shift;

  $fscale=1.0 unless (defined $fscale);

  $self->_sendCommand("cons hmcm clear");
  
  my $ref=SICHO::new();
  $ref->readChain($cons->{reffile});

  foreach my $c (@{$cons->{list}}) {
    my $force=(defined $c->{force}?$c->{force}:1.0);
    $force*=$fscale;

    my $ch=$self->{molecule}->getChain($c->{chain});
    die "invalid restraint list: cannot find chain $c->{chain}"
      if (!defined $ch);

    my $rfirst=$self->{molecule}->firstResNum($ch);
    my $rlast=$self->{molecule}->lastResNum($ch);

    if ($c->{to}>=$rfirst && $c->{from}<=$rlast) {
      $c->{from}=$rfirst if ($c->{from}<$rfirst);
      $c->{to}=$rlast    if ($c->{to}>$rlast);
    }

    for (my $i=$c->{from}; $i<=$c->{to}; $i++) {
      my $chain=$ref->{sidechain}->[$i];
      my $x=($chain->{xcoor}-$ref->{offset}->{xcoor})*$ref->{resolution};
      my $y=($chain->{ycoor}-$ref->{offset}->{ycoor})*$ref->{resolution};
      my $z=($chain->{zcoor}-$ref->{offset}->{zcoor})*$ref->{resolution};
      $self->_sendCommand("cons hmcm force $force -\n refx $x refy $y refz $z -\n select resid $i .and. .not. -\n ( type h* .or. type n .or. type c .or. -\n type o .or. type oct* .or. type ot* ) end");
    }
  }
  $self->{_havehmcm}=1;
}

## method: clearRestraints()
## clears all harmonic restraints

sub clearRestraints {
  my $self=shift;
  $self->_sendCommand("cons harm clear");
  $self->_sendCommand("cons hmcm clear")
    if (defined $self->{_havehmcm} && $self->{_havehmcm});
  $self->_sendCommand("rgyr reset")
    if (defined $self->{_havermsd} && $self->{_havermsd});
  $self->{_keepRestraint}=0;
  $self->{_havehmcm}=0;
  $self->{_havermsd}=0;
}


## method: applyBias(parameters)
## imposes a harmonic restraint of the type described
## in the <mark>biastype</mark> data structure.
## This data structure must contain the key "type",
## while other type-dependent keys may also be present.

sub applyBias {
  my $self=shift;
  my %par=@_;

  die "need to setup protein structure first"
    if (!defined $self->{molecule});

die "Bias type is not defined\n" if (!defined $par{type});

  if ($par{type} eq "rg") {
    die "rg bias must be cleared before being reused"
      if ((defined $self->{_biasstatus}->{rg}) && 
          ($self->{_biasstatus}->{rg} eq "on"));
    my $cmd="rgyr";
    $cmd .= " force " . $par{force};
    $cmd .= " refe " . $par{target};
    $cmd .= " sele type ca end";
    $self->_sendCommand($cmd);
    $self->{_biasstatus}->{rg}="on";
  } elsif ($par{type} eq "rmsd") {
    die "rg bias must be cleared before being reused"
      if ((defined $self->{_biasstatus}->{rg}) &&
          ($self->{_biasstatus}->{rg} eq "on"));
    if (!$self->{_loadedref}) {
      $self->loadReference($par{pdb});
      $self->{_loadedref}=1;
    }
    my $cmd="rgyr";
    $cmd .= " force " . $par{force};
    $cmd .= " refe " . $par{target};
    $cmd .= " rmsd orie comp ";
    $cmd .= " sele type ca .or. type cb end";
    $self->_sendCommand($cmd);
    $self->{_biasstatus}->{rg}="on";
  } elsif ($par{type} eq "rho") {
    die "rho bias must be cleared before being reused"
      if ((defined $self->{_biasstatus}->{rho}) && 
          ($self->{_biasstatus}->{rho} eq "on"));
    my $cmd="dmco";
    my $def;
    if ((defined $par{sidechains}) && ($par{sidechains} eq "on")) {
      $def="define bb sele (type ca .or. type c .or. type n .or. type o ) end";
      $self->_sendCommand($def);
      $def="define sd sele .not. resname tip3 .and. .not. (bb .or. hydrogen) end";
      $self->_sendCommand($def);
    } else {
      $def="define aa sele .not. resname tip3 .and. .not. hydrogen end";
      $self->_sendCommand($def);
    }
    $cmd .= " force " . $par{force};
    $cmd .= " refe " . $par{target};
    $par{qcut} = 6.5 if (! defined $par{qcut});
    $cmd .= " cutoff " . $par{qcut};
    my $contacts=();
    my $fname;
    if (defined $par{qfile}) {
      if (-r $par{qfile}) {
        $fname=$par{qfile};
      }
    }
    die "Cannot find file defining native contacts" if (! defined $fname);

    my $qin=&GenUtil::getInputFile($fname);
    while (<$qin>) {
      my $line=$_;
      chomp $line;
      push(@{$contacts}, $line);
    }
    close $qin;
    $cmd .= " ncontact " . ($#{$contacts}+1) . "\n";
    while ($#{$contacts} > -1) {
      $cmd .= pop(@{$contacts}) . "\n";
    }
    $self->_sendCommand($cmd);
    $self->{_biasstatus}->{rho}="on";
  } elsif ($par{type} eq "cons") {
    my $cmd=sprintf("scalar cons set %f select prop cons .gt. 0.00000001 end",$par{force});
    $self->_sendCommand($cmd);    
  } else {
    die "this bias type is not yet supported";
  }

}

## method: clearBias(parameters)
## clear the harmonic restraint of the type described
## in the <mark>biastype</mark> data structure.

sub clearBias {
  my $self=shift;
  my %par=@_;

  die "Bias type is not defined\n" if (!defined $par{type});

  if ($par{type} eq "rg" || $par{type} eq "rmsd") {
    die "no rg bias was found to clear"
      if ((! defined $self->{_biasstatus}->{rg}) || 
          ($self->{_biasstatus}->{rg} eq "off"));
    $self->_sendCommand("rgyr rese");
    $self->{_biasstatus}->{rg}="off";
  } elsif ($par{type} eq "rho") {
    die "no rho bias was found to clear"
      if ((! defined $self->{_biasstatus}->{rho}) || 
          ($self->{_biasstatus}->{rho} eq "off"));
    $self->_sendCommand("dmco force 0");
    $self->{_biasstatus}->{rho}="off";
  } elsif ($par{type} eq "cons") {
  } else {
    die "this bias type is not yet supported";
  }

  return 1;
}

## method: val = getBiasVal(parameters)
## return the current value of the variable to which
## the restraint was applied

sub getBiasVal {
  my $self=shift;
  my %par=@_;

  die "Bias type is not defined\n" if (!defined $par{type});

  if ($par{type} eq "rg") {
    return $self->getRg();
  } elsif ($par{type} eq "rmsd") {
    if (!$self->{_loadedref}) {
      $self->loadReference($par{pdb});
      $self->{_loadedref}=1;
    }
    return $self->getRMSD();
  } elsif ($par{type} eq "rho") {
    return $self->getRho();
  } elsif ($par{type} eq "cons") {
    return 0;
  } else {
    die "this bias type is not yet supported";
  }

  return 1;
}


## method: verbose(command)
## runs an arbitrary CHARMM command. Please note, that 
## for commands requiring multiple lines all lines have
## to be given at once as a single argument.

sub verbose {
  my $self=shift;
  my $txt=join("",@_);
  $self->_sendCommand("$txt");
}

## method: stream(command)
## sends a multiple command lines to CHARMM. It returns
## only after the last command is finished

sub stream {
  my $self=shift;
  my $commands=shift;
  $self->_sendCommand($commands);
  $self->_getCHARMMOutput("PERLDONE\n  \n");
}

## method: writePDB(file)
## has CHARMM write out the current coordinates to a file
## in PDB format.

sub writePDB {
  my $self=shift;
  my $pdbfile=shift;
  $self->_sendCommand("open unit 10 write form name \"$pdbfile\"");
  $self->_sendCommand("write coor pdb unit 10\n*");
  $self->_sendCommand("close unit 10");
  $self->_getCHARMMOutput("PERLDONE\n  \n");
}

## method: writeCRD(file)
## has CHARMM write out the current coordinates to a file
## in CRD format.

sub writeCRD {
  my $self=shift;
  my $crdfile=shift;
  $self->_sendCommand("open unit 10 write form name \"$crdfile\"");
  $self->_sendCommand("write coor card unit 10\n*");
  $self->_sendCommand("close unit 10");
  $self->_getCHARMMOutput("PERLDONE\n  \n");
}

## method: orient()
## orients the current structure centered at the origin
## and with the principal moments of inertia aligned with
## the principal axes.

sub orient {
  my $self=shift;

  $self->_sendCommand("coor orie");
}

## method: $stat = coorStats()
## gets coordinate statistics for the current structure
## and returns a data structure with the following
## elements: <mark>atoms</mark> (number of atoms),
## <mark>xmin</mark>, <mark>xmax</mark>, <mark>xcenter</mark>,
## <mark>ymin</mark>, <mark>ymax</mark>, <mark>ycenter</mark>,
## <mark>zmin</mark>, <mark>zmax</mark>, <mark>zcenter</mark>
## (cartesian dimensions).

sub coorStats {
  my $self=shift;

  $self->_sendCommand("coor stats");
  $self->_getCHARMMOutput("PERLDONE\n  \n");

  my $rec={};
  foreach my $l ( split(/\n/,$self->{_lastOutput}) ) {
    if ($l=~/STATISTICS FOR +([0-9]+)/) {
      $rec->{atoms}=$1;
    } elsif ($l=~/XMIN = +([0-9\.\-]+) +XMAX = +([0-9\.\-]+) +XAVE = +([0-9\.\-]+)/) {
      $rec->{xmin}=$1;
      $rec->{xmax}=$2;
      $rec->{xcenter}=$3;
    } elsif ($l=~/YMIN = +([0-9\.\-]+) +YMAX = +([0-9\.\-]+) +YAVE = +([0-9\.\-]+)/) {
      $rec->{ymin}=$1;
      $rec->{ymax}=$2;
      $rec->{ycenter}=$3;
    } elsif ($l=~/ZMIN = +([0-9\.\-]+) +ZMAX = +([0-9\.\-]+) +ZAVE = +([0-9\.\-]+)/) {
      $rec->{zmin}=$1;
      $rec->{zmax}=$2;
      $rec->{zcenter}=$3;
    }      
  }

  return $rec;
}

## method: $energy = poissonBoltzmann()
## calculates the non-hydrophobic contribution to the 
## solvation free energy by solving the 
## Poisson-Boltzmann equation. The only argument
## is the grid size used for the finite difference
## solution grid.

sub poissonBoltzmann {
  my $self=shift;
  $self->_getpar(@_);

  my $cstat=$self->coorStats();
  
  my $maxx=($cstat->{xmax}-$cstat->{xmin})/2.0;
  $cstat->{xcenter}=($cstat->{xmax}-$maxx);
 
  my $maxy=($cstat->{ymax}-$cstat->{ymin})/2.0;
  $cstat->{ycenter}=($cstat->{ymax}-$maxy);

  my $maxz=($cstat->{zmax}-$cstat->{zmin})/2.0;
  $cstat->{zcenter}=($cstat->{zmax}-$maxz);

  my $nxcel=2*(int(($maxx+$self->{par}->{pbdelta})/$self->{par}->{dcel}))+1;
  my $nycel=2*(int(($maxy+$self->{par}->{pbdelta})/$self->{par}->{dcel}))+1;
  my $nzcel=2*(int(($maxz+$self->{par}->{pbdelta})/$self->{par}->{dcel}))+1;

  $self->_sendCommand("PBEQ");
  $self->_sendCommand("scalar wmain = radius");
  $self->_scaleRadii("wmain");

  my $surface=($self->{par}->{smooth})?"smooth swin 0.3":"watr $self->{par}->{proberad} reen";

  if ($self->{par}->{smooth}) {
    $self->_sendCommand("scalar wmain add 0.3");
    $self->_sendCommand("scalar wmain mult 0.952");
  }

  my $epspr=(defined $self->{par}->{epspr})?
    $self->{par}->{epspr}:$self->{par}->{epsp};

  $self->_sendCommand(sprintf("SOLVE epsp $self->{par}->{epsp} epsw $self->{par}->{epsw} nclx %d ncly %d nclz %d dcel %1.4f xbce %1.4f ybce %1.4f zbce %1.4f maxi 10000 intbp %s",$nxcel,$nycel,$nzcel,$self->{par}->{dcel},$cstat->{xcenter},$cstat->{ycenter},$cstat->{zcenter},$surface));
  my $pboutsolv=$self->_processPBOutput();
  $self->_sendCommand(sprintf("SOLVE epsp $epspr epsw $self->{par}->{epsr} nclx %d ncly %d nclz %d dcel %1.4f xbce %1.4f ybce %1.4f zbce %1.4f maxi 10000 intbp %s",$nxcel,$nycel,$nzcel,$self->{par}->{dcel},$cstat->{xcenter},$cstat->{ycenter},$cstat->{zcenter},$surface));
  my $pboutvac=$self->_processPBOutput();
  $self->_sendCommand("RESET");
  $self->_sendCommand("END");

  return $pboutsolv-$pboutvac;
}


## method: $energy = atomPoissonBoltzmann(index)
## calculates the Poisson-Boltzmann energy for a single atom
## with unity charge while all other charges are set to zero
## The atom index and the finite difference solution grid size
## are required as arguments. Please note: This command does not
## restore the original charges from the force field.

sub atomPoissonBoltzmann {
  my $self=shift;
  my $inx=shift;
  my $keepcharge=shift;

  $self->_getpar(@_);

  if (defined $keepcharge && $keepcharge) {
    if (!defined $self->{_charge}) {
       $self->{_charge}=$self->getScalar("charge");
    }
  }

  $self->_sendCommand("scalar charge set 0.0");

  foreach my $k ( split(/=/,$inx) ) {
    my ($chain,$resid,$aname,$ainx)=split(/:/,$k);
    my $charge=1.0;
    if (defined $keepcharge && $keepcharge) {
      $charge=$self->{_charge}->[$ainx];
    }
    my $segid=$self->{molecule}->getResidue($resid,$chain)->{seg};
    $self->_sendCommand("scalar charge set $charge select atom $segid $resid $aname end");
  }

  return $self->poissonBoltzmann();
}

sub _processPBOutput {
  my $self=shift;

  my $eval;
  $self->_getCHARMMOutput("PERLDONE\n  \n") 
    if ($self->{_lastOutput} eq "");
  foreach my $l ( split(/\n/,$self->{_lastOutput}) ) {
    if ($l=~/Electrostatic energy \[KCAL\/MOL\] = +([0-9\.]+)/) {
      $eval=$1;
    }
  }
  return $eval;
}

## method: $list = getScalar(name[,selection])
## obtains the list of values for a scalar
## quantity in CHARMM.

sub getScalar {
  my $self=shift;
  my $name=shift;
  my $sele=shift;

  my $res=();

  $sele="all" if (!defined $sele);
  $self->_sendCommand("scalar $name show sele $sele end");
  $self->_getCHARMMOutput("PERLDONE\n  \n");
  foreach my $l ( split(/\n/,$self->{_lastOutput}) ) {
    if ($l=~/\(.*\) *(.*)$/) {
      push(@{$res},$1);
    }
  }

  return $res;
}

## method: $val = reportVariable(name)
## returns the value of a CHARMM ?variable.

sub reportVariable {
  my $self=shift;
  my $name=shift;

  my $val;

  $self->_sendCommand("set tmp$name ?$name");
  $self->_getCHARMMOutput("PERLDONE\n  \n") 
    if ($self->{_lastOutput} eq "");

  foreach my $l ( split(/\n/,$self->{_lastOutput}) ) {
    if ($l=~/Parameter: .+ <- "(.+)"/) {
      $val=$1;
    }
  }

  return $val;
}

## method: solvAccessSurf([accu[,rprobe]])
## calculates the solvent accessible surface area.
## A desired accuracy and probe radius may be
## given as argument.

sub solvAccessSurf {
  my $self=shift;
  my $accu=shift;
  my $rprobe=shift;
  
  $accu=0.005 if (!defined $accu);
  $rprobe=1.4 if (!defined $rprobe);

  $self->_sendCommand("scalar wmain = radius");
  if (defined $self->{par}->{scalerad} && $self->{par}->{scalerad}) {
    $self->_scaleRadii("wmain");
  }

  $self->_sendCommand("coor surface accu $accu rprobe $rprobe weight");
}  

## method: $sasa = getSASAOutput() 
## extracts solvent accessible surface area output after
## <mark>solvAccessSurf</mark> has been called and
## returns a data structure with the surface area in <mark>area</mark>
## and the associated hydrophobic contribution to the free energy of
## solvation in <mark>energy</mark>

sub getSASAOutput {
  my $self=shift;

  my $sasa={};

  $self->_getCHARMMOutput("PERLDONE\n  \n") 
    if ($self->{_lastOutput} eq "");
  foreach my $l ( split(/\n/,$self->{_lastOutput}) ) {
    if ($l=~/SURFAC: TOTAL = +([0-9\.]+)/) {
      $sasa->{area}=$1;
      $sasa->{energy}=$1*$self->{par}->{sasagamma}+$self->{par}->{sasadelta};
    }
  }
  return $sasa;
}

## method: logSASA()
## writes out solvent accessible surface area and the hydrophobic contribution
## to the solvation free energy to the energy log file after
## <mark>solvAccessSurf</mark> has been called.

sub logSASA {
  my $self=shift;

  if (defined $self->{handle}->{enerlog} ) {
    my $sasa=$self->getSASAOutput();
    $self->{handle}->{enerlog}->printf("SASA: %f %f\n",$sasa->{area},$sasa->{energy});
  }
}

## method: $ener = getEnergy()
## evaluates the energy of the current conformation
## and returns a data structure with the energy components

sub getEnergy {
  my $self=shift;
  my $sellist=shift;
  my $selfe=shift;

  if (!defined $sellist) {
    $self->_sendCommand("energy");
  } else {
    my $selarr=();
    foreach my $c ( @{$sellist} ) {
      my $ch=$self->{molecule}->getChain($c->{chain});
      die "invalid restraint list: cannot find chain $c->{chain}"
	if (!defined $ch);
      
      my $rfirst=$self->{molecule}->firstResNum($ch);
      my $rlast=$self->{molecule}->lastResNum($ch);

      if ($c->{to}>=$rfirst && $c->{from}<=$rlast) {
	$c->{from}=$rfirst if ($c->{from}<$rfirst);
	$c->{to}=$rlast    if ($c->{to}>$rlast);

	my $rfrom=$self->{molecule}->getResidueInChain($c->{from},$ch);
	die "invalid restraint list: cannot find first residue $c->{chain}$c->{from} in molecule"
	  if (!defined $rfrom);
    
	my $r=$ch->{res};

	my $rfrominx=$ch->{resinx}->{$c->{from}};
	my $rtoinx=$ch->{resinx}->{$c->{to}};

	die "invalid restraint list: cannot find last residue $c->{chain}$c->{to} in molecule"
	  if (!defined $rtoinx);
    
	my $sel;
	my $first=$rfrominx;
	my $segid=$rfrom->{seg};
	for (my $ir=$rfrominx+1; $ir<=$rtoinx; $ir++) {
	  if ($r->[$ir]->{seg} ne $segid) {
	    my $n1=$r->[$first]->{num};
	    my $n2=$r->[$ir-1]->{num};
	    $sel="( resid ".(($first!=$ir-1)?"$n1:$n2":"$n1")." .and. segid ".$segid." )";
	    push(@{$selarr},$sel);
	    $first=$ir;
	    $segid=$r->[$ir]->{seg};
	  }
	}

	my $n1=$r->[$first]->{num};
	$sel="( resid ".(($first!=$rtoinx)?"$n1:$c->{to}":"$n1")." .and. segid ".$segid." )";
	push(@{$selarr},$sel);
      }
    }

    my $cmd="inte select ( -\n";
    $cmd.=join(" .or. - \n",@{$selarr});
    if (defined $selfe && $selfe) {
      $cmd.=" ) end";
    } else {
      $cmd.=" ) end select all end";
    }
    $self->_sendCommand($cmd);
  }

  my $earr=$self->_processEneOutput();
  return $earr->[0];
}

## method: $rg = getRg()
## evaluates and returns the radius of gyration
## of the current conformation

sub getRg {
  my $self=shift;
  $self->_sendCommand("coor rgyr sele type ca end");

  my $currRg;

  $self->_getCHARMMOutput("PERLDONE\n  \n")
    if ($self->{_lastOutput} eq "");

  foreach my $l ( split(/\n/,$self->{_lastOutput}) ) {
    if ($l=~/Radius of gyration= *(.*) *Net/) {
      $currRg=$1;
    }
  }
  return $currRg;
}

## method: $rmsd = getRMSD()
## evaluates and returns the RMSD with respect to
## the structure in COMP

sub getRMSD { 
  my $self=shift;
  $self->_sendCommand("scalar sca1 copy x select all end");
  $self->_sendCommand("scalar sca2 copy y select all end");
  $self->_sendCommand("scalar sca3 copy z select all end");
  $self->_sendCommand("coor orie rms select type ca .or. type cb end");
  
  my $currRMSD;

  $self->_getCHARMMOutput("PERLDONE\n  \n")
    if ($self->{_lastOutput} eq "");

  foreach my $l ( split(/\n/,$self->{_lastOutput}) ) {
    if ($l=~/THUS RMS DIFF IS +([0-9\.]+)/) {
      $currRMSD=$1;
    }
  }

  $self->_sendCommand("scalar x copy sca1 select all end");
  $self->_sendCommand("scalar y copy sca2 select all end");
  $self->_sendCommand("scalar z copy sca3 select all end");

  return $currRMSD;
}

## method: $rho = getRho()
## evaluates and returns the "continuous" fraction
## of native contacts of the current conformation

sub getRho {
  my $self=shift;

  my $rho = "N/A";
  $rho = $self->reportVariable("rho")
    if (defined $self->{_biasstatus}->{rho});

  return $rho;
}

## function: pdbformat = getConvType(paramset)
## return the PDB format needed for a given a force field
## parameter set

sub getConvType {
  my $par=shift;

  my $convtype;
  if ($par eq "19" || $par eq "eef1") {
    $convtype="CHARMM19";
  } elsif ($par eq "a94") {
    $convtype="CHAMBER";
  } else {
    $convtype="CHARMM22";
  }
  return $convtype;
}


## function: initTrajectory(file)

sub initTrajectory {
  my $self=shift;
  my $fname=shift;
  my $wname=shift;
  my $title=shift;

  if (-r $fname) {
    $self->_sendCommand("open read unform unit 22 name \"$fname\"");
    if (defined $wname) {
      $self->_sendCommand("open write unform unit 77 name \"$wname\"");
      $self->_sendCommand("traj iread 22 nread 1 iwrite 77\n* $title\n*");
    } else {
      $self->_sendCommand("traj iread 22 nread 1");
    }
    
    $self->{_trajio}=22;
    $self->{_trajframes}=-1;
    $self->{_readframes}=0;
    $self->{_trajfreq}=0;
  } else {
    printf STDERR "file $fname cannot be read\n";
  }
}

## function: nextFrame 

sub nextFrame {
  my $self=shift;
  die "need to initialize trajectory first" 
    if (!defined $self->{_trajio});
  return 0 if ($self->{_readframes}++>=$self->{_trajframes} && $self->{_trajframes}>0);
  $self->_sendCommand("traj read");
  if ($self->{_trajframes}<=0) {
    $self->_getCHARMMOutput("PERLDONE\n  \n")
     if ($self->{_lastOutput} eq "");
    foreach my $l ( split(/\n/,$self->{_lastOutput}) ) {
     if ($l=~/NUMBER OF COORDINATE SETS IN FILE: +([0-9]+)/) {
        $self->{_trajframes}=$1;
     } elsif ($l=~/FREQUENCY FOR SAVING COORDINATES: +([0-9]+)/) {
        $self->{_trajfreq}=$1;
     } elsif ($l=~/NUMBER OF PREVIOUS DYNAMICS STEPS: +([0-9]+)/) {
        $self->{_trajprev}=$1;
     }
    } 
    $self->{_firstframe}=$self->{_trajprev}/$self->{_trajfreq};
  }
  return 1;
}
 
## function: closeTrajectory

sub closeTrajectory {
  my $self=shift;
   
  die "need to initialize trajectory first" 
    if (!defined $self->{_trajio});
  $self->_sendCommand("close unit $self->{_trajio}");
  undef $self->{_trajio};
}

## function: ($values)=parseOutput($rex)

sub parseOutput {
  my $self=shift;
  my $rex=shift;

  $self->_getCHARMMOutput("PERLDONE\n  \n") if ($self->{_lastOutput} eq "");
 
  if ($self->{_lastOutput}=~/$rex/) {
    return ($1,$2,$3,$4,$5,$6,$7,$8,$9);
  } else {
    return (undef,undef,undef,undef,undef,undef,undef,undef,undef);
  }
}

## function: ($ener,$sasa) = readEnergy(file)
## reads an energy log file and returns the
## data in two data structures. The first
## is a list of entries with the individual
## energy components . The second one
## contains the solvent accesible surface area
## and hydrophobic free energy of solvation,
## if available from the log file.

sub readEnergy {
  my $self=();

  my $elog=&GenUtil::getInputFile(shift);

  my ($sasa,$esasa)=(0.0,0.0);

  while (<$elog>) {
    chomp;
    s/^ +//;

    my ($tag,$rest)=split(/:/);
      $rest=~s/^ +//;

    if ($tag eq "SASA") {
      ($sasa,$esasa)=split(/ +/,$rest);
    } else {
      my $rec={};
      ($rec->{step},$rec->{total},$rec->{vdwaals},
       $rec->{elec},$rec->{gb},$rec->{asp},$rec->{constr})=split(/ +/,$rest);
      push (@{$self},$rec);
    }
  }

  undef $elog;

  return ($self,$esasa);
}

### process MD energy output ######

sub _processMDEneOutput {
  my $self=shift;
  
  my $earr=();
  
  $self->_getCHARMMOutput("PERLDONE\n  \n")
    if ($self->{_lastOutput} eq "");

  my $erec;
  my $save=0;

  foreach my $l ( split(/\n/,$self->{_lastOutput}) ) {
    if ($l=~/^DYNA>/) {
      $save=1;
      $erec={};
      $erec->{step}=substr($l,8,6);
      $erec->{time}=substr($l,14,13);
      $erec->{total}=substr($l,27,13);
      $erec->{kine}=substr($l,40,13);
      $erec->{pot}=substr($l,53,13);
      $erec->{temp}=substr($l,66,13);

      foreach my $n ( qw ( bonds angles cmap ureyb dihedrals impropers vdwaals elec gb asp volume constr ) ) {
	$erec->{$n}=0.0;
      }
    } elsif ($l=~/^DYNA INTERN>/) {
      $erec->{bonds}=substr($l,14,13);
      $erec->{angles}=substr($l,27,13);
      $erec->{ureyb}=substr($l,40,13);
      $erec->{dihedrals}=substr($l,53,13);
      $erec->{impropers}=substr($l,66,13);
    } elsif ($l=~/^DYNA CROSS>/) {
      $erec->{cmap}=substr($l,14,13);
    } elsif ($l=~/^DYNA EXTERN>/) {
      $erec->{vdwaals}=substr($l,14,13);
      $erec->{elec}=substr($l,27,13);
      $erec->{asp}=substr($l,53,13);
    } elsif ($l=~/^DYNA PBEQ>/) {
      $erec->{gb}=substr($l,40,13);
    } elsif ($l=~/^DYNA PRESS>/) {
      $erec->{volume}=substr($l,66,13);
    } elsif ($l=~/^DYNA CONSTR>/) {
      $erec->{constr}=substr($l,14,13);
    } elsif ($l=~/^ +DYNA A     =/) {
      $erec->{boxsize}=substr($l,19,10);
    } else {
      if ($save) {
	foreach my $n ( keys %{$erec} ) {
	  $erec->{$n}=0.0 if ($erec->{$n}=~/\*/);
	}
	push(@{$earr},$erec);
	$save=0;
      }
    }
  }
  return $earr;
}

### process PERT output ######

sub _processPertOutput {
  my $self=shift;
  
  my $earr=();
  
  $self->_getCHARMMOutput("PERLDONE\n  \n")
    if ($self->{_lastOutput} eq "");

  my $erec;
  my $save=0;

  foreach my $l ( split(/\n/,$self->{_lastOutput}) ) {
    if ($l=~/^ PERTURBATION> results,/) {
      $save=1;
      $erec={};
      $erec->{lstart}=substr($l,31,12);
      $erec->{lstop}=substr($l,51,12);
      $erec->{llast}=substr($l,71,12);
      $erec->{steps}=substr($l,105,10);
    } elsif ($l=~/^ PERTURBATION> result:/) {
      $erec->{expave1}=substr($l,30,12);
      $erec->{expave2}=substr($l,43,12);
      $erec->{expflc1}=substr($l,63,12);
      $erec->{expflc2}=substr($l,76,12);
      $erec->{difave}=substr($l,96,12);
      $erec->{difflc}=substr($l,116,12);
    } elsif ($l=~/^ PERTURBATION> TP/) {
      $erec->{tptot}=substr($l,43,12);
      $erec->{tpforward}=substr($l,65,12);
      $erec->{tpbackward}=substr($l,88,12);
    } elsif ($l=~/^ PERTURBATION> TI/) {
      $erec->{titot}=substr($l,43,12);
      $erec->{tiforward}=substr($l,65,12);
      $erec->{tibackward}=substr($l,88,12);
    } else {
      if ($save) {
	foreach my $n ( keys %{$erec} ) {
	  $erec->{$n}=0.0 if ($erec->{$n}=~/\*/);
	}
	push(@{$earr},$erec);
	$save=0;
      }
    }
  }
  return $earr;
}


### process energy output ######

sub _processEneOutput {
  my $self=shift;
  
  my $earr=();

  $self->_getCHARMMOutput("PERLDONE\n  \n") 
    if ($self->{_lastOutput} eq "");

  my $erec;
  my $save=0;

  foreach my $l ( split(/\n/,$self->{_lastOutput}) ) {
    if ($l=~/^(MINI|ENER|INTE)>/) {
      $save=1;
      $erec={};
      $erec->{step}=substr($l,8,6);
      $erec->{total}=substr($l,14,13);
      foreach my $n ( qw ( bonds cmap angles ureyb dihedrals impropers vdwaals asp elec gb constr ) ) {
	$erec->{$n}=0.0;
      }
    } elsif ($l=~/^(MINI|ENER|INTE) INTERN>/) {
      $erec->{bonds}=substr($l,14,13);
      $erec->{angles}=substr($l,27,13);
      $erec->{ureyb}=substr($l,40,13);
      $erec->{dihedrals}=substr($l,53,13);
      $erec->{impropers}=substr($l,66,13);
    } elsif ($l=~/^(MINI|ENER|INTE) CROSS>/) {
      $erec->{cmap}=substr($l,14,13);
    } elsif ($l=~/^(MINI|ENER|INTE) EXTERN>/) {
      $erec->{vdwaals}=substr($l,14,13);
      $erec->{elec}=substr($l,27,13);
      $erec->{asp}=substr($l,53,13);
      $erec->{user}=substr($l,66,13);
    } elsif ($l=~/^(MINI|ENER|INTE) PBEQ>/) {
      $erec->{gb}=substr($l,40,13);
    } elsif ($l=~/^(MINI|ENER|INTE) CONSTR>/) {
      $erec->{constr}=substr($l,14,13);
    } elsif ($l=~/^(MINI|ENER|INTE) SOLVAT>/) {
      $erec->{solvation}=substr($l,14,13);
    } elsif ($l=~/^(MINI|ENER|INTE)   ACE1>/) {
      $erec->{acehyd}=substr($l,14,13);
      $erec->{aceself}=substr($l,27,13);
      $erec->{acescreen}=substr($l,40,13);
      $erec->{acecoulomb}=substr($l,53,13);
    } elsif ($l=~/^(MINI|ENER|INTE)   ACE2>/) {
      $erec->{acesolv}=substr($l,14,13);
      $erec->{aceinter}=substr($l,27,13);
    } elsif ($l=~/^(MINI|ENER|INTE) EWALD>/) {
      $erec->{ewksum}=substr($l,14,13);
      $erec->{ewself}=substr($l,27,13);
      $erec->{ewexcl}=substr($l,40,13);
    } else {
      if ($save) {
	foreach my $n ( keys %{$erec} ) {
	  $erec->{$n}=0.0 if ($erec->{$n}=~/\*/);
	}
	push(@{$earr},$erec);
	$save=0;
      }
    }
  }
  return $earr;
}

### get CHARMM output ######

sub _getCHARMMOutput {
  my $self=shift;
  my $expect=shift;

  $self->{_lastOutput}="";

  my $han='';
  vec($han,fileno($self->{handle}->{fromcharmm}),1)=1;
  
  do {
    my $rout;
    select($rout=$han,undef,undef,undef);

    my $line="";
    sysread($self->{handle}->{fromcharmm},$line,100000);

    $self->{_lastOutput}.=$line;
    $self->{handle}->{outlog}->print($line)
      if (defined $self->{handle}->{outlog});

    if ($self->{_lastOutput} =~ /ABNORMAL TERMINATION/ || 
       $self->{_lastOutput} =~ /TERMINATING/) { 
      printf STDERR "*** CHARMM terminated abnormally\n";
      printf STDERR "*** Last command: $self->{_lastCommand}\n";
      printf STDERR "*** CHARMM output follows:\n";
      my $shortoutput=$self->{_lastOutput};
      $shortoutput=~s/[ \t]+\/---------.*//s;
      printf STDERR $shortoutput;
      $self->_closeAll();
      exit (1);
    }
  } until (!defined $expect || $expect eq "" || $self->{_lastOutput} =~ /$expect/);
}

### wait for CHARMM prompt and send command ######

sub _sendCommand { 
  my $self=shift;
  my $cmd=shift;
  my $dontsendlabel=shift;

  $self->_getCHARMMOutput("PERLDONE\n  \n") 
    if ($self->{_lastOutput} eq "");

  $self->{handle}->{tocharmm}->print($cmd,"\n");
  $self->{handle}->{cmdlog}->print($cmd,"\n")
    if (defined $self->{handle}->{cmdlog});

  $self->{_lastCommand}=$cmd;
  $self->{_lastOutput}="";

  $self->{handle}->{tocharmm}->print("label PERLDONE\n")
    unless (defined $dontsendlabel && $dontsendlabel);
}

### setParameter ######

sub setParameter {
  my $self=shift;
  $self->_getpar(@_);
}

### _getpar ######

sub _getpar {
  my $self=shift;
  my %arg=@_;
  foreach my $n ( keys %arg ) {
    if (exists $self->{par}->{$n}) {
      $self->{par}->{$n}=$arg{$n};
      if (($n eq "eef1" || $n eq "sasa") && $arg{$n}) {
        $self->{par}->{param}="eef1";
      }
    } else {
      printf STDERR "Unknown CHARMM parameter $n will be ignored!\n";
    }
  }
}

### _getSel ######

sub _getSel {
  my $a=lc shift;
  
  if ($a eq "ca") {
    return "atom * * ca";
  } elsif ($a eq "cb") {
    return "atom * * cb";
  } elsif ($a eq "cab") {
    return "(atom * * ca .or. atom * * cb)";
  } elsif ($a eq "heavy") {
    return ".not. hydrogen";
  } else {
    return "all";
  }
}

### _closeAll ######

sub _closeAll {
  my $self=shift;

  close $self->{handle}->{tocharmm};
  close $self->{handle}->{fromcharmm};

  waitpid($self->{_charmmpid},0);

  undef $self->{handle}->{cmdlog} if (defined $self->{handle}->{cmdlog});
  undef $self->{handle}->{outlog} if (defined $self->{handle}->{outlog});
  undef $self->{handle}->{enerlog} if (defined $self->{handle}->{enerlog});
}


1;

