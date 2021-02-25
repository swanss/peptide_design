# Author: Gevorg Grigoryan
# Start date: August 25, 2003
# This module defines all the common paths and variables
# needed by scripts and other packages written in our lab.
package DEFINITIONS;
#use strict;
use File::Spec;
use Exporter ();
@ISA = qw(Exporter);

@EXPORT = qw($ROOT_DEF $BIN_DEF $SCRIPT_DEF $LIB_DEF $PERLLIB_DEF $PARAM_DEF $CUSTOM_PARAM_DEF $DELPHI_DEF $CHARMM_DEF $DEF_VERSION
             $ROTLIB_DEF $FS_DEF $MULTE_DEF $DEF_VER_CONTROL $STRIDE_DEF $NACCESS_DEF $PRANK_DEF $NPROC_DEF $MPI_COMM_WORLD_DEF
             $ASURF_DEF $PEP21_DEF $PEP11_DEF $PEP5_DEF $LOCAL_DEF $LOCAL_DEF_ROOT $LOCAL_TIMEOUT $NFS_DELAY);

# Version control mechanism:
# DEFINITIONS must be included first in the root script, so that
# it sets up the version hash table. Each module is responsible for
# modifying the entries in the table corresponding to it when it is included.
my %DEF_VER_CONTROL;
$DEF_VER_CONTROL = \%DEF_VER_CONTROL;
#my %ver;
#$DEF_VER_CONTROL->{ver} = \%ver;
my %date;
$DEF_VER_CONTROL->{date} = \%date;
#$DEF_VER_CONTROL->{ver}->{DEFINITIONS} = "1.0";
$DEF_VER_CONTROL->{date}->{DEFINITIONS} = "06/17/04";

# First define all global paths
# -- find root directory
#$ROOT_DEF = "/home/gevorg/work/ProCEDe/trunk/";
$ROOT_DEF = $INC{"DEFINITIONS.pm"}; $ROOT_DEF =~ s/\/\//\//g;
(undef, $ROOT_DEF, undef) = File::Spec->splitpath($ROOT_DEF);
$ROOT_DEF = `cd $ROOT_DEF../; pwd`; chomp($ROOT_DEF); $ROOT_DEF .= "/";

$BIN_DEF = "$ROOT_DEF/bin";
$SCRIPT_DEF = "$ROOT_DEF/scripts/";
$LIB_DEF = "$ROOT_DEF/lib/";
$PERLLIB_DEF = "$ROOT_DEF/modules";
# ..the default parameter file path
$PARAM_DEF = "$ROOT_DEF/param";
# ..if the user specifies a custom path, this variable will change
$CUSTOM_PARAM_DEF = "/dev/null";
$ROTLIB_DEF = undef;
$LOCAL_DEF = undef;
$LOCAL_DEF_ROOT = undef;
$LOCAL_TIMEOUT = undef;

# Then define paths to commonly used executables
$DELPHI_DEF = "$BIN_DEF/delphi";
$CHARMM_DEF = "$BIN_DEF/charmm.33b1";
$MULTE_DEF = "$BIN_DEF/multe";
$STRIDE_DEF = "$BIN_DEF/stride";
$NACCESS_DEF = "$BIN_DEF/naccess/naccess";
$ASURF_DEF = "$BIN_DEF/pep/asurf";
$PEP21_DEF = "$BIN_DEF/pep/pep21";
$PEP11_DEF = "$BIN_DEF/pep/pep11";
$PEP5_DEF = "$BIN_DEF/pep/pep5";

# File space hash table
my %FS_DEF;
$FS_DEF = \%FS_DEF;

# Constants relating to parallel runs
$PRANK_DEF = 0; # the rank of the process
$NPROC_DEF = 1; # the number of processes
$MPI_COMM_WORLD_DEF = 0; # global world communicator
$NFS_DELAY = 10; # number of seconds to wait for NFS to "catch" up with changes after extensive writing



=head2

 Title   : getParmFile
 Usage   : DEFINITIONS::getParmFile($filename);
 Function: Returns the full path to the parameter file.
           First checks under the custom parameter directory (if
           the user has specified any) and then in the standard
           parameter directory.
 Returns : Full path to parameter file. Exits with an error if the
           file is not found.
 Args    : Parameter file name.

=cut

sub getParmFile {
  my $fname = shift;
  my $sentinel = shift;
  
  if ((defined($LOCAL_DEF)) && (-e "$LOCAL_DEF/$fname")) {
    return File::Spec->rel2abs("$LOCAL_DEF/$fname");
  } elsif (-e "$CUSTOM_PARAM_DEF/$fname") {
    if (defined($LOCAL_DEF) && (-d "$LOCAL_DEF")) {
      GENERAL::csystem("cp $CUSTOM_PARAM_DEF/$fname $LOCAL_DEF/$fname");
      return File::Spec->rel2abs("$LOCAL_DEF/$fname");
    } else {
      return "$CUSTOM_PARAM_DEF/$fname";
    }
  } elsif (-e "$PARAM_DEF/$fname") {
    if (defined($LOCAL_DEF) && (-d "$LOCAL_DEF")) {
      GENERAL::csystem("cp $PARAM_DEF/$fname $LOCAL_DEF/$fname");
      return File::Spec->rel2abs("$LOCAL_DEF/$fname");
    } else {
      return "$PARAM_DEF/$fname";
    }
  } else {
    return $sentinel if (defined($sentinel));
    GENERAL::error("Could not find parameter file \"$fname\"!");
  }
}


=head2

 Title   : getExec
 Usage   : DEFINITIONS::getExec($binfilename);
 Function: Returns the full path to the executable file.
           First checks if a local copy of the file already exists, and if not,
           then checks under the standard binary directory and copies the file locally.
 Returns : Full path to executable file. Exits with an error if the
           file is not found.
 Args    : Executable file name.

=cut

sub getExec {
  my $fname = shift;
  my $sentinel = shift;
  
  if ((defined($LOCAL_DEF)) && (-e "$LOCAL_DEF/$fname")) {
    return File::Spec->rel2abs("$LOCAL_DEF/$fname");
  } elsif (-e "$BIN_DEF/$fname") {
    if (defined($LOCAL_DEF) && (-d "$LOCAL_DEF")) {
      GENERAL::csystem("cp $BIN_DEF/$fname $LOCAL_DEF/$fname");
      return File::Spec->rel2abs("$LOCAL_DEF/$fname");
    } else {
      return "$BIN_DEF/$fname";
    }
  } else {
    return $sentinel if (defined($sentinel));
    GENERAL::error("Could not find executable file \"$fname\"!");
  }
}

1;
