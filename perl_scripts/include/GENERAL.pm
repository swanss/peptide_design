=head1 NAME GENERAL

GENERAL - a collection of general utility 

=head1 SYNOPSIS

someone needs to write it

=head1 DESCRIPTION

Functions: 
(1) 

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




package GENERAL;
#$DEF_VER_CONTROL->{ver}->{GENERAL} = "1.0";
$DEF_VER_CONTROL->{date}->{GENERAL} = "06/17/04";

use strict;
use IO::File;
use IO::Handle;
#use DESIGN;
use Data::Dumper;
use DEFINITIONS;
use Text::Wrap;
use Math::Trig;
#use Parallel::MPI::Simple;
use POSIX ":sys_wait_h";
use Fcntl;
use POSIX qw(:unistd_h :errno_h);
#use Term::ReadKey;

#--------------------------------------------------------------
#-------------------Section I: File Operation------------------
#--------------------------------------------------------------

=head2

 Title   :  usage
 Usage   :  GENERAL::usage("Performs Magic", "Required:", "", "-m", "amount of magic", "-l", "length", "Optional:", "-d", "dissappearing rabits");
 Function:  Returns usage string properly formatted.
 Returns :  Returns usage string properly formatted.
 Args    :  1. Program title (a summary of what it does)
            2. Option name, option description combinations
            3. If option description is empty, treated as a new category name

=cut
sub usage {
  my $title = shift;
  GENERAL::requireArgs($title);
  my @opts = @_;
  if (scalar(@opts)%2 != 0) {
    GENERAL::error("Opts array must have even number of elements");
  }
  
  my $tab = "   ";
  my $sep = " - ";
#  my ($wchar, $hchar, $wpixels, $hpixels) = Term::ReadKey::GetTerminalSize();
#  $wchar = $wchar - 2;
  my $cols = $Text::Wrap::columns;
#  $Text::Wrap::columns = $wchar;
  my $usage = Text::Wrap::wrap("", "", "\n$title\n");
  for (my $i = 0; $i < scalar(@opts); $i = $i + 2) {
    my $opn = $opts[$i];
    my $opt = $opts[$i+1];
    if ($opt =~ /^\s*$/) {
      $usage .= Text::Wrap::wrap("", "", "\n$opn\n");
    } elsif ($opn =~ /^\s*$/) {
      $opt = Text::Wrap::wrap("$tab", "$tab", $opt);
      $usage .= "$opt\n";
    } else {
#      $Text::Wrap::columns = $wchar - length($tab) - length($sep);
      $Text::Wrap::columns = $cols - length($tab) - length($sep) - length($opn);
      $opt = Text::Wrap::wrap("", "", $opt);
      my $par = "$tab" . (" " x (length($sep) + length($opn)));
      $opt =~ s/\n/\n$par/g;
      $usage .= "$tab$opn$sep$opt\n";
#      $Text::Wrap::columns = $wchar;
      $Text::Wrap::columns = $cols;
    }
  }
  $usage .= "\n";
  $Text::Wrap::columns = $cols;
  return $usage;
}

# function: GetBase(string)
# description: Given a file name, get the base for that file name.
#              That is, remove the last extension of the file name.
# example:     blah.pdb  ---> blah
# example:     blah.pdb.jpg ---> blah.pdb
# usage: &GetBase("blah.pdb");

sub GetBase {
    my $fname = shift;
    my $ext = shift;
    my $fn = shift;
    
    (my $basename = $fname) =~ s/\.([^\.\/]*)$//;
    if (defined($ext)) {
      if (defined($1)) { $$ext = $1; }
      else { $$ext = ""; }
    }
    if (defined($fn)) {
      $basename =~ s/^.+\/([^\/]+)$/$1/g;
    }
    return $basename;
}


#-----------------------------------------------------------------------------------------
# Author: Gevorg Grigoryan
# Date: 09/08/03
# Description: A set of "careful" system commands.
# Careful means that they check for errors and if any occur, they
# halt the program. These should always be called when the result
# of the system command is at all important for the program.

# A careful general system call. Halts the program if an error occurs.
sub csystem {
  my $cmd = shift;
  my $errmsg = shift;
  my $success = shift;
  if (!defined($success)) { $success = 0; }
  my $strict = shift;
  if (!defined($strict)) { $strict = 1; }

#  undef $!;
  my $ret = system($cmd);
  if ($ret != $success) {
    if ($strict > 0) {
      if (defined($errmsg)) {
        GENERAL::error($errmsg);
      } else {
        GENERAL::error("System error occured while trying to execute command \"$cmd\". Return code is $ret.");
      }
    } elsif ($strict == 0) {
      if (defined($errmsg)) {
        GENERAL::warning($errmsg);
      } else {
        GENERAL::warning("System error occured while trying to execute command \"$cmd\". Return code is $ret.");
      }
      return $ret;
    } else {
      return $ret;
    }
  }
  return $ret;
}

=head2

 Title   :  sys
 Usage   :  my ($ret, $out) = GENERAL::sys($cmd);
 Function:  Executes an arbitrary command and returns both the error code and the standard output.
 Returns :  Returns the error code and the standard output.
 Args    :  1. command string

=cut
sub sys {
  my $cmd = shift;
  my $out = `$cmd`;
  return ($?, $out);
}

=head2

 Title   :  bashExecute
 Usage   :  my $ret = GENERAL::bashExecute($cmd);
 Function:  Executes an arbitrary command in a shell after executing the user's .bashrc file. Thus,
            most commands that are valid in an interactive bash shell can be executed with this function.
 Returns :  Returns the standard output of the command.
 Args    :  1. command to execute in bash
            2. optional: expected return code (upon success). Default is 0.

=cut
sub bashExecute {
  my $com = shift;
  my $suc = shift;
  my $out = shift;
  GENERAL::requireArgs($com);
  $com = GENERAL::escape($com, "\"");
  $suc = 0 if (!defined($suc));
  my $delay = shift;
  $delay = 0 if (!defined($delay));
  my $mach = GENERAL::GetMachine();
  my $tmpf = "__123__bashrc_tmp__$$\_$mach";

  my $user = GENERAL::GetUser();
  sleep($delay) if ($delay);
  system("cp /u/$user/.bashrc $tmpf");
  sleep($delay) if ($delay);
  system("echo \"" . GENERAL::escape($com) . "\" >> $tmpf");
  sleep($delay) if ($delay);
  my $ret = undef;
  # if the output is asked, collect it; otherwise send it to stdout
  if (defined(wantarray)) { $ret = `bash $tmpf`; }
  else { system("bash $tmpf"); }
  if ($? !~ /^$suc$/) { GENERAL::error("error while trying to execute command '$com' through bash (error code = $?" . (defined($ret) ? ", return value $ret)" : ")")); }
  chomp($ret) if (defined($ret));
  GENERAL::csystem("rm $tmpf");
  sleep($delay) if ($delay);
  return $ret;
}


=head2

 Title   :  crm
 Usage   :  GENERAL::crm($filename);
            GENERAL::crm(@manyfiles);
 Function:  A carefull system call to remove files (uses unlink)
 Returns :  
 Args    :  1. either a file name or an array of file names

=cut
sub crm {
#  $! = "";
  my $num = unlink(@_);
  if ($num != scalar(@_)) {
    GENERAL::error("System error occured while trying to delete file(s): \'" . join("', '", @_) . "'");
    exit(-1);
  }
}


=head2

 Title   :  crmdir
 Usage   :  GENERAL::crmdir($dirname);
            GENERAL::crmdir(@manydirs);
 Function:  A carefull system call to remove directories (uses rmdir)
 Returns :
 Args    :  1. either a file name or an array of file names

=cut
sub crmdir {
#  $! = "";
  foreach my $dir (@_) {
    if (!rmdir($dir)) {
      GENERAL::error("Error removing directory '$dir'.");
    }
  }
}


# A careful current directory change. Halts the program if an error occurs.
# (Sometimes you want to cd into a subdirectory and remove everything there
# silently. It's bad news if you misstype the directory name and it does not
# exist! This way we can avoid that.)
sub cchdir {
  my $dir = shift;
#  $! = "";
  my $ret = chdir($dir);
  if ($ret == 0) {
    GENERAL::error("System error occured while trying to cd to \"$dir\"");
  }
}

# A careful make direcotry. Halts the program if an error occurs.
sub cmkdir {
  my $dir = shift;
  my $rec = shift;
  my $strict = shift;
  
#  $! = "";
  if (-d $dir) {
    if (defined($strict) && ($strict > 1)) {
      GENERAL::error("Error in cmkdir: tried to create directory \"$dir\" while it already exists.\nIgnoring...");
    } elsif (defined($strict) && ($strict == 1)) {
      GENERAL::warning("Warning in cmkdir: tried to create directory \"$dir\" while it already exists.\nIgnoring...");
    } else {
      return;
    }
  }

  if (defined($rec) && ($rec)) {
    $dir = GENERAL::Trim($dir);
    my @dirs = split(/\//, $dir);
    my $path = ($dir =~ /^\// ? "/" : "");
    foreach my $d (@dirs) {
      if ($d ne "") {
        if ($path eq "") { $path = $d; }
        else { $path = $path . "/$d"; }
        if (! -d $path) {
          my $ret = mkdir($path);
          if ($ret == 0) {
            GENERAL::error("System error occured while trying to create directory \"$dir\"");
            exit(-1);
          }
        }
      }
    }
  } else {
    my $ret = mkdir($dir);
    if ($ret == 0) {
      GENERAL::error("System error occured while trying to create directory \"$dir\"");
      exit(-1);
    }
  }
}

# Saved the specified structure into the specified file name.
# Expects a reference and a file name.
sub SaveStruct {
  my $struct = shift;
  my $fname = shift;

  $Data::Dumper::Purity = 1;
  open (OUT, ">$fname") or die "Cannot open \"$fname\": $!";
  print OUT Data::Dumper->Dump([$struct], ['*___123___struct']);
  close OUT;
}


# Load the data structure saved in the specified file.
# Expects file name and optional flag. If set to 1, will
# return references rather than the structures themselves.
sub LoadStruct {
  my $fname = shift;
  my $ref = shift;
  if (!defined $ref) {$ref=0;}

  # Read in code
  my $fh = GENERAL::GetInFH($fname);
  my $tmp = $/;
  undef $/;  # read in all text of the dump file as one string
  my $code = <$fh>;
  close($fh);
  $/ = $tmp;

  # Remove comments and blank lines
  $code =~ s/#.*$//mg; # remove comments
  $code =~ s/^\s*$//mg; # make blank lines empty
  $code =~ s/\n\n+/\n/mg; # remove empty lines
  $code = GENERAL::Trim($code); # to remove first empty line if any

  # Check what was saved
  if ($code =~ /^\s*(\$\S+)[\s=]/) {
#    my $struct = eval("return my $code");
    my $struct = eval("my $code\nreturn $1;\n");
    GENERAL::error("Cannot create data structure from \"$fname\": $@") if ($@);
    if ($ref) { return \$struct; }
    else { return $struct; }
  } elsif ($code =~ /^\s*(%\S+)[\s=]/) {
#    my %struct = eval("return my $code");
    my %struct = eval("my $code\nreturn $1;\n");
    GENERAL::error("Cannot create data structure from \"$fname\": $@") if ($@);
    if ($ref) { return \%struct; }
    else { return %struct; }
  } elsif ($code =~ /^\s*(@\S+)[\s=]/) {
#    my @struct = eval("return my $code");
    my @struct = eval("my $code\nreturn $1;\n");
    GENERAL::error("Cannot create data structure from \"$fname\": $@") if ($@);
    if ($ref) { return \@struct; }
    else { return @struct; }
  } else {
    GENERAL::error("File $fname does not seem to contain a saved structure! The extracted code was:\n$code\n");
  }

}

=head2

 Title   :  redirectSTDOUT
 Usage   :  my $old = GENERAL::redirectSTDOUT($filename);
            GENERAL::redirectSTDOUT($old, 1);
 Function:  Redirects standard out into the given file
 Returns :  Data structure for going back to outputting to standard out.
 Args    :  1. File name to redict standard output into.
            2. optional: if set, redirects back to normal stadrd out,
               interprets the first argument as the output from the first call.

=cut

sub redirectSTDOUT {
  my $par = shift;
  my $rev = shift;
  $rev = 0 if (!defined($rev));

  if ($rev == 0) {
    my $oldSO;
    open($oldSO, ">&STDOUT") or die "Can't dup STDOUT: $!";
    open(STDOUT, '>', $par) or die "Can't redirect STDOUT to $par: $!";
    select STDOUT; $| = 1;  # make unbuffered
    return $oldSO;
  } else {
    close(STDOUT);
    open(STDOUT, ">&", $par) or die "Can't dup \$oldSO: $!";
  }
}


=head2

 Title   :  isInteger
 Usage   :  GENERAL::isInteger($num);
 Function:  Returns 1 if the passed string is an integer, 0 otherwise.
 Returns :  1/0.
 Args    :  string or number

=cut
sub isInteger {
  my $num = shift;
  return ($num =~ /(^\s*[-+]?[1-9][0-9]*\s*$|^\s*[-+]?0\s*$)/) ? 1:0;
}

=head2

 Title   :  isNumeric
 Usage   :  GENERAL::isNumeric($num);
 Function:  Returns 1 if the argument is a valid number (as per perl), 0 otherwise.
 Returns :  1/0.
 Args    :  string or number

=cut
sub isNumeric {
  my($cand) = shift;

  my $not_num = 0;
  local $^W = 1;
  local $SIG{__WARN__} = sub {
    $not_num = $_[0] =~ /^Argument ".*?" isn't numeric/;
  };
  () = $cand + 0;
  return !$not_num;
}

=head2

 Title   :  isNumericArray
 Usage   :  GENERAL::isNumericArray($num);
 Function:  Returns 1 if the argument is a valid number (as per perl), 0 otherwise.
 Returns :  1/0.
 Args    :  string or number

=cut
sub isNumericArray {
  my $arr = shift;
  foreach my $val (@$arr) {
    return 0 if (!GENERAL::isNumeric($val));
  }
  return 1;
}

=head2

 Title   :  isReNumeric
 Usage   :  GENERAL::isReNumeric($num);
 Function:  Returns 1 if the argument is a real number, 0 otherwise.
 Returns :  1/0.
 Args    :  string or number

=cut
sub isReNumeric {
  my $cand = shift;
  return (isNumeric($cand) && ($cand !~ /nan/) && ($cand !~ /inf/));
}



=head2

 Title   :  round
 Usage   :  my $int = GENERAL::round($double);
 Function:  Rounds to the nearest integer. Works for negative and positive numbers.
 Returns :  Integer value.
 Args    :  (1) Number to round.
=cut

sub round {
  my $num = shift;
  return int($num + .5 * ($num <=> 0));
}


=head2

 Title   :  floor
 Usage   :  my $int = GENERAL::floor($double);
 Function:  Rounds to the nearest integer in the direction of -inf.
 Returns :  Integer value.
 Args    :  (1) Number to round.
=cut

sub floor {
  my $num = shift;
  if (GENERAL::isInteger($num)) { return $num; }
  elsif ($num > 0) { return int($num); }
  else { return int($num) - 1; }
}


=head2

 Title   :  ceil
 Usage   :  my $int = GENERAL::ceil($double);
 Function:  Rounds to the nearest integer in the direction of +inf.
 Returns :  Integer value.
 Args    :  (1) Number to round.
=cut

sub ceil {
  my $num = shift;
  if (GENERAL::isInteger($num)) { return $num; }
  elsif ($num < 0) { return int($num); }
  else { return int($num) + 1; }
}


=head2

 Title   :  min
 Usage   :  GENERAL::min(\@array);
 Function:  Finds the smallest value in the array and its index.
 Returns :  the smallest value
 Args    :  (1) A reference to an array values.
            (2) Optional: a reference to a variable. If specified,
                the source variable will be set to the index of the
                first occurance of the smallest value of the array.
=cut

sub min {
  my $arr;
  my $ind;
  if (ref($_[0]) eq "ARRAY") {
    $arr = shift;
    $ind = shift;
  } else {
    $arr = \@_;
    $ind = undef;
  }

  my $min = $arr->[0];
  if (defined($ind)) { $$ind = 0; }
  for (my $i = 1; $i < scalar(@$arr); $i++) {
    if ($arr->[$i] < $min) {
      $min = $arr->[$i];
      if (defined($ind)) { $$ind = $i; }
    }
  }
  return $min;
}

=head2

 Title   :  max
 Usage   :  GENERAL::max(\@array);
 Function:  Finds the largest value in the array and its index.
 Returns :  the smallest value
 Args    :  (1) A reference to an array values.
            (2) Optional: a reference to a variable. If specified,
                the source variable will be set to the index of the
                first occurance of the largest value of the array.
=cut

sub max {
  my $arr;
  my $ind;
  if (ref($_[0]) eq "ARRAY") {
    $arr = shift;
    $ind = shift;
  } else {
    $arr = \@_;
    $ind = undef;
  }

  my $max = $arr->[0];
  if (defined($ind)) { $$ind = 0; }
  for (my $i = 1; $i < scalar(@$arr); $i++) {
    if ($arr->[$i] > $max) {
      $max = $arr->[$i];
      if (defined($ind)) { $$ind = $i; }
    }
  }
  return $max;
}

=head2

 Title   :  corrCoef
 Usage   :  GENERAL::corrCoef(\@array1, \@array2);
 Function:  Finds the the correlation coefficient between two equally sized arrays
            of numbers. The correlation coefficient is defined as cov(x, y)/(std(x)*std(y)).
            The range of values are from -1 (perfect anti-correlation) to 1 (perfect correlation).
 Returns :  Correlation coefficient.
 Args    :  (1) First array.
            (2) Second array.
=cut

sub corrCoef {
  my $arr1 = shift;
  my $arr2 = shift;
  GENERAL::requireArgs($arr1, $arr2);
  
  my $num = GENERAL::cov($arr1, $arr2);
  my $den = ((GENERAL::stdev($arr1))*(GENERAL::stdev($arr2)));
  if ($num*$den == 0) {
    return -100000;
  } else {
    return $num/$den;
  }
}


=head2

 Title   :  cov
 Usage   :  GENERAL::cov(\@array1, \@array2);
 Function:  Finds the the covariance between two equally sized arrays of numbers. Covariance
            is defined as cov(x, y) = mean((x - mean(x)) .* (y - mean(y))).
 Returns :  Covariance.
 Args    :  (1) First array.
            (2) Second array.
=cut

sub cov {
  my $arr1 = shift;
  my $arr2 = shift;
  GENERAL::requireArgs($arr1, $arr2);
  
  if (scalar(@$arr1) != scalar(@$arr2)) {
    GENERAL::error(sprintf("Arrays are of different size: %d and %d", scalar(@$arr1), scalar(@$arr2)));
  }
  
  my $cov = 0;
  my $mean1 = GENERAL::mean($arr1);
  my $mean2 = GENERAL::mean($arr2);
  for (my $i = 0; $i < scalar(@$arr1); $i++) {
    $cov += ($arr1->[$i] - $mean1)*($arr2->[$i] - $mean2);
  }
  $cov /= scalar(@$arr1);
  return $cov;
}


=head2

 Title   :  rmsd
 Usage   :  GENERAL::rmsd(\@array1, \@array2);
 Function:  Calculates the room mean square deviation between two equally sized arrays of numbers.
 Returns :  RMSD
 Args    :  (1) First array.
            (2) Second array.
=cut

sub rmsd {
  my $arr1 = shift;
  my $arr2 = shift;
  GENERAL::requireArgs($arr1, $arr2);
  
  if (scalar(@$arr1) != scalar(@$arr2)) {
    GENERAL::error(sprintf("Arrays are of different size: %d and %d", scalar(@$arr1), scalar(@$arr2)));
  }
  
  my $rmsd = 0;
  for (my $i = 0; $i < scalar(@$arr1); $i++) {
    $rmsd += ($arr1->[$i] - $arr2->[$i])**2;
  }
  return sqrt($rmsd/scalar(@$arr1));
}


=head2

 Title   :  mean
 Usage   :  GENERAL::mean(\@array1);
 Function:  Calculates the mean of an array of numbers.
 Returns :  Mean.
 Args    :  (1) Array.
=cut

sub mean {
  my $arr = shift;
  GENERAL::requireArgs($arr);
  return GENERAL::sum($arr)/scalar(@$arr);
}


=head2

 Title   :  median
 Usage   :  GENERAL::median(\@array1);
 Function:  Calculates the median of an array of numbers.
 Returns :  Median
 Args    :  (1) Array referrence.
=cut

sub median {
  my $arr = shift;
  GENERAL::requireArgs($arr);

  my @sarr = sort { $a <=> $b } @$arr;
  if (scalar(@sarr) % 2) {
    return $sarr[(scalar(@sarr)-1)/2];
  } else {
    return 0.5*($sarr[scalar(@sarr)/2 - 1] + $sarr[scalar(@sarr)/2]);
  }
}


=head2

 Title   :  mad
 Usage   :  GENERAL::mad(\@array1);
 Function:  Calculates the median absolute deviation.
 Returns :  Median absolute deviation
 Args    :  (1) Array referrence.
=cut

sub mad {
  my $arr = shift;
  GENERAL::requireArgs($arr);

  my $median = GENERAL::median($arr);
  my @adev = @$arr;
  for (my $i = 0; $i < scalar(@$arr); $i++) {
    $adev[$i] = abs($arr->[$i] - $median);
  }
  return GENERAL::median(\@adev);
}


=head2

 Title   :  stdev
 Usage   :  GENERAL::stdev(\@array1);
 Function:  Calculates the standard deviation of an array of numbers.
 Returns :  Standard deviation.
 Args    :  (1) Array.
=cut

sub stdev {
  my $arr = shift;
  GENERAL::requireArgs($arr);
  
  my $ss = 0;
  foreach my $v (@$arr) { $ss += $v**2; }
  my $var = $ss/scalar(@$arr) - (GENERAL::mean($arr))**2;
  return sqrt(GENERAL::max($var, 0));
}

=head2

 Title   :  sum
 Usage   :  my $sum = GENERAL::sum(\@array);
 Function:  Returns the sum of all elements in the array.
 Returns :  Returns the sum of all elements in the array.
 Args    :  (1) An array values.

=cut

sub sum {
  my @arr;
  if (ref($_[0])) {
    my $arr = shift;
    @arr = @$arr;
  } else {
    @arr = @_;
  }
  my $sum = 0;
  foreach (@arr) { $sum += $_; }
  return $sum;
}

=head2

 Title   :  prod
 Usage   :  my $prod = GENERAL::prod(\@array);
 Function:  Returns the product of all elements in the array.
 Returns :  Returns the product of all elements in the array.
 Args    :  (1) An array values.

=cut

sub prod {
  my @arr;
  if (ref($_[0])) {
    my $arr = shift;
    @arr = @$arr;
  } else {
    @arr = @_;
  }
  my $prod = 1;
  foreach (@arr) { $prod *= $_; }
  return $prod;
}



=head2

 Title   :  shuffle
 Usage   :  GENERAL::shuffle(\@array);
 Function:  Randomly permutes the entries of the array. Taken from the Perl cookbook (FisherYatesShuffle)
 Returns :  Nothing
 Args    :  (1) Reference to an array

=cut

sub shuffle {
  my $array = shift;
  return 0 if (scalar(@$array) == 0);
  my $i;
  for ($i = @$array; --$i; ) {
      my $j = int rand ($i+1);
      next if $i == $j;
      @$array[$i,$j] = @$array[$j,$i];
  }
}

=head2

 Title   :  informationContent
 Usage   :  GENERAL::informationContent(\@seqs, 20);
            GENERAL::informationContent(\@seqs, 20, 'pos', \@positions);
 Function:  Calculates the information content in the given sequence alignment.
 Returns :  Information content in bits.
 Args    :  (1) Reference to an array of sequences. Each array entry is a reference to an array of
                string, each string being a "letter" (either protein or DNA or whatever). All
                sequences must be of the same length (must have the same number of letters).
            (2) Total number of letters in the relevant alphabet (e.g. 20 for proteins, 4 for DNA, etc)
             *  Additional parameters are specified in name/value pairs. Possibilities are:
             ** 'pos' -- a list of a subset of sequence positions to calculate information content over
             ** 'lowCountCorrection' -- controls whether the low-cont correction is applied and what kind.
                 0 - means do not do anything. 1 - means apply the standard low-count correction
                 (ala Basharin, Basharin. G. P. (1959). Theory Probability Appl. 4, 333-336;
                  what WebLogo uses), which is default. 2 - means apply the new and improved
                 low-count correction SPECIFIC FOR PROTEINS (polynomial fit to empirical observation of
                 information content in random alignments with an alphabet of 20 characters).
             ** 'prevSaved' -- a data structure (hash reference) from a previous call to the function. This assumes
                that the new call is only different by one sequence, i.e. the last one. Will reuse the save
                information to make this a constant-time operation in the number of sequences and will update
                the data structure. The first time this is used, a reference to an empty hash should be passed
                to initialize the structure correctly.
            
=cut

sub informationContent {
  my $seqs = shift;
  my $S = shift;
  my $N = scalar(@$seqs);
  GENERAL::requireArgs($seqs, $S);
  GENERAL::assert(scalar(@_) % 2 == 0, "expected an even number of inputs to form name/value pairs");
  GENERAL::assert($N != 0, "empty alignment passed!");
  my $W = scalar(@{$seqs->[0]});
  my %opts = @_;
  foreach my $name (keys(%opts)) { GENERAL::assert($name =~ /^(pos|lowCountCorrection|prevSaved)$/ ? 1 : 0, "unknown option '$name'"); }
  if (!defined($opts{pos})) {
    my @pos; for (my $i = 0; $i < $W; $i++) { push(@pos, $i+1); }
    $opts{pos} = \@pos;
  }
  $opts{lowCountCorrection} = 1 if (!defined($opts{lowCountCorrection}));
  my ($P, $usePrev);
  if (defined($opts{prevSaved})) {
    $P = $opts{prevSaved};
    $usePrev = 1 if (scalar(keys(%$P)) > 0);
  } else {
    my %tmp; $P = \%tmp;
    $usePrev = 0;
  }
  if (!$usePrev) {
    foreach my $posi (@{$opts{pos}}) {
      my %tmp; $P->{$posi-1} = \%tmp;
    }
  }

  # compute probabilities of all observed letters at all required positions
  my $I = 0; my $Is = 0;
  foreach my $posi (@{$opts{pos}}) {
    my $i = $posi - 1;
    for (my $k = defined($opts{prevSaved}) ? scalar(@$seqs) - 1 : 0; $k < scalar(@$seqs); $k++) {
      my $seq = $seqs->[$k];
      GENERAL::assert(scalar(@$seq) > $i, "position $i is past the end of sequence '" . join(" ", @$seq) . "'");
      $P->{$i}->{$seq->[$i]} = 0 if (!defined($P->{$i}->{$seq->[$i]}));
      $P->{$i}->{$seq->[$i]}++;
    }
    my $Hmax = log($S)/log(2);
    my $H = 0;
    foreach my $let (keys(%{$P->{$i}})) {
      my $p = $P->{$i}->{$let} / scalar(@$seqs);
      $H -= $p*log($p)/log(2);
    }
    $I += $Hmax - $H;
    if ($opts{lowCountCorrection} == 1) {
      # standard low-count correction (ala Basharin)
      $I -= (log(exp(1))/log(2))*($S-1)/(2*$N);
    } elsif ($opts{lowCountCorrection} == 2) {
      # improved low-count correction (fit to empirical data for proteins; performs much better for very low counts)
      $I -= -0.9372/sqrt($N) + 23.8502/$N - 84.8337/($N**2) + 196.7665/($N**3) - 230.3894/($N**4) + 99.8655/($N**5);
      $Is += 1/(0.001045 * $N**2 + 0.08307 * $N + 3.3563);
    }
  }

  if (($opts{lowCountCorrection} == 2) && wantarray) {
    return ($I, $Is);
  }
  return $I;
}


=head2

 Title   :  getUnit
 Usage   :  my @uv = GENERAL::getUnit(\@v);
 Function:  Returns a unit vector in the direction of the given vector
 Returns :  
 Args    :  

=cut

sub getUnitVector {
  my $v = shift;
  my $L = 0;
  foreach my $e (@$v) { $L += $e**2; }
  $L = sqrt($L);
  my @v = @$v;
  for (my $i = 0; $i < scalar(@v); $i++) {
    $v[$i] = $v[$i]/$L;
  }
  return @v;
}

=head2

 Title   :  scale
 Usage   :  my @v = GENERAL::scale(\@v, $sf);
 Function:  scales the given vector by the given value
 Returns :  
 Args    :  

=cut

sub scale {
  my $v = shift;
  my $sf = shift;
  for (my $i = 0; $i < scalar(@$v); $i++) { $v->[$i] *= $sf; }
}

=head2

 Title   :  vecScale
 Usage   :  my @v = GENERAL::scale(\@v, $sf);
 Function:  scales the given vector by the given value
 Returns :  
 Args    :  

=cut

sub vecScale {
  my $v = shift;
  my $sf = shift;
  my @v = @$v;
  for (my $i = 0; $i < scalar(@v); $i++) { $v[$i] *= $sf; }
  return @v;
}

=head2

 Title   :  vecSum
 Usage   :  my @v = GENERAL::vecSum(\@v1, \@v2);
 Function:  vector sum
 Returns :  
 Args    :  

=cut

sub vecSum {
  my $v1 = shift;
  my $v2 = shift;
  GENERAL::assert(scalar(@$v1) == scalar(@$v2));
  my @v = @$v1;
  for (my $i = 0; $i < scalar(@$v2); $i++) { $v[$i] += $v2->[$i]; }
  return @v;
}

=head2

 Title   :  vecDiff
 Usage   :  my @v = GENERAL::vecDiff(\@v1, \@v2);
 Function:  vector difference
 Returns :  
 Args    :  

=cut

sub vecDiff {
  my $v1 = shift;
  my $v2 = shift;
  GENERAL::assert(scalar(@$v1) == scalar(@$v2));
  my @v = @$v1;
  for (my $i = 0; $i < scalar(@$v2); $i++) { $v[$i] -= $v2->[$i]; }
  return @v;
}

=head2

 Title   :  vecNorm
 Usage   :  my @v = GENERAL::vecNorm(\@v1);
 Function:  vector norm
 Returns :  
 Args    :  

=cut

sub vecNorm {
  my $v = shift;
  my $n = 0;
  for (my $i = 0; $i < scalar(@$v); $i++) { $n += ($v->[$i])**2; }
  return sqrt($n);
}

=head2

 Title   :  cross
 Usage   :  my @cp = GENERAL::cross(\@vec1, \@vec2);
 Function:  Returns the cross product between two 3D vectors.
 Returns :  
 Args    :  1. Vector 1 (reference to an array)
            2. Vector 2 (reference to an array)
            3. Reference flag. Will return an array reference if set.

=cut

sub cross {
  my $v1 = shift;
  my $v2 = shift;
  GENERAL::requireArgs($v1, $v2);
  my $ptr = shift;
  $ptr = 0 if (!defined($ptr));

  if ((scalar(@$v1) != 3) || (scalar(@$v2) != 3)) { GENERAL::error("Cross product implemented only for 3D vectors"); }
  my @cp;
  $cp[0] = $v1->[1]*$v2->[2] - $v1->[2]*$v2->[1];
  $cp[1] = $v1->[2]*$v2->[0] - $v1->[0]*$v2->[2];
  $cp[2] = $v1->[0]*$v2->[1] - $v1->[1]*$v2->[0];
  if ($ptr) { return \@cp; }
  else { return @cp; }
}


=head2

 Title   :  dot
 Usage   :  my $dp = GENERAL::dot(\@vec1, \@vec2);
 Function:  Returns the dot product between two vectors.
 Returns :  
 Args    :  1. Vector 1 (reference to an array)
            2. Vector 2 (reference to an array)

=cut

sub dot {
  my $v1 = shift;
  my $v2 = shift;
  GENERAL::requireArgs($v1, $v2);

  if (scalar(@$v1) != scalar(@$v2)) { GENERAL::error("Vectors must be of the same length!"); }
  my $dp = 0;
  for (my $i = 0; $i < scalar(@$v1); $i++) {
    $dp += $v1->[$i]*$v2->[$i];
  }

  return $dp;
}


=head2

 Title   :  mod
 Usage   :  my $dp = GENERAL::mod($n, $base);
 Function:  Returns the proper mod of the first number base second number
 Returns :
 Args    :  1. number to mod
            2. base of the mod

=cut

sub mod {
  return $_[0] - floor($_[0] / $_[1]) * $_[1];
}

=head2

 Title   :  arrayRef
 Usage   :  my $ref = GENERAL::arrayRef(1, 2, 3);
 Function:  Returns a reference to the array of parameters
 Returns :
 Args    :  1. parameter 1
            ....

=cut

sub arrayRef {
  my @arr = @_;
  return \@arr;
}


=head2

 Title   :  angleDiff
 Usage   :  my $dp = GENERAL::angleDiff($a, $b);
 Function:  Calculates the difference between two angles between -360 and 360. The difference is guaranteed to be between
            -180 and 180, where the negative sign means "clockwise" difference and positive "counter-clockwise".
 Returns :  angle 1 - angle 2
 Args    :  1. angle 1
            2. angle 2

=cut

sub angleDiff {
  my $a = shift;
  my $b = shift;
  GENERAL::requireArgs($a, $b);
#  GENERAL::assert((abs($a) <= 180) && (abs($b) <= 180), "Angles must be within [-180; 180], but are $a and $b!");

  my $da = mod((mod($a, 360) - mod($b, 360)), 360);
  $da -= 360 if ($da > 180);

  return $da;
}


=head2

 Title   :  angleMean
 Usage   :  my $dp = GENERAL::angleMean(@arr);
 Function:  Calculates the mean angle, on a circle.
 Returns :  Mean angle, in degrees
 Args    :  1. array of angles in degrees

=cut

sub angleMean {
  my @arr = @_;

  my $avcos = 0;
  my $avsin = 0;
  foreach my $a (@arr) {
    $avcos += cos($a*pi()/180);
    $avsin += sin($a*pi()/180);
  }
  $avcos /= scalar(@arr);
  $avsin /= scalar(@arr);
  my $n = sqrt($avcos**2 + $avsin**2);
  $avcos /= $n;
  $avsin /= $n;

  my $angle = acos($avcos);
  $angle = -$angle if ($avsin < 0);

  return $angle*180/pi();
}


=head2

 Title   :  linspace
 Usage   :  my @arr = GENERAL::linspace($beg, $end, $inc);
 Function:  Returns an array with values from $beg to $end in increments of $inc
            (equivalent to Matlab's linspace).
 Returns :  An array.
 Args    :  (1) Beginning value.
            (2) End value.
            (3) Increment.

=cut

sub linspace {
  my $beg = shift;
  my $end = shift;
  my $inc = shift;
  GENERAL::requireArgs($beg, $end, $inc);

  # Prevent infinite loops
  if (($inc == 0) || (($beg < $end) && ($inc < 0)) || (($beg > $end) && ($inc > 0))) {
    GENERAL::error("From $beg to $end with increments of $inc will cause an infinite loop!\n");
  }

  # Make the array
  my @arr;
  my $val = $beg;
  while ($val <= $end) {
    push(@arr, $val);
    $val += $inc;
  }

  return @arr;
}


=head2

 Title   :  ones
 Usage   :  my @arr = GENERAL::ones($N, $val);
 Function:  Returns an array filled with $val (default 1's) of specified length.
 Returns :  An array.
 Args    :  (1) Array length.
            (2) Optional: value to fill the array with
            (3) Optional: if set to 1, will return a reference to the created array. 0 by default.

=cut

sub ones {
  my $N = shift;
  GENERAL::requireArgs($N);
  my $val = 1;
  $val = shift if (scalar(@_) != 0);
  my $ptr = shift;

  if ($N < 0) { GENERAL::error("Array length specified as $N!\n"); }

  my @arr;
  $#arr = ($N-1);
  for (my $i = 0; $i < $N; $i++) { $arr[$i] = $val; }

  if (defined($ptr) && $ptr) {
    return \@arr;
  } else {
    return @arr;
  }
}


=head2

 Title   :  stdiff
 Usage   :  GENERAL::stdiff(\@array1, \@array2);
 Function:  Returns the difference between two sets (set 1 - set 2).
 Returns :  Array
 Args    :  (1) Set 1 (array).
            (2) Set 2 (array).
=cut

sub setdiff {
  my $set1 = shift;
  my $set2 = shift;
  GENERAL::requireArgs($set1, $set2);
  
  my (@diff, %h2);
  foreach my $e2 (@$set2) { $h2{$e2} = 1; }
  foreach my $e1 (@$set1) {
    if (!defined($h2{$e1})) { push(@diff, $e1); }
  }
  return @diff;
}

=head2

 Title   :  intersect
 Usage   :  GENERAL::intersect(\@array1, \@array2);
 Function:  Returns the intersection between two sets (set 1 - set 2).
 Returns :  Array
 Args    :  (1) Set 1 (array).
            (2) Set 2 (array).
=cut

sub intersect {
  my $set1 = shift;
  my $set2 = shift;
  GENERAL::requireArgs($set1, $set2);
  
  my (@diff, %h2);
  foreach my $e2 (@$set2) { $h2{$e2} = 1; }
  foreach my $e1 (@$set1) {
    if (defined($h2{$e1})) { push(@diff, $e1); }
  }
  return @diff;
}

=head2

 Title   :  union
 Usage   :  GENERAL::union(\@array1, \@array2);
 Function:  Returns the union of two sets.
 Returns :  Array
 Args    :  (1) Set 1 (array).
            (2) Set 2 (array).
=cut

sub union {
  my $set1 = shift;
  my $set2 = shift;
  GENERAL::requireArgs($set1, $set2);
  
  my %h;
  foreach my $e (@$set1) { $h{$e} = 1; }
  foreach my $e (@$set2) { $h{$e} = 1; }
  return keys(%h);
}


=head2   FUNCTION leastSquares3D

 Title   :  superimpose
 Usage   :  PDB::leastSquares3D(\@X, \@Y, \@Z)
 Function:  Finds the least-squares line, in 3D, that optimally fits the give set of points, in 3D.
 Returns :  An array of six values - first three are an origin of the best-fitting axis, and the next three are the direction (normalized vector)
 Args    :  (1) reference to an array of x coordinates
            (2) reference to an array of y coordinates
            {3} reference to an array of z coordinates
=cut

sub leastSquares3D {
  my $x = shift;
  my $y = shift;
  my $z = shift;
  GENERAL::requireArgs($x, $y, $z);
  GENERAL::assert((scalar(@$x) == scalar(@$y)) && (scalar(@$x) == scalar(@$z)), "coordinate vectors are of different length!");
  my $pi = 4*atan(1);

  my $Xm = 0; my $Ym = 0; my $Zm = 0;
  my $Sxx = 0; my $Syy = 0; my $Szz = 0; my $Sxy = 0; my $Sxz = 0; my $Syz = 0;
  for (my $i = 0; $i < scalar(@$x); $i++) {
    $Xm += $x->[$i];
    $Ym += $y->[$i];
    $Zm += $z->[$i];
    $Sxx += ($x->[$i])**2;
    $Syy += ($y->[$i])**2;
    $Szz += ($z->[$i])**2;
    $Sxy += $x->[$i]*$y->[$i];
    $Sxz += $x->[$i]*$z->[$i];
    $Syz += $y->[$i]*$z->[$i];
  }
  my $n = scalar(@$x);
  $Xm /= $n;
  $Ym /= $n;
  $Zm /= $n;

  $Sxx = $Sxx/$n - $Xm*$Xm;
  $Syy = $Syy/$n - $Ym*$Ym;
  $Szz = $Szz/$n - $Zm*$Zm;
  $Sxy = $Sxy/$n - $Xm*$Ym;
  $Sxz = $Sxz/$n - $Xm*$Zm;
  $Syz = $Syz/$n - $Ym*$Zm;

  my $c0 = $Sxy*$Sxz*($Szz-$Syy)+$Syz*($Sxy*$Sxy-$Sxz*$Sxz);
  my $c1 = $Sxz*$Syz*(2*$Sxx-$Syy-$Szz)+$Sxy*(2*$Syz*$Syz-$Sxy*$Sxy-$Sxz*$Sxz)+$Sxy*($Szz-$Sxx)*($Szz-$Syy);
  my $c2 = $Sxy*$Sxz*($Sxx+$Syy-2*$Szz)+$Syz*($Sxz*$Sxz+$Syz*$Syz-2*$Sxy*$Sxy)+$Syz*($Sxx-$Syy)*($Szz-$Sxx);
  my $c3 = $Sxz*$Syz*($Syy-$Sxx)+$Sxy*($Sxz*$Sxz-$Syz*$Syz);

  my $r = $c2/$c3;
  my $s = $c1/$c3;
  my $t = $c0/$c3;

  my $p = $s - $r*$r/3;
  my $q = 2*$r*$r*$r/27 - $r*$s/3 + $t;
  my $R = $q*$q/4 + $p*$p*$p/27;

  # go through all solutions and pick the best one
  my @sol;
  for (my $i = 0; $i < 3; $i++) {
    # choose solution and find a
    my $a;
    if ($R > 0) {
      $a = -$r/3 + (-$q/2+$R**(1/2))**(1/3) + (-$q/2-$R**(1/2))**(1/3);
    } else {
      my $rho = (-$p*$p*$p/27)**(1/2);
      my $phi = acos(-$q/2/$rho);

      $a = -$r/3 + 2 * $rho**(1/3) * cos(($phi + $i*2*$pi)/3);
    }

    # find b and compute the least-squares distance for this solution
    my $b = ( $a*($Szz-$Sxx)+(1-$a*$a)*$Sxz ) / ($Sxy+$a*$Syz) ;
    my $sum = 1 + $a*$a + $b*$b;
    my $u = ((1+$b*$b)*$Xm-$a*$b*$Ym+$a*$Zm)  / $sum;
    my $v = (-$a*$b*$Xm+(1+$a*$a)*$Ym+$b*$Zm) / $sum;
    my $w = ($a*$Xm+$b*$Ym+($a*$a+$b*$b)*$Zm) / $sum;
    my $dist = 0;
    for (my $i = 0; $i < $n; $i++) {
      my $xk = ((1+$b*$b)*$x->[$i] - $a*$b*$y->[$i] + $a*$z->[$i])/$sum;
      my $yk = (-$a*$b*$x->[$i] + (1+$a*$a)*$y->[$i] + $b*$z->[$i])/$sum;
      my $zk = ($a*$x->[$i] + $b*$y->[$i] + ($a*$a+$b*$b)*$z->[$i])/$sum;
      $dist += ($xk-$u)**2 + ($yk-$v)**2 + ($zk-$w)**2;
    }
    if ((scalar(@sol) == 0) || ($sol[-1] > $dist)) {
      @sol = ($a, $b, $u, $v, $w, $dist);
    }
    last if ($R > 0);
  }

  my ($a, $b, $u, $v, $w, $dist) = @sol;
  my $rl = sqrt($a**2 + $b**2 + 1);
  return ($Xm, $Ym, $Zm, $a/$rl, $b/$rl, -1/$rl);
}

=head2

 Title   :  irand
 Usage   :  my $random_integer = GENERAL::irand($beg, $end);
 Function:  Returns a random integer in the range $beg to $end (including the ends).
 Returns :  Returns a random integer in the range $beg to $end (including the ends).
 Args    :  (1) Beginning of the range (integer).
            (2) End of the range (integer).

=cut
sub irand {
  my $n = shift;
  my $m = shift;
  return int($n + ($m - $n +1)*rand());
}

=head2

 Title   :  exp_eval
 Usage   :  GENERAL::exp_eval($expr, $hash);
 Function:  Evaluates the string expression with values from the
            given hash table. For example, expression '(1 - s*num)**p'
            will be evaluated as '(1 - $hash->{s}*$hash->{num})**$hash->{p}'
 Returns :  The result of the evaluation.
 Args    :  (1) string expression
            (2) reference to hash table with values

=cut

sub expr_eval {
  my $oexpr = shift;
  my $hash = shift;
  GENERAL::requireArgs($oexpr, $hash);

  # Pad expression with spaces so that each variable name has a "previous" and "next" characters
  my $expr = " $oexpr ";
  foreach my $key (keys(%$hash)) {
    $expr =~ s/([^A-Za-z0-9_])$key([^A-Za-z0-9_])/$1$hash->{$key}$2/g;
  }
  my $ans = eval("return $expr");
  if ($@) {
    GENERAL::error("Could not evaluate expression \"$oexpr\", which parsed into \"$expr\"!");
  }
  return $ans;
}


=head2
 
 Title   :  LambertW
 Usage   :  my $W = GENERAL::LambertW($val);
 Function:  Calculates the value of the Lambert's W function (the solution to the
            equation x*exp(x) = $val). This algorithm was converted from a Matlab function
            found on MATLAB Central File Exchage written by Pascal Getreuer. The original
            M-file:
            
            function w = lambertw(b,z)
            %LAMBERTW  Lambert W-Function.
            %   W = LAMBERTW(Z) computes the principal value of the Lambert 
            %   W-Function, the inverse of Z = W*exp(W).  Z may be a 
            %   complex scalar or array.  For real Z, the result is real on
            %   the principal branch for Z >= -1/e.
            %
            %   W = LAMBERTW(B,Z) specifies which branch of the Lambert 
            %   W-Function to compute.  If Z is an array, B may either be an
            %   integer array of the same size as Z or an integer scalar.  
            %   If Z is a scalar, B may be an array of any size.
            %   
            %   The algorithm uses series approximations as initializations
            %   and Halley's method as developed in Corless, Gonnet, Hare,
            %   Jeffrey, Knuth, "On the Lambert W Function", Advances in
            %   Computational Mathematics, volume 5, 1996, pp. 329-359.
            
            % Pascal Getreuer 2005
            
            if nargin == 1
              z = b;
              b = 0;
            end
            
            w = (1 - 2*abs(b)).*sqrt(2*(exp(1)*z + 1)) - 1;
            j = find(abs(z + exp(-1)) > 1.5 - abs(b) | (b.*imag(z) > 0));
            tmp = log((z == 0) + z) + (i*2*pi)*b;
            w(j) = tmp(j) - log(tmp(j));
            
            for k = 1:16
              c1 = exp(w);
              c2 = w.*c1 - z;
              w1 = w + (w ~= 1);
              dw = c2./(c1.*w1 - 0.5*((w + 2).*c2./w1));
              w = w - dw;
              
              if abs(dw) < eps*(2+abs(w))
                  break;
              end      
            end
             
 Returns :  Value of Lambert's W function on the principal branch.
 Args    :  (1) value

=cut

sub LambertW {
  my $z = shift;
  my ($w, $c1, $c2, $w1, $j, $tmp, $dw);
  my $eps = 2.2*10**-16; # machine precision
    
  $w = sqrt(2*(exp(1)*$z + 1)) - 1;
  $j = (abs($z + exp(-1)) > 1.5) ? 1:0;
  if ($j) {
    $tmp = log(($z == 0) + $z);
    $w = $tmp - log($tmp);
  }
  for (my $k = 1; $k <= 16; $k++) {
    $c1 = exp($w);
    $c2 = $w*$c1 - $z;
    $w1 = $w + (($w != 1) ? 1 : 0);
    $dw = $c2/($c1*$w1 - 0.5*(($w + 2)*$c2/$w1));
    $w = $w - $dw;
    
    last if (abs($dw) < $eps*(2 + abs($w)));
  }
  
  return $w;
}


=head2

 Title   :  my_MPI_Barrier
 Usage   :  GENERAL::my_MPI_Barrier(MPI_COMM_WORLD);
            GENERAL::my_MPI_Barrier(MPI_COMM_WORLD, 5);
 Function:  Calls MPI barrier and then sleep for a specified amount of
            seconds. This is useful if one wants to make sure that NFS
            synchronizes before processes go on and start reading the
            same files.
 Returns :  nothing
 Args    :  (1) world reference
            (2) number of seconds to sleep after the barrier (default is 5).

=cut

sub my_MPI_Barrier {
  my $world = shift;
  my $stime = shift;
  # By default sleep for 0 second
  if (!defined($stime)) { $stime = 0; }
  
  MPI_Barrier($world);
  if ($stime > 0) {
    if ($PRANK_DEF == 0) {
      print "Waiting after MPI_Barrier ($stime seconds)...\n";
      flush STDOUT;
    }
    sleep($stime);
  }
}


=head2

 Title   :  readBinFile
 Usage   :  GENERAL::readBinFile($filename);
            GENERAL::readBinFile($filename, 'd');
 Function:  Reads in a binary file and returns an array of values.
 Returns :  An array of read values.
 Args    :  (1) file name
            (2) optional: packing format to pass to the unpack
                function. The default is 'd' (double precision
                floating point value).

=cut

sub readBinFile {
  my $file = shift;
  my $format = shift;
  
  if (!defined($file)) {
    GENERAL::error("Binary file name not specified");
  }

  if (!defined($format)) { $format = 'd'; }

  my $fh = GENERAL::GetInFH($file);
#  my $fh;
#  open($fh, "< $file") or die "Error in readBinFile: unable to open file $file\n";
  binmode($fh);

  my $block = "";
  my @unbuf = ();
  while (sysread($fh, $block, 1024)) {
    push(@unbuf, unpack("$format*", $block));
  }
  close($fh);

  return \@unbuf;
}

=head2

 Title   :  readFasta
 Usage   :  my $S = GENERAL::readFasta($filename);
 Function:  Returns a reference to an array of read sequences. Each is a hash
            with entries 'name' and 'seq', corresponding to the name and the
            sequence, respectively.
 Returns :  An array of sequence entries.
 Args    :  (1) fasta file name

=cut

sub readFasta {
  my $fname = shift;
  my (@S, $line);
  my $ifh = GENERAL::GetInFH($fname);
  my $name = "";
  GENERAL::skipTo($ifh, '^\s*>', \$name);
  while ($name ne "") {
    $name =~ s/^>//;
    my %seq; $seq{name} = $name; $seq{seq} = "";
    $name = "";
    while ($line = <$ifh>) {
      $line = GENERAL::Trim($line);
      if ($line =~ /^>/) { $name = $line; last; }
      $seq{seq} .= $line;
    }
    $seq{seqArr} = GENERAL::arrayRef(split("", $seq{seq}));
    push(@S, \%seq);
  }
  close($ifh);

  return \@S;
}


=head2

 Title   :  hashEntry
 Usage   :  my $val = GENERAL::hashEntry(\%hash, $key, $defaultValue);
 Function:  Returns the value at the given key of the given hash, unless that
            key does not exist, in which case returns the specified default value
            if given or undef otherwise.
 Returns :  value at the key or default value or undef.
 Args    :  (1) hash reference.
            (2) key.
            (3) optional: default value. If not given, will return undef on non-existing hash entry.

=cut

sub hashEntry {
  my $hash = shift;
  my $value = shift;
  my $default = shift;

  if (defined($hash->{$value})) {
    return $hash->{$value};
  } else {
    return $default;
  }
}

=head2

 Title   :  findin
 Usage   :  my $ind = GENERAL::findin($val, @arr);
 Function:  Returns the first index of the array at which the given value exists.
            If the value is not in the array, returns -1.
 Returns :  Index or -1.
 Args    :  (1) Value.
            (2) Array.

=cut

sub findin {
  my ($val, @ar) = @_;

  for (my $i = 0; $i < @ar; $i++) {
    if ($ar[$i] eq $val) {return $i;}
  }
  return -1;
}


=head2

 Title   :  findin_regexp
 Usage   :  my $ind = GENERAL::findin_regexp($regexp, @arr);
 Function:  Returns the first index of the array at which the given regular expression produces a match.
            If the value is not in the array, returns -1.
 Returns :  Index or -1.
 Args    :  (1) Regular expression string.
            (2) Array.

=cut

sub findin_regexp {
  my ($re, @ar) = @_;

  for (my $i = 0; $i < @ar; $i++) {
    if ($ar[$i] =~ /$re/) {return $i;}
  }
  return -1;
}


=head2

 Title   :  findinSorted
 Usage   :  my $ind = GENERAL::findinSorted($val, @arr);
 Function:  Returns the index of the array at which the given value exists.
            Assumes the array is sorted in ascending order.
            If the value is not in the array, returns -1.
 Returns :  Index or -1.
 Args    :  (1) Value.
            (2) Array.

=cut

sub findinSorted {
  my ($val, @ar) = @_;
  my $b = 0;
  my $e = scalar(@ar)-1;
  
  while ($b <= $e) {
    my $i = int(($b+$e)/2);
    if ($ar[$i] < $val) {
      $b = $i + 1;
    } elsif ($ar[$i] > $val) {
      $e = $i - 1;
    } else {
      return $i;
    }
  }
  
  return -1;
}

=head2

 Title   :  requireArgs
 Usage   :  GENERAL::requireArgs($arg1, $arg2, @arg3, \@arg4, $arg5, ...);
 Function:  Halts the code if any of its arguments are not defined. Prints
            an informative message (including the name of the caller function).
 Returns :  nothing
 Args    :  (1) A list of required arguments.

=cut

sub requireArgs {
  my $i = 1;
  foreach my $arg (@_) {
    if (!defined($arg)) {
      GENERAL::error("Error in function " . (caller(1))[3] . ": required argument number $i (out of " . scalar(@_) . ") not defined!!!");
    }
    $i++;
  }
}


=head2

 Title   :  GetFileInfo
 Usage   :  my $ans = GENERAL::GetFileInfo($file_name, $info_type);
 Function:  Returns useful information about the file.
 Returns :  Returns useful information about the file.
 Args    :  (1) File name.
            (2) Type of information to be returned.

=cut

sub GetFileInfo {
  my $fname = shift;
  my $infot = uc(shift);
  GENERAL::requireArgs($fname, $infot);

  if ($infot =~ /LINES/) {
#    my $str = `wc -l $fname`;
#    chomp($str);
#    if ($str =~ /No such file/) {
#      die "Error in GetFileInfo: invalid file name \"$fname\"!";
#    }
#    my @ar = split(" ", $str);
#    return $ar[0];
    my $fh = GENERAL::GetInFH($fname);
#    open FILE, "< $fname" or GENERAL::error("Could not open file $fname!");
    my $count = 0;
    $count += tr/\n/\n/ while sysread($fh, $_, 2 ** 16);
    close($fh);
    return $count;
  }

  GENERAL::error("invalid information type \"$infot\"!");
}

=head2

 Title   :  error
 Usage   :  GENERAL::error("error message");
 Function:  Prints out the specified error message along with the name of the caller
            function and halts the program.  Wraps the text of the message.
 Returns :  nothing
 Args    :  (1) Error message string.
            (2) Optional logical flag (default true). If true, the error message will
                also be logged in addition to being printed to stdout and stderr. If
                false, the error message will not be logged.

=cut

sub error {
  my $emsg = shift;
  my $log = shift;
  if (!defined($log)) { $log = 1; }

  my $fc; my $i=0;
  while (1) {
    if (!defined(caller($i))) { last; }
    if (defined(caller($i+1))) {
      $fc .= (caller($i))[3] . " (line " . (caller($i))[2] . ") <= ";
    } else {
      $fc .= (caller($i))[3] . " (line " . (caller($i))[2] . ")";
    }
    $i++;
  }
  my $msg;
  if (defined((caller(1))[3])) {
    if ($NPROC_DEF > 1) {
      $msg = "Error in " . (caller(1))[3] . " (process $PRANK_DEF):\n" . $emsg . "\n" . "Trace: $fc\n";
    } else {
      $msg = "Error in " . (caller(1))[3] . ":\n" . $emsg . "\n" . "Trace: $fc\n";
    }
  } else {
    if ($NPROC_DEF > 1) {
      $msg = "Error in main (process $PRANK_DEF):\n" . $emsg . "\n" . "Trace: $fc\n";
    } else {
      $msg = "Error in main:\n" . $emsg . "\n" . "Trace: $fc\n";
    }
  }
  if (defined($!) && !($! eq "")) {
    $msg = $msg . "Error message = $!\n";
  }
  $msg .= "Directory = " . GENERAL::GetDir() . "\n";
  $msg .= "Machine = " . GENERAL::GetMachine() . "\n";
  print(wrap("", "", $msg));
  if ($log) {
    GENERAL::GLog(wrap("", "", $msg));
  }
#  if (defined($NPROC_DEF) && ($NPROC_DEF > 1)) {
#    MPI_Finalize();
#  }
  die(wrap("", "", $msg));
}


=head2

 Title   :  assert
 Usage   :  GENERAL::assert($a == 1, "a is not equal to 1!")
 Function:  Checks if the specified condition is correct. If not, prints an error
            using the specified error message.
 Returns :  nothing
 Args    :  (1) Logical statement
            (2) Optional: error message.

=cut

sub assert {
  my $cond = shift;
  GENERAL::requireArgs($cond);
  my $emes = shift;
  if (!defined($emes)) { $emes = "Assertion failed!"; }

  if (!defined($cond) || !$cond) { GENERAL::error($emes); }
  return $cond;
}



=head2

 Title   :  warning
 Usage   :  GENERAL::warning("error message");
 Function:  Throws a warning with the specified message along with the name of the caller
            function. Wraps the text of the message. Also, writes the warning to the global
            log file.
 Returns :  nothing
 Args    :  (1) Warning message string.
            (2) Optional logical flag (default true). If true, the warning message will
                also be logged in addition to being printed to stdout and stderr. If
                false, the warning message will not be logged.

=cut

sub warning {
  my $wmsg = shift;
  my $log = shift;
  if (!defined($log)) { $log = 1; }

  my $fc; my $i=0;
  while (1) {
    if (!defined(caller($i))) { last; }
    if (defined(caller($i+1))) {
      $fc .= (caller($i))[3] . " <= ";
    } else {
      $fc .= (caller($i))[3];
    }
    $i++;
  }
  my $msg;
  if (defined((caller(1))[3])) {
    if ($NPROC_DEF > 1) {
      $msg = "Warning in " . (caller(1))[3] . " (process $PRANK_DEF):\n" . $wmsg . "\n" . "Trace: $fc\n";
    } else {
      $msg = "Warning in " . (caller(1))[3] . ":\n" . $wmsg . "\n" . "Trace: $fc\n";
    }
  } else {
    if ($NPROC_DEF > 1) {
      $msg = "Warning in main (process $PRANK_DEF):\n" . $wmsg . "\n" . "Trace: $fc\n";
    } else {
      $msg = "Warning in main:\n" . $wmsg . "\n" . "Trace: $fc\n";
    }
  }
  print(wrap("", "", $msg));
  if ($log) {
    GENERAL::GLog(wrap("", "", $msg));
  }
}


=head2

 Title   :  splitTasks
 Usage   :  my $a = GENERAL::splitTasks($number_of_tasks, $number_of_processes);
 Function:  Splits the specified number of tasks into the specified number of
            processes trying to keep things as balanced as possible. It is guaranteed
            that the difference in number of tasks assigned to any pair of processes
            is at most 1.
 Returns :  A reference to an array of hashes containing the first and last task to
            be performed by any process. So $a->[$i]->{beg} is the index of the first
            task to be performed by process with rank $i and $a->[$i]->{end} is the
            last task to be performed by that process. Process ranks start with 0 and
            task indices start with 1. Note, if ($a->[$i]->{beg} > $a->[$i]->{end}) the
            process is assigned 0 tasks. In addition, populates $a->[$i]->{num} with
            the total number of tasks to be performed by process with rank $i.
 Args    :  (1) Number of tasks.
            (2) Number of processes.
=cut

sub splitTasks {
  my $t = shift;
  my $p = shift;
  GENERAL::requireArgs($t, $p);

  my @a;
  if ($t >= $p) {
    my $n = ($p == 0) ? 0 : int($t/$p);
    my $rem = $t - $n*$p;
    my $k = 1;
    # So the first $rem processes do $n+1 tasks
    for (my $i=0; $i < $rem; $i++) {
      $a[$i] = {};
      $a[$i]->{beg} = $k;
      $a[$i]->{end} = $k+$n;
      $k += $n+1;
    }
    # While the remaining ones do $n tasks
    for (my $i = $rem; $i < $p; $i++) {
      $a[$i] = {};
      $a[$i]->{beg} = $k;
      $a[$i]->{end} = $k+$n-1;
      $k += $n;
    }
  } else {
    my $k = 1;
    # So the first $t processes do 1 task
    for (my $i=0; $i < $t; $i++) {
      $a[$i] = {};
      $a[$i]->{beg} = $k;
      $a[$i]->{end} = $k;
      $k += 1;
    }
    # While the remaining ones do nothing
    for (my $i = $t; $i < $p; $i++) {
      $a[$i] = {};
      $a[$i]->{beg} = $k;
      $a[$i]->{end} = $k-1;
    }
  }

  for (my $i = 0; $i < scalar(@a); $i++) {
    $a[$i]->{num} = $a[$i]->{end} - $a[$i]->{beg} + 1;
  }
  return \@a;
}


=head2

 Title   :  splitTasksNonContig
 Usage   :  my $a = GENERAL::splitTasksNonContig($number_of_tasks, $number_of_processes, \@aray_of_complete_tasks);
 Function:  Splits the specified number of tasks into the specified number of
            processes trying to keep things as balanced as possible. It is guaranteed
            that the difference in number of tasks assigned to any pair of processes
            is at most 1. The difference between this function and splitTasks is that this
            one can handle non-contiguous splitting as well as takes into accounts that
            some tasks may already be finished.
 Returns :  A reference to an array of arrays containing the task indices for each process to perform.
 Args    :  (1) Number of tasks.
            (2) Number of processes.
            (3) Optional: an array of taks indices for complete tasks.
=cut

sub splitTasksNonContig {
  my $nt = shift;
  my $p = shift;
  GENERAL::requireArgs($nt, $p);
  my $comp = shift;
  my @comp = ();
  # ascending sort
  if (defined($comp)) { @comp = sort {$a <=> $b } @$comp; }
  
  if (defined($comp) && (($comp[scalar(@comp)-1] > $nt) || ($comp[0] < 1))) {
    GENERAL::error("Array of completed tasks has task indices out of bounds - should be from 1 to $nt, but is from $comp[0] to ".$comp[scalar(@comp)-1].".");
  }
  if (defined($comp) && (scalar(@comp) > $nt)) {
    GENERAL::error(sprintf("Array of completed tasks is longer than the number of tasks (%d versus $nt)!", scalar(@comp)));
  }

  my $t; # number of tasks left to complete
  if (defined($comp)) { $t = $nt - scalar(@$comp); }
  else { $t = $nt; }
  my @tlist;
  if ($t >= $p) {
    my $n = int($t/$p);
    my $rem = $t - $n*$p;
    my $k = 1;
    # So the first $rem processes do $n+1 tasks
    my $ti = 1; my $ci = 0;
    for (my $i = 0; $i < $rem; $i++) {
      $tlist[$i] = ();
      # proc $i does $n+1 tasks
      for (my $ii = 1; $ii <= $n+1; $ii++) {
        while (($ci < scalar(@comp)) && ($ti == $comp[$ci])) { $ti++; $ci++; }
        push(@{$tlist[$i]}, $ti);
        $ti++;
      }
    }
    # While the remaining ones do $n tasks
    for (my $i = $rem; $i < $p; $i++) {
      $tlist[$i] = ();
      # proc $i does $n tasks
      for (my $ii = 1; $ii <= $n; $ii++) {
        while (($ci < scalar(@comp)) && ($ti == $comp[$ci])) { $ti++; $ci++; }
        push(@{$tlist[$i]}, $ti);
        $ti++;
      }
    }
    return \@tlist;
  } else {
    my $k = 1;
    # So the first $t processes do 1 task
    my $ti = 1; my $ci = 0;
    for (my $i = 0; $i < $t; $i++) {
      $tlist[$i] = ();
      # proc $i does 1 task
      while (($ci < scalar(@comp)) && ($ti == $comp[$ci])) { $ti++; $ci++; }
      push(@tlist, $ti);
      $ti++;
    }
    # While the remaining ones do nothing
    for (my $i = $t; $i < $p; $i++) {
      $tlist[$i] = ();
    }
    return \@tlist;
  }
}


=head2

 Title   :  limitAccess
 Usage   :  GENERAL::limitAccess($user_regexp, $message);
 Function:  Limits access by allowing only certain users (defined by the regular
            expression to be applied to the user name) to continue, while causes
            a halt for all others.
 Returns :  Nothing.
 Args    :  1. Regular expression defining user names, for which access is allowed.
            2. Optional: message to be printed before halting when access is denied.

=cut
sub limitAccess {
  my $uregexp = shift;
  my $mes = shift;

  my $user = GENERAL::GetUser();
  if ($user !~ /$uregexp/) {
    if (!defined($mes)) { $mes = "Access to code denied for user '$user'"; }
    GENERAL::error($mes);
  }
}

#-----------------------------------------------------------------------------------------



# function: MakeDir
# description: Give directory, recursively build up directory
# usage:    &MakeDir(/home/people/jchen/test/blah)

sub MakeDir {
    my $dir=shift;
    my @f=split(/\/+/,$dir);

    my $t="";
    while (@f) { 
        $t.=(shift @f)."/";
        system "mkdir $t 2>/dev/null" if (!-d $t);
    }
}

# function: GetSubDirList
# description: Get directory list in curreent directory
# usage:   &GENERAL::GetSubDirList()


sub GetSubDirList {
     my @filelists = `ls -l`;
     chomp @filelists;
     my @dirlists = grep { /^d/ } @filelists;
     
     my @results;
     foreach (@dirlists) {
         chomp($_);
         my @line=split(" ", $_);
         push (@results, $line[8]);
     }
     return \@results;
}


=head2

 Title   :  fileCheck
 Usage   :  my $bool = GENERAL::fileCheck(file, $operator_string);
 Function:  Given a file name applies the specified perl operator to it and returns the
            result. Optionally, also looks for the compressed version of the file if 
            the operator returns false with the file itself.
 Returns :  0 or 1.
 Args    :  1. File name.
            2. Operator string (such as "e" or "r")
            2. Optional: logical flag. If set, will look for the compressed version in case
               the test on the file itself returns false. Otherwise, will not look for the 
               compressed version. The defaul is to look.

=cut

sub fileCheck {
  my $fname = shift;
  my $op = shift;
  my $uncom = shift;
  if (!defined($uncom)) { $uncom = 1; }

  my $r = eval("if (-$op \"$fname\") { return 1; } else { return 0; }");
  if (($r == 0) && $uncom) {
    $r = eval("if ((-$op \"$fname.bz2\") || (-$op \"$fname.gz\")) { return 1; } else { return 0; }");
  }
  return $r;
}


=head2

 Title   :  GetInFH
 Usage   :  my $fh = GENERAL::GetInFH($file);
 Function:  Given a file name returns a read handle to it. If the appropriate flag is given,
            will search for a compressed form of the file (in case the specified name does
            not exist) and return a handle to an uncompression stream.
 Returns :  File handle.
 Args    :  1. File name.
            2. Optional: logical flag. If set, will look for the compressed version in case
               the file does not exist. Otherwise, will not look for the compressed version.
               The defaul is to look.

=cut

sub GetInFH {
  my $fname = shift;
  my $uncom = shift;
  if (!defined($uncom)) { $uncom = 1; }

  if (!defined $fname || $fname eq "") { GENERAL::error("Undefined or empty file name passed!"); }

  if (ref $fname) {
    return $fname;
  }
  
  my $newfile;
  if (! -r $fname) {
    if (!$uncom) {
      GENERAL::error("File \"$fname\" not openable");
    }
    if (-r "$fname.bz2") {
      $newfile = new IO::File("bzip2 -dc $fname.bz2 |");
      if (!defined($newfile)) { GENERAL::error("Open on $fname.bz2 failed"); }
    } elsif (-r "$fname.gz") {
      $newfile = new IO::File("gzip -dc $fname.gz |");
      if (!defined($newfile)) { GENERAL::error("Open on $fname.gz failed"); }
    } else {
      GENERAL::error("Could not open file '$fname' or its compressed version");
    }
  } else {
    $newfile = new IO::File($fname);
    if (!defined($newfile)) { GENERAL::error("Open on $fname failed"); }
  }

  return $newfile;
}


=head2

 Title   :  GetOutFH
 Usage   :  my $fh = GENERAL::GetOutFH($file);
 Function:  Given a file name returns a write handle to it. If the appropriate flag is given,
            will pipe the output through a compressor and write the compressed version of the
            file.
 Returns :  File handle.
 Args    :  1. File name.
            2. Optional: compression switch. If "bzip2", bzip2 will be used for compression.
               If "gzip", gzip will be used for compression. Actual file will be written to
               file.bz2 or file.gz in these cases.

=cut

sub GetOutFH {
  my $fname = shift;
  if (!defined $fname || $fname eq "") { die "Error in GetOutFH: undefined or empty file name passed!\n"; }
  my $comp = shift;
  if (!defined($comp)) { $comp = ""; }
  
  if (ref $fname) {
    return $fname;
  } else {
    my $newfile;
    if (uc($comp) =~ /BZIP2/) {
      $newfile = new IO::File("| bzip2 -zcf > $fname.bz2");
      if (!defined($newfile)) { GENERAL::error("Open for '| bzip2 -zcf > $fname.bz2' failed"); }
    } elsif (uc($comp) =~ /GZIP/) {
      $newfile = new IO::File("| gzip -cf > $fname.gz");
      if (!defined($newfile)) { GENERAL::error("Open for '| gzip -cf > $fname.gz' failed"); }
    } else {
      $newfile = new IO::File(">$fname");
      if (!defined($newfile)) { GENERAL::error("Open for '> $fname' failed"); }
    }
    return $newfile;
  }
}

=head2

 Title   :  GetFH
 Usage   :  my $fh = GENERAL::GetOutFH($file, $mode);
 Function:  Given a file name and a mode for opening returns a handle to it.
 Returns :  File handle.
 Args    :  1. File name.
            2. Opening mode.

=cut

sub GetFH {
  my $fname = shift;
  my $mode = shift;
  GENERAL::requireArgs($fname, $mode);
  
  if (ref $fname) {
    return $fname;
  } else {
    my $newfile = new IO::File $fname, $mode;
    if (!defined($newfile)) {
      GENERAL::error("Failed to open file '$fname' in mode '$mode'");
    }
    return $newfile;
  }
}


# function: GetAppendOutFH
# description: Given a file name, generates file handle to write to the end of that file
#              If a reference is passed, returns the reference and does nothing.
# usage: &GetInFH("test.out")
sub GetAppendOutFH {
    my $fname = shift;
      
    if (!defined $fname || $fname eq "") { GENERAL::error("undefined or empty file name passed!"); }

    if (ref $fname) {
      return $fname;
    } else {
      my $outputfh;
      $outputfh = IO::File->new(">> $fname") || GENERAL::error("can't create file $fname");
      return $outputfh;
    }
}    
    
    
=head2

 Title   :  GetInDH
 Usage   :  my $fd = GENERAL::GetInDH($file);
 Function:  Given a directory name returns a read handle to it.
 Returns :  Directory handle.
 Args    :  1. Directory name.

=cut

sub GetInDH {
  my $dname = shift;

  if (!defined $dname || $dname eq "") { GENERAL::error("Undefined or empty directory name passed!"); }

  if (ref $dname) {
    return $dname;
  }
  
  my $dh;
  opendir($dh, $dname) || GENERAL::error("Could not open directory $dname!");

  return $dh;
}


# function: createFileSpace
# description: creates the file space (directory structure) and assigns
#              common file names, so that they can be used later in the code.
#              All the information is stored in the hash table FS_DEF,
#              the pointer to which is the first argument.
#              Appropriate logging is done. Note, since one may have different
#              file spaces for different applications (or different modes), the
#              functions expects an integer as the second parameter to tell it which file space
#              to create. Also, it is possible to specify which files are to be
#              removed after use and which should stay. This way one may have a
#              debug file space.
# inputs:      0. pointer to the hash table in which to store the file space info
#              1. integer specifying which file space to create
#              2. flag specifying whether the file space is to be created as opposed
#                 to just initializing the file space hash table with the right paths.
#                 This parameter is optional and assumed to be 1 by default.
#              3. optional: home directory. If not specified, current directory will be used
# usage:
sub setFileSpace {
  my $FS_DEF = shift;
  my $c = shift;
  my $create = shift;
  if (!defined($create)) { $create = 1; }
  my $home = shift;

  if ($c == 1) {
    # Directory names
    if (!defined($home)) {
      $home = `pwd`; chomp($home);
    }
    $FS_DEF->{dhome} = $home;
    $FS_DEF->{designd} = "$FS_DEF->{dhome}/design";
    $FS_DEF->{solsd} = "$FS_DEF->{designd}/sols";
    $FS_DEF->{drots} = "$FS_DEF->{dhome}/rots";
    $FS_DEF->{detable} = "$FS_DEF->{designd}/energytable";
    $FS_DEF->{dself} = "$FS_DEF->{dhome}/self";
    $FS_DEF->{dselfeef} = "$FS_DEF->{dself}/eef";
    $FS_DEF->{dselfmulte} = "$FS_DEF->{dself}/multe";
    $FS_DEF->{dunfold} = "$FS_DEF->{dhome}/unfold";
    $FS_DEF->{dpair} = "$FS_DEF->{dhome}/pair";
    $FS_DEF->{dpair_multe} = "$FS_DEF->{dpair}/multe";
    $FS_DEF->{dpair_eef} = "$FS_DEF->{dpair}/eef";

    # File names
    $FS_DEF->{rot_statesf} = "$FS_DEF->{designd}/rot_states.self";
    $FS_DEF->{rotinpf} = "$FS_DEF->{dhome}/rot.inp";
    $FS_DEF->{econf} = "$FS_DEF->{dhome}/energy.conf";
    $FS_DEF->{rotinp1f} = "$FS_DEF->{designd}/rot_input";
    $FS_DEF->{rotinp2f} = "$FS_DEF->{designd}/rot_input2";
    $FS_DEF->{rotinp3f} = "$FS_DEF->{designd}/rot_input3";
    $FS_DEF->{aainputf} = "$FS_DEF->{designd}/aa_input";

#    $FS_DEF->{unfold_eef_tabf} = "$FS_DEF->{dunfold}/%%%_eef1.tab";
#    $FS_DEF->{unfold_multe_tabf} = "$FS_DEF->{dunfold}/%%%_multe.tab";
    $FS_DEF->{unfold_eef_tabf} = "$FS_DEF->{dunfold}/eef1_unfold.tab";
    $FS_DEF->{unfold_multe_tabf} = "$FS_DEF->{dunfold}/multe_unfold.tab";
#    $FS_DEF->{unfold_SA_tabf} = "$FS_DEF->{dunfold}/%%%_sa.tab";
#    $FS_DEF->{unfold_GB_tabf} = "$FS_DEF->{dunfold}/%%%_gb.tab";
    $FS_DEF->{unfold_SA_tabf} = "$FS_DEF->{dunfold}/sa_unfold.tab";
    $FS_DEF->{unfold_GB_tabf} = "$FS_DEF->{dunfold}/gb_unfold.tab";

    $FS_DEF->{self_eef_tabf} = "$FS_DEF->{dselfeef}/self_eef1.tab";
    $FS_DEF->{self_multe_tabf} = "$FS_DEF->{dselfmulte}/self_multe.tab";
    $FS_DEF->{self_SA_tabf} = "$FS_DEF->{dself}/self_sa.tab";
    $FS_DEF->{self_GB_tabf} = "$FS_DEF->{dself}/self_gb.tab";
    $FS_DEF->{born_radf} = "$FS_DEF->{dself}/born_radii.dat";

    $FS_DEF->{self_tabf} = "$FS_DEF->{detable}/SelfEnergyTable";
    $FS_DEF->{unfold_tabf} = "$FS_DEF->{detable}/UnfoldEnergyTable";
    $FS_DEF->{pair_tabf} = "$FS_DEF->{detable}/PairEnergyTable";

    $FS_DEF->{unf_binf} = "$FS_DEF->{detable}/energy_uf.bin";
    $FS_DEF->{ener_tot_binf} = "$FS_DEF->{detable}/energy_f_minus_uf.bin";
    $FS_DEF->{ener_nounf_binf} = "$FS_DEF->{detable}/energy_f.bin";

    $FS_DEF->{unf_logf} = "$FS_DEF->{detable}/energy_uf.log";
    $FS_DEF->{ener_tot_logf} = "$FS_DEF->{detable}/energy_f_minus_uf.log";
    $FS_DEF->{ener_nounf_logf} = "$FS_DEF->{detable}/energy_f.log";

    $FS_DEF->{rot_pdbf} = "$FS_DEF->{drots}/rot%%-%.pdb";
    $FS_DEF->{rot_crdf} = "$FS_DEF->{drots}/rot%%-%.crd";
    $FS_DEF->{rot_outf} = "$FS_DEF->{drots}/rot%%-%.out";
    $FS_DEF->{rot_ipdbf} = "$FS_DEF->{drots}/rot%%-%.%%%.pdb";
    $FS_DEF->{rot_icrdf} = "$FS_DEF->{drots}/rot%%-%.%%%.crd";
    $FS_DEF->{rot_ioutf} = "$FS_DEF->{drots}/rot%%-%.%%%.out";

#    $FS_DEF->{pair_eef_tabf} = "$FS_DEF->{dpair_eef}/paireef%.tab";
#    $FS_DEF->{pair_multe_tabf} = "$FS_DEF->{dpair_multe}/pairmulte%.tab";
    $FS_DEF->{pair_eef_tot_tabf} = "$FS_DEF->{dpair_eef}/paireef.tab";
    $FS_DEF->{pair_multe_tot_tabf} = "$FS_DEF->{dpair_multe}/pairmulte.tab";
    $FS_DEF->{pair_epic_tabf} = "$FS_DEF->{dpair}/pairepic.tab";
    $FS_DEF->{pair_SA_tabf} = "$FS_DEF->{dpair}/pair_sa.tab";
    $FS_DEF->{pair_GB_tabf} = "$FS_DEF->{dpair}/pair_gb.tab";
    $FS_DEF->{stat_pot_tabf} = "$FS_DEF->{detable}/statpot.tab";
    $FS_DEF->{eef_solv_tabf} = "$FS_DEF->{dself}/solv.tab";

    $FS_DEF->{term_trackerf} = "$FS_DEF->{dhome}/.energyterm.track";
    $FS_DEF->{self_ctrlf} = "$FS_DEF->{detable}/.selftotal.track";
    $FS_DEF->{pair_ctrlf} = "$FS_DEF->{detable}/.pairtotal.track";
    $FS_DEF->{unf_ctrlf} = "$FS_DEF->{detable}/.unfoldtotal.track";

    # Flags
    $FS_DEF->{del_charmm_inp} = 1; # Do we delete charmm input files?
    $FS_DEF->{del_pair_multe_tabf} = 1;
    $FS_DEF->{del_pair_eef_tabf} = 1;

    if ($create) {
      # Create directory structure
      cmkdir($FS_DEF->{designd});
      cmkdir($FS_DEF->{solsd});
      cmkdir($FS_DEF->{drots});
      cmkdir($FS_DEF->{detable});
      cmkdir($FS_DEF->{dself});
      cmkdir($FS_DEF->{dselfeef});
      cmkdir($FS_DEF->{dselfmulte});
      cmkdir($FS_DEF->{dunfold});
      cmkdir($FS_DEF->{dpair});
      cmkdir($FS_DEF->{dpair_multe});
      cmkdir($FS_DEF->{dpair_eef});
    }
  } else {
    die "Error in createFileSpace: I do not know how to create file space $c\n";
  }
}


=head2

 Title   : file2array
 Usage   : GENERAL::file2array($filanem)
 Function: Returns the contents of a file as an array of lines.
 Returns : A reference to an array of strings.
 Args    : 1. File name or file handle.
           2. Optional: first line index (discard all lines before this one)
           3. Optional: last line index (discard all lines after this one)
=cut

sub file2array {
  my $fh;
  my $cl = 0;
  if (ref($_[0])) {
    $fh = shift;
  } else {
    $fh = GENERAL::GetInFH(shift);
    $cl = 1;
  }
  my $start = shift;
  my $end = shift;
  
  my @lines = <$fh>;
  chomp(@lines);
  
  splice (@lines, 0, $start) if (defined $start);
  splice (@lines, -$end) if (defined $end);
  
  if ($cl) { close($fh); }
  
  return \@lines;
}


=head2

 Title   : skipTo
 Usage   : my $suc = GENERAL::skipTo($fh, '^\s*BEGIN')
 Function: Forwards the file handle until a line matching the given regular expression is encountered.
           Returns true if such a line is found, returns false if end of file is reached without finding such a line.
 Returns : Returns true if such a line is found, returns false if end of file is reached without finding such a line.
 Args    : 1. File handle
           2. Regular expression to skip to
           3. optional: pointer to a string that will be set to the matching line read
           4. optional: pointer to a string that will be filled with skipped lines
=cut

sub skipTo {
  my $fh = shift;
  my $re = shift;
  GENERAL::requireArgs($fh, $re);
  my $linep = shift;
  my $skipped = shift;
  $$skipped = "" if (defined($skipped));

  while (<$fh>) {
    if (/$re/) {
      $$linep = $_;
      return 1;
    } elsif (defined($skipped)) {
      $$skipped .= $_;
    }
  }
  return 0;
}

# function: Compress
# description: compress a file given file name
# usage: &Compress("test.out");

sub Compress {
    my $f=shift;
    system "gzip -f -$f" if (-r $f);
}


# function: Remove
# description: delete  a file given file name
# usage: &Remove("test.out");

sub Remove {
    my $f=shift;
    unlink $f if (-r $f);
}

sub Removes {
    my @f=@_;
    foreach my $i (@f) {
         unlink $i if (-r $i);
    }
}


# function: SetLogFile
# description: sets a log file for debug messages
# usage: &SetLogFile("test.log");

sub SetLogFile {
  my $logfile=&GetOutFH(shift);
  $logfile->autoflush(1);
}

# function: File2String
# description: This function takes a file handle and spit out the file content
# as a big string.
# usage: File2String($file_handle)

sub File2String {

    my $fh = shift;
    local $/ = undef;
    my $bigstring = <$fh>;
    return $bigstring;

}



# function: ChopArray
# description: this function takes array reference and number of pieces and size as input
# it return the reference to array after chopping
# usage: ChopArray($array, 10, 4); it will chop the araay into 10 pieces
#        each pieces with size of 4.


sub ChopArray {
    my $arf=shift;
    my $pieces=shift;
    my $psize=shift;

    if ( ! (defined $pieces) || ( $pieces == 0)) {
    print " You cann't chop array into 0 elements! \n";
    return;
    }
    
    my @chopped;
    for (my $j=0; $j< $pieces; $j++) {
    push (@chopped, ShiftArray($arf, $psize)); 
    }

    #now return the array of references pointed to chopped arrays
    return @chopped;
}


# function: GetTime 
# description: give the local time: by Month, Day, Hour and Min
# usage:  &GetTime();
sub GetTime {
    my $mode = shift;
    my ($week, $mon, $day, $hours, $zone, $year);
#    $! = "";
    my $dates = `date`;
    ($week, $mon, $day, $hours, $zone, $year) = split(" ", $dates);
    my ($hour, $min, $second) = split (/:/, $hours);
    my $st = "$mon"."-"."$day"."-"."$year";
    if (defined($mode) && ($mode =~ /1/)) {
      $st = "$st-$hour".":"."$min";
    } else {
      $st = $st."  "."$hour".":"."$min";
    }
    if (defined($mode) && ($mode =~ /sec/)) {
      $st = $st . ":$second";
    }
    return $st;
} 


# Author: Gevorg Grigoryan
# Date: 09/08/03
# function: GetUser 
# description: returns the current username
sub GetUser {
    my $user = `whoami`;
    chomp $user;
    return $user;
} 


# Author: Gevorg Grigoryan
# Date: 09/08/03
# function: GetMachine
# description: returns the current username
sub GetMachine {
#  $! = "";
  my $machine = `hostname`;
  chomp $machine;
  return $machine;
}

=head2

 Title   : getArch
 Usage   : GENERAL::getArch()
 Function: Returns a string representing the architecture of the machine.
 Returns : A string representing the architecture of the machine.
 Args    : None

=cut

sub getArch {
  my $arch = `uname -m`;
  return GENERAL::Trim($arch);
}

# Author: Gevorg Grigoryan
# Date: 12/12/03
# function: GetDir
# description: returns the current username
sub GetDir {
#  $! = "";
  my $dir = `pwd`;
  chomp $dir;
  return $dir;
}


# Author: Gevorg Grigoryan
# Date: 09/08/03
# function: listFiles(filename_pattern)
# description: returns the result of doing 'ls' on the specified pattern as
#              an array of file names. If no matches are found, an empty array
#              is returned.
sub listFiles {
  my $filename = shift;

  my $files = `ls -a`;
  my @files = split(" ", $files);
  my @matches = ();
  foreach my $file (@files) {
    if ($file =~ m/$filename/) {
      push(@matches, $file);
    }
  }
  return @matches;
} 

# Author: Gevorg Grigoryan
# Date: 09/08/03
# function: Log(string, mode)
# description: Logging function. Appends things to the file name specified
# in environmental variable LOG_DEF. If the environmental variable does not exist,
# does nothing. The second parameter is optional.
# If the second parameter is 2, the string is treated as the name
# of a file which is to be copied to the end of the long file. In this case,
# a string saying "The contents of file $$$: " is inserted before
# the contents of the file.
# If the second parameter is omitted or if it is 1, the specified string is 
# treated as a value. If it is not a reference, it is treated as a string and
# appended as is to the end of the file. If it is a reference, the Data::Dumper
# routine is used to recuresively traverse the data structure and the resulting
# string is appended to the end of the file.
# Each long entry is separated by "separator" lines. 
# Each separator contains the time (to the second) of the log entry.
sub GLog {
  my $string = shift;
  my $mode = shift;

  if (!defined($ENV{GLOG_FILE})) {
    return;
  }
  if (!defined($mode)) { $mode = 1; }

  # If more than one process, make sure only one is writing to log file at a time
  if ($NPROC_DEF > 1) {
    #  Get a lock on the lock file
    my $lockfile = "$FS_DEF->{dhome}/.loglock";
    open (LOCK, ">$lockfile") or GENERAL::error("Process $PRANK_DEF unable to open lock file to start logging!", 0);
    flock LOCK, 2;
  }


  my $separator = "*"x60;
  my $mes = "GLog: $separator\nGLog (process $PRANK_DEF):" . GENERAL::GetTime("sec") . "\n";
  # Simply append a value
  if ($mode == 1) {
    if (ref($string)) {
      my $out = Dumper($string);
      $mes .= $out . "\n";
    } else {
      $mes .= $string . "\n";
    }
    $mes = GENERAL::escape($mes, "\"");
    if (system("echo \"$mes\" >> $ENV{GLOG_FILE}") != 0) {
      GENERAL::error("Error logging message \"$mes\" to file $ENV{GLOG_FILE}!", 0);
    }
  # String is a file
  } elsif ($mode == 2) {
    $mes .= "GLog: contents of file $string:\n";
    $mes = GENERAL::escape($mes, "\"");
    $string = GENERAL::escape($string, "\"");
    if ((system("echo \"$mes\" >> $ENV{GLOG_FILE}") != 0) || (system("cat \"$string\" >> $ENV{GLOG_FILE}") != 0)) {
      GENERAL::error("Error logging the contents of file \"$string\" with header message \"$mes\" to file $ENV{GLOG_FILE}!", 0);
    }
  } else {
    GENERAL::error("Error while logging: unknown mode of logging - $mode!", 0);
  }

  # If more than one process, make sure only one is writing to log file at a time
  if ($NPROC_DEF > 1) {
    #  Release lock
    flock LOCK, 8;
    close LOCK;
  }
} 

=head2

 Title   : getFileLock
 Usage   : GENERAL::getFileLock($filanem)
 Function: Gets a lock on a file with the given file name. If the file does
           not exist, it is created. Only returns when a lock is obtained.
 Returns : File handle.
 Args    : 1. The name of the lock file.

=cut

sub getFileLock {
  my $lockfile = shift;

  #  Get a lock on the lock file
  my $fh;
  open ($fh, ">$lockfile") or GENERAL::error("Process $PRANK_DEF unable to open lock file $lockfile", 0);
  flock $fh, 2;
  return $fh;
}

=head2

 Title   : releaseFileLock
 Usage   : GENERAL::releaseFileLock($filanem)
 Function: Releases a lock on a file and closes the file.
 Returns : Nothing.
 Args    : 1. File handle of the locked file.

=cut

sub releaseFileLock {
  my $fh = shift;

  flock $fh, 8;
  close $fh or GENERAL::error("Process $PRANK_DEF unable to close lock file", 0);
}


=head2

 Title   : fcntl_lock
 Usage   : GENERAL::fcntl_lock($fh, $lock_type)
 Function: Implements fcntl locking. Taken from recipe 7.22 Perl Cookbook with small modifications.
 Returns : 0 on success. If times out, returns the PID of the blocking process. On other errors, returns -1.
 Args    : 1. File handle of the locked file.
           2. Type of lock ("sh" meaning F_RDLCK - shared or "ex" meaning F_WRLCK - exclusive)
           3. Timeout - how long to wait (nagative means indefinitely). Note, waiting only happens if
              another process is holding a conflicting lock on the file.

=cut

sub fcntl_lock {
  my ($fh, $type, $timeout) = @_;
  GENERAL::requireArgs($fh, $type, $timeout);
  if (uc($type) =~ /EX/) { $type = F_WRLCK; }
  elsif (uc($type) =~ /SH/) { $type = F_RDLCK; }
  else { GENERAL::error("Unknown lock type '$type'"); }
  
  ##print "$$: Locking $start, $till\n";
  my $lock = struct_flock($type, SEEK_SET, 0, 0, 0);
  my $tim = time(); my $d = $tim; my $blocker = 0;
  while (!fcntl($fh, F_SETLK, $lock)) {
    if ($! == EAGAIN || $! == EACCES) {
#      $blocker = (struct_flock($lock))[-1]; # this does not work for some reason!
      if (($timeout >= 0) && ($tim + $timeout < time())) { $blocker = 1; last; }
      if ($d + 10 < time()) { printf("fcntl_lock: Waiting for an exclusive lock...\n"); $d = time(); }
    } else {
      return -1;
    }
  }
  return $blocker;
}  

=head2

 Title   : fcntl_unlock
 Usage   : GENERAL::fcntl_unlock($fh)
 Function: Implements fcntl unlocking. Taken from recipe 7.22 Perl Cookbook with small modifications.
 Returns : 0 on success. Non-zero otherwise.
 Args    : 1. File handle of the be un locked file.

=cut

sub fcntl_unlock {
  my $fh = shift;
  GENERAL::requireArgs($fh);
  my $lock = struct_flock(F_UNLCK, SEEK_SET, 0, 0, 0);
  if (!fcntl($fh, F_SETLK, $lock)) { return 0; }
  return 1;
}

=head2

 Title   : struct_flock
 Usage   : GENERAL::struct_flock($type, $whence, $start, $len, $pid)
 Function: Packs up the struct flock structure for use with fcntl (see man page for fcntl(2)). If used in list context,
           unpacks the packed struct.
 Returns : The packed struct (or a list of unpacked values)
 Args    : $type, $whence, $start, $len, $pid (se eman page fcntl(2)) OR the packed structure.

=cut

sub struct_flock {
  # c2ph says: typedef='s2 l2 i', sizeof=16
  my $FLOCK_STRUCT = 's s l l i';
  
  if (wantarray) {
    my ($type, $whence, $start, $len, $pid) = unpack($FLOCK_STRUCT, $_[0]);
    return ($type, $whence, $start, $len, $pid);
  } else {
    my ($type, $whence, $start, $len, $pid) = @_;
    return pack($FLOCK_STRUCT, $type, $whence, $start, $len, $pid);
  }
}


=head2

 Title   : createLocalSpace
 Usage   : GENERAL::createLocalSpace()
 Function: Creates a local directory (on the local drive) for fast temporary file I/O.
           Every process (if there are more than one) makes its own directory because
           it is possible that the processes are on different machines, so different
           local disks will need to be accessed. Global variable $LOCAL_DEF is modified.
 Returns : Nothing.
 Args    : None

=cut

sub createLocalSpace {
  my $rdir = shift;
  if (!defined($rdir)) { $rdir = "/tmp"; }
  my $ver = shift;
  $ver = 1 if (!defined($ver));

  my $cdir = GENERAL::GetDir();
  $LOCAL_DEF_ROOT = $rdir;
  GENERAL::cchdir($LOCAL_DEF_ROOT);
  # number of minutes of idle time before local directory can be reused
  $LOCAL_TIMEOUT = 120;

  my $username = GENERAL::GetUser();
  # Get a lock on a file so that only one process operates on the same
  # hard drive at a time
  my $lfh = GENERAL::getFileLock(".energy\_$username.lock");

  # Create username directory, if needed
  if (! -e "$username/") {
    GENERAL::cmkdir($username);
  }

  # If our temp directory does not exist on this hard drive, create it.
  # Otherwise, we'll overwrite it - this way no need to clean up after a failure.
  my $tmp_root = "$username/energy_temp";
  if (! -e $tmp_root) {
    GENERAL::cmkdir($tmp_root);
    my $local_dir = "$tmp_root/run_1/";
    GENERAL::cmkdir($local_dir);
    $LOCAL_DEF = "$LOCAL_DEF_ROOT/$local_dir";
  } else {
#    my $mindex = 0;
    my $found = 0; # found an old one to reuse?
    my %used; # a hash of used directory names
    # If the directory does exist, see if there are any left-over directories
    # from previous runs, which are too old and re-use them.
    my $str = `ls $tmp_root`; chomp($str);
    my @list = split(" ", $str);
    foreach my $dir (@list) {
      if (-d "$tmp_root/$dir") {
        # was any file in this directory modified recently?
        my $new_files = `find $tmp_root/$dir/ -mmin -$LOCAL_TIMEOUT -type f`;
        $new_files .= `find $tmp_root/$dir/ -amin -$LOCAL_TIMEOUT -type f`;
        $new_files = GENERAL::Trim($new_files);
        if ($new_files eq "") {
          # if no new files, see if the process that owned this directory is still running
          my $takenf = "$LOCAL_DEF_ROOT/$tmp_root/$dir/.taken";
          if (-e $takenf) {
            my $owner = `cat $takenf`; $owner = GENERAL::Trim($owner);
            if (GENERAL::isInteger($owner)) {
              next if (kill 0, $owner);
            }
          }
          # if no, use it
          $LOCAL_DEF = "$LOCAL_DEF_ROOT/$tmp_root/$dir";
          GENERAL::csystem("rm -fr $LOCAL_DEF/*");
          $found = 1;
          last;
        }
        $used{$dir} = 1;
      }
    }
    if (!$found) {
      my $i = 0;
      while (1) {
        my $dir = "run\_$i";
        $i++;
        next if (defined($used{$dir}));
        my $local_dir = "$tmp_root/$dir";
        GENERAL::cmkdir($local_dir);
        $LOCAL_DEF = "$LOCAL_DEF_ROOT/$local_dir";
        last;
      }
    }
  }

  GENERAL::csystem("echo '$$' > $LOCAL_DEF/.taken"); # before lock is released, create a new file so that other processes know this directory is taken
  GENERAL::releaseFileLock($lfh);
  print "Local temp directory for process $PRANK_DEF is $LOCAL_DEF on host " . GENERAL::GetMachine() . "\n" if ($ver);
  GENERAL::cchdir($cdir);

  return $LOCAL_DEF;
}

=head2

 Title   : destroyLocalSpace
 Usage   : GENERAL::destroyLocalSpace()
 Function: Frees local file space. Makes sure parallel processes go in order.
 Returns : Nothing.
 Args    : None

=cut

sub destroyLocalSpace {
  my $username = GENERAL::GetUser();
  my $lfh = GENERAL::getFileLock("$LOCAL_DEF_ROOT/.energy\_$username.lock");
  GENERAL::csystem("rm -rf $LOCAL_DEF");
  GENERAL::releaseFileLock($lfh);
  undef $LOCAL_DEF;
}

=head2

 Title   : touchLocalSpace
 Usage   : GENERAL::touchLocalSpace()
 Function: This is useful for long functions that do not do any local I/O - from time
           to time they should execute this function to make ensure owenership of the
           local space.
 Returns : Nothing.
 Args    : None

=cut

sub touchLocalSpace {
  if (defined($LOCAL_DEF)) {
    GENERAL::csystem("touch $LOCAL_DEF");
  }
}


=head2

 Title   : Trim
 Usage   : GENERAL::Trim()
 Function: Delete leading and trailing white space.
 Returns : 
 Args    : None

=cut

sub Trim {
  my @inp = @_;
  my $length = $#inp;

  # chop of leading and trailing blankspacing.
  foreach (@inp) {
    s/^\s+//;
    s/\s+$//;
  }

  if ($length == 0) {
    return $inp[0];
  } else {
    return @inp;
  }
}


=head2

 Title   : escape
 Usage   : GENERAL::escape("Some string possibly dangerous for the shell!")
 Function: Escapes all the characters possibly dangerous for the shell.
 Returns : The string with escaped characters.
 Args    : 1. the scring in question.
           2. optional: a list of characters to be escape. If not 
              specified, all characters "potentially dangerous" for the 
              shell will be escaped.

=cut

sub escape {
  my $str = shift;
  my $chars = shift;
  if (!defined($chars)) { $chars = ";<>\*\|`&\$!#\(\)\[\]\{\}:'\""; }

  $str =~ s/([$chars])/\\$1/g;
  return $str;
}


=head2

 Title   : REAPER
 Usage   : $SIG{CHLD} = \&GENERAL::REAPER;
 Function: A general dead child process reaper - avoids zombies. Taken from recipe 16.9 of Perl Cookbook.
 Returns : Nothing
 Args    : 
 
=cut

sub REAPER {
  my $stiff;

#  # First, get info about all child processes
#  my $chps = `ps --ppid $$ | grep -v "PID"`;
#  my @chps = split("\n", GENERAL::Trim($chps));
#  my %chps;
#  foreach my $chp (@chps) {
#    my @arr = split(" ", GENERAL::Trim($chp));
#    $chps{$arr[0]}{"name"} = substr($chp, 23);
#  }

  while (($stiff = waitpid(-1, &WNOHANG)) > 0) {
    GENERAL::warning("Child process with pid $stiff seems to haved died (status $?).");
  }
#  $SIG{CHLD} = \&GENERAL::REAPER;
}

# Merges two hash tables (A and B).
# Input: pointer to hash A.
#        pointer to hash B.
#        pointer to resulting hash C.
# C is an optional parameter. If not defined, a new hash is created
# and a pointer to it returned. If C is the same as either A or B,
# the function copies only the other hash into C (to save time).
sub mergeHashes {
  my $A = shift;
  my $B = shift;
  my $C = shift;
  my %C;
  
  if (!defined($C)) {
    $C = \%C;
  }
  
  my ($k, $v);
  if ($C != $A) {
    while ( ($k,$v) = each(%$A)) {
      $C->{$k} = $v;
    }
  }
  
  if ($C != $B) {
    while ( ($k,$v) = each(%$B)) {
      $C->{$k} = $v;
    }
  }
  
  return $C;
}

# function: AskInputString
# description: give a question, get back the answer
# usage: my $answer = &AskInputString($question);

sub AskForInputString {

    my $question=shift;
    print  $question, "\n";
    my $hint = shift;
    print $hint, "\n" if (  defined $hint);
    my $answer = <STDIN>;
    chomp $answer;
    return $answer;

}


# function: AskForManyLines
# description: give a question, get back the answer
# usage: my $answer = &AskForManyLines($question);

sub AskForManyLines {

    my $question=shift;
    print  $question;
    my $hint = shift;
    print $hint  if (  defined $hint);
    my $finalanswer="";
 
    READINPUT: 
      while (<>) {
            last READINPUT if ($_ =~ /END|end|End/);
            $finalanswer .=$_; 
      }

    return $finalanswer;

}



## function: $fullpath = FindExe(name)
## attempts to locate the full path of a executable
## with the given name

sub FindExe {

    my $name=shift;
  
    foreach my $p ( split(/:/,$ENV{PATH}) ) {

    return "$p/$name"
        if (-d $p && -x "$p/$name");
    }
    return undef;
}

# function: Log(tag,text)
# writes a debug message to a log file
# if one has been set previously with setLogFile
# The message consists of a tag and the actual text.

sub Log {
#    my $tag=shift;
#    my $s=shift;
#    printf $logfile "## %s ##  %s\n",$tag,$s
#    if (defined $logfile);
  GENERAL::error("This function has been deprecated!");
}

## function: $list = ReadFormat(handle,reclen,numrec)
## reads a formatted line with fields of a specific record length

sub ReadFormat {
    my $handle=shift;
    my $reclen=shift;
    my $numrec=shift;
    
    my @arr=();
    
  FORMATINP:
    while (<$handle>) {
    for (my $i=0; $i+$reclen<=length($_); $i+=$reclen) {
        push(@arr,substr($_,$i,$reclen));
        last FORMATINP if ($#arr+1>=$numrec);
    }      
    }
    printf STDERR "read only %d elements (%d requested)\n",$#arr+1,$numrec
    if ($#arr+1<$numrec);
    
    return @arr;
}

## function: $buffer = ReadData(handle,pattern)
## reads data from the given handle into a buffer until 
## a pattern is found.

sub ReadData {
    
    my $handle=shift;
    my $pattern=shift;
    
    my $buffer="";
    while (<$handle>) {
    return $buffer 
        if ($pattern ne "" && $_=~/$pattern/);
    $buffer.=$_;
    }
    return $buffer;
}


## function: $status = ValidGzip(filename)
## checks whether a file is a valid compressed gzip file.

sub ValidGzip {
    my $f=shift;
    
    my $ret=`gzip -t $f`;
    
    return ($ret!~/not/);
}

## function: $status = CheckFile(filename)
## checks whether a file is available and readable

sub CheckFile {
    my $f=shift;
    return ((-r $f && !-z $f) || (-r "$f.gz" && validGzip("$f.gz")));
}

## function: $value = dihedral(coor1,coor2,coor3,coor4)
## calculates the dihedral value for four points given
## by the coor* arguments. 

sub Dihedral {
    my $c1 = shift;
    my $c2 = shift;
    my $c3 = shift;
    my $c4 = shift;
    
    my $dx12 = $c1->{xcoor}-$c2->{xcoor};
    my $dy12 = $c1->{ycoor}-$c2->{ycoor};
    my $dz12 = $c1->{zcoor}-$c2->{zcoor};
    
    my $dx23 = $c2->{xcoor}-$c3->{xcoor};
    my $dy23 = $c2->{ycoor}-$c3->{ycoor};
    my $dz23 = $c2->{zcoor}-$c3->{zcoor};
    
    my $dx43 = $c4->{xcoor}-$c3->{xcoor};
    my $dy43 = $c4->{ycoor}-$c3->{ycoor};
    my $dz43 = $c4->{zcoor}-$c3->{zcoor};
    
    my $px1 = $dy12*$dz23-$dy23*$dz12;
    my $py1 = $dz12*$dx23-$dz23*$dx12;
    my $pz1 = $dx12*$dy23-$dx23*$dy12;
    
    my $np1 = sqrt($px1*$px1+$py1*$py1+$pz1*$pz1);
    
    $px1/=$np1;
    $py1/=$np1;
    $pz1/=$np1;
    
    my $px2 = $dy43*$dz23-$dy23*$dz43;
    my $py2 = $dz43*$dx23-$dz23*$dx43;
    my $pz2 = $dx43*$dy23-$dx23*$dy43;
    my $np2 = sqrt($px2*$px2+$py2*$py2+$pz2*$pz2);
    
    $px2/=$np2;
    $py2/=$np2;
    $pz2/=$np2;
    
    my $dp12 = $px1*$px2+$py1*$py2+$pz1*$pz2;

    # sometimes there is some sort of a strange numerical instability and
    # $dp12*$dp12 ends up being 1 + 10^-16.
    my $sin2 = 1 - $dp12*$dp12;
    if ($sin2 < 0) {
      if ($sin2 < -10**-4) {
        GENERAL::error("The square of a sine was too small ($sin2 - below 0)!");
      }
      $sin2 = 0;
    }
    
    my $angle = pi()/2.0-atan2($dp12,sqrt($sin2)); # /
    
    my $px3 = $py1*$pz2-$py2*$pz1;
    my $py3 = $pz1*$px2-$pz2*$px1;
    my $pz3 = $px1*$py2-$px2*$py1;
    
    my $dp233 = $px3*$dx23+$py3*$dy23+$pz3*$dz23;
    
    if ($dp233 > 0.0) {
      $angle=-$angle;
    }
    
    return ($angle/pi())*180.0;
}

sub angle {
  my $c1 = shift;
  my $c2 = shift;
  my $c3 = shift;

  my @v21 = ($c1->{xcoor} - $c2->{xcoor}, $c1->{ycoor} - $c2->{ycoor}, $c1->{zcoor} - $c2->{zcoor});
  @v21 = GENERAL::getUnitVector(\@v21);
  my @v23 = ($c3->{xcoor} - $c2->{xcoor}, $c3->{ycoor} - $c2->{ycoor}, $c3->{zcoor} - $c2->{zcoor});
  @v23 = GENERAL::getUnitVector(\@v23);
  my $cos = GENERAL::dot(\@v21, \@v23);
  return (180/pi()) * atan2(sqrt(1 - $cos**2), $cos);
}

## function: ($avg,stddev,$n) = Average(list)
## calculates average and standard deviation of
## values given in the list argument

sub Average {
    
    my $list=shift;
    
    my $sum=0.0;
    my $sumsq=0.0;
    my $n=0;
    
    die "list not defined"
    if (!defined $list || $#{$list}<0);
    
    return ($list->[0],0.0,1)
    if ($#{$list}==0);
    
    foreach my $s ( @{$list} ) {
    my $val=$s;
    $sum+=$val;
    $sumsq+=$val*$val;
    $n++;
    }
    
    my $avg=$sum/$n;
    my $stddev=sqrt(($sumsq-$n*($avg*$avg))/($n-1));
    
    return ($avg,$stddev,$n);
}


## function: $status = DataAvailable(handle,timeout)
## checks whether input data is available at the given
## handle.

sub DataAvailable {
    my $handle=shift;
    my $timeout=shift;

    $timeout=0.5 if ($timeout <= 0);
    
    my ($nfound,$rmask,$nread,$line);
    
    $rmask="";
    vec($rmask,fileno($handle),1)=1;
    ($nfound,$rmask)=select($rmask,undef,undef,$timeout);
    
    return $nfound;
}





#  function: $list = FragListFromOption(string)
#  converts a fragment list command line argument
#  into a list data structure

sub FragListFromOption {

    my $fragoption=shift;
    my $fraglist=();
    
    die "invalid list specification: >$fragoption<"
    if ($fragoption !~ /^[A-Za-z0-9:\.=_\-]+$/);

    my $fdef=undef;
    foreach my $f (split(/=/,$fragoption)) {
    my @fl=split(/_/,$f);
    my @l=split(/:/,$fl[0]);
    my $rec={};
    my ($tchain,$tfrom)=($l[0]=~/([A-Za-z]*)([0-9]*)/);
    $rec->{chain}=$tchain if ($tchain ne "");
    $rec->{from}=$tfrom if ($tfrom ne "");
    if (defined $l[1]) {
        ($rec->{to}=$l[1])=~s/^[A-Za-z]+//;
    } else {
        $rec->{to}=$rec->{from};
    }
    if (defined $fl[1]) {
        $rec->{force}=$fl[1];
        $fdef=$fl[1];
    } else {
        $rec->{force}=$fdef;
    }
    push(@{$fraglist},$rec);
    }
    return $fraglist;
}


## function: $string = FragOptionFromList(list[,force])
## generates a command line option string from
## a fragment list data structure. A force constant may be
## given separately to be included in the output 
## if the force field is not defined in the input list.

sub FragOptionFromList {
    my $fraglist=shift;
    my $force=shift;
    
    my $lastforce;
    
  if (defined $fraglist) {
      my @list=();
      foreach my $f (@{$fraglist}) {
      my $str=($f->{from}==$f->{to})?"$f->{from}":"$f->{from}:$f->{to}";
      
      $str=(uc $f->{chain}).$str if (defined $f->{chain} && $f->{chain} ne "");
      
      my $fval=(defined $force)?
          $force:((defined $f->{force})?$f->{force}:undef);
      
      if (defined $fval) {
          if (!defined $lastforce || $fval ne $lastforce) {
          $str.="_$fval";
          $lastforce=$fval;
          }
      }
      push(@list,$str);
      }
      return join("=",@list); 
  } else {
      return "";
  }
}



## function: $list = FragListFromArray(arr[,chain])
## generates a fragment list data structure
## from an array of residues

sub FragListFromArray {
    my $arr=shift;
    my $chain=shift;
    
    my $list=();
  
    my $trec={};
    $trec->{chain}=$chain;
    $trec->{from}=$arr->[0];
    for (my $i=1; $i<=$#{$arr}; $i++) {
    if ($arr->[$i]-$arr->[$i-1]>1) {
        $trec->{to}=$arr->[$i-1];
        push(@{$list},$trec);
        $trec={};
        $trec->{chain}=$chain;
        $trec->{from}=$arr->[$i];
    }
    }
    $trec->{to}=$arr->[$#{$arr}];
    push(@{$list},$trec);
    
  return $list;
}



## function: $list = LimRange(list,min,max)
## generates a sublist from the list given
## as argument with values that lie within
## the range given by min and
## max.

sub LimRange {
  my $list=shift;
  my $min=shift;
  my $max=shift;

  my $outlist=();

  foreach my $s ( @{$list} ) {
    my $val=$s->{val};
    push(@{$outlist},$s)
      if ((!defined $min || $val>=$min) && 
      (!defined $max || $val<=$max));
  }
  
  return $outlist;
}


## function: ($list,$lastlimit,$average,$stddev,$n) = 
## function:  limCore(list[,parameters])
## generates a sublist by excluding data points that
## are further away from the mean than a multiple of the standard
## deviation. Parameters in hash-style key=>value pair form can be given 
## to select whether all values below (lower)
## and above (upper) are cropped, to change
## the multiple mult (default: 3) and to change a 
## tolerance limit tol for terminating the iterative 
## solution.<BR>
## Return values are the reduced list, the last limit value used,
## and the average, standard deviation and number of points for the
## reduced list.

sub LimCore {
    my $list=shift;
    
    die "list not defined"
    if (!defined $list || $#{$list}<0);
    
    return ($list,0.0,$list->[0]->{val},0.0,1)
    if ($#{$list}==0);
    
    my %par=@_;
    
    my $lower=(defined $par{lower})?$par{lower}:1;
    my $upper=(defined $par{upper})?$par{upper}:1;
    my $mult=(defined $par{mult})?$par{mult}:3.0;
    my $tol=(defined $par{tol})?$par{tol}:0.05;
    
    my $avg;
    my $stddev;
    my $limit=1.0E20;
    my $lastlimit;
    
    my $n; 
    
    do {
    $lastlimit=$limit;
    
    $n=0;
    my $sum=0.0;
    my $sumsq=0.0;
    
    foreach my $vl ( @{$list} ) {
        my $v=$vl->{val};
        if ((!$upper || $v<($avg+$limit)) && 
        (!$lower || $v>($avg-$limit))) {
        $sum+=$v;
        $sumsq+=$v*$v;
        $n++;
        }
    }
    $avg=$sum/$n;
    $stddev=sqrt(($sumsq-$n*($avg*$avg))/($n-1));
    $limit=$mult*$stddev;
    } while($lastlimit-$limit>$tol);
    
    my $outlist=&GenUtil::limRange($list,
     ($lower)?($avg-$lastlimit):undef,($upper)?($avg+$lastlimit):undef);
    
    return ($outlist,$lastlimit,$avg,$stddev,$n);
}



# Function: Ensemble
# Description: Take a list of energy and optional Temperature. 
#              give the emsemble averaged energy data structure

sub Ensemble {
    my $elist = shift;
    my $temperature = shift;
    $temperature=300.0 if (!defined $temperature);

    # define boltzmann constant in kcal/(mol K)
    my $boltzmann=1.987191E-3;
    my $Constant= (-1.0)/($boltzmann * $temperature);


    my %ensemble;
    my @templist = ();
    my $sum = 0.0;
    # To avoid all exponents being zero due to large energies, find
    # the minimal energy and subtract it from exponents in both the numerator
    # and the denamenator (equivalent to multiplying by exp(Emin/kt)/exp(Emin/kt)
    my $men = GENERAL::min($elist);

    for (my $i=0; $i <= $#{$elist}; $i++) {
      $templist[$i]=exp($Constant*($elist->[$i]-$men));
      push (@{$ensemble{'energy'}},  $elist->[$i]);
      $sum += $templist[$i];
    }

    # Calculate ensemble average

    my @prob;
    my $average;
    # we want to find out waht is most probable one
    my $best=0.0;
    my $indexofbest;

    for (my $i=0; $i <= $#{$elist}; $i++) {
      $prob[$i]=$templist[$i]/$sum;
      if ( $prob[$i] >= $best) {
           $best = $prob[$i];
           $indexofbest = $i;
      }
      push (@{$ensemble{'probability'}}, $prob[$i]);
      $average +=  $prob[$i]*$elist->[$i];
    }

    $ensemble{'average'} = $average;
    $ensemble{'best'} = $best;
    $ensemble{'indexofbest'} = $indexofbest;
    $ensemble{'partition'} = $sum;
 
    return \%ensemble;
}
  

# Function: EnsembleProp
# Description: Take a list of property we want to average and ensemble reference.
#              give the emsemble averaged property

sub EnsembleProp {

    my $proplist=shift;
    my $ensemble=shift;
    my $average;

    for (my $i=0; $i <= $#{$proplist}; $i++) {
    $average += $proplist->[$i] * $ensemble->{probability}->[$i];
    } 

    return $average;
}


# Function: RoundIt
# Description: round the number into integer
#              used in the backbone-dependent library and 
#              phi, psi value

sub RoundIt {

    my $i = shift;

    my $j = $i/10.0;
    my $intf = sprintf("%3d", $j);
    my $newj = $intf*10;
    
    # rounding number is done in follwoing fashion: 
    #    -15.2  --> -20
    #    -14.9  --> -10
    #     15.2  -->  20
    #     14.9  -->  10
    if ( ($i-$newj) > 5.0 && $i> 0.0 ) {
      return $newj+10;
    } elsif  ( ($i-$newj) <5.0 && $i> 0.0 ) {
      return $newj;
    } elsif  ( abs($i-$newj) > 5.0 && $i<0.0) {
      return $newj-10;
    } elsif  ( abs($i-$newj) <5.0  && $i<0.0) {
      return $newj;
    }

}
            
sub FixCRDTotalNum {

    my $fname=shift;
    my $inp = &GENERAL::GetInFH($fname);
    my $newfname=shift;
    my $newfile = &GENERAL::GetOutFH($newfname);

    my $linecount=0;
    my $total = `cat $fname | wc -l`;
    chomp $total;
    $total -=1;

    while ( <$inp> ) {
       $linecount++;
       print $newfile $_ if ($linecount >=2);
       if ($linecount == 1) {
            printf $newfile "%5d\n", $total;
       }
    }
 
    system (" mv  $newfname $fname");

}



=head2 

 Title   :  MergeAddArray
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut



sub MergeAddArray {
    
    my $arr1=shift;
    my $arr2=shift;
    my $lamda1=shift;
    my $lamda2=shift;

    $lamda1=1.0 if (!defined $lamda1);
    $lamda2=1.0 if (!defined $lamda2);

    if ($#{$arr1} != $#{$arr2} ) {
      die "Error in MergeAddArray: the two arrays don't have same size (" . scalar(@{$arr1}) . " and " . scalar(@{$arr2}) . ")!\n";
    }

    my $newarr=[];
    for (my $i=0; $i<=$#{$arr1}; $i++) {
      $newarr->[$i]=$lamda1*$arr1->[$i]+ $lamda2*$arr2->[$i];
    }

    return $newarr;
}



=head2 

 Title   : MergeSubtractArray
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut




sub MergeSubtractArray {
    
    my $arr1=shift;
    my $arr2=shift;

    if ($#{$arr1} != $#{$arr2} ) {
    die " two arrays don't have same size!\n";
    }

    my $newarr=[];

    for (my $i=0; $i<=$#{$arr1}; $i++) {
    $newarr->[$i]= $arr1->[$i]- $arr2->[$i];
    }

    return $newarr;
}



=head2 

 Title   : MergeAppendArray
 Usage   : 
 Function: This function is designed for merging tables with the same
           number of rows. A table in this context is an array of strings,
           with each string referring to a row and containing a concatination
           of cells in that row (separated by a space " "). This function
           accepts an array of references to tables and returns a reference to
           a bigger table. The rows of this bigger table are formed by appending
           rows of the smaller tables. Naturally, all tables thus have to have
           the same number of rows. Rows of the output table are actually not
           strings, but arrays of values, which are obtained by braking up the
           rows of the input tables. So the output table is really a reference
           to an array of references to arrays.
 Returns : a reference to the merged table.
 Args    : 

=cut

sub MergeAppendArray {

    my @inp=@_;
    my $size=$#{$inp[0]};

    # check the size of array
    foreach my $i (@inp) {
      die "Error in MergeAppendArray: sizes don't agree with each other.\n" if ( $size != $#{$i});
    }

    my $newref=[];
    for (my $i=0; $i<=$#{$inp[0]}; $i++) {
      # split each line and append it.
      my @lines;
      foreach my $j (@inp){
          my @line=split(" ", $j->[$i]);
          push (@lines, @line);
      }
      push (@{$newref}, \@lines);
    }
    return $newref;
}



=head2 

 Title   :  gaussianRand
 Usage   :  gaussianRand()
 Function:  generates a gaussian distributed random number with mean 0 and standard deviation of 1
 Example :  To generate a gaussian randomly distributed number with 
            mean value of A, standard deviation of B, do
            A + gaussianRand()*B;
 Returns :  
 Args    : 

=cut

sub gaussianRand {
    my ($u1, $u2);  # uniformly distributed random numbers
    my $w;          # variance, then a weight
    my ($g1, $g2);  # gaussian-distributed numbers

    do {
        $u1 = 2 * rand() - 1;
        $u2 = 2 * rand() - 1;
        $w = $u1*$u1 + $u2*$u2;
    } while ( $w >= 1 );

    $w = sqrt( (-2 * log($w))  / $w );
    $g2 = $u1 * $w;
    $g1 = $u2 * $w;
    # return both if wanted, else just one
    return wantarray ? ($g1, $g2) : $g1;
}


=head2

 Title   :  runCharmm2
 Usage   :  runCharmm2
 Function:
  Example:
 Returns :
 Args    :

=cut

sub runCharmm2 {

   my $charmm =shift;
   my $type=shift;
   my $pairnum=shift;
   my $dir=shift;
   $dir='.' if (!defined $dir);
   my $currentdir=`pwd`;
   chomp($currentdir);

   chdir($dir);

   for(my $j=1; $j<=$pairnum; $j++)  {
        my $inp=$type.$j.'.inp';
        my $output=$type.$j.'.out';
        `$charmm < $inp 2>&1 1>& $output`;
        my $charmmline = "$charmm < $inp 2>&1 1>& $output";
        print $charmmline, "\n";
        &GENERAL::Remove("$output");
        &GENERAL::Remove("$inp");

  }

  chdir($currentdir);

}




=head2 

 Title   : pairResMin 
 Usage   : 
 Function: This is a perl version of 
           calculate the minimal distance between 
           a pair of rotamer. It was implemented for 
           the purpose of convenience, not for performace.  
           When you do real calculation. you should always using 
           c version code.  

 Returns : 
 Args    : 

=cut



sub pairResMin {


    my $site1=shift;
    my $res1=shift;
    my $site2=shift;
    my $res2=shift;
    my $pairnum=shift;
    my $fname1='rotmulte'.$site1."-".$res1.'.dat';
    my $fname2='rotmulte'.$site2."-".$res2.'.dat';
    my $out='pairmin'.$pairnum.'.dat';
    my $outf=&GENERAL::GetOutFH($out);

    my $f1=&GENERAL::GetInFH($fname1);
    my $f2=&GENERAL::GetInFH($fname2);
    
    my @a1=<$f1>;
    my @a2=<$f2>;

    chomp @a1;
    chomp @a2;
    
    my $head1=shift @a1;
    my $head2=shift @a2;
    my ($atomno1, $rotnum1) = split(" ", $head1);
    my ($atomno2, $rotnum2) = split(" ", $head2);

    my $s1=$#a1;
    my $s2=$#a2;

    my $counter=0; 
    # loop over all pair of rotamers

    for (my $i=1; $i <=$rotnum1; $i++) {
    for (my $j=1; $j <=$rotnum2; $j++) {
        
        $counter++;
        my $mind=9999.0;
        # loop over all the atom pairs
        for (my $k=($i-1)*$atomno1; $k< ($i*$atomno1); $k++) {
        for (my $l=($j-1)*$atomno2; $l< ($j*$atomno2); $l++) {
            # calculate mini distance 
            my @xyz1=split(" ", $a1[$k]);
            my @xyz2=split(" ", $a2[$l]);
            my $dx=$xyz1[0]-$xyz2[0];
            my $dy=$xyz1[1]-$xyz2[1];
            my $dz=$xyz1[2]-$xyz2[2];
            my $dd=$dx*$dx+ $dy*$dy +$dz*$dz;
            $mind=$dd if ($dd<$mind);
        }
        }
        
        $outf->printf("%12.3f\n", sqrt($mind));
    }
    }
} 


=head2 

 Title   : 
 Usage   : 
 Function:  A utilities to clean up the directory
 Returns : 
 Args    : 

=cut



sub ClearPairGZ {

    my $target=shift;
    my $dir=shift; 
    $dir='./' if (!defined $dir);

    my $pairmulte;
    my $paireef;

    if (!defined $target) {
        $pairmulte=1; 
        $paireef=1;
    } elsif ( uc($target) =~ /MULTE/ ) {
        $pairmulte=1;
        $paireef=0;
    } elsif ( uc($target) =~ /EEF/ ) {
        $pairmulte=0;
        $paireef=1;
    }

    if ($pairmulte == 1) {
        my $multe=$dir."pair/Multe";
        chdir($multe);
        system("/bin/rm -rf *.gz");
        chdir("../../");
    }

    if ($paireef ==1 ) {
        my $eef=$dir."pair/EEF";
        chdir($eef);
        system("/bin/rm -rf *.gz");
        chdir("../../");
    }
} 


=head2 

 Title   : readDesignOut  
 Usage   : $config=&readDesignOut($parm);
 Function: go to the design directory grab the output
           get the final configuration 
 Returns : reference to array of hash $config
 Args    : reference to parameter hash

=cut


sub readDesignOut {
  
    # read in all the parameters related with design 
    my $parm=shift;
    my $method=$parm->{'packing_method'};
    my $parentdir=$parm->{'design_parent_dir'};   
    my $subdir=$parm->{'design_sub_dir'};
    my $ebin=$parm->{'energy_bin'};
    my $inp =$parm->{'design_input'};

    my $design=DESIGN::new($method, $inp, $ebin);
    my $direct=$parentdir.'/'.$subdir;
    $design->parseOutput($direct);

    return $design->{'solution'}->{'results'};

}


=head2 

 Title   : ReadQuickRotinp 
 Usage   : $sites=&ReadQuickRotinp($rotinp);
 Function: gives the information of chain id and resid information
           by quickly look up the rot.inp file. instead of going 
           through rotamer object methods
 Returns : $sites: reference to array of array 
 Args    : reference to parameter hash

=cut

sub ReadQuickRotinp {
  
    my $rotinp=shift;
    my $site=&GENERAL::GetInFH($rotinp);
    my @sites=<$site>;
    chomp @sites;
    my $sites=[];
    for (@sites) {
        chomp $_;
        my @lines=split(" ", $_);
        my @newlines = splice(@lines,0,2);
        push (@{$sites}, \@newlines);
    }
    close $site;
    return $sites;
}



=head2

 Title   : ReadQuickRotinp2
 Usage   : $sites=&ReadQuickRotinp2($rotinp);
 Function: gives the information of chain id and resid information
           by quickly look up the rot.inp file. instead of going
           through rotamer object methods
           similar to ReadQuickRotinp method. But has the additional information 
           about residue name
 Returns : $sites: reference to array of array
 Args    : reference to parameter hash

=cut


sub ReadQuickRotinp2 {
 
    my $rotinp=shift;
    my $site=&GENERAL::GetInFH($rotinp);
    my @sites=<$site>;
    chomp @sites;
    my $sites=[];
    for (@sites) {
        chomp $_;
        my @lines=split(" ", $_);
        splice(@lines,2,1);
        push (@{$sites}, \@lines);
    }
    return $sites;
}




=head2

 Title   : PrintParm
 Usage   : &printParm($parm);
 Function: print Parameter on the screen
 Returns : none
 Args    : 

=cut

sub PrintParm {
   my $parm = shift;
   foreach my $key (keys %{$parm}) {
       printf("%s -> '%s'\n", $key, $parm->{$key});
   }

}


=head2 

 Title   : ReadParm  
 Usage   : my $parm = GENERAL::ReadParm($filename);
 Function: Reads a parameter file formatted as:
           parameter_name = parameter_value
           ...
           Optionally, a separator character can be specified, in which 
           case it does not have to be "=".
 Returns : Reference to param hash
 Args    : 1. File name
           2. Optional: separator string that separates names from values.
           3. Optional: array reference. If specified, parameter names
              will be pushed onto this array in the order read.

=cut


sub ReadParm {
  my $fname = shift;
  my $sep = shift;
  if (!defined($sep)) { $sep = "="; }
  my $names = shift;

  my $inp = GENERAL::GetInFH($fname);
  my %parm;

  while (<$inp>) {
    next if ( /^\#/ );
    next if ( /^\s+$/);
    next if (!/$sep/);
    my @line = split(/$sep+/, $_, 2);
    $parm{GENERAL::Trim($line[0])} = GENERAL::Trim($line[1]);
    if (defined($names)) { push(@$names, GENERAL::Trim($line[0])); }
  }

  return \%parm;
}


=head2 

 Title   : SetParm
 Usage   : &SetParm($parmref, $key, $value)
 Function: reset the parameter value
 Returns : none
 Args    : (1) reference to parm hash
           (2) the key enetry
           (3) the value entry 
=cut

sub SetParm {
  
    my $parm=shift; 
    my $key=shift;
    my $value=shift;
    $parm->{$key}=$value;

}


=head2 

 Title   : BuildXcoorHash 
 Usage   : 
 Function: Take a pdb file, build a look up table 
           key is x,y,z-coord, value is  
 Returns : 
 Args    : 

=cut


sub BuildXcoorHash {

     my $fname=shift;
     my $fh=&GENERAL::GetInFH($fname);

     my %xcoorhash;
     while (<$fh>) {
        if ( /^ATOM/) {
            my @line = split(" ", $_);
            my $xkey=sprintf("%10.4f", $line[6]);
            my $ykey=sprintf("%10.4f", $line[7]);
            my $zkey=sprintf("%10.4f", $line[8]);
            my $key=$xkey.$ykey.$zkey;
            my $value=$line[4].$line[5].$line[3];
            $xcoorhash{$key}=$value;
        }
     }

     return \%xcoorhash;
}




sub ReadCharge {

    my $fname=shift;
    my $fh=&GENERAL::GetInFH($fname);

    my %parmhash;
    my $key;

    while (<$fh>) {
        next if (/^!/);
        next if (/^$/);
        my @lines=split(" ", $_);


        if ( /^RESI /) {
            $key=$lines[1];
            $parmhash{$key}={};
        }

        if ( /^ATOM /) {
            $parmhash{$key}->{$lines[1]}=$lines[2];
        }

    }

    return \%parmhash;
}


=head2

 Title   :
 Usage   :
 Function:
 Returns :
 Args    :

=cut

sub WriteChargeFile {

    my $outf=shift;
    my $inf=shift;
    my $chargeparm=shift;
    my $unfoldbit=shift;
    $unfoldbit=0 if (!defined $unfoldbit);

    my $infh=&GENERAL::GetInFH($inf);
    my $outfh=&GENERAL::GetOutFH($outf);

    $outfh->print('atom__resnumbc_charge_', "\n");

    my $lcount=0;
    while (<$infh>) {
        if ( /^ATOM/ ) {
            $lcount++;
            my ($atomname, $resname, $chain, $resnum)
                = &GENERAL::Trim(substr($_, 12,4), substr($_,17,4),substr($_,21, 1),substr($_, 22,4));
            print $lcount, " Atom Name\n" if ( !defined $atomname);
            print $lcount, " Res Name\n" if (!defined $resname);
            print $lcount, " Chain\n" if (!defined $chain);
            print $lcount, " Res Number\n" if (!defined $resnum);
            print $lcount, " Charge \n" if (!defined  $chargeparm->{$resname}->{$atomname});
            if ($unfoldbit==1) {  
        $chain='A';
            }
            $outfh->printf("%-4s%5s%4d%1s%8.5f\n", $atomname, $resname, $resnum, $chain, $chargeparm->{$resname}->{$atomname});
        }

    }

}



=head2 

 Title   : ReadDesignConfig   
 Usage   : my $rec=&ReadDesignConfig($parm);
 Function: read in the solution file from design output           
 Returns : reference to array of array 
 Args    : (1) param hash

=cut

sub ReadDesignConfig {

    my $parm=shift;
    my $parentdir=$parm->{'design_parent_dir'};
    my $subdir=$parm->{'design_sub_dir'};
    my $solution = $parm->{'solution'};
  
    my $cur=`pwd`;
    chomp $cur;
    
    my $direct=$parentdir.'/'.$subdir;
    chdir($direct);
    my $fh=&GENERAL::GetInFH($solution);

    my $rec=[];    
    while (<$fh>) {
      next if (/^#/);
      my @lines=split(" ", $_);
      push (@$rec, \@lines);
    }

    chdir($cur);
    return $rec;
}


=head2

 Title   : WriteTemplatePDB
 Usage   :
 Function:
 Returns :
 Args    :

=cut


sub WriteTemplatePDB {
    
    my $parm=shift;
    my $workpdb=shift;
    my $isdesignres=shift;

    my $workdir=$parm->{'pdb_self_rot'};
    my $pdb=&GENERAL::GetInFH("../$workpdb");   
    my $out=&GENERAL::GetOutFH("$workdir/backbone.pdb");

    
    while (<$pdb>) {
    if (/^ATOM/) {
           next if ( / H /);
           my ($id, $name, $resnum) = &GENERAL::Trim(substr($_,21, 1), substr($_,17,4),substr($_, 22,4));
           my $key2= $id.$resnum.$name;
       if (!exists $isdesignres->{$key2} || / N | CA | C | O | H / ) {
        $out->print($_);
       } else {
            next;
           }
        } elsif ( /^TER/ ) {
           $out->print($_);
    } else {
       next;
        }
    }
}



=head2 

 Title   : WriteSelfPDB
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut


sub WriteSelfPDB {
    my $parm=shift;
    my $solutionkey=shift;
    my $xraykey=shift;
    my $i=shift;
    my $xraybit=$parm->{'xbit'};
    my $tbit=$parm->{'self_tbit'};
#    my $dbit= $parm->{'self_dbit'};
 
    $xraybit=0 if (!defined $xraybit);
    
    my $key;
    if ($xraybit == 0) {
    $key=$solutionkey->[$i];
    } else {
        $key=$xraykey->[$i];
    }

    my $prefix;
    my $workdir;
    if ($tbit ==1) {
       $prefix='self-';
       $workdir=$parm->{'pdb_self_asp'};
    } else {
       $prefix ='rotself-';
       $workdir=$parm->{'pdb_self_rot'};
    }

     #  print $key, "\n";
    my $pdbfile= $parm->{'mixed_pdb_file'};
    my $newpdb=$prefix.$i.'.pdb';
    my $out=&GENERAL::GetOutFH("$workdir/$newpdb");  
    my $pdb=&GENERAL::GetInFH("../$pdbfile");
    printf "prepare file:....   %15s\n", $newpdb;

    my $charge=lc($parm->{'charge_marker'});

    while (<$pdb>) {
       if (/^ATOM/) { 
            my ($id, $name, $resnum, $rotnum)
                = &GENERAL::Trim(substr($_,21, 1), substr($_,17,4),substr($_, 22,4), substr($_,26,2));

             my $key2= $id.$resnum.$name.$rotnum;
             if ( $key eq $key2 )  {
               $out->print($_);
             } elsif ( $tbit ==1 && substr($_, 27,1) eq " " ) {
                   substr($_, 21, 1) = sprintf("%1s", 'X') if ($charge eq 'on');
                  # print $_;
        #           print $charge, "\n";
                   $out->print($_);     
             } else { 
                   next;
             }
       } elsif (/^TER/) {
             $out->print($_);
       } else {
             next;
       }
     
    }
}



=head2 

 Title   : WriteSelfDelphiPDB 
 Usage   :  
 Function: write self pdb file for delphi input
 Returns : 
 Args    : 

=cut


sub WriteSelfDelphiPDB {

    my $parm=shift;
    my $workpdb=shift;
    my $designkey=shift;
    # index to charge the side chain 
    my $i=shift;
    my $workdir=$parm->{'pdb_self_delphi'};
        my $prefix='delphi-';
    my $newpdb=$prefix.$i.'.pdb';   

    my $key=$designkey->[$i];
        print "prepare the self file for delphi ....  $newpdb \n";
    my $out=&GENERAL::GetOutFH("$workdir/$newpdb");
        my $pdb=&GENERAL::GetInFH("../$workpdb");   

    while (<$pdb>) {
        if (/^ATOM/) {
             my ($id, $name, $resnum) = &GENERAL::Trim(substr($_,21, 1), substr($_,17,4),substr($_, 22,4));
         my $key2= $id.$resnum.$name;
                 print $key2, "\n";
         if ( $key eq $key2 && ! / CA | N | H | C | O / )  {
                   $out->print($_);
         }  else {
           substr($_, 21, 1) =sprintf("%1s", 'X');
           $out->print($_);
                 }
             } elsif (/^TER/) {
                   $out->print($_);
             } else {
                   next;
             }
        }

        close $pdb;
        close $out;
}




=head2 

 Title   : WritePairPDB
 Usage   : &WritePairPDB($parm, \@solutionkeys, \@xraykeys, $l, $xbit); 
 Function: 
 Returns : 
 Args    : 

=cut



sub WritePairPDB {

    my $parm=shift;
    my $solutionkey=shift;
    my $xraykey=shift;
    my $i=shift;
    my $j=shift;
    my $n=shift;
    my $tbit =$parm->{'pair_tbit'};
    my $xbit=$parm->{'xbit'};

    my $key1;
    my $key2;

    if ($xbit == 0) {
    $key1=$solutionkey->[$i];
    $key2=$solutionkey->[$j];
    } else {
        $key1=$xraykey->[$i];
    $key2=$xraykey->[$j];
    }


    my $prefix;
    my $workdir;
    if ($tbit ==1) {
       $prefix='pair-';
       $workdir=$parm->{'pdb_pair_bb'};
    } else {
       $prefix ='rotpair-';
       $workdir=$parm->{'pdb_pair_rot'};
    }


    my $pdbfile= $parm->{'mixed_pdb_file'};
    my $pdb=&GENERAL::GetInFH("../$pdbfile");
    my $newpdb=$prefix.$n.'.pdb';
    my $out=&GENERAL::GetOutFH("$workdir/$newpdb");  
    printf "prepare file:....  %15s\n", $newpdb;
    my $charge=lc($parm->{'charge_pdb'}); 

    while (<$pdb>) {
       if (/^ATOM/) { 
            my ($id, $name, $resnum, $rotnum)
                = &GENERAL::Trim(substr($_,21, 1), substr($_,17,4),substr($_, 22,4), substr($_,26,2));

        my $key3= $id.$resnum.$name.$rotnum;
        if ($key3 eq $key1 ) {
        $out->print($_);
            } elsif ($key3 eq $key2) {
#                substr($_, 21, 1) =sprintf("%1s", 'X') if ($charge eq 'on');
#                substr($_, 26, 1) =sprintf("%1s", 'S') if ($surface eq 'on');
        $out->print($_);
        } elsif ( $tbit ==1 && substr($_, 27,1) eq " " ) {
                substr($_, 21, 1) =sprintf("%1s", 'X') if ($charge eq 'on');
                $out->print($_);
        } else {
        next;
        }

    } elsif (/^TER/) {
        $out->print($_);
    } else {
        next;
    }
   }
}




sub SplitWorkingPDB {

     my $fname=shift;
     my $ishash=shift;
     my $fh=&GENERAL::GetInFH($fname);
     my $tfh = &GENERAL::GetOutFH('template.pdb');
     my $rfh = &GENERAL::GetOutFH('rotamer.pdb');
  
     while ( <$fh>) {
          if ( /^ATOM/) {
              my @line = split(" ", $_);
              my $key= $line[4].$line[5].$line[3];
              if ( !exists $ishash->{$key} || /CA | N | O | H | C / ) {
                $tfh->print($_);
              } else {
                    $rfh->print($_);
              }
          }
     }
}         



sub RunSurface  {

    my $InputPDBFile = shift;
    my $parm=shift;
    my $surfaceexe=$parm->{'exec_parent_dir'}.$parm->{'surface_exe'};
    my $surface=$parm->{'exec_parent_dir'}.'/surface';

    my $AccallInput = &GENERAL::GetOutFH("accall.input");

    print $AccallInput <<"EOFAI";
PDBFILE $InputPDBFile
VDWFILE $surface/vdw.radii
STDFILE $surface/standard.data
PROBE 1.40
ZSLICE 0.05
EOFAI

    system("$surfaceexe < accall.input");

    unlink("accall.input") or die "can not delete file: accall.input \n" ;

}


sub ParseBackbone {
    my $fname=shift;
    my $outfh=shift;
    my $fh= &GENERAL::GetInFH($fname);

   while (<$fh>) {
        if ( /^RES/) {
           my ($non, $pol)=&GENERAL::Trim(substr($_,28,7 ), substr($_,41,7));
           $outfh->printf("%8.3f %8.3f\n", $non, $pol);
        }
   }
    unlink($fname);
}




sub ParseRotself {
    my $fname=shift;
    my $areafh=shift;
    my $fh= &GENERAL::GetInFH($fname);

   while (<$fh>) {
        if ( /^RES/) {
           my ($non, $pol)=&GENERAL::Trim(substr($_,28,7 ), substr($_,41,7));
           $areafh->printf("%8.3f %8.3f\n", $non, $pol);
        }
   }
    unlink($fname);
}


sub ParseTemplate {
    my $fname=shift;
    my $areafh=shift;
    my $fh= &GENERAL::GetInFH($fname);
    
   while (<$fh>) {
        if ( /^TOTAL/) {
           my @line=split(" ", $_);  
           my ($non, $pol)=($line[6], $line[7]);
           $areafh->printf("%8.3f %8.3f\n", $non, $pol);
        }
   }
    unlink($fname);
}



sub ParseSelf {
    my $fname=shift;
    # for final area calculations   
    my $areafh=shift;
    # for the asp calculations
    my $outfh=shift;

    # for the isolated rotamer    
    my $rot=shift;
    # for the empty template area
    my $bb=shift;
    # Asp parameters 
    my $par1=shift;

    # open the .rsa file
    my $fh= &GENERAL::GetInFH($fname);

   my ($t1, $t2);
   my $sum1=0.0;
   my $sumn1=0.0;
   my $sump1=0.0;
   my $sum2=0.0;
   my $sumn2=0.0;
   my $sump2=0.0;
   while (<$fh>) {
        # for the self part, one needs to count desolvation energy of both 
        # rotamer and backbone. desolvation energy of rotamer is isolate rotamer-
        # exposed area of rotamer in the rotamer+template
        # desolvation energy of backbone is template area (without any rotamer) -
        # exposed area of template in the rotamer+template 

        if ( /^RES/ && substr($_, 8,1) ne "X" ) {
           my ($non, $pol)=&GENERAL::Trim(substr($_,28,7 ), substr($_,41,7));
           print $non, "   ", $pol, "\n";
           $t1=$non;
           $t2=$pol;
           $non -=$rot->[0];
           $pol -=$rot->[1];
       # applying scaling factors
           $non *=0.62;
           $pol *=0.62;

           $sum1=$sum1+$non+$pol;
           $sumn1 +=$non;
           $sump1 +=$pol;
           print $non, "   ", $pol, "\n\n";
           $areafh->printf("%8.3f %8.3f ", $non, $pol);
           $non *=$par1->[0];
           $pol *=$par1->[1];
           $sumn2 +=$non;
           $sump2 +=$pol;
           $outfh->printf("%8.3f %8.3f ", $non, $pol);
           $sum2=$sum2+$non+$pol;
        } else {
        }   
 
        if ( /^TOTAL/) {
           my @line=split(" ", $_);
           my ($non, $pol)=($line[6], $line[7]);
           print $non, "   ", $pol, "\n";
           print $bb->[0], "   ", $bb->[1], "\n";

       # first get exposed area of template in the rotamer+template 
           $non -=$t1;
           $pol -=$t2;
   
           print $non, "   ", $pol, "\n";
       # then get the buried area of template in the rotamer+template
           $non -=$bb->[0];
           $pol -=$bb->[1];
           print $non, "   ", $pol, "\n";

       # then get the scaled buried areas
       $non *=0.62;
           $pol *=0.62;
           $sumn1 +=$non;
           $sump1 +=$pol;
           $sum1= $sum1+$non+$pol;
           $areafh->printf(" %9.3f %9.3f", $non, $pol);


           $non *=$par1->[0];
           $pol *=$par1->[1];
           $sum2=$sum2+$non+$pol;
           $sumn2 +=$non;
           $sump2 +=$pol;
           $outfh->printf(" %9.3f %9.3f", $non, $pol);
           $areafh->printf(" %9.3f %9.3f %9.3f\n",$sumn1, $sump1, $sum1);
           $outfh->printf(" %9.3f %9.3f %9.3f\n", $sumn2, $sump2, $sum2);
        } else {
        }
   } 

    unlink($fname);
}



sub ParsePair {
    my $fname=shift;
    my $areafh=shift;
    my $outfh=shift;
    my $roti=shift;
    my $rotj=shift;
    my $par1=shift;
    my $par2=shift;


    my $fh= &GENERAL::GetInFH($fname);

    my $counter=0;
    my $sum1=0.0;
    my $sumn1=0.0;
    my $sump1=0.0;
    my $sum2=0.0;
    my $sumn2=0.0;
    my $sump2=0.0;


   while (<$fh>) {
        if ( /^RES/) {
           $counter++;
           my ($non, $pol)=&GENERAL::Trim(substr($_,28,7 ), substr($_,41,7));
           my ($non1, $pol1);
           if ( $counter == 1) {
                 $non -=$roti->[0];
                 $pol -=$roti->[1];
         $non *=0.62;
                 $pol *=0.62;
                 $sum1=$sum1+$non+$pol;
                 $sumn1 +=$non;
                 $sump1 +=$pol;
                 $areafh->printf("%8.3f %8.3f ", $non, $pol);
                 $non1 = $non * $par1->[0];
                 $pol1 = $pol * $par1->[1];
                 $sum2 = $sum2+$non1+$pol1;
                 $sumn2 +=$non1;
                 $sump2 +=$pol1;
                 $outfh->printf("%8.3f %8.3f ", $non1, $pol1);
           } else {
                 $non -=$rotj->[0];
                 $pol -=$rotj->[1];
                 $non *=0.62;
                 $pol *=0.62;
                 $sumn1 +=$non;
                 $sump1 +=$pol;
                 $sum1= $sum1+$non+$pol;
                 $areafh->printf("%8.3f %8.3f ", $non, $pol);
                 $areafh->printf("%8.3f %8.3f %8.3f \n", $sumn1, $sump1, $sum1);
                 $non1 = $non * $par2->[0];
                 $pol1 = $pol * $par2->[1];
                 $sum2= $sum2+$non1+$pol1;
                 $sumn2 +=$non1;
                 $sump2 +=$pol1;
                 $outfh->printf("%8.3f %8.3f ", $non1, $pol1);
                 $outfh->printf("%8.3f %8.3f %8.3f \n", $sumn2, $sump2, $sum2);
      }

#      $areafh->printf("%8.3f %8.3f ", $non, $pol);
#           $outfh->printf("%8.3f %8.3f ", $non1, $pol1);
        }
   }

#   $areafh->printf("\n");
#   $outfh->printf("\n");
   unlink($fname);

}


sub WriteDelphiInp {
    my %arg=@_;
    my $fhandle=&GENERAL::GetOutFH('fort.10');

#   'gsize', 'default', 'focus', 'high', 'diel_in', 2, 'diel_out', 1

    if (lc($arg{'gsize'}) =~ /default/ ) {
        $fhandle->print('!gsize=65',  "\n");
    } else {
        $fhandle->print('gsize=',$arg{'gsize'}, "\n");
    }
   
    if (lc($arg{'focus'}) eq 'high') {
        $fhandle->print('perfil=92.',  "\n");
    $fhandle->print('bndcon=3',  "\n");
    } else {
        $fhandle->print('perfil=23',  "\n");
        $fhandle->print('bndcon=4',  "\n");
    }

    $fhandle->print('center(0.,0.,0)',  "\n");

    $fhandle->print('indi=', $arg{'diel_in'}, "\n");
    $fhandle->print('exdi=', $arg{'diel_out'}, "\n");

    if ( $arg{'diel_out'} != 1) { 
       # $fhandle->print("scale=4.0\n");
        $fhandle->print('salt=0.145',  "\n");
        $fhandle->print('ionrad=2',  "\n");
    }   

    $fhandle->print('nonit=0',  "\n");
    $fhandle->print("energy\(G\,C\,S\)",  "\n");
    $fhandle->print("out\(modpdb\)",  "\n");
    $fhandle->print("out\(frc\)",  "\n");

    if (lc($arg{'focus'}) eq 'high') {
        $fhandle->print('in', "\(phi, file=\'temp.phi\'\)",  "\n");
    } else {
        $fhandle->print('out', "\(phi, file=\'temp.phi\'\)",  "\n");
    }

    close  $fhandle; 
}


sub RunDeeSpec {
    
    my $e=shift;
    my $verbose=shift;
    $verbose=1 if (!defined $verbose);
    die "you must provide energy file name!\n" if (!defined $e); 
    my $design=DESIGN::new('dee-astar', 'dee-astar-spec.inp', $e);
    print "using combination of Dead-End Elimination and A-star Algorithm ...\n" if ($verbose==1);
    print "the input file is dee-astar-spec.inp ...\n" if ($verbose==1);
    print "create design input file......\n" if ($verbose==1);
    $design->writeDesignInput(1, 0.0001, 0, 'rot_subsets');
    print "run design......\n" if ($verbose==1);
    $design->runDesign('design');
    print "parse design output files...\n" if ($verbose==1);
    $design->parseOutput();
    return $design->{solution}->{results}->[0];
}


sub RunDee {

    my $e=shift;
    my $verbose=shift;
    $verbose=1 if (!defined $verbose);
    die "you must provide energy file name!\n" if (!defined $e);
    my $design=DESIGN::new('dee-astar', 'dee-astar.inp', $e);
    print "using combination of Dead-End Elimination and A-star Algorithm ...\n" if ($verbose==1);
    print "the input file is dee-astar.inp ...\n" if ($verbose==1);
    print "create design input file......\n" if ($verbose==1);
    $design->writeDesignInput(0, 0.0001, 0);
    print "run design......\n" if ($verbose==1);
    $design->runDesign('design');
    print "parse design output files...\n" if ($verbose==1);
    $design->parseOutput();
    return $design->{solution}->{results}->[0];
}


sub BuildSequence {

    my $rotamer=shift;
    my $parm=shift;
    my $res=shift;
    my @seq;

    foreach my $site (@{$rotamer->{dslist}}) {
    if ($site->{chain} eq $parm->{chain} && $site->{iresnum} eq $parm->{iresnum} ) {
        push (@seq, $res);
    } else {
        push (@seq, $site->{reslist}->[0]->{resname});
    }
    }
    
    return \@seq;
}

1;
