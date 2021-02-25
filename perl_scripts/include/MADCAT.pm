# By Gevorg Grigoryan

package MADCAT;
$DEF_VER_CONTROL->{date}->{MADCAT} = "04/12/2012";

use strict;
use GENERAL;
use PDB;
use Math::Libm ':all';

=head2

 Title   :  new
 Usage   :  my $madcat = MADCAT::new("niceRun")
            my $madcat = MADCAT::new('home', '/home/grigoryanlab/library/MaDCaT/', 'mapList', '/home/grigoryanlab/library/MaDCaT/lists/bc-30.list')
 Function:  MADCAT object constructor.
 Returns :  none
 Args    :  1. name of run
            *  optionalm 'home': root directory for desired MADCAT installation (should have bin/ and scripts/)
            *  optional, 'mapList': path to map list for use in searches 
            *  optional, 'keepOut': if set (default), will keep all of the disk-written output; otherwise, will read
               and return search result, destroying disk-written data 
            *  optional, 'wdir':    working directory for MaDCaT (defaults to current directory if keepOut is set, and a temporary directory otherwise)
            *  optional, 'silent':  if set (not by default), will hide MaDCaT-related stdout
            *  optional, 'profile': name of a canned profile for a run; sets many options, making much of the above unnecessary. Choices are:
              ** 'anthillIronFS' - MaDCaT will be run on the anthill cluster, using map a list referring to the network-stored maps on IronFS. The job will be
                 submitted only to nodes that can access this FS. In this case 'list' (if specified) gets interpreted as the name of the list, not the path.
                 If not specified, defaults to bc-30.list.
              ** 'anthillLocal'  - MaDCaT will be run on the anthill cluster, using map a list referring to locally-stored maps on the running node. The job will be
                 submitted only to nodes that are known to store maps locally. In this case 'list' (if specified) gets interpreted as the name of the list, not the path.
                 If not specified, defaults to bc-30.list.
              ** 'anthill'       - MaDCaT will be run on the anthill cluster, but the map list can be explicitly specified via a path.
              For all cluster profiles the property 'time' can be set (wall time in minutes); if not set, defaults to 2 hrs 59 minutes and 59 seconds.

=cut

sub new {
  my $self = {};
  my $name = shift;
  GENERAL::requireArgs($name);
  my %par = @_;
  bless $self;

  # set input options specifying general run information
  $self->{base} = $name;
  $self->{home} = defined($par{home}) ? $par{home} : "/home/grigoryanlab/library/MaDCaT/";
  $self->{cluster}->{run} = 0;
  if (defined($par{profile})) {
    if ($par{profile} =~ /\s*anthillIronFS\s*/) {
      $self->{cluster}->{run} = 1;
      $self->{cluster}->{server} = "anthill";
#      $self->{cluster}->{nodeList} = MADCAT::readNodes("/home/grigoryanlab/library/MaDCaT/lists/ironfs/nodes");
      push(@{$self->{cluster}->{attributes}}, "-l ironfs");
      $self->{cluster}->{time} = defined($par{time}) ? MADCAT::convertRunTime($par{time}) : "2:59:59";
      if (defined($par{list})) {
        $self->{list} = "/home/grigoryanlab/library/MaDCaT/lists/ironfs/$par{list}";
      } else {
        $self->{list} = "/home/grigoryanlab/library/MaDCaT/lists/ironfs/bc-30.list";
      }
    } elsif ($par{profile} =~ /\s*anthillLocal\s*/) {
      $self->{cluster}->{run} = 1;
      $self->{cluster}->{server} = "anthill";
      $self->{cluster}->{time} = defined($par{time}) ? MADCAT::convertRunTime($par{time}) : "2:59:59";
      $self->{cluster}->{nodeList} = MADCAT::readNodes("/home/grigoryanlab/library/MaDCaT/lists/nodeTemp/nodes");
      if (defined($par{list})) {
        $self->{list} = "/home/grigoryanlab/library/MaDCaT/lists/nodeTemp/$par{list}";
      } else {
        $self->{list} = "/home/grigoryanlab/library/MaDCaT/lists/nodeTemp/bc-30.list";
      }
    } elsif ($par{profile} =~ /\s*anthill\s*/) {
      $self->{cluster}->{run} = 1;
      $self->{cluster}->{server} = "anthill";
      $self->{cluster}->{time} = defined($par{time}) ? MADCAT::convertRunTime($par{time}) : "2:59:59";
      if (defined($par{list})) {
        $self->{list} = $par{list};
      } else {
        $self->{list} = "/home/grigoryanlab/library/MaDCaT/lists/bc-30.list";
      }
    } else {
      GENERAL::error("Unknown profile name '$par{profile}");
    }
  } else {
    $self->{list} = defined($par{list}) ? $par{list} : "/home/grigoryanlab/library/MaDCaT/lists/bc-30.list";
  }
  $self->{keepOut} = defined($par{keepOut}) ? $par{keepOut} : 1;
  if (defined($par{wdir})) {
    $self->{wdir} = $par{wdir};
  } else {
    if (($self->{keepOut} == 0) && ($self->{cluster}->{run} == 0)) {
      $self->{wdir} = GENERAL::createLocalSpace();
      $self->{wdirIsTemp} = 1;
    } else {
      $self->{wdir} = "./";
    }
  }
  $self->{silent} = (defined($par{silent}) && $par{silent}) ? " > /dev/null": "";
  GENERAL::assert((-e $self->{list} ? 1 : 0), "MaDCaT plans to work on list $self->{list}, but it does not exist!");

  # set default detailed MaDCaT options
  $self->setDefaultOptions();

  # set the binary paths
  $self->{bins}->{createDM} = $self->{home} . "/bin/createDM";
  $self->{bins}->{madcat} = $self->{home} . "/bin/madcat";
  $self->{bins}->{moreInfo} = $self->{home} . "/bin/structsFromMatches";

  # initialize hashes of created files and directories
  $self->{generatedFiles} = ();
  $self->{generatedDirs} = ();

  return $self;
}

DESTROY {
  my $self = shift;

#  $self->finish() if (!defined $self->{_finished} || !$self->{_finished});
}


=head2

 Title   :  finish
 Usage   :  $madcat->finish();
 Function:  Cleans up the MaDCaT run.
 Returns :  Nothing
 Args    :  None.

=cut

sub finish {
  my $self = shift;
  
  if (defined($self->{wdirIsTemp})) {
    GENERAL::destroyLocalSpace($self->{wdir});
  } else {
    foreach my $file (keys(%{$self->{generatedFiles}})) {
      if (($self->{generatedFiles}->{$file} > $self->{keepOut}) && (-e "$self->{wdir}/$file")) {
        GENERAL::crm("$self->{wdir}/$file");
      }
    }
    foreach my $dir (keys(%{$self->{generatedDirs}})) {
      if (($self->{generatedDirs}->{$dir} > $self->{keepOut}) && (-d "$self->{wdir}/$dir")) {
        GENERAL::csystem("rm -r $self->{wdir}/$dir");
      }
    }
  }
  $self->{_finished} = 1;
}

sub setDefaultOptions {
  my $self = shift;
  my %opts;
  $self->{opts}->{diag} = 1;
  $self->{opts}->{topN} = 5000;
  $self->{opts}->{matchOut} = $self->{base} . ".match";
  $self->{opts}->{seqOut} = $self->{base} . ".seq";
  $self->{opts}->{structOut} = "";
  $self->{opts}->{structOutType} = "match";
  $self->{opts}->{linkerSeq} = "";
  $self->{opts}->{linkerStruct} = "";
  $self->{opts}->{cutsc} = "";
  $self->{opts}->{greed} = "";
  $self->{opts}->{cutsz} = "";
}

sub setOption {
  my $self = shift;
  my $name = shift;
  my $value = shift;
  GENERAL::warning("unknown option name '$name'") if (!defined($self->{opts}->{$name}));
  $self->{opts}->{$name} = $value;
}

sub setClusterOption {
  my $self = shift;
  my $name = shift;
  my $value = shift;
  GENERAL::warning("unknown option name '$name'") if (!defined($self->{cluster}->{$name}));
  $self->{cluster}->{$name} = $value;
}

=head2

 Title   :  runMaDCaT
 Usage   :  $madcat->runMaDCaT('ipdb', "frag.pdb");
 Function:  Runs MaDCaT on the given input.
 Returns :  Nothing
 Args    :  Argument name/value pairs:
            * 'ipdb' or 'imap' -- one of these required (but can't use both!)

=cut

sub runMaDCaT {
  my $self = shift;
  GENERAL::assert(scalar(@_) % 2 == 0, "Expected an even numer of inputs for name/value pairs!");
  my %runOpts = @_;
  GENERAL::assert(defined($runOpts{ipdb}) || defined($runOpts{imap}), "must define either an query structure or map!");

  # prepare query map
  my $imap;
  if (defined($runOpts{imap})) {
    $imap = $runOpts{imap};
  } else {
    # create the map
    $self->{opts}->{ipdb} = $runOpts{ipdb};
    $imap = $self->{base} . ".map";
    GENERAL::csystem("cd $self->{wdir}; $self->{bins}->{createDM} --p $self->{opts}->{ipdb} --o $imap $self->{silent}");
    $self->{generatedFiles}->{$imap} = 1;
  }
  $self->{opts}->{imap} = $imap;

  # run MaDCaT
  my $cmd = sprintf("$self->{bins}->{madcat} --map $self->{opts}->{imap} --compList $self->{list} --diag %d --topN %d --matchOut %s --seqOut %s", $self->{opts}->{diag}, $self->{opts}->{topN}, $self->{opts}->{matchOut}, $self->{opts}->{seqOut});
  $self->{generatedFiles}->{$self->{opts}->{matchOut}} = 1;
  $self->{generatedFiles}->{$self->{opts}->{seqOut}} = 1;

  if ($self->{opts}->{structOutType} ne "") {
    $cmd .= " --structOutType $self->{opts}->{structOutType}";
  }
  if ($self->{opts}->{structOut} ne "") {
    $cmd .= " --structOut $self->{opts}->{structOut}";
    if ($self->{opts}->{structOutType} =~ /file/) {
      $self->{generatedFiles}->{$self->{opts}->{structOut}} = 1;
    } else {
      $self->{generatedDirs}->{$self->{opts}->{structOut}} = 1;
    }
  }
  if ($self->{opts}->{linkerSeq} ne "") {
    $cmd .= " --linkerSeq $self->{opts}->{linkerSeq}";
    $self->{generatedFiles}->{$self->{opts}->{linkerSeq}} = 1 if ($self->{opts}->{linkerSeq} !~ /^\s*together\s*$/);
  }
  if ($self->{opts}->{linkerStruct} ne "") {
    $cmd .= " --linkerStruct $self->{opts}->{linkerStruct}";
    $self->{generatedDirs}->{$self->{opts}->{linkerStruct}} = 1 if ($self->{opts}->{linkerStruct} !~ /^\s*together\s*$/);
  }
  if ($self->{opts}->{cutsc} ne "") {
    $cmd .= " --cutsc $self->{opts}->{cutsc}";
  }
  if ($self->{opts}->{greed} ne "") {
    $cmd .= " --greed $self->{opts}->{greed}";
  }
  if ($self->{opts}->{cutsz} ne "") {
    $cmd .= " --cutsz $self->{opts}->{cutsz}";
  }
  if ($self->{cluster}->{run}) {
    $self->{cluster}->{attempt} = 0 if (!defined($self->{cluster}->{attempt})); 
    if ($self->{cluster}->{server} eq "anthill") {
      my $odir = GENERAL::GetDir();
      GENERAL::cchdir($self->{wdir});
      my $shf = "$self->{base}.sh";
      my $ofh = GENERAL::GetOutFH($shf);
      $ofh->printf("#!/bin/bash\n#\$ -j y\n#\$ -cwd\n#\$ -V\n");
      $ofh->printf("#\$ -l hostname=\'$self->{cluster}->{nodeList}\'\n") if (defined($self->{cluster}->{nodeList}));
      $ofh->printf("#\$ -l h_rt=$self->{cluster}->{time}\n");
      $ofh->printf("#\$ -l vf=2G\n\n");
      if (defined($self->{cluster}->{attributes})) {
        foreach my $attr (@{$self->{cluster}->{attributes}}) {
          $ofh->printf("#\$ $attr\n\n");
        }
      }
      $ofh->printf("$cmd $self->{silent}\n\n");
      close($ofh);
      $self->{generatedFiles}->{$shf} = 1;
      my $out = `qsub $shf`;
      if ($out !~ /Your job (\d+) .+ has been submitted/) {
        GENERAL::error("Could not submit MaDCaT job to cluster - got message '$out'")
      }
      $self->{cluster}->{jobId} = $1;
      $self->{generatedFiles}->{"$shf.o$1"} = 2;
      GENERAL::cchdir($odir);
    } else {
      GENERAL::error("Don't know how to run jobs on cluster '$self->{cluster}->{server}'")
    }
    $self->{cluster}->{attempt}++;
  } else {
    GENERAL::csystem("cd $self->{wdir}; $cmd $self->{silent}");
  }
}


=head2

 Title   :  waitForJob
 Usage   :  $madcat->waitForJob();
 Function:  Block waits until the submitted job finishes.
 Returns :  Did the job finish (1) or did it timeout (0).
 Args    :  1. optional: interval in seconds between cluster polls (every 10 seconds by default).
            2. optional: timeout in minutes. Will return 0 if timeout expires without job finishing.
               By default, there is no timeout. 

=cut

sub waitForJob {
  my $self = shift;
  my $inter = shift;
  my $timeout = shift;
  $inter = 10 if (!defined($inter));
  my $t = time();
  while (1) {
    last if ($self->isJobFinished() || (defined($timeout) && ((time() - $t)/60 > $timeout)));
    sleep($inter);
  }
}

=head2

 Title   :  isJobFinished
 Usage   :  $madcat->isJobFinished();
 Function:  Check whether cluster-submitted job has finished.
 Returns :  1 if finished, 0 otherwise.
 Args    :  None 

=cut

sub isJobFinished {
  my $self = shift;

  GENERAL::assert(defined($self->{cluster}->{jobId}), "it seems no job was submitted to the cluster!");
  my $out = `qstat | grep "^$self->{cluster}->{jobId}"`;
  $out = GENERAL::Trim($out);
  return ($out eq "");
}

=head2

 Title   :  readMatchInfo
 Usage   :  $madcat->readMatchInfo();
 Function:  Reads the sequences of the matches from the last run.
 Returns :  Array of match seqs and RMSDs/map scores.
 Args    :  Pairs of name/value pairs. Possibilities:
            * 'readMatchFile' - 0 (default) for not reading the match file, 1 - for reading it

=cut

sub readMatchInfo {
  my $self = shift;
  GENERAL::assert(scalar(@_) % 2 == 0, "expected an even numer of inputs for name/value pairs!");
  my %opts = @_;
  $opts{readMatchFile} = 0 if (!defined($opts{readMatchFile}));

  # read sequence file  
  GENERAL::assert((-e "$self->{wdir}/$self->{opts}->{seqOut}" ? 1 : 0), "could not find output file $self->{opts}->{seqOut} -- perhaps MaDCaT did not run yet?");
  my $lines = GENERAL::file2array("$self->{wdir}/$self->{opts}->{seqOut}");
  my $mlines = $opts{readMatchFile} ? GENERAL::file2array("$self->{wdir}/$self->{opts}->{matchOut}") : undef;
  my (@seqs, @seqStrs, @rmsds, @mapSc, @matches, @matchFiles, @matchNames, @matchBases);
  for (my $i = 0; $i < scalar(@$lines); $i++) {
    my $line = $lines->[$i];
    my @line = split(" ", $line);
    GENERAL::assert(scalar(@line) >= 3, "Line '$line' of sequence file could not be parsed!");
    my $rmsd = shift @line;
    $rmsd = 0 if ($rmsd =~ /nan/);
    push(@rmsds, $rmsd);
    push(@mapSc, shift @line);
    push(@seqs, \@line);
    push(@seqStrs, join(" ", @line));
    if (defined($mlines)) {
      my $match = $mlines->[$i];
      push(@matches, $match);
      my @match = split(" ", $match);
      GENERAL::assert(scalar(@match) > 3, "could not parse match '$match'");
      my $file = $match[2]; $file =~ s/'//g;
      push(@matchFiles, $file);
      $file =~ s/^.*\/([^\/]+)$/$1/;
      push(@matchNames, $file);
      $file =~ s/^([^\.]+)\..*$/$1/;
      push(@matchBases, $file);
    }
  }
  my %res;
  $res{seqs} = \@seqs;
  $res{seqStrings} = \@seqStrs;
  $res{RMSDs} = \@rmsds;
  $res{mapScores} = \@mapSc;
  if ($opts{readMatchFile}) {
    $res{matches} = \@matches;
    $res{matchFiles} = \@matchFiles;
    $res{matchNames} = \@matchNames;
    $res{matchBases} = \@matchBases;
  }
  return \%res;
}

=head2

 Title   :  contextProbability
 Usage   :  $madcat->contextProbability();
 Function:  Assesses the chances of the context in which the given query motif is encountered.
 Returns :  A probability-like metric.
 Args    :  1. PDB object containing the query structure. 
            2. PDB object with the query plus its context (the original structure from which the query came).
            * name/value combinations. Choices are:
            ** 'pos' -- a subset of query positions to consider the context for. By default, will do all.
            ** 'matchFile' -- a file with matches to the query that will be considered in assessing the proability
               of the original context. By default, all matches from the previous run will be used.
            ** 'imap' -- input map that was used for the search. By default, will use the map used for the last
               search (this option is useful when there was no explicit search in this object).
=cut

sub contextProbability {
  my $self = shift;
  my $qpdb = shift;
  my $opdb = shift;
  GENERAL::requireArgs($self, $qpdb, $opdb);
  GENERAL::assert(scalar(@_) % 2 == 0, "expected an even numer of inputs for name/value pairs!");
  my %opts = @_;
  if (!defined($opts{pos})) {
    my @pos = GENERAL::linspace(1, scalar($qpdb->ConRes()), 1);
    $opts{pos} = \@pos;
  }
  GENERAL::assert(defined($opts{imap}) || defined($self->{opts}->{imap}), "no input map specified, and there appears to have been no explicit search prior to the call!");
  my $imap = defined($opts{imap}) ? $opts{imap} : $self->{opts}->{imap};

  # create temporary directory for dumping structures
  my $ldir = GENERAL::createLocalSpace(undef, 0);
  my $mdir = "$ldir/matchStructs";
  my $fdir = "$ldir/fullStructs";

  # dump full and match-region structures for all matches
  my $matchFile = defined($opts{matchFile}) ? $opts{matchFile} : $self->{opts}->{matchOut};
  GENERAL::csystem("$self->{bins}->{moreInfo} --map $imap --diag $self->{opts}->{diag} --matchOut $matchFile --structOut $mdir --structOutType match $self->{silent}");
  GENERAL::csystem("$self->{bins}->{moreInfo} --map $imap --diag $self->{opts}->{diag} --matchOut $matchFile --structOut $fdir --structOutType full $self->{silent}");

  # read each match and characterize the environment around the residues in question
  my $n = GENERAL::GetFileInfo($matchFile, 'lines');
  my $qlen = scalar($qpdb->ConRes());
  my @C;
  for (my $i = 1; $i <= $n; $i++) {
    my $matchpdb = PDB::new(sprintf("$mdir/match%05d.pdb", $i), "PDB", undef, undef, 'norename', 1);
    next if (scalar($matchpdb->ConRes()) != $qlen);
    push(@C, MADCAT::matchContext($matchpdb, PDB::new(sprintf("$fdir/match%05d.pdb", $i), "PDB", undef, undef, 'norename', 1), $opts{pos}));
  }
  my $mC = GENERAL::mean(\@C);
  my $sC = GENERAL::stdev(\@C);

  # compare to the environment of the query and compute probability
  my $c = MADCAT::matchContext($qpdb, $opdb, $opts{pos});
  my $z = abs($c - $mC)/($sC);
  my $p = 1 - erf($z/sqrt(2));

  # clean up
  GENERAL::destroyLocalSpace($ldir);

  return $p;
}


=head2

 Title   :  matchContext
 Usage   :  MADCAT::matchContext();
 Function:  Characterizes the context of a given match.
 Returns :  The number of CA atoms (in the full match structure) within 12 A of CA atoms of specified match residues.
 Args    :  1. PDB object corresponding to the matching region.
            2. PDB object corresponding to the full entry of the match.
            3. Reference to an array of match positions positions (1-initiated indices) to consider.
=cut

sub matchContext {
  my $mpdb = shift;
  my $fpdb = shift;
  my $pos = shift;
  my $dcut = 12;

  # calculate number of CA atoms within the cutoff distance of the given residues, not counting
  # residues within the matching regions
  my @allRes = $mpdb->ConRes();
  my (@mres, %mres);
  foreach my $pi (@$pos) {
    GENERAL::assert($pi <= scalar(@allRes), sprintf("position $pi requested from match that has only %d residues!", scalar(@allRes)));
    my $res = $allRes[$pi-1];
    my $mres = $fpdb->getResByInd($res->{chain}->{id}, $res->{resnum}, 1, 3);
    push(@mres, $mres);
    $mres{$mres} = 1;
  }
  my $n = 0;
  foreach my $mres (@mres) {
    my @ca = (PDB::getAtomInRes($mres, "CA", -1));
    @ca = (PDB::_newAtom(1, 1, "CA", PDB::centerOfMass($mres->{atom}), 1, 1, 0)) if ($ca[0] eq -1);
    my $atoms = $fpdb->atomsWithin(\@ca, $dcut, 0);
    foreach my $a (@$atoms) {
      next if ($a->{atomname} ne "CA");
      next if (defined($mres{$a->{residue}}));
      $n++;
    }
  }
  return $n/scalar(@mres);
}

=head2

 Title   :  mostInformativeMatchSubset
 Usage   :  $madcat->mostInformativeMatchSubset();
 Function:  Reads the sequences of the matches from the last run.
 Returns :  Array of match seqs and RMSDs/map scores.
 Args    :  1. summary of search results - the output from MADCAT::readMatchInfo
            * name/value combinations. Choices are:
            ** 'pos' -- a subset of alignment positions to consider
            ** 'width' -- the length that match sequences are expected to have. If specified,
               matches with sequences of a different length will be discarded.
            ** 'id'    -- percent identity to use as the sequence redundancy filter (100% by default)
            ** 'writeUnique' -- specifies the base name of .seq and .match files to write with
               just the "unique" sequences considered from all of the matches (using whatever the
               percent identity cutoff). By default, these files are not written.
            ** 'writeUsed' -- specifies the base name of .seq and .match files to write with
               just the sequences comprising the most informative sub-alignment. By default,
               these files are not written.
            ** 'writeRedundant' -- specifies the base name for .seq and .match files that include
               all hits (redundancy not removed) with RMSDs below the worst RMSD included in the
               most informative subset.
            ** 'nativeSeq' -- gives the "native" sequence of the query (as a space-separated, three-letter code
               string), such that matches with sequence too close to it will be discarded (too close is
               defined via 'id')
            ** 'nativeID' -- the PDB ID of the entry from which the query originated; matches
               from this entry will be discarded.
            ** 'onePerID' -- if set to 1 (default is 0), will allow at most one match per PDB ID (will
               take the best)

=cut

sub mostInformativeMatchSubset {
  my $srch = shift;
  GENERAL::requireArgs($srch);
  GENERAL::assert(scalar(@_) % 2 == 0, "expected an even numer of inputs for name/value pairs!");
  GENERAL::assert(scalar(@{$srch->{seqs}}) > 0, "empty search results passed!");
  my ($sofh, $mofh, $rsofh, $rmofh);

  my %opts = @_;
  if (defined($opts{id})) {
    GENERAL::assert(($opts{id} > 0) && ($opts{id} <= 100), "redundancy ID cutoff should be between 0 and 100");
  } else { $opts{id} = 100; }
  if (!defined($opts{pos})) {
    my @pos;
    for (my $i = 0; $i < scalar(@{$srch->{seqs}->[0]}); $i++) {
      push(@pos, $i+1);
    }
    $opts{pos} = \@pos;
  }
  $opts{onePerID} = 0 if (!defined($opts{onePerID}));
  GENERAL::assert(defined($srch->{matches}), "search results object must contain match-file lines to enable exclusion based on database IDs") if (defined($opts{nativeID}) || $opts{onePerID});

  # mark as invalid all sequences not of the same length as needed
  my @ok = GENERAL::ones(scalar(@{$srch->{seqs}}));
  if (defined($opts{width})) {
    for (my $i = 0; $i < scalar(@{$srch->{seqs}}); $i++) {
      if (scalar(@{$srch->{seqs}->[$i]}) != $opts{width}) {
        $ok[$i] = 0;
      }
    }
  }
  GENERAL::assert(GENERAL::sum(@ok) > 0, "no match sequences had the requested length!");

  # first, sort all matches by RMSD
  my @sinds = GENERAL::linspace(0, scalar(@{$srch->{seqs}})-1, 1);
  @sinds = sort { $srch->{RMSDs}->[$a] <=> $srch->{RMSDs}->[$b] } @sinds;

  # now build a new array of sequence in order of low-to-high RMSD, removing redundancy
  my (@seqs, @inds, %h, %hID);
  if (defined($opts{nativeSeq})) {
    $opts{nativeSeq} = GENERAL::Trim($opts{nativeSeq});
    $h{$opts{nativeSeq}} = 1;
    my @tmp = split(" ", $opts{nativeSeq});
    $opts{nativeSeq} = \@tmp;
  }
  for (my $i = 0; $i < scalar(@sinds); $i++) {
    my $ii = $sinds[$i];
    next if ($ok[$ii] == 0);
    next if (defined($opts{nativeID}) && ($srch->{matchBases}->[$ii] =~ /$opts{nativeID}/i));
    next if ($opts{onePerID} && defined($hID{$srch->{matchBases}->[$ii]}));

    # is the new sequence within the given cutoff of any of the previously admitted ones? 
    if ($opts{id} == 100) {
      next if (defined($h{$srch->{seqStrings}->[$ii]}));
      $h{$srch->{seqStrings}->[$ii]} = 1;
    } else {
      my $skip = 0;
      foreach my $seq (defined($opts{nativeSeq}) ? ($opts{nativeSeq}, @seqs) : @seqs) {
        my $n = 0;
        for (my $k = 0; $k < scalar(@$seq); $k++) {
          $n += ($seq->[$k] eq $srch->{seqs}->[$ii]->[$k]);
        }
        if ($n > $opts{id}*scalar(@{$srch->{seqs}->[$ii]})/100) {
          $skip = 1; last;
        }
      }
      next if ($skip);
    }

    # if no, admit this one
    push(@seqs, $srch->{seqs}->[$ii]);
    push(@inds, $ii);
    $hID{$srch->{matchBases}->[$ii]} = 1;
  }

  # now add the collected non-redundant sequence into the alignment in order of low-to-high RMSD and compute information content
  my (@I, @subAlgn, %prevSaved);
  if (defined($opts{writeUnique})) { $sofh = GENERAL::GetOutFH($opts{writeUnique} . ".seq"); $mofh = GENERAL::GetOutFH($opts{writeUnique} . ".match"); }
  for (my $i = 0; $i < scalar(@seqs); $i++) {
    push(@subAlgn, $seqs[$i]);
#    $I[$i] = GENERAL::informationContent(\@subAlgn, 20, 'pos', $opts{pos}, 'lowCountCorrection', 2, 'prevSaved', \%prevSaved);
#    $I[$i] = 0 if ($I[$i] < 0);
    my ($I, $Is) = GENERAL::informationContent(\@subAlgn, 20, 'pos', $opts{pos}, 'lowCountCorrection', 2, 'prevSaved', \%prevSaved);
    $I[$i] = GENERAL::max(0, $I - 1*$Is); # if the information content is not at least one standard deviation away from random expectation, consider it zero
    $I[$i] = ($i < 5 ? 0 : $I[$i]);
    # if asked, write files with unique hits
    if (defined($opts{writeUnique})) {
      my $ii = $inds[$i];
      # in this sequence file, instead of map scores write the information content with up to this sequence
      $sofh->printf("%f %f %s\n", $srch->{RMSDs}->[$ii], $I[$i], $srch->{seqStrings}->[$ii]);
      $mofh->printf("%s\n", $srch->{matches}->[$ii]);
    }
  }
  if (defined($opts{writeUnique})) { close($sofh); close($mofh); }

  # find where the information content peaks - that constitutes the most informative sub-alignment
  my $mi; GENERAL::max(\@I, \$mi);

  # package, calculate stats, and return
  my %res; my $k = 0;
  if (defined($opts{writeUsed})) { $sofh = GENERAL::GetOutFH($opts{writeUsed} . ".seq"); $mofh = GENERAL::GetOutFH($opts{writeUsed} . ".match"); } 
  if (defined($opts{writeRedundant})) { $rsofh = GENERAL::GetOutFH($opts{writeRedundant} . ".seq"); $rmofh = GENERAL::GetOutFH($opts{writeRedundant} . ".match"); } 
  for (my $i = 0; $i <= $mi; $i++) {
    my $ii = $inds[$i];
    push(@{$res{seqs}}, $srch->{seqs}->[$ii]);
    push(@{$res{seqStrings}}, $srch->{seqStrings}->[$ii]);
    push(@{$res{RMSDs}}, $srch->{RMSDs}->[$ii]);
    push(@{$res{mapScores}}, $srch->{mapScores}->[$ii]);
    push(@{$res{matches}}, $srch->{matches}->[$ii]) if (defined($srch->{matches}));
    # if asked, write sub-alignment files
    if (defined($opts{writeUsed})) {
      $sofh->printf("%f %f %s\n", $srch->{RMSDs}->[$ii], $srch->{mapScores}->[$ii], $srch->{seqStrings}->[$ii]);
      $mofh->printf("%s\n", $srch->{matches}->[$ii]);
    }
    if (defined($opts{writeRedundant})) {
      for (; $k < scalar(@sinds); $k++) {
        my $jj = $sinds[$k];
        if ($srch->{RMSDs}->[$jj] <= $srch->{RMSDs}->[$ii]) {
          $rsofh->printf("%f %f %s\n", $srch->{RMSDs}->[$jj], $srch->{mapScores}->[$jj], $srch->{seqStrings}->[$jj]);
          $rmofh->printf("%s\n", $srch->{matches}->[$jj]);
        } else { last; }
      }
    }
  }
  if (defined($opts{writeUsed})) { close($sofh); close($mofh); } 
  $res{posI} = GENERAL::ones(scalar(@{$res{seqs}->[0]}), 0, 1);
  $res{posIstd} = GENERAL::ones(scalar(@{$res{seqs}->[0]}), 0, 1);
  $res{prob} = GENERAL::ones(scalar(@{$res{seqs}->[0]}), 0, 1);
  $res{count} = GENERAL::ones(scalar(@{$res{seqs}->[0]}), 0, 1);
  for (my $i = 0; $i < scalar(@{$res{seqs}->[0]}); $i++) {
    my @pos = ($i+1);
    ($res{posI}->[$i], $res{posIstd}->[$i]) = GENERAL::informationContent($res{seqs}, 20, 'pos', \@pos, 'lowCountCorrection', 2);
    my %tmp; $res{prob}->[$i] = \%tmp;
    my %tmp1; $res{count}->[$i] = \%tmp1;
    for (my $j = 0; $j < scalar(@{$res{seqs}}); $j++) {
      my $lett = $res{seqs}->[$j]->[$i];
      $res{count}->[$i]->{$lett} = 0 if (!defined($res{count}->[$i]->{$lett}));
      $res{count}->[$i]->{$lett}++;
    }
    foreach my $lett (keys(%{$res{count}->[$i]})) {
      $res{prob}->[$i]->{$lett} = $res{count}->[$i]->{$lett}/scalar(@{$res{seqs}});
    }
  }

  # estimate motif abundance based on the initial slope of the log of CDF of RMSD
  my $ni = 100;
  if (scalar(@inds) > $ni) {
    $res{abundance} = (scalar(@{$srch->{seqs}->[0]}) - 1)*log($ni)/($srch->{RMSDs}->[$inds[$ni]] - $srch->{RMSDs}->[$inds[0]]);
  } else {
    $res{abundance} = "unknown";
  }

  return \%res;
}


=head2

 Title   :  wasRunSuccessful
 Usage   :  $madcat->wasRunSuccessful();
 Function:  Determines whether the last run was successful or not.
 Returns :  1/0 for successful/not.
 Args    :  None

=cut

sub wasRunSuccessful {
  my $self = shift;

  # check sequence and match files
  return 0 if (! -e "$self->{wdir}/$self->{opts}->{seqOut}");
  return 0 if (GENERAL::GetFileInfo("$self->{wdir}/$self->{opts}->{seqOut}", 'lines') != $self->{opts}->{topN});
  return 0 if (! -e "$self->{wdir}/$self->{opts}->{matchOut}");
  return 0 if (GENERAL::GetFileInfo("$self->{wdir}/$self->{opts}->{matchOut}", 'lines') != $self->{opts}->{topN});

  # check structure out files, if those were requested
  if ($self->{opts}->{structOut} ne "") {
    return 0 if ((($self->{opts}->{structOutType} =~ /file/) && (! -f $self->{opts}->{structOut})) ||
                (($self->{opts}->{structOutType} !~ /file/) && (! -d $self->{opts}->{structOut})));
  }

  # check linker sequence and structure files, if those were requested
  if (($self->{opts}->{linkerSeq} ne "") && ($self->{opts}->{linkerSeq} !~ /together/)) {
    return 0 if (! -f $self->{opts}->{linkerSeq});
  }
  if (($self->{opts}->{linkerStruct} ne "") && ($self->{opts}->{linkerStruct} !~ /together/)) {
    return 0 if (! -d $self->{opts}->{linkerStruct});
  }
  
  return 1;
}

=head2

 Title   :  readNodes
 Usage   :  $madcat->readNodes("nodes");
 Function:  Reads the list of allowed cluster nodes.
 Returns :  Reference to an array of allowed cluster nodes.
 Args    :  None

=cut

sub readNodes {
  my $nlf = shift;
  my $ifh = GENERAL::GetInFH($nlf);
  my @nodes;
  foreach my $line (<$ifh>) {
    $line =~ s/\#.*$//g;
    $line = GENERAL::Trim($line);
    next if ($line eq "");
    push(@nodes, $line);
  }
  return join("|", @nodes);
}


=head2

 Title   :  convertRunTime
 Usage   :  $madcat->convertRunTime(179);
 Function:  Converts time in minutes to a property-formatted wall-time string for the cluster.
 Returns :  Wall-time string.
 Args    :  None

=cut

sub convertRunTime {
  my $rem = shift;
  my $hrs = GENERAL::floor($rem/60);
  $rem = $rem - $hrs*60;
  my $mins = GENERAL::floor($rem);
  $rem = $rem - $mins;
  my $secs = GENERAL::floor($rem*60);
  return "$hrs:$mins:$secs";
}
1;
