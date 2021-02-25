# module for submitting jobs to a cluster
# author: Gevorg Grigoryan
# start date: April 15, 2014

package CLUSTER;

use GENERAL;

=head2

 Title   :  new
 Usage   :  my $c = CLUSTER::new("anthill", "discovery")
 Function:  CLUSTER object constructor.
 Returns :  new CLUSTER object.
 Args    :  1-n. names of one or more cluster to use

=cut

sub new {
  my $self = {};
  if (scalar(@_) == 0) { @_ = ("anthill"); }
  bless $self;

  foreach my $cname (@_) {
    GENERAL::error("don't know how to use cluster '$cname'") if ($cname !~ /^(anthill|discovery)$/);
    my %c;
    $c{name} = $cname;
    $c{maxPause} = 60;
    $c{minPause} = 0.01;
    $c{pauses} = ();
    if ($cname eq "anthill") {
      if (-e "/data/scratch") { $c{tmpdir} = "/data/scratch/" . GENERAL::GetUser() . "/jobscripts"; }
      else { $c{tmpdir} = "/tmp/" . GENERAL::GetUser() . "/jobscripts"; }
      GENERAL::cmkdir($c{tmpdir}, 1);
    } else {
      $c{tmpdir} = "/tmp/jobscripts";
      GENERAL::cmkdir($c{tmpdir}, 1);
    }

    push(@{$self->{clusters}}, \%c);
    $self->{clustByName}{$cname} = \%c;
    push(@{$self->{clustNames}}, $cname);
  }

  return $self;
}

=head2

 Title   :  submit
 Usage   :  $c->submit("cluster" => "anthill", "cmd" => "ls -l", "time" => "1:00:00");
            $c->submit("cmd" => "ls -l", "time" => "1:00:00"); # when dealing with only one cluster
 Function:  Submits a job consisting of the given command(s) to the cluster
 Returns :  nothing
 Args    :  1. command
            2. various attributes

=cut

sub submit {
  my $self = shift;
  my %opts = @_;
  if (!(defined($opts{'cmd'}) && defined($opts{'time'}))) {
    GENERAL::error("options 'cmd' and 'time' must be specified at a minimum");
  }
  GENERAL::assert(defined($opts{cluster}), "need to specify which cluster waiting on") if (scalar(@{$self->{clusters}}) > 1);

  my $clusterName = $self->{clustNames}->[0];
  $clusterName = $opts{cluster} if (defined($opts{cluster}));
  my $cluster = $self->{clustByName}->{$clusterName};
  $opts{m} = 1.0 if (!defined($opts{m}));
  GENERAL::assert(GENERAL::isNumeric($opts{m}) && ($opts{m} > 0), "memory does not make sense!");

  my $sname;
  if (defined($opts{name})) {
    $sname = "$opts{name}.sh";
  } else {
    do {
      $sname = "$cluster->{tmpdir}/GG$$." . GENERAL::irand(0, 10000) . ".sh";
      # very infrequently, clean old script files
      if (rand() < 0.001) {
        GENERAL::csystem("find $cluster->{tmpdir}/*.sh -mtime +7 -exec rm {} \\;");
      }
    } while (-e $sname);
  }
  my $ofh = GENERAL::GetOutFH($sname);
  if ($clusterName eq "anthill") {
    $ofh->print("#!/bin/bash\n");
    $ofh->print("#\$ -j y\n");     # combine stdout/stderr
    $ofh->print("#\$ -cwd\n");     # start job in current directory
    $ofh->print("#\$ -V\n");       # pass along environmental variables
    $ofh->print("#\$ -l vf=$opts{m}G\n");       # pass along environmental variables
    $ofh->print("#\$ -l h_rt=$opts{time}\n");
    $ofh->print("#\$ -N $opts{n}\n") if (defined($opts{n}));
    $ofh->print("#\$ -l hostname='tanto-0-*|katana-0-*|gridiron-*'\n") if (defined($opts{fast}) && $opts{fast});
    if (defined($opts{attr})) {
      if (ref($opts{attr}) eq "ARRAY") {
        foreach my $attr (@{$opts{attr}}) { $ofh->print("#\$ $attr\n"); }
      } else {
        $ofh->print("#\$ $opts{attr}\n");
      }
    }
    $ofh->print("#\$ -q '$opts{q}'\n") if (defined($opts{"q"}));
    $ofh->print("#\$ -pe $opts{pe}\n") if (defined($opts{pe}));
    if (defined($opts{cwd})) {
      $ofh->printf("cd $opts{cwd}\n");
    }
    $ofh->printf($opts{cmd}. "\n");

  } elsif ($clusterName eq "discovery") {
    $ofh->print("#!/bin/bash -l\n");
    $ofh->print("#PBS -N $opts{n}\n") if (defined($opts{n}));
    $ofh->print("#PBS -q largeq\n");
    $ofh->print("#PBS -l nodes=1:ppn=1\n");
    $ofh->print("#PBS -l walltime=$opts{time}\n");
    if (defined($opts{attr})) {
      if (ref($opts{attr}) eq "ARRAY") {
        foreach my $attr (@{$opts{attr}}) { $ofh->print("#PBS $attr\n"); }
      } else {
        $ofh->print("#PBS $opts{attr}\n");
      }
    }
    if (defined($opts{cwd})) {
      $ofh->printf("cd $opts{cwd}\n");
    } else {
      $ofh->printf("cd \$PBS_O_WORKDIR\n");
    }
    $ofh->print("\nshopt -s expand_aliases\n");
    $ofh->printf($opts{cmd}. "\n");
  }
  close($ofh);

  # submit script
  if (!defined($opts{dry})) {
    my $good = successfulSubmitOutput($clusterName);
    my $ans = $self->keepTrying($clusterName, "qsub $sname", $good);
    $ans =~ /$good/;
    if (wantarray()) {
      return ($1, $sname);
    } else {
      return $1;
    }
  } else {
    return $sname;
  }
}


=head2

 Title   :  successfulSubmitOutput
 Usage   :  my $outStr = $c->successfulSubmitOutput("anthill")
 Function:  returns a regular expression characterizing the output expected from a successful qsub command
 Returns :  string
 Args    :  1. cluster name

=cut

sub successfulSubmitOutput {
  my $cname = shift;

  if ($cname eq "anthill") {
    return 'Your job (\d+) .+ has been submitted';
  } elsif ($cname eq "discovery") {
    return '(\d+)\..+\.dartmouth\.edu';
  } else {
    GENERAL::error("unknown cluster name '$cname'");
  }
}

=head2

 Title   :  numberOfWaitingJobs
 Usage   :  $c->waitOnWaitingJobs(10)
            $c->waitOnWaitingJobs(10, "anthill")
            $c->waitOnWaitingJobs(10, "anthill", 1)
 Function:  returns when the number of waiting jobs on the cluster falls below the given limit
 Returns :  
 Args    :  1. limit number of jobs -- return when the number of jobs falls below this limit
            2. cluster name
            3. optional: set to 1 if should wait on the total number of jobs, not just those waiting

=cut

sub waitOnWaitingJobs {
  my $self = shift;
  my $ncut = shift;
  my $cname = shift;
  my $totalJobs = shift; $totalJobs = 0 if (!defined($totalJobs));
  if (!defined($cname)) {
    GENERAL::error("need to specify which cluster waiting on") if (scalar(@{$self->{clusters}}) > 1);
    $cname = $self->{clustNames}->[0];
  }

  $self->resetPause($cname, "waitOnWaitingJobs");
  while (1) {
    my $n = $totalJobs ? scalar($self->queuedJobs($cname)) : $self->numberOfWaitingJobs($cname);
    return $n if ($n < $ncut);
    $self->pause($cname, "waitOnWaitingJobs");
  }
}

=head2

 Title   :  numberOfWaitingJobs
 Usage   :  my $n = $c->numberOfWaitingJobs()
            my $n = $c->numberOfWaitingJobs("anthill")
 Function:  gets the number of jobs waiting to run on the given cluster or all clusters
 Returns :  number of jobs
 Args    :  1-n. names of clusters (no arguments to do all clusters)

=cut

sub numberOfWaitingJobs {
  my $self = shift;
  my @cnames;
  if (scalar(@_) > 0) {
    @cnames = @_;
  } else {
    @cnames = @{$self->{clustNames}};
  }

  foreach my $cname (@cnames) {
    my $ans = $self->qstat($cname, $cname);
    if ($cname eq "anthill") {
      my @ans = split("\n", GENERAL::Trim($ans));
      my $n = 0;
      foreach my $line (@ans) {
        $n++ if ($line =~ / qw /);
      }
      return $n;
    }
  }
}

=head2

 Title   :  queuedJobs
 Usage   :  my $jobIDs = $c->queuedJobs()
            my $jobIDs = $c->queuedJobs("anthill")
 Function:  returns a list of job IDs that are currently on the queued
            for the current user
 Returns :  reference to a list of job IDs
 Args    :  1. name of cluster

=cut

sub queuedJobs {
  my $self = shift;
  my $cname = shift; $cname = "anthill" if (!defined($cname));

  my $ans = $self->qstat($cname, $cname);
  if ($cname eq "anthill") {
    my @ans = split("\n", GENERAL::Trim($ans));
    @ans = @ans[2 .. scalar(@ans)-1];
    my @jobIDs;
    foreach my $line (@ans) {
      my @line = split(" ", $line);
      GENERAL::assert(GENERAL::isInteger($line[0]), "could not parse qstat line '$line'");
      push(@jobIDs, $line[0]);
    }
    return \@jobIDs;
  }
}

=head2

 Title   :  finishedJobs
 Usage   :  my ($done, $remaining) = $c->queuedJobs(\@submittedJobIDs)
            my ($done, $remaining) = $c->queuedJobs(\@submittedJobIDs, "anthill")
 Function:  given a list of job IDs that have been submitted, figures out
            which ones are still running and which ones have apparently finished
 Returns :  references to lists of job IDs that are done and those that remain
 Args    :  1. reference to a list of submitted jobs
            2. name of cluster

=cut

sub finishedJobs {
  my $self = shift;
  my $subJobs = shift;
  GENERAL::requireArgs($self, $subJobs);
  my $cname = shift; $cname = "anthill" if (!defined($cname));

  my $inQ = $self->queuedJobs();
  my %inQ;
  foreach my $id (@$inQ) { $inQ{$id} = 1; }
  my (@done, @rem);
  foreach my $id (@$subJobs) {
    if (defined($inQ{$id})) { push(@rem, $id); }
    else { push(@done, $id); }
  }
  return (\@done, \@rem);
}

=head2

 Title   :  qstat
 Usage   :  my $ans = $c->qstat("anthill")
 Function:  runs qstat and returns the result as a string (stdout and stderr are combined)
 Returns :  the resut of qstat with stderr and stdout combined
 Args    :  1. cluster name

=cut

sub qstat {
  my $self = shift;
  my $cname = shift;

  $self->resetPause($cname, "qstat");
  if ($cname eq "anthill") {
    for (my $i = 0; 1; $i++) {
      my $ans = `qstat 2>&1`;
      if ($ans =~ /error/im) {
        $self->pause($cname, "qstat");
      } else {
        return $ans;
      }
      if ($i % 10 == 0) {
        GENERAL::warning("have been failing to qstat for a while, reply:\n---------------\n$ans\n---------------\nWill keep trying and hoping for the best...\n");
      }
    }
  }
}

=head2

 Title   :  keepTrying
 Usage   :  my $ans = $c->keepTrying("anthill", "qstat")
            my $ans = $c->keepTrying("anthill", "qsub jon.sh", "Your job \d+ (.+) has been submitted")
 Function:  will keep trying to run the given command until it succeeds (pausing as necessary to go easy on the queue manager)
 Returns :  the resut of the command with stderr and stdout combined
 Args    :  1. cluster name
            2. command
            3. optional: the expected "correct" patter (not checked by default)
            4. optional: the expected erroneous pattern (default is /error/i)

=cut

sub keepTrying {
  my $self = shift;
  my $cname = shift;
  my $cmd = shift;
  my $good = shift;
  my $bad = shift;
  $bad = "error" if (!defined($bad));

  $self->resetPause($cname, "keepTrying $cmd");
  if (($cname eq "anthill") || ($cname eq "discovery")) {
    for (my $i = 0; 1; $i++) {
      my $ans = `$cmd`;
      # not seeing the good message (if it is defined) is sufficient to think there was an error
      if ((defined($good) && ($ans =~ /$good/)) || (!defined($good) && ($ans !~ /$bad/i))) {
        return $ans;
      } else {
        $self->pause($cname, "keepTrying $cmd");
      }
      if ($i % 10 == 0) {
        GENERAL::warning("have been failing to run '$cmd' for a while, reply:---------------\n$ans\n---------------\nWill keep trying and hoping for the best...\n");
      }
    }
  }
}

=head2

 Title   :  pause
 Usage   :  $c->pause("anthill", "my pause")
 Function:  sleeps for a duration of time given by the pause and increases the pause for future sleeps
 Returns :  nothing
 Args    :  1. cluster name
            2. pause name

=cut

sub pause {
  my $self = shift;
  my $cname = shift;
  my $pause = shift;
  sleep($self->{clustByName}->{$cname}->{pauses}->{$pause});
  $self->{clustByName}->{$cname}->{pauses}->{$pause} = GENERAL::min($self->{clustByName}->{$cname}->{pauses}->{$pause} + 10, $self->{clustByName}->{$cname}->{maxPause});
}

=head2

 Title   :  resetPause
 Usage   :  $c->resetPause("anthill")
 Function:  resets the given pause duration of the given cluster
 Returns :  nothing
 Args    :  1. cluster name
            2. pause name

=cut

sub resetPause {
  my $self = shift;
  my $cname = shift;
  my $pause = shift;
  $self->{clustByName}->{$cname}->{pauses}->{$pause} = $self->{clustByName}->{$cname}->{minPause};
}

1;
