package ETABLE2;
#$DEF_VER_CONTROL->{ver}->{ETABLE2} = "1.0";
$DEF_VER_CONTROL->{date}->{ETABLE2} = "06/17/04";

use strict;
use GENERAL;
use DEFINITIONS;
#use Devel::Size;

# Global variables for the module
my $INF = "(nan|inf|\\*+)"; # these strings commonly end up in energy tables if values are too large

=head2 

 Title   : writeSolvation
 Usage   : 
 Function: write out solvation term for each residues
 Returns : 
 Args    : 

=cut

sub processSolvationTable {
   my $rotamer=shift;
   my $fh = GENERAL::GetInFH(DEFINITIONS::getParmFile("reference_eef1"));
   my %solv;

   # Hash reference solvation energies by aa name
   while (<$fh>) {
     my @line=split(" ", $_);
     $solv{$line[0]}=$line[1];
   }

   # Write out reference solvation energies for each rotamer
   my $outfh = GENERAL::GetOutFH($FS_DEF->{eef_solv_tabf});
   my $sitelist=$rotamer->{dslist};
   my $counter=0;
   foreach my $s (@{$sitelist}) {
     # loop over all the residues
     foreach my $r (@{$s->{reslist}}) {
       my $resname=$r->{resname};
       # for each residue, get the energies for all the  rotamers
       foreach (@{$r->{rotamerlist}}) {
         # print the energy value into the file
         #print $solv->{$resname}, "\n";
         $outfh->printf("%14.5f\n", $solv{$resname});
       }
     }
   }
   close($outfh);
}



=head2

 Title   : writeDesignEnergyTable
 Usage   :
 Function: Creates binary energy table(s) for design. It simply
           checks the energy tracker structure to see which terms
           are to be included with the current energy scheme.
 Returns : nothing
 Args    : 1. rotamer object
           2. energy configuration parameter hash table
           3. This is an optional parameter designed to allow the user
              to save the energies corresponding to a smaller design
              problem (with only a subset of design sites/residues/rotamers).
              This parameter is a string tag, which is going to be used
              to check whether every rotamer belongs to the smaller subset
              or not (0 means it does not, and 1 means it does). So if the
              parameter is "subset_valid" lets say, then $rot->{subset_valid}
              will be used to decide whether this rotamer (or rather the
              energies involving this rotamer) are to be included in the tables
              or not.
           4. Optional: self table. If not specified, it will be read from
              disk according to the entries in the energy tracker.
           5. Optional: unfold table. If not specified, it will be read from
              disk according to the entries in the energy tracker. Note, that if
              it is specified, it will be assumed that the Boltzmann averaging (if
              it is neccessary) is already performed.
           6. Optional: pair table. If not specified, it will be read from
              disk according to the entries in the energy tracker.
           7. Optional: self valid flag. If set, it means that the self energies
              passed to this function contain only the energies for valid rotamers.
              Otherwise, it is assumed that all self energies are passed. The latter
              is the default.

=cut

sub writeDesignEnergyTable {
  my $rot = shift;
  my $p = shift;
  GENERAL::requireArgs($rot, $p);
  my $svalid = shift;
  if (!defined($svalid)) { $svalid = "valid"; }
  my $st = shift;
  my $ut = shift;
  my $pt = shift;
  my ($du, $ds, $dp) = (0, 0, 0);
  if (defined($ut)) { $du = 1; }
  if (defined($st)) { $ds = 1; }
  if (defined($pt)) { $dp = 1; }
  my $self_valid = shift;
  if (!defined($self_valid)) { $self_valid = 0; }

  # ----- First get the total self energy
  my $ns;
  if ($self_valid) { $ns = $rot->countValidSelfTotal(); }
  else { $ns = $rot->countSelfTotal(); }
  if (!$ds) {
    print "Reading self term... " . GENERAL::GetTime() . "\n";
    $st = getSelfTotal($p, 2);
    print "Done reading self term. " . GENERAL::GetTime() . "\n";
  }
  # Sanity check
  if (scalar(@$st) != $ns) {
    GENERAL::error("Error in writeDesignEnergyTable: total number of self terms is " .
    scalar(@$st) . ", while the number of rotamers is $ns!");
  }

  # ----- Now get the total unfolded energy
  if ($p->{unfolded_state_calc} =~ /yes/i) {
    if (!$du) {
      print "Reading unfold term... " . GENERAL::GetTime() . "\n";
      $ut = getUnfoldTotal($p, 1);
      print "Done reading unfold term. " . GENERAL::GetTime() . "\n";
    }
    # Sanity check
    if (scalar(@$ut) != $ns) {
      GENERAL::error("Error in writeDesignEnergyTable: total number of unfold terms is " .
      scalar(@$ut) . ", while the number of rotamers is $ns!");
    }
    # The unfolded state is special, in that we may want to use the
    # ensemble average of all rotamers for each residue.
    if ((!$du) && ($p->{escheme} eq $p->{escheme})) {
      print "Calculating ea for unfold term... " . GENERAL::GetTime() . "\n";
      ETABLE2::stateEnsembleAverage($rot, $ut);
      print "Done calculating ea for unfold term... " . GENERAL::GetTime() . "\n";
    }
  }

  # ----- And finally, the the total pair energy
  my $np = $rot->countValidPairTotal();
  my $nel; # number of actual elements (either read of passed)
  if (!$dp) {
    print "Reading pair term... " . GENERAL::GetTime() . "\n";
    ($pt, $nel) = getPairTotal($p, $np, 1);
    print "Done reading pair term. " . GENERAL::GetTime() . "\n";
  } else {
    $nel = scalar(@$pt);
  }
  # Sanity check
  if ($nel != $np) {
    GENERAL::error("Error in writeDesignEnergyTable: total number of pair terms is $nel, while the total number of pairs of valid rotamers is $np!");
  }

  print "Writing... " . GENERAL::GetTime() . "\n";
  my $etf;
  # Do we want the unfolded energy separately or added to self?
  if ($p->{unfolded_state_sep}) {
    $etf = GENERAL::GetOutFH($FS_DEF->{ener_nounf_binf});
    if ($p->{unfolded_state_calc} =~ /yes/i) {
      my $utf = GENERAL::GetOutFH($FS_DEF->{unf_binf});
      my $i = 0;
      foreach my $site (@{$rot->{dslist}}) {
        foreach my $res (@{$site->{reslist}}) {
          foreach my $r (@{$res->{rotamerlist}}) {
            if ($r->{valid} && $r->{$svalid}) { $utf->print(pack("d",  $ut->[$i])); }
            next if ($self_valid && ($r->{valid} == 0));
            $i++;
          }
        }
      }
      # We need to add a list of zeros for pairwise terms in in the unfolded energy table
      # loop over all pairs of valid rotamers
      my $sitelist = $rot->{dslist};
      for (my $s1=0; $s1 <=$#{$sitelist}; $s1++) {
        for (my $res1=0; $res1<=$#{$sitelist->[$s1]->{reslist}}; $res1++) {
          for(my $rot1=0; $rot1<=$#{$sitelist->[$s1]->{reslist}->[$res1]->{rotamerlist}}; $rot1++) {
            my $r1 = $sitelist->[$s1]->{reslist}->[$res1]->{rotamerlist}->[$rot1];
            next if ($r1->{valid} == 0);
            for (my $s2=$s1+1; $s2 <=$#{$sitelist}; $s2++) {
              for (my $res2=0; $res2<=$#{$sitelist->[$s2]->{reslist}}; $res2++) {
                for(my $rot2=0; $rot2<=$#{$sitelist->[$s2]->{reslist}->[$res2]->{rotamerlist}}; $rot2++) {
                  my $r2 = $sitelist->[$s2]->{reslist}->[$res2]->{rotamerlist}->[$rot2];
                  next if ($r2->{valid} == 0);
                  if ($r1->{$svalid} && $r2->{$svalid}) {
                    $utf->print(pack("d", 0.0));
                  }
                }
              }
            }
          }
        }
      }
      # Write a log for the user of what went into the unfold energy file
      LogEnergyBin($FS_DEF->{unf_ctrlf}, $FS_DEF->{unf_logf}, $FS_DEF->{unf_binf}, "Unfold State Energy:\n");
      close($utf);
    }

    my $i = 0;
    foreach my $site (@{$rot->{dslist}}) {
      foreach my $res (@{$site->{reslist}}) {
        foreach my $r (@{$res->{rotamerlist}}) {
          if ($r->{valid} && $r->{$svalid}) { $etf->print(pack("d",  $st->[$i])); }
          next if ($self_valid && ($r->{valid} == 0));
          $i++;
        }
      }
    }
    # Write a log for the user of what went into the total energy file
    LogEnergyBin($FS_DEF->{self_ctrlf}, $FS_DEF->{ener_nounf_logf}, $FS_DEF->{ener_nounf_binf}, "Self Energy:\n");
  } else {
    $etf = GENERAL::GetOutFH($FS_DEF->{ener_tot_binf});
    my $i = 0;
    foreach my $site (@{$rot->{dslist}}) {
      foreach my $res (@{$site->{reslist}}) {
        foreach my $r (@{$res->{rotamerlist}}) {
          if ($r->{valid} && $r->{$svalid}) { $etf->print(pack("d",  $st->[$i] - (($p->{unfolded_state_calc} =~ /yes/i) ? $ut->[$i] : 0)));}
          next if ($self_valid && ($r->{valid} == 0));
          $i++;
        }
      }
    }
    # Write a log for the user of what went into the total energy file
    LogEnergyBin($FS_DEF->{unf_ctrlf}, $FS_DEF->{ener_tot_logf}, $FS_DEF->{ener_tot_binf}, "Unfold State Energy:\n") if ($p->{unfolded_state_calc} =~ /yes/i);
    LogEnergyBin($FS_DEF->{self_ctrlf}, $FS_DEF->{ener_tot_logf}, $FS_DEF->{ener_tot_binf}, "Self Energy:\n", 1);
  }

  # And now, complete the energy table by including the pairwise terms
  my $sitelist = $rot->{dslist};
  my $i = 0;
  for (my $s1=0; $s1 <=$#{$sitelist}; $s1++) {
    for (my $res1=0; $res1<=$#{$sitelist->[$s1]->{reslist}}; $res1++) {
      for(my $rot1=0; $rot1<=$#{$sitelist->[$s1]->{reslist}->[$res1]->{rotamerlist}}; $rot1++) {
        my $r1 = $sitelist->[$s1]->{reslist}->[$res1]->{rotamerlist}->[$rot1];
        next if ($r1->{valid} == 0);
        for (my $s2=$s1+1; $s2 <=$#{$sitelist}; $s2++) {
          for (my $res2=0; $res2<=$#{$sitelist->[$s2]->{reslist}}; $res2++) {
            for(my $rot2=0; $rot2<=$#{$sitelist->[$s2]->{reslist}->[$res2]->{rotamerlist}}; $rot2++) {
              my $r2 = $sitelist->[$s2]->{reslist}->[$res2]->{rotamerlist}->[$rot2];
              next if ($r2->{valid} == 0);
              if ($r1->{$svalid} && $r2->{$svalid}) {
                $etf->print(pack("d", $pt->[$i]));
              }
              $i++;
            }
          }
        }
      }
    }
  }
  close($etf);
  print "Done writing... " . GENERAL::GetTime() . "\n";
  # Continue log for the user of what went into the total energy file
  if ($p->{unfolded_state_sep}) {
    LogEnergyBin($FS_DEF->{pair_ctrlf}, $FS_DEF->{ener_nounf_logf}, $FS_DEF->{ener_tot_binf}, "Pairwise Energy:\n", 1);
  } else {
    LogEnergyBin($FS_DEF->{pair_ctrlf}, $FS_DEF->{ener_tot_logf}, $FS_DEF->{ener_tot_binf}, "Pairwise Energy:\n", 1);
  }

  # Return the energy tables
  return ($st, $ut, $pt);
}


=head2

 Title   : stateEnsembleAverage
 Usage   : ETABLE2::stateEnsembleAverage($rotamer, \@ener);
 Function: Ensemble averages rotamer energies. I.e. for each rotamer, sets its energy
           to the ensemble average energy of all rotamers of that residue.
 Returns : Nothing.
 Args    : 1. Rotamer object
           2. Reference to energy array. There must be as many entries in the energy
              array as there are rotamers in the rotamer object.

=cut

sub stateEnsembleAverage {
  my $rot = shift;
  my $ener = shift;
  GENERAL::requireArgs($rot, $ener);

  if ($rot->countSelfTotal() != scalar(@$ener)) {
    GENERAL::error("Total number of rotamers is " . ($rot->countSelfTotal()) . " while total number of terms is " . scalar(@$ener));
  }
    
  my $i=0;
  foreach my $site (@{$rot->{dslist}}) {
    foreach my $res (@{$site->{reslist}}) {
      my $rote = [];
      my $j = $i;
      foreach my $r (@{$res->{rotamerlist}}) {
        push(@$rote, $ener->[$j]);
        $j++;
      }
      my $ensemble = GENERAL::Ensemble($rote);
      foreach my $r (@{$res->{rotamerlist}}) {
        $ener->[$i] = $ensemble->{average};
        $i++;
      }
    }
  }
  return $ener;
}

=head2

 Title   : LogEnergyBin
 Usage   :
 Function: Creates a log of the content of the given energy bin file.
 Returns : nothing
 Args    : 1. control file name
           2. log file name
           3. energy bin file name
           4. optional: overwrite or append

=cut

sub LogEnergyBin {
  my $ctrlf = shift;
  my $logf = shift;
  my $binf = shift;
  my $title = shift;
  my $app = shift;
  if (!defined $app) { $app = 0; }

  my $pc = GENERAL::LoadStruct($ctrlf, 1);
  my $log;
  if ($app) {
    $log = GENERAL::GetAppendOutFH($logf);
  } else {
    $log = GENERAL::GetOutFH($logf);
  }

  if (!$app) {
    print $log "File $binf was created with the following terms and scale factors:\n\n";
  }

  if (defined($title)) {
    print $log $title;
  }

  for (my $i=0; $i < scalar(@{$pc->{current}->{col_names}}); $i++) {
    print $log $pc->{current}->{col_names}->[$i] . "\t:\t" . $pc->{current}->{sf}->[$i] ."\n";
  }
  print $log "\n";
  close($log);
}

=head2

 Title   : processSelfTable
 Usage   :
 Function: Creates a large informative table for the user of
           all the calculated self energy terms. Does this by
           reading the energy progress file and simply combining
           all the terms relating to self energy.
 Returns : nothing
 Args    : 1. rotamer object
           2. delimiter string to use when creating the table
           3. optional - flag to specify whether the table is to be
              written to disk or simply returned. Default is 1 (write).

=cut

sub processSelfTable {
  my $rot = shift;
  my $delim = shift;
  if (!defined $delim) {$delim = ",";}
  my $write = shift;
  if (!defined($write)) { $write = 1; }

  # Get the term tracker structure and fish out the right entry
  my $th = GENERAL::LoadStruct($FS_DEF->{term_trackerf}, 1);
  my $en = $th->{self};

  # Now actually read in all tables and combine them
  my @arr;
  my $table = \@arr;
  my $rows = -1;
  my $pterm = -1;

  # Open large table for output
  my $outf;
  if ($write) { $outf = GENERAL::GetOutFH($FS_DEF->{self_tabf}); }

  # Go through each term
  my @names;
  foreach my $tkey (sort(keys(%{$en}))) {
    my $term = $en->{$tkey};
    my $ttab = ();

    # Have we calculated this term?
    if (!defined($term->{files})) { next; }

    # The header for this term
    foreach my $col (@{$term->{col_names}}) {
      push(@names, $col);
    }

    foreach my $file (@{$term->{files}}) {
      my $tab = readEnergyTable($file);
      # Append to the end since the same term
      push(@{$ttab}, @$tab);
    }
    if ($rows < 0) { $rows = scalar(@$ttab); }
    elsif ($rows != scalar(@$ttab)) {
      die "Error in processSelfTable: term \"$term->{name}\" has an inconsistent number of rows (" .
      scalar(@$ttab) . ") as compared to term \"$pterm->{name}\" ($rows)!\n";
    }
    # Merge rows since different terms
    for (my $j=0; $j < scalar(@$ttab); $j++) {
      push(@{$table->[$j]}, @{$ttab->[$j]});
    }
    $pterm = $term;
  }

  # Make sure we have the same number of rows as rotamers
  if (scalar(@$table) != $rot->countSelfTotal()) {
    die "Error in processSelfTable: total number of terms is " .
    scalar(@$table) . " while the number of rotamers is " .
    $rot->countSelfTotal() . "\n";
  }

  # And finally output the combined table
  if ($write) {
    print $outf "\"Site\"$delim\"Residue Name\"$delim\"Rotamer Number\"$delim\"Eliminated?\"$delim\"" . join("\"$delim\"", @names) . "\"\n";
    # Then the data
    foreach my $site (@{$rot->{dslist}}) {
      foreach my $res (@{$site->{reslist}}) {
        my $i = 1;
        foreach my $rot (@{$res->{rotamerlist}}) {
          print $outf "$site->{chain}$site->{iresnum}$delim$res->{resname}$delim$i$delim$rot->{valid}";
          $i++;
          my $row = shift(@{$table});
          foreach my $val (@{$row}) {
            print $outf "$delim$val";
          }
          print $outf "\n";
        }
      }
    }
    close($outf);
  } else {
    return (\@names, $table);
  }
  
}


=head2 

 Title   :  processUnfoldTable
 Usage   : 
 Function: Creates a large informative table for the user of
           all the calculated unfold energy terms. Does this by
           reading the energy progress file and simply combining
           all the terms relating to unfolded energy.
 Returns : nothing
 Args    : 1. rotamer object
           2. optional - delimiter string to use when creating the table
           3. optional - flag to specify whether the table is to be
              written to disk or simply returned. Default is 1 (write).

=cut

sub processUnfoldTable {
  my $rot = shift;
  my $delim = shift;
  if (!defined $delim) {$delim = ",";}
  my $write = shift;
  if (!defined($write)) { $write = 1; }

  # Get the term tracker structure and fish out the right entry
  my $th = GENERAL::LoadStruct($FS_DEF->{term_trackerf}, 1);
  my $en = $th->{unfold};

  # Now actually read in all tables and combine them
  my $table = ();
  my $rows = -1;
  my $pterm = -1;

  # Open the large energy table file
  my $outf;
  if ($write) { $outf = GENERAL::GetOutFH($FS_DEF->{unfold_tabf}); }

  # Go through each term
  my @names;
  foreach my $tkey (sort(keys(%{$en}))) {
    my $term = $en->{$tkey};
    my $ttab = ();

    # Have we calculated this term?
    if (!defined($term->{files})) { next; }

    # The header for this term
    foreach my $col (@{$term->{col_names}}) {
      push(@names, $col);
    }

    foreach my $file (@{$term->{files}}) {
      my $tab = readEnergyTable($file);
      # Append to the end since the same term
      push(@{$ttab}, @$tab);
    }
    if ($rows < 0) { $rows = scalar(@$ttab); }
    elsif ($rows != scalar(@$ttab)) {
      die "Error in processUnfoldTable: term \"$term->{name}\" has an inconsistent number of rows (" .
      scalar(@$ttab) . ") as compared to term \"$pterm->{name}\" ($rows)!\n";
    }
    # Merge rows since different terms
    for (my $j=0; $j < scalar(@$ttab); $j++) {
      push(@{$table->[$j]}, @{$ttab->[$j]});
    }
    $pterm = $term;
  }

  # Make sure we have the same number of rows as rotamers
  if (scalar(@$table) != $rot->countSelfTotal()) {
    die "Error in processUnfoldTable: total number of terms is " .
    scalar(@$table) . " while the number of rotamers is " .
    $rot->countSelfTotal() . "\n";
  }

  # And finally output the combined table
  if ($write) {
    print $outf "\"Site\"$delim\"Residue Name\"$delim\"Rotamer Number\"$delim\"Eliminated?\"$delim\"" . join("\"$delim\"", @names) . "\"\n";
    # Then the data
    foreach my $site (@{$rot->{dslist}}) {
      foreach my $res (@{$site->{reslist}}) {
        my $i = 1;
        foreach my $rot (@{$res->{rotamerlist}}) {
          print $outf "$site->{chain}$site->{iresnum}$delim$res->{resname}$delim$i$delim$rot->{valid}";
          $i++;
          my $row = shift(@{$table});
          foreach my $val (@{$row}) {
            print $outf "$delim$val";
          }
          print $outf "\n";
        }
      }
    }
    close($outf);
  } else {
    return (\@names, $table);
  }
}




=head2 

 Title   : processPairTable
 Usage   : 
 Function: Process the energy files in pair directory and 
           report total energy table
 Returns : 
 Args    : 1. rotamer object
           2. delimiter string to use when creating the table
           3. optional - flag to specify whether the table is to be
              written to disk or simply returned. Default is 1 (write).

=cut

sub processPairTable {
  my $rot = shift;
  my $delim = shift;
  if (!defined $delim) {$delim = ",";}
  my $write = shift;
  if (!defined($write)) { $write = 1; }

  # Get the term tracker structure and fish out the right entry
  my $th = GENERAL::LoadStruct($FS_DEF->{term_trackerf}, 1);
  my $en = $th->{pair};

  # Now actually read in all tables and combine them
  my @table;
  my $rows = -1;
#  my $pterm = -1;

  # Open the large energy table file
  my $outf;
  if ($write) { $outf = GENERAL::GetOutFH($FS_DEF->{pair_tabf}); }
  
  # Go through each term
  my @names;
  foreach my $tkey (sort(keys(%{$en}))) {
    my $term = $en->{$tkey};
    my $ttab = ();

    # Have we calculated this term?
    if (!defined($term->{files})) { next; }

    # The header for this term
    foreach my $col (@{$term->{col_names}}) {
      push(@names, $col);
    }
    ETABLE2::readEnergyTerm($term, undef, \@table);
  }

  # Make sure we have the same number of rows as pairs of rotamers
  if (scalar(@table) != $rot->countValidPairTotal()) {
    die "Error in processPairTable: total number of terms is " .
    scalar(@table) . " while the number of valid rotamers is " .
    $rot->countValidPairTotal() . "\n";
  }

  # And finally output the combined table
  if ($write) {
    print $outf "\"Site 1\"$delim\"Residue Name 1\"$delim\"Rotamer Number 1\"$delim\"Site 2\"$delim\"Residue Name 2\"$delim\"Rotamer Number 2\"$delim\"".
                join("\"$delim\"", @names) . "\"\n";
    # Then the data
    my $sitelist = $rot->{dslist};
    for (my $s1=0; $s1 <=$#{$sitelist}; $s1++) {
      my $site1 = $sitelist->[$s1];
      for (my $res1=0; $res1<=$#{$sitelist->[$s1]->{reslist}}; $res1++) {
        my $r1 = $sitelist->[$s1]->{reslist}->[$res1];
        my $i = 0;
        for(my $rot1=0; $rot1<=$#{$sitelist->[$s1]->{reslist}->[$res1]->{rotamerlist}}; $rot1++) {
          $i++;
          next if ($sitelist->[$s1]->{reslist}->[$res1]->{rotamerlist}->[$rot1]->{valid} == 0);
          for (my $s2=$s1+1; $s2 <=$#{$sitelist}; $s2++) {
            my $site2 = $sitelist->[$s2];
            for (my $res2=0; $res2<=$#{$sitelist->[$s2]->{reslist}}; $res2++) {
              my $r2 = $sitelist->[$s2]->{reslist}->[$res2];
              my $j = 0;
              for(my $rot2=0; $rot2<=$#{$sitelist->[$s2]->{reslist}->[$res2]->{rotamerlist}}; $rot2++) {
                $j++;
                next if ($sitelist->[$s2]->{reslist}->[$res2]->{rotamerlist}->[$rot2]->{valid} == 0);
                print $outf "$site1->{chain}$site1->{iresnum}$delim$r1->{resname}$delim$i$delim";
                print $outf "$site2->{chain}$site2->{iresnum}$delim$r2->{resname}$delim$j";
                my $row = shift(@table);
                foreach my $val (@{$row}) {
                  print $outf "$delim$val";
                }
                print $outf "\n";
              }
            }
          }
        }
      }
    }
    close($outf);
  } else {
    return (\@names, \@table);
  }
}



=head2 

 Titlm   : 
 Usage   : 
 Function:  process the energy set any term above 
            the ecut value ( criterion for clash ) 
            as a million 
 Returns : 
 Args    : 

=cut


sub getSelfCrash {
    my $arr=shift;
    my $selfecut=shift;
    my $bige=1000000.0;
    foreach my $j (@{$arr}) {
      foreach my $k (@{$j}) {
        #print $k, "\n";
        $k=$bige if ( $k =~ /\*/ );
        $k=$bige if ( $k =~ /inf/ );
        $k=$bige if ( $k =~ /nan/ );
        $k=$bige if ($k >= $selfecut);
      }
    }
    return $arr;
} 



=head2 

 Title   : 
 Usage   : 
 Function:  process the energy of pair term. 
            set the energy * as 1 million 
            * indicate an infinite energy or clash
 Returns : 
 Args    : 

=cut



sub getPairCrash {

    my $arr=shift;
#    my $selfecut=shift;
    my $bige=1000000.0;
    foreach my $j (@{$arr}) {
        my $col=0;
        foreach my $k (@{$j}) {
           #print $k, "\n";
           $k=$bige if ( $k =~ /\*/ );
           $k=$bige if ( $k =~ /inf/ );
           $k=$bige if ( $k =~ /nan/ );
           $col++;
           if ($col == 2 && $k >=100.0) {
		$k=100.0;
           }
           if ($col == 2 && $k <=-100.0) {
                $k=-100.0;
           }
	   # i decide to not put the energy above ecut to ceiling 
	   # sounds like it gives the problems in scmf calculations
        }
    }
   
    return $arr;
}


sub getPairCrash2 {
     

     my $arr=shift;
#    my $selfecut=shift;
   
     my $bige=500.0;

    foreach my $j (@{$arr}) {
           #print $k, "\n";
           $j=$bige if ($j =~ /inf/ );
           $j=$bige if ( $j =~ /nan/ );
           $j=$bige if ( $j > $bige );
           # i decide to not put the energy above ecut to ceiling
           # sounds like it gives the problems in scmf calculations
    }
  
    return $arr;
}



sub getSelfTotalEnergy {

     my $suffix=shift; 
     $suffix='vdw' if (!defined $suffix);
     GENERAL::cchdir("EnergyTable");

     my $total='selftotal-'."$suffix".'.dat';
     my $energy=&GENERAL::GetInFH($total);

     my @arr=<$energy>;
     foreach my $j (@arr) {
        chomp $j;
        $j=&GENERAL::Trim($j);
     }

     GENERAL::cchdir("../");
     return \@arr;
}


=head2 

 Title   : pickEnergyTerms
 Usage   : 
 Function: given a total energy table array and enery  
           terms we want to combine. produce a new 
           total energy array.  
 Returns : 
 Args    : 

=cut


sub pickEnergyTerms {

    my $totalenergy=shift;
    my $eterms=shift;

    my @newarray;
    for(my $i=0; $i<=$#{$totalenergy}; $i++) {
      my $line=$totalenergy->[$i];
      my $totale=0.0;
      # loop over all the energy terms we want to add up
      foreach (@{$eterms}) {
        if ($line->[$_] =~ /(nan|inf|\*)/) {
          $totale = "nan";
          last;
        }
        $totale += $line->[$_];
      }
      push (@newarray, $totale);
    }

    return \@newarray;

}




=head2

 Title   :
 Usage   :
 Function: a function to decide which columns to combine in energy table
           to produce new energy array
 Returns :
 Args    :

=cut



sub getEterms {

    my $option=shift;
    my $pick=shift;
    my @eterms;
    if (!defined $option) { die "Error in getEterms: parameter \"option\" undefined!\n"; }
    if (!defined $pick) { die "Error in getEterms: parameter \"pick\" undefined!\n"; }


    if ( $pick==1 ) {
        if ($option eq 'unfold' || $option eq 'self') {
            @eterms=(1,3,5,7);
            return \@eterms;
         }
        if ($option eq 'pair') {
            @eterms=(0,2,4);
            return \@eterms;
         }
    }

    if ( $pick==2) {
        if ($option eq 'unfold' || $option eq 'self') {
            @eterms=(1,2,3,5);
            return \@eterms;
         }
        if ($option eq 'pair') {
            @eterms=(0,1,2);
            return \@eterms;
         }
    }

    if ( $pick==3) {
        if ($option eq 'unfold' || $option eq 'self') {
            @eterms=(1,3,5);
            return \@eterms;
         }
        if ($option eq 'pair') {
            @eterms=(0,2);
            return \@eterms;
         }
    }

    if ( $pick==4) {
        if ($option eq 'unfold' || $option eq 'self') {
            @eterms=(1,2,3);
            return \@eterms;
         }
        if ($option eq 'pair') {
            @eterms=(0,1);
            return \@eterms;
         }
    }


    if ( $pick==5) {
        if ($option eq 'unfold' || $option eq 'self') {
            @eterms=(1,3,7);
            return \@eterms;
         }
        if ($option eq 'pair') {
            @eterms=(0,4);
            return \@eterms;
         }
    }

   if ( $pick==6) {
         if ($option eq 'unfold' || $option eq 'self') {
            @eterms=(1,3);
            return \@eterms;
         }
        if ($option eq 'pair') {
            @eterms=(0);
            return \@eterms;
         }
    }

}



=head2 

 Title   : readEnergyTable
 Usage   : my $table = ETABLE2::readEnergyTable($file_name)
 Function: Reads energy table file, returns a reference
           to a list of lists. If asked, checks to make sure that
           all values are numeric or fixes some common problems (i.e. inf,
           nan, ****** all get substituted by 10^6) or both.
 Returns : An array reference.
 Args    : 1. File name of the energy table.
           2. Optional strictness level:
              0 - do nothing (just read and split)
              1 - fix common problems
              2 - throw a warning if values are not numeric
              3 - 1+2
              4 - die with a message if values are not numeric
              5 - 4+1 - default
=cut

sub readEnergyTable {

    my $fname = shift;
    my $strict = shift;
    if (!defined $strict) {$strict = 5;}

    my $etable = GENERAL::GetInFH($fname, 1);
    my $earray=[];
    my $row;
    while ($row = ETABLE2::readEnergyTableRow($etable, ($strict =~ /^[135]$/))) {
      for (my $i=0; ($i < scalar(@$row)) && ($strict =~ /^[2345]$/); $i++) {
        if (!GENERAL::isReNumeric($row->[$i])) {
          my $msg = "Value $row->[$i] is not numeric while reading $fname (part of line \"" . join(" ", @$row) . "\")";
          if ($strict =~ /^[23]$/) {
            GENERAL::error($msg);
          } elsif ($strict =~ /^[45]$/) {
            GENERAL::warning($msg);
          }
        } else {
          $row->[$i] = $row->[$i] + 0.0; # convert to number -> save tons of space
        }
      }
      push (@{$earray}, $row);
    }
    close $etable;
    return $earray;
}



=head2

 Title   : readEnergyTableRow
 Usage   : my $row = ETABLE2::readEnergyTableRow($file_handle)
 Function: Reads the next row in the energy table given the file
           handle. Returns a reference to a list of values. Also,
           if asked, fixes some common problems (i.e. inf,
           nan, ****** all get substituted by 10^6).
           Returns 0 on end of file.
 Returns : An array reference or 0 on EOF.
 Args    : 1. File handle of the energy table.
           2. Optional strictness flag: if set, common problems are
           fixed (strings indicating infinity are set to 10^6).
           The default for the flag is to be set.
=cut

sub readEnergyTableRow {
  my $fh = shift;
  my $strict = shift;
  if (!defined($strict)) { $strict = 1; }

  while (<$fh>) {
    next if (/^#/);
    next if ( /^\s+$/);
    chomp;
    if ($strict) { s/(\s$INF\s|^$INF\s|\s$INF$|^$INF$)/ 1000000 /g; }
    my @line = split(" ", $_);
    return \@line;
  }
  return 0;
}

=head2 

 Title   :  
 Usage   : 
 Function:  write unfold energy term: 
           (1) Multe VDW + EEF + torsion + EEF-electrostatics: 1,3,5,7
           (2) Multe VDW + EEF + torsion + Mutle-elect: 1,2,3,5
           (3) Multe VDW + EEF solvation + torsion:  1,3,5
           (4) Multe only:   1,2,3
           (5) Multe VDW + EEF-elec + Torsion: 1,3,7  
           (6) Multe VDW only
 Returns : 
 Args    : 

=cut

sub writeUnfoldData {

    my $rotamer=shift;
    my $suffix=shift;
    $suffix='default' if (!defined $suffix);
    my $choice=shift;
    $choice=1 if (!defined $choice);
    my $eterms=&ETABLE2::getEterms('unfold', $choice);
    my $bige=1000000.0;
    my $UnfoldEnergyTable=shift;
    $UnfoldEnergyTable='UnfoldEnergyTable' if (!defined $UnfoldEnergyTable);
    my $diel=shift;
    $diel="" if (!defined $diel || $diel ==4);

    # now determine if we have solvation term
    my $sol=map {/5/} @{$eterms};

    # now loop over all the design site
    my $sitelist=$rotamer->{dslist};

    &GENERAL::MakeDir("EnergyTable");
    &GENERAL::MakeDir("UPS");
    chdir("EnergyTable");

    # readin energy Table
    my $earray=&ETABLE2::readEnergyTable($UnfoldEnergyTable);
    my $energy=&ETABLE2::pickEnergyTerms($earray, $eterms);

    system("cp $PARAM_DEF/reference_eef1 .");
    my $fh=&GENERAL::GetInFH('reference_eef1');
    my $solv={};

    while (<$fh>) {
         my @line=split(" ", $_);
         $solv->{$line[0]}=$line[1];
    }
    close $fh;

    my $dat='unfold'.$diel.'-'."$suffix".'.dat';
    my $datf=&GENERAL::GetOutFH($dat);
    my $unfoldbin='unfold'.$diel.'-'."$suffix".'.bin';
    my $unfoldbinf=&GENERAL::GetOutFH($unfoldbin);

    my $unfoldsum='../UPS/summary'.$diel.'-'."$suffix".'.res';
    my $unfolddist='../UPS/summary'.$diel.'-'."$suffix".'.dist';
    my $out1=&GENERAL::GetOutFH($unfoldsum);
    my $out2=&GENERAL::GetOutFH($unfolddist);

 
    my $counter=0;
    foreach my $s (@{$sitelist}) {
	
	my $chain=lc($s->{chain});
        my $cchain=$s->{chain};
	my $iresnum=$s->{iresnum};


	# loop over all the residues

	foreach my $r (@{$s->{reslist}}) {
  
	    # for each residue create a new array to hold the energy
	    my @temp;
	    
	    my $resname=$r->{resname};

	    # for each residue, get the energies for all the  rotamers 
	    for (my $rot=0; $rot<=$#{$r->{rotamerlist}}; $rot++) {
		#push the energy below the crash value into @temp
		my $tempe=shift @{$energy};
		push (@temp,$tempe);
	    }

	    # now do ensemble average and print it into the file 
	    my $ensemble=&GENERAL::Ensemble(\@temp);
            if ($sol) {
               print $resname, "\n"  if (!defined $solv->{$resname});
               $ensemble->{'average'} +=$solv->{$resname};
            }
	    $out1->printf("%1s%5d%5s%12.5f\n",$cchain, $iresnum, $resname, $ensemble->{'average'});

	    my $index=-1;
	    foreach (@{$r->{rotamerlist}}) {
		# we only want to print out the number of rotamers after backbone 
		# clashing elimination. the number of 0.000 should be the same as the number of 
                # pair energy terms
		if ($_->{valid} ==1 ) {
		    # print the energy value into the file
		    $counter++;
		    $datf->printf("%8d%12.5f\n", $counter, $ensemble->{'average'});
                    my $average= $ensemble->{'average'}; 
                    $unfoldbinf->print(pack("d",  $average));
                }
                    $index++;
		    my $state;
		    # check g- or trans or g+
	            if ( $resname eq "ALA" ) {
			$state='g+';	
		    } elsif ($_->{chi}->[0] < 120.0 &&  $_->{chi}->[0] >= 0.0 ) {
			$state='g+';
		    } elsif (  $_->{chi}->[0] >= -120.0 && $_->{chi}->[0] < 0.0 ) {
			$state='g-';
		    } else {
			$state='t';
		    }
		    # now print out: chain, resnum, resname, state, probablitiy
		    $out2->printf("%1s%5d%5s%3s%6.2f\n", $cchain, $iresnum, $resname, $state, $ensemble->{'probability'}->[$index]);
#                    print $cchain, $iresnum,"  ", $resname, "\n";
		
	    }
	}
    }

    # finally, we need to add a list of zero for pairwise term in unfoldenergy.dat file
    
    # loop over all the pairwise term
    for (my $s1=0; $s1 <=$#{$sitelist}; $s1++) {
	for (my $res1=0; $res1<=$#{$sitelist->[$s1]->{reslist}}; $res1++) {
	    for (my $s2=$s1+1; $s2 <=$#{$sitelist}; $s2++) {
		for (my $res2=0; $res2<=$#{$sitelist->[$s2]->{reslist}}; $res2++) {
		    for(my $rot1=0; $rot1<=$#{$sitelist->[$s1]->{reslist}->[$res1]->{rotamerlist}}; $rot1++) {
                        next if ($sitelist->[$s1]->{reslist}->[$res1]->{rotamerlist}->[$rot1]->{valid} ==0 );
                        for(my $rot2=0; $rot2<=$#{$sitelist->[$s2]->{reslist}->[$res2]->{rotamerlist}}; $rot2++) {
                            next if ($sitelist->[$s2]->{reslist}->[$res2]->{rotamerlist}->[$rot2]->{valid} ==0 );
			    # every residue pair has a pair.inp
			    $counter++;
			    $datf->printf("%8d%14.5f\n", $counter, 0.0);
                            $unfoldbinf->print(pack("d", 0.0));
			}
		    }
		}
	    }
	}
    }

   close $unfoldbinf;
   close $datf;
   close $out2;
   close $out1;
   chdir("../");
}



=head2
 Title   :  writeUnfoldData2
 Usage   :  writeUnfoldData2($rotamer, 'UnfoldEnergyTable', $vdw, $diel, $eef);

=cut


sub writeUnfoldData2 {

  my $rotamer=shift;
  my $parm=shift;

  my $vdw = $parm->{vdw};
  my $eef = $parm->{eef};
  my $diel = $parm->{dde};

  if (!defined $vdw || !defined $diel || !defined $eef) {
    die "Error in writeUnfoldData2: energy components not defined!\n";
  }
  my $multe = 4.0/$diel;
  my $suffix = sprintf("%.1f_%s_%.1f", $vdw, $diel, $eef);
  # now loop over all the design site
  my $sitelist = $rotamer->{dslist};

  # readin energy Table
  my $earray = ETABLE2::readEnergyTable($FS_DEF->{unfold_tabf});
  my $energy = ETABLE2::pickEnergyTerms2($earray, $vdw, $multe, $eef);
  my $fh = GENERAL::GetInFH(DEFINITIONS::getParmFile("reference_eef1"));
  my $solv={};

  while (<$fh>) {
    my @line=split(" ", $_);
    $solv->{$line[0]}=$line[1];
  }
  close($fh);

  my $dat = $FS_DEF->{unfold_datf};
  $dat =~ s/%/$suffix/;
  my $datf = GENERAL::GetOutFH($dat);
  my $unfoldbin = $FS_DEF->{unfold_binf};
  $unfoldbin =~ s/%/$suffix/;
  my $unfoldbinf = GENERAL::GetOutFH($unfoldbin);


  my $counter=0;
  foreach my $s (@{$sitelist}) {
    my $chain=lc($s->{chain});
    my $cchain=$s->{chain};

    # loop over all the residues
    foreach my $r (@{$s->{reslist}}) {
      # for each residue create a new array to hold the energy
      my @temp;
      my $resname=$r->{resname};

      # for each residue, get the energies for all the  rotamers
      for (my $rot=0; $rot<=$#{$r->{rotamerlist}}; $rot++) {
        #push the energy below the crash value into @temp
        my $tempe=shift @{$energy};
        push (@temp,$tempe);
      }

      # now do ensemble average and print it into the file
      my $ensemble = GENERAL::Ensemble(\@temp);
      print $resname, "\n"  if (!defined $solv->{$resname});
      $ensemble->{'average'} += $solv->{$resname};

      my $index=-1;
      foreach (@{$r->{rotamerlist}}) {
        # we only want to print out the number of rotamers after backbone
        # clashing elimination. the number of 0.000 should be the same as the number of
        # pair energy terms
        if ($_->{valid} ==1 ) {
          # print the energy value into the file
          $counter++;
          $datf->printf("%8d%12.5f\n", $counter, $ensemble->{'average'});
          my $average= $ensemble->{'average'};
          $unfoldbinf->print(pack("d",  $average));
        }
      }
    }
  }

  # finally, we need to add a list of zero for pairwise term in unfoldenergy.dat file

  # loop over all the pairwise term
  for (my $s1=0; $s1 <=$#{$sitelist}; $s1++) {
    for (my $res1=0; $res1<=$#{$sitelist->[$s1]->{reslist}}; $res1++) {
      for (my $s2=$s1+1; $s2 <=$#{$sitelist}; $s2++) {
        for (my $res2=0; $res2<=$#{$sitelist->[$s2]->{reslist}}; $res2++) {
          for(my $rot1=0; $rot1<=$#{$sitelist->[$s1]->{reslist}->[$res1]->{rotamerlist}}; $rot1++) {
            next if ($sitelist->[$s1]->{reslist}->[$res1]->{rotamerlist}->[$rot1]->{valid} ==0 );
            for(my $rot2=0; $rot2<=$#{$sitelist->[$s2]->{reslist}->[$res2]->{rotamerlist}}; $rot2++) {
              next if ($sitelist->[$s2]->{reslist}->[$res2]->{rotamerlist}->[$rot2]->{valid} ==0 );
              # every residue pair has a pair.inp
              $counter++;
              #			    $datf->printf("%8d%14.5f\n", $counter, 0.0);
              $unfoldbinf->print(pack("d", 0.0));
            }
          }
        }
      }
    }
  }
  close($datf);
  close($unfoldbinf);

}


sub writeUnfoldData3 {
      
    my $rotamer=shift;
    my $propen_hash=shift;
    # now loop over all the design site
    my $sitelist=$rotamer->{dslist};
    &GENERAL::MakeDir("EnergyTable");
    chdir("EnergyTable");

    my $unfold='unfold-jeep'.'.dat';
    my $unfoldf=&GENERAL::GetOutFH($unfold);

    my $counter=0;
    #print "current counter: $counter ...\n"; 
    foreach my $s (@{$sitelist}) {
	my $chain=lc($s->{chain});
        my $cchain=$s->{chain};

	# loop over all the residues

	foreach my $r (@{$s->{reslist}}) {
  
	    # for each residue create a new array to hold the energy
	    my $resname=uc($r->{resname});

                foreach (@{$r->{rotamerlist}}) {
		if ($_->{valid} ==1 ) {
		    # print the energy value into the file
                    $counter++;
                    my $average= $propen_hash->{$resname}; 
                    print "$resname  value is not available!\n" if (!defined $average);
                    $unfoldf->printf("%8d%12.5f\n", $counter, $average);
                }

	    }
	}
    }

    close $unfoldf;
    chdir("../");

}



=head2

 Title   : readEnergyTerm
 Usage   : my $tab = ETABLE2::readEnergyTerm($term_track->{self}->{self_sa}, $param);
           ETABLE2::readEnergyTerm($term_track->{self}->{self_sa}, $param, $old_tab);
 Function: Reads the specified energy term. If a previous table is passed as the last
           parameter and it is not empty, the new term is either added (algebraically)
           or appended to the rows of the old table.
           If the specified table is empty, it is expanded to accomodate the newly read
           term. If no old table is specified, the new term is read into an array, a
           reference to which is returned.
           The actual values of the scale factors are taken from the param hash table
           passed to the function.
 Returns : an array of energies, 1 or 0 (see above).
 Args    : 1. reference to the term from the appropriate term tracker structure.
           2. Optional: the param hash table (where the values of the scale factors are given).
              If not specified, entries in rows will not be combined into one value and
              the entire table corresponding to the term will be returned.
           3. optional. Old energy table. If specified, the new term is added to the
              old values (unless the old table is empty, in which case it is treated
              as a table of zeros of the right size). This is convenient to use when
              adding large terms together (such as pair terms) because it is not
              neccessary to ever keep two terms in memory simultaneously. If old
              table is specified and no parameter hash table is given, the rows of the
              new term are appended to the old table instead of being added.

=cut

sub readEnergyTerm {
  my $term = shift;
  GENERAL::requireArgs($term);
  my $p = shift;
  my $ttab = shift;
  my $c;
  my $plen; # length of the previous table
  if (!defined($ttab)) {
    $c = 1;
    $ttab = ();
  } elsif (scalar(@$ttab) == 0) {
    $c = 2;
    $plen = 0;
  } else {
    $c = 3;
    $plen = scalar(@$ttab);
  }
  # Set up a local warning signal handler for non-numeric warnings
  local $^W = 1; # turn all warnings (the same as -w)
  local $SIG{__WARN__} = sub {
    if ($_[0] =~ /^Argument "(.*?)" isn't numeric/) {
      GENERAL::error("Value $1 is not numeric!");
    }
  };

  if ((defined($p)) && (!defined($term->{eschemes}->{$p->{escheme}}))) {
    GENERAL::warning("Warning in readEnergyTerm: it looks like term \"$term->{name}\" is not needed under " .
                     "the current energy scheme ($p->{escheme})!");
    return;
  }
  if ((!defined($term->{files})) || (scalar(@{$term->{files}}) == 0)) {
    GENERAL::error("It looks like term \"$term->{name}\" has not been calculated!");
  }

  my $sch;
  if (defined($p)) { $sch = $term->{eschemes}->{$p->{escheme}}; }
  my $i = 0; # Element index
  foreach my $file (@{$term->{files}}) {
    # All paths are defined relative to the root directory
    $file = File::Spec->rel2abs($file, $FS_DEF->{dhome});
    my $fh = GENERAL::GetInFH($file, 1);
    my $r;
    # Get scale factors for columns only once
    my @sfs;
    if (defined($p)) {
      for (my $k=0; $k < scalar(@{$sch->{cols_used}}); $k++) {
        # Get the right scale factor for this column
        push(@sfs, GENERAL::expr_eval($sch->{sf}->[$k], $p));
      }
    }
    while ($r = ETABLE2::readEnergyTableRow($fh, 1)) {
      my $val = 0;
      if (defined($sch->{fun})) {
        for (my $k=0; $k < scalar(@{$sch->{cols_used}}); $k++) {
          # Add the contribution of this column after it is scaled and the appropriate function is applied
          my $f = $sch->{fun}->[$k]; # function
          my $sval = $r->[$sch->{cols_used}->[$k]] * $sfs[$k]; # scaled value
          $f =~ s/\[arg\]/$sval/g; # use the scaled value as the function argument
          $val += eval($f); # compute the value of the function
        }
      } elsif (defined($p)) {
        for (my $k=0; $k < scalar(@{$sch->{cols_used}}); $k++) {
          # Add the contribution of this column after it is scaled
          $val += $r->[$sch->{cols_used}->[$k]] * $sfs[$k];
        }
      } else {
        # Convert strings to numbers
        for (my $k=0; $k < scalar(@$r); $k++) { $r->[$k] += 0.0; }
        $val = $r;
      }
      if ($c == 3) {
        if ($plen < $i+1) {
          GENERAL::error("Term \"$term->{name}\" has more rows (at least " . ($i+1) . ") than the previous term (" . scalar($plen) . " rows)!");
        }
        if (defined($p)) {
          $ttab->[$i] += $val;
        } else {
          push(@{$ttab->[$i]}, @$val);
        }
      } else {
        push(@{$ttab}, $val);
      }
      $i++;
#      if ($i % 2**18 == 0) {
#        print "$i: size = " . Devel::Size::size($ttab) . ", total size = " .  Devel::Size::total_size($ttab) . "\n";
#      }
    }
    close($fh);
  }
  if (($c == 3) && ($i < scalar(@$ttab))) {
    GENERAL::error("Term \"$term->{name}\" has fewer rows ($i) than the previous term (" . scalar(@$ttab) . " rows)!");
  }

  if ($c == 1) {
    return $ttab;
  } else {
    return $i;
  }

}



=head2

 Title   : getSelfTotal
 Usage   : ETABLE2::getSelfTotal()
 Function: Returns an array containing total self energies of all
           rotamers. The definition of the self energy (which terms
           go into it with which scale factors) is taken from the
           energy progress file. The actual values of the scale factors
           are taken from the param hash table passed to the function.
           This function also leaves information about which terms were
           used in a file such that external functions can log this for the
           user and also such that if scale factors change in stage 3 (after
           rotamers are eliminated) this function warns the user of the problem.
 Returns : an array of rotamer total self energies.
 Args    : 1. energy configuration hash table.
           2. optional. If set to 1, the function is going to assume that this
              is the first call, and it will create a tracking object file from scratch
              (overwriting the previous one, if it existed).
              If set to 2, the function performs a check to make sure
              that scale factors used now are the same as the ones used in the
              first run. If not, reports informative warnings. Also, it will update
              the "current" part of the tracking object file.
              If set to 0, the function does not bother doing anything related to the
              tracking structure. The default is 0.
           3. optional: reference to the energy track structure to use.

=cut

sub getSelfTotal {
  my $p = shift;
  GENERAL::requireArgs($p);
  my $s = shift;
  if (!defined($s)) { $s = 0; }
  my $th = shift;

  # Do we have anything to check?
  my $pc = {};
  if ($s == 2) {
    $pc = GENERAL::LoadStruct($FS_DEF->{self_ctrlf}, 1);
    # Make sure that all the scale factors previously used for rotamer elimination have not changed
    for (my $i=0; $i < scalar(@{$pc->{first}->{sf_names}}); $i++) {
      my $sf_name = $pc->{first}->{sf_names}->[$i];
      my $col_name = $pc->{first}->{col_names}->[$i];
      if (GENERAL::expr_eval($sf_name, $p) != $pc->{first}->{sf}->[$i]) {
        warn "\nWARNING\tWARNING\tWARNING\tWARNING\nWarning from getSelfTotal: you have changed the value of\n".
        "scale factor $sf_name (for column named \"$col_name\")\nfrom " .
        $pc->{first}->{sf}->[$i] . " to " . GENERAL::expr_eval($sf_name, $p) .
        ".\nThis change will be reflected upon the total calculated self energy.\n".
        "HOWEVER, it will NOT be reflected upon which rotamers are eliminated!!!\nWARNING\tWARNING\tWARNING\tWARNING\n\n";
      }
    }
  }

  # Get the term tracker structure and get the right entry
  if (!defined($th)) {
    $th = GENERAL::LoadStruct($FS_DEF->{term_trackerf}, 1);
  }
  my $en = $th->{self};

  # Now actually read in all tables and combine them
  my $table = [];
  my $rows = -1;
  my $pterm = -1;
  # Set up the "previous call" structure
  if (($s == 2) || ($s == 1)) {
    $pc->{current} = {};
    $pc->{current}->{col_names} = [];
    $pc->{current}->{sf} = [];
    $pc->{current}->{sf_names} = [];
    if ($s == 1) {
      $pc->{first}->{col_names} = [];
      $pc->{first}->{sf} = [];
      $pc->{first}->{sf_names} = [];
    }
  }
  # Go through each term
  foreach my $tkey (keys(%{$en})) {
    my $term = $en->{$tkey};
    # Is this term necessary for the current energy scheme?
    if (defined($term->{eschemes}->{$p->{escheme}})) {
#      my $ttab = ();
      my $sch = $term->{eschemes}->{$p->{escheme}};
      # Has the term been calculated?
      if (!defined($term->{files}) || (scalar(@{$term->{files}}) == 0)) {
        GENERAL::error("It looks like term \"$term->{name}\" has not been calculated although it is needed under the current ".
        "energy scheme ($p->{escheme})!");
      }
      # ------------ Document how this term is used --------------------
      if (($s == 2) || ($s == 1)) {
        for (my $ii=0; $ii < scalar(@{$sch->{cols_used}}); $ii++) {
          my $ci = $sch->{cols_used}->[$ii];
          my $sf = GENERAL::expr_eval($sch->{sf}->[$ii], $p);
          push(@{$pc->{current}->{col_names}}, $term->{col_names}->[$ci]);
          push(@{$pc->{current}->{sf}}, $sf);
          push(@{$pc->{current}->{sf_names}}, $sch->{sf}->[$ii]);
          if ($s == 1) {
            push(@{$pc->{first}->{col_names}}, $term->{col_names}->[$ci]);
            push(@{$pc->{first}->{sf}}, $sf);
            push(@{$pc->{first}->{sf_names}}, $sch->{sf}->[$ii]);
          }
        }
      }
      # ----------------------------------------------------------------

      # Read this term and add it to the total
      ETABLE2::readEnergyTerm($term, $p, $table);

    }
  }

  # Save information about this call
  if (($s == 2) || ($s == 1)) {
    GENERAL::SaveStruct($pc, $FS_DEF->{self_ctrlf});
  }
  return $table;
}


=head2

 Title   : getUnfoldTotal
 Usage   : ETABLE2::getUnfoldTotal()
 Function: Returns an array containing total unfold energies of all
           rotamers. The definition of the unfold energy (which terms
           go into it with which scale factors) is taken from the
           energy progress file. The actual values of the scale factors
           are taken from the param hash table passed to the function.
           This function also leaves information about which terms were
           used in a file such that external functions can log this for the
           user.
 Args    : 1. energy configuration hash table.
           2. optional. If set to 1, the function is going to create a tracking 
              object file (overwriting the previous one, if it existed).
              If set to 0, the function does not bother doing anything related to the
              tracking structure. The default is 0.
           3. optional: reference to the energy track structure to use.

=cut

sub getUnfoldTotal {
  my $p = shift;
  GENERAL::requireArgs($p);
  my $s = shift;
  if (!defined($s)) { $s = 0; }
  my $th = shift;

  my $pc = {};
  # Get the term tracker structure and get the right entry
  if (!defined($th)) {
    $th = GENERAL::LoadStruct($FS_DEF->{term_trackerf}, 1);
  }
  my $en = $th->{unfold};

  # Now actually read in all tables and combine them
  my $table = [];
  my $rows = -1;
  my $pterm = -1;
  # Set up the "previous call" structure
  $pc->{current} = {};
  $pc->{current}->{col_names} = [];
  $pc->{current}->{sf} = [];
  # Go through each term
  foreach my $tkey (keys(%{$en})) {
    my $term = $en->{$tkey};
    # Is this term necessary for the current energy scheme?
    if (defined($term->{eschemes}->{$p->{escheme}})) {
#      my $ttab = ();
      my $sch = $term->{eschemes}->{$p->{escheme}};
      # Has the term been calculated?
      if (scalar(@{$term->{files}}) == 0) {
        die "Error in getUnfoldTotal: it looks like term \"$term->{name}\" ".
        "has not been calculated although it is needed under the current ".
        "energy scheme ($p->{escheme})!\n";
      }
      # ------------ Document how this term is used --------------------
      for (my $ii=0; $ii < scalar(@{$sch->{cols_used}}); $ii++) {
        my $ci = $sch->{cols_used}->[$ii];
        my $sf = GENERAL::expr_eval($sch->{sf}->[$ii], $p);
        push(@{$pc->{current}->{col_names}}, $term->{col_names}->[$ci]);
        push(@{$pc->{current}->{sf}}, $sf);
      }
      # ----------------------------------------------------------------

      # Read this term and add it to the total
      ETABLE2::readEnergyTerm($term, $p, $table);

    }
  }

  # Save information about this call
  if ($s) {
    GENERAL::SaveStruct($pc, $FS_DEF->{unf_ctrlf});
  }

  return $table;
}


=head2

 Title   : getPairTotal
 Usage   : ETABLE2::getPairTotal()
 Function: Returns an array containing total pair energies of all valid
           rotamer pairs. The definition of the pair energy (which terms
           go into it with which scale factors) is taken from the
           energy progress file. The actual values of the scale factors
           are taken from the param hash table passed to the function.
           This function also leaves information about which terms were
           used in a file such that external functions can log this for the
           user.
 Args    : 1. energy configuration hash table.
           2. optinal: the number of elements to be read. If specified,
              array of this size is pre-allocated, which should save time for
              large tables.
           3. optional: If set to 1, the function is going to create a tracking 
              object file (overwriting the previous one, if it existed).
              If set to 0, the function does not bother doing anything related to the
              tracking structure. The default is 0.
           4. optional: verbose flag. Default is true (be verbose).
           5. optional: reference to the energy track structure to use.

=cut

sub getPairTotal {
  my $p = shift;
  GENERAL::requireArgs($p);
  my $size = shift;
  my $s = shift;
  if (!defined($s)) { $s = 0; }
  my $verbose = shift;
  if (!defined($verbose)) { $verbose = 1; }
  my $th = shift;

  my $pc = {};

  # Get the term tracker structure and get the right entry
  if (!defined($th)) {
    $th = GENERAL::LoadStruct($FS_DEF->{term_trackerf}, 1);
  }
  my $en = $th->{pair};

  # Now actually read in all tables and combine them
  # Set up the "previous call" structure
  $pc->{current} = {};
  $pc->{current}->{col_names} = [];
  $pc->{current}->{sf} = [];
  my @table = ();
  if (defined($size)) {
    print "Allocating memory for table...\n" unless(!$verbose);
    $#table = ($size-1);
    print "Done allocating memory for table...\n" unless(!$verbose);
  }
  my $table = \@table;
  my $N = 0; # number of rows read
  # Go through each term
  foreach my $tkey (keys(%{$en})) {
    my $term = $en->{$tkey};
    # Is this term necessary for the current energy scheme?
    if (defined($term->{eschemes}->{$p->{escheme}})) {
#      my $ttab = ();
      my $sch = $term->{eschemes}->{$p->{escheme}};
      # Has the term been calculated?
      if (!defined($term->{files}) || (scalar(@{$term->{files}}) == 0)) {
        GENERAL::error("It looks like term '$term->{name}' has not been calculated although it is needed under the current energy scheme ($p->{escheme})!");
      }
      # ------------ Document how this term is used --------------------
      for (my $ii=0; $ii < scalar(@{$sch->{cols_used}}); $ii++) {
        my $ci = $sch->{cols_used}->[$ii];
        my $sf = GENERAL::expr_eval($sch->{sf}->[$ii], $p);
        push(@{$pc->{current}->{col_names}}, $term->{col_names}->[$ci]);
        push(@{$pc->{current}->{sf}}, $sf);
      }
      # ----------------------------------------------------------------

      # Read this term and add it to the total
      $N = ETABLE2::readEnergyTerm($term, $p, $table);
    }
  }

  # Save information about this call
  if ($s) {
    GENERAL::SaveStruct($pc, $FS_DEF->{pair_ctrlf});
  }
  return ($table, $N);
}


sub writeSelfData {

    my $ecut=shift;
    $ecut=100.0 if (!defined $ecut);
    my $suffix=shift;
    $suffix='default' if (!defined $suffix);
    my $choice=shift;
    $choice=1 if (!defined $choice);
    my $eterms=&ETABLE2::getEterms('self', $choice);
    my $SelfEnergyTable=shift;
    $SelfEnergyTable = 'SelfEnergyTable' if (!defined $SelfEnergyTable);
    my $diel=shift;
    $diel="" if (!defined $diel || $diel == 4); 

    my $bige=1000000.0;
    chdir("EnergyTable");

    my $fname1="self".$diel."-"."$suffix".'.dat';
    my $edat=&GENERAL::GetOutFH($fname1);
    my $fname2="selftotal".$diel."-"."$suffix".'.dat';
    my $total=&GENERAL::GetOutFH("$fname2");

    my $sol=map {/5/} @{$eterms};

    # readin energy Table
    my $earray=&ETABLE2::readEnergyTable($SelfEnergyTable);
    my $energy=&ETABLE2::pickEnergyTerms($earray, $eterms);
    my $newenergy=$energy;
   
    my @solv=[]; 
    if ($sol) {
        my $fh=&GENERAL::GetInFH("solv.dat");
        @solv=<$fh>;
        chomp @solv;
        $newenergy=&GENERAL::MergeAddArray($energy, \@solv);
     }

    my $fname4="../Design/rot_states.self";
    my $rotstate=&GENERAL::GetInFH($fname4);
    my @states=();

    # collect rotamer state info
    while (<$rotstate>) {
        chomp $_;
        if ( $_ ne "") {
	    my @lines=split(" ", $_); 
            push (@states, @lines);
        }
    } 

  
    # now print total energy table and after clash check table    
    foreach (@{$newenergy}) {
         my $totale= $_;
	 my $sflag=shift @states;
         if ($sflag == 1) {
            $edat->printf("%14.5f\n", $totale);
         }

         $total->printf("%3d %14.5f\n", $sflag, $totale);
    }
   close $total;
   close $edat;

    chdir("../");
}


=head2 

 Title   : pickEnergyTerms2
 Usage   : 
 Function: given a total energy table array and enery  
           terms we want to combine. produce a new 
           total energy array. same as pickEnergyTerms
           except we add some coeffciency for eef and diel 
           only works for choice 2:  multe + vdw + eef
 Returns : 
 Args    : (1) reference to energy table
           (2) column to pick 
           (4) value of coeffciency:
                     vdw
                     multe
                     eef   
=cut


sub pickEnergyTerms2 {
  my $totalenergy=shift;
  my $vdw=shift;
  my $multe=shift;
  my $eef=shift;

  my @newarray;
  for(my $i=0; $i<=$#{$totalenergy}; $i++) {
    my $line=$totalenergy->[$i];
    my $totale=0.0;
    # loop over all the energy terms we want to add up
    $totale = $vdw*$line->[1] + $multe*$line->[2] + $line->[3] + $eef*$line->[5];
    push (@newarray, $totale);
  }

  return \@newarray;
}

=head2 
 Title   :  writeSelfData2
 Usage   :  writeSelfData2($ecut, 'SelfEnergyTable', $vdw, $diel, $eef);

=cut



sub writeSelfData2 {

    my $rot = shift;
    my $parm = shift;

    my $eterms = ETABLE2::getEterms('self', 2);
    my $vdw=$parm->{vdw};
    my $eef=$parm->{eef};
    my $diel=$parm->{dde};

    my $suffix = sprintf("%.1f_%s_%.1f", $vdw, $diel, $eef);

    $vdw=1.0 if (!defined $vdw);
    $diel=4.0 if (!defined $diel);
    $eef=1.0 if (!defined $eef);
    my $multe=4.0/$diel;

    my $fname1 = $FS_DEF->{self_datf};
    $fname1 =~ s/%/$suffix/;
    my $edat = GENERAL::GetOutFH($fname1);

    # readin energy Table
    my $earray = ETABLE2::readEnergyTable($FS_DEF->{self_tabf});
    my $energy = ETABLE2::pickEnergyTerms2($earray, $vdw, $multe, $eef);
   
    my @solv=[]; 
    my $fh = GENERAL::GetInFH($FS_DEF->{eef_solv_tabf});
    @solv = <$fh>;
    chomp @solv;
    my $newenergy = GENERAL::MergeAddArray($energy, \@solv);

    # now print total energy table and after clash check table
    my $i = 0;
    foreach my $site ($rot->{dslist}) {
      foreach my $res ($site->{reslist}) {
        foreach my $rot ($res->{rotamerlist}) {
          if ($rot->{valid}) {
            $edat->printf("%14.5f\n", $newenergy->[$i]);
          }
          $i++;
        }
      }
    }
    close($edat);
}




sub writePairData {
  my $choice = shift;
  my $eterms = ETABLE2::getEterms('pair', $choice);

  my $fname = $FS_DEF->{pair_datf};
  $fname =~ s/%/$choice/;
  my $datf = GENERAL::GetOutFH($fname);

  # readin energy Table
  my $etable = GENERAL::GetInFH($FS_DEF->{pair_tabf});

  while (<$etable>) {
    next if (/^\#/);
    next if ( /^\s+$/);
    chomp $_;
    my @lines=split(" ", $_);
    my $totale=0.0;

    foreach my $index (@{$eterms}) {
      $totale += $lines[$index];
    }

    $datf->printf("%14.5f\n", $totale);
  
  }
  close $datf;
  close $etable;
}


=head2

 Usage:  writePairData2('pairEnergyTable', $vdw, $multe, $eef);


=cut

sub writePairData2 { 

    my $tablefname=shift;
    $tablefname='pairEnergyTable' if (!defined $tablefname);
    my $parm=shift;
    my $vdw=$parm->{vdw};
    my $eef=$parm->{eef};
    my $diel=$parm->{dde};
    my $suffix=$parm->{suffix};

    $vdw=1.0 if (!defined $vdw);
    $diel=4.0 if (!defined $diel);
    $eef=1.0 if (!defined $eef);
    my $multe=1.0;

 
    chdir("EnergyTable");

    my $fname='pair'.'-'.$suffix.'.dat';
    my $datf=&GENERAL::GetOutFH($fname);

    # readin energy Table

    $tablefname='pairEnergyTable' if (!defined $tablefname);
    my $etable=&GENERAL::GetInFH($tablefname);
 
    while (<$etable>) {
        next if (/^\#/);
        next if ( /^\s+$/);
        chomp $_;
        my @lines=split(" ", $_);
        my $totale=0.0;

        $totale = $vdw*$lines[0] + $multe*$lines[1] + $eef*$lines[2];

        # Assign total energy  upper limit as 500.0
        if ($totale >= 500.0) {
            $totale=500.0;
         }

         $datf->printf("%14.5f\n", $totale);

    }

    chdir("../");
}



=head2 

 Title   : buildHybridPotential
 Usage   : ETABLE2::buildHybridPotential('coupled',0.2)
           or ETABLE2::buildHybridPotential('uncoupled',0.2)
 Function: two choices to build hybrid potential 
           (A) coupled potential:  (1-lamda)*physical + lamda)* stat
           (B) uncoupled potential: physical + lamda*stat
 Returns : 
 Args    : 

=cut


sub  buildHybridPotential {

    my $choice=shift;
    my $lamda=shift;
    my $suffix1=shift;
    my $suffix2=shift;

    $suffix1='default' if (!defined $suffix1);
    my $newlamda=int(100*$lamda);
    $suffix2=$choice.'-'.$newlamda if (!defined $suffix2);

    my $selfdata='self-'."$suffix1".'.dat';
    my $fname='self-'."$suffix2".'.dat';

    chdir("EnergyTable");
    my $selfe=&GENERAL::GetInFH($selfdata);
    my $stat=&GENERAL::GetInFH("StatEnergy.tab");
    my $out=&GENERAL::GetOutFH($fname);


    # dump the data into array
    my @selfe=<$selfe>;
    chomp (@selfe);  
    my @stat=<$stat>;
    chomp(@stat);

    my $newarray;
    if (lc($choice) eq 'coupled') {
	my $lamda1=1.0-$lamda;
	$newarray=&GENERAL::MergeAddArray(\@selfe,\@stat, $lamda1, $lamda);
    } elsif ( lc($choice) eq 'uncoupled') {
	$newarray=&GENERAL::MergeAddArray(\@selfe,\@stat, 1.0, $lamda);
    }

    		    
    foreach (@{$newarray}) { 
	$out->printf("%14.5f\n", $_);
    }

    chdir("../");
}



=head2 

 Title   : packEnergyBin
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub packEnergyBin {

    my $suffix1=shift;
    my $suffix2=shift;
    my $suffix3=shift;
    my $diel=shift;

    $suffix1='default' if (!defined $suffix1);
    $suffix2='default' if (!defined $suffix2);
    $suffix3='default' if (!defined $suffix3);
    $diel="" if (!defined $diel || $diel==4);

    chdir("EnergyTable");
    my $fname1='self'.$diel.'-'.$suffix1.'.dat';
    my $fname2='pair'.$diel.'-'.$suffix2.'.dat';
    my $fname3='energy'.$diel.'-'.$suffix3.'.dat';
    my $fname4='energy'.$diel.'-'.$suffix3.'.bin';

    my $fh1=&GENERAL::GetInFH($fname1);
    my $fh2=&GENERAL::GetInFH($fname2);
    my $out=&GENERAL::GetOutFH($fname3);
    my $out2=&GENERAL::GetOutFH($fname4);


    my $num=0;
     while (<$fh1>) {
	$num++;
        $_=&GENERAL::Trim($_);
	$out->printf("%5d %14.5f\n",$num, $_);
	$out2->print(pack("d", $_));
   }

    while (<$fh2>) {
	$num++;
	$num=$num % 100000;
        $_=&GENERAL::Trim($_);
	$out->printf("%5d %14.5f\n",$num, $_);
	$out2->print(pack("d", $_));
    }

   close $fh1;
   close $fh2;
   close $out2;
   close $out;
    chdir("../");

}


=head2

 Title   : packScaledEnergyBin 
 Usage   : packScaledEnergyBin(1.0, 16, 0.8);
 Function: 
 Returns :
 Args    :

=cut




sub packScaledEnergyBin {

    my $parm=shift;
    my $vdw=$parm->{vdw};
    my $eef=$parm->{eef};
    my $diel=$parm->{dde};
    my $suffix=$parm->{suffix};

    $vdw=1.0 if (!defined $vdw);
    $diel=4.0 if (!defined $diel);
    $eef=1.0 if (!defined $eef);
    my $multe=4.0/$diel;

    my $unfoldbit;
    if (lc($parm->{'unfolded_state_calc'}) ne  'no') {
        $unfoldbit=1;
    } else {
        $unfoldbit=0;
    }
 
    my $dbit=shift;
    $dbit=0 if (!defined $dbit);

    chdir("EnergyTable");
    my $fname1;
    if ($unfoldbit==0 ) {
    	$fname1="self"."-".$suffix.'.dat';
    } elsif ($unfoldbit==1) {
#        $fname1="selfunfold"."-".$suffix.'.dat';
         $fname1="self"."-".$suffix.'.dat';
    } else {
    }

    my $fname2='pair'."-".$suffix.'.dat';
    my $fname4='energy'."-".$suffix.'.bin';

    my $fh1=&GENERAL::GetInFH($fname1);
    my $fh2=&GENERAL::GetInFH($fname2);
    my $fh3=&GENERAL::GetInFH('BetaDistance') if ($dbit==1);
    my $out2=&GENERAL::GetOutFH($fname4);


     while (<$fh1>) {
        $_=&GENERAL::Trim($_);
        $out2->print(pack("d", $_));
   }

    while (<$fh2>) {
        $_=&GENERAL::Trim($_);
        $out2->print(pack("d", $_));
    }


    if ($dbit==1) {
         $out2->print(pack("CCC", 255, 255, 255));

        while (<$fh3>) {
            $_=&GENERAL::Trim($_);
            $out2->print(pack("d", $_));
       }

    }

   close $fh1;
   close $fh2;
   close $out2;
   close $fh3 if ($dbit==1);
   chdir("../");

}







=head2 

 Title   :  packEnergyBin2
 Usage   : 
 Function:  pack the energy bin file with ditance
 Returns : 
 Args    : 

=cut



sub packEnergyDistanceBin {

    my $suffix1=shift;
    my $suffix2=shift;
    my $suffix3=shift;
    my $diel=shift;

    $suffix1='default' if (!defined $suffix1);
    $suffix2='default' if (!defined $suffix2);
    $suffix3='default' if (!defined $suffix3);
    $diel="" if (!defined $diel || $diel==4);

    chdir("EnergyTable");
    my $fname1='self'.$diel.'-'.$suffix1.'.dat';
    my $fname2='pair'.$diel.'-'.$suffix2.'.dat';
    my $fname3='energy'.$diel.'-'.$suffix3.'_distance.dat';
    my $fname4='energy'.$diel.'-'.$suffix3.'_distance.bin';

    my $fh1=&GENERAL::GetInFH($fname1);
    my $fh2=&GENERAL::GetInFH($fname2);
    my $fh3=&GENERAL::GetInFH('BetaDistance');
    my $out=&GENERAL::GetOutFH($fname3);
    my $out2=&GENERAL::GetOutFH($fname4);


    my $num=0;
     while (<$fh1>) {
        $num++;
        $_=&GENERAL::Trim($_);
        $out->printf("%5d %14.5f\n",$num, $_);
        $out2->print(pack("d", $_));
   }

    while (<$fh2>) {
        $num++;
        $num=$num % 100000;
        $_=&GENERAL::Trim($_);
        $out->printf("%5d %14.5f\n",$num, $_);
        $out2->print(pack("d", $_));
    }

#    $out2->print(pack("B", "111111111111111111111111"));     
    $out2->print(pack("CCC", 255, 255, 255));

    while (<$fh3>) {
        $num++;
        $num=$num % 100000;
        $_=&GENERAL::Trim($_);
        $out->printf("%5d %14.5f\n",$num, $_);
        $out2->print(pack("d", $_));
    }


    chdir("../");

}



sub packMan {

    # self file	
    my $fname1=shift;
    # pair file
    my $fname2=shift;
    # binary file
    my $fname3=shift;
	
   chdir("EnergyTable");

   my $fh1=&GENERAL::GetInFH($fname1);
   my $fh2=&GENERAL::GetInFH($fname2);
   my $out=&GENERAL::GetOutFH($fname3);

   while (<$fh1>) {
        $_=&GENERAL::Trim($_);
	$out->print(pack("d", $_));
   }

   while (<$fh2>) {
        $_=&GENERAL::Trim($_);
        $out->print(pack("d", $_));
    }

   close $fh1;
   close $fh2;
   close $out;
   chdir("../");

}

1;
