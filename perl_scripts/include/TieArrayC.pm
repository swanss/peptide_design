package TieArrayC;

use ArrayC;
use strict;
use warnings;
use Carp;

# we don't need Tie::Array anymore -- all methods implemented!
# use base 'Tie::Array';

sub DESTROY {
  freeArrayC($_[0]->{ptr}) if ($_[0]->{ptr} ne 0);
}

sub TIEARRAY {
  my $class = shift;
  my $type  = shift;
  my $size = shift; $size = 0 if (!defined($size));
  croak $@ if $@;
  my ($itype, $get_function_ptr, $set_function_ptr, $null_value);

  # parameter error check
  croak("Unknown type '$type'") if ($type !~ /^int|double$/i);
  if ($type =~ /^int$/i) {
    $itype = 2;
    $get_function_ptr = \&accessIntArrayElementC;
    $set_function_ptr = \&setIntArrayElementC;
    $null_value = 0;
  } elsif ($type =~ /^double$/i) {
    $itype = 1;
    $get_function_ptr = \&accessDoubleArrayElementC;
    $set_function_ptr = \&setDoubleArrayElementC;
    $null_value = 0.0;
  } else {
    croak("Unrecognized array type '$type'");
  }
  croak("Size must be integer, but is $size") if (!isInteger($size));

  # allocate C array if size was specified
  my $ptr = 0;
  $ptr = allocateArrayC($size, $itype) if ($size > 0);

  bless {
      ptr    => $ptr,
      type   => $type,
      itype  => $itype,
      size   => $size,
      get    => $get_function_ptr,
      set    => $set_function_ptr,
      undefv => $null_value,
  };
}

sub FETCH {
  return $_[0]->{get}->($_[0]->{ptr}, $_[1]);
}

sub FETCHSIZE {
  return $_[0]->{size};
}

sub STORE {
  my ($this, $index, $value) = @_;

  # resize if needed
  if ($this->FETCHSIZE() - $index < 1) {
    $this->STORESIZE($index + 1);
  }

  # store
  $this->{set}->($this->{ptr}, $index, $value);
}

sub STORESIZE {
  my ($this, $count) = @_;

  if ($this->FETCHSIZE() != $count) {
    if ($count == 0) {
      freeArrayC($this->{ptr});
      $this->{ptr} = 0;
    } else {
      $this->{ptr} = resizeArrayC($this->{ptr}, $count, $this->{itype});
    }
    $this->{size} = $count;
  }
  return $count;
}

sub EXISTS {
  my ($this, $key) = @_;
  return ($this->FETCHSIZE() > $key);
}

sub EXTEND {
}

sub DELETE {
  my ($this, $index) = @_;

  # zero down
  $this->{set}->($this->{ptr}, $index, $this->{undefv});
}

sub CLEAR {
  my $this = shift;
  $this->STORESIZE(0);
}

sub PUSH {    # append
  my $this = shift;

  my $oldsize = $this->{size};
  $this->STORESIZE($this->{size} + scalar(@_));
  for (my $i = 0; $i < scalar(@_); $i++) {
    $this->{set}->($this->{ptr}, $oldsize + $i, $_[$i]);
  }
}

sub UNSHIFT {    # prepend
  my $this = shift;

  if (scalar(@_) > 0) {
    my $oldsize = $this->{size};
    $this->STORESIZE($this->{size} + scalar(@_));
    for (my $i = $oldsize-1; $i >= 0; $i--) {
      $this->{set}->($this->{ptr}, $i + scalar(@_), $this->{get}->($this->{ptr}, $i));
    }
    for (my $i = 0; $i < scalar(@_); $i++) {
      $this->{set}->($this->{ptr}, $i, $_[$i]);
    }
  }
}

sub POP {
  my $this = shift;

  my $val;
  my $newsize = $this->FETCHSIZE() - 1;
  if ($newsize >= 0) {
    $val = $this->{get}->($this->{ptr}, $this->{size}-1);
    $this->STORESIZE($newsize);
  }
  return $val;
}

sub SHIFT {
  my $this = shift;

  my $val;
  my $newsize = $this->FETCHSIZE() - 1;
  if ($newsize >= 0) {
    $val = $this->{get}->($this->{ptr}, $this->{size}-1);
    for (my $i = 0; $i < $this->{size}-1; $i++) {
      $this->{set}->($this->{ptr}, $i, $this->{get}->($this->{ptr}, $i+1));
    }
    $this->STORESIZE($newsize);
  }
  return $val;
}

sub isInteger {
  my $num = shift;
  return ($num =~ /(^\s*[-+]?[1-9][0-9]*\s*$|^\s*[-+]?0\s*$)/) ? 1:0;
}


1;
__END__

