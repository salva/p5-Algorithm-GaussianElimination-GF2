package Algorithm::GaussianElimination::GF2;

our $VERSION = '0.01';

use strict;
use warnings;

sub new {
    my $class = shift;
    my $self = { eqs => [] };
    bless $self, $class;
}

sub add_equation {
    my ($self, $eq) = @_;
    push @{$self->{eqs}}, $eq;
}

sub new_equation {
    my $self = shift;
    my $eq = Algorithm::GaussianElimination::GF2::Equation->new(@_);
    $self->add_equation($eq);
    $eq;
}

sub _first_1 {
    $_[0] =~ /[^\0]/g or return undef;
    my $end = pos($_[0]) * 8 - 1;
    for my $i (($end - 7) .. $end) {
        return $i if vec($_[0], $i, 1);
    }
}

sub dump {
    my $self = shift;
    my $eqs = $self->{eqs};
    my $len = 0;
    for (@$eqs) {
        $len = $_->[2] if $_->[2] > $len;
    }
    printf "GF(2) system of %d equations and %d variables\n", scalar(@$eqs), $len;
    for (@$eqs) {
        $_->[2] = $len;
        $_->dump;
    }
    print "\n";
}

sub solve {
    my $self = shift;
    my $eqs = $self->{eqs};
    my $len = 0;
    for my $i (0..$#$eqs) {
        my $pivot = $eqs->[$i];
        $len = $pivot->[2] if $pivot->[2] > $len;
        my $v = $pivot->[0];
        my $ix = _first_1($v);
        if (defined $ix) {
            for my $j (($i+1)..$#$eqs) {
                my $target = $eqs->[$j];
                if (vec($target->[0], $ix, 1)) {
                    $target->[0] ^= $v;
                    $target->[1] ^= $pivot->[1];
                }
            }
        }
        else {
            # unsolvable!
            return if $pivot->[1];
        }
    }
    $_->[2] = $len for @$eqs;

    my @free;
    my @sol;
    for my $eq (reverse @$eqs) {
        my $v = $eq->[0];
        my $ix = _first_1($v);
        if (defined $ix) {
            my $b = $eq->[1];
            for my $i (($ix + 1) .. ($len - 1)) {
                if (vec $v, $i, 1) {
                    if (defined $sol[$i]) {
                        $b ^= $sol[$i];
                    }
                    else {
                        push @free, $i;
                        $sol[$i] = 0;
                    }
                }
            }
            $sol[$ix] = $b;
        }
    }
    for my $i (0..$len - 1) {
        unless (defined $sol[$i]) {
            push $i, @free;
            $sol[$i] = 0;
        }
    }
    return \@sol unless wantarray;
    my @base0;
    for my $free (@free) {
        my @sol0;
        $sol0[$_] = 0 for @free;
        $sol0[$free] = 1;
        for my $eq (reverse @$eqs) {
            my $v = $eq->[0];
            my $ix = _first_1($v);
            if (defined $ix) {
                my $b = 0;
                for my $i (($ix + 1) .. ($len - 1)) {
                    if (vec $v, $i, 1) {
                        defined $sol0[$i] or die "internal error: unexpected free var found";
                        $b ^= $sol0[$i];
                    }
                }
                $sol0[$ix] = $b;
            }
        }
        push @base0, \@sol0;
    }
    return \@sol, @base0;
}

package Algorithm::GaussianElimination::GF2::Equation;

sub new {
    my $class = shift;
    my $self = ['', 0, 0];
    bless $self, $class;
    if (@_) {
        $self->[1] = pop @_;
        for my $ix (0..$#_) {
            vec($self->[0], $ix, 1) = $_[$ix]
        }
        $self->[2] = @_;
    }
    $self
}

sub a {
    my ($self, $ix, $v) = @_;
    if (defined $v) {
        $self->[2] = $ix + 1 if $self->[2] <= $ix;
        return vec($self->[0], $ix, 1) = $v;
    }
    return vec($self->[0], $ix, 1);
}

sub b {
    my ($self, $v) = @_;
    if (defined $v) {
        return $self->[1] = ($v ? 1 : 0);
    }
    return $self->[1];
}

sub dump {
    my $self = shift;
    my $last = $self->[2] - 1;
    my @a = map vec($self->[0], $_, 1), 0.. $last;
    print "@a | $self->[1]\n";
}

1;
__END__

=head1 NAME

Algorithm::GaussianElimination::GF2 - Solve linear systems of equations on GF(2)

=head1 SYNOPSIS

  use Algorithm::GaussianElimination::GF2;

  my $age = Algorithm::GaussianElimination::GF2->new;
  $age->add_equation(1, 0, 0, 1 => 1);
  $age->add_equation(0, 0, 1, 1 => 0);
  my ($sol, @base0) = $age->solve;

  # or you can also create the equations setting elements at given
  # positions:

  my $age = Algorithm::GaussianElimination::GF2->new;
  my $eq1 = $age->add_equation;
  $eq1->a(0, 1);
  $eq1->a(3, 1);
  $eq1->b(1);
  my $eq2 = $age->add_equation;
  $eq2->a(3, 1);
  $eq2->a(4, 1);
  $eq2->b(0);
  my ($sol, @base0) = $age->solve;


=head1 DESCRIPTION

This module implements a variation of the Gaussian Elimination
algorithm over GF(2).

=head2 API

Those are the interesting methods:

=over 4

=item ($sol, @base0) = $age->solve

=item $sol = $age->solve

This method solves the system of equations.

If it is unsolvable returns an empty list.

If is is solvable and deterministic returns and array reference containing the solution

If it is solvable and undeterministic it returns one solution and a base to the vector space of the null system.

=back

=head2 EXPORT

None by default.



=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 AUTHOR

Salvador Fandino, E<lt>salva@E<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2011 by Salvador Fandino

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.14.2 or,
at your option, any later version of Perl 5 you may have available.


=cut
