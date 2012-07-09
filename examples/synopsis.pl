use Algorithm::GaussianElimination::GF2;

my $age = Algorithm::GaussianElimination::GF2->new;
$age->new_equation(1, 0, 0, 1 => 1);
$age->new_equation(0, 0, 1, 1 => 0);
my ($sol, @base0) = $age->solve;
print "sol: @$sol\n";

# or you can also create the equations setting elements at given
# positions:

my $age = Algorithm::GaussianElimination::GF2->new;
my $eq1 = $age->new_equation;
$eq1->a(0, 1);
$eq1->a(3, 1);
$eq1->b(1);
my $eq2 = $age->new_equation;
$eq2->a(2, 1);
$eq2->a(3, 1);
$eq2->b(0);
my ($sol, @base0) = $age->solve;

print "sol: @$sol\n";
