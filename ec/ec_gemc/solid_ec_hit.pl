use strict;
use warnings;

our %configuration;
our %parameters;

sub define_solidec_hit
{
	# uploading the hit definition
	my %hit = init_hit();
	$hit{"name"}            = "solidec";
	$hit{"description"}     = "solid ec hit definition";
	$hit{"identifiers"}     = "id";
	$hit{"signalThreshold"} = "0.5*MeV";
	$hit{"timeWindow"}      = "4*ns";
	$hit{"prodThreshold"}   = "1*mm";
	$hit{"maxStep"}         = "1*cm";
	$hit{"delay"}           = "10*ns";
	$hit{"riseTime"}        = "1*ns";
	$hit{"fallTime"}        = "1*ns";
	$hit{"mvToMeV"}         = 100;
	$hit{"pedestal"}        = -20;
	print_hit(\%configuration, \%hit);
}


1;