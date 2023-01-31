#!/usr/bin/perl -w
use strict;
use warnings;
our %detector;
our %configuration;
our %parameters;

use Getopt::Long;
use Math::Trig;

my $DetectorName = 'solid_DDVCS_muon_largeangle';

my $DetectorMother="root";

sub solid_DDVCS_muon_largeangle
{
make_solid_DDVCS_muon_largeangle();
}

 my $color="00ffff";
 my $material="G4_Fe";

sub make_solid_DDVCS_muon_largeangle
{
    my $z=(209+570+50)/2;
    my %detector=init_det();
    $detector{"name"}        = "$DetectorName";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $z*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = $color;
    $detector{"type"}       = "Cons";
    my $Rmin1 = 290;
    my $Rmax1 = 326;
    my $Rmin2 = 290;
    my $Rmax2 = 326;
    my $Dz    = (570+50-209)/2;
    my $Sphi  = 0;
    my $Dphi  = 360;
    $detector{"dimensions"}  = "$Rmin1*cm $Rmax1*cm $Rmin2*cm $Rmax2*cm $Dz*cm $Sphi*deg $Dphi*deg";
    $detector{"material"}   = "$material";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "no";
    $detector{"hit_type"}    = "no";
    $detector{"identifiers"} = "no";
    print_det(\%configuration, \%detector);
}
