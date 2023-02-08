#!/usr/bin/perl -w
use strict;
use warnings;
our %detector;
our %configuration;
our %parameters;

use Getopt::Long;
use Math::Trig;

my $DetectorName = 'solid_DDVCS_muon_barrel';

my $DetectorMother="root";

sub solid_DDVCS_muon_barrel
{
make_solid_DDVCS_muon_barrel();
}

 my $color="00ffff";
 my $material="G4_Fe";
 my $color_scint="ff0000"; 
 my $material_scint="G4_PLASTIC_SC_VINYLTOLUENE"; 

sub make_solid_DDVCS_muon_barrel
{
    my %detector=init_det();
    
    $detector{"name"}        = "$DetectorName\_1_virt";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm 0*cm";
    $detector{"rotation"}   = "0*deg 0*deg 22.5*deg";
    $detector{"color"}      = "CC6633";
    $detector{"type"}       = "Pgon";
    $detector{"dimensions"} = "0*deg 360*deg 8*counts 2*counts 280*cm 280*cm 280.1*cm 280.1*cm -266*cm 182*cm";
    $detector{"material"}   = "G4_Galactic";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "flux";
    $detector{"hit_type"}    = "flux";
    $detector{"identifiers"} = "id manual 6310000";
    print_det(\%configuration, \%detector);
    
    $detector{"name"}        = "$DetectorName\_1_scint";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm 0*cm";
    $detector{"rotation"}   = "0*deg 0*deg 22.5*deg";
    $detector{"color"}      = $color_scint;
    $detector{"type"}       = "Pgon";
    $detector{"dimensions"} = "0*deg 360*deg 8*counts 2*counts 285*cm 285*cm 290*cm 290*cm -266*cm 182*cm";
    $detector{"material"}   = "$material_scint";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "flux";
    $detector{"hit_type"}    = "flux";
    $detector{"identifiers"} = "id manual 6301000";
    print_det(\%configuration, \%detector);
}
