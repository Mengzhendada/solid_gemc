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
 my $color_scint="ff0000"; 
 my $material_scint="G4_PLASTIC_SC_VINYLTOLUENE"; 

sub make_solid_DDVCS_muon_largeangle
{
    my $z=(209+570+50)/2;
    my $Dz = (570+50-209)/2;
    
    my %detector=init_det();
    
    $detector{"name"}        = "$DetectorName\_1_virt";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $z*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "CC6633";
    $detector{"type"}       = "Tube";
    $detector{"dimensions"}  = "287*cm 287.1*cm $Dz*cm 0*deg 360*deg";
    $detector{"material"}   = "G4_Galactic";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "flux";
    $detector{"hit_type"}    = "flux";
    $detector{"identifiers"} = "id manual 6210000";
    print_det(\%configuration, \%detector);
    
    $detector{"name"}        = "$DetectorName\_1_scint";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $z*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = $color_scint;
    $detector{"type"}       = "Tube";
    $detector{"dimensions"}  = "290*cm 295*cm $Dz*cm 0*deg 360*deg";
    $detector{"material"}   = "$material_scint";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "flux";
    $detector{"hit_type"}    = "flux";
    $detector{"identifiers"} = "id manual 6201000";
    print_det(\%configuration, \%detector);
}
