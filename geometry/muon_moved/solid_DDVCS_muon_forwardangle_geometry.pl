#!/usr/bin/perl -w
use strict;
use warnings;
our %detector;
our %configuration;
our %parameters;

use Getopt::Long;
use Math::Trig;

my $DetectorName = 'solid_DDVCS_muon_forwardangle';

my $DetectorMother="root";

# 1 of 8 side of 3rd layer of CLEO is roughly 98.55"x14.21"x210" (250x36x534cm)
# echa side is made of 2 pieces of 62" and 37"
# we have 13 pieces: 7 at 62" x 14" x 210" and 6 at 37" x 14" x 210", while Cornell 3 for their use
# see drawings at
# https://solid.jlab.org/files/cleo_manual/20120314%20Accelerator%20SCANS/Roll%202%20-%20Laboratory%20of%20Nuclear%20Studies/6052-303%20-Sh%2003.pdf
# https://solid.jlab.org/files/cleo_manual/20120314%20Accelerator%20SCANS/Roll%202%20-%20Laboratory%20of%20Nuclear%20Studies/6052-303%20-Sh%2006.pdf
#  https://solid.jlab.org/files/cleo_manual/20120314%20Accelerator%20SCANS/Roll%202%20-%20Laboratory%20of%20Nuclear%20Studies/
# 
#  photos at
#  https://solid.jlab.org/files/cleo_manual/cleo_scan1.pdf
#  https://solid.jlab.org/files/cleo_manual/cleo_scan2.pdf

sub solid_DDVCS_muon_forwardangle
{
make_solid_DDVCS_muon_forwardangle_1();
make_solid_DDVCS_muon_forwardangle_2();
make_solid_DDVCS_muon_forwardangle_3();
}

 my $color="00ffff";
 my $material="G4_Fe";

sub make_solid_DDVCS_muon_forwardangle_1
{
    my $z=570+50+18;    
    my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_11";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "160*cm 0*cm $z*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = $color;
    $detector{"type"}       = "Box";
    $detector{"dimensions"}  = "125*cm 267*cm 18*cm";
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
    
    $detector{"name"}        = "$DetectorName\_12";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "-160*cm 0*cm $z*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = $color;
    $detector{"type"}       = "Box";
    $detector{"dimensions"}  = "125*cm 267*cm 18*cm";
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


sub make_solid_DDVCS_muon_forwardangle_2
{
    my $z=570+50+18+18+10+18;    
    my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_21";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 160*cm $z*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = $color;
    $detector{"type"}       = "Box";
    $detector{"dimensions"}  = "267*cm 125*cm 18*cm";
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
    
    $detector{"name"}        = "$DetectorName\_22";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm -160*cm $z*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = $color;
    $detector{"type"}       = "Box";
    $detector{"dimensions"}  = "267*cm 125*cm 18*cm";
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


sub make_solid_DDVCS_muon_forwardangle_3
{
    my $z=570+50+18+18+10+18+18+10+18;    
    my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_31";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "160*cm 0*cm $z*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = $color;
    $detector{"type"}       = "Box";
    $detector{"dimensions"}  = "125*cm 267*cm 18*cm";
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
    
    $detector{"name"}        = "$DetectorName\_32";
    $detector{"mother"}      = "$DetectorMother";
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "-160*cm 0*cm $z*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = $color;
    $detector{"type"}       = "Box";
    $detector{"dimensions"}  = "125*cm 267*cm 18*cm";
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
