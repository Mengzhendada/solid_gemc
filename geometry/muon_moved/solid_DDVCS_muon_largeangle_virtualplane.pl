#!/usr/bin/perl -w
use strict;
use warnings;
our %detector;
our %configuration;
our %parameters;

use Getopt::Long;
use Math::Trig;

my $DetectorName = 'solid_DDVCS_muon_largeangle_virtualplane';

my $DetectorMother="root";

sub solid_DDVCS_muon_largeangle_virtualplane
{
make_solid_DDVCS_muon_largeangle_virtualplane_endcapdonut_1();
# make_solid_DDVCS_muon_largeangle_virtualplane_endcapdonut_2();
# make_solid_DDVCS_muon_largeangle_virtualplane_barrel_1();
make_solid_DDVCS_muon_largeangle_virtualplane_barrel_2();
}

sub make_solid_DDVCS_muon_largeangle_virtualplane_endcapdonut_1
{
 my $z=(209+570+50)/2;
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_endcapdonut_1";
 $detector{"mother"}      = "root";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm $z*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "CC6633";
 $detector{"type"}        = "Cons";
  my $Rmin1 = 286;
  my $Rmax1 = 286.1;
  my $Rmin2 = 286;
  my $Rmax2 = 286.1;
  my $Dz    = (570+50-209)/2;
  my $Sphi  = 0;
  my $Dphi  = 360;
 $detector{"dimensions"}  = "$Rmin1*cm $Rmax1*cm $Rmin2*cm $Rmax2*cm $Dz*cm $Sphi*deg $Dphi*deg";
 $detector{"material"}    = "G4_Galactic";
 $detector{"mfield"}      = "no";
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 1;
 $detector{"style"}       = 0;
 $detector{"sensitivity"} = "flux";
 $detector{"hit_type"}    = "flux";
 $detector{"identifiers"} = "id manual 6210000";
 print_det(\%configuration, \%detector);
 }
 
sub make_solid_DDVCS_muon_largeangle_virtualplane_endcapdonut_2
{
 my $z=(209+570+50)/2;
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_endcapdonut_2";
 $detector{"mother"}      = "root";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm $z*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "CC6633";
 $detector{"type"}        = "Cons";
  my $Rmin1 = 327;
  my $Rmax1 = 327.1;
  my $Rmin2 = 327;
  my $Rmax2 = 327.1;
  my $Dz    = (570+50-209)/2;
  my $Sphi  = 0;
  my $Dphi  = 360;
 $detector{"dimensions"}  = "$Rmin1*cm $Rmax1*cm $Rmin2*cm $Rmax2*cm $Dz*cm $Sphi*deg $Dphi*deg";
 $detector{"material"}    = "G4_Galactic";
 $detector{"mfield"}      = "no";
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 1;
 $detector{"style"}       = 0;
 $detector{"sensitivity"} = "flux";
 $detector{"hit_type"}    = "flux";
 $detector{"identifiers"} = "id manual 6220000";
 print_det(\%configuration, \%detector);
 }
 
sub make_solid_DDVCS_muon_largeangle_virtualplane_barrel_1
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_barrel_1";
 $detector{"mother"}      = "root";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm 0*cm";
 $detector{"rotation"}   = "0*deg 0*deg 22.5*deg";
 $detector{"color"}       = "CC6633";
 $detector{"type"}       = "Pgon";
# match solenoid_v4 
 $detector{"dimensions"} = "0*deg 360*deg 8*counts 2*counts 235*cm 235*cm 235.1*cm 235.1*cm -241*cm 182*cm";
 $detector{"material"}    = "G4_Galactic";
 $detector{"mfield"}      = "no";
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 1;
 $detector{"style"}       = 0;
 $detector{"sensitivity"} = "flux";
 $detector{"hit_type"}    = "flux";
 $detector{"identifiers"} = "id manual 6310000";
 print_det(\%configuration, \%detector);
}


sub make_solid_DDVCS_muon_largeangle_virtualplane_barrel_2
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_barrel_2";
 $detector{"mother"}      = "root";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm 0*cm";
 $detector{"rotation"}   = "0*deg 0*deg 22.5*deg";
 $detector{"color"}       = "CC6633"; 
 $detector{"type"}       = "Pgon";
# match solenoid_v4
 $detector{"dimensions"} = "0*deg 360*deg 8*counts 2*counts 280*cm 280*cm 280.1*cm 280.1*cm -266*cm 182*cm";
 $detector{"material"}    = "G4_Galactic";
 $detector{"mfield"}      = "no";
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 1;
 $detector{"style"}       = 0;
 $detector{"sensitivity"} = "flux";
 $detector{"hit_type"}    = "flux";
 $detector{"identifiers"} = "id manual 6320000";
 print_det(\%configuration, \%detector);
}
