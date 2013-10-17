#!/usr/bin/perl -w

use strict;

use lib ("$ENV{GEMC}/database_io");
use geo;
use geo qw($pi);

my $DetectorName = 'CLEO_SIDIS_ec_forwardangle_kill';

my $envelope = "solid_$DetectorName";
my $file     = "solid_$DetectorName.txt";

my $rmin      = 1;
my $rmax      = 1000000;

my %detector = ();    # hash (map) that defines the gemc detector
$detector{"rmin"} = $rmin;
$detector{"rmax"} = $rmax;

use Getopt::Long;
use Math::Trig;

my $DetectorMother="root";

# HOD2SLATS       2
# HOD2MEDIUM    654
# HOD2MOTHER   'HALL' 
# HOD2IDTYPE       42
# HOD2GATE       80.
# HOD2THRES       0
# HOD2SHAP      'CONE' 'CONE'
# HOD2SIZE1      9.0  108.  197. 111. 202.  
# # HOD2SIZE2      2.0  107.00  196. 108. 197.  
# HOD2TYPE        1 2
# HOD2POSX        2*0.
# HOD2POSY        2*0.
# HOD2POSZ        740. 729. 

 my $color="0000ff";
 my $material="Kryptonite";
 my $material_shield="Kryptonite";

sub make_ec_forwardangle
{
#  my $z=775-350;
#   my $z=765-350;
#   my $z=430;
  my $z=440;

    $detector{"name"}        = "$DetectorName";
    $detector{"mother"}      = "$DetectorMother" ;
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $z*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = $color;
    $detector{"type"}       = "Cons";
    my $Rmin1 = 98;
    my $Rmax1 = 230;
    my $Rmin2 = 98;
    my $Rmax2 = 230;
    my $Dz    = 25;
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
    $detector{"identifiers"} = "";
    print_det(\%detector, $file);
}
make_ec_forwardangle();

sub make_ec_forwardangle_shield
{
 $detector{"name"}        = "$DetectorName\_shield";
 $detector{"mother"}      = "root";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "0*cm 0*cm 413.8*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = $color;
 $detector{"type"}        = "Cons";
  my $Rmin1 = 98;
  my $Rmax1 = 230;
  my $Rmin2 = 98;
  my $Rmax2 = 230;
  my $Dz    = 0.56;
  my $Sphi  = 0;
  my $Dphi  = 360;
 $detector{"dimensions"}  = "$Rmin1*cm $Rmax1*cm $Rmin2*cm $Rmax2*cm $Dz*cm $Sphi*deg $Dphi*deg";
 $detector{"material"}    = $material_shield;
 $detector{"mfield"}      = "no";
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 1;
 $detector{"style"}       = 0;
 $detector{"sensitivity"} = "no";
 $detector{"hit_type"}    = "";
 $detector{"identifiers"} = "";
 print_det(\%detector, $file);
}
make_ec_forwardangle_shield();
