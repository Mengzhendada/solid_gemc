use strict;
use warnings;
our %detector;
our %configuration;
our %parameters;

use Getopt::Long;
use Math::Trig;

my $DetectorName = 'solid_SIDIS_target_He3_collimator';

my $DetectorMother="root";

my $material_collimator = "SL_target_He3_TungstenPowder";

sub solid_SIDIS_target_He3_collimator
{
make_collimator_upstream();
make_collimator_downstream();
}


sub make_collimator_upstream
{
    my $z_collimator_upstream	= -355;

    my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_upstream";
    $detector{"mother"}      = "$DetectorMother" ;
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm $z_collimator_upstream*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "003300";
    $detector{"type"}       = "Cons";
    #cover 7.5-15deg
#     my $Rmin1 = 1.3;
#     my $Rmax1 = 2.7;
#     my $Rmin2 = 2.6;
#     my $Rmax2 = 5.4;
    #cover 7.0-15deg
    my $Rmin1 = 1.23;	#block target upstream window z=-370 at 7deg
#     my $Rmax1 = 4.00;	#block beamline upstream window z=-375 at 15deg
    my $Rmax1 = 6.35;	#block beamline upstream window z=-375 at 23deg
    my $Rmin2 = 2.46;	#block target upstream window z=-370 at 7deg
#     my $Rmax2 = 6.70;	#block beamline upstream window z=-375 at 15deg
    my $Rmax2 = 10.60;	#block beamline upstream window z=-375 at 23deg
    my $Dz    = 5;
    $detector{"dimensions"}  = "$Rmin1*cm $Rmax1*cm $Rmin2*cm $Rmax2*cm $Dz*cm 0*deg 360*deg";
    $detector{"material"}   = "$material_collimator";
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

# sub make_collimator_downstream
# {
#     my $z_collimator_downstream	= -315;
#     
#     my %detector=init_det();
#     $detector{"name"}        = "$DetectorName\_downstream";
#     $detector{"mother"}      = "$DetectorMother" ;
#     $detector{"description"} = $detector{"name"};
#     $detector{"pos"}        = "0*cm 0*cm $z_collimator_downstream*cm";
#     $detector{"rotation"}   = "0*deg 0*deg 0*deg";
#     $detector{"color"}      = "003300";
#     $detector{"type"}       = "Cons";
#     #cover 7.5-15deg
# #     my $Rmin1 = 1.3;
# #     my $Rmax1 = 2.7;
# #     my $Rmin2 = 2.6;
# #     my $Rmax2 = 5.4;
# #     my $Dz    = 5;
#     #cover 7.0-15deg
#     my $Rmin1 = 1.23;	#block target downstream window z=-330 at 7deg
#     my $Rmax1 = 2.70;	#block target downstream window z=-330 at 15deg
#     my $Rmin2 = 1.85;	#block beamline downstream window z=-325 at 7deg
#     my $Rmax2 = 5.40;	#block target downstream window z=-330 at 15deg
#     my $Dz    = 5;
#     my $Rmin1 = 1.23;	#block target downstream window z=-330 at 7deg
#     $detector{"dimensions"}  = "$Rmin1*cm $Rmax1*cm $Rmin2*cm $Rmax2*cm $Dz*cm 0*deg 360*deg";
#     $detector{"material"}   = "$material_collimator";
#     $detector{"mfield"}     = "no";
#     $detector{"ncopy"}      = 1;
#     $detector{"pMany"}       = 1;
#     $detector{"exist"}       = 1;
#     $detector{"visible"}     = 1;
#     $detector{"style"}       = 1;
#     $detector{"sensitivity"} = "no";
#     $detector{"hit_type"}    = "no";
#     $detector{"identifiers"} = "no";
#     print_det(\%configuration, \%detector);
# }

sub make_collimator_downstream
{ 
    my %detector=init_det();
    $detector{"name"}        = "$DetectorName\_downstream";
    $detector{"mother"}      = "$DetectorMother" ;
    $detector{"description"} = $detector{"name"};
    $detector{"pos"}        = "0*cm 0*cm 0*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "003300";
    $detector{"type"}        = "Polycone";
    #cover 7.0-15deg
    my $Rmin1 = 1.23;	#block target downstream window z=-330 at 7deg
    my $Rmax1 = 2.70;	#block target downstream window z=-330 at 15deg
    my $Rmin2 = 1.71;	#block beamline downstream window z=-325 at 6deg
    my $Rmax2 = 5.40;	#block target downstream window z=-330 at 15deg
    my $Rmin3 = 2.85;	#block beamline downstream window z=-325 at 6deg
    my $Rmax3 = 5.40;	#flat out
    my $Z1 = -320;
    my $Z2 = -310;    
    my $Z3 = -300;    
#     my $Rmin1 = 1.23;	#block target downstream window z=-330 at 7deg
#     my $Rmin1 = 1.13;	#block target downstream window z=-330 at 6.5deg    
#     my $Rmax1 = 2.70;	#block target downstream window z=-330 at 15deg
#     my $Rmin2 = 2.85;	#block beamline downstream window z=-325 at 6.5deg
#     my $Rmax2 = 8.04;	#block target downstream window z=-330 at 15deg
#     my $Dz    = 10;
    $detector{"dimensions"}  = "0*deg 360*deg 3*counts $Rmin1*cm $Rmin2*cm $Rmin3*cm $Rmax1*cm $Rmax2*cm $Rmax3*cm $Z1*cm $Z2*cm $Z3*cm";
    $detector{"material"}   = "$material_collimator";
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

