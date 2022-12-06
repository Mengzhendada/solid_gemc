
my $cleo_pvdis_fw_zpos = 194.0 + $cleo_fw_hthk;
# my $cleo_pvdis_fw_zpos = 189.2 + $cleo_fw_hthk;
#my $cleo_bw_zpos = 301.0 - $cleo_bw_hthk;
my $cleo_bw_zpos = 312.0 - $cleo_bw_hthk;
my $cleo_sidis_fw_zpos = 97.0 + $cleo_fw_hthk;
my $bisec_hthk = 0.5;


sub make_mother
{
    $detector{"name"}        = $detMom;
    $detector{"mother"}      = "root";
    $detector{"description"} = "SoLID Hall container";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"color"}       = "ff0000";
    $detector{"type"}        = "Tube";
    $detector{"dimensions"}  = "0.0*cm 500.0*cm 800*cm 0*deg 360*deg";
    $detector{"material"}    = "Vacuum";
#    $detector{"material"}    = "Air_Opt";
    $detector{"mfield"}      = "no";
    $detector{"ncopy"}       = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 0;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "no";
    $detector{"hit_type"}    = "no";
    $detector{"identifiers"} = "no";

    print_det(\%configuration, \%detector);

}1;

sub make_tank_BaBar_pvdis
{
    $detector{"name"} = $namePre."Tank_1";
    $detector{"mother"} = $detMom ;
    $detector{"description"} = "Tank part 1";
    $detector{"pos"} = "0.0*cm 0.0*cm 0.0*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "FF9900";
    $detector{"type"}       = "Polycone";
    $detector{"dimensions"} = "0.0*deg 360.0*deg 2*counts 65.25*cm 70.65*cm 150.54*cm 169.57*cm 200.00*cm 225.0*cm";
    $detector{"material"}   = "Component";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "no";
    $detector{"hit_type"}    = "no";
    $detector{"identifiers"} = "Tank_1";

    print_det(\%configuration, \%detector);

    $detector{"name"} = $namePre."Tank_2";
    $detector{"mother"} = $detMom ;
    $detector{"description"} = "Tank part 2";
    $detector{"pos"} = "0.0*cm 0.0*cm 0.0*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "FF9900";
    $detector{"type"}       = "Polycone";
    $detector{"dimensions"} = "0.0*deg 360.0*deg 4*counts  240.0*cm  190.0*cm  70.65*cm  87.05*cm 257.0*cm  265.0*cm  265.0*cm  265.0*cm 209.0*cm  225.0*cm  225.0*cm  301.0*cm";
    $detector{"material"}   = "Component";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "no";
    $detector{"hit_type"}    = "no";
    $detector{"identifiers"} = "Tank_2";

    print_det(\%configuration, \%detector);

    $detector{"name"} = $namePre."Tank";
    $detector{"mother"} = $detMom ;
    $detector{"description"} = "Tank total";
    $detector{"pos"} = "0.0*cm 0.0*cm 0.0*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "FF9900";
    $detector{"type"}       = "Operation:@ ".$namePre."Tank_1 + ".$namePre."Tank_2";
    #$detector{"material"}   = "PVDIS_C4F8O_N2_Gas";
    $detector{"material"}   = "SL_LGCCgas_SIDIS";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "no";
    $detector{"hit_type"}    = "no";
    $detector{"identifiers"} = "Tank";

    print_det(\%configuration, \%detector);
}1;

sub make_tank_CLEO_pvdis
{
    $detector{"name"}        = $namePre."Tank";
    $detector{"mother"}      = $detMom ;
    $detector{"description"} = "tank";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"color"}       = "FF9900";
    $detector{"type"}        = "Polycone";
	#older settings:
  #  $detector{"dimensions"}  = "0*deg 360*deg 4*counts 71*cm 73*cm 73*cm 85*cm 144*cm 155*cm 265*cm 265*cm 194*cm 209.01*cm 209.01*cm 301*cm";
 #   $detector{"dimensions"}  = "0*deg 180*deg 4*counts 65*cm 67*cm 85*cm 85*cm 144*cm 155*cm 265*cm 265*cm 194*cm 209.01*cm 209.01*cm 301*cm";
 #   $detector{"dimensions"}  = "0*deg 180*deg 4*counts 71*cm 73*cm 73*cm 85*cm 144*cm 155*cm 265*cm 265*cm 194*cm 209.01*cm 209.01*cm 301*cm"; #for cross-section of tank

    #$detector{"dimensions"}  = "0*deg 360*deg 4*counts 65*cm 67*cm 67*cm 85*cm 144*cm 155*cm 265*cm 265*cm 194*cm 209.01*cm 209.01*cm 301*cm";
	
# 	$detector{"dimensions"}  = "0*deg 360*deg 5*counts  66.5*cm 69.5*cm 71.5*cm 71.5*cm 85*cm 127.0*cm 150*cm 163.0*cm 263*cm 263*cm 189.2*cm 209.0*cm 225.0*cm 225.0*cm 312*cm";
	#$detector{"dimensions"}  = "0*deg 360*deg 4*counts  69.5*cm 71.5*cm 71.5*cm 85*cm 150*cm 163.0*cm 263*cm 263*cm 209.0*cm 225.0*cm 225.0*cm 312*cm";
#    $detector{"dimensions"}  = "0*deg 180*deg 4*counts 65*cm 67*cm 67*cm 85*cm 144*cm 155*cm 265*cm 265*cm 194*cm 209.01*cm 209.01*cm 301*cm";

#current, standard nose
$detector{"dimensions"}  = "0*deg 360*deg 5*counts 63.2*cm 69.3*cm 71*cm 71*cm 83*cm 135*cm 150*cm 163*cm 263*cm 263*cm 194*cm 210*cm 226*cm 226*cm 312*cm";
if($use_modNose){
	#current, alt-nose (7 deg):
	#$detector{"dimensions"}  = "0*deg 360*deg 6*counts 52.43*cm 64.34*cm 66.3*cm 68.27*cm 68.27*cm 78.83*cm 126*cm 135*cm 150*cm 163*cm 263*cm 263*cm 97*cm 194*cm 210*cm 226*cm 226*cm 312*cm";
	#current, alt-nose (6.5 deg):
	#$detector{"dimensions"}  = "0*deg 360*deg 6*counts 48.65*cm 59.70*cm 61.5*cm 63.35*cm 63.35*cm 73.14*cm 126*cm 135*cm 150*cm 163*cm 263*cm 263*cm 97*cm 194*cm 210*cm 226*cm 226*cm 312*cm";
	#including cable space:
	#$detector{"dimensions"}  = "0*deg 360*deg 5*counts 59.70*cm 61.5*cm 63.35*cm 63.35*cm 73.14*cm 130*cm 145*cm 163*cm 263*cm 263*cm 194*cm 210*cm 232*cm 226*cm 312*cm";
	#add a 0.1 mm to inner diameter:
	$detector{"dimensions"}  = "0*deg 360*deg 5*counts 59.75*cm 61.55*cm 63.4*cm 63.4*cm 73.2*cm 130*cm 145*cm 163*cm 263*cm 263*cm 194*cm 210*cm 226*cm 226*cm 312*cm";
}	


#    if($use_pvdis){
#	$detector{"material"}    = "PVDIS_C4F8O_N2_Gas";
	$detector{"material"}   = "SL_LGCCgas_SIDIS";
#    }else{
#	$detector{"material"}    = "SIDIS_CO2_gas";
#	$detector{"material"}    = "Vacuum";
#	$detector{"material"}   = "SL_LGCCgas_SIDIS";
#    }
    $detector{"mfield"}      = "no";
    $detector{"ncopy"}       = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
#    $detector{"visible"}     = 0;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "no";
    $detector{"hit_type"}    = "no";
    $detector{"identifiers"} = "no";
    print_det(\%configuration, \%detector);

#back window:
    $detector{"name"}        = $namePre."Tank_window_back";
    $detector{"mother"}      = $namePre."Tank";
    $detector{"description"} = "tank window back";
    $detector{"pos"}         = "0*cm 0*cm ".$cleo_bw_zpos."*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"color"}       = "FF9900";
    $detector{"type"}        = "Tube";
    $detector{"dimensions"}  = "85.1*cm 264.9*cm ".$cleo_bw_hthk."*cm 0*deg 360*deg";
 #   $detector{"material"}    = "Aluminum";
    $detector{"material"}    = "SL_PVF";
    $detector{"mfield"}      = "no";
    $detector{"ncopy"}       = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "no";
    $detector{"hit_type"}    = "no";
    $detector{"identifiers"} = "no";
    print_det(\%configuration, \%detector);

#bisection

    $detector{"name"}        = $namePre."Tank_bisection1";
    $detector{"mother"}      = $detMom ;
    $detector{"description"} = "tank bisection1";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"color"}       = "FF9900";
    $detector{"type"}        = "Polycone";
    $detector{"dimensions"}  = "0*deg 360*deg 4*counts 71*cm 73*cm 73*cm 85*cm 144*cm 155*cm 265*cm 265*cm 194*cm 209.01*cm 209.01*cm 301*cm";
    $detector{"mfield"}      = "no";
    $detector{"material"}   = "Component";
    $detector{"ncopy"}       = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 0;
 #   $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "no";
    $detector{"hit_type"}    = "no";
    $detector{"identifiers"} = "no";
#    print_det(\%configuration, \%detector);

    $detector{"name"}        = $namePre."Tank_bisection2";
    $detector{"mother"}      = $detMom ;
    $detector{"description"} = "tank bisection1";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"color"}       = "FF9900";
    $detector{"type"}        = "Box";
    $detector{"dimensions"}  = "301*cm ".$bisec_hthk."*cm 301*cm";
    $detector{"mfield"}      = "no";
    $detector{"material"}   = "Component";
    $detector{"ncopy"}       = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 0;
 #   $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "no";
    $detector{"hit_type"}    = "no";
    $detector{"identifiers"} = "no";
 #   print_det(\%configuration, \%detector);

    $detector{"name"} = $namePre."Tank_bisection";
    $detector{"mother"} = $detMom ;
    $detector{"description"} = "Tank Bisection";
    $detector{"pos"} = "0.0*cm 0.0*cm 0.0*cm";
    $detector{"rotation"}   = "0*deg 0*deg 0*deg";
    $detector{"color"}      = "FF9900";
    $detector{"type"}       = "Operation:@ ".$namePre."Tank_bisection1 * ".$namePre."Tank_bisection2";
    #$detector{"material"}   = "PVDIS_C4F8O_N2_Gas";
    $detector{"material"}   = "SL_LGCCgas_SIDIS";
    $detector{"mfield"}     = "no";
    $detector{"ncopy"}      = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "no";
    $detector{"hit_type"}    = "no";
    $detector{"identifiers"} = "Tank";

#    print_det(\%configuration, \%detector);



}1;

sub make_tank_CLEO_frontwindow
{
#front window:
    $detector{"name"}        = $namePre."Tank_window_front";
    $detector{"mother"}      = $namePre."Tank";
    $detector{"description"} = "tank window front";
    $detector{"pos"}         = "0*cm 0*cm ".$cleo_pvdis_fw_zpos."*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"color"}       = "FF9900";
    $detector{"type"}        = "Tube";
#    $detector{"dimensions"}  = "65.1*cm 143.9*cm ".$cleo_fw_hthk."*cm 0*deg 360*deg";
 	$detector{"dimensions"}  = "66.5*cm 127.0*cm ".$cleo_fw_hthk."*cm 0*deg 360*deg";
    $detector{"material"}    = "SL_PVF";
    $detector{"mfield"}      = "no";
    $detector{"ncopy"}       = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 0;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "no";
    $detector{"hit_type"}    = "no";
    $detector{"identifiers"} = "no";
    print_det(\%configuration, \%detector);

#For LEM detector plane:

    $detector{"name"}        = $namePre."Tank_LEM_window";
    $detector{"mother"}      = $namePre."Tank";
    $detector{"description"} = "tank window detector plane";
    $detector{"pos"}         = "0*cm 0*cm ".($cleo_pvdis_fw_zpos + 0.5)."*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"color"}       = "FF9900";
    $detector{"type"}        = "Tube";
    $detector{"dimensions"}  = "65.1*cm 143.9*cm 0.001*mm 0*deg 360*deg";
    if($use_pvdis){
	#$detector{"material"}    = "PVDIS_C4F8O_N2_Gas";
	$detector{"material"}   = "SL_LGCCgas_SIDIS";
    }else{
	#$detector{"material"}    = "SIDIS_CO2_gas";
	$detector{"material"}   = "SL_LGCCgas_SIDIS";
    }
    $detector{"mfield"}      = "no";
    $detector{"ncopy"}       = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 0;
    $detector{"style"}       = 1;
#    $detector{"sensitivity"} = "FLUX";
#    $detector{"hit_type"}    = "FLUX";
    $detector{"sensitivity"} = "no";
    $detector{"hit_type"}    = "no";
    $detector{"identifiers"} = "no";
    print_det(\%configuration, \%detector);


}1;

sub make_tank_CLEO_sidis
{
	print "!!!!!!!!!!!!   Making world and tank!!!!\n";
	# this is just the snout.  should be removed once everything is working
#    $detector{"name"}        = $namePre."Tank2";
#    $detector{"mother"}      = $detMom ;
#    $detector{"description"} = "tank";
#    $detector{"pos"}         = "0*cm 0*cm 0*cm";
#    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
#    $detector{"color"}       = "FF9900";
#    $detector{"type"}        = "Polycone";
#    #$detector{"dimensions"}  = "0*deg 360*deg 2 58*cm 71*cm 127*cm 144*cm 97*cm 194*cm";
#    $detector{"dimensions"}  = "0*deg 360*deg 2*counts 58*cm 65*cm 127*cm 144*cm 97*cm 194*cm";
#    #$detector{"material"}    = "SIDIS_CO2_gas";
#    $detector{"material"}   = "SL_LGCCgas_SIDIS";
#    $detector{"mfield"}      = "no";
#    $detector{"ncopy"}       = 1;
#    $detector{"pMany"}       = 1;
#    $detector{"exist"}       = 1;
#    $detector{"visible"}     = 1;
#    $detector{"style"}       = 1;
#    $detector{"sensitivity"} = "no";
#    $detector{"hit_type"}    = "no";
#    $detector{"identifiers"} = "no";
#    print_det(\%configuration, \%detector);

    $detector{"name"}        = $namePre."Tank";
    $detector{"mother"}      = $detMom ;
    $detector{"description"} = "tank";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"color"}       = "FF9900";
    $detector{"type"}        = "Polycone";
	#older settings:
#    $detector{"dimensions"}  = "0*deg 360*deg 5*counts 58*cm 65*cm 67*cm 67*cm 85*cm 127*cm 144*cm 155*cm 265*cm 265*cm 97*cm 194*cm 209.01*cm 209.01*cm 301*cm";
#  $detector{"dimensions"}  = "0*deg 360*deg 6*counts  58*cm 66.5*cm 69.5*cm 70.5*cm 70.5*cm 86*cm 127*cm 128.5*cm 153*cm 164.0*cm 264*cm 264*cm 97*cm 189.2*cm 210.0*cm 225.0*cm 225.0*cm 312*cm";	
#     $detector{"dimensions"}  = "0*deg 360*deg 6*counts  53*cm 66.5*cm 69.5*cm 71.5*cm 71.5*cm 85*cm 127*cm 127.0*cm 150*cm 163.0*cm 263*cm 263*cm 97*cm 189.2*cm 209.0*cm 225.0*cm 225.0*cm 312*cm";

#current, standard nose:

	$detector{"dimensions"}  = "0*deg 360*deg 6*counts 50*cm 63.2*cm 69.3*cm 71*cm 71*cm 83*cm 126*cm 135*cm 150*cm 163*cm 263*cm 263*cm 97*cm 194*cm 210*cm 226*cm 226*cm 312*cm";    
	
	
	if($use_modNose){
		#current, alt-nose (7 deg):
		#$detector{"dimensions"}  = "0*deg 360*deg 6*counts 52.43*cm 64.34*cm 66.3*cm 68.27*cm 68.27*cm 78.83*cm 126*cm 135*cm 150*cm 163*cm 263*cm 263*cm 97*cm 194*cm 210*cm 226*cm 226*cm 312*cm";
		#current, alt-nose (6.5 deg):
		#$detector{"dimensions"}  = "0*deg 360*deg 6*counts 48.65*cm 59.70*cm 61.5*cm 63.35*cm 63.35*cm 73.14*cm 126*cm 135*cm 150*cm 163*cm 263*cm 263*cm 97*cm 194*cm 210*cm 226*cm 226*cm 312*cm";
		#including cable space:
		#include a mm to inner radius:
		#$detector{"dimensions"}  = "0*deg 360*deg 6*counts 48.65*cm 59.70*cm 61.5*cm 63.35*cm 63.35*cm 73.14*cm 130*cm 130*cm 145*cm 163*cm 263*cm 263*cm 97*cm 194*cm 210*cm 232*cm 226*cm 312*cm";
		$detector{"dimensions"}  = "0*deg 360*deg 6*counts 48.7*cm 59.75*cm 61.55*cm 63.4*cm 63.4*cm 73.2*cm 130*cm 130*cm 145*cm 163*cm 263*cm 263*cm 97*cm 194*cm 210*cm 232*cm 232*cm 312*cm";
	}	


	$detector{"material"}   = "SL_LGCCgas_SIDIS";
#	$detector{"material"}    = "SIDIS_CO2_gas";
    $detector{"mfield"}      = "no";
    $detector{"ncopy"}       = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "no";
    $detector{"hit_type"}    = "no";
    $detector{"identifiers"} = "no";
    print_det(\%configuration, \%detector);



#front window:
    $detector{"name"}        = $namePre."Tank_window_front";
    $detector{"mother"}      = $namePre."Tank";
    $detector{"description"} = "tank window front";
    $detector{"pos"}         = "0*cm 0*cm ".$cleo_sidis_fw_zpos."*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"color"}       = "FF9900";
    $detector{"type"}        = "Tube";
    $detector{"dimensions"}  = "58.1*cm 126.9*cm ".$cleo_fw_hthk."*cm 0*deg 360*deg";
    $detector{"material"}    = "SL_PVF";
    $detector{"mfield"}      = "no";
    $detector{"ncopy"}       = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "no";
    $detector{"hit_type"}    = "no";
    $detector{"identifiers"} = "no";
    print_det(\%configuration, \%detector);

}1;

sub make_snout
{
    $detector{"name"}        = $namePre."_Snout";
    $detector{"mother"}      = $detMom ;
    $detector{"description"} = "tank";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"color"}       = "FF9900";
    $detector{"type"}        = "Polycone";
    #$detector{"dimensions"}  = "0*deg 360*deg 4 71*cm 73*cm 73*cm 85*cm 144*cm 155*cm 265*cm 265*cm 194*cm 209.01*cm 209.01*cm 301*cm";
    $detector{"dimensions"}  = "0*deg 360*deg 4*counts 65*cm 67*cm 85*cm 85*cm 144*cm 155*cm 265*cm 265*cm 194*cm 209.01*cm 209.01*cm 301*cm";
#    $detector{"dimensions"}  = "0*deg 180*deg 4*counts 71*cm 73*cm 73*cm 85*cm 144*cm 155*cm 265*cm 265*cm 194*cm 209.01*cm 209.01*cm 301*cm"; #for cross-section of tank
    $detector{"material"}    = "Lead";
    $detector{"mfield"}      = "no";
    $detector{"ncopy"}       = 1;
    $detector{"pMany"}       = 1;
    $detector{"exist"}       = 1;
#     $detector{"visible"}     = 0;
   $detector{"visible"}     = 1;
    $detector{"style"}       = 1;
    $detector{"sensitivity"} = "no";
    $detector{"hit_type"}    = "no";
    $detector{"identifiers"} = "no";
    print_det(\%configuration, \%detector);
}1;
