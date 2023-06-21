#!/bin/csh -f

source /jlab/2.5/ce/jlab.csh
setenv SoLID_GEMC /solid_gemc
setenv LD_LIBRARY_PATH ${GEMC}:${LD_LIBRARY_PATH}
setenv PATH ${SoLID_GEMC}/source/${GEMC_VERSION}:${PATH}

#run simulation
solid_gemc solid.gcard -FILTER_HITS=1 -BEAM_P="e-,5.3*GeV,18*deg,180*deg" -SPREAD_P="5.3*GeV,0*deg,0*deg" -BEAM_V="(0,0,0)cm" -SPREAD_V="(0,0)cm" -N=1e4 | cat >> out.log
#convert evio to root
#need more option for banks other than flux
# which  evio2root
evio2root -INPUTF=out.evio -B="/work/halla/solid/tianye/HallCbeamtest/solid_gemc/prototype/ec_pro/HallC2022/geometry/cc_pro   /work/halla/solid/tianye/HallCbeamtest/solid_gemc/prototype/ec_pro/HallC2022/geometry/solid_ec  /work/halla/solid/tianye/HallCbeamtest/solid_gemc/prototype/ec_pro/HallC2022/geometry/solid_gem /work/halla/solid/tianye/HallCbeamtest/solid_gemc/prototype/ec_pro/HallC2022/geometry/solid_spd " -R="flux" | cat >> out.log

# cp /solid_gemc/analysis/SoLIDFileReduce.C ./
# root -b -q 'SoLIDFileReduce.C+("SIDIS","out.root","2")'

# mv -f out.root_skim out.root


