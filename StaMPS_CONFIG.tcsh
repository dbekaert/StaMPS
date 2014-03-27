setenv STAMPS       "/home/ahooper/StaMPS_v3.2"
setenv SAR          "/home/ahooper/software/ROI_PAC_3_0"
setenv GETORB_BIN   "/home/ahooper/software/getorb/bin"
setenv SAR_ODR_DIR  "/home/ahooper/software/SAR_FILES/ODR"
#setenv SAR_PRC_DIR  "/home/ahooper/software/SAR_FILES/PRC"
setenv VOR_DIR      "/home/ahooper/software/SAR_FILES/VOR"
setenv INS_DIR      "/home/ahooper/software/SAR_FILES/INS"
setenv DORIS_BIN    "/home/ahooper/software/doris_v4.02/bin"
setenv TRIANGLE_BIN "/home/ahooper/software/triangle/bin"
setenv SNAPHU_BIN   "/home/ahooper/software/snaphu-v1.4.2/bin"

setenv ROI_PAC      "$SAR/ROI_PAC"
#####################################
# ROI_PAC VERSION 3 
#####################################
setenv INT_BIN      "$ROI_PAC/INT_BIN"
setenv INT_SCR      "$ROI_PAC/INT_SCR"
#####################################

#####################################
# ROI_PAC VERSION 2.3 and before 
#####################################
#set MACH=`uname -s`
#if ($MACH == "HP-UX") then
#  setenv ARCHC      HP
#else if ($MACH == "IRIX") then
#  setenv ARCHC      SGI
#else if ($MACH == "SunOS") then
#  setenv ARCHC      SUN
#else if ($MACH == "Linux") then
#  setenv ARCHC      LIN
#else if ($MACH == "Darwin") then
#  setenv ARCHC      MAC
#endif
#setenv INT_LIB      "$ROI_PAC/LIB/$ARCHC"
#setenv INT_BIN      "$ROI_PAC/BIN/$ARCHC"
#setenv FFTW_LIB     "$SAR/FFTW/$ARCHC""_fftw_lib"
#####################################

#####################################
# shouldn't need to change below here
#####################################

setenv MY_BIN       "$INT_BIN"
setenv MATLABPATH   $STAMPS/matlab:`echo $MATLABPATH`
setenv DORIS_SCR    "$STAMPS/DORIS_SCR"

# Needed for ROI_PAC (a bit different to standard)

### use points not commas for decimals, and give dates in US english
setenv LC_NUMERIC "en_US.UTF-8"
setenv LC_TIME "en_US.UTF-8"


setenv MY_SAR       "$SAR"
setenv OUR_SCR      "$MY_SAR/OUR_SCR"
setenv MY_SCR       "$STAMPS/ROI_PAC_SCR"

setenv SAR_TAPE     "/dev/rmt/0mn"

set path = ( $path $STAMPS/bin $MY_SCR $INT_BIN $INT_SCR $OUR_SCR $DORIS_SCR $GETORB_BIN $DORIS_BIN $TRIANGLE_BIN $SNAPHU_BIN )


