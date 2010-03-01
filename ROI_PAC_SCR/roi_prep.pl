#!/usr/bin/perl
### roi_prep.pl
###########################################################
# Modified Version for StaMPS. 2007/12/01 Andy Hooper
#
# 07/2009 AH Azimuth resolution generalized for all sensors
###########################################################
#
#
# $Id: roi_prep.pl,v 1.2 2007/08/22 18:22:32 ericf Exp $
#
# $Log: roi_prep.pl,v $
# Revision 1.2  2007/08/22 18:22:32  ericf
# changed default patch size to 16K for L-band data
#
# Revision 1.1.1.1  2007/04/17 22:28:01  ericmg
# initial import into erda insarcvs
#
# Revision 1.1.1.1  2007/04/03 03:40:10  ericmg
# initial import into erda insarcvs
#
# Revision 1.5  2007/01/26 01:26:30  ericf
# fixed slc_ymax sometimes not an integer, courtesy Paul Lundgren
#
# Revision 1.4  2006/12/06 21:17:45  ericf
# changed variable names and equations to clarify calculations
# modified calculation of the last usable data in SLC (SLC_YMAX)
#
# Revision 1.3  2005/12/16 18:45:39  bswift
# RSC keyword documentation changes.
#

$] >= 5.004 or die "Perl version must be >= 5.004 (Currently $]).\n";

use Env qw(INT_SCR INT_BIN);
use lib "$INT_SCR";  #### Location of Generic.pm
use Generic;

###Usage info/check

sub Usage{

`$INT_SCR/pod2man.pl  $INT_SCR/roi_prep.pl`;
exit 1;
}
@ARGV ==2 or Usage();
@args = @ARGV;

$prefix       = shift;
$orbit_type   = shift;

### Check for input files:
if    (-r "$prefix.proc"){$processfile = "$prefix.proc";}
elsif (-r "../roi.proc") {$processfile = "../roi.proc";}
else  {$processfile = "";}

if ($processfile){
##########################
  Message "Read the inputfile";
##########################
  open IN, "$processfile" or die "Can't read $processfile\n";
  while (chomp($line = <IN>)){
    $line =~ /=/ or next;
    $line =~ s/\s//g; ###Remove all whitespace
    ($keyname, $value) = split /=/, $line;
    $$keyname = $value;
  }
  close(IN);
}

# changed these to compute below EJF Jan. 2005
# --set to zero if not in .proc files so calc is done below
$before_z_ext       or $before_z_ext       = 0;
$after_z_ext        or $after_z_ext        = 0;
$near_rng_ext       or $near_rng_ext       = 0;
$far_rng_ext        or $far_rng_ext        = 0;
$valid_samples      or $valid_samples      = 0;

# moved to calculate below EJF 2007/8/21 $patch_size or $patch_size = 4096;

$number_of_patches  or $number_of_patches  = 0;
###Can also set number_of_patches and out_pixel

### For StaMPS AJH
$mean_pixel_rng     or $mean_pixel_rng     = 3156;
################## 

$az_frac = 0.8; # default amount of azimuth extension (fraction of synthetic aperture)
$rg_frac = 0.5; # default amount of range extension (fraction of chirp)
$az_round = 0.55; # fraction of last patch to add for rounding--decreased from 0.8 EJF 06/11/28

### For StaMPS AJH
#$dopfile    = "roi.dop";
$dopfile    = "$prefix.raw";
################## 

$infile     = "$prefix.raw";
$outfile    = "$prefix.slc";
$roiin      = "$prefix"."_roi.in";
$roiout     = "$prefix"."_roi.out";

##########################
Message "Checking I/O";
##########################
@Infiles    = ($infile, "$infile.rsc", "$dopfile.rsc");
@Outfiles   = ($roiin, "$outfile.rsc");
IOcheck(\@Infiles, \@Outfiles);
Log("roi_prep.pl", @args);


###########################
Message "setting roi default";
###########################
$speed_of_light = 299792458;

##########################################
Message "Reading resource file: $infile.rsc";
##########################################

$orbit_number    = Use_rsc "$infile read ORBIT_NUMBER";
$width           = Use_rsc "$infile read WIDTH";
$xmin            = Use_rsc "$infile read XMIN";
$xmax            = Use_rsc "$infile read XMAX";
### For StaMPS AJH
$ymin or $ymin            = Use_rsc "$infile read YMIN";
$ymax or $ymax            = Use_rsc "$infile read YMAX";
##################
$velocity        = Use_rsc "$infile read VELOCITY";
$height          = Use_rsc "$infile read HEIGHT";
$earth_radius    = Use_rsc "$infile read EARTH_RADIUS";
$starting_rng    = Use_rsc "$infile read STARTING_RANGE";
$prf             = Use_rsc "$infile read PRF";
$wavelength      = Use_rsc "$infile read WAVELENGTH";
$pulse_length    = Use_rsc "$infile read PULSE_LENGTH";
$chirp_slope     = Use_rsc "$infile read CHIRP_SLOPE";
$sampling_freq   = Use_rsc "$infile read RANGE_SAMPLING_FREQUENCY";
$rng_offset      = Use_rsc "$infile read RANGE_OFFSET";
$i_bias          = Use_rsc "$infile read I_BIAS";
$q_bias          = Use_rsc "$infile read Q_BIAS";
$first_line_utc  = Use_rsc "$infile read FIRST_LINE_UTC";
$orbit_direction = Use_rsc "$infile read ORBIT_DIRECTION";
$PLANET_GM       = Use_rsc "$infile read PLANET_GM";

$AntennaSide     = Use_rsc "$infile read ANTENNA_SIDE"; 
if(    $AntennaSide == -1 ){ $AntennaSide = "Right"; }
elsif( $AntennaSide == 1 ){  $AntennaSide = "Left"; }
else{                        $AntennaSide = "Unknown"; }

$AntennaLen      = Use_rsc "$infile read ANTENNA_LENGTH";

$dop_rng0        = Use_rsc "$dopfile read DOPPLER_RANGE0";
$dop_rng1        = Use_rsc "$dopfile read DOPPLER_RANGE1";
$dop_rng2        = Use_rsc "$dopfile read DOPPLER_RANGE2";
$dop_rng3        = Use_rsc "$dopfile read DOPPLER_RANGE3";

### For StaMPS AJH
#$sl_azim_res     = Use_rsc "$dopfile read SL_AZIMUT_RESOL";
$sl_azim_res  = $AntennaLen/2;
##################
#
$squint		 = Use_rsc "$dopfile read SQUINT";

if($use1dopp){
$dop_rng0        = Use_rsc "$infile read DOPPLER_RANGE0";
$dop_rng1        = Use_rsc "$infile read DOPPLER_RANGE1";
$dop_rng2        = Use_rsc "$infile read DOPPLER_RANGE2";
$dop_rng3        = Use_rsc "$infile read DOPPLER_RANGE3";
$squint		 = Use_rsc "$infile read SQUINT";
}    

### For StaMPS AJH
$dop_rng0=$dop_rng0+$dop_rng1*$mean_pixel_rng;
$dop_rng0=$dop_rng0+$dop_rng2*$mean_pixel_rng*$mean_pixel_rng;
$dop_rng1=0;
$dop_rng2=0;
$dop_rng3=0;
##################

############################
Message "Computing parameters";
############################
$chirp_samps = int($sampling_freq * $pulse_length + 0.5);

Message "chirp length in samples $chirp_samps";
# set default range extensions to 1/2 chirp
$near_rng_ext or $near_rng_ext = int($chirp_samps * $rg_frac );
$far_rng_ext or $far_rng_ext = int($chirp_samps * $rg_frac );  # this is now the real far range extension

$rng_pixel_size          = $speed_of_light/$sampling_freq/2;
$out_pixel or $out_pixel = int($xmax/2-$xmin/2-$chirp_samps+$near_rng_ext+$far_rng_ext);
$slc_starting_rng        = $starting_rng - $rng_pixel_size*$near_rng_ext;

# use far range
$slc_far_range = $starting_rng + $out_pixel*$rng_pixel_size;

# use specified if present or change default based on wavelength EJF 2007/8/20
$patch_size or $patch_size = 1;

if ($patch_size == 1) {  # ugly Perl, but could not add the if to the "or" above
  if ( $wavelength < 0.07 )
    {$patch_size = 4096;} # e.g., C-band
  else
    {$patch_size = 16384;}   # need much longer FFT for L-band
}

Message "Azimuth patch length $patch_size";

$synth_apert_samps = int($wavelength*$slc_far_range*$prf/($AntennaLen*$velocity) +0.5);

Message "Synthetic aperture length in samples $synth_apert_samps";
# default azimuth extensions
$before_z_ext or $before_z_ext = int($az_frac * $synth_apert_samps + 0.5);
$after_z_ext or $after_z_ext = int($az_frac * $synth_apert_samps + 0.5);

# valid azimuth samples per patch
$valid_samples or $valid_samples = $patch_size - $synth_apert_samps;
# need to make sure this is an even number
$valid_samples = (2* int($valid_samples/2));

# processing region relative to raw data
$proc_start   = int($ymin-$before_z_ext +0.5);
$proc_end     = int($ymax+$after_z_ext +0.5);

# slc start line is offset due to FFT wrap-around trimming
$slc_first_line_offset = ($patch_size-$valid_samples)/2+$proc_start;

$delta_line_utc          = 1/$prf;
$azimuth_pixel_size      = $velocity/$prf;
# data trimmed on both ends by ($patch_size-$valid_samples)/2, 
#   add $az_round to make sure we get end when taking int()
$number_of_patches  
or $number_of_patches    = int(($proc_end-$proc_start-($patch_size-$valid_samples))/$valid_samples + $az_round );

$slc_first_line_utc      = $first_line_utc + $slc_first_line_offset/$prf;
$slc_length              = $number_of_patches*$valid_samples;
$slc_last_line_utc           = $slc_first_line_utc + $slc_length/$prf;
$slc_center_line_utc         = ($slc_first_line_utc + $slc_last_line_utc)/2;

# last remotely useful SLC is $az_frac * $synth_apert_samps past last raw data
# calc SLC end relative to raw data
$slc_end = $slc_length-1 - $slc_first_line_offset;

if ( $slc_end < $ymax + $az_frac*$synth_apert_samps ) 
  { $slc_ymax = $slc_length-1; }  # all SLC good data
else
 { $slc_ymax                = $ymax + int($az_frac*$synth_apert_samps); }  # part of SLC not good

$dop0                    = $dop_rng0;
$dop1                    = $dop_rng1;
$dop2                    = $dop_rng2;
$dop3                    = $dop_rng3;

##########################################
Message "Writing resource file: $outfile.rsc";
##########################################
`cp $infile.rsc $outfile.rsc`;

Use_rsc "$outfile write FILE_LENGTH              $slc_length";
Use_rsc "$outfile write ORBIT_NUMBER             $orbit_number"; 
Use_rsc "$outfile write WIDTH                    $out_pixel";
Use_rsc "$outfile write FILE_START               1";   
Doc_rsc(
 RSC_Tip => 'Starting line (record) number',
 RSC_Doc => q[
   Starting line (record) number of file.
   Is combined with YMIN to determine actual starting line number.
   ],
 RSC_Derivation => q[
   Set to constant '1' in roi_prep.pl.
   ],
 RSC_Comment => q[
   Only used as 'start' parameter of trees (in new_cut.pl),
   corr_flag and grass (in unwrap.pl).

   In source code first record read (offs) is calculated as follows -
     offs=start+ymin-1;    /* first line of file to start reading/writing */
   ],
 RSC_Type => Int,
 RSC_Unit => 'record'
);


Use_rsc "$outfile write XMIN                     0";  
Use_rsc "$outfile write XMAX                     $out_pixel";  
Use_rsc "$outfile write YMIN                     0";   
Use_rsc "$outfile write YMAX                     $slc_ymax";   
Use_rsc "$outfile write DOPPLER_RANGE0           $dop_rng0";  
Use_rsc "$outfile write DOPPLER_RANGE1           $dop_rng1";    
Use_rsc "$outfile write DOPPLER_RANGE2           $dop_rng2";   
Use_rsc "$outfile write DOPPLER_RANGE3           $dop_rng3";   
Use_rsc "$outfile write EARTH_RADIUS             $earth_radius";
Use_rsc "$outfile write STARTING_RANGE           $slc_starting_rng"; 
Use_rsc "$outfile write RAW_DATA_RANGE           $starting_rng"; 
Doc_rsc(
 RSC_Tip => 'Range of first sample in raw data file',
 RSC_Doc => q[
   See make_raw.pl STARTING_RANGE.
   ],
 RSC_Derivation => q[
   Value is copy of .raw STARTING_RANGE keyword.
   ],
 RSC_Type => Real,
 RSC_Unit => 'SI:meter',
);


Use_rsc "$outfile write PRF                      $prf";  
Use_rsc "$outfile write WAVELENGTH               $wavelength";   
Use_rsc "$outfile write PULSE_LENGTH             $pulse_length";  
Use_rsc "$outfile write CHIRP_SLOPE              $chirp_slope";
Use_rsc "$outfile write RANGE_SAMPLING_FREQUENCY $sampling_freq"; 
Use_rsc "$outfile write RANGE_PIXEL_SIZE         $rng_pixel_size";
Use_rsc "$outfile write AZIMUTH_PIXEL_SIZE       $azimuth_pixel_size";  
Doc_rsc(
 RSC_Derivation => q[
   # $velocity and $prf are read from raw resource file in roi_prep.pl
   $azimuth_pixel_size      = $velocity/$prf;
   ],
 RSC_Type => Real,
 RSC_Unit => 'SI:meter',
);


Use_rsc "$outfile write DELTA_LINE_UTC           $delta_line_utc";
Doc_rsc(
 RSC_Tip => 'Time between lines (records)',
 RSC_Derivation => q[
   $delta_line_utc          = 1/$prf;
   ],
 RSC_Comment => q[
   Initial value is IPP.
   ],
 RSC_Type => Real,
 RSC_Unit => 'SI:second/record',
);


Use_rsc "$outfile write SQUINT                   $squint";
Use_rsc "$outfile write RANGE_OFFSET             $rng_offset";
Doc_rsc(
 RSC_Derivation => q[
   Keyword RANGE_OFFSET doesn't exist in *.raw.rsc,
   so a value of 0 is returned from Use_rsc.
   ],
 RSC_Comment => q[
   Does not appear to be used anywhere.
   ],
);


Use_rsc "$outfile write FIRST_LINE_UTC           $slc_first_line_utc";
Use_rsc "$outfile write CENTER_LINE_UTC          $slc_center_line_utc";
Use_rsc "$outfile write LAST_LINE_UTC            $slc_last_line_utc";
Use_rsc "$outfile write RLOOKS                   1";
Doc_rsc(
 RSC_Tip => 'Range Looks',
 RSC_Doc => q[
   'Number of Range Looks' based on RDF usage in resamp.pl.
   ],
 RSC_Derivation => q[
   Set to constant '1' in roi_prep.pl.
   Updated in diffnsim.pl
     $new_rlks    = $rlks*$nlook;
   Updated in look.pl
     $newRlooks          = $oldRlooks*$Rlooks;
   ],
 RSC_Comment => q[
   Does anybody ever take fractional looks?
   ],
 RSC_Type => Int,
 RSC_Unit => 'element/element',
);


Use_rsc "$outfile write ALOOKS                   1";
Doc_rsc(
 RSC_Tip => 'Azimuth Looks',
 RSC_Doc => q[
   'Number of Azimuth Looks' based on RDF usage in resamp.pl.
   ],
 RSC_Derivation => q[
   Set to constant '1' in roi_prep.pl.
   Updated in diffnsim.pl
     $new_alks    = $alks*$nlook;
   Updated in look.pl
     $newAlooks          = $oldAlooks*$Alooks;
   ],
 RSC_Comment => q[
   Does anybody ever take fractional looks?
   ],
 RSC_Type => Int,
 RSC_Unit => 'record/record',
);



###############################################
# Get peg info and update geometry parameters #
# now extracted for SLC, not raw data         #
###############################################

system "cp ${infile}.rsc debug.rsc";

system "$INT_SCR/GetPeg.pl $outfile $orbit_type";

($name = $outfile) =~ s/\.[^.]*$//; #strip the last extension from the input file name

open PEGOUT, "$name.peg.out";
while( defined( $Line = <PEGOUT> ) ){
	if( $Line =~ /Peg Lat\/Lon , H =\s+(\S+)\s+(\S+)\s+(\S+)/ ){
		$Latd   = $1;
		$Lond   = $2;
		$PegHgt = $3;
	}
	if( $Line =~ /Peg Heading =\s+(\S+)/ ){
		$hdgd = $1;
	}
	if( $Line =~ /Vertical Fit:\s+(\S+)\s+(\S+)\s+(\S+)/ ){
		@Height_poly = ( $1, $2, $3 );
	}
	if( $Line =~ /Horizontal Fit:\s+(\S+)\s+(\S+)\s+(\S+)/ ){
		@CrossT_poly = ( $1, $2, $3 );
	}
	if( $Line =~ /Vertical Velocity Fit:\s+(\S+)\s+(\S+)/ ){
		@Vert_V_poly = ( $1, $2 );
	}
	if( $Line =~ /Cross-Track Velocity Fit:\s+(\S+)\s+(\S+)/ ){
		@CrossT_V_poly = ( $1, $2 );
	}
	if( $Line =~ /Along-Track Velocity Fit:\s+(\S+)\s+(\S+)/ ){
		@AlongT_V_poly = ( $1, $2 );
	}
	if( $Line =~ /Platform SCH Velocity \(m\/s\):\s+(\S+)\s+(\S+)\s+(\S+)/ ){
		@VelocitySCH = ( $1, $2, $3 );
		$velocity_mid = Norm( @VelocitySCH );
	}
	if( $Line =~ /Platform SCH Acceleration \(m\/s\^2\):\s+(\S+)\s+(\S+)\s+(\S+)/ ){
		@AccelerationSCH = ( $1, $2, $3 );
	}
	if( $Line =~ /Time to first\/middle scene:\s+\S+\s+(\S+)/ ){
		$PegUtc = $1;
	}
}

close PEGOUT;

$HgtDt = $Height_poly[1] * $VelocitySCH[0];

Use_rsc "$outfile write HEIGHT                   $Height_poly[0]";
Doc_rsc(
 RSC_Tip => 'Platform Altitude at Peg Point',
 RSC_Doc => q[
   "Platform Altitude" based on RDF usage in diffnsim.pl

   First coefficient (constant term) of orbit "Vertical Fit" polynomial
   output by get_peg_info run by GetPeg.pl.

   Polynomial is function of SCH 'S' coordinate.
   Value is in 'H' direction, which is height above SCH reference sphere.

   SCH Coordinate system.
   ],
 RSC_Derivation => q[
   See baseline/get_peg_info.f
   ],
 RSC_Type => Real,
 RSC_Unit => 'SI:meter',
);

Use_rsc "$outfile write HEIGHT_DS                $Height_poly[1]";
Doc_rsc(
 RSC_Tip => 'Platform Altitude Rate at Peg Point',
 RSC_Doc => q[
   "Platform Altitude Rate" based on RDF usage in diffnsim.pl

   Second coefficient (linear term) of orbit "Vertical Fit" polynomial
   output by get_peg_info run by GetPeg.pl.

   Polynomial is function of SCH 'S' coordinate.
   Value is in 'H' direction, which is height above SCH reference sphere.

   SCH Coordinate system.
   ],
 RSC_Derivation => q[
   See baseline/get_peg_info.f
   ],
 RSC_Type => Real,
 RSC_Unit => 'SI:meter/SI:meter',
);

Use_rsc "$outfile write HEIGHT_DDS               $Height_poly[2]";
Doc_rsc(
 RSC_Tip => 'Platform Altitude Acceleration at Peg Point',
 RSC_Doc => q[
   "Platform Altitude Acceleration" based on RDF usage in diffnsim.pl

   Third coefficient (quadratic term) of orbit "Vertical Fit" polynomial
   output by get_peg_info run by GetPeg.pl.

   Polynomial is function of SCH 'S' coordinate.
   Value is in 'H' direction, which is height above SCH reference sphere.

   SCH Coordinate system.
   ],
 RSC_Derivation => q[
   See baseline/get_peg_info.f
   ],
 RSC_Type => Real,
 RSC_Unit => 'SI:meter/SI:meter**2',
);

Use_rsc "$outfile write HEIGHT_DT                $HgtDt";
Doc_rsc(
 RSC_Tip => 'Platform Altitude change w.r.t. time at Peg Point',
 RSC_Derivation => q[
   $HgtDt = $Height_poly[1] * $VelocitySCH[0];
   ],
 RSC_Comment => q[
   Does not appear to ever be used.
   ],
 RSC_Type => Real,
 RSC_Unit => 'SI:second',
);

Use_rsc "$outfile write CROSSTRACK_POS           $CrossT_poly[0]";
Use_rsc "$outfile write CROSSTRACK_POS_DS        $CrossT_poly[1]";
Use_rsc "$outfile write CROSSTRACK_POS_DDS       $CrossT_poly[2]";
Use_rsc "$outfile write VELOCITY                 $velocity_mid";
Doc_rsc(
 RSC_Tip => 'Norm of Platform SCH Velocity at Peg Point',
 RSC_Doc => q[
   "Body fixed S/C velocities" based on RDF usage in roi_prep.pl and autofocus.pl.
   "Spacecraft Along Track Velocity" based on RDF usage in inverse3d.pl
   "Platform Velocity" based on RDF usage in diffnsim.pl and phase2base.pl
   ],
 RSC_Derivation => q[
   $velocity_mid = Norm( @VelocitySCH );
   ],
 RSC_Type => Real,
 RSC_Unit => 'SI:meter/SI:second',
);

Use_rsc "$outfile write VELOCITY_S               $VelocitySCH[0]";
Doc_rsc(
 RSC_Tip => 'Platform Velocity S Component',
 RSC_Doc => q[
   'S' Component of 'Platform SCH Velocity'
   produced by get_peg_info run by GetPeg.pl.

   SCH Coordinate system.
   ],
 RSC_Derivation => q[
   See baseline/get_peg_info.f
   ],
 RSC_Comment => q[
   Appears to only be used by roi_prep.pl and autofocus.pl
   ],
 RSC_Type => Real,
 RSC_Unit => 'SI:meter/SI:second',
);

Use_rsc "$outfile write VELOCITY_C               $VelocitySCH[1]";
Doc_rsc(
 RSC_Tip => 'Platform Velocity C Component',
 RSC_Doc => q[
   'C' Component of 'Platform SCH Velocity'
   produced by get_peg_info run by GetPeg.pl.

   SCH Coordinate system.
   ],
 RSC_Derivation => q[
   See baseline/get_peg_info.f
   ],
 RSC_Comment => q[
   Appears to only be used by roi_prep.pl and autofocus.pl
   Note - this is not the speed of light.
   ],
 RSC_Type => Real,
 RSC_Unit => 'SI:meter/SI:second',
);
Use_rsc "$outfile write VELOCITY_H               $VelocitySCH[2]";
Doc_rsc(
 RSC_Tip => 'Platform Velocity H Component',
 RSC_Doc => q[
   'H' Component of 'Platform SCH Velocity'
   produced by get_peg_info run by GetPeg.pl.

   SCH Coordinate system.
   ],
 RSC_Derivation => q[
   See baseline/get_peg_info.f
   ],
 RSC_Comment => q[
   Appears to only be used by roi_prep.pl and autofocus.pl
   ],
 RSC_Type => Real,
 RSC_Unit => 'SI:meter/SI:second',
);

Use_rsc "$outfile write ACCELERATION_S           $AccelerationSCH[0]";
Doc_rsc(
 RSC_Tip => 'Platform Acceleration S Component',
 RSC_Doc => q[
   'S' Component of 'Platform SCH Acceleration'
   produced by get_peg_info run by GetPeg.pl.

   SCH Coordinate system.
   ],
 RSC_Derivation => q[
   See baseline/get_peg_info.f
   ],
 RSC_Comment => q[
   Appears to only be used by roi_prep.pl and autofocus.pl
   ],
 RSC_Type => Real,
 RSC_Unit => 'SI:meter/SI:second**2',
);

Use_rsc "$outfile write ACCELERATION_C           $AccelerationSCH[1]";
Doc_rsc(
 RSC_Tip => 'Platform Acceleration C Component',
 RSC_Doc => q[
   'C' Component of 'Platform SCH Acceleration'
   produced by get_peg_info run by GetPeg.pl.

   SCH Coordinate system.
   ],
 RSC_Derivation => q[
   See baseline/get_peg_info.f
   ],
 RSC_Comment => q[
   Appears to only be used by roi_prep.pl and autofocus.pl
   ],
 RSC_Type => Real,
 RSC_Unit => 'SI:meter/SI:second**2',
);

Use_rsc "$outfile write ACCELERATION_H           $AccelerationSCH[2]";
Doc_rsc(
 RSC_Tip => 'Platform Acceleration H Component',
 RSC_Doc => q[
   'H' Component of 'Platform SCH Acceleration'
   produced by get_peg_info run by GetPeg.pl.

   SCH Coordinate system.
   ],
 RSC_Derivation => q[
   See baseline/get_peg_info.f
   ],
 RSC_Comment => q[
   Appears to only be used by roi_prep.pl and autofocus.pl
   ],
 RSC_Type => Real,
 RSC_Unit => 'SI:meter/SI:second**2',
);

Use_rsc "$outfile write VERT_VELOCITY            $Vert_V_poly[0]";
Use_rsc "$outfile write VERT_VELOCITY_DS         $Vert_V_poly[1]";
Use_rsc "$outfile write CROSSTRACK_VELOCITY      $CrossT_V_poly[0]";
Use_rsc "$outfile write CROSSTRACK_VELOCITY_DS   $CrossT_V_poly[1]";
Use_rsc "$outfile write ALONGTRACK_VELOCITY      $AlongT_V_poly[0]";
Use_rsc "$outfile write ALONGTRACK_VELOCITY_DS   $AlongT_V_poly[1]";
Use_rsc "$outfile write LATITUDE     $Latd";
Doc_rsc(
 RSC_Tip => 'Latitude of SCH Peg Point',
 RSC_Doc => q[
   'Lat' value of 'Peg Lat/Lon , H'
   produced by get_peg_info run by GetPeg.pl.

   ],
 RSC_Derivation => q[
   See baseline/get_peg_info.f
   ],
 RSC_Comment => q[
   Appears to only be used in inverse3d.pl to specify "Peg Point Data"
   ],
 RSC_Type => Real,
 RSC_Unit => 'SI:degree',
);

Use_rsc "$outfile write LONGITUDE    $Lond";
Doc_rsc(
 RSC_Tip => 'Longitude of SCH Peg Point',
 RSC_Doc => q[
   'Lon' value of 'Peg Lat/Lon , H'
   produced by get_peg_info run by GetPeg.pl.

   ],
 RSC_Derivation => q[
   See baseline/get_peg_info.f
   ],
 RSC_Comment => q[
   Appears to only be used in inverse3d.pl to specify "Peg Point Data"
   ],
 RSC_Type => Real,
 RSC_Unit => 'SI:degree',
);

Use_rsc "$outfile write HEADING      $hdgd";
Doc_rsc(
 RSC_Tip => 'Heading of SCH Peg Point',
 RSC_Doc => q[
   'Peg Heading' value
   produced by get_peg_info run by GetPeg.pl.
   ],
 RSC_Derivation => q[
   See baseline/get_peg_info.f
   ],
 RSC_Comment => q[
   Appears to only be used in inverse3d.pl to specify "Peg Point Data"
   ],
 RSC_Type => Real,
 RSC_Unit => 'SI:degree',
);

Use_rsc "$outfile write PEG_UTC      $PegUtc";
Doc_rsc(
 RSC_Tip => 'Scene start time-of-day',
 RSC_Doc => q[
   'first' value of 'Time to first/middle scene'
   produced by get_peg_info run by GetPeg.pl.
   ],
 RSC_Derivation => q[
   See baseline/get_peg_info.f
   ],
 RSC_Comment => q[
   Appears to only be used in inverse3d.pl to calculate
     $PegLine = int( ( $PegUtc - $slc_first_line_utc ) / $delta_line_utc );
   which is used for
     'Reference Line for SCH Coordinates' RDF value
   ],
 RSC_Type => Real,
 RSC_Unit => 'SI:second',
);

#########################################
Message "Writing roi input file: ${prefix}_roi.in";
#########################################

if ($near_rng_ext < 0){ 
  $first_range_bin = 1-$near_rng_ext;
  $near_rng_ext=0;
}
else{
  $first_range_bin = 1;
}

open (ROI,">$roiin") or die "Can't write to $roiin\n";

$xmin = int(($xmin/2)+0.5);

print ROI <<END;

RDF File for ROI           
Master input file                                (-) = $infile                         ! Raw data for channel 1
Second input file                                (-) = /dev/null  		       ! Raw data for channel 2
Master SLC                                       (-) = $outfile                        ! SLC channel 1
Second SLC                                       (-) = /dev/null                       ! SLC channel 2
Eight look image file                            (-) = 8lk                             ! Eight look image file
Debug flag                                       (-) = 0                               ! 
Bytes per input line                             (-) = $width    $width                ! files 1 and 2
Good bytes per input line                        (-) = $xmax $xmax                     ! including header, files 1 and 2
First line to read                               (-) = $proc_start                           ! 
Number of patches                                (-) = $number_of_patches              ! org 40, 1 for testing
First sample pair to use                         (-) = $xmin                           ! 
Number of valid pulses                           (-) = $valid_samples                  ! 
Azimuth Patch Size                               (-) = $patch_size                     !
### For StaMPS AJH
Deskew?                                          (-) = y                               ! 
##################

Caltone location                                 (-) = 0 0                             !  % of sample rate, files 1 and 2
Start range bin, number to process               (-) = $first_range_bin $out_pixel     ! org 6626, test with 1600
Delta azimuth, range pixels for second file      (-) = 0 0                    ! 
Master file Doppler centroid coefs               (-) = $dop0 $dop1 $dop2 $dop3         ! (Hz/prf) quadratic polynomial
Second file Doppler centroid coefs               (-) = $dop0 $dop1 $dop2 $dop3         ! (Hz/prf) quadratic polynomial
Doppler average method                           (-) = 1                               ! 1: use file 1 dopp, 2: file 2, 3: avg
Earth radius                                     (m) = $earth_radius                   ! 
Body fixed S/C velocities                      (m/s) = $velocity $velocity             ! files 1, 2
Spacecraft height                                (m) = $height $height                 ! files 1, 2
Planet GM					 (-) = $PLANET_GM
Left, Right or Unknown Pointing 		 (-) = $AntennaSide
SCH Velocity Vector 1			       (m/s,m/s,m/s) = $VelocitySCH[0] $VelocitySCH[1] $VelocitySCH[2]
SCH Velocity Vector 2			       (m/s,m/s,m/s) = $VelocitySCH[0] $VelocitySCH[1] $VelocitySCH[2]
SCH Acceleration Vector 1		       (-,-,-) = $AccelerationSCH[0] $AccelerationSCH[1] $AccelerationSCH[2] 
SCH Acceleration Vector 2		       (-,-,-) = $AccelerationSCH[0] $AccelerationSCH[1] $AccelerationSCH[2]



Range of first sample in data file               (m) = $starting_rng $starting_rng     ! files 1, 2
PRF                                              (-) = $prf $prf                       ! pps, files 1, 2
i/q means                                        (-) = $i_bias $q_bias $i_bias $q_bias ! file 1 (i,q), file 2 (i,q)
Flip i/q?                                        (-) = n                               ! 
Azimuth resolution                               (m) = $sl_azim_res
                               !  
Number of azimuth looks                          (-) = 4                               ! 
Range sampling rate                             (Hz) = $sampling_freq                  ! 
Range chirp slope                                (-) = $chirp_slope                    ! Hz/s
Range pulse duration                             (s) = $pulse_length                   ! 
Range chirp extension points                     (-) = $near_rng_ext                   ! 
Secondary range migration correction?            (-) = n                               ! 
Radar wavelength                                 (m) = $wavelength                     ! 
Range spectral weighting                         (-) = 1.                              ! 
Spectral shift fraction                          (-) = 0. 0.                           ! % of range/azimuth bandwidth to remove
Linear resampling coefs                          (-) = 0 0 0 0                         ! sloper,  intr,  slopea,  inta 
Linear resampling deltas                         (-) = 0 0 0 0                         ! dsloper, dintr, dslopea, dinta
AGC file                                         (-) = $infile.agcf                    ! AGC file name
DWP file                                         (-) = $infile.dwpf                    ! DWP file name
END
close(ROI);
###############################################################################
###############################################################################

exit 0;

=pod


=head1 USAGE

B<roi_prep.pl> I<date> OrbitType

date: radical of the input file I<date>.raw

Looks for I<date>.proc, then roi.proc in the next directory up, then uses default values

 defaults:
 before_z_ext = 0.8 * synthetic aperture
 after_z_ext  = 0.8 * synthetic aperture
 near_rng_ext = 1/2 chirp samples
 far_rng_ext  = 1/2 chirp samples
 patch_size   = 4096 (16384 for L-band)

Can also set number_of_patches, valid_samples (lines saved per patch) and out_pixel

=head1 FUNCTION

Creates input files for ROI

=head1 ROUTINES CALLED

none

=head1 CALLED BY

process.pl

=head1 FILES USED

I<date>.raw.rsc

roi.dop (generated by dopav.pl)

I<date>.proc [optional]

roi.proc [optional]


=head1 FILES CREATED

I<date>.slc.rsc

I<date>_roi.in

I<date>_roi.out

=head1 HISTORY

Shell Script : Francois ROGEZ 96/98

Perl  Script : Rowena LOHMAN 04/18/98

Frederic Crampe, Jan 10, 2000

Eric Fielding, Jan. 23, 2005

Eric Fielding, Nov. 29, 2006

=head1 LAST UPDATE

Eric Fielding, Aug. 22, 2007

=cut
