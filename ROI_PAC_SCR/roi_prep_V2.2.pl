#!/usr/bin/perl
### roi_prep.pl

$] >= 5.004 or die "Perl version must be >= 5.004 (Currently $]).\n";

use Env qw(INT_SCR INT_BIN);
use lib "$INT_SCR";  #### Location of Generic.pm
use Generic;

###Usage info/check

sub Usage{

`$INT_SCR/pod2man.pl  $INT_SCR/roi_prep.pl`;
exit 1;
}
@ARGV ==1 or Usage();
@
rgs = @ARGV;

$prefix       = shift;

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
    Message "$keyname = $value";
  }
  close(IN);
}

$before_z_ext       or $before_z_ext       = 1500;
$after_z_ext        or $after_z_ext        = 1500;
$near_rng_ext       or $near_rng_ext       = 700;
$far_rng_ext        or $far_rng_ext        = 0;
$valid_samples      or $valid_samples      = 3072;
$patch_size         or $patch_size         = 4096;
$number_of_patches  or $number_of_patches  = 0;
$mean_pixel_rng     or $mean_pixel_rng     = 3156;
###Can also set number_of_patches

$infile     = "$prefix.raw";
$outfile    = "$prefix.slc";
$roiin      = "$prefix"."_roi.in";
$roiout     = "$prefix"."_roi.out";

##########################
Message "Checking I/O";
##########################
@Infiles    = ($infile, "$infile.rsc");
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
$VelocitySCH[0]  = Use_rsc "$infile read VELOCITY_S";
$VelocitySCH[1]  = Use_rsc "$infile read VELOCITY_C";
$VelocitySCH[2]  = Use_rsc "$infile read VELOCITY_H";
$AccelerationSCH[0]  = Use_rsc "$infile read ACCELERATION_S";
$AccelerationSCH[1]  = Use_rsc "$infile read ACCELERATION_C";
$AccelerationSCH[2]  = Use_rsc "$infile read ACCELERATION_H";

$AntennaSide     = Use_rsc "$infile read ANTENNA_SIDE"; 
if(    $AntennaSide == -1 ){ $AntennaSide = "Right"; }
elsif( $AntennaSide == 1 ){  $AntennaSide = "Left"; }
else{                        $AntennaSide = "Unknown"; }

$dop_rng0        = Use_rsc "$infile read DOPPLER_RANGE0";
$dop_rng1        = Use_rsc "$infile read DOPPLER_RANGE1";
$dop_rng2        = Use_rsc "$infile read DOPPLER_RANGE2";
$dop_rng3        = Use_rsc "$infile read DOPPLER_RANGE3";
$squint		 = Use_rsc "$infile read SQUINT";


$dop_rng0=$dop_rng0+$dop_rng1*$mean_pixel_rng;
$dop_rng0=$dop_rng0+$dop_rng2*$mean_pixel_rng*$mean_pixel_rng;
$dop_rng1=0;
$dop_rng2=0;
$dop_rng3=0;

############################
Message "Computing parameters";
############################
$ymin     = $ymin-$before_z_ext;
$ymax     = $ymax+$after_z_ext;


#if ($orbit_direction eq "ascending" ){
#  $squint   = 0.222; 
#  $dop_rng0 = -0.15;
#  $dop_rng1 = 0; 
#  $dop_rng2 = 0; 
#}
#else squint=0.222 ; dop_rng0=0.25  ; dop_rng1=0 ; dop_rng2=0 
#else{
#  $squint   = "0.1";
#  $dop_rng0 = "0.25";
#  $dop_rng1 = "0";
#  $dop_rng2 = "0";
#} 

$delta_line_utc          = 1/$prf;
$azimuth_pixel_size      = $velocity/$prf;
$rng_pixel_size          = $speed_of_light/$sampling_freq/2;
$out_pixel or $out_pixel = int($xmax/2-$xmin/2+$near_rng_ext+$far_rng_ext);
$number_of_patches  
or $number_of_patches    = int(($ymax-$ymin-($patch_size-$valid_samples))/$valid_samples+0.99);
$first_line_utc          = $first_line_utc + (($patch_size-$valid_samples)/2+$ymin)/$prf;
$slc_starting_rng        = $starting_rng - $rng_pixel_size*$near_rng_ext;
$slc_length              = $number_of_patches*$valid_samples;
$last_line_utc           = $first_line_utc + $slc_length/$prf;
$center_line_utc         = ($first_line_utc + $last_line_utc)/2;
$slc_ymax                = $slc_length -1;
#$dop0                   = $dop_rng0 + ($dop_rng1 + $dop_rng2*$starting_rng)*$starting_rng;
#$dop1                   = ($dop_rng1+2*$dop_rng2*$slc_starting_rng*$rng_pixel_size)*$rng_pixel_size;
#$dop2                   = $dop_rng2 * ($rng_pixel_size*$rng_pixel_size);
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
Use_rsc "$outfile write XMIN                     0";  
Use_rsc "$outfile write XMAX                     $out_pixel";  
Use_rsc "$outfile write YMIN                     0";   
Use_rsc "$outfile write YMAX                     $slc_ymax";   
Use_rsc "$outfile write DOPPLER_RANGE0           $dop_rng0";  
Use_rsc "$outfile write DOPPLER_RANGE1           $dop_rng1";    
Use_rsc "$outfile write DOPPLER_RANGE2           $dop_rng2";   
Use_rsc "$outfile write DOPPLER_RANGE3           $dop_rng3";   
Use_rsc "$outfile write VELOCITY                 $velocity"; 
Use_rsc "$outfile write HEIGHT                   $height";  
Use_rsc "$outfile write EARTH_RADIUS             $earth_radius";
Use_rsc "$outfile write STARTING_RANGE           $slc_starting_rng"; 
Use_rsc "$outfile write RAW_DATA_RANGE           $starting_rng"; 
Use_rsc "$outfile write PRF                      $prf";  
Use_rsc "$outfile write WAVELENGTH               $wavelength";   
Use_rsc "$outfile write PULSE_LENGTH             $pulse_length";  
Use_rsc "$outfile write CHIRP_SLOPE              $chirp_slope";
Use_rsc "$outfile write RANGE_SAMPLING_FREQUENCY $sampling_freq"; 
Use_rsc "$outfile write RANGE_PIXEL_SIZE         $rng_pixel_size";
Use_rsc "$outfile write AZIMUTH_PIXEL_SIZE       $azimuth_pixel_size";  
Use_rsc "$outfile write DELTA_LINE_UTC           $delta_line_utc";
Use_rsc "$outfile write SQUINT                   $squint";
Use_rsc "$outfile write RANGE_OFFSET             $rng_offset";
Use_rsc "$outfile write FIRST_LINE_UTC           $first_line_utc";
Use_rsc "$outfile write CENTER_LINE_UTC          $center_line_utc";
Use_rsc "$outfile write LAST_LINE_UTC            $last_line_utc";
Use_rsc "$outfile write RLOOKS                   1";
Use_rsc "$outfile write ALOOKS                   1";
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

$xmin = $xmin/2;

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
First line to read                               (-) = $ymin                           ! 
Number of patches                                (-) = $number_of_patches              ! org 40, 1 for testing
First sample pair to use                         (-) = $xmin                           ! 
Number of valid pulses                           (-) = $valid_samples                  ! 
Azimuth Patch Size                               (-) = $patch_size                     !
Deskew?                                          (-) = y                               ! 
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
Azimuth resolution                               (m) = 5                    ! 
Number of azimuth looks                          (-) = 1                               ! 
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

B<roi_prep.pl> I<date>

date: radical of the input file I<date>.raw

Looks for I<date>.proc, then roi.proc in the next directory up, then uses default values

 defaults:
 before_z_ext = 1500
 after_z_ext  = 1500
 near_rng_ext = 700
 far_rng_ext  = 0

Can also set number_of_patches

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

=head1 LAST UPDATE

Frederic Crampe, Jan 10, 2000

=cut
