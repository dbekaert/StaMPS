#!/usr/bin/perl
### make_raw_ers.pl
## 28-Oct-2010 Adapted for Stamps by MCC
#
use Env qw(INT_SCR INT_BIN);
use lib "$INT_SCR";  #### Location of Generic.pm
use Generic;
use POSIX qw(ceil floor);

###Usage info/check
sub Usage{

`$INT_SCR/pod2man.pl  $INT_SCR/make_raw.pl`;
exit 1;
}
@ARGV >= 3  or Usage();
@args = @ARGV;

$orbit_type        = shift;
$leader_file       = shift;
$outname           = shift ;
$reference_counter = shift;

# get list of image file inputs
@imagery=split /\s+/, `ls $leader_file*` or die "No Imagery files\n";
print STDERR "@imagery \n";

#################
Message "Checking I/O";
#################
@Infiles  = ($leader_file, @imagery);
@Outfiles = ("$outname.raw",  "$outname.raw.rsc");
&IOcheck(\@Infiles, \@Outfiles);
Log ("make_raw.pl", @args);

###########################
Message "General definitions";
###########################

$C                        = 299792458;
# new value
$ONE_WAY_DELAY            = 3.4449923715353358e-06;
# old value
#$ONE_WAY_DELAY            = -2.0168e-06;
$RANGE_SAMPLING_FREQUENCY = 18.950420e6;
$ANTENNA_SIDE             = -1;
$ANTENNA_LENGTH           = 12;
$SEC_PER_PRI_COUNT        = 210.94e-09;

###############
# FR 20020903
$ANTENNA_LENGTH           = 10;

###############
#Soren's paper
$RANGE_SAMPLING_FREQUENCY = 18.962468e6;
$SEC_PER_PRI_COUNT        = 210.943006e-9;
$ONE_WAY_DELAY            = 6.622e-6 / 2;
#chirp_length 37.12e-6
#chirp bandwidth 15.50829e6;

$PLANET_GM       = 3.98600448073E+14;
$PLANET_SPINRATE = 7.29211573052E-05;

$RPVersion = 2.3;  #ROI_PAC version

#############################
Message "Getting facility name";
#############################

$inputnum       = @imagery;
print "number of input imagery files, $inputnum\n";
for ($j = 0; $j < $inputnum; ++$j){
 open(RAW,$imagery[$j]) or die "Can't open Imagery file\n";
 for ($i=1; $i<79; $i++) {
   ($name,$value) = split /\=/,<RAW>;
   chomp($value);
   #Message "Name: $name";
   #Message "Value: $value";
   if ($j==0 && $name eq "PRODUCT") {($foo,$time) = split /\"/,$value;}
   if ($name eq "SOFTWARE_VER") {($foo,$proc_ver) = split /\"/,$value;}
   if ($name eq "SWATH") {($foo,$swath)=split /\"/,$value;}
   if ($name eq "TX_RX_POLAR") {($foo,$polar)=split /\"/,$value;}
   if ($j==0 && $name eq "SENSING_START") {($foo,$stime)=split /\"/,$value;}
   if ($j==($inputnum-1) && $name eq "SENSING_STOP") {($foo,$sense_sp)=split /\"/,$value;}
   if ($name eq "PROC_CENTER") {($foo,$facility)=split /\"/,$value;}
   if ($name eq "REL_ORBIT") {($foo,$track)=split /\+/,$value;}
   if ($name eq "ABS_ORBIT") {($foo,$orbit_num)=split /\+/,$value;}
 }
 close(RAW);
}

print "orbit number $orbit_num, track $track\n";
print "sensing start $stime, sensing stop $sense_sp\n";

# ERS in Envisat format header does not contain frame number
$frame=0;

$dt=substr($time,14,8);
$yr=substr($time,14,4);
$mo=substr($time,18,2);
$da=substr($time,20,2);
$hr=substr($time,23,2);
$mn=substr($time,25,2);
$sc=substr($time,27,2);
$ms=substr($stime,-6)/1000;
$outname or $outname=$dt;  # use date if not set in call
print STDERR " $yr $mo $da $hr $mn $sc $ms\n";
print STDERR " $proc_ver \n";
print STDERR "Facility: $facility \n";
print STDERR "Acqsn date: $dt \n";

#$version    = ByteRead ($leader_file, 1782, 8);

$sat_num = substr $time, 61, 1 ;
$sat = "ERS$sat_num";
print "satellite $sat\n";

$width = 11498;
$xmin = 266;
$xmax = 11498;
$LineCounterFirstByte=54;
$MPHpSPH=3203;
$SWSToffset=58;
$PRIoffset=60;

##################################
  Message "Finding reference counter, pri and prf";
##################################
$num_bytes    = 2;     #both SWST and PRI are 2-byte integers.

$imagery_file = $imagery[0]; 

$pri_count    = ByteRead ($imagery_file, $MPHpSPH+$PRIoffset+$width,  $num_bytes);
$swst_counter = ByteRead ($imagery_file, $MPHpSPH+$SWSToffset+$width, $num_bytes);

$swst_counter = unpack("n",$swst_counter);
$pri_count    = unpack("n",$pri_count);

print STDERR "swst counter read: $swst_counter\n";
print STDERR "swst offset read: $SWSToffset\n";

$pri = ($pri_count + 2.0) * $SEC_PER_PRI_COUNT;
$prf = 1/$pri;
print STDERR "PRI: $pri\n";
print STDERR "PRF: $prf\n";

unless ($reference_counter){
  $reference_counter = $swst_counter;
}
Message "Reference counter = $reference_counter";
print STDERR "Reference counter = $reference_counter";

############################
Message "Checking line number";
############################
$LineCounterLength = 4;
$RemoveInputFiles  = 1;

$parse_args="$width $xmin $LineCounterFirstByte $LineCounterLength $MPHpSPH $RemoveInputFiles";
$call_parse="$INT_BIN/new_parse_ers $parse_args @imagery tmp_IMAGERY.raw" ;

Message "$call_parse";
`$call_parse`;

################################
Message "Writing the imagery resource file";
##############################################
###Calculate file length
$size         = -s "tmp_IMAGERY.raw" or die "tmp_IMAGERY.raw has zero size\n";
$file_length = $size/$width;
print STDERR " File length : $file_length";
$clength=0;  # Envisat provides first sample line time, not scene center, so we fib here
$rms=floor($ms);

$chirp = 0.419137466e12;
$wavelength = 0.0565646;
$pulse_length = 37.10e-06;

Use_rsc "tmp_IMAGERY.raw write FIRST_FRAME                 $frame";
Use_rsc "tmp_IMAGERY.raw write FIRST_FRAME_SCENE_CENTER_TIME            $yr$mo$da$hr$mn$sc$rms";
Use_rsc "tmp_IMAGERY.raw write FIRST_FRAME_SCENE_CENTER_LINE            $clength";
Use_rsc "tmp_IMAGERY.raw write DATE                                     $dt";
Use_rsc "tmp_IMAGERY.raw write FIRST_LINE_YEAR                          $yr";
Use_rsc "tmp_IMAGERY.raw write FIRST_LINE_MONTH_OF_YEAR                 $mo";
Use_rsc "tmp_IMAGERY.raw write FIRST_LINE_DAY_OF_MONTH                  $da";
Use_rsc "tmp_IMAGERY.raw write FIRST_CENTER_HOUR_OF_DAY                 $hr";
Use_rsc "tmp_IMAGERY.raw write FIRST_CENTER_MN_OF_HOUR                  $mn";
Use_rsc "tmp_IMAGERY.raw write FIRST_CENTER_S_OF_MN                     $sc";
Use_rsc "tmp_IMAGERY.raw write FIRST_CENTER_MS_OF_S                     $ms";
Use_rsc "tmp_IMAGERY.raw write PROCESSING_FACILITY                      $facility";
Use_rsc "tmp_IMAGERY.raw write PROCESSING_SYSTEM                        $proc_ver";
#Use_rsc "tmp_IMAGERY.raw write PROCESSING_VERSION                      $sat";
Use_rsc "tmp_IMAGERY.raw write WAVELENGTH                        $wavelength";
Use_rsc "tmp_IMAGERY.raw write PULSE_LENGTH                      $pulse_length";
Use_rsc "tmp_IMAGERY.raw write CHIRP_SLOPE                        $chirp";
Use_rsc "tmp_IMAGERY.raw write I_BIAS                        15.5";
Use_rsc "tmp_IMAGERY.raw write Q_BIAS                        15.5";

###############################
#Message "Reading the leader file";
###############################
#$call_leader="$INT_BIN/leader2rsc $leader_file $INT_SCR/format_leaderfile_$facility tmp_IMAGERY.raw.rsc"; 
#
#Message "$call_leader";
#`$call_leader`;

#Add some lines at top and bottom
$pad_top    = 0;
$pad_bottom = 0;
$call_shift="$INT_BIN/delay_shift tmp_IMAGERY.raw tmp_IMAGERY.new shift.out $SWSToffset $reference_counter $width $xmin $xmax $pad_top $pad_bottom";

Message   "$call_shift";
$status = `$call_shift`;
if ($status =~ /No good range counter found/){
  die "No good range counter found using reference counter = $reference_counter\n";
}


`mv tmp_IMAGERY.new tmp_IMAGERY.raw`;
`cp shift.out shift.out.rsc`;
### Changed shift.out to shift.out.rsc so that Use_rsc can read it

$swst_counter = Use_rsc "shift.out read SWST_COUNTER";
$width        = Use_rsc "shift.out read WIDTH";
$xmax         = $width;
###Find true length using new width
$size         = -s "tmp_IMAGERY.raw" or die "tmp_IMAGERY.raw has zero size\n";
$file_length  = $size/$width;

$first_line_utc   = ($hr*60+$mn)*60+$sc+$ms/1000.;
$last_line_utc    = $first_line_utc+$pri*$file_length;
$center_utc       = ($first_line_utc+$last_line_utc)/2;

################################
Message "Writing the imagery resource file";
################################
$range_pixel_size = $C / $RANGE_SAMPLING_FREQUENCY /2;
$swst             = 9*$pri+$swst_counter*$SEC_PER_PRI_COUNT;
$starting_range   = ($swst/2-($ONE_WAY_DELAY))*$C;

Use_rsc "tmp_IMAGERY.raw write PLATFORM                 $sat";
Use_rsc "tmp_IMAGERY.raw write BEAM                     $swath";
Use_rsc "tmp_IMAGERY.raw write POLARIZATION             $polar";
Use_rsc "tmp_IMAGERY.raw write ORBIT_NUMBER             $orbit_num";
Use_rsc "tmp_IMAGERY.raw write RANGE_BIAS               $range_bias";
Use_rsc "tmp_IMAGERY.raw write RANGE_PIXEL_SIZE         $range_pixel_size";
Use_rsc "tmp_IMAGERY.raw write PRF                      $prf";
Use_rsc "tmp_IMAGERY.raw write ANTENNA_SIDE             $ANTENNA_SIDE";
Use_rsc "tmp_IMAGERY.raw write ANTENNA_LENGTH           $ANTENNA_LENGTH";
Use_rsc "tmp_IMAGERY.raw write FILE_LENGTH              $file_length";
Use_rsc "tmp_IMAGERY.raw write WIDTH                    $width";
Use_rsc "tmp_IMAGERY.raw write YMIN                     0";
Use_rsc "tmp_IMAGERY.raw write YMAX                     $file_length";
Use_rsc "tmp_IMAGERY.raw write RANGE_SAMPLING_FREQUENCY $RANGE_SAMPLING_FREQUENCY";
Use_rsc "tmp_IMAGERY.raw write PLANET_GM                $PLANET_GM";
Use_rsc "tmp_IMAGERY.raw write PLANET_SPINRATE          $PLANET_SPINRATE";

Use_rsc "tmp_IMAGERY.raw write FIRST_LINE_UTC  $first_line_utc";
Use_rsc "tmp_IMAGERY.raw write CENTER_LINE_UTC $center_utc";
Use_rsc "tmp_IMAGERY.raw write LAST_LINE_UTC   $last_line_utc";
Use_rsc "tmp_IMAGERY.raw write ONE_WAY_DELAY            $ONE_WAY_DELAY"; 
Use_rsc "tmp_IMAGERY.raw write STARTING_RANGE           $starting_range"; 
Use_rsc "tmp_IMAGERY.raw write XMIN                     $xmin";
Use_rsc "tmp_IMAGERY.raw write XMAX                     $xmax";
Use_rsc "tmp_IMAGERY.raw write WIDTH                    $xmax";

######################################################################################
Message "Reading state vectors in SARLEADER's header, Building hdr_data_points_$outname file"; 
######################################################################################
$day   = Use_rsc "tmp_IMAGERY.raw read FIRST_LINE_DAY_OF_MONTH"; 
$month = Use_rsc "tmp_IMAGERY.raw read FIRST_LINE_MONTH_OF_YEAR"; 
$year  = Use_rsc "tmp_IMAGERY.raw read FIRST_LINE_YEAR";

open HDR, ">hdr_data_points_$outname.rsc" or die "Can't write to hdr_data_points_$outname.rsc\n";
 
if ($orbit_type eq "HDR"){  
# Attention : HDR ne marche pas !
 $numofsarl=$#ld_file+1;
 $countsarl=0;
 foreach $imagery (@imagery){
  $countsarl=$countsarl+1;
 open(RAW,$imagery[$countsarl]) or die "Can't open Imagery file\n";
 for ($i=1; $i<79; $i++) {
   ($name,$value) = split /\=/,<RAW>;
   chomp($value);
   #Message "Name: $name";
   #Message "Value: $value";
   if ($name eq "STATE_VECTOR_TIME") {($foo,$time_data1) = split /\"/,$value;}
   if ($name eq "X_POSITION") {($foo,$xx[0])=split /\"/,$value;}
   if ($name eq "Y_POSITION") {($foo,$yy[0])=split /\"/,$value;}
   if ($name eq "Z_POSITION") {($foo,$zz[0])=split /\"/,$value;}
   if ($name eq "X_VELOCITY") {($foo,$vvx[0])=split /\"/,$value;}
   if ($name eq "Y_VELOCITY") {($foo,$vvy[0])=split /\"/,$value;}
   if ($name eq "Z_VELOCITY") {($foo,$vvz[0])=split /\+/,$value;}
 }
 close(RAW);
 print HDR "$time_data[$i] $xx[$i] $yy[$i] $zz[$i] $vvx[$i] $vvy[$i] $vvz[$i]\n";
  }
}else {
 for ($i=0;$i<5;$i++){
  $time=($last_line_utc-$first_line_utc)*$i/4+$first_line_utc;
  ($q1,$q2,$q3,$q4,$q5, $x, $y, $z, $vx, $vy,$vz) = split /\s+/,
    `$INT_SCR/state_vector.pl $year$month$day $time $sat $orbit_type $date`;
  Status "state_vector.pl";
  print HDR "$time $x $y $z $vx $vy $vz\n";
 }
}

close(HDR);


###############################
Message "Using Orbit Information"; 
###############################
#Stamps : Lat and Lon changed to Latd and Lond, MCC
($q1,$q2,$Latd,$Lond,$height_mid, $x0, $y0, $z0, $vx0, $vy0,$vz0) = split /\s+/,
    `$INT_SCR/state_vector.pl $year$month$day $center_utc $sat $orbit_type $date`;
Status "state_vector.pl";

$pi   = atan2(1,1)*4;

$Lat=$Latd*$pi/180.; # convert to radians
$Lon=$Lond*$pi/180.; # MCC


if ($orbit_type eq "HDR"){
 $ae    = 6378137;             #GRS80 reference ellipsoid
 $flat  = 1/298.257223563;
 $r     = sqrt($x0**2+$y0**2+$z0**2);
 $r1    = sqrt($x0**2+$y0**2);
 $Lat   = atan2($z0,$r1);
 $Lon   = atan2($y0,$x0);
 $H     = $r-$ae;
 for ($i=1; $i<7; $i++){
  $N      = $ae/(sqrt(1-$flat*(2-$flat)*sin($Lat)**2));
  $TanLat = $z0/$r1/(1-(2-$flat)*$flat*$N/($N+$H));
  $Lat    = atan2($TanLat,1);
  $H      = $r1/cos($Lat)-$N;
 }
 $height_mid=$H; 
}

$ae   = 6378137;                        #WGS84 reference ellipsoid
$flat = 1./298.257223563;
$N    = $ae/sqrt(1-$flat*(2-$flat)*sin($Lat)**2);
$re_mid=$N;

$ve=-sin($Lon)*$vx0+cos($Lon)*$vy0;
$vn=-sin($Lat)*cos($Lon)*$vx0-sin($Lat)*sin($Lon)*$vy0+cos($Lat)*$vz0;
$hdg = atan2($ve,$vn);
$e2 = $flat*(2-$flat);
$M = $ae*(1-$e2)/(sqrt(1-$e2*sin($Lat)**2))**3;
$earth_radius_mid = $N*$M/($N*(cos($hdg))**2+$M*(sin($hdg))**2);

($q1,$q2,$q3,$q4,$height_top, $x0, $y0, $z0, $vx, $vy,$vz) = split /\s+/,
    `$INT_SCR/state_vector.pl $year$month$day $first_line_utc $sat $orbit_type $date`;
Status "state_vector.pl";

if ($orbit_type eq "HDR"){
 $ae    = 6378137;             #GRS80 reference ellipsoid
 $flat  = 1/298.257223563;
 $r     = sqrt($x0**2+$y0**2+$z0**2);
 $r1    = sqrt($x0**2+$y0**2);
 $Lat   = atan2($z0,$r1);
 $Lon   = atan2($y0,$x0);
 $H     = $r-$ae;
 for ($i=1; $i<7; $i++){
  $N      = $ae/(sqrt(1-$flat*(2-$flat)*sin($Lat)**2));
  $TanLat = $z0/$r1/(1-(2-$flat)*$flat*$N/($N+$H));
  $Lat    = atan2($TanLat,1);
  $H      = $r1/cos($Lat)-$N;
 }
 $height_top=$H; 
}

$height_dt=($height_mid-$height_top)/($center_utc-$first_line_utc);
if ($vz0 > 0) {$orbit_direction =  "ascending";}
else          {$orbit_direction = "descending";}
$velocity_mid=sqrt($vx0**2 + $vy0**2 + $vz0**2);

#$Lat=$Latd*180./$pi;
#$Lon=$Lond*180./$pi;
#$hdg=$hdgd*180./$pi;

$Latd=$Lat*180./$pi;
$Lond=$Lon*180./$pi;
$hdgd=$hdg*180./$pi;


Use_rsc "tmp_IMAGERY.raw write HEIGHT       $height_top";
Use_rsc "tmp_IMAGERY.raw write HEIGHT_DT    $height_dt";
Use_rsc "tmp_IMAGERY.raw write VELOCITY     $velocity_mid";
Use_rsc "tmp_IMAGERY.raw write LATITUDE     $Latd";
Use_rsc "tmp_IMAGERY.raw write LONGITUDE    $Lond";
Use_rsc "tmp_IMAGERY.raw write HEADING      $hdgd";
Use_rsc "tmp_IMAGERY.raw write EQUATORIAL_RADIUS   $ae";
Use_rsc "tmp_IMAGERY.raw write ECCENTRICITY_SQUARED $e2";
#Use_rsc "tmp_IMAGERY.raw write EARTH_EAST_RADIUS $N";
#Use_rsc "tmp_IMAGERY.raw write EARTH_NORTH_RADIUS $M";
Use_rsc "tmp_IMAGERY.raw write EARTH_RADIUS $earth_radius_mid";
Use_rsc "tmp_IMAGERY.raw write ORBIT_DIRECTION $orbit_direction";

#####################################################
$default = "$INT_SCR/default.raw.rsc";
Use_rsc    "tmp_IMAGERY.raw.rsc merge $default";
#####################################################

###############################
Message "Doppler Computation"; 
###############################

$line_0=100;
$line_1=$file_length - 200;
$i_bias = Use_rsc "tmp_IMAGERY.raw read I_BIAS";

$dop_rng1=0;
$dop_rng2=0;
$dop_rng3=0;

if ($xmin > 0) {
  $line_header = $xmin - 1;
}
else {$line_header = 0;}
$lastbyte=$xmax;  # the right edge seems to be all data

open DOP, ">dopiq.in" or die "Can't write to dopiq.in\n";

print DOP <<END;
tmp_IMAGERY.raw
$xmax,$line_header,$lastbyte
$line_0,$line_1
$i_bias,$prf
END
close(DOP);

#`$INT_BIN/dopiq-new <<end
#tmp_IMAGERY.raw
#$xmax,$line_header,$lastbyte
#$line_0,$line_1
#$i_bias,$prf
#end
#`;


#dopiq-new estimated wrong doppler, changed  to dopiq MCC
Message "using dopiq instead dopiq-new";

`$INT_BIN/dopiq <<end
tmp_IMAGERY.raw
$xmax,$line_0,$line_1
$i_bias,$prf
end
`;



$DopUnwStr = `$INT_SCR/DopUnw.pl dop.out`;
$DopUnwStr =~ /Quadratic doppler:\s+(\S+)\s+(\S+)\s+(\S+)/;
@QuadDopCoeff = ( $1, $2, $3 );

$wavelength  = Use_rsc "tmp_IMAGERY.raw read WAVELENGTH";
$sin_theta = sqrt( 1 - ($height_top / $starting_range)**2 );
$fd = $QuadDopCoeff[0] * $prf;
$sin_squint = $fd / ( 2 * $velocity_mid * $sin_theta ) * $wavelength; # neglect vertical velocity

$squint_rad = atan2( $sin_squint, sqrt( 1 - $sin_squint ** 2 ) );
$squint_deg = 180 * $squint_rad / $pi;

Use_rsc "tmp_IMAGERY.raw write DOPPLER_RANGE0  $QuadDopCoeff[0]";
Use_rsc "tmp_IMAGERY.raw write DOPPLER_RANGE1  $QuadDopCoeff[1]";
Use_rsc "tmp_IMAGERY.raw write DOPPLER_RANGE2  $QuadDopCoeff[2]";
Use_rsc "tmp_IMAGERY.raw write DOPPLER_RANGE3  0.";
Use_rsc "tmp_IMAGERY.raw write SQUINT          $squint_deg";
Use_rsc "tmp_IMAGERY.raw write ROI_PAC_VERSION     $RPVersion";


`mv tmp_IMAGERY.raw_parse_line.out  ${outname}_parse_line.out`;
`mv tmp_IMAGERY.raw.rsc             ${outname}.raw.rsc`;
`mv tmp_IMAGERY.raw                 ${outname}.raw`;

#########################
Message "Raw data ready for processing";
#########################

exit 0;

=pod

=head1 USAGE

B<make_raw.pl> I<orbit_type leader_file outname [facility] [reference_counter] >

reference_counter: default is read from imagery file

=head1 FUNCTION

Creates I<outname>.raw and I<outname>.raw.rsc from imagery files

=head1 ROUTINES CALLED

delay_shift

leader2rsc

new_parse_ers

state_vector.pl

dd

dopiq

=head1 CALLED BY

none

=head1 FILES USED

IMAGERY*

SARLEADER*

=head1 FILES CREATED

I<outname>.raw

I<outname>.raw.rsc

I<outname>_parse_line.out

shift.out

shift.out.rsc

=head1 HISTORY

Shell Script : Francois ROGEZ 96/98
Perl  Script : Rowena LOHMAN 04/18/98
Modifications: Frederic CRAMPE, Oct 21, 1998

=head1 LAST UPDATE

Frederic CRAMPE, Oct 06, 1999

=cut

