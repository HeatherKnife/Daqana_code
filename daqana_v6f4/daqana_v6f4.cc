// 2019/10/25 (DT). The new version daqana_v6f1 is created by Diego, starting from the previous version (daqana_v6e4) that was used by Diego for the analysis of PPACs, and that is stored in the Dropbox folder of the Medley team. This present version (daqana_v6f1) is a modification of that one to try to implement some of the changes that were done by Bernardo, and that were saved, unfortunately, under the same name "daqana_v6e4" in the VERDI folder. Therefore, this new version does not use the code by Bernardo, except for a few lines of code that can be implemented here.


//2019/07/01 (DT). I modify it to analyze one Silicon detector (channel A) and three PPACs (ch. B, C, D).
// The original version for that purpose (daqana_v6e3) has been lost when the hard disk of the computer at the lab crashed on 2019/06/17. No backup existed.
// Therefore, I create the present version (daqana_v6e4) starting from the daqana_v6e and modifying it to use 3 PPACs and 1 Silicon. (The "3" in daqana_v6e3 meant, originally, "3 PPACs", so that there were no versions "v6e2" or similar;  daqana_v6e3 was a direct modification of daqana_v6e).

//2018/10/15 (DT). New variables "chargeRaw" and "chargeFiltered" to store the integral of the pulse (i.e., the charge in case of a signal from a current-sensitive amplifier).
#include <cassert>
#include <cmath>
#include <cstring>
#include <fftw3.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <map>
#include <stdexcept>
#include <stdint.h>
#include <TGraph.h>
#include <TF1.h>
#include <TString.h>
#include <vector>
#include <snappy.h>
#include "RootFileManager.hh"

//#define VERSION 6
//const uint32_t header_size[VERSION+1] = {0,29,30,30,38,46,50};

// DT 2019/10/25. We have to write VERSION 7 because the data acquisition was done using the version 7 from Benjaminas's laptop.
 #define VERSION 7  //is still versione 6, I just changed the version in order to make it work on my computer (BB 28/08/2019)
 const uint32_t header_size[VERSION+1] = {0,29,30,30,38,46,50,52};
//
// End of DT 2019/10/25.



/* File format
// VERSION 1:
  Header (29 bytes)
  {
    UInt32 Version
    UInt64 Date (in decimal: yyyymmddhhmmss)
    UInt8  Active channels
    UInt8  Active triggers
    SInt16 Trigger level
    UInt32 Number of records
    UInt32 Samples per record
    UInt32 Pre trigger samples
    UInt8  Interleaving setting
  }
  Records[Number of records]
  {
    UInt64 Timestamp
    channels[Number of active channels]
    {
      SInt16 samples[Samples per record] // 12 bit value stored with 16 bit
    }
  }
// VERSION 2:
  Header (30 bytes)
  {
    UInt32 Version
    UInt64 Date (in decimal: yyyymmddhhmmss)
    UInt8  Active channels
    UInt8  Active triggers
    SInt16 Trigger level
    UInt32 Number of records
    UInt32 Samples per record
    UInt32 Pre trigger samples
    UInt8  Interleaving setting
    SInt8  Trigger Edge
  }
  Records[Number of records]
  {
    UInt64 Timestamp
    channels[Number of active channels]
    {
      SInt16 samples[Samples per record] // 12 bit value stored with 16 bit
    }
  }
// VERSION 3: (Identical to version 2)
  Header (30 bytes)
  {
    UInt32 Version
    UInt64 Date (in decimal: yyyymmddhhmmss)
    UInt8  Active channels
    UInt8  Active triggers
    SInt16 Trigger level
    UInt32 Number of records
    UInt32 Samples per record
    UInt32 Pre trigger samples
    UInt8  Interleaving setting
    SInt8  Trigger Edge
  }
  Records[Number of records]
  {
    UInt64 Timestamp
    channels[Number of active channels]
    {
      SInt16 samples[Samples per record] // 12 bit value stored with 16 bit
    }
  }
// VERSION 4: (Adds DC-offsets)
  Header (38 bytes)
  {
    UInt32    Version
    UInt64    Date (in decimal: yyyymmddhhmmss)
    UInt8     Active channels
    UInt8     Active triggers
    SInt16    Trigger level
    UInt32    Number of records
    UInt32    Samples per record
    UInt32    Pre trigger samples
    UInt8     Interleaving setting
    UInt8     Trigger Edge
    SInt16[4] DC offsets for channel A, B, C, D
  }
  Records[Number of records]
  {
    UInt64 Timestamp
    channels[Number of active channels]
    {
      SInt16 samples[Samples per record] // 12 bit value stored with 16 bit
    }
  }
// VERSION 5: (Adds DC-offsets)
  Header (46 bytes)
  {
    UInt32    Version
    UInt64    Date (in decimal: yyyymmddhhmmss)
    UInt64    Acquisition duration in microsecs
    UInt8     Active channels
    UInt8     Active triggers
    SInt16    Trigger level
    UInt32    Number of records
    UInt32    Samples per record
    UInt32    Pre trigger samples
    UInt8     Interleaving setting
    UInt8     Trigger Edge
    SInt16[4] DC offsets for channel A, B, C, D
  }
  Records[Number of records]
  {
    UInt64 Timestamp
    channels[Number of active channels]
    {
      SInt16 samples[Samples per record] // 12 bit value stored with 16 bit
    }
  }
VERSION 6: (Added Compression and Trigger offset capability)
  Header (50 bytes)
  {
    UInt32    Version
    UInt64    Start time (decimal value: yyyymmddhhmmss)
    UInt64    Acquisition duration in microsecs
    UInt8     Active channels
    UInt8     Active triggers
    SInt16    Trigger level
    UInt16    Tfrigger level offset
    UInt32    Number of records
    UInt32    Samples per record
    UInt32    Pre trigger samples
    UInt8     Interleaving setting
    UInt8     Trigger edge setting
    UInt16     1 compress /  0 do not compress 
    SInt16[4] DC-offset for channel A, B, C, D
  }
  Records[Number of records]
  {
    UInt64 Timestamp
    channels[Number of active channels]
    {
      SInt16 samples[Samples per record] // 12 bit value stored with 16 bit
    }
  }
*/

/* Header variables */
uint32_t file_format_version;
uint64_t date;
uint64_t acq_duration;
uint8_t active_ch;
uint8_t active_tr;
int16_t tr_lvl;
int16_t tr_lvl_off;
uint8_t tr_edge;
uint32_t n_rec;
uint32_t smp_per_rec;
uint32_t pre_tr_smp;
uint8_t interleaving;
uint16_t compress;
double ns_per_smp;
int16_t dc_offset[4];

/* Other globals */
uint32_t do_store;
const uint32_t max_stored_events = 10;
const uint8_t ch_mask[4] = {1,2,4,8};
/* the _only_ arrays and plans to use with FFTW */
fftw_complex* waveform_t;
fftw_complex* waveform_f;
fftw_plan fftw_forward_plan;
fftw_plan fftw_backward_plan;

/* Function declarations */
void    read_header( std::ifstream& );
bool    get_active_ch( uint8_t ch ) { return active_ch & ch_mask[ch]; }
bool    get_active_tr( uint8_t ch ) { return active_tr & ch_mask[ch]; }
template<class A> double baseline( const A* );
template<class A> double baseline( const A*, uint32_t, uint32_t );
template<class A> double baseline_var( const A*, double );
template<class A> double check_pileup( const A*, double, std::vector< std::pair<double,double> >&, int32_t polarity );
template<class A> double get_fwhm( const A*, uint32_t, double );
template<class A> double fit_line( const A*, double&, double&, int32_t, int32_t );
template<class A> double get_min(  const A*, double& );
template<class A> double get_min(  const A*, double&, uint32_t, uint32_t );
template<class A> double get_extremum_polyfit( const A*, double&, uint32_t, uint32_t );
template<class A> double get_max(  const A*, double& );
template<class A> double get_max(  const A*, double&, uint32_t, uint32_t );
template<class A> double left_threshold( const A*, uint32_t, double );
template<class A> double right_threshold( const A*, uint32_t, double );
template<class A> double threshold( const A*, double );
template<class A> double integratePeak( const A*, uint32_t, uint32_t, double);

template<class A> double CFD( const A*, double*, double, uint32_t , double, int );  // DT 2019/10/25 
template<class A> double Leading_edge(const A*, int); // DT 2019/10/25

double linear_interpolation( double, double, double, double, double );
double quadratic_extremum( double, double, double, double& );

// void    handle_record( uint32_t, int16_t**, double**, double** );
void    handle_record( uint32_t, int16_t**, double**, double**, double ** ); // DT 2019/10/25. To compare two ways of doing CFD.  

template<class A> void store(    const A*, const TString );
                  void store(    const fftw_complex*, const TString );
template<class A> void rectangle_smooth( const A*, double*, uint32_t );
template<class A> void triangle_smooth( const A*, double*, uint32_t );
template<class A> void savitzky_golay_quad_smooth( const A*, double*, uint32_t );
template<class A> void subtract( const A*, const double*, double* );
template<class A> void diff(     const A*, double* );
template<class A> void MWD(const A* , double* , int);
template<class A> void Trapezoidal_Filter(const A* , double* , int,int);
template<class A> void moving_average(const A* , double* , int);
fftw_complex* lowpass(    fftw_complex*, double );
fftw_complex* highpass(   fftw_complex*, double );
fftw_complex* bandpass(   fftw_complex*, double, double );
fftw_complex* bandfilter( fftw_complex*, double, double );


/* FFT functions */
void fftw_init();
fftw_complex* fftw_forward( fftw_complex* );
fftw_complex* fftw_backward( fftw_complex* );
fftw_complex* to_freq( const int16_t* );
fftw_complex* to_freq( const double* );
fftw_complex* to_freq( const fftw_complex* );
fftw_complex* to_time( const int16_t* );
fftw_complex* to_time( const fftw_complex* );

// Handles for root file manager
handle h_ts;
handle h_pileup;
handle h_dist;
std::map< int, handle > h_time;
//std::map< int, handle > h_time_cfd;
std::map< int, handle > h_time_cfd_raw; // DT 2019/10/25. To compare two ways of doing CFD.
std::map< int, handle > h_time_cfd_smooth;  // DT 2019/10/25. To compare two ways of doing CFD.
std::map< int, handle > h_fwhm;
std::map< int, handle > h_ph;
std::map< int, handle > h_max;
std::map< int, handle > h_bl_pre;
std::map< int, handle > h_bl_post;
std::map< int, handle > h_bl_k;
std::map< int, handle > h_bl_m;
std::map< int, handle > h_bl_chi2;
std::map< int, handle > h_chargefiltered;
std::map< int, handle > h_chargeraw;
std::map< int, handle > h_leftlimit;
std::map< int, handle > h_rightlimit;


int main(int argc, char** argv)
{
  assert( argc > 2 );

  /* Open the output waveform file */
  TString fname(argv[argc-1]);
  RootFileManager* rm = RootFileManager::GetInstance();
 // std::cout << "Writing to file: " << fname << std::endl;
  rm->OpenFile( fname, "RECREATE" );
  rm->AddFolder("Waveforms");
  handle h_dtree = rm->NewTree("EventData");

  std::ifstream input(argv[1], std::ifstream::binary);
  /* Read data header of first file
   * it is assumed all other files have the same settings! */
  read_header( input );
  int16_t* raw[4];
  static char* temporary  = new char[snappy::MaxCompressedLength(smp_per_rec*sizeof(int16_t))];  // temporary buffer if reading compressed data
  double* smooth[4];
//  double* bipolar[4]; // DT 2019/10/25
  double* bipolar_raw[4]; // DT 2019/10/25 To compare two ways of doing CFD.
  double* bipolar_smooth[4];  // DT 2019/10/25 To compare two ways of doing CFD.
  h_ts     = rm->NewInt64( h_dtree, "Timestamp" );
  h_pileup = rm->NewFloat( h_dtree, "Pileup" );
  h_dist   = rm->NewFloat( h_dtree, "PileupDistance" );
  TString ch_name[4];
  ch_name[0] = "A";
  ch_name[1] = "B";
  ch_name[2] = "C";
  ch_name[3] = "D";
 // std::cout << "max_stored_events: " << max_stored_events << std::endl;

  for( int ch = 0; ch < 4; ++ch )
    if( get_active_ch(ch) )
    {
      raw[ch]       = new int16_t[smp_per_rec];
      smooth[ch]    = new double[smp_per_rec];

//    bipolar[ch]    = new double[smp_per_rec]; 
      bipolar_raw[ch]    = new double[smp_per_rec];  // DT 2019/10/25. To compare two ways of doing CFD.
      bipolar_smooth[ch]    = new double[smp_per_rec];  // DT 2019/10/25. To compare two ways of doing CFD.

      h_time[ch]    = rm->NewFloat( h_dtree, "Time_"         + ch_name[ch] );
//      h_time_cfd[ch]    = rm->NewFloat( h_dtree, "Time_CFD_"         + ch_name[ch] );
      h_time_cfd_raw[ch]    = rm->NewFloat( h_dtree, "Time_CFD_raw_"         + ch_name[ch] ); // DT 2019/10/25. To compare two ways of doing CFD.
      h_time_cfd_smooth[ch]    = rm->NewFloat( h_dtree, "Time_CFD_smooth_"         + ch_name[ch] ); // DT 2019/10/25. To compare two ways of doing CFD.
      h_fwhm[ch]    = rm->NewFloat( h_dtree, "FWHM_"         + ch_name[ch] );
      h_ph[ch]      = rm->NewFloat( h_dtree, "PH_"           + ch_name[ch] );
      h_max[ch]     = rm->NewFloat( h_dtree, "max_"          + ch_name[ch] );
      h_chargefiltered[ch]  = rm->NewFloat( h_dtree, "ChargeFiltered_"       + ch_name[ch] );
      h_chargeraw[ch]       = rm->NewFloat( h_dtree, "ChargeRaw_"       + ch_name[ch] );
      h_leftlimit[ch]       = rm->NewFloat( h_dtree, "leftLimit_"       + ch_name[ch] );
      h_rightlimit[ch]        = rm->NewFloat( h_dtree, "rightLimit_"       + ch_name[ch] );
      h_bl_pre[ch]  = rm->NewFloat( h_dtree, "BL_pre_"       + ch_name[ch] );
      h_bl_post[ch] = rm->NewFloat( h_dtree, "BL_post_"      + ch_name[ch] );
      h_bl_k[ch]    = rm->NewFloat( h_dtree, "BL_Slope_"     + ch_name[ch] );
      h_bl_m[ch]    = rm->NewFloat( h_dtree, "BL_Intercept_" + ch_name[ch] );
      h_bl_chi2[ch] = rm->NewFloat( h_dtree, "BL_RedChi2_"   + ch_name[ch] );
    }
  
  for( int file = 0; file < argc-2; ++file )
  {
  //  std::cout << "Reading file: " << argv[file+1] << std::endl;
    std::ifstream input(argv[file+1], std::ifstream::binary);
    input.ignore(header_size[file_format_version]); // Skip header

    /* Read data file to the end */

    for( unsigned int rec = 0; rec < n_rec; ++rec )
    {
     uint64_t len;
      /* Read one full record */
      if( !input.good() )
        throw std::runtime_error("Error in input stream.");
      uint64_t timestamp;

      input.read( reinterpret_cast<char*>(&timestamp), sizeof(timestamp) );
      rm->FillInt64( h_ts, timestamp );

      for( int ch = 0; ch < 4; ++ch ){
        if( get_active_ch(ch) && compress ==0 ){
            
          input.read( reinterpret_cast<char*>(raw[ch]),
                      sizeof(raw[ch][0])*smp_per_rec );
        }
        else if ( get_active_ch(ch) && compress ==1 ){
          input.read( reinterpret_cast<char*>(&len), sizeof(len) );  
            
          input.read( temporary,    len);  // reading the length of the compressed data i will read later

        //
        // Uncompressing the vlaues

        snappy::RawUncompress(temporary, len, reinterpret_cast<char*>(raw[ch]));    
        
        }
      }

      do_store = false;
      for( int ch = 0; ch < 4; ++ch )
        if( get_active_ch(ch) )
        {
          /* Handle everything that concerns only this record */
          try{
//            handle_record(ch,raw,smooth,bipolar);
  handle_record(ch,raw,smooth,bipolar_raw,bipolar_smooth);// DT 2019/10/25. To compare two ways of doing CFD.  
        }
          catch( const char* e )
          {
            std::cerr << e << std::endl;
            return -1;
          }
        }

      /* Store waveforms */
      if( do_store )
      {
        static uint32_t n_stored = 0;
        if( n_stored < max_stored_events )
        {
          TString rec_name("file_");
          rec_name += file;
          rec_name += "_rec_";
          rec_name += rec;
          rm->Cd();
          rm->Cd("Waveforms");
          rm->AddAndCd(rec_name);
          for( int ch = 0; ch < 4; ++ch )
            if( get_active_ch(ch) )
            {
              TString name = rec_name + "_ch_";
              name += ch;
              store( raw[ch], "raw_"+name );
              store( smooth[ch], "smooth_"+name );
           //   store( bipolar[ch], "bipolar_"+name );// DT 2019/10/25. To compare two ways of doing CFD.  
              store( bipolar_raw[ch], "bipolar_raw_"+name );// DT 2019/10/25. To compare two ways of doing CFD.  
              store( bipolar_smooth[ch], "bipolar_smooth_"+name );// DT 2019/10/25. To compare two ways of doing CFD.  
            }
          ++n_stored;
        }
      }
      rm->CloseEntry(h_dtree);
    }
  }

  /* Close waveform file file */
  rm->CloseFile();
//  std::cout << "Analysis done!" << std::endl << std::endl;
}

//void handle_record( uint32_t ch, int16_t** raw, double** smooth, double** bipolar)
void handle_record( uint32_t ch, int16_t** raw, double** smooth, double** bipolar_raw, double** bipolar_smooth) // DT 2019/10/25. To compare two ways of doing CFD.  
{
  RootFileManager* rm = RootFileManager::GetInstance();

  static double* temp = new double[smp_per_rec];

  /* Channel specific actions */
  // Channel A
  // Time signal, peak position is important
/* 
  if( ch > 10 ) // Dummy value to avoid this case. Change when needed;
  {
    rectangle_smooth( raw[ch], smooth[ch], 11 );
    triangle_smooth( smooth[ch], temp, 11 );
    triangle_smooth( temp, smooth[ch], 11 );
  }
  // Channel B
  // Channel C
  // Channel D
*/
  // Energy signals, peak height is important
  if(ch<4)
    {
    rectangle_smooth( raw[ch], smooth[ch], 11 );
    savitzky_golay_quad_smooth( smooth[ch], temp, 11 );
    savitzky_golay_quad_smooth( temp, smooth[ch], 11 );
 //   Trapezoidal_Filter( smooth[ch], , 11 );

    }

  // Common to all channels
    // The function baseline receives the pulse in the first entrance, the starting channel to start to calculate the average in the second entrance, and the number of channels for the average to be calculated in the thrird entrance.
  double bl_k, bl_m;
  double pos_min = 0; //This variable is associated with the x value corresponding to the y min value and is started in zero so that returns the real value of x associated with y_min = min
  double min = get_min(smooth[ch], pos_min, 0, pre_tr_smp); //Take care of the pre_trg_smp for every pulse!!! care about how much of the pulse is going through.
  //The program returns in some way the value of x, which is pos_min through a pointer or a direction of memory
  double pos_max = 0;
  double max = get_max(smooth[ch], pos_max, pre_tr_smp/2, smp_per_rec);
  double bl_chi2 = fit_line( raw[ch], bl_k, bl_m, 0, 3*pre_tr_smp/10);
  double bl_pre  = baseline( smooth[ch], 0, 3*pre_tr_smp/10);  //First average for the beggining of the somooth pulse
  double bl_preRaw  = baseline( raw[ch], 0, 3*pre_tr_smp/10);  //First average for the beggining of the raw pulse
  double start_bl_post = smp_per_rec - 3*pre_tr_smp/10;
  double bl_post = baseline(smooth[ch], start_bl_post, smp_per_rec - start_bl_post);
  //double bl_post = baseline(smooth[ch], smp_per_rec-pre_tr_smp/2, pre_tr_smp/2 - 10); //Second average for the end of the smooth pulse
  //double pulse_h_end = baseline(smooth[ch], smp_per_rec - pre_tr_smp/2, sup_limit_PH - (smp_per_rec - pre_tr_smp/2)); // Third average for the end of the pulse minus 10 channels
  double pos;
  double val;
  double time;
//  double time_cfd;
  double time_cfd_raw; // DT 2019/10/25. To compare two ways of doing CFD.
  double time_cfd_smooth;  // DT 2019/10/25. To compare two ways of doing CFD.

  static uint32_t pileup = 0;



if(ch==0) // Silicon signals after PreAmp. Positive signals.
  {
  val = get_max( smooth[ch], pos);

//std::cout << "I am here" << std::endl;
  

// time = left_threshold( smooth[ch], pos, bl_pre + 0.3 * (val - bl_pre) );
//  time_cfd = CFD( smooth[ch], bipolar[ch], 0.3, 20., -1800., 1);  // DT 2019/10/25; positive signal -> pol=1
 time_cfd_raw = CFD( raw[ch], bipolar_raw[ch], 0.3, 20., -1300., 1); // DT 2019/10/25. To compare two ways of doing CFD.
 time_cfd_smooth = CFD( smooth[ch], bipolar_smooth[ch], 0.3, 20., -1300., 1);  // DT 2019/10/25. To compare two ways of doing CFD.


//    time = left_threshold( smooth[ch], pos, 0.5 * val );
//    val = get_extremum_polyfit( raw[ch], pos, 
//        std::max( 0, int32_t(pos / ns_per_smp + 0.5) - 200 ), 
//        std::min( int32_t(smp_per_rec), int32_t(pos / ns_per_smp + 0.5) + 200 ) );

//std::cout << "ch: " << ch << " max= " << val << " pos= " << pos << " time= " << time << std::endl;
  }


// For MCP, use Raw signal for timing;
else if(ch==1||ch==2||ch==3) // MCP signals. Negative signals.
  {
  val = get_min( smooth[ch], pos);
  time = left_threshold( raw[ch], pos, bl_pre + 0.2 * (val - bl_pre) );
//  time_cfd = CFD( smooth[ch], bipolar[ch], 0.3, 1.,  -1800., 0);  // DT 2019/10/25; negative signal -> pol=0

//  time_cfd_raw = CFD( raw[ch], bipolar_raw[ch], 0.2, 300.,  -1800., 0);   // DT 2019/10/25. To compare two ways of doing CFD.
//  time_cfd_smooth = CFD( smooth[ch], bipolar_smooth[ch], 0.1, 300.,  -1800., 0);  // DT 2019/10/25. To compare two ways of doing CFD.
  time_cfd_raw = CFD( raw[ch], bipolar_raw[ch], 0.3, 1.,  -1800., 0);  
  time_cfd_smooth = CFD( smooth[ch], bipolar_smooth[ch], 0.3, 1.,  -1800., 0);  
 
  std::cout << "time_cfd_smooth" << time_cfd_smooth << std::endl;


//    val = get_extremum_polyfit( raw[ch], pos, 
//        std::max( 0, int32_t(pos / ns_per_smp + 0.5) - 200 ), 
//        std::min( int32_t(smp_per_rec), int32_t(pos / ns_per_smp + 0.5) + 200 ) );

//  std::cout << "ch: " << ch << " min= " << val << " pos= " << pos << " time= " << time << std::endl; 
  }




std::vector< std::pair<double, double> > peaks;
pileup = check_pileup( smooth[ch], std::min(val/2, 1.0*tr_lvl), peaks, -1 ); // Trying to avoid getting triggered by the reflection
rm->FillFloat( h_pileup, pileup );
double distance = 0;
if( peaks.size() > 1 )
  distance = ( peaks[1].first - peaks[0].first );
rm->FillFloat( h_dist, distance );


double ph = fabs(bl_pre - bl_post);
//double ph_Gauss = fabs(max - bl_pre);

//double ph=val;

 //std::cout << " bl_pre " << bl_pre << " ph " << ph << " max " << max << " pos_max " << pos_max << " bl_post " << bl_post << std::endl;

double fwhm = get_fwhm( smooth[ch], int32_t(pos / ns_per_smp + 0.5), bl_pre );

double leftLimit  =  left_threshold( smooth[ch],  int32_t(pos/ ns_per_smp + 0.5), bl_pre + 0.3*(val-bl_pre) );
double rightLimit = right_threshold( smooth[ch],  int32_t(pos/ ns_per_smp + 0.5), bl_post + 0.3*(val-bl_post) );
//double rightLimit = leftLimit + 50;
double chargeFiltered = integratePeak(smooth[ch], int32_t(leftLimit / ns_per_smp + 0.5), int32_t(rightLimit / ns_per_smp + 0.5), bl_pre); // Diego 2018/10/15;
double chargeRaw      = integratePeak(   raw[ch], int32_t(leftLimit / ns_per_smp + 0.5), int32_t(rightLimit / ns_per_smp + 0.5), bl_preRaw);  // Diego 2018/10/15;


// std::cout << "ch: " << ch << " time= " << time << " ph= " << ph << std::endl;

  rm->FillFloat( h_time[ch],    time );

//  rm->FillFloat( h_time_cfd[ch],    time_cfd );   // DT 2019/10/25;
  rm->FillFloat( h_time_cfd_raw[ch],    time_cfd_raw );   // DT 2019/10/25. To compare two ways of doing CFD.
  rm->FillFloat( h_time_cfd_smooth[ch],    time_cfd_smooth );   // DT 2019/10/25. To compare two ways of doing CFD.

  rm->FillFloat( h_fwhm[ch],    fwhm );
  rm->FillFloat( h_ph[ch],      ph );
  rm->FillFloat( h_max[ch],  max);
  rm->FillFloat( h_chargefiltered[ch],  chargeFiltered );
  rm->FillFloat( h_chargeraw[ch],  chargeRaw );
  rm->FillFloat( h_leftlimit[ch],  leftLimit );
  rm->FillFloat( h_rightlimit[ch],  rightLimit );
  rm->FillFloat( h_bl_pre[ch],  bl_pre  );
  rm->FillFloat( h_bl_post[ch], bl_post );
  rm->FillFloat( h_bl_k[ch],    bl_k );
  rm->FillFloat( h_bl_m[ch],    bl_m );
  rm->FillFloat( h_bl_chi2[ch], bl_chi2 );
  
  do_store = true;
}

void read_header( std::ifstream& input )
{
  std::cout << "*** BEGIN HEADER ***" << std::endl;
  input.read( reinterpret_cast<char*>(&file_format_version), sizeof(file_format_version) );
  std::cout << "Detected version:    " << file_format_version << std::endl;
  if( VERSION < file_format_version )
    std::cerr << "Warning: file format version exceeds program version!" << std::endl;
  input.read( reinterpret_cast<char*>(&date), sizeof(date) );
  if( file_format_version >= 5 )
    input.read( reinterpret_cast<char*>(&acq_duration), sizeof(acq_duration) );
  input.read( reinterpret_cast<char*>(&active_ch), sizeof(active_ch) );
  input.read( reinterpret_cast<char*>(&active_tr), sizeof(active_tr) );
  input.read( reinterpret_cast<char*>(&tr_lvl), sizeof(tr_lvl) );
  if( file_format_version >= 6 )
    input.read( reinterpret_cast<char*>(&tr_lvl_off), sizeof(tr_lvl_off) );
  input.read( reinterpret_cast<char*>(&n_rec), sizeof(n_rec) );
  input.read( reinterpret_cast<char*>(&smp_per_rec), sizeof(smp_per_rec) );
  input.read( reinterpret_cast<char*>(&pre_tr_smp), sizeof(pre_tr_smp) );
  input.read( reinterpret_cast<char*>(&interleaving), sizeof(interleaving) );
  if( file_format_version >= 2 )
    input.read( reinterpret_cast<char*>(&tr_edge), sizeof(tr_edge) );
  if( file_format_version >= 6 )
    input.read( reinterpret_cast<char*>(&compress), sizeof(compress) );
  if( file_format_version >= 4 )
    for( int i = 0; i < 4; ++i )
      input.read( reinterpret_cast<char*>(dc_offset+i), sizeof(dc_offset[i]) );
  ns_per_smp = interleaving == 0 ? 1.0 : 0.5;

  std::cout << std::setfill('0');
  std::cout << "File recorded:       " << date/10000000000 << "-" <<
    std::setw(2) << (date/100000000)%100 << "-" <<
    std::setw(2) << (date/1000000)  %100 << " " <<
    std::setw(2) << (date/10000)    %100 << ":" <<
    std::setw(2) << (date/100)      %100 << ":" <<
    std::setw(2) <<  date           %100 << std::endl;
  std::cout << "Acq duration:        " << acq_duration << std::endl;
  std::cout << "Active channels:     " << int(active_ch) << std::endl;
  std::cout << "Trigger channels:    " << int(active_tr) << std::endl;
  std::cout << "Trigger level:       " << tr_lvl << std::endl;
  std::cout << "Trigger level offset:" << tr_lvl_off << std::endl;
  if( file_format_version >= 4 ) 
  {
    std::cout << "DC offset channel A: " << dc_offset[0] << std::endl;
    std::cout << "DC offset channel B: " << dc_offset[1] << std::endl;
    std::cout << "DC offset channel C: " << dc_offset[2] << std::endl;
    std::cout << "DC offset channel D: " << dc_offset[3] << std::endl;
  }
  if( file_format_version >= 2 ) 
    std::cout << "Trigger edge:        " << int(tr_edge) << std::endl;
  std::cout << "Number of records:   " << n_rec << std::endl;
  std::cout << "Samples per record:  " << smp_per_rec << std::endl;
  std::cout << "Pre trigger samples: " << pre_tr_smp << std::endl;
  std::cout << "Interleaving:        " << int(interleaving) << std::endl;
  std::cout << "Compression:         " << compress << std::endl;
  std::cout << "*** END HEADER ***" << std::endl;
}

template<class A> double baseline( const A* wave )
{
  return baseline( wave, 0, pre_tr_smp/2 );
}

template<class A> double baseline( const A* wave, uint32_t start_smp, uint32_t n_smp )
{
  double bline = 0;
  double end_smp = start_smp + n_smp;
  for( uint32_t smp = 0; smp < n_smp; ++smp )
  {
    bline += wave[start_smp + smp];
  //  std::cout << " bline " << bline << " wave " << wave[smp] << " smp "<< smp << " n_smp " << n_smp << std::endl; 
  } 
  // std::cout << "-------------------------------------------------------Next----------------------------------------------------------" << std::endl;
  bline /= n_smp;
  return bline;
}

template<class A> double baseline_var( const A* wave, double bl )
{
  uint32_t n_smp = pre_tr_smp/2;
  double var = 0;
  for( uint32_t smp = 0; smp < n_smp; ++smp )
  {
    double tmp = (wave[smp] - bl);
    var += tmp*tmp;
  }
  var /= n_smp;
  return var;
}

template<class A> double fit_line( const A* a, double& k, double& m, int32_t xi, int32_t xf )
{
  if( xi < 0 ) xi = 0;
  if( xf > int32_t(smp_per_rec) ) xf = smp_per_rec;
  double N = xf - xi;
  double y_mean = 0;
  double xy_sum = 0;
  double x2_sum = 0;
  for( int32_t i = xi; i < xf; ++i )
  {
    y_mean += a[i];
    xy_sum += i*a[i];
    x2_sum += i*i;
  }
  y_mean /= N;
  double x_mean = 0.5 * ( xf + xi - 1 );

  k = ( xy_sum - N * x_mean * y_mean ) / ( x2_sum - N * x_mean * x_mean );
  m = y_mean - k * x_mean;

  double chi2_red = 0;
  for( int32_t i = xi; i < xf; ++i )
    chi2_red = std::pow( a[i] - k*i - m, 2 );
  chi2_red /= N;
  
  k /= ns_per_smp;
  return chi2_red;
}

template<class A> double check_pileup( const A* wave, double threshold, std::vector< std::pair<double, double> >& peaks, int32_t polarity )
{
  peaks.clear();
  polarity  = polarity < 0 ? -1 : 1;

  uint32_t smp = 0;
  while( smp < smp_per_rec )
  {
    while( smp < smp_per_rec && polarity*wave[smp] < polarity*threshold ) ++smp;
    // wave is now above threshold for positive polarity
    uint32_t start = smp++;
    if( start >= smp_per_rec ) break;
    while( smp < smp_per_rec && polarity*wave[smp] >= polarity*threshold ) ++smp;
    double pos;
    double value = polarity < 0 ? get_min( wave, pos, start, smp ) : get_max( wave, pos, start, smp );

    peaks.push_back( std::pair<double,double>(pos,value) );
  }
  return peaks.size();
}

template<class A> double get_fwhm( const A* wave, uint32_t peak_pos, double zero_level )
{
  double value = wave[peak_pos];
  double half = ( value + zero_level ) / 2;
  double wl = left_threshold( wave, peak_pos, half ); //left bound
  double wr = right_threshold( wave, peak_pos, half ); // right bound
  return (wr - wl) * ns_per_smp;
}

template<class A> double left_threshold( const A* wave, uint32_t pos, double limit )
{
  double value = wave[pos];
  //std::cout << "limit " << limit << std::endl;
  //std::cout << "value " << value << std::endl;
  int32_t i;
  for( i = pos; i >= 0; --i )
    if( ( value > limit && wave[i] < limit ) || ( value < limit && wave[i] > limit ) )
    //	std::cout << "wave pos " << wave[i] << std::endl;
      return linear_interpolation( i, i+1, wave[i], wave[i+1], limit );
  return 0;
}

template<class A> double right_threshold( const A* wave, uint32_t pos, double limit )
{
  double value = wave[pos];
  for( uint32_t i = pos; i < smp_per_rec; ++i )
    if( ( value > limit && wave[i] < limit ) || ( value < limit && wave[i] > limit ) )
      return linear_interpolation( i-1, i, wave[i-1], wave[i], limit );
  return smp_per_rec - 1;
}

template<class A> double get_min( const A* wave, double& pos)
{
  return get_min( wave, pos, 0, smp_per_rec );
}

template<class A> double get_min( const A* wave, double& pos, uint32_t start_bin, uint32_t end_bin )
{
  if( start_bin >= end_bin )
    throw "get_min(): Invalid interval";
  uint32_t min_pos = 0;
  end_bin = std::min( end_bin, smp_per_rec );
  A min_val = std::numeric_limits<A>::max(); // Initialize with the largest positive number; The use of std::numeric_limits<A>::min() would give the smallest positive number (that is, near to zero);


  for( uint32_t smp = start_bin; smp < end_bin; ++smp )
  {
// std::cout << "smp= " << smp << " wave[smp]= " << wave[smp] << std::endl;

    if( min_val > wave[smp] )
    {
      min_val = wave[smp];
      min_pos = smp;
    }
  }
  if( min_pos <= start_bin || min_pos >= end_bin - 1 )
  {
    pos = min_pos * ns_per_smp;
    return min_val;
  }
  double value = quadratic_extremum( wave[min_pos-1], min_val, wave[min_pos+1], pos );
  pos = (min_pos + pos) * ns_per_smp;
  return value;
}

// pos contains initial guess
template<class A> double get_extremum_polyfit( const A* wave, double& pos, uint32_t start_bin, uint32_t end_bin )
{
  if( start_bin >= end_bin )
    throw "get_min(): Invalid interval";
  uint32_t n_bins = end_bin - start_bin;
  TGraph points( n_bins );
  for( uint32_t i = 0; i < n_bins; ++i )
    points.SetPoint( i, (start_bin+i) * ns_per_smp, wave[start_bin + i] );
  static TF1 polynomial("PolyFit", "pol2", start_bin * ns_per_smp, end_bin * ns_per_smp );
  points.Fit( &polynomial, "Q" );
  double fit_pos = -polynomial.GetParameter(1) / ( 2 * polynomial.GetParameter(2) );
  if( fit_pos > start_bin * ns_per_smp && fit_pos < end_bin * ns_per_smp )
    pos = fit_pos;
  return polynomial.Eval(pos);
}

template<class A> double get_max( const A* wave, double& pos )
{
  return get_max( wave, pos, 0, smp_per_rec );
}

template<class A> double get_max( const A* wave, double& pos, uint32_t start_bin, uint32_t end_bin )
{
  if( start_bin >= end_bin )
    throw "get_max(): Invalid interval";
  uint32_t max_pos = 0;
  end_bin = std::min( end_bin, smp_per_rec );
  A max_val = -std::numeric_limits<A>::max(); // Initialize with the largest negative number; The use of std::numeric_limits<A>::min() would give the smallest positive number (that is, near to zero);

  for( uint32_t smp = start_bin; smp < end_bin; ++smp )
  {
//std::cout << "smp= " << smp << " wave[smp]= " << wave[smp] << " temp.max= " << max_val << std::endl;
    if( max_val < wave[smp] )
    {
      max_val = wave[smp];
      max_pos = smp;
//std::cout << "max_val= " << max_val << " max_pos= " << smp << std::endl;
    }
  }
  if( max_pos <= start_bin || max_pos >= end_bin - 1 )
  {
    pos = max_pos * ns_per_smp;
    return max_val;
  }
  double value = quadratic_extremum( wave[max_pos-1], max_val, wave[max_pos+1], pos );
  pos = (max_pos + pos) * ns_per_smp;
  return value;
}

template<class A> double threshold( const A* wave, double limit )
{
  return right_threshold( wave, 0, limit );
}

double linear_interpolation( double x1, double x2, double y1, double y2, double y )
{
  double delta_x = ( y - y1 ) * ( x2 - x1 ) / ( y2 - y1 );
 /* std::cout << "delta x " << delta_x << std::endl;
  std::cout << "y " << y << std::endl;
  std::cout << "y1 " << y1 << std::endl;
  std::cout << "y2 " << y2 << std::endl;
  std::cout << "x1 " << x1 << std::endl;
  std::cout << "x2 " << x2 << std::endl;
  std::cout << " " << std::endl;*/
  return x1 + delta_x;
}

// Assumes equidistant points and returns difference from centre value
double quadratic_extremum( double y1, double y2, double y3, double& delta_x )
{
  delta_x = ( y3 - y1 ) / ( 2 * ( y1 + y3 - 2 * y2 ) );
  return (y1 + y3 - 2*y2)/2 * delta_x*delta_x + (y3 - y1)/2 * delta_x + y2;
}

template<class A> void store( const A* waveform, const TString name )
{
  TGraph* wf_graph = new TGraph(smp_per_rec);
  for( unsigned int smp = 0; smp < smp_per_rec; ++smp )
    wf_graph->SetPoint( smp, smp*ns_per_smp, waveform[smp] );
  wf_graph->SetName( name );
  wf_graph->SetTitle( name );

  RootFileManager* rm = RootFileManager::GetInstance();
  rm->Adopt( wf_graph );
}

void store( const fftw_complex* waveform, const TString name )
{
  TGraph* wf_graph_r = new TGraph(smp_per_rec);
  TGraph* wf_graph_i = new TGraph(smp_per_rec);
  for( unsigned int smp = 0; smp < smp_per_rec; ++smp )
  {
    wf_graph_r->SetPoint( smp, smp*ns_per_smp, waveform[smp][0] );
    wf_graph_i->SetPoint( smp, smp*ns_per_smp, waveform[smp][1] );
  }
  wf_graph_r->SetName(  name+"_real" );
  wf_graph_i->SetName(  name+"_imag" );
  wf_graph_r->SetTitle( name+"_real" );
  wf_graph_i->SetTitle( name+"_imag" );

  RootFileManager* rm = RootFileManager::GetInstance();
  rm->Adopt( wf_graph_r );
  rm->Adopt( wf_graph_i );
}

template<class A> void savitzky_golay_quad_smooth( const A* in, double* out, uint32_t width )
{
  const uint32_t reach = ( width - 1 )/2;
  switch( width )
  {
    case 11:
      for( uint32_t i = 0; i < reach; ++i )
      {
        out[i] = in[i];
        out[smp_per_rec - 1 - i] = in[smp_per_rec - 1 - i];
      }
      for( uint32_t i = reach; i < smp_per_rec - reach; ++i )
      {
        out[i] = - 36 * ( in[i-5] + in[i+5] ) 
                 + 9  * ( in[i-4] + in[i+4] )
                 + 44 * ( in[i-3] + in[i+3] ) 
                 + 69 * ( in[i-2] + in[i+2] ) 
                 + 84 * ( in[i-1] + in[i+1] ) 
                 + 89 * in[i];
        out[i] /= 429;
      }
      break;
    case 9:
      for( uint32_t i = 0; i < reach; ++i )
      {
        out[i] = in[i];
        out[smp_per_rec - 1 - i] = in[smp_per_rec - 1 - i];
      }
      for( uint32_t i = reach; i < smp_per_rec - reach; ++i )
      {
        out[i] = - 21 * ( in[i-4] + in[i+4] )
                 + 14 * ( in[i-3] + in[i+3] )
                 + 39 * ( in[i-2] + in [i+2] )
                 + 54 * ( in[i-1] + in[i+1] ) 
                 + 59 * in[i];
        out[i] /= 231;
      }
      break;
    case 7:
      for( uint32_t i = 0; i < reach; ++i )
      {
        out[i] = in[i];
        out[smp_per_rec - 1 - i] = in[smp_per_rec - 1 - i];
      }
      for( uint32_t i = reach; i < smp_per_rec - reach; ++i )
      {
        out[i] = - 2 * ( in[i-3] + in[i+3] )
                 + 3 * ( in[i-2] + in [i+2] ) 
                 + 6 * ( in[i-1] + in[i+1] ) 
                 + 7 * in[i];
        out[i] /= 21;
      }
      break;
    default:
      std::cerr << "Savitzky-Golay Quadractic smoothing is not implemented for a smoothing width of " <<
        width << " defaulting to 5." << std::endl;
    case 5:
      for( uint32_t i = 0; i < reach; ++i )
      {
        out[i] = in[i];
        out[smp_per_rec - 1 - i] = in[smp_per_rec - 1 - i];
      }
      for( uint32_t i = reach; i < smp_per_rec - reach; ++i )
      {
        out[i] = - 3  * ( in[i-2] + in [i+2] ) 
                 + 12 * ( in[i-1] + in[i+1] ) 
                 + 17 * in[i];
        out[i] /= 35;
      }
      break;
  }
}

template<class A> void rectangle_smooth( const A* in, double* out, uint32_t width )
{
  const uint32_t reach = ( width - 1 )/2;
  for( uint32_t i = reach; i < smp_per_rec - reach; ++i )
  {
    double sum = 0;
    for( uint32_t j = i - reach; j < i + reach + 1; j++ )
      sum += in[j];
    out[i] = sum/(2*reach+1);
  }
  for( uint32_t i = 0; i < reach; ++i )
  {
    out[i] = in[i];
    out[smp_per_rec - 1 - i] = in[smp_per_rec - 1 - i];
  }
}

template<class A> void triangle_smooth( const A* in, double* out, uint32_t width )
{
  const int32_t reach = ( width - 1 )/2;
  const double norm = (reach+1)*(reach+1);
  for( uint32_t i = reach; i < smp_per_rec - reach; ++i )
  {
    out[i] = 0;
    for( int32_t j = -reach; j <= reach; ++j )
      out[i] += in[i-j] * (reach + 1 - TMath::Abs(j) );
    out[i] /= norm;
  }
  for( int32_t i = 0; i < reach; ++i )
  {
    out[i] = in[i];
    out[smp_per_rec - 1 - i] = in[smp_per_rec - 1 - i];
  }
}

template<class A> void subtract( const A* in1, const double* in2, double* out )
{
  for( unsigned int i = 0; i < smp_per_rec; ++i )
    out[i] = in1[i] - in2[i];
}

template<class A> void diff( const A* in, double* out )
{
  for( unsigned int i = 0; i < smp_per_rec - 1; ++i )
    out[i] = ( in[i+1] - in[i] ) / ns_per_smp;
}

fftw_complex* lowpass( fftw_complex* freq, double high_freq )
{
  /* Removes any frequencies higher than high_freq */
  return bandpass( freq, 0, high_freq );
}

fftw_complex* highpass( fftw_complex* freq, double low_freq )
{
  /* Removes any frequencies lower than low_freq */
  return bandpass( freq, low_freq, smp_per_rec );
}

fftw_complex* bandpass( fftw_complex* freq, double low_freq, double high_freq )
{
  uint32_t  low_k =  low_freq * smp_per_rec * ns_per_smp*1e-9 + 0.5;
  uint32_t high_k = high_freq * smp_per_rec * ns_per_smp*1e-9 + 0.5;
  /* Removes any frequencies lower than high_freq
   * but higher than low_frew */
  for( uint32_t smp = 0; smp < std::max(low_k, uint32_t(0)); ++smp )
    freq[smp][0] = freq[smp_per_rec - smp - 1][0] =
      freq[smp][1] = freq[smp_per_rec - smp - 1][1] = 0;
  for( uint32_t smp = high_k + 1; smp < smp_per_rec; ++smp )
    freq[smp][0] = freq[smp_per_rec - smp - 1][0] =
      freq[smp][1] = freq[smp_per_rec - smp - 1][1] = 0;
  return freq;
}

fftw_complex* bandfilter( fftw_complex* freq, double low_freq, double high_freq )
{
  uint32_t  low_k =  low_freq * smp_per_rec * ns_per_smp*1e-9 + 0.5;
  uint32_t high_k = high_freq * smp_per_rec * ns_per_smp*1e-9 + 0.5;
  /* Removes any frequencies between low_freq and high_freq */
  for( uint32_t smp = std::max(uint32_t(0),low_k);
      smp < std::min(high_k+1, smp_per_rec); ++smp )
    freq[smp][0] = freq[smp_per_rec - smp - 1][0] =
      freq[smp][1] = freq[smp_per_rec - smp - 1][1] = 0;
  return freq;
}


// time_cfd_raw = CFD( raw[ch], bipolar_raw[ch], 0.3, 1.,  -1800., 0);   // DT 2019/10/25. To compare two ways of doing CFD.
//  time_cfd_smooth = CFD( smooth[ch], bipolar_smooth[ch], 0.3, 1.,  -1800., 0);  // DT 2019/10/25. To compare two ways of doing CFD.




// DT 2019/10/25. CFD function as written by BB;
// 
// *************
//
// CFD function for both positive (pol = 1) and negative (pol = 0) signals
template<class A> double CFD(const A* wave, double* out, double frac, uint32_t delay, double Threshold, int pol)
//  time_cfd_raw = CFD( raw[ch], bipolar_raw[ch], 0.3, 1.,  -1800., 0);  
//  time_cfd_smooth = CFD( smooth[ch], bipolar_smooth[ch], 0.3, 1.,  -1800., 0);  
 {
  double trigger = 0;
  double pos = 0;
  double val_max;
  double val_min;

  if (pol == 1) {

    val_max =  get_max( wave, pos );
 //   std::cout << "val max " << val_max << std::endl;
  //  if( val_max < Threshold) {
  //        std::cout<< ": CFD below threshold"<< std::endl;
  //    return 2;       }

  // Bipolar pulse 
    for( uint32_t i = 0; i < smp_per_rec; ++i ) {
    if( i <= delay ) { 
      out[i] =  (-1)*frac*(wave[i]);
               }
    else {
      out[i] = wave[i-delay] + (-1)*frac*wave[i];
      } 
   // std::cout << "out[i] " << out[i] << std::endl;
    }

  // Find zero crossing
    uint32_t t = delay;
 //   std::cout << "t " << t << std::endl;
    while( ( out[t] < Threshold ) && ( t < smp_per_rec ) ) {
      val_max = 0.;
      t++;
 //     std::cout << "t aumentado " << t << std::endl;
  		}

    for ( uint32_t i = 0; i < smp_per_rec; i++) {
      if (val_max  < out[i]) {
       //std::cout<< "wl[" << i << "] = " << w1[i]<< std::endl;
        val_max = out[i];
        t = i;
        //std::cout<< t << " T first if" << std::endl;   
      }
       // std::cout<< t << " t del max " << std::endl;
    //    std::cout<< " w1[" << i << "] = " << w1[i] << std::endl;    
    }
    while( (out[t] > 0.) && (t > delay) )  {
      trigger = 0.;
      if(out[t + 1] == out[t]) {
        trigger = 0.5*(t + (t + 1));
        //std::cout<< "trigger first case:  " << trigger  <<std::endl;   
         }
      else {
        trigger = t - (t + 1 - t) * out[t] / (out[t + 1] - out[t]);            
         // trigger = linear_interpolation (t-1, t, w1[t-1], w1[t], w1[t+1]); } 
        }
      t = t-1;                 
    }
    }
  else if (pol == 0) {
      val_min =  get_min(wave, pos);
  //  std::cout << "val min " << val_min << std::endl;
  //  if( val_min > Threshold) {
  //        std::cout<< ": CFD below threshold"<< std::endl;
  //    return -2;       }
  // Bipolar pulse 
    for(uint32_t i = 0; i < smp_per_rec; i++) {
    	if( i <= delay ) {
      	out[i] =  (-1)*frac*(wave[i]);
    	    			 }
    	else {
      	out[i] = wave[i-delay] + (-1)*frac*wave[i] ;
     		 } 
    	}
  // Find zero crossing
    uint32_t t = delay;
    while( (out[t] > Threshold) && (t < smp_per_rec) ) {
      val_min = 0.;
      t++;                                            }

    for (uint32_t i = 0; i < smp_per_rec; i++) {
      if (val_min  > out[i]) {
         //std::cout<< "wl[" << i << "] = " << w1[i]<< std::endl;
        val_min = out[i];
        t = i;
       //    std::cout<< t << " T primo if" << std::endl;   
      }
  //    std::cout<< t << " t del max " << std::endl;
  //    std::cout<< " w1[" << i << "] = " << w1[i] << std::endl;    
    }
    while( (out[t] < 0.) && (t > delay) )  {
      trigger = 0.;
      if(out[t+1] == out[t]) {
        trigger = 0.5*(t + (t+1));
        //std::cout<< "trigger primo caso:  " << trigger  <<std::endl;  
          }
      else {
        trigger= t - (t + 1 - t)* out[t]/ (out[t + 1] - out [t]);            
         // trigger = linear_interpolation (t-1, t, w1[t-1], w1[t], w1[t+1]); }
        }
      t = t-1;                 
    }
    }

  //std::cout<< "trigger:  " << trigger  <<std::endl;
  return trigger;
  //std::cout << "-------------------------------------------------------Next----------------------------------------------------------" << std::endl;
}



// 
// 
// *************
//
// DT 2019/10/25. END of CFD function as written by BB;


















/* * * * * * * * * * * * * * * *
 * Fourier transform functions *
 * * * * * * * * * * * * * * * */ 
void fftw_init()
{
  fftw_forward_plan = fftw_plan_dft_1d(smp_per_rec, waveform_t, waveform_f, FFTW_FORWARD, FFTW_MEASURE);
  fftw_backward_plan = fftw_plan_dft_1d(smp_per_rec, waveform_f, waveform_t, FFTW_BACKWARD, FFTW_MEASURE);
}

fftw_complex* to_freq( const int16_t* waveform )
{
  /* Copy data into complex array */
  for( unsigned int i = 0; i < smp_per_rec; ++i)
  {
    waveform_t[i][0] = waveform[i];
    waveform_t[i][1] = 0;
  }
  return fftw_forward(waveform_t);
}

fftw_complex* to_freq( const double* waveform )
{
  /* Copy data into complex array */
  for( unsigned int i = 0; i < smp_per_rec; ++i)
  {
    waveform_t[i][0] = waveform[i];
    waveform_t[i][1] = 0;
  }
  return fftw_forward(waveform_t);
}

fftw_complex* to_freq( const fftw_complex* waveform )
{
  if( waveform != waveform_t )
    /* Copy data into complex array */
    for( unsigned int i = 0; i < smp_per_rec; ++i)
    {
      waveform_t[i][0] = waveform[i][0];
      waveform_t[i][1] = waveform[i][1];
    }
  return fftw_forward(waveform_t);
}

fftw_complex* to_time( const int16_t* waveform )
{
  /* Copy data into complex array */
  for( unsigned int i = 0; i < smp_per_rec; ++i)
    waveform_f[i][0] = waveform[i];
  return fftw_backward(waveform_f);
}

fftw_complex* to_time( const fftw_complex* waveform )
{
  if( waveform != waveform_f )
    /* Copy data into complex array */
    for( unsigned int i = 0; i < smp_per_rec; ++i)
    {
      waveform_f[i][0] = waveform[i][0];
      waveform_f[i][1] = waveform[i][1];
    }
  return fftw_backward(waveform_f);
}

fftw_complex* fftw_forward( fftw_complex* waveform)
{
  assert( waveform == waveform_t);
  for( unsigned int i = 0; i < smp_per_rec; ++i)
  {
    waveform_f[i][0] = 0;
    waveform_f[i][1] = 0;
  }
  fftw_execute(fftw_forward_plan);
  return waveform_f;
}

fftw_complex* fftw_backward( fftw_complex* waveform)
{
  assert( waveform == waveform_f);
  for( unsigned int i = 0; i < smp_per_rec; ++i)
  {
    waveform_t[i][0] = 0;
    waveform_t[i][1] = 0;
  }
  fftw_execute(fftw_backward_plan);
  /* Normalise */
  for( unsigned int i = 0; i < smp_per_rec; ++i)
  {
    waveform_t[i][0] = waveform_t[i][0]/smp_per_rec;
    waveform_t[i][1] = waveform_t[i][1]/smp_per_rec;
  }
  return waveform_t;
}



template<class A>  void moving_average(const A* in, double* out, int widthofmoving)
{
    
    /*  Moving average for numer of points defined in "width"  for whole waveform
    
    */
  int  n;
  int WindowSum = 0;
  int size = smp_per_rec;
  int width = widthofmoving;
  // where do we start smoothing no point take all pre trigger sample

      //
      // presumming so later i would not need to call loop again
      //
      
      for (int i = 0; i < n; ++i)
      {
          WindowSum = WindowSum + in[i];
      }

      for (int i=0; i<smp_per_rec; ++i)
    {
         if (i >= size - n){
            WindowSum  = WindowSum - in[i-n-1];
            width = width -1;
            out[i] = ( WindowSum ) / double(width); 
            }
         else if (i < n){
            WindowSum  = WindowSum + in[i+n];
            out[i] = ( WindowSum ) / double( n + i+1);}
         else{
            WindowSum  = WindowSum + in[i+n] - in[i-n-1];
            out[i] = ( WindowSum ) / double(width);}
    }
    
}


// Diego 2018/10/15.
// To calculate the integral of a signal;
template<class A> double integratePeak(const A* in, uint32_t low, uint32_t high, double bl_value)
{
  double sum=0;
  double integral=0;
  for( uint32_t i = low; i <= high; i++ )
  sum=sum+(in[i]-bl_value);

  integral=TMath::Abs(sum*ns_per_smp);

//  std::cout << "integral= " << integral << std::endl;
  return integral;
}



//
// Special filters for pileup separation 
//

template<class A>  void  Trapezoidal_Filter(const A* in, double* out, int k, int l)
{ /*  This helps to separate the pileups and determine the amplitude
    
 Nuclear Instruments and Methods in Physics Research A 345 (1994) 337-345  */
    double p[2] = {0,0};
    double d =0;
    /*  Testing whether this filter can be applied with given parameters */
     if (l < k)
         throw std::runtime_error(" error: l < k in trapezoidal filter ");
        


     double M = 1 / (exp(5.0e-5 ) - 1 );
         
     for (int n=0; n<smp_per_rec; ++n)
        {

         if (n - l - k >= 0)
             d = (in[n] - in[n-k] - in[n-l]  + in[n-k-l]);
         
         else if (n- l >= 0)
             d = (in[n]   - in[n-k] - in[n-l]);
         
         else if (n-k >= 0 )
             d =( in[n] - in[n-k]);
         
         else if (n>=0 )
             d = (in[n] );
        
        
         p[1]  =  p[0] + d;
         out[n]  =  out[n-1] + M*d + p[1];
         p[0] = p[1];
        }

    
    
}

            
            
template<class A>  void MWD(const A* in, double* out, int width)
{ /* This is moving average deconvolution   .  can be used to seperate the pile ups and determine pulse height
     taken from reference : " Extraction of energy and time from pile-up pulses
        with fast sampling ADC analysis techniques " Anton Roth,Nuclear Physics Division, June 16, 2016*/

  int ii;
  double Pre_amp_dec = 0.000065;         // Specific For Preamp. Decay constant
  int WindowSum = 0;
  
     
      for (int i=0; i<smp_per_rec; ++i)
    {
     ii = i - width;
     
     
     if (ii < 0){
        out[i] = in[i] + Pre_amp_dec* WindowSum/width;
        WindowSum  = WindowSum  + in[i];
                }
    else {
        out[i] = in[i] - in[ii] + Pre_amp_dec* WindowSum/width;
        WindowSum  = WindowSum - in[ii] + in[i];
         }
     
     
    }

}

