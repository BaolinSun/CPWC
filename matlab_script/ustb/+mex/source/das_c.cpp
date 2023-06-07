#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <vector>
#include <concurrent_vector.h>
#include <tuple>
#include <valarray>
#include <set>
#if defined(_WIN_)
#include <ppl.h>           // Requires VS2010+
#elif defined (_UNIX_)
#include <parallel_for.h>  // Requires Intel tbb
#endif

// compulsory input
#define	M_P         prhs[0]	// channel_data [time, channel, frame]
#define	M_FS		prhs[1] // sampling frequency (Hz)
#define M_T0		prhs[2]	// initial time (s)

#define	M_APO_TX	prhs[3]	// transmit apodization [pixel, wave]
#define	M_APO_RX	prhs[4]	// receive apodization [pixel, channel]

#define	M_DELAY_TX  prhs[5]	// transmit delay [pixel, wave]
#define	M_DELAY_RX  prhs[6]	// receive delay [pixel, channel]

#define	M_FD		prhs[7] // modulation frequency (Hz)
#define	M_SUM		prhs[8] // sum mode 0 -> NONE, 1->RX, 2->TX, 3->BOTH

// optional input
#define	M_VERBOSE	prhs[9] // verbose flag [Optional]

// output
#define	M_D			plhs[0] // delayed data [pixel, channel, frame]

// constants
#define EPS 1e-6
#define PI 3.14159265359

// constants
#define NONE 0
#define RX 1
#define TX 2
#define BOTH 3

#define VERSION "1.1.2"

// types
typedef std::vector<float*> vec_p_float; // change to single

//-----------------------------------------------------------------------------

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    ///////////////////////////////
    // CHECKING ARGUMENTS
    ///////////////////////////////////////
    // number of inputs
    if (nrhs<9) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:nrhs", "Too few input arguments");
    if (nrhs>10) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:nrhs", "Too many input arguments");
    if (nlhs>1) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:nlhs", "Too many output arguments");
    
    ///////////////////////////////////////
    // VERBOSE
    ///////////////////////////////////////
    bool verbose = false;
    if (nrhs == 8)  verbose = (*((int*)mxGetData(M_VERBOSE)) > 0);
    
    if (verbose) {
        mexPrintf("---------------------------------------------------------------\n");
        mexPrintf(" USTB General beamformer\n");
        mexPrintf("---------------------------------------------------------------\n");
        mexPrintf(" Single precision\n");
        mexPrintf(" Vers:  %s\n",VERSION);
        mexPrintf(" Auth:  Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>\n");
        mexPrintf(" Date:  2017/10/05\n");
        mexPrintf("---------------------------------------------------------------\n");
        mexEvalString("drawnow;");
    }
    
    ///////////////////////////////////////
    // SUM (NONE/RX/TX/BOTH)
    ///////////////////////////////////////
    int sum_mode = (*((int*)mxGetData(M_SUM)));
    if((sum_mode<0)||(sum_mode>3)) mexErrMsgTxt("Unknown sum mode. Available: 0 -> NONE, 1->RX, 2->TX, 3->BOTH");
    
    ///////////////////////////////////////
    // CHANNEL DATA
    ///////////////////////////////////////
    // dimensions
    int ndim = (int)mxGetNumberOfDimensions(M_P);
    if (ndim<2 || ndim>4) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Unknown channel data format. Expected from 2 to 4 dimensions: [time, channel, wave, frame]");
    // size
    const mwSize* p_dim = mxGetDimensions(M_P);
    int L = (int)p_dim[0];                              // number of time samples
    int N = (int)p_dim[1];                              // number of channels
    int W = (int)1; if (ndim > 2) W = (int)p_dim[2];	// number of waves
    int F = (int)1;	if (ndim > 3) F = (int)p_dim[3];	// number of frames
    // single or double precision
    if (mxIsDouble(M_P)) mexErrMsgTxt("The channel data should be single precision");
    // complex or real
    bool complex_data=mxIsComplex(M_P);
    
    if (verbose) {
        if (complex_data) mexPrintf("Data Type                       Complex\n");
        else mexPrintf("Data Type                    Real\n");
        switch (sum_mode) {
            case NONE:
                mexPrintf("Sum                             None\n");
                break;
            case RX:
                mexPrintf("Sum                             Channels\n");
                break;
            case TX:
                mexPrintf("Sum                             Waves\n");
                break;
            case BOTH:
                mexPrintf("Sum                             Both channels and waves\n");
                break;
        }
        mexPrintf("Time Samples                    %i\n", L);
        mexPrintf("Channels						%i\n", N);
        mexPrintf("Waves                           %i\n", W);
        mexPrintf("Frames							%i\n", F);
        mexEvalString("drawnow;");
    }
    
    ///////////////////////////////////////
    // DELAY TX
    ///////////////////////////////////////
    // dimensions
    ndim = (int)mxGetNumberOfDimensions(M_DELAY_TX);
    if (ndim>2) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Unknown transmit delay format. Expected 2 dimensions: [pixel, wave]");
    // size
    p_dim = mxGetDimensions(M_DELAY_TX);
    int P = (int)p_dim[0];                  // number of pixels
    int W_check = (int)1;
    if (ndim == 2) W_check = (int)p_dim[1];	// number of waves
    if (W!=W_check) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Number of waves in channel data & transmit delays do not match.");
    // single or double precision
    if (mxIsDouble(M_DELAY_TX)) mexErrMsgTxt("The transmiy delays should be single precision");
    
    ///////////////////////////////////////
    // DELAY RX
    ///////////////////////////////////////
    // dimensions
    ndim = (int)mxGetNumberOfDimensions(M_DELAY_RX);
    if (ndim>2) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Unknown receive delay format. Expected 2 dimensions: [pixel, channel]");
    // size
    p_dim = mxGetDimensions(M_DELAY_RX);
    int P_check = (int)p_dim[0];                  // number of pixels
    if (P!=P_check) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Number of pixels in transmit and receive delays do not match.");
    int N_check = (int)1;
    if (ndim == 2) N_check = (int)p_dim[1];       // number of channels
    if (N!=N_check) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Number of channels in channel data & receive delays do not match.");
    // single or double precision
    if (mxIsDouble(M_DELAY_RX)) mexErrMsgTxt("The transmiy delays should be single precision");
    
    if (verbose) {
        mexPrintf("Pixels							%i\n", P);
        mexEvalString("drawnow;");
    }
    
    ///////////////////////////////////////
    // APO TX
    ///////////////////////////////////////
    // dimensions
    ndim = (int)mxGetNumberOfDimensions(M_APO_TX);
    if (ndim>2) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Unknown transmit apodization format. Expected 2 dimensions: [pixel, wave]");
    // size
    p_dim = mxGetDimensions(M_APO_TX);
    P_check = (int)p_dim[0];                  // number of pixels
    if (P!=P_check) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Number of pixels in transmit apodization and transmit delays do not match.");
    W_check = (int)1;
    if (ndim == 2) W_check = (int)p_dim[1];	// number of waves
    if (W!=W_check) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Number of waves in channel data & transmit apodization do not match.");
    // single or double precision
    if (mxIsDouble(M_APO_TX)) mexErrMsgTxt("The transmiy apodization should be single precision");
    
    ///////////////////////////////////////
    // APO RX
    ///////////////////////////////////////
    // dimensions
    ndim = (int)mxGetNumberOfDimensions(M_APO_RX);
    if (ndim>2) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Unknown receive apodization format. Expected 2 dimensions: [pixel, channel]");
    // size
    p_dim = mxGetDimensions(M_APO_RX);
    P_check = (int)p_dim[0];                  // number of pixels
    if (P!=P_check) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Number of pixels in receive apodization and receive delays do not match.");
    N_check = (int)1;
    if (ndim == 2) N_check = (int)p_dim[1];       // number of channels
    if (N!=N_check) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Number of channels in channel data & receive apodization do not match.");
    // single or double precision
    if (mxIsDouble(M_APO_RX)) mexErrMsgTxt("The transmiy delays should be single precision");
    
    ///////////////////////////////////////
    // SAMPLING FREQUENCY
    ///////////////////////////////////////
    // dimensions
    const size_t len_Fs = mxGetNumberOfElements(M_FS);
    if (len_Fs != 1) mexErrMsgTxt("The sampling frequency should be an escalar");
    // single or double precision
    if (mxIsDouble(M_FS)) mexErrMsgTxt("The sampling frequency should be single precision");
    // read data
    const float Fs = *((float*)mxGetData(M_FS));
    const float dt = 1 / Fs;
    if (verbose) mexPrintf("Sampling frequency				%0.2f MHz\n", Fs / 1e6);
    
    ///////////////////////////////////////
    // INITIAL TIME
    ///////////////////////////////////////
    // dimensions
    const size_t len_t0 = mxGetNumberOfElements(M_T0);
    if (len_t0 != 1) mexErrMsgTxt("The initial time should be an escalar");
    // single or double precision
    if (mxIsDouble(M_T0)) mexErrMsgTxt("The initial time should be single precision");
    // read data
    const float t0 = *((float*)mxGetData(M_T0));
    
    if (verbose) mexPrintf("Initial time					%0.2f us\n", t0*1e6);
    
    ///////////////////////////////////////
    // MODULATION FREQUENCY
    ///////////////////////////////////////
    // dimension
    const size_t len_fd = mxGetNumberOfElements(M_FD);
    if (len_fd != 1) mexErrMsgTxt("Modulation frequency 'fd' should be an escalar");
    // single or double precision
    if (mxIsDouble(M_FD)) mexErrMsgTxt("The modulation frequency should be single precision");
    // read data
    const float fd = *((float*)mxGetData(M_FD));
    
    float wd = 0;
    bool IQ_version = false;
    if (std::abs(fd) > EPS) {
        if (!complex_data) mexErrMsgTxt("The modulation frequency > 0 but the input data is real. Check inputs.");
        wd =  2 * (float)PI * fd;
        IQ_version = true;
    }
    
    if (verbose) {
        mexPrintf ("Modulation frequency:           %0.2f MHz\n", fd / 1e6);
        if(IQ_version) mexPrintf("IQ data:                        true\n");
        else mexPrintf("IQ data:                        false\n");
        mexPrintf("---------------------------------------------------------------\n");
        mexEvalString("drawnow;");
    }
    
    ///////////////////////////////////////
    // OUTPUT MATRIX
    ///////////////////////////////////////
    mwSize out_size2[4];
    out_size2[0] = P;  // pixels
    out_size2[3] = F;  // frames
    
    switch(sum_mode) {
        case NONE:
            out_size2[1] = N;  // channels
            out_size2[2] = W;  // waves
            if (verbose) {
                mexPrintf("Output data size:               %d x %d x %d x %d\n", P,N,W,F);
                mexPrintf("---------------------------------------------------------------\n");
            }
            break;
        case RX:
            out_size2[1] = 1;  // channels
            out_size2[2] = W;  // waves
            if (verbose) {
                mexPrintf("Output data size:               %d x %d x %d x %d\n", P,1,W,F);
                mexPrintf("---------------------------------------------------------------\n");
            }
            break;
        case TX:
            out_size2[1] = N;  // channels
            out_size2[2] = 1;  // waves
            if (verbose) {
                mexPrintf("Output data size:               %d x %d x %d x %d\n", P,N,1,F);
                mexPrintf("---------------------------------------------------------------\n");
            }
            break;
        case BOTH:
            out_size2[1] = 1;  // channels
            out_size2[2] = 1;  // waves
            if (verbose) {
                mexPrintf("Output data size:               %d x %d x %d x %d\n", P,1,1,F);
                mexPrintf("---------------------------------------------------------------\n");
            }
            break;
    }
    // define array
    M_D = mxCreateNumericArray(4, (const mwSize*)&out_size2, mxSINGLE_CLASS, mxCOMPLEX);
    float* Dr = (float*)mxGetData(M_D);
    float* Di = (float*)mxGetImagData(M_D);
    
    // size variables
    const unsigned int LN=L*N;
    const unsigned int LNW=LN*W;
    const unsigned int PN=P*N;
    const unsigned int PW=P*W;
    const unsigned int PNW=PN*W;
    const unsigned int PF=P*F;
    const unsigned int PFN=PF*N;
    
    // Pointers to channel data (real & imaginary part)
    float* Pr = (float*)mxGetData(M_P);                     // real part of channel data
    float* Pi = NULL;
    if (complex_data) Pi = (float*)mxGetImagData(M_P);		// imaginary part of channel data
    
    // Pointers to delay & apodization matrices
    float* p_tx_delay  = (float*)mxGetData(M_DELAY_TX);
    float* p_tx_apo  = (float*)mxGetData(M_APO_TX);
    float* p_rx_delay  = (float*)mxGetData(M_DELAY_RX);
    float* p_rx_apo  = (float*)mxGetData(M_APO_RX);
    
    //////////////////////////////////////////////////////
    // Beamforming loop
    mexPrintf("USTB General beamformer MEX v%s ........", VERSION);
    mexEvalString("drawnow;");
    
#if defined (_WIN_)
    Concurrency::parallel_for(0, P, [&](int pp) { // pixel loop -> WIN
#elif defined(_UNIX_)
        tbb::strict_ppl::parallel_for(0, P, [&](int pp) { // pixel loop -> UNIX
#endif
            
            for (int rx = 0; rx<N; rx++) { // channel loop
                // references
                float& c_rx_delay = *(p_rx_delay + P*rx + pp);      // delay
                float& c_rx_apo = *(p_rx_apo + P*rx + pp);          // apodization
                
                if (c_rx_apo>1e-3) for (int w = 0; w<W; w++) { // wave loop
                    // references
                    float& c_tx_delay = *(p_tx_delay + P*w + pp);   // delay
                    float& c_tx_apo = *(p_tx_apo + P*w + pp);       // apodization
                    
                    if (c_tx_apo>1e-3) { // we skip any apodization value smaller than -60 dB
                    
                        float c_delay = c_tx_delay + c_rx_delay ;   // delay
                        float c_apo = c_tx_apo * c_rx_apo;          // apodization
                        
                        float denay = (c_delay - t0)*Fs;		// untruncated sample number
                        int n0 = (int)floor(denay);             // truncated sample number
                        float b = denay - n0;                   // linear interpolation coefficient 2
                        float a = 1 - b;                        // linear interpolation coefficient 1
                        
                        if (n0>0 && n0<(L - 1)) {
                            if (IQ_version) { // data must be complex
                                float phase = wd*c_delay;
                                float coswt = cos(phase);
                                float sinwt = sin(phase);
                                for (int f = 0; f<F; f++) { // frame loop
                                    // reference to data
                                    float& prn0=*(Pr + L*rx + LN*w + LNW*f + n0);
                                    float& prn1=*(Pr + L*rx + LN*w + LNW*f + n0 + 1);
                                    float& pin0=*(Pi + L*rx + LN*w + LNW*f + n0);
                                    float& pin1=*(Pi + L*rx + LN*w + LNW*f + n0 + 1);
                                    
                                    // fractional part of the delay -> linear interpolation
                                    float re = a*prn0 + b*prn1;
                                    float im = a*pin0 + b*pin1;
                                    
                                    // apply phase change and delay
                                    switch(sum_mode) {
                                        case NONE:
                                            Dr[pp + P*rx + PN*w + PNW*f] = c_apo*(re*coswt - im*sinwt);
                                            Di[pp + P*rx + PN*w + PNW*f] = c_apo*(im*coswt + re*sinwt);
                                            break;
                                        case RX:
                                            Dr[pp + P*w + PW*f] += c_apo*(re*coswt - im*sinwt);
                                            Di[pp + P*w + PW*f] += c_apo*(im*coswt + re*sinwt);
                                            break;
                                        case TX:
                                            Dr[pp + P*rx + PN*f] += c_apo*(re*coswt - im*sinwt);
                                            Di[pp + P*rx + PN*f] += c_apo*(im*coswt + re*sinwt);
                                            break;
                                        case BOTH:
                                            Dr[pp + P*f] += c_apo*(re*coswt - im*sinwt);
                                            Di[pp + P*f] += c_apo*(im*coswt + re*sinwt);
                                            break;
                                    }
                                }
                            } else if (complex_data){ // no modulation but complex data -> hilbert transformed RF data
                                for (int f = 0; f<F; f++) { // frame loop
                                    // reference to data
                                    float& prn0=*(Pr + L*rx + LN*w + LNW*f + n0);
                                    float& prn1=*(Pr + L*rx + LN*w + LNW*f + n0 + 1);
                                    float& pin0=*(Pi + L*rx + LN*w + LNW*f + n0);
                                    float& pin1=*(Pi + L*rx + LN*w + LNW*f + n0 + 1);
                                    
                                    switch(sum_mode) {
                                        case NONE:
                                            Dr[pp + P*rx + PN*w + PNW*f] = c_apo*(a*prn0 + b*prn1);
                                            Di[pp + P*rx + PN*w + PNW*f] = c_apo*(a*pin0 + b*pin1);
                                            break;
                                        case RX:
                                            Dr[pp + P*w + PW*f] += c_apo*(a*prn0 + b*prn1);
                                            Di[pp + P*w + PW*f] += c_apo*(a*pin0 + b*pin1);
                                            break;
                                        case TX:
                                            Dr[pp + P*rx + PN*f] += c_apo*(a*prn0 + b*prn1);
                                            Di[pp + P*rx + PN*f] += c_apo*(a*pin0 + b*pin1);
                                            break;
                                        case BOTH:
                                            Dr[pp + P*f] += c_apo*(a*prn0 + b*prn1);
                                            Di[pp + P*f] += c_apo*(a*pin0 + b*pin1);
                                            break;
                                    }
                                }
                            } else { // channel data is real -> conventional RF beamforming: problem with pixel-based approach
                                for (int f = 0; f<F; f++) { // frame loop
                                    // reference to data
                                    float& prn0=*(Pr + L*rx + LN*w + LNW*f + n0);
                                    float& prn1=*(Pr + L*rx + LN*w + LNW*f + n0 + 1);
                                    
                                    switch(sum_mode) {
                                        case NONE:
                                            Dr[pp + P*rx + PN*w + PNW*f] = c_apo*(a*prn0 + b*prn1);
                                            break;
                                        case RX:
                                            Dr[pp + P*w + PW*f] += c_apo*(a*prn0 + b*prn1);
                                            break;
                                        case TX:
                                            Dr[pp + P*rx + PN*f] += c_apo*(a*prn0 + b*prn1);
                                            break;
                                        case BOTH:
                                            Dr[pp + P*f] += c_apo*(a*prn0 + b*prn1);
                                            break;
                                    }
                                }
                            }
                        }
                    }
                }   // end wave loop
            }   // end channel
        }); // end pixel loop
        
        mexPrintf(".....done!\n");
        mexEvalString("drawnow;");
        
        if (verbose) {
            mexPrintf("---------------------------------------------------\n");
            mexEvalString("drawnow;");
        }
        
        return;
        
    }
    
    
    
    
    
    
