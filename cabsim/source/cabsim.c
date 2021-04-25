#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "fftw3.h"
#include "lv2/lv2plug.in/ns/lv2core/lv2.h"
#include "IRdata.h"
#include "sndfile.h"

/**********************************************************************************************************************************************************/
#define REAL 0
#define IMAG 1

#define SIZE 512

//plugin URI
#define PLUGIN_URI "http://VeJaPlugins.com/plugins/Release/cabsim"

//macro for Volume in DB to a coefficient
#define DB_CO(g) ((g) > -90.0f ? powf(10.0f, (g) * 0.05f) : 0.0f)

typedef enum {IN, OUT, ATTENUATE, MODEL}PortIndex;

fftwf_complex *outComplex;
fftwf_complex *IRout;
fftwf_complex *convolved;

fftwf_plan fft;
fftwf_plan ifft;
fftwf_plan IRfft;

/**********************************************************************************************************************************************************/

            // TOEDOE:: make a worker thread for non realtime IR.wav loading. might be used in the future for this plugin :) everything else works

/**********************************************************************************************************************************************************/

typedef struct{
    float const *in;
    float *out;
    float *model;
    float *outbuf;
    float *inbuf;
    float *IR;
    float *overlap;
    float *oA;
    float *oB;
    float *oC;
    //float *impulseResponse;
    const float *attenuation;
} Cabsim;
/**********************************************************************************************************************************************************/
//functions

/**********************************************************************************************************************************************************/
static LV2_Handle
instantiate(const LV2_Descriptor*   descriptor,
double                              samplerate,
const char*                         bundle_path,
const LV2_Feature* const* features)
{
    Cabsim* cabsim = (Cabsim*)malloc(sizeof(Cabsim));
    return (LV2_Handle)cabsim;
}
/**********************************************************************************************************************************************************/
static void connect_port(LV2_Handle instance, uint32_t port, void *data)
{
    Cabsim* cabsim = (Cabsim*)instance;

    switch ((PortIndex)port)
    {
        case IN:
            cabsim->in = (float*) data;
            break;
        case OUT:
            cabsim->out = (float*) data;
            break;
        case ATTENUATE:
            cabsim->attenuation = (const float*) data;
            break;
        case MODEL:
            cabsim->model = (float*) data;
    }
}
/**********************************************************************************************************************************************************/
void activate(LV2_Handle instance)
{
    Cabsim* const cabsim = (Cabsim*)instance;
/*
    cabsim->impulseResponse = (float *) malloc((SIZE)*sizeof(float)); 

    SF_INFO sndInfo;
    SNDFILE *sndFile = sf_open("/root/.lv2/cabsim.lv2/test2_48000.wav", SFM_READ, &sndInfo);
    int channelsCount = sndInfo.channels;

    int bufferSize = 256;
    float buffer[bufferSize * channelsCount];
    
    cabsim->impulseResponse = (float*) malloc(
        sndInfo.frames * channelsCount * sizeof(float)
    );

    int length = sndInfo.frames;
    int offset = 0;

    while (length) {
        int n = (length > bufferSize) ? bufferSize : length;

        n = sf_readf_float(sndFile, buffer, n);

        for (int i = 0; i < n * channelsCount; i++) {
            cabsim->impulseResponse[offset + i] = buffer[i];
        }

        offset += n * channelsCount;
        length -= n;
    }
    sf_close(sndFile);
*/
    outComplex = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*(SIZE));
    IRout =  (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*(SIZE));
    convolved =  (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*(SIZE));

    cabsim->overlap = (float *) malloc((SIZE)*sizeof(float));
    cabsim->outbuf = (float *) malloc((SIZE)*sizeof(float));
    cabsim->inbuf = (float *) malloc((SIZE)*sizeof(float));
    cabsim->IR = (float *) malloc((SIZE)*sizeof(float));
    cabsim->oA = (float *) malloc((SIZE)*sizeof(float));
    cabsim->oB = (float *) malloc((SIZE)*sizeof(float));
    cabsim->oC = (float *) malloc((SIZE)*sizeof(float));

    fft = fftwf_plan_dft_r2c_1d(SIZE, cabsim->inbuf, outComplex, FFTW_ESTIMATE); 
    IRfft = fftwf_plan_dft_r2c_1d(SIZE, cabsim->IR, IRout, FFTW_ESTIMATE);
    ifft = fftwf_plan_dft_c2r_1d(SIZE, convolved, cabsim->outbuf, FFTW_ESTIMATE);
}

/**********************************************************************************************************************************************************/
void run(LV2_Handle instance, uint32_t n_samples)
{
    Cabsim* const cabsim = (Cabsim*)instance;    

    const float const *in = cabsim->in;
    float *out = cabsim->out;
    float *outbuf = cabsim->outbuf;
    float *inbuf = cabsim->inbuf;
    float *IR = cabsim->IR;
    float *overlap = cabsim->overlap;
    float *oA = cabsim->oA;
    float *oB = cabsim->oB;
    float *oC = cabsim->oC;
    //float *impulseResponse = cabsim->impulseResponse;
    int model = (int)(*(cabsim->model));
    const float attenuation = *cabsim->attenuation; 

    const float coef = DB_CO(attenuation);

    uint32_t i, j, m;
    
    int multiplier = 1;

    if(n_samples == 128)
        multiplier = 4;
    else if(n_samples == 256)
        multiplier = 2;

    //copy inputbuffer and IR buffer with zero padding.
    for ( i = 0; i < n_samples * multiplier; i++)
    {
        inbuf[i] = (i < n_samples) ? (in[i] * coef * 0.2f): 0.0f;
        IR[i] = (i < n_samples) ? convolK(model,i) : 0.0f;
    }

    fftwf_execute(fft);
    fftwf_execute(IRfft);

    //complex multiplication
    for(m = 0; m < ((n_samples / 2) * multiplier) ;m++)
    {
        //real component
        convolved[m][REAL] = outComplex[m][REAL] * IRout[m][REAL] - outComplex[m][IMAG] * IRout[m][IMAG];
        //imaginary component
        convolved[m][IMAG] = outComplex[m][REAL] * IRout[m][IMAG] + outComplex[m][IMAG] * IRout[m][REAL];
    }

    fftwf_execute(ifft);

    //normalize output with overlap add.
    if(n_samples == 256)
    {
        for ( j = 0; j < n_samples * multiplier; j++)
        {
            if(j < n_samples)
            {
                out[j] = ((outbuf[j] / (n_samples * multiplier)) + overlap[j]);
            }
            else
            {
                overlap[j - n_samples] = outbuf[j]  / (n_samples * multiplier);
            }
        }
    }
    else if (n_samples == 128)      //HIER VERDER GAAN!!!!!! oA, oB, oC changed malloc to calloc. (initiate buffer with all zeroes)
    {
        for ( j = 0; j < n_samples * multiplier; j++)
        {
            if(j < n_samples)   //runs 128 times filling the output buffer with overap add
            {
                out[j] = (outbuf[j] / (n_samples * multiplier) + oA[j] + oB[j] + oC[j]);
            }
            else   
            {
                oC[j - n_samples] = oB[j]; // 128 samples of usefull data
                oB[j - n_samples] = oA[j];  //filled with samples 128 to 255 of usefull data
                oA[j - n_samples] = (outbuf[j] / (n_samples * multiplier)); //filled with 384 samples
            }
        }
    }
    memset(IR, 0, sizeof(SIZE));

}

/**********************************************************************************************************************************************************/
void deactivate(LV2_Handle instance)
{
    // TODO: include the deactivate function code here
}
/**********************************************************************************************************************************************************/
void cleanup(LV2_Handle instance)
{
    fftwf_destroy_plan(fft);
    fftwf_destroy_plan(ifft);
    fftwf_destroy_plan(IRfft);
    //free fft memory
    fftwf_free(outComplex);
    fftwf_free(IRout);
    fftwf_free(convolved);
    //free allocated memory
    free(instance);
}
/**********************************************************************************************************************************************************/
const void* extension_data(const char* uri)
{
    return NULL;
}
/**********************************************************************************************************************************************************/
static const LV2_Descriptor Descriptor = {
    PLUGIN_URI,
    instantiate,
    connect_port,
    activate,
    run,
    deactivate,
    cleanup,
    extension_data
};
/**********************************************************************************************************************************************************/
LV2_SYMBOL_EXPORT
const LV2_Descriptor* lv2_descriptor(uint32_t index)
{
    if (index == 0) return &Descriptor;
    else return NULL;
}
/**********************************************************************************************************************************************************/
