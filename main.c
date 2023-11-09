
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <complex.h>
#include <assert.h>

#include "include/raylib.h"
#include "realfft.h"

#define Float_Complex float complex
#define TWO_PI 6.28318530717959
#define N (1<<11)
#define TARGET_FREQ_SIZE 9

float input_raw_Data[N] = {0.0};
float input_with_window_function[N] = {0.0};
float output_raw_Data[N];
float frequency_magnitude[N/2] = {0.0};

// coeeficient for Flat top window.
const float A[5] = {0.21557895, 0.41663158, 0.277263158, 0.083578947, 0.006947368};

//--------------------------------------------------------------------------------------
float target_frequencies[TARGET_FREQ_SIZE] = {20.0, 60.0, 250.0, 500.0, 2000.0, 4000.0, 6000.0, 16000.0, 22050.0};

float smooth_spectrum[TARGET_FREQ_SIZE-1] = {0.0};

void pushSample(float sample)
{
    memmove(input_raw_Data, input_raw_Data + 1, (N - 1)*sizeof(input_raw_Data[0])); // Moving history to the left
    input_raw_Data[N-1] = sample;                                                   // Adding last average value
}

void ProcessAudioStream(void *bufferData, unsigned int frames) 
{
    /**
     * https://cdecl.org/?q=float+%28*fs%29%5B2%5D
     * Samples internally stored as <float>s
    */
    float (*samples)[2] = bufferData;
    
    for (size_t frame = 0; frame < frames; frame++)
    {
        /* code */
        pushSample(samples[frame][0]);                   // We only interested in the left audio channel.
    }
    return;
}

float getAmp(float a, float b)
{
    // float mag = sqrtf(a*a + b*b);
    // return powf(mag, 2.0);
    return fabsf(a*a + b*b);
}

int main(void)
{
    // Initialization
    //--------------------------------------------------------------------------------------
    const int screenWidth  = 512;
    const int screenHeight = 512;
    float smoothingFactor  = 20.0;

    InitWindow(screenWidth, screenHeight, "Mckenzie - FFT AUDIO VISUALIZER");

    InitAudioDevice();              // Initialize audio device

    Music music = LoadMusicStream("resources/Cloud phonk.mp3");

    AttachAudioStreamProcessor(music.stream, ProcessAudioStream);

    PlayMusicStream(music);

    music.looping = false;

    unsigned int fs = music.stream.sampleRate;

    SetTargetFPS(60);               // Set to render at 60 frames-per-second
    //--------------------------------------------------------------------------------------

    // Main game loop
    while (!WindowShouldClose())    // Detect window close button or ESC key
    {
        // Get Framerate in seconds.
        float dt = GetFrameTime();
        printf("%f\n", dt);

        // Update
        //----------------------------------------------------------------------------------
        UpdateMusicStream(music);   // Update music buffer with new stream data

        //----------------------------------------------------------------------------------
        if(IsMusicReady(music))
        {
            // Apply the hanning window.
            for (size_t n = 0; n < N; n++)
            {
                /* code */
                float t = (float) n / (N-1);
                float h = (0.5 * (1.0 - cosf(TWO_PI * t)));
                input_with_window_function[n] = input_raw_Data[n] * h * 1.44;
            }
            
            // Do FFT
            memcpy(output_raw_Data, input_with_window_function, N*sizeof(float));
            realft(output_raw_Data, N, 1);

            // Calculate amplitudes.
            for (size_t c = 0; c < (N/2); c++)
            {
                /* code */
                float amp = getAmp(output_raw_Data[2*c], output_raw_Data[2*c+1]);      // even indexes are real values and odd indexes are complex value.
                frequency_magnitude[c] = amp;
            }       
                   
        }
        //----------------------------------------------------------------------------------
        if (!IsMusicStreamPlaying(music)) {
            memset(input_raw_Data, 0, N*sizeof(float));
            memset(input_with_window_function, 0, N*sizeof(float));
            memset(output_raw_Data, 0, N*sizeof(float));
            memset(frequency_magnitude, 0, (N/2)*sizeof(float));
            memset(smooth_spectrum, 0, (TARGET_FREQ_SIZE-1)*sizeof(float));
            break;
        }
        //----------------------------------------------------------------------------------

        float spectrum[TARGET_FREQ_SIZE-1] = {0.0};
        float n_freq[TARGET_FREQ_SIZE-1]   = {0.0};
        //----------------------------------------------------------------------------------

        /******************************************************************************************************************
         * Get the amplitude of the specific set of bins using Perseval's theorem.                                        *
         * https://www.ni.com/docs/en-US/bundle/labwindows-cvi/page/advancedanalysisconcepts/lvac_parseval_s_theorem.html *
         * https://www.gaussianwaves.com/2015/07/significance-of-rms-root-mean-square-value/                              *
         * https://stackoverflow.com/questions/15455698/extract-treble-and-bass-from-audio-in-ios                         *
         *                                                                                                                *
         ******************************************************************************************************************/
        for (size_t bin_index = 0; bin_index < (N/2); bin_index++)
        {
            /* code */
            float freq = bin_index * (fs / N);
            for (size_t j = 0; j < TARGET_FREQ_SIZE-1; j++)
            {
                /* code */
                if (freq >= target_frequencies[j] && freq < target_frequencies[j+1])
                {
                    spectrum[j] = spectrum[j] + frequency_magnitude[bin_index];
                    n_freq[j] = n_freq[j] + 1;
                }
            }
        }
        
        for (size_t i = 0; i < TARGET_FREQ_SIZE-1; i++)
        {
            /* code */
            float rms = sqrtf(spectrum[i]/n_freq[i]);        // Calculate the RMS value.
            float dBFS = 20.0*log10f(rms);                   // Convert RMS value to Decibel Full Scale (dBFS) value.
            if (dBFS < 0.0) { dBFS = 0.0; }
            spectrum[i] = dBFS;
            smooth_spectrum[i] += (spectrum[i] - smooth_spectrum[i]) * dt * smoothingFactor;    // Smoothing the spectrum output value.
        }
        //----------------------------------------------------------------------------------
        // Draw
        //----------------------------------------------------------------------------------
        BeginDrawing();
            ClearBackground(RAYWHITE);
            float h = GetRenderHeight() / 2;
            for (size_t i = 0; i < TARGET_FREQ_SIZE-1; i++)
            {
                /* code */
                DrawRectangleLines(128.0 + (i * 32), h - (3.0*smooth_spectrum[i]), 30, (3.0*smooth_spectrum[i]), BLACK);
            }         
        EndDrawing();
        
        //----------------------------------------------------------------------------------
    }

    // De-Initialization
    //--------------------------------------------------------------------------------------
    UnloadMusicStream(music);   // Unload music stream buffers from RAM

    DetachAudioStreamProcessor(music.stream, ProcessAudioStream);  // Disconnect audio stream processor

    CloseAudioDevice();         // Close audio device (music streaming is automatically stopped)

    CloseWindow();              // Close window and OpenGL context
    //--------------------------------------------------------------------------------------

    return 0;
}
