
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <complex.h>
#include <assert.h>

#include "include/raylib.h"
#include "realfft.h"

#define SCREEN_HEIGHT 512
#define SCREEN_WIDTH 512

#define TWO_PI 6.28318530717959
#define N (1<<11)
#define TARGET_FREQ_SIZE 9

#define BG_COLOR (Color) {237, 244, 242, 255}
#define TEXT_COLOR (Color) {115, 93, 51, 255}
#define SPECTRUM_COLOR (Color) {124, 131, 99, 255}
#define PROGRESS_BAR_COLOR (Color) {49, 71, 58, 255}

typedef struct
{
    /* data */
    float input_raw_Data[N];
    float input_data_with_windowing_function[N];
    float output_raw_Data[N];
    float frequency_magnitude[N/2];
    float smooth_spectrum[TARGET_FREQ_SIZE-1];
} Data;

Data data;

void pushSample(float sample)
{
    memmove(data.input_raw_Data, data.input_raw_Data + 1, (N - 1)*sizeof(data.input_raw_Data[0])); // Moving history to the left
    data.input_raw_Data[N-1] = sample;                                                   // Adding last average value
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

void applyHanningWindow()
{
    for (size_t n = 0; n < N; n++)
    {
        /* code */
        float t = (float) n / (N-1);
        float h = (0.5 * (1.0 - cosf(TWO_PI * t)));
        data.input_data_with_windowing_function[n] = data.input_raw_Data[n] * h * 1.44;
    }
}

void performFFT()
{
    memcpy(data.output_raw_Data, data.input_data_with_windowing_function, N*sizeof(float));
    realft(data.output_raw_Data, N, 1);
}

float getAmp(float a, float b)
{
    // float mag = sqrtf(a*a + b*b);
    // return powf(mag, 2.0);
    return fabsf(a*a + b*b);
}

void calculateAmplitudes()
{
    for (size_t c = 0; c < (N/2); c++)
    {
        /* code */
        float amp = getAmp(data.output_raw_Data[2*c], data.output_raw_Data[2*c+1]);      // even indexes are real values and odd indexes are complex value.
        data.frequency_magnitude[c] = amp;
    }  
}

void performPersevalTheorem(float spectrum[], float n_freq[], float target_frequencies[], unsigned int fs)
{

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
                spectrum[j] = spectrum[j] + data.frequency_magnitude[bin_index];
                n_freq[j] = n_freq[j] + 1;
            }
        }
    }
}

void RMS_TO_DBFS(float spectrum[], float n_freq[], float dt, float smoothingFactor)
{
    for (size_t i = 0; i < TARGET_FREQ_SIZE-1; i++)
    {
        /* code */
        float rms = sqrtf(spectrum[i]/n_freq[i]);        // Calculate the RMS value.
        float dBFS = 20.0*log10f(rms);                   // Convert RMS value to Decibel Full Scale (dBFS) value.
        if (dBFS < 0.0) { dBFS = 0.0; }
        spectrum[i] = dBFS;
        data.smooth_spectrum[i] += (spectrum[i] - data.smooth_spectrum[i]) * dt * smoothingFactor;    // Smoothing the spectrum output value.
    }
}

void visualizeSpectrum()
{
    float h = SCREEN_HEIGHT / 2;
    for (size_t i = 0; i < TARGET_FREQ_SIZE-1; i++)
    {
        /* code */
        DrawRectangleLines(128.0 + (i * 32), h - (3.0*data.smooth_spectrum[i]), 30, (3.0*data.smooth_spectrum[i]), SPECTRUM_COLOR);
    }  
}

void displayProgressBar(int time_played)
{
    DrawRectangle(128, (SCREEN_HEIGHT / 2) + 32 - 12, (int)time_played, 12, PROGRESS_BAR_COLOR); 
}

void cleanUp()
{
    memset(data.input_raw_Data, 0, N*sizeof(float));
    memset(data.input_data_with_windowing_function, 0, N*sizeof(float));
    memset(data.output_raw_Data, 0, N*sizeof(float));
    memset(data.frequency_magnitude, 0, (N/2)*sizeof(float));
    memset(data.smooth_spectrum, 0, (TARGET_FREQ_SIZE-1)*sizeof(float));
}

int main(void)
{
    // Initialization

    InitWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "Mckenzie - FFT AUDIO VISUALIZER");
    InitAudioDevice();              // Initialize audio device
    Music music = LoadMusicStream("resources/KIRA-vynth.mp3");
    AttachAudioStreamProcessor(music.stream, ProcessAudioStream);
    PlayMusicStream(music);
    music.looping = false;

    //--------------------------------------------------------------------------------------
    float target_frequencies[TARGET_FREQ_SIZE] = {20.0, 60.0, 250.0, 500.0, 2000.0, 4000.0, 6000.0, 16000.0, 22050.0};
    float smoothingFactor  = 20.0;
    unsigned int fs = music.stream.sampleRate;

    //--------------------------------------------------------------------------------------
    float durations = GetMusicTimeLength(music);
    float time_played = 0.0f;
    char *music_title = "KIRA - Vynth";

    SetTargetFPS(60);               // Set to render at 60 frames-per-second

    // Main game loop
    while (!WindowShouldClose())    // Detect window close button or ESC key
    {
        // Get Framerate in seconds.
        float dt = GetFrameTime();

        // Update
        //----------------------------------------------------------------------------------
        UpdateMusicStream(music);   // Update music buffer with new stream data

        //----------------------------------------------------------------------------------
        if(IsMusicReady(music))
        {
            // Apply the hanning window.
            applyHanningWindow();
            
            // Do FFT
            performFFT();

            // Calculate amplitudes.
            calculateAmplitudes();     
                   
        }
        //----------------------------------------------------------------------------------
        if (!IsMusicStreamPlaying(music)) {
            cleanUp();
            CloseWindow();
        }
        //----------------------------------------------------------------------------------
        float spectrum[TARGET_FREQ_SIZE-1] = {0.0};
        float n_freq[TARGET_FREQ_SIZE-1]   = {0.0};

        performPersevalTheorem(spectrum, n_freq, target_frequencies, fs);
        
        RMS_TO_DBFS(spectrum, n_freq, dt, smoothingFactor);

        //----------------------------------------------------------------------------------
        time_played = GetMusicTimePlayed(music)/durations*(SCREEN_WIDTH - 256);

        //----------------------------------------------------------------------------------
        // Draw
        //----------------------------------------------------------------------------------
        BeginDrawing();
            ClearBackground(BG_COLOR);
            visualizeSpectrum();
            displayProgressBar((int)time_played);
            DrawText(music_title, 128, (SCREEN_HEIGHT/2) + 64 - 20, 20, TEXT_COLOR);    
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
