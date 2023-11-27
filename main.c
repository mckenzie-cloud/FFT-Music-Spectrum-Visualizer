#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

#include "include/raylib.h"
#include "realfft.h"

#define SCREEN_HEIGHT 512
#define SCREEN_WIDTH 512

#define TWO_PI 6.28318530717959
#define N (1<<12)
#define TARGET_FREQ_SIZE 10

#define BG_COLOR (Color) {0, 0, 0, 255}
#define TEXT_COLOR (Color) {244, 234, 224, 255}
#define SPECTRUM_COLOR (Color) {250, 246, 240, 255}
#define PROGRESS_BAR_COLOR (Color) {244, 223, 200, 255}

typedef struct
{
    /* data */
    float input_raw_Data[N];
    float input_data_with_windowing_function[N];
    float output_raw_Data[N];
    float frequency_sqr_mag[N/2];
    float smooth_spectrum[TARGET_FREQ_SIZE-1];
} Data;

Data data;

void pushSample(float sample)
{
    memmove(data.input_raw_Data, data.input_raw_Data + 1, (N - 1)*sizeof(data.input_raw_Data[0])); // Moving history to the left
    data.input_raw_Data[N-1] = sample;                                                   // Adding last average value
}

void ProcessAudioStreamCallback(void *bufferData, unsigned int frames) 
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

float getSqrMag(float a, float b)
{
    return a*a + b*b;
}

void calculateAmplitudes()
{
    for (size_t c = 0; c < (N/2); c++)
    {
        /* code */
        float sqr_mag = getSqrMag(data.output_raw_Data[2*c], data.output_raw_Data[2*c+1]);      // even indexes are real values and odd indexes are complex value.
        data.frequency_sqr_mag[c] = sqr_mag;
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
        float freq = (bin_index+1) * (fs / N);
        for (size_t j = 0; j < TARGET_FREQ_SIZE-1; j++)
        {
            /* code */
            if (freq >= target_frequencies[j] && freq < target_frequencies[j+1])
            {
                spectrum[j] = spectrum[j] + data.frequency_sqr_mag[bin_index];
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
        if (rms > 0.0)
        {
            float dBFS = 10.0*log10f(rms);               // Convert RMS value to Decibel Full Scale (dBFS) value.
            if (dBFS < 0.0) { dBFS = 0.0; }
            spectrum[i] = fabsf(dBFS);
        }
        data.smooth_spectrum[i] += (spectrum[i] - data.smooth_spectrum[i]) * dt * smoothingFactor;    // Smoothing the spectrum output value.
    }
}

void visualizeSpectrum(float spectrum_scaling_factor)
{
    float h = SCREEN_HEIGHT / 2;
    for (size_t i = 0; i < TARGET_FREQ_SIZE-1; i++)
    {
        /* code */
        DrawRectangleLines(184.0 + (i * 16), h - (spectrum_scaling_factor*data.smooth_spectrum[i]), 15, (spectrum_scaling_factor*data.smooth_spectrum[i]), SPECTRUM_COLOR);
    }  
}

void displayProgressBar(int time_played)
{
    DrawRectangle(32, (SCREEN_HEIGHT / 2) + 32 - 10, (int)time_played, 10, PROGRESS_BAR_COLOR); 
}

void cleanUp()
{
    memset(data.input_raw_Data, 0, N*sizeof(float));
    memset(data.input_data_with_windowing_function, 0, N*sizeof(float));
    memset(data.output_raw_Data, 0, N*sizeof(float));
    memset(data.frequency_sqr_mag, 0, (N/2)*sizeof(float));
    memset(data.smooth_spectrum, 0, (TARGET_FREQ_SIZE-1)*sizeof(float));
}

int main(void)
{
    // Initialization

    InitWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "Mckenzie - FFT AUDIO VISUALIZER");
    InitAudioDevice();              // Initialize audio device

    //--------------------------------------------------------------------------------------
    const char *music_file_path = "";
    Music music_stream;
    bool  music_loaded = false;

    //--------------------------------------------------------------------------------------
    float target_frequencies[TARGET_FREQ_SIZE] = {20.0, 40.0, 80.0, 160.0, 300.0, 600.0, 1200.0, 5000.0, 10000.0, 22050.0};
    float smoothingFactor  = 20.0;
    float spectrum_scaling_factor = 3.0;

    //--------------------------------------------------------------------------------------
    float durations   = 0.0;
    float time_played = 0.0f;
    unsigned int fs   = 0;
    const char *music_title = "";

    SetTargetFPS(60);               // Set to render at 60 frames-per-second

    // Main game loop
    while (!WindowShouldClose())    // Detect window close button or ESC key
    {
        // Get Framerate in seconds.
        float dt = GetFrameTime();

        // Update
        //----------------------------------------------------------------------------------
        if (music_loaded)
        {
            UpdateMusicStream(music_stream);   // Update music buffer with new stream data
        }

        // Handle drag & drop file.
        //----------------------------------------------------------------------------------
        if (IsFileDropped())
        {
            FilePathList droppedFiles = LoadDroppedFiles();     // Load dropped filepaths
            music_file_path = droppedFiles.paths[0];            // We are only interested in a single file.
            if (IsFileExtension(music_file_path, ".mp3"))       // Check file extension (including point: .png, .wav)
            {
                /**
                 * Before loading a new music file, make sure to unload any existing ones.
                */
                if (music_loaded)
                {
                    if (IsMusicStreamPlaying(music_stream))         // Check if music is playing
                    {
                        StopMusicStream(music_stream);              // Stop music playing
                    }
                    UnloadMusicStream(music_stream);                // Unload music stream buffers from RAM
                    DetachAudioStreamProcessor(music_stream.stream, ProcessAudioStreamCallback);  // Disconnect audio stream processor
                }

                music_stream = LoadMusicStream(music_file_path);          // Load music stream from file

                if (IsMusicReady(music_stream))                           // Checks if a music stream is ready
                {
                    music_loaded = true;
                    music_title  = GetFileNameWithoutExt(music_file_path);
                    durations    = GetMusicTimeLength(music_stream);
                    fs           = music_stream.stream.sampleRate;
                    //----------------------------------------------------------------------------------
                    cleanUp();
                    AttachAudioStreamProcessor(music_stream.stream, ProcessAudioStreamCallback);    // Attach audio stream processor to stream, receives the samples as <float>s
                    PlayMusicStream(music_stream);                        // Start music playing
                }
            }

            UnloadDroppedFiles(droppedFiles);                         // Unload dropped filepaths
        }

        //----------------------------------------------------------------------------------
        if(IsMusicStreamPlaying(music_stream))
        {
            // Apply the hanning window.
            applyHanningWindow();
            
            // Do FFT
            performFFT();

            // Calculate amplitudes.
            calculateAmplitudes();     

            //----------------------------------------------------------------------------------
            float spectrum[TARGET_FREQ_SIZE-1] = {0.0};
            float n_freq[TARGET_FREQ_SIZE-1]   = {0.0};

            performPersevalTheorem(spectrum, n_freq, target_frequencies, fs);
            
            RMS_TO_DBFS(spectrum, n_freq, dt, smoothingFactor);
            
            time_played = GetMusicTimePlayed(music_stream)/durations*(SCREEN_WIDTH - 64);
        }

        //----------------------------------------------------------------------------------
        // Draw
        //----------------------------------------------------------------------------------
        BeginDrawing();
            ClearBackground(BG_COLOR);
            DrawText("Drag and Drop a music file | MP3 FILE ONLY!", 5, 5, 12, TEXT_COLOR);
            if (IsMusicStreamPlaying(music_stream))
            {
                visualizeSpectrum(spectrum_scaling_factor);
                displayProgressBar((int)time_played);
                DrawText(music_title, 112, (SCREEN_HEIGHT/2) + 64 - 12, 12, TEXT_COLOR);  
            } 
        EndDrawing();
        
        //----------------------------------------------------------------------------------
    }

    // De-Initialization
    //--------------------------------------------------------------------------------------
    if (music_loaded)
    {
        UnloadMusicStream(music_stream);                                      // Unload music stream buffers from RAM
        DetachAudioStreamProcessor(music_stream.stream, ProcessAudioStreamCallback);  // Disconnect audio stream processor
    }

    CloseAudioDevice();         // Close audio device (music streaming is automatically stopped)

    CloseWindow();              // Close window and OpenGL context
    //--------------------------------------------------------------------------------------

    return EXIT_SUCCESS;
}
