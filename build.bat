@echo off
gcc -Wall -O3 -std=c11 fourier1.c realfft.c main.c -o run -I include/ -L lib/ -lraylib -lopengl32 -lgdi32 -lwinmm