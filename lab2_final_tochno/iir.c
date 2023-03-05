#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "data1.h"

int N = 20;
float Pi=3.141592653589793;

float in_var1[100]={A_ARRAY};
float window_rect[100];
float window_bartlett[100];
float window_hann[100];
float in_var1_x_window_rect[100];
float in_var1_x_window_bartlett[100];
float in_var1_x_window_hann[100];

float fft_in_var1[100];
float fft_window_rect[100];
float fft_window_bartlett[100];
float fft_window_hann[100];
float fft_in_var1_x_window_rect[100];
float fft_in_var1_x_window_bartlett[100];
float fft_in_var1_x_window_hann[100];

void window_one(){
	int i=0;
	for(i=0;i<100;i++){
		if(i>=0 && i<=(N-1)){
			window_rect[i]=1;
		}
		else{
			window_rect[i]=0;
		}
	}
}

void window_two(){
	int i=0;
	for(i=0;i<100;i++){
	    if(i<21){
    		if(i<11){
    			window_bartlett[i]=i;
    		}
    		else{
    			window_bartlett[i]=20-i;
    		}
	    }
	    else{
	        window_bartlett[i]=0;
	    }
	}
}

void window_three(){
    int i=0;
	for(i=0;i<100;i++){
		if(i<21){
			window_hann[i]=0.5-0.5*cos(2*Pi*i/N);
		}
		else{
			window_hann[i]=0;
		}
	}
}

void multiply(){
    int i = 0;
    for(i=0;i<100;i++){
        in_var1_x_window_rect[i]=in_var1[i]*window_rect[i];
        in_var1_x_window_bartlett[i]=in_var1[i]*window_bartlett[i];
        in_var1_x_window_hann[i]=in_var1[i]*window_hann[i];
    }
}


void calculateDFT()
{
    float Xr[100];
    float Xi[100];
    int k, n = 0;
    
	for (k = 0; k < 100; k++) {
        Xr[k] = 0;
        Xi[k] = 0;
        for (n = 0; n < 100; n++) {
            Xr[k] = (Xr[k] + in_var1[n] * cos(2 * 3.141592 * k * n / 100));
            Xi[k] = (Xi[k]- in_var1[n] * sin(2 * 3.141592 * k * n / 100));
        }
        fft_in_var1[k]=pow(Xr[k]*Xr[k]+Xi[k]*Xi[k],0.5);
    }

	for (k = 0; k < 100; k++) {
        Xr[k] = 0;
        Xi[k] = 0;
        for (n = 0; n < 100; n++) {
            Xr[k] = (Xr[k] + window_rect[n] * cos(2 * 3.141592 * k * n / 100));
            Xi[k] = (Xi[k]- window_rect[n] * sin(2 * 3.141592 * k * n / 100));
        }
        fft_window_rect[k]=pow(Xr[k]*Xr[k]+Xi[k]*Xi[k],0.5);
    }

	for (k = 0; k < 100; k++) {
        Xr[k] = 0;
        Xi[k] = 0;
        for (n = 0; n < 100; n++) {
            Xr[k] = (Xr[k] + window_bartlett[n] * cos(2 * 3.141592 * k * n / 100));
            Xi[k] = (Xi[k]- window_bartlett[n] * sin(2 * 3.141592 * k * n / 100));
        }
        fft_window_bartlett[k]=pow(Xr[k]*Xr[k]+Xi[k]*Xi[k],0.5);
    }

	for (k = 0; k < 100; k++) {
        Xr[k] = 0;
        Xi[k] = 0;
        for (n = 0; n < 100; n++) {
            Xr[k] = (Xr[k] + window_hann[n] * cos(2 * 3.141592 * k * n / 100));
            Xi[k] = (Xi[k]- window_hann[n] * sin(2 * 3.141592 * k * n / 100));
        }
        fft_window_hann[k]=pow(Xr[k]*Xr[k]+Xi[k]*Xi[k],0.5);
    }

    for (k = 0; k < 100; k++) {
        Xr[k] = 0;
        Xi[k] = 0;
        for (n = 0; n < 100; n++) {
            Xr[k] = (Xr[k] + in_var1_x_window_rect[n] * cos(2 * 3.141592 * k * n / 100));
            Xi[k] = (Xi[k]- in_var1_x_window_rect[n] * sin(2 * 3.141592 * k * n / 100));
        }
        fft_in_var1_x_window_rect[k]=pow(Xr[k]*Xr[k]+Xi[k]*Xi[k],0.5);
    }
    
    for (k = 0; k < 100; k++) {
        Xr[k] = 0;
        Xi[k] = 0;
        for (n = 0; n < 100; n++) {
            Xr[k] = (Xr[k] + in_var1_x_window_bartlett[n] * cos(2 * 3.141592 * k * n / 100));
            Xi[k] = (Xi[k]- in_var1_x_window_bartlett[n] * sin(2 * 3.141592 * k * n / 100));
        }
        fft_in_var1_x_window_bartlett[k]=pow(Xr[k]*Xr[k]+Xi[k]*Xi[k],0.5);
    }
    
    for (k = 0; k < 100; k++) {
        Xr[k] = 0;
        Xi[k] = 0;
        for (n = 0; n < 100; n++) {
            Xr[k] = (Xr[k] + in_var1_x_window_hann[n] * cos(2 * 3.141592 * k * n / 100));
            Xi[k] = (Xi[k]- in_var1_x_window_hann[n] * sin(2 * 3.141592 * k * n / 100));
        }
        fft_in_var1_x_window_hann[k]=pow(Xr[k]*Xr[k]+Xi[k]*Xi[k],0.5);
    }
}

void main()
{
	window_one();
	window_two();
	window_three();
	multiply();
	calculateDFT();
}

