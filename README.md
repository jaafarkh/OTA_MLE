# OTA_MLE
Orer tracking analysis by Maximum Likelihood Estimator
This repository includes Matlab code for the Order Tracking Analysis using three methods:
(*) Our developed Maximum Likelihood Estimator
(*) Time Variant Discrete Fourier Transform of Blough
(*) Vold-Kalman Order Tracking of Vold
Please cite as:
Order Tracking Analysis Using Maximum Likeihood Estimator in the Presence of Crossing Orders and Low Resolution Tacho Signal
(to be published)

testOTA_Sig1.m: this code applies OTA for simulated data or data obtained from simulation of constant-amp orders
testOTA_Sig2.m:this code applies OTA for simulated data or data obtained from simulation of 2-DOF system
testOTA_Sig3.m: this code applies OTA for simulated data obtained from numerical simulation of a rectangular SSSS plate using Matlab code inside /Plate
testOTA_Sig4.m: this code applies OTA for actual vibration data obtained from experimental results, vibration data is located in /otaP/input/plateViv.txt

The code is supplied as is without any warranty for the purpose of reproduction of the results and/or be used by other scholars
