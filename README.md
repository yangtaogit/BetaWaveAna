# BetaWaveAna
Beta test wave analysis example

```cpp
string W7_350V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-II-F5-L1-15_70_beta/220V_350V/C2BiasV";  //Reference channel 
string w7_350V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-II-F5-L1-15_70_beta/220V_350V/C3BiasV";  // DUT     
WaveAnalysis W7_350V(W7_350V_t1,w7_350V_t2);
W7_350V.FullAnalysis(); 
```
