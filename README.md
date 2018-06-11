# matlabbing
### Various MATLAB scripts I have created over the years:
##### (I also did a lot of MATLAB for my Master's and Bachelor's theses, which can be found in their respective repos...)                     


#### Sobol Hemisphere (sobol_hemisphere.m)

Checking how [Sobol sequence](https://en.wikipedia.org/wiki/Sobol_sequence) samples map to a unit hemisphere with a particular FOV cone using the [Shirley-Chiu](https://doi.org/10.1080/10867651.1997.10487479) square to disk mapping. I did this while coding my own [path tracer](https://en.wikipedia.org/wiki/Path_tracing) as part of the _Advanced Computer Graphics_ course at Aalto University.        


#### Serato DJ Control Signal Analysis (serato.m)

Investigating how the [Serato DJ / Scractch Live](https://serato.com/) DJ software's timecode [control signal](https://en.wikipedia.org/wiki/Vinyl_emulation_software) works. More on my [webpage](http://www.esgrove.fi/analysing-the-serato-dj-timecode-signal/).              


#### Pi value approximation using Monte Carlo integration (montecarlo_pi.m)

I implemented the classic example of [Monte Carlo integration](https://en.wikipedia.org/wiki/Monte_Carlo_integration) for approximating the value of Pi for myself when studying Monte Carlo methods. A really stupid and slow method for approximating the value of Pi but it illustrates nicely the basic idea behind Monte Carlo methods.


#### True-peak Limiting (true_peak.m)

Investigating digital audio sampling and [true-peaks](https://techblog.izotope.com/2015/08/24/true-peak-detection/) (_inter-sample distortion_), especially how excessively most tracks tend to exceed 0 dBTP in actuality though they are seemingly staying at or below 0 dBFS. When these tracks are converted from the digital samples to continuous-time waveforms (i.e., analog signals), over 0 dBTP peaks can and will result in distortion. This is why all audio should be mastered with true-peak limiting, which is sadly not the case with many commercial audio releases even today. A great scientific paper regarding this topic is [Stop Counting Samples](http://www.aes.org/e-lib/browse.cfm?elib=13806) by Thomas Lund.


#### Dithering (dither.m)

Investigating dithering in practice, utilizing sound samples (sine waves) dithered with the [iZotope Ozone 7 Maximizer plugin](https://www.izotope.com/). Most audio people quote the common "_-96 dB for 16-bit audio_" and the "_6 dB per bit_" rule when talking about the [dynamic range of digital audio](https://en.wikipedia.org/wiki/Dynamic_range#Audio), but these simplifications are not the actual truth, which a lot of people don't seem to be aware of. It is quite possible have an audible -120 dB peak or even quieter frequency component in a 16-bit audio file even though noise floor is understood to be at -96 dB. In practice, the effective dynamic range of 16-bit audio is much higher (typically around 120 dB) with the help of noiseshaping, which is demostrated in this script.      


#### Pixel filtering (pixel_filter.m)

Simple script for checking out the shape of common filtering kernels for images (in this case, a gaussian and the Mitchell-Netravali filter) depending on the parameters and implemention style. Done while implementing antialiasing filters to my own 3D-rendering programs.   