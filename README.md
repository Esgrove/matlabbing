# matlabbing
#### Various MATLAB scripts I have created over the years:                   



### Sobol Hemisphere (sobol_hemisphere.m)

Checking how [Sobol sequence](https://en.wikipedia.org/wiki/Sobol_sequence) samples map to a unit hemisphere with a particular FOV cone using the [Shirley-Chiu](https://doi.org/10.1080/10867651.1997.10487479) square to disk mapping. I did this while coding my own [path tracer](https://en.wikipedia.org/wiki/Path_tracing) as part of the _Advanced Computer Graphics_ course at Aalto University.        


### Serato DJ Control Signal Analysis (serato.m)

Investigating and understanding how the [Serato DJ / Scractch Live](https://serato.com/) DJ-software's timecode [control signal](https://en.wikipedia.org/wiki/Vinyl_emulation_software) works. More on my [webpage](http://www.esgrove.fi/analysing-the-serato-dj-timecode-signal/).              


### Pi value approximation using Monte Carlo integration (montecarlo_pi.m)

I implemented the classic example of [Monte Carlo integration](https://en.wikipedia.org/wiki/Monte_Carlo_integration) for approximating the value of pi for myself when studying Monte Carlo methods.       


### True-peak Limiting (true_peak.m)

Investigating digital audio sampling and [true-peaks](https://techblog.izotope.com/2015/08/24/true-peak-detection/) (_inter-sample distortion_), especially how excessively most tracks tend to exceed 0 dBTP in actuality though they are seemingly staying at or below 0 dBFS.


### Dithering (dither.m)

Investigating dithering in practice, utilizing sound samples (sine waves) processed (dithered) with the [iZotope Ozone 7 Maximizer plugin](https://www.izotope.com/). Most audio people quote the common "_-96 dB for 16-bit audio_" and the "_6 dB per bit_" rule when talking about the [dynamic range of digital audio](https://en.wikipedia.org/wiki/Dynamic_range#Audio), but these simplifications are not the actual truth, which a lot of people don't seem to understand. It is quite possible have a audible -120 dB peak or even quieter audio component in a 16-bit audio file even though noise floor is understood to be at -96 dB. In practice, the effective dynamic range of 16-bit audio is much higher (closer to 120 dB) through the help of noiseshaping, which is demostrated in this script.                


