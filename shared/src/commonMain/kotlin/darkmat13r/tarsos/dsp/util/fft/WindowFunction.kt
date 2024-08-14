/*
*      _______                       _____   _____ _____  
*     |__   __|                     |  __ \ / ____|  __ \ 
*        | | __ _ _ __ ___  ___  ___| |  | | (___ | |__) |
*        | |/ _` | '__/ __|/ _ \/ __| |  | |\___ \|  ___/ 
*        | | (_| | |  \__ \ (_) \__ \ |__| |____) | |     
*        |_|\__,_|_|  |___/\___/|___/_____/|_____/|_|     
*                                                         
* -------------------------------------------------------------
*
* TarsosDSP is developed by Joren Six at IPEM, University Ghent
*  
* -------------------------------------------------------------
*
*  Info: http://0110.be/tag/TarsosDSP
*  Github: https://github.com/JorenSix/TarsosDSP
*  Releases: http://0110.be/releases/TarsosDSP/
*  
*  TarsosDSP includes modified source code by various authors,
*  for credits and info, see README.
* 
*/
/*
 *  Copyright (c) 2007 - 2008 by Damien Di Fede <ddf@compartmental.net>
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU Library General Public License as published
 *   by the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Library General Public License for more details.
 *
 *   You should have received a copy of the GNU Library General Public
 *   License along with this program; if not, write to the Free Software
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */
package darkmat13r.tarsos.dsp.util.fft

import kotlin.math.PI

/**
 * A Window function represents a curve which is applied to a sample buffer to
 * reduce the introduction of spectral leakage in the Fourier transform.
 *
 *
 *
 * **Windowing**
 *
 *
 * Windowing is the process of shaping the audio samples before transforming
 * them to the frequency domain. The Fourier Transform assumes the sample buffer
 * is is a repetitive signal, if a sample buffer is not truly periodic within
 * the measured interval sharp discontinuities may arise that can introduce
 * spectral leakage. Spectral leakage is the speading of signal energy across
 * multiple FFT bins. This "spreading" can drown out narrow band signals and
 * hinder detection.
 *
 *
 * A [windowing
 * function](http://en.wikipedia.org/wiki/Window_function) attempts to reduce spectral leakage by attenuating the measured
 * sample buffer at its end points to eliminate discontinuities. If you call the
 * `window()` function with an appropriate WindowFunction, such as
 * `HammingWindow()`, the sample buffers passed to the object for
 * analysis will be shaped by the current window before being transformed. The
 * result of using a window is to reduce the leakage in the spectrum somewhat.
 *
 *
 * `WindowFunction` handles work associated with various window
 * functions such as the Hamming window. To create your own window function you
 * must extend `WindowFunction` and implement the
 * [value][.value] method which defines the shape of the window
 * at a given offset. `WindowFunction` will call this method to apply
 * the window to a sample buffer. The number passed to the method is an offset
 * within the length of the window curve.
 *
 * @author Damien Di Fede
 * @author Corban Brook
 */
abstract class WindowFunction
/**
 * Construct a new window.
 */
protected constructor() {
    protected var length: Int = 0


    /**
     * Apply the window function to a sample buffer.
     *
     * @param samples
     * a sample buffer
     */
    fun apply(samples: FloatArray) {
        this.length = samples.size

        for (n in samples.indices) {
            samples[n] *= value(samples.size, n)
        }
    }

    /**
     * Generates the curve of the window function.
     *
     * @param length
     * the length of the window
     * @return the shape of the window function
     */
    fun generateCurve(length: Int): FloatArray {
        val samples = FloatArray(length)
        for (n in 0 until length) {
            samples[n] = 1f * value(length, n)
        }
        return samples
    }

    /**
     * The value of the window function
     * @param length with the lengt of the window (in samples)
     * @param index at index
     * @return The value of the window function at the requested index.
     */
    protected abstract fun value(length: Int, index: Int): Float

    companion object {


    }
}
