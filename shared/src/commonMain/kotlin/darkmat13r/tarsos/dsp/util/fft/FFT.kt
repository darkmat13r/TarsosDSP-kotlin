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
package darkmat13r.tarsos.dsp.util.fft

import kotlin.jvm.JvmOverloads
import kotlin.math.PI
import kotlin.math.atan2
import kotlin.math.sqrt


/**
 * Wrapper for calling a hopefully Fast Fourier transform. Makes it easy to
 * switch FFT algorithm with minimal overhead.
 * Support for window functions is also present.
 *
 * @author Joren Six
 */
class FFT constructor(size: Int, private val windowFunction: WindowFunction? = null) {
    /**
     * Forward FFT.
     */
    private val fft: FloatFFT = FloatFFT(size)
    private val fftSize: Int = size
    private val window: FloatArray? = windowFunction?.generateCurve(size)

    /**
     * Create a new fft of the specified size. Apply the specified window on the samples before a forward transform.
     * arning: the window is not applied in reverse when a backwards transform is requested.
     * @param size The size of the fft.
     * @param windowFunction Apply the specified window on the samples before a forward transform.
     * arning: the window is not applied in reverse when a backwards transform is requested.
     */

    /**
     * Computes forward DFT.
     *
     * @param data
     * data to transform.
     */
    suspend fun forwardTransform(data: FloatArray) {
        if (windowFunction != null) {
            for (i in data.indices) {
                data[i] = data[i] * window!![i]
            }
            //windowFunction.apply(data);
        }
        fft.realForward(data)
    }

    /**
     * do a complex forward transform
     * @param data do a complex forward transform on these complex numbers
     */
    suspend fun complexForwardTransform(data: FloatArray) {
        if (windowFunction != null) {
            for (i in data.indices) {
                data[i] = data[i] * window!![i]
            }
            //windowFunction.apply(data);
        }
        fft.complexForward(data)
    }

    /**
     * Computes inverse DFT.
     * Warning, does not reverse the window function.
     * @param data
     * data to transform
     */
    suspend fun backwardsTransform(data: FloatArray) {
        fft.realInverse(data, true)
    }

    /**
     * Calculate the frequency of the bin.
     * @param binIndex The index of the bin.
     * @param sampleRate The sample rate of the audio.
     * @return The frequency in Hz of the bin.
     */
    fun binToHz(binIndex: Int, sampleRate: Float): Double {
        return binIndex * sampleRate / fftSize.toDouble()
    }

    /**
     * The size of the fft.
     * @return The size of the fft.
     */
    fun size(): Int {
        return fftSize
    }

    /**
     * Returns the modulus of the element at index bufferCount. The modulus,
     * magnitude or absolute value is (a²+b²) ^ 0.5 with a being the real part
     * and b the imaginary part of a complex number.
     *
     * @param data
     * The FFT transformed data.
     * @param index
     * The index of the element.
     * @return The modulus, magnitude or absolute value of the element at index
     * bufferCount
     */
    fun modulus(data: FloatArray, index: Int): Float {
        val realIndex = 2 * index
        val imgIndex = 2 * index + 1
        val modulus = data[realIndex] * data[realIndex] + data[imgIndex] * data[imgIndex]
        return sqrt(modulus.toDouble()).toFloat()
    }

    /**
     * Calculates the the modulus for each element in data and stores the result
     * in amplitudes.
     *
     * @param data
     * The input data.
     * @param amplitudes
     * The output modulus info or amplitude.
     */
    fun modulus(data: FloatArray, amplitudes: FloatArray) {
        check(data.size / 2 == amplitudes.size)
        for (i in amplitudes.indices) {
            amplitudes[i] = modulus(data, i)
        }
    }

    /**
     * Computes an FFT and converts the results to polar coordinates (power and
     * phase). Both the power and phase arrays must be the same length, data
     * should be double the length.
     *
     * @param data
     * The input audio signal.
     * @param power
     * The power (modulus) of the data.
     * @param phase
     * The phase of the data
     */
    suspend fun powerPhaseFFT(data: FloatArray, power: FloatArray, phase: FloatArray) {
        check(data.size / 2 == power.size)
        check(data.size / 2 == phase.size)
        if (windowFunction != null) {
            windowFunction.apply(data)
        }
        fft.realForward(data)
        powerAndPhaseFromFFT(data, power, phase)
    }


    /**
     * Returns magnitude (or power) and phase for the FFT transformed data.
     * @param data The FFT transformed data.
     * @param power The array where the magnitudes or powers are going to be stored. It is half the length of data (FFT size).
     * @param phase The array where the phases are going to be stored. It is half the length of data (FFT size).
     */
    fun powerAndPhaseFromFFT(data: FloatArray, power: FloatArray, phase: FloatArray) {
        phase[0] = PI.toFloat()
        power[0] = -data[0]
        for (i in 1 until power.size) {
            val realIndex = 2 * i
            val imgIndex = 2 * i + 1
            power[i] = sqrt(
                (data[realIndex] * data[realIndex] + data[imgIndex] * data[imgIndex]).toDouble()
            ).toFloat()
            phase[i] = atan2(data[imgIndex].toDouble(), data[realIndex].toDouble())
                .toFloat()
        }
    }

    /**
     * Beatroot expects a different first power element
     * @param data The FFT transformed data.
     * @param power The array where the magnitudes or powers are going to be stored. It is half the length of data (FFT size).
     * @param phase The array where the phases are going to be stored. It is half the length of data (FFT size).
     */
    suspend fun powerPhaseFFTBeatRootOnset(data: FloatArray, power: FloatArray, phase: FloatArray) {
        powerPhaseFFT(data, power, phase)
        power[0] = sqrt((data[0] * data[0] + data[1] * data[1]).toDouble())
            .toFloat()
    }

    /**
     * Multiplies to arrays containing imaginary numbers. The data in the first argument
     * is modified! The real part is stored at `2*i`, the imaginary part `2*i+i`
     * @param data The array with imaginary numbers that is modified.
     * @param other The array with imaginary numbers that is not modified.
     * Data and other need to be the same length.
     */
    fun multiply(data: FloatArray, other: FloatArray) {
        check(data.size == other.size)
        if (data.size != other.size) {
            throw IllegalArgumentException("Both arrays with imaginary numbers shouldb e of equal length")
        }
        var i = 1
        while (i < data.size - 1) {
            val realIndex = i
            val imgIndex = i + 1
            val tempReal =
                data[realIndex] * other[realIndex] + -1 * data[imgIndex] * other[imgIndex]
            val tempImg = data[realIndex] * other[imgIndex] + data[imgIndex] * other[realIndex]
            data[realIndex] = tempReal
            data[imgIndex] = tempImg
            i += 2
        }
    }
}
