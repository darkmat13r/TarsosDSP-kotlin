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
*
*  I took Joren's Code and changed it so that 
*  it uses the FFT to calculate the difference function.
*  TarsosDSP is developed by Joren Six at 
*  The Royal Academy of Fine Arts & Royal Conservatory,
*  University College Ghent,
*  Hoogpoort 64, 9000 Ghent - Belgium
*  
*  http://tarsos.0110.be/tag/TarsosDSP
*
*/
package darkmat13r.tarsos.dsp.pitch

import darkmat13r.tarsos.dsp.util.fft.FloatFFT


/**
 * An implementation of the YIN pitch tracking algorithm which uses an FFT to
 * calculate the difference function. This makes calculating the difference
 * function more performant. See [the YIN paper.](http://recherche.ircam.fr/equipes/pcm/cheveign/ps/2002_JASA_YIN_proof.pdf) This implementation is done by [Matthias Mauch](mailto:matthias.mauch@elec.qmul.ac.uk) and is
 * based on [Yin] which is based on the implementation found in [aubio](http://aubio.org) by Paul Brossier.
 *
 * @author Matthias Mauch
 * @author Joren Six
 * @author Paul Brossier
 */
class FastYin  (
    /**
     * The audio sample rate. Most audio has a sample rate of 44.1kHz.
     */
    private val sampleRate: Float, bufferSize: Int,
    /**
     * The actual YIN threshold.
     */
    private val threshold: Double = DEFAULT_THRESHOLD
) : PitchDetector() {
    /**
     * The buffer that stores the calculated values. It is exactly half the size
     * of the input buffer.
     */
    private val yinBuffer = FloatArray(bufferSize / 2)

    /**
     * The result of the pitch detection iteration.
     */
    private val result: PitchResult.Pitch

    //------------------------ FFT instance members
    //Initializations for FFT difference step
    /**
     * Holds the FFT data, twice the length of the audio buffer.
     */
    private val audioBufferFFT = FloatArray(2 * bufferSize)

    /**
     * Half of the data, disguised as a convolution kernel.
     */
    private val kernel = FloatArray(2 * bufferSize)

    /**
     * Buffer to allow convolution via complex multiplication. It calculates the auto correlation function (ACF).
     */
    private val yinStyleACF = FloatArray(2 * bufferSize)

    /**
     * An FFT object to quickly calculate the difference function.
     */
    private val fft: FloatFFT

    /**
     * Create a new pitch detector for a stream with the defined sample rate.
     * Processes the audio in blocks of the defined size.
     *
     * @param sampleRate
     * The sample rate of the audio stream. E.g. 44.1 kHz.
     * @param bufferSize
     * The size of a buffer. E.g. 1024.
     * @param threshold
     * The parameter that defines which peaks are kept as possible
     * pitch candidates. See the YIN paper for more details.
     */
    /**
     * Create a new pitch detector for a stream with the defined sample rate.
     * Processes the audio in blocks of the defined size.
     *
     * @param audioSampleRate
     * The sample rate of the audio stream. E.g. 44.1 kHz.
     * @param bufferSize
     * The size of a buffer. E.g. 1024.
     */
    init {
        fft = FloatFFT(bufferSize)
        result = PitchResult.Pitch()
    }

    /**
     * The main flow of the YIN algorithm. Returns a pitch value in Hz or -1 if
     * no pitch is detected.
     *
     * @return a pitch value in Hz or -1 if no pitch is detected.
     */
    override suspend fun detect(audioBuffer: FloatArray): PitchResult {
        val pitchInHertz: Float

        // step 2
        difference(audioBuffer)

        // step 3
        cumulativeMeanNormalizedDifference()

        // step 4
        val tauEstimate = absoluteThreshold()

        // step 5
        if (tauEstimate != -1) {
            val betterTau = parabolicInterpolation(tauEstimate)

            // step 6
            // TODO Implement optimization for the AUBIO_YIN algorithm.
            // 0.77% => 0.5% error rate,
            // using the data of the YIN paper
            // bestLocalEstimate()

            // conversion to Hz
            pitchInHertz = sampleRate / betterTau
        } else {
            // no pitch found
            pitchInHertz = -1f
        }

        result.hz = pitchInHertz

        return result
    }

    /**
     * Implements the difference function as described in step 2 of the YIN
     * paper with an FFT to reduce the number of operations.
     */
    private suspend fun difference(audioBuffer: FloatArray) {
        // POWER TERM CALCULATION
        // ... for the power terms in equation (7) in the Yin paper
        val powerTerms = FloatArray(yinBuffer.size)
        for (j in yinBuffer.indices) {
            powerTerms[0] += audioBuffer[j] * audioBuffer[j]
        }
        // now iteratively calculate all others (saves a few multiplications)
        for (tau in 1 until yinBuffer.size) {
            powerTerms[tau] =
                powerTerms[tau - 1] - audioBuffer[tau - 1] * audioBuffer[tau - 1] + audioBuffer[tau + yinBuffer.size] * audioBuffer[tau + yinBuffer.size]
        }

        // YIN-STYLE AUTOCORRELATION via FFT
        // 1. data
        for (j in audioBuffer.indices) {
            audioBufferFFT[2 * j] = audioBuffer[j]
            audioBufferFFT[2 * j + 1] = 0f
        }
        fft.complexForward(audioBufferFFT)


        // 2. half of the data, disguised as a convolution kernel
        for (j in yinBuffer.indices) {
            kernel[2 * j] = audioBuffer[yinBuffer.size - 1 - j]
            kernel[2 * j + 1] = 0f
            kernel[2 * j + audioBuffer.size] = 0f
            kernel[2 * j + audioBuffer.size + 1] = 0f
        }
        fft.complexForward(kernel)

        // 3. convolution via complex multiplication
        for (j in audioBuffer.indices) {
            yinStyleACF[2 * j] =
                audioBufferFFT[2 * j] * kernel[2 * j] - audioBufferFFT[2 * j + 1] * kernel[2 * j + 1] // real
            yinStyleACF[2 * j + 1] =
                audioBufferFFT[2 * j + 1] * kernel[2 * j] + audioBufferFFT[2 * j] * kernel[2 * j + 1] // imaginary
        }
        fft.complexInverse(yinStyleACF, true)


        // CALCULATION OF difference function
        // ... according to (7) in the Yin paper.
        for (j in yinBuffer.indices) {
            // taking only the real part
            yinBuffer[j] =
                powerTerms[0] + powerTerms[j] - 2 * yinStyleACF[2 * (yinBuffer.size - 1 + j)]
        }
    }

    /**
     * The cumulative mean normalized difference function as described in step 3
     * of the YIN paper. <br></br>
     * `
     * yinBuffer[0] == yinBuffer[1] = 1
    ` *
     */
    private fun cumulativeMeanNormalizedDifference() {
        yinBuffer[0] = 1f
        var runningSum = 0f
        var tau = 1
        while (tau < yinBuffer.size) {
            runningSum += yinBuffer[tau]
            yinBuffer[tau] *= tau / runningSum
            tau++
        }
    }

    /**
     * Implements step 4 of the AUBIO_YIN paper.
     */
    private fun absoluteThreshold(): Int {
        // Uses another loop construct
        // than the AUBIO implementation
        var tau: Int
        // first two positions in yinBuffer are always 1
        // So start at the third (index 2)
        tau = 2
        while (tau < yinBuffer.size) {
            if (yinBuffer[tau] < threshold) {
                while (tau + 1 < yinBuffer.size && yinBuffer[tau + 1] < yinBuffer[tau]) {
                    tau++
                }
                // found tau, exit loop and return
                // store the probability
                // From the YIN paper: The threshold determines the list of
                // candidates admitted to the set, and can be interpreted as the
                // proportion of aperiodic power tolerated
                // within a periodic signal.
                //
                // Since we want the periodicity and and not aperiodicity:
                // periodicity = 1 - aperiodicity
                result.probability = 1 - yinBuffer[tau]
                break
            }
            tau++
        }


        // if no pitch found, tau => -1
        if (tau == yinBuffer.size || yinBuffer[tau] >= threshold || result.probability > 1.0) {
            tau = -1
            result.probability = 0f
            result.pitched = false
        } else {
            result.pitched = true
        }

        return tau
    }

    /**
     * Implements step 5 of the AUBIO_YIN paper. It refines the estimated tau
     * value using parabolic interpolation. This is needed to detect higher
     * frequencies more precisely. See http://fizyka.umk.pl/nrbook/c10-2.pdf and
     * for more background
     * http://fedc.wiwi.hu-berlin.de/xplore/tutorials/xegbohtmlnode62.html
     *
     * @param tauEstimate
     * The estimated tau value.
     * @return A better, more precise tau value.
     */
    private fun parabolicInterpolation(tauEstimate: Int): Float {
        val betterTau: Float

        val x0 = if (tauEstimate < 1) {
            tauEstimate
        } else {
            tauEstimate - 1
        }
        val x2 = if (tauEstimate + 1 < yinBuffer.size) {
            tauEstimate + 1
        } else {
            tauEstimate
        }
        if (x0 == tauEstimate) {
            betterTau = if (yinBuffer[tauEstimate] <= yinBuffer[x2]) {
                tauEstimate.toFloat()
            } else {
                x2.toFloat()
            }
        } else if (x2 == tauEstimate) {
            betterTau = if (yinBuffer[tauEstimate] <= yinBuffer[x0]) {
                tauEstimate.toFloat()
            } else {
                x0.toFloat()
            }
        } else {
            val s0 = yinBuffer[x0]
            val s1 = yinBuffer[tauEstimate]
            val s2 = yinBuffer[x2]
            // fixed AUBIO implementation, thanks to Karl Helgason:
            // (2.0f * s1 - s2 - s0) was incorrectly multiplied with -1
            betterTau = tauEstimate + (s2 - s0) / (2 * (2 * s1 - s2 - s0))
        }
        return betterTau
    }

    companion object {
        /**
         * The default YIN threshold value. Should be around 0.10~0.15. See YIN
         * paper for more information.
         */
        private const val DEFAULT_THRESHOLD = 0.20

        /**
         * The default size of an audio buffer (in samples).
         */
        const val DEFAULT_BUFFER_SIZE: Int = 2048

        /**
         * The default overlap of two consecutive audio buffers (in samples).
         */
        const val DEFAULT_OVERLAP: Int = 1536
    }
}
