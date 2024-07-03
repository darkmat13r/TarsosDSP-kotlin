package darkmat13r.tarsos.dsp.pitch

import kotlin.math.max

class McLeodPitchMethod(
    /**
     * The audio sample rate. Most audio has a sample rate of 44.1kHz.
     */
    val audioSampleRate: Float,
    private val audioBufferSize: Int = DEFAULT_BUFFER_SIZE,
    /**
     * Defines the relative size the chosen peak (pitch) has.
     */
    val cutOff: Double = DEFAULT_CUTOFF

) : PitchDetector() {
    companion object {
        /**
         * The expected size of an audio buffer (in samples).
         */
        private const val DEFAULT_BUFFER_SIZE = 1024

        /**
         * Overlap defines how much two audio buffers following each other should
         * overlap (in samples). 75% overlap is advised in the MPM article.
         */
        private const val DEFAULT_OVERLAP = 768

        /**
         * Defines the relative size the chosen peak (pitch) has. 0.93 means: choose
         * the first peak that is higher than 93% of the highest peak detected. 93%
         * is the default value used in the Tartini user interface.
         */
        private const val DEFAULT_CUTOFF: Double = 0.97

        /**
         * For performance reasons, peaks below this cutoff are not even considered.
         */
        private const val SMALL_CUTOFF: Double = 0.5

        /**
         * Pitch annotations below this threshold are considered invalid, they are
         * ignored.
         */
        private const val LOWER_PITCH_CUTOFF: Double = 80.0 // Hz

    }

    /**
     * Contains a normalized square difference function value for each delay
     * (tau).
     */
    private val nsdf = FloatArray(audioBufferSize)

    /**
     * The probability of the last detected pitch.
     */
    private var probability = 0f

    /**
     * The x and y coordinate of the top of the curve (nsdf).
     */
    private var turningPointX: Float = 0f

    /**
     * The x and y coordinate of the top of the curve (nsdf).
     */
    private var turningPointY: Float = 0f

    /**
     * A list with minimum and maximum values of the nsdf curve.
     */
    private val maxPositions: ArrayList<Int> = ArrayList()

    /**
     * A list of estimates of the period of the signal (in samples).
     */
    private val periodEstimates: ArrayList<Float> = ArrayList()

    /**
     * A list of estimates of the amplitudes corresponding with the period
     * estimates.
     */
    private val ampEstimates: ArrayList<Float> = ArrayList()


    /**
     * Implements the normalized square difference function. See section 4 (and
     * the explanation before) in the MPM article. This calculation can be
     * optimized by using an FFT. The results should remain the same.
     *
     * @param audioBuffer
     *            The buffer with audio information.
     */
    private fun normalizedSquareDifference(audioBuffer: FloatArray) {
        for (tau in audioBuffer.indices) {
            var acf = 0f
            var divisorM = 0f
            for (i in 0 until audioBuffer.size - tau) {
                acf += audioBuffer[i] * audioBuffer[i + tau]
                divisorM += audioBuffer[i] * audioBuffer[i] + audioBuffer[i + tau] * audioBuffer[i + tau]
            }
            nsdf[tau] = 2 * acf / divisorM
        }
    }

    override fun detect(audioBuffer: FloatArray): PitchResult {
        var pitch = 0f

        // 0. Clear previous results (Is this faster than initializing a list
        // again and again?)
        maxPositions.clear()
        periodEstimates.clear()
        ampEstimates.clear()

        // 1. Calculate the normalized square difference for each Tau value.
        normalizedSquareDifference(audioBuffer)

        // 2. Peak picking time: time to pick some peaks.
        peakPicking()

        var highestAmplitude = Double.NEGATIVE_INFINITY

        for (tau in maxPositions) {
            // make sure every annotation has a probability attached
            highestAmplitude = max(highestAmplitude, nsdf[tau].toDouble())

            if (nsdf[tau] > SMALL_CUTOFF) {
                // Calculates turningPointX and Y
                prabolicInterpolation(tau)
                //store the turning points
                ampEstimates.add(turningPointY)
                periodEstimates.add(turningPointX)

                //Remember the highest amplitude
                highestAmplitude = max(highestAmplitude, turningPointY.toDouble())
            }
        }

        if (periodEstimates.isEmpty()) {
            pitch = -1f
        } else {
            // Use the overall maximum to calculate a cutoff.
            // The cutoff value is based on the highest value and a relative
            // threshold.
            val actualCutoff = cutOff * highestAmplitude

            // Find first period above or equal to cutoff
            var periodIndex = 0
            for (i in 0 until ampEstimates.size) {
                if(ampEstimates[i] > actualCutoff){
                    periodIndex = i
                    break
                }
            }

            val period = periodEstimates[periodIndex]
            val pitchEstimate = audioSampleRate / period
            if(pitchEstimate > LOWER_PITCH_CUTOFF){
                pitch = pitchEstimate
            }else{
                pitch = -1f
            }
        }
        probability = highestAmplitude.toFloat()

        return PitchResult.Pitch(
            hz = pitch,
            probability = probability
        )
    }

    /**
     * <p>
     * Finds the x value corresponding with the peak of a parabola.
     * </p>
     * <p>
     * a,b,c are three samples that follow each other. E.g. a is at 511, b at
     * 512 and c at 513; f(a), f(b) and f(c) are the normalized square
     * difference values for those samples; x is the peak of the parabola and is
     * what we are looking for. Because the samples follow each other
     * <code>b - a = 1</code> the formula for <a
     * href="http://fizyka.umk.pl/nrbook/c10-2.pdf">parabolic interpolation</a>
     * can be simplified a lot.
     * </p>
     * <p>
     * The following ASCII ART shows it a bit more clear, imagine this to be a
     * bit more curvaceous.
     * </p>
     *
     * <pre>
     *     nsdf(x)
     *       ^
     *       |
     * f(x)  |------ ^
     * f(b)  |     / |\
     * f(a)  |    /  | \
     *       |   /   |  \
     *       |  /    |   \
     * f(c)  | /     |    \
     *       |_____________________> x
     *            a  x b  c
     * </pre>
     *
     * @param tau
     *            The delay tau, b value in the drawing is the tau value.
     */
    fun prabolicInterpolation(tau: Int) {
        val nsdfa = nsdf[tau - 1]
        val nsdfb = nsdf[tau]
        val nsdfc = nsdf[tau + 1]
        val bValue = tau
        val bottom = nsdfc + nsdfa - 2 * nsdfb
        if (bottom == 0f) {
            turningPointX = bValue.toFloat()
            turningPointY = nsdfb
        } else {
            val delta = nsdfa - nsdfc
            turningPointX = bValue + delta / (2 * bottom)
            turningPointY = nsdfb - delta * delta / (8 * bottom)
        }
    }

    /**
     * <p>
     * Implementation based on the GPL'ED code of <a
     * href="http://tartini.net">Tartini</a> This code can be found in the file
     * <code>general/mytransforms.cpp</code>.
     * </p>
     * <p>
     * Finds the highest value between each pair of positive zero crossings.
     * Including the highest value between the last positive zero crossing and
     * the end (if any). Ignoring the first maximum (which is at zero). In this
     * diagram the desired values are marked with a +
     * </p>
     *
     * <pre>
     *  f(x)
     *   ^
     *   |
     *  1|               +
     *   | \      +     /\      +     /\
     *  0| _\____/\____/__\/\__/\____/_______> x
     *   |   \  /  \  /      \/  \  /
     * -1|    \/    \/            \/
     *   |
     * </pre>
     *
     * @param nsdf
     *            The array to look for maximum values in. It should contain
     *            values between -1 and 1
     * @author Phillip McLeod
     */
    fun peakPicking() {
        var pos = 0
        var curMaxPos = 0

        // find the first negative zero crossing
        while (pos < (nsdf.size - 1) / 3 && nsdf[pos] > 0) {
            pos++
        }

        //Loop over all the values below zero
        while (pos < nsdf.size - 1 && nsdf[pos] <= 0.0) {
            pos++
        }

        //Can happen if output[0] is NaN
        if (pos == 0) {
            pos = 1
        }

        while (pos < nsdf.size - 1) {
            check(nsdf[pos] >= 0)
            if (nsdf[pos] > nsdf[pos - 1] && nsdf[pos] >= nsdf[pos + 1]) {
                if (curMaxPos == 0) {
                    // The first max (between zero crossing)
                    curMaxPos = pos
                } else if (nsdf[pos] > nsdf[curMaxPos]) {
                    // A higher max (between the zero crossings)
                    curMaxPos = pos
                }
                pos++

                // A negative zero crossing
                if (pos < nsdf.size - 1 && nsdf[pos] <= 0) {
                    // If there was a maximum add it to the list of maxima
                    if (curMaxPos > 0) {
                        maxPositions.add(curMaxPos)
                        curMaxPos = 0 // Clear the maximum position, so we
                        //start looking for a new ones
                    }
                    while (pos < nsdf.size - 1 && nsdf[pos] <= 0f) {
                        pos++ //loop over all the values below zero
                    }
                }
                if (curMaxPos > 0) { // If there was a maximum in the last part
                    maxPositions.add(curMaxPos) // add it to the vector of maxima
                }
            }
        }
    }

}