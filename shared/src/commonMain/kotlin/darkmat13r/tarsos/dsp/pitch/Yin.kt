package darkmat13r.tarsos.dsp.pitch

class Yin(
    private val audioSampleRate: Float,
    private val bufferSize: Int,
    private val threshold: Float = DEFAULT_THRESHOLD
) : PitchDetector() {

    private val yinBuffer: FloatArray = FloatArray(bufferSize / 2)
    private var result = PitchResult.Pitch()


    companion object {
        private const val DEFAULT_THRESHOLD = 0.20f
        private const val DEFAULT_OVERLAP = 512
        private const val DEFAULT_BUFFER_SIZE = 1024
    }

    override fun detect(audioBuffer: FloatArray): PitchResult {
        var tauEstimate = 0

        // Step 2
        difference(audioBuffer)

        //Step 3
        cumulativeMeanNormalizedDifference()

        //Step 4
        tauEstimate = absoluteThreshold()

        //Step 5
        if (tauEstimate != -1) {
            var betterTau: Int = parabolicInterpolation(tauEstimate)


            // step 6
            // TODO Implement optimization for the AUBIO_YIN algorithm.
            // 0.77% => 0.5% error rate,
            // using the data of the YIN paper
            // bestLocalEstimate()

            // conversion to Hz
            result.hz = audioSampleRate / betterTau
        } else {
            // no pitch found
            result.hz = -1f
        }
        return result
    }

    /**
     * Implements the difference function as described in step 2 of the YIN
     * paper.
     */
    private fun difference(audioBuffer: FloatArray) {
        var delta = 0f
        for (tau in yinBuffer.indices) {
            yinBuffer[tau] = 0f
        }
        for (tau in 1 until yinBuffer.size) {
            for (i in yinBuffer.indices) {
                delta = audioBuffer[i] - audioBuffer[i + tau]
                yinBuffer[tau] += delta * delta
            }
        }
    }

    /**
     * The cumulative mean normalized difference function as described in step 3
     * of the YIN paper. <br>
     * <code>
     * yinBuffer[0] == yinBuffer[1] = 1
     * </code>
     */
    private fun cumulativeMeanNormalizedDifference() {
        yinBuffer[0] = 1f
        var runningSum = 0f
        for (tau in 1 until yinBuffer.size) {
            runningSum += yinBuffer[tau]
            yinBuffer[tau] *= tau / runningSum
        }
    }

    private fun absoluteThreshold(): Int {
        var tau = 2
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
                result.probability = 1f - yinBuffer[tau]
                break
            }
            tau++
        }


        // if no pitch found, tau => -1
        if (tau == yinBuffer.size || yinBuffer[tau] >= threshold) {
            tau = -1
            result.probability = 0f
            result.pitched = false
        }else{
            result.pitched = true
        }
        return tau
    }


    private fun parabolicInterpolation(tauEstimate: Int): Int {
        var betterTau = 0
        var x0 = 0
        var x2 = 0
        if (tauEstimate < 1) {
            x0 = tauEstimate
        } else {
            x0 = tauEstimate - 1
        }

        if (tauEstimate + 1 < yinBuffer.size) {
            x2 = tauEstimate + 1
        } else {
            x2 = tauEstimate
        }

        if (x0 == tauEstimate) {
            if (yinBuffer[tauEstimate] <= yinBuffer[x2]) {
                betterTau = tauEstimate
            } else {
                betterTau = x2
            }
        } else if (x2 == tauEstimate) {
            if (yinBuffer[tauEstimate] <= yinBuffer[x0]) {
                betterTau = tauEstimate
            } else {
                betterTau = x0
            }
        } else {
            val s0 = yinBuffer[x0]
            val s1 = yinBuffer[tauEstimate]
            val s2 = yinBuffer[x2]


            // fixed AUBIO implementation, thanks to Karl Helgason:
            // (2.0f * s1 - s2 - s0) was incorrectly multiplied with -1
            betterTau = (tauEstimate + (s2 - s0) / (2 * (2 * s1 - s2 - s0))).toInt()


        }
        return betterTau
    }
}