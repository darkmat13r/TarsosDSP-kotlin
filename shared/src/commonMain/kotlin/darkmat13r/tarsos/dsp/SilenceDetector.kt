package darkmat13r.tarsos.dsp

import io.ktor.utils.io.ByteReadChannel
import kotlin.math.log10
import kotlin.math.pow

/**
 * The silence detector breaks the audio processing pipeline when silence is detected, when there is sound
 * (no silence) it hands over the audio data to the next AudioProcessor.
 */
class SilenceDetector(val threshold: Double = DEFAULT_SILENCE_THRESHOLD) : AudioProcessor {
    companion object {
        private const val DEFAULT_SILENCE_THRESHOLD = -70.0  //db
    }

    /**
     * Calculates the local (linear) energy of an audio buffer.
     *
     * @param buffer
     *            The audio buffer.
     * @return The local (linear) energy of an audio buffer.
     */
    private fun localEnergy(buffer: FloatArray): Double {
        var power: Double = 0.0
        for (element in buffer) {
            power += element * element
        }
        return power
    }

    /**
     * Returns the dBSPL for a buffer.
     *
     * @param buffer
     *            The buffer with audio information.
     * @return The dBSPL level for the buffer.
     */
    private fun soundPressureLeve(buffer: FloatArray): Double {
        var value = localEnergy(buffer).pow(0.5)
        value /= buffer.size
        return linearToDecibel(value)
    }

    /**
     * Converts a linear to a dB value.
     *
     * @param value
     *            The value to convert.
     * @return The converted value.
     */
    private fun linearToDecibel(value: Double): Double = 20 * log10(value)

    /**
     * Checks if the dBSPL level in the buffer falls below a certain threshold.
     *
     * @param buffer
     *            The buffer with audio information.
     * @param silenceThreshold
     *            The threshold in dBSPL
     * @return True if the audio information in buffer corresponds with silence,
     *         false otherwise.
     */
    fun isSilence(buffer: FloatArray, silenceThreshold: Double = threshold): Boolean =
        soundPressureLeve(buffer) < silenceThreshold


    override fun processFull(audioFloatBuffer: FloatArray, audioByteBuffer: ByteArray): Boolean {
        return !isSilence(audioFloatBuffer)
    }

    override fun processOverlapping(
        audioFloatBuffer: FloatArray,
        audioByteBuffer: ByteArray
    ): Boolean {
        return processFull(audioFloatBuffer, audioByteBuffer)
    }

    override fun processingFinished() {
        //DO NOTHING
    }
}