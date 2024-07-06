package darkmat13r.tarsos.dsp

import darkmat13r.tarsos.dsp.io.AudioFloatConverter
import darkmat13r.tarsos.dsp.io.AudioFormat
import kotlin.math.log10
import kotlin.math.sqrt

/**
 * An audio event flows through the processing pipeline. The object is reused for performance reasons.
 * The arrays with audio information are also reused, so watch out when using the buffer getter and setters.
 *
 * @author Joren Six
 */
class AudioEvent(private val format: AudioFormat) {

    private val converter: AudioFloatConverter =
        AudioFloatConverter.Factory.create(format)

    /**
     * The audio block in floats
     * @return The float representation of the audio block.
     */
    /**
     * Set a new audio block.
     * @param floatBuffer The audio block that is passed to the next processor.
     */
    /**
     * The audio data encoded in floats from -1.0 to 1.0.
     */
    lateinit var floatBuffer: FloatArray
        private set

    /**
     * The audio data encoded in bytes according to format.
     */
    private var byteBuffer: ByteArray? = null

    /**
     * The overlap in samples between blocks of audio
     * @return The overlap in samples between blocks of audio
     */
    /**
     * Change the default overlap (in samples)
     * @param newOverlap The new overlap between audio blocks.
     */
    /**
     * The overlap in samples.
     */
    var overlap: Int = 0

    /**
     * The length of the stream, expressed in sample frames rather than byte
     * @return  The length of the stream, expressed in sample frames rather than bytes
     */
    /**
     * The length of the stream, expressed in sample frames rather than bytes
     */
    private val frameLength: Long = 0

    /**
     * The number of bytes processed before this event. It can be used to calculate the time stamp for when this event started.
     */
    private var bytesProcessed: Long = 0

    private var bytesProcessing = 0

    val sampleRate: Float
        /**
         * The audio sample rate in Hz
         * @return  The audio sample rate in Hz
         */
        get() = format.sampleRate

    val bufferSize: Int
        /**
         * The size of the buffer in samples.
         * @return The size of the buffer in samples.
         */
        get() = floatBuffer.size ?: -1

    /**
     * Change the number of bytes processed
     *
     * @param bytesProcessed The number of bytes processed.
     */
    fun setBytesProcessed(bytesProcessed: Long) {
        this.bytesProcessed = bytesProcessed
    }

    val timeStamp: Float
        /**
         * Calculates and returns the time stamp at the beginning of this audio event.
         * @return The time stamp at the beginning of the event in seconds.
         */
        get() = bytesProcessed / format.frameSize / format.sampleRate

    val endTimeStamp: Float
        /**
         * The timestamp at the end of the buffer (in seconds)
         * @return The timestamp at the end of the buffer (in seconds)
         */
        get() = (bytesProcessed + bytesProcessing) / format.frameSize / format.sampleRate

    val samplesProcessed: Long
        /**
         * The number of samples processed.
         * @return The number of samples processed.
         */
        get() = bytesProcessed / format.frameSize

    val progress: Double
        /**
         * Calculate the progress in percentage of the total number of frames.
         *
         * @return a percentage of processed frames or a negative number if the
         * number of frames is not known beforehand.
         */
        get() = bytesProcessed / format.frameSize / frameLength.toDouble()

    /**
     * Return a byte array with the audio data in bytes.
     * A conversion is done from float, cache accordingly on the other side...
     *
     * @return a byte array with the audio data in bytes.
     */
    fun getByteBuffer(): ByteArray {
        val length: Int = bufferSize * format.frameSize
        if (byteBuffer == null || byteBuffer!!.size != length) {
            byteBuffer = ByteArray(length)
        }
        converter.toByteArray(inBuff = floatBuffer, outBuff = byteBuffer!!)
        return byteBuffer!!
    }

    fun setAudioBuffer(buffer: FloatArray) {
        this.floatBuffer = buffer
    }

    val rMS: Double
        /**
         * Calculates and returns the root mean square of the signal. Please
         * cache the result since it is calculated every time.
         * @return The [RMS](http://en.wikipedia.org/wiki/Root_mean_square) of
         * the signal present in the current buffer.
         */
        get() = calculateRMS(floatBuffer)


    /**
     * Returns the dBSPL for a buffer.
     *
     * @return The dBSPL level for the buffer.
     */
    fun getdBSPL(): Double {
        return soundPressureLevel(floatBuffer)
    }

    /**
     * Set all sample values to zero.
     */
    fun clearFloatBuffer() {
        floatBuffer.fill(0f, 0, floatBuffer.size)
    }

    /**
     * Checks whether this block of audio is silent
     * @param silenceThreshold the threshold in spl to use.
     * @return True if SPL is below the threshold. False otherwise.
     */
    fun isSilence(silenceThreshold: Double): Boolean {
        return soundPressureLevel(floatBuffer) < silenceThreshold
    }

    /**
     * The number of bytes being processed.
     * @param bytesProcessing Sets the number of bytes being processed.
     */
    fun setBytesProcessing(bytesProcessing: Int) {
        this.bytesProcessing = bytesProcessing
    }

    companion object {
        /**
         * Calculates and returns the root mean square of the signal. Please
         * cache the result since it is calculated every time.
         * @param floatBuffer The audio buffer to calculate the RMS for.
         * @return The [RMS](http://en.wikipedia.org/wiki/Root_mean_square) of
         * the signal present in the current buffer.
         */
        fun calculateRMS(floatBuffer: FloatArray): Double {
            var rms = 0.0
            for (i in floatBuffer.indices) {
                rms += (floatBuffer[i] * floatBuffer[i]).toDouble()
            }
            rms /= floatBuffer.size
            rms = sqrt(rms)
            return rms
        }

        /**
         * Returns the dBSPL for a buffer.
         *
         * @param buffer
         * The buffer with audio information.
         * @return The dBSPL level for the buffer.
         */
        private fun soundPressureLevel(buffer: FloatArray): Double {
            val rms = calculateRMS(buffer)
            return linearToDecibel(rms)
        }

        /**
         * Converts a linear to a dB value.
         *
         * @param value
         * The value to convert.
         * @return The converted value.
         */
        private fun linearToDecibel(value: Double): Double {
            return 20.0 * log10(value)
        }
    }
}
