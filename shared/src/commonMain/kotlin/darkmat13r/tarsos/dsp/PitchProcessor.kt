package darkmat13r.tarsos.dsp

import darkmat13r.tarsos.dsp.pitch.McLeodPitchMethod
import darkmat13r.tarsos.dsp.pitch.PitchDetector
import darkmat13r.tarsos.dsp.pitch.PitchResult
import darkmat13r.tarsos.dsp.pitch.Yin

/**
 * Initialize a new pitch processor.
 *
 * @param algorithm
 *            An enum defining the algorithm.
 * @param sampleRate
 *            The sample rate of the buffer (Hz).
 * @param bufferSize
 *            The size of the buffer in samples.
 * @param bufferOverlap
 *            The size of the overlap between two consecutive buffers (in
 *            samples).
 * @param totalLengthInSamples
 *            The total length of the stream (in samples).
 * @param handler
 *            The handler handles detected pitch.
 */
class PitchProcessor(
    private val algorithm: PitchEstimationAlgorithm,
    private val sampleRate: Float,
    private val bufferSize: Int,
    private val bufferOverlap: Int,
    private val totalLengthInSamples: Long,
    private val handler: DetectedPitchHandler
) : AudioProcessor {

    private val detector : PitchDetector

    private var processedSamples: Int = 0

    init {
       detector = when(algorithm){
            PitchEstimationAlgorithm.YIN -> McLeodPitchMethod(sampleRate, bufferSize)
            PitchEstimationAlgorithm.MPM -> Yin(sampleRate, bufferSize)
        }
    }

    /**
     * A list of pitch estimation algorithms.
     */
    enum class PitchEstimationAlgorithm {
        /**
         * See {@link Yin} for the implementation. Or see <a href=
         * "http://recherche.ircam.fr/equipes/pcm/cheveign/ps/2002_JASA_YIN_proof.pdf"
         * >the YIN article</a>.
         */
        YIN,

        /**
         * See {@link McLeodPitchMethod}. It is described in the article "<a
         * href=
         * "http://miracle.otago.ac.nz/postgrads/tartini/papers/A_Smarter_Way_to_Find_Pitch.pdf"
         * >A Smarter Way to Find Pitch</a>".
         */
        MPM
    }

    /**
     * An interface to handle detected pitch.
     *
     */
    interface DetectedPitchHandler {
        /**
         * Handle a detected pitch.
         *
         * @param pitch
         *            The pitch in Hz. -1 is returned if no pitch is detected.
         * @param timestamp
         *            A time stamp associated with the detection.
         * @param progress
         *            If the length of the stream is known beforehand (a file) a
         *            progress indication is possible. It is a percentage. If a
         *            stream is analyzed a negative value is returned.
         */
        fun handlePitch(pitch: PitchResult, timestamp: Float, progress: Float)
    }

    override fun processFull(audioFloatBuffer: FloatArray, audioByteBuffer: ByteArray): Boolean {
        processedSamples += audioByteBuffer.size
        val pitch = detector.detect(audioFloatBuffer)
        val timestamp = processedSamples / sampleRate
        val progress = processedSamples / totalLengthInSamples
        handler.handlePitch(pitch, timestamp, progress.toFloat())
        return true
    }

    override fun processOverlapping(
        audioFloatBuffer: FloatArray,
        audioByteBuffer: ByteArray
    ): Boolean {
        processedSamples -= bufferOverlap
        return processFull(audioFloatBuffer, audioByteBuffer)
    }

    override fun processingFinished() {
        //DO NOTHING
    }
}