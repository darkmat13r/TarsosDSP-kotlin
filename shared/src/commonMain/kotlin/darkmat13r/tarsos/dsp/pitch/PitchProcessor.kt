package darkmat13r.tarsos.dsp.pitch

import darkmat13r.tarsos.dsp.AudioEvent
import darkmat13r.tarsos.dsp.AudioProcessor

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
    private val sampleRate: Float,
    private val bufferSize: Int,
    private val handler: DetectedPitchHandler,
    private val algorithm: PitchEstimationAlgorithm? = null,
    pitchDetector: PitchDetector? = null,
) : AudioProcessor {

    private val detector: PitchDetector

    private var processedSamples: Int = 0

    init {
        detector = when (algorithm) {
            PitchEstimationAlgorithm.YIN -> Yin(sampleRate, bufferSize)
            PitchEstimationAlgorithm.MPM -> McLeodPitchMethod(sampleRate, bufferSize)
            else -> pitchDetector ?: throw Exception("Either set algorithm or set pitch detector")
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
        fun handlePitch(pitch: PitchResult, audioEvent: AudioEvent)
    }

    override suspend fun process(audioEvent: AudioEvent): Boolean {
        val audioFloatBuffer: FloatArray = audioEvent.floatBuffer
        val pitch = detector.detect(audioFloatBuffer)
        handler.handlePitch(pitch, audioEvent)
        return true
    }


    override fun processingFinished() {
        //DO NOTHING
    }
}