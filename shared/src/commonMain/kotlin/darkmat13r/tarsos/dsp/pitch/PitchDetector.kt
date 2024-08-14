package darkmat13r.tarsos.dsp.pitch

/**
 * A pitch detector is capable of analyzing a buffer with audio information
 * and return a pitch estimation in Hz.
 */
abstract class PitchDetector {

    /**
     * Analyzes a buffer with audio information and estimates a pitch in Hz.
     * Currently, this interface only allows one pitch per buffer.
     *
     * The implementation of this function should process the audio data provided
     * in the buffer and return an estimation of the pitch. The pitch detection
     * algorithm may vary, but the result should include the detected pitch frequency
     * and, optionally, a probability indicating the confidence or clarity of the detected pitch.
     *
     * The audio buffer contains raw audio samples that are not modified by this function,
     * allowing the buffer to be reused for other types of audio analysis, such as FFT (Fast Fourier Transform).
     *
     * @param audioBuffer
     *            The buffer containing audio information. This buffer is a FloatArray
     *            where each element represents a single audio sample. The samples are
     *            typically normalized to the range [-1.0, 1.0].
     *
     * @return  An estimation of the pitch result. The returned PitchResult should contain:
     *          - A `freq` field representing the estimated pitch frequency in Hz.
     *          - An optional `probability` field which, if calculated by the algorithm, indicates
     *            the confidence or clarity of the detected pitch. This value can help distinguish
     *            between voiced and unvoiced segments or assess the reliability of the pitch detection.
     *
     * Example:
     *
     * val pitchResult = pitchDetector.detect(audioBuffer)
     * when (pitchResult) {
     *     is PitchResult.Pitch -> {
     *         val frequency = pitchResult.freq
     *         val confidence = pitchResult.probability
     *         // Use the frequency and confidence values as needed
     *     }
     *     else -> {
     *         // Handle other types of PitchResult if applicable
     *     }
     * }
     */
    abstract suspend fun detect(audioBuffer: FloatArray) : PitchResult
}