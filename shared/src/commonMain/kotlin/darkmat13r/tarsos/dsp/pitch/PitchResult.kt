package darkmat13r.tarsos.dsp.pitch

sealed interface PitchResult{
    data class Pitch(
        /**
         * pitch in Hz
         */
        val hz : Float,
        /**
         * Some algorithms can calculate a probability (noisiness, (a)periodicity,
         * salience, voicedness or clarity measure) for the detected pitch. This is
         * somewhat similar to the term voiced which is used in speech recognition.
         * This probability should be calculated together with the pitch but is
         * returned using a call to this method. So if you want the probability of a
         * buffer: first call getPitch(buffer) and then getProbability().
         * @return A probability
         */
        val probability : Float = -1f
    ) : PitchResult
}