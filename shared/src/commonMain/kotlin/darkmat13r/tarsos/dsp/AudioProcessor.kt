package darkmat13r.tarsos.dsp

interface AudioProcessor {
    /**
     * Process the first (complete) buffer. Once the first complete buffer is
     * processed the remaining buffers are overlapping buffers and processed
     * using the processOverlapping method (Even if overlap is zero).
     * @param audioEvent
     *            The audio event that contains audio data.
     * @return False if the chain needs to stop here, true otherwise. This can be used to implement e.g. a silence detector.
     */
    suspend fun process(audioEvent: AudioEvent) : Boolean


    /**
     * Notify the AudioProcessor that no more data is available and processing
     * has finished. Can be used to deallocate resources or cleanup.
     */
    fun processingFinished()
}