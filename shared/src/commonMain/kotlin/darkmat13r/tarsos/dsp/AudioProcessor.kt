package darkmat13r.tarsos.dsp

interface AudioProcessor {
    /**
     * Process the first (complete) buffer. Once the first complete buffer is
     * processed the remaining buffers are overlapping buffers and processed
     * using the processOverlapping method (Even if overlap is zero).
     * @param audioFloatBuffer
     *            The buffer to process using the float data type.
     * @param audioByteBuffer
     *            The buffer to process using raw bytes.
     * @return False if the chain needs to stop here, true otherwise. This can be used to implement e.g. a silence detector.
     */
    fun processFull(audioFloatBuffer : FloatArray, audioByteBuffer : ByteArray) : Boolean

    /**
     * Do the actual signal processing on an overlapping buffer. Once the
     * first complete buffer is processed the remaining buffers are
     * overlapping buffers and are processed using the processOverlapping
     * method. Even if overlap is zero.
     * @param audioFloatBuffer
     *            The buffer to process using the float data type.
     * @param audioByteBuffer
     *            The buffer to process using raw bytes.
     */
    fun processOverlapping(audioFloatBuffer: FloatArray, audioByteBuffer : ByteArray) : Boolean

    /**
     * Notify the AudioProcessor that no more data is available and processing
     * has finished. Can be used to deallocate resources or cleanup.
     */
    fun processingFinished()
}