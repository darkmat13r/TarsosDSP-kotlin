package darkmat13r.tarsos.dsp.io

import darkmat13r.tarsos.dsp.io.exception.IOException

/**
 * Decouples the audio input stream
 * @author Joren Six
 */
interface AudioInputStream {
    /**
     * Skip a number of bytes before reading the remaining bytes.
     * @param bytesToSkip The number of bytes to skip.
     * @return The number of bytes skipped.
     * @throws IOException If the underlying  if an input or output error occurs
     * #see read
     */
    @Throws(IOException::class)
    fun skip(bytesToSkip: Long): Long

    /**
     * Reads up to a specified maximum number of bytes of data from the audio
     * stream, putting them into the given byte array.
     *
     * This method will always read an integral number of frames.
     * If `len` does not specify an integral number
     * of frames, a maximum of `len - (len % frameSize)
    ` *  bytes will be read.
     *
     * @param b the buffer into which the data is read
     * @param off the offset, from the beginning of array `b`, at which
     * the data will be written
     * @param len the maximum number of bytes to read
     * @return the total number of bytes read into the buffer, or -1 if there
     * is no more data because the end of the stream has been reached
     * @throws IOException if an input or output error occurs
     * @see .skip
     */
    @Throws(IOException::class)
    fun read(b: ByteArray, off: Int, len: Int): Int

    /**
     * Closes this audio input stream and releases any system resources associated
     * with the stream.
     * @throws IOException if an input or output error occurs
     */
    @Throws(IOException::class)
    fun close()

    /**
     *
     * @return The format of the underlying audio
     */
    var format: AudioFormat

    var frameLength: Long
}
