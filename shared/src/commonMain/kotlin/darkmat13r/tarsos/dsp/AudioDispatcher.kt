package darkmat13r.tarsos.dsp

import co.touchlab.kermit.Logger
import darkmat13r.tarsos.dsp.io.AudioFloatConverter
import darkmat13r.tarsos.dsp.io.AudioFormat
import darkmat13r.tarsos.dsp.io.AudioInputStream
import darkmat13r.tarsos.dsp.io.exception.IOException
import kotlin.math.round

class AudioDispatcher(
    /**
     * The audio stream (in bytes), conversion to float happens at the last
     * moment.
     */
    private val audioInputStream: AudioInputStream,
    private val audioBufferSize: Int,
    private val bufferOverlap: Int
) {

    /**
     * This buffer is reused again and again to store audio data using the float
     * data type.
     */
    private lateinit var audioFloatBuffer: FloatArray

    /**
     * This buffer is reused again and again to store audio data using the byte
     * data type.
     */
    private lateinit var audioByteBuffer: ByteArray

    /**
     * A list of registered audio processors. The audio processors are
     * responsible for actually doing the digital signal processing
     */
    private val audioProcessors: ArrayList<AudioProcessor> = arrayListOf()

    /**
     * Converter converts an array of floats to an array of bytes (and vice
     * versa).
     */
    private val converter: AudioFloatConverter

    private val format: AudioFormat = audioInputStream.format

    /**
     * The floatOverlap: the number of elements that are copied in the buffer
     * from the previous buffer. Overlap should be smaller (strict) than the
     * buffer size and can be zero. Defined in number of samples.
     */
    private var floatOverlap = 0
    private var floatStepSize: kotlin.Int = 0

    /**
     * The overlap and stepsize defined not in samples but in bytes. So it
     * depends on the bit depth. Since the int datatype is used only 8,16,24,...
     * bits or 1,2,3,... bytes are supported.
     */
    private var byteOverlap = 0
    private var byteStepSize: kotlin.Int = 0

    /**
     * The number of bytes to skip before processing starts.
     */
    private var bytesToSkip: Long = 0

    /**
     * Position in the stream in bytes. e.g. if 44100 bytes are processed and 16
     * bits per frame are used then you are 0.5 seconds into the stream.
     */
    private var bytesProcessed: Long = 0

    /**
     * The audio event that is send through the processing chain.
     */
    private val audioEvent: AudioEvent

    /**
     * If true the dispatcher stops dispatching audio.
     */
    private var stopped = false

    /**
     * If true then the first buffer is only filled up to buffer size - hop size
     * E.g. if the buffer is 2048 and the hop size is 48 then you get 2000 times
     * zero 0 and 48 actual audio samples. During the next iteration you get
     * mostly zeros and 96 samples.
     */
    private var zeroPadFirstBuffer = false

    /**
     * If true then the last buffer is zero padded. Otherwise the buffer is
     * shortened to the remaining number of samples. If false then the audio
     * processors must be prepared to handle shorter audio buffers.
     */
    private var zeroPadLastBuffer = false


    init {

        setStepSizeAndOverlap(audioBufferSize, bufferOverlap)
        audioEvent = AudioEvent(format)
        audioEvent.setAudioBuffer(audioFloatBuffer)
        audioEvent.overlap = bufferOverlap

        converter = AudioFloatConverter.Factory.create(format)

        stopped = false
        bytesToSkip = 0
        zeroPadLastBuffer = true
    }

    /**
     * Set a new step size and overlap size. Both in number of samples. Watch
     * out with this method: it should be called after a batch of samples is
     * processed, not during.
     *
     * @param audioBufferSize
     * The size of the buffer defines how much samples are processed
     * in one step. Common values are 1024,2048.
     * @param bufferOverlap
     * How much consecutive buffers overlap (in samples). Half of the
     * AudioBufferSize is common (512, 1024) for an FFT.
     */
    fun setStepSizeAndOverlap(audioBufferSize: Int, bufferOverlap: Int) {
        audioFloatBuffer = FloatArray(audioBufferSize)

        floatOverlap = bufferOverlap
        floatStepSize = audioFloatBuffer.size - floatOverlap

        audioByteBuffer = ByteArray(audioFloatBuffer.size * format.frameSize)
        byteOverlap = floatOverlap * format.frameSize
        byteStepSize = floatStepSize * format.frameSize
    }

    /**
     * Skip a number of seconds before processing the stream.
     * @param seconds The number of seconds to skip
     */
    fun skip(seconds: Double) {
        bytesToSkip = (round(seconds * format.sampleRate) * format.frameSize).toLong()
    }

    /**
     * if zero pad is true then the first buffer is only filled up to  buffer size - hop size
     * E.g. if the buffer is 2048 and the hop size is 48 then you get 2000x0 and 48 filled audio samples
     * @param zeroPadFirstBuffer true if the buffer should be zeroPadFirstBuffer, false otherwise.
     */
    fun setZeroPadFirstBuffer(zeroPadFirstBuffer: Boolean) {
        this.zeroPadFirstBuffer = zeroPadFirstBuffer
    }

    /**
     * If zero pad last buffer is true then the last buffer is filled with zeros until the normal amount
     * of elements are present in the buffer. Otherwise the buffer only contains the last elements and no zeros.
     * By default it is set to true.
     *
     * @param zeroPadLastBuffer A boolean to control whether the last buffer is zero-padded.
     */
    fun setZeroPadLastBuffer(zeroPadLastBuffer: Boolean) {
        this.zeroPadLastBuffer = zeroPadLastBuffer
    }

    /**
     * Adds an AudioProcessor to the chain of processors.
     *
     * @param audioProcessor
     * The AudioProcessor to add.
     */
    fun addAudioProcessor(audioProcessor: AudioProcessor) {
        audioProcessors.add(audioProcessor)
        Logger.i("Added an audioprocessor to the list of processors: $audioProcessor")
    }

    /**
     * Removes an AudioProcessor to the chain of processors and calls its `processingFinished` method.
     *
     * @param audioProcessor
     * The AudioProcessor to remove.
     */
    fun removeAudioProcessor(audioProcessor: AudioProcessor) {
        audioProcessors.remove(audioProcessor)
        audioProcessor.processingFinished()
        Logger.i("Remove an audioprocessor to the list of processors: $audioProcessor")
    }

    fun run() {
        var bytesToRead = 0

        if (bytesToSkip != 0L) {
            skipToStart()
        }

        //Read the first (and in some cases last) audio block
        try {
            //Need to get correct time info when skipping first x seconds
            audioEvent.setBytesProcessed(bytesProcessed)
            bytesToRead = readNextAudioBlock()
        } catch (e: IOException) {
            val message = "Error while reading audio input stream: " + e.message
            Logger.w(message)
            throw Error(message)
        }

        Logger.i("bytesToRead ${bytesToRead}")

        // As long as the stream has not ended
        while (bytesToRead != 0 && !stopped) {

            Logger.i("bytesToRead iteration ${bytesToRead}")
            //Makes sure the right buffers are processed, they can be changed by audio processors.
            try{
                for (processor in audioProcessors) {
                    if (!processor.process(audioEvent)) {
                        //skip to the next audio processors if false is returned.
                        break
                    }
                }
            }catch (e: Exception){
                e.printStackTrace()
            }

            if (!stopped) {
                //Update the number of bytes processed
                bytesProcessed += bytesToRead
                audioEvent.setBytesProcessed(bytesProcessed)

                // Read, convert and process consecutive  overlapping buffers
                // Slide the buffer.
                try {
                    bytesToRead = readNextAudioBlock()
                    Logger.i("bytesToRead  readNextAudioBlock ${bytesToRead}")
                    audioEvent.overlap = floatOverlap
                } catch (e: IOException) {
                    val message = "Error while reading audio input stream: " + e.message
                    Logger.w(message)
                    throw Error(message)
                }
            }
        }
        Logger.i("stopped ${stopped}")

        // Notify all processors that no more data is available.
        // when stop() is called processingFinished is called explicitly, no need to do this again.
        // The explicit call is to prevent timing issues.
        if (!stopped) {
            stop()
        }
    }

    /**
     * Reads the next audio block. It tries to read the number of bytes defined
     * by the audio buffer size minus the overlap. If the expected number of
     * bytes could not be read either the end of the stream is reached or
     * something went wrong.
     *
     * The behavior for the first and last buffer is defined by their corresponding the zero pad settings. The method also handles the case if
     * the first buffer is also the last.
     *
     * @return The number of bytes read.
     * @throws IOException
     *             When something goes wrong while reading the stream. In
     *             particular, an IOException is thrown if the input stream has
     *             been closed.
     */
    private fun readNextAudioBlock(): Int {
        Logger.i("readNextAudioBlock ${floatOverlap} < ${audioFloatBuffer.size}")
        check(floatOverlap < audioFloatBuffer.size)

        // Is this the first buffer?
        val isFirstBuffer = (bytesProcessed == 0L || bytesProcessed == bytesToSkip)

        var offsetInBytes: Int = 0
        var offsetInSamples: Int = 0
        var bytesToRead: Int = 0

        // Determine the amount of bytes to read from the stream
        if (isFirstBuffer && zeroPadFirstBuffer) {
            //If this is the first buffer and we do not want to zero page the
            //first buffer then read a full buffer
            bytesToRead = audioByteBuffer.size
            //With an offset in bytes of zero
            offsetInBytes = 0
            offsetInSamples = 0
        } else {
            // In all other cases read the amount of bytes defined by the step size
            bytesToRead = byteStepSize
            offsetInBytes = byteOverlap
            offsetInSamples = floatOverlap
        }

        // Shift the audio information using copyTo since it is probably faster than
        // manually shifting it
        if (!isFirstBuffer && audioFloatBuffer.size == floatOverlap + floatStepSize) {
            audioFloatBuffer.copyInto(
                audioFloatBuffer,
                destinationOffset = 0,
                startIndex = floatStepSize,
                endIndex = floatStepSize + floatOverlap
            )
        }

        //Total amount of bytes read
        var totalBytesRead = 0

        // The amount of bytes read from the stream during one iteration
        var bytesRead = 0

        // Is the end of the stream reached?
        var endOfStream = false

        Logger.i("Always try to read the `bytesToRead` amount of bytes")
        // Always try to read the `bytesToRead` amount of bytes
        // Unless the stream is closed (stopped is true) or no bytes could be read during
        //on iteration
        Logger.i("stopped ${stopped},endOfStream ${endOfStream} , totalBytesRead ${totalBytesRead}bytesToRead ${bytesToRead}")
        while (!stopped && !endOfStream && totalBytesRead < bytesToRead) {
            /*try {
                bytesRead = audioInputStream.read(
                    audioByteBuffer,
                    offsetInBytes + totalBytesRead,
                    bytesToRead - totalBytesRead
                )
            } catch (ex: IndexOutOfBoundsException) {
                // The pipe decoder generates an out of bounds if end
                // of stream is reached. Ugly hack...
                bytesRead = -1
                ex.printStackTrace()
            }*/
            bytesRead = audioInputStream.read(
                audioByteBuffer,
                offsetInBytes + totalBytesRead,
                bytesToRead - totalBytesRead
            )

            if (bytesRead == -1) {
                // The end of the stream is reached if the number of bytes read
                // during this iteration equals -1
                endOfStream = true
            } else {
                //Otherwise add the number of bytes read to the total
                totalBytesRead += bytesRead
            }
        }

        Logger.i("endOfStream ${endOfStream}")
        if (endOfStream) {
            // Could not read a full buffer from the stream , there are two options:
            if (zeroPadLastBuffer) {
                for (i in offsetInBytes + totalBytesRead until audioByteBuffer.size) {
                    audioByteBuffer[i] = 0
                }
                converter.toFloatArray(
                    inBuff = audioByteBuffer,
                    inOffset = offsetInBytes,
                    outBuff = audioFloatBuffer,
                    outOffset = offsetInSamples,
                    outLen = floatStepSize
                )
            } else {
                // Send a smaller buffer through the chain
                val audioByteBufferContent = audioByteBuffer
                audioByteBuffer = ByteArray(offsetInBytes + totalBytesRead)
                for (i in audioFloatBuffer.indices) {
                    audioByteBuffer[i] = audioByteBufferContent[i]
                }

                val totalSamplesRead = totalBytesRead / format.frameSize
                audioFloatBuffer = FloatArray(offsetInSamples + totalBytesRead / format.frameSize)
                converter.toFloatArray(
                    inBuff = audioByteBuffer,
                    inOffset = offsetInBytes,
                    outBuff = audioFloatBuffer,
                    outOffset = offsetInSamples,
                    outLen = totalSamplesRead
                )
            }
        } else if (bytesToRead == totalBytesRead) {
            // The expected amount of bytes have been read from the stream
            if (isFirstBuffer && !zeroPadFirstBuffer) {
                converter.toFloatArray(
                    inBuff = audioByteBuffer,
                    outBuff = audioFloatBuffer,
                )
            } else {
                converter.toFloatArray(
                    inBuff = audioByteBuffer,
                    inOffset = offsetInBytes,
                    outBuff = audioFloatBuffer,
                    outOffset = offsetInSamples,
                    outLen = floatStepSize
                )
            }
        } else if (!stopped) {
            // If the end of the stream has not been reached and the number of bytes read is not the
            // Expected amount of bytes ,  then we are in an invalid state;
            throw IOException(
                "The end of the audio stream has not been reached and the number of bytes read ($totalBytesRead) is not equal "
                        + "to the expected amount of bytes($bytesToRead)."
            )
        }

        //Makes sure AudioEvent contains correct info.
        audioEvent.setAudioBuffer(audioFloatBuffer)
        audioEvent.overlap = offsetInSamples
        Logger.i("totalBytesRead ${totalBytesRead}")
        return totalBytesRead
    }

    private fun skipToStart() {
        var skipped = 0L
        try {
            skipped = audioInputStream.skip(bytesToSkip)
            if (skipped != bytesToSkip) {
                throw IOException()
            }
            bytesProcessed += bytesToSkip
        } catch (e: IOException) {
            val message: String =
                "Did not skip the expected amount of bytes, " +
                        " $skipped skipped, $bytesToSkip expected!"
            Logger.w(message)
            throw Error(message)
        }
    }

    /**
     * Stops dispatching audio data.
     */
    fun stop() {
        stopped = true
        for (processor in audioProcessors) {
            processor.processingFinished()
        }
        try {
            audioInputStream.close()
        } catch (e: IOException) {
            Logger.e(
                "Closing audio stream error.",
                e
            )
        }
    }


    /**
     *
     * @return The currently processed number of seconds.
     */
    fun secondsProcessed(): Float {
        return bytesProcessed / (format.sampleSizeInBits / 8) / format.sampleRate / format.channels
    }

    /**
     * Set a new audio buffer
     * @param audioBuffer The audio buffer to use.
     */
    fun setAudioFloatBuffer(audioBuffer: FloatArray) {
        audioFloatBuffer = audioBuffer
    }

    /**
     * @return True if the dispatcher is stopped or the end of stream has been reached.
     */
    fun isStopped(): Boolean {
        return stopped
    }
}