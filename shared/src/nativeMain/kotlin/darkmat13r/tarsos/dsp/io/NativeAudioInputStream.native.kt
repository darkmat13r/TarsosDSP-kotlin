package darkmat13r.tarsos.dsp.io

import co.touchlab.kermit.Logger
import darkmat13r.tarsos.dsp.io.exception.IOException
import kotlinx.cinterop.ByteVar
import kotlinx.cinterop.ExperimentalForeignApi
import kotlinx.cinterop.convert
import kotlinx.cinterop.pointed
import kotlinx.cinterop.refTo
import kotlinx.coroutines.delay
import kotlinx.coroutines.runBlocking
import platform.AVFAudio.AVAudioEngine
import platform.AVFAudio.AVAudioFormat
import platform.AVFAudio.AVAudioSession
import platform.AVFAudio.AVAudioSessionCategoryPlayAndRecord
import platform.AVFAudio.sampleRate
import platform.AVFAudio.setActive
import platform.posix.memcpy
import kotlin.concurrent.AtomicReference

actual class NativeAudioInputStream(
    override var format: AudioFormat,
) : AudioInputStream {

    private val audioEngine = AVAudioEngine()
    private val inputNode = audioEngine.inputNode
    @OptIn(ExperimentalForeignApi::class)
    private val audioFormat = AVAudioFormat(standardFormatWithSampleRate = format.sampleRate.toDouble(), channels =  format.channels.convert()) ?: error("Unable to create AVAudioFormat")

    private val bufferQueue = AtomicReference<ArrayDeque<ByteArray>>(ArrayDeque())
    private var isRunning = true

    init {
        setup()
    }

    @OptIn(ExperimentalForeignApi::class)
    private fun setup(){
        val bufferSize = getBufferSize()
        val session = AVAudioSession.sharedInstance()
        session.setCategory(AVAudioSessionCategoryPlayAndRecord, error = null)
        session.setActive(true, error = null)
        inputNode.installTapOnBus(0u, bufferSize.convert(), audioFormat) { buffer, time ->
            if (!isRunning) return@installTapOnBus
            if (buffer == null) return@installTapOnBus

            val channelData = buffer.audioBufferList?.pointed?.mBuffers?.pointed?.mData
            val dataLength = buffer.audioBufferList?.pointed?.mBuffers?.pointed?.mDataByteSize?.toInt() ?: 0

            if (channelData != null && dataLength > 0) {
                val byteArray = ByteArray(dataLength)
                memcpy(byteArray.refTo(0), channelData, dataLength.convert())
                bufferQueue.value = bufferQueue.value.apply { add(byteArray) }
            }
        }

        audioEngine.prepare()
        audioEngine.startAndReturnError(null)
    }

    override fun skip(bytesToSkip: Long): Long {
        throw IOException("Cannot skip in audio stream")
    }

    override fun read(b: ByteArray, off: Int, len: Int): Int {
        while (bufferQueue.value.isEmpty() && isRunning) {
            runBlocking { delay(10) }
        }

        val buffer = bufferQueue.value.removeFirstOrNull() ?: return 0
        val bytesRead = minOf(buffer.size, len)
        buffer.copyInto(b, off, 0, bytesRead)
        return bytesRead
    }

    override fun close() {
        isRunning = false
        inputNode.removeTapOnBus(0u)
        audioEngine.stop()
    }

    override var frameLength: Long
        get() = -1
        set(value) {}

    fun getBufferSize() = format.frameSize * format.sampleSizeInBits
}