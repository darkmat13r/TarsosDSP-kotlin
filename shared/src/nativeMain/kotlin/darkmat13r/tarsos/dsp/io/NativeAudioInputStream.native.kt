package darkmat13r.tarsos.dsp.io

import darkmat13r.tarsos.dsp.io.exception.IOException
import platform.AVFAudio.AVAudioEngine

actual class NativeAudioInputStream(
    val audioEngine : AVAudioEngine,
    override var format: AudioFormat,
) : AudioInputStream {
    override fun skip(bytesToSkip: Long): Long {
        throw IOException("Can not skip in audio stream")
    }

    override fun read(b: ByteArray, off: Int, len: Int): Int {
        TODO("Not yet implemented")
    }

    override fun close() {
        TODO("Not yet implemented")
    }

    override var frameLength: Long
        get() = -1
        set(value) {}
}