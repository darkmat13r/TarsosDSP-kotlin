package darkmat13r.tarsos.dsp.io

import android.media.AudioRecord
import java.io.IOException

actual class NativeAudioInputStream(
    private val audioRecord : AudioRecord,
    override var format: AudioFormat,
) : AudioInputStream {

    override var frameLength: Long
        get() = -1
        set(value) {}

    override fun skip(bytesToSkip: Long): Long {
        throw IOException("Can not skip in audio stream")
    }

    override fun read(b: ByteArray, off: Int, len: Int): Int {
       return audioRecord.read(b, off, len)
    }

    override fun close() {
        this.audioRecord.stop()
        this.audioRecord.release()
    }
}