package darkmat13r.tarsos.dsp.io

import android.media.AudioRecord
import android.util.Log
import java.io.IOException

actual class NativeAudioInputStream(
    private val audioRecord : AudioRecord,
    override var format: AudioFormat,
) : AudioInputStream {

    companion object{
        private val TAG = NativeAudioInputStream::class.simpleName
    }

    override var frameLength: Long
        get() = -1
        set(value) {}

    override fun skip(bytesToSkip: Long): Long {
        throw IOException("Can not skip in audio stream")
    }

    override fun read(b: ByteArray, off: Int, len: Int): Int {
        Log.i(TAG, "Audio buffer reading")
       return audioRecord.read(b, off, len)
    }

    override fun close() {
        this.audioRecord.stop()
        this.audioRecord.release()
    }
}