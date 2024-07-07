package audio

import darkmat13r.tarsos.dsp.AudioDispatcher
import darkmat13r.tarsos.dsp.io.NativeAudioInputStream

actual class AudioManager {
    private  var audioStream: NativeAudioInputStream? = null

    actual fun createAudioDispatcher(): AudioDispatcher {
        val format = darkmat13r.tarsos.dsp.io.AudioFormat(
            AudioConfig.SAMPLE_RATE.toFloat(),
            AudioConfig.SAMPLE_RATE_BITS,
            AudioConfig.CHANNEL_COUNT,
            signed = true,
            bigEndian = false
        )
        audioStream = NativeAudioInputStream(format)
        val bufferSize = audioStream!!.getBufferSize()
        return AudioDispatcher(audioStream!!, bufferSize, AudioConfig.OVERLAP)
    }

    actual fun getBufferSize(): Int {
        return audioStream?.getBufferSize() ?: 0
    }

    actual fun stop() {
        audioStream?.close()
    }

}