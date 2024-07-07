package audio

import co.touchlab.kermit.Logger
import darkmat13r.tarsos.dsp.AudioDispatcher
import darkmat13r.tarsos.dsp.io.NativeAudioInputStream

actual class AudioManager {
    private val format = darkmat13r.tarsos.dsp.io.AudioFormat(
        AudioConfig.SAMPLE_RATE.toFloat(),
        AudioConfig.SAMPLE_RATE_BITS,
        AudioConfig.CHANNEL_COUNT,
        signed = true,
        bigEndian = false
    )
    private val audioStream: NativeAudioInputStream by lazy {
        NativeAudioInputStream(format)
    }

    actual fun createAudioDispatcher(): AudioDispatcher {
        val bufferSize = audioStream.getBufferSize()
        return AudioDispatcher(audioStream, bufferSize, AudioConfig.OVERLAP)
    }

    actual fun getBufferSize(): Int {
        return audioStream.getBufferSize()
    }

    actual fun stop() {
        audioStream.close()
    }

}