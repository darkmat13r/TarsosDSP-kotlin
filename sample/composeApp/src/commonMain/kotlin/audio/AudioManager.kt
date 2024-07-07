package audio

import darkmat13r.tarsos.dsp.AudioDispatcher

object AudioConfig{
    const val SAMPLE_RATE = 22050
    const val SAMPLE_RATE_BITS = 16
    const val CHANNEL_COUNT = 1
    const val OVERLAP = 0
    const val BUFFER_SIZE = 2048
}

expect class AudioManager() {
    fun getBufferSize() : Int
    fun createAudioDispatcher(): AudioDispatcher
    fun stop()
}