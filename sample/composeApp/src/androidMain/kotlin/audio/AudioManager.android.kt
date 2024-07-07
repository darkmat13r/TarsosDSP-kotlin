package audio

import android.Manifest
import android.annotation.SuppressLint
import android.content.pm.PackageManager
import android.media.AudioFormat
import android.media.AudioManager
import android.media.AudioRecord
import android.media.MediaRecorder
import android.media.audiofx.NoiseSuppressor
import android.util.Log
import androidx.core.app.ActivityCompat
import darkmat13r.tarsos.dsp.AudioDispatcher
import darkmat13r.tarsos.dsp.io.NativeAudioInputStream
import org.koin.core.logger.Logger

actual class AudioManager {
    private var noiseSuppressor: NoiseSuppressor? = null

    private companion object {
        private val TAG = AudioManager::class.simpleName

    }

    actual fun createAudioDispatcher(): AudioDispatcher {
        val bufferSize = getBufferSize()
        val audioRecord = getAudioRecord(bufferSize).apply {
            startRecording()
        }

        if (NoiseSuppressor.isAvailable()) {
            startNoiseSuppressor(audioRecord.audioSessionId)
        }
        return getAudioDispatcher(audioRecord, bufferSize)
    }

    private fun getAudioDispatcher(audioRecord: AudioRecord, bufferSize: Int): AudioDispatcher {
        val format = darkmat13r.tarsos.dsp.io.AudioFormat(
            AudioConfig.SAMPLE_RATE.toFloat(),
            AudioConfig.SAMPLE_RATE_BITS,
            AudioConfig.CHANNEL_COUNT,
            signed = true,
            bigEndian = false
        )
        val audioStream = NativeAudioInputStream(audioRecord, format)

        return AudioDispatcher(audioStream, bufferSize, AudioConfig.OVERLAP)
    }

    actual fun stop() {
        noiseSuppressor?.apply {
            enabled = false
            release()
        }
    }

    private fun startNoiseSuppressor(audioSessionId: Int) {
        runCatching {
            noiseSuppressor = NoiseSuppressor.create(audioSessionId)
                .apply { enabled = true }
        }.onFailure(::logError)
    }

    private fun logError(throwable: Throwable) {
        throwable.printStackTrace()
    }

    @SuppressLint("MissingPermission")
    private fun getAudioRecord(bufferSize: Int) =
        AudioRecord(
            MediaRecorder.AudioSource.MIC,
            AudioConfig.SAMPLE_RATE,
            AudioFormat.CHANNEL_IN_MONO,
            AudioFormat.ENCODING_PCM_16BIT,
            bufferSize * 2
        )


    actual fun getBufferSize(): Int {
        val minBufferSize = AudioRecord.getMinBufferSize(
            AudioConfig.SAMPLE_RATE,
            AudioFormat.CHANNEL_IN_MONO,
            AudioFormat.ENCODING_PCM_16BIT
        )
        val minAudioBufferSizeInSamples = minBufferSize / 2

        return if (minAudioBufferSizeInSamples > AudioConfig.BUFFER_SIZE) minAudioBufferSizeInSamples else AudioConfig.BUFFER_SIZE
    }
}