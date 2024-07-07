package tuner

import androidx.lifecycle.Lifecycle
import androidx.lifecycle.LifecycleEventObserver
import androidx.lifecycle.LifecycleOwner
import androidx.lifecycle.lifecycleScope
import audio.AudioConfig
import audio.AudioManager
import coroutines.IODispatcher
import darkmat13r.tarsos.dsp.AudioDispatcher
import darkmat13r.tarsos.dsp.AudioEvent
import darkmat13r.tarsos.dsp.pitch.PitchProcessor
import darkmat13r.tarsos.dsp.pitch.PitchResult
import data.tuner.ChromaticScale
import data.tuner.Tuning
import data.tuner.TuningDeviationPrecision
import data.tuner.TuningDeviationResult
import kotlinx.coroutines.CoroutineScope
import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.flow.MutableStateFlow
import kotlinx.coroutines.flow.asStateFlow
import kotlinx.coroutines.launch
import kotlinx.coroutines.withContext
import kotlinx.coroutines.yield
import utils.logError
import kotlin.math.absoluteValue
import kotlin.math.log2
import kotlin.math.roundToInt


class Tuner(private val audioManager: AudioManager) :
    PitchProcessor.DetectedPitchHandler {

    private val _state by lazy { MutableStateFlow<Tuning>(Tuning()) }
    val state by lazy { _state.asStateFlow() }

    private var audioDispatcher: AudioDispatcher? = null
    private var pitchProcessor: PitchProcessor? = null


    suspend fun startRecording() {
        withContext(IODispatcher) {
            runCatching {
                val bufferSize = audioManager.getBufferSize()

                pitchProcessor = PitchProcessor(
                    PitchProcessor.PitchEstimationAlgorithm.YIN,
                    AudioConfig.SAMPLE_RATE.toFloat(),
                    bufferSize,
                    this@Tuner
                )
                audioDispatcher = audioManager.createAudioDispatcher().apply {
                    addAudioProcessor(pitchProcessor!!)
                    run()
                }
            }.onFailure(::logError)
        }
    }


    suspend fun stopRecording() {
        withContext(IODispatcher) {
            runCatching {
                audioManager.stop()
                audioDispatcher?.apply {
                    pitchProcessor?.let { removeAudioProcessor(it) }
                    stop()
                }
                pitchProcessor = null
                audioDispatcher = null
            }.onFailure(::logError)
        }
    }

    override fun handlePitch(pitch: PitchResult, audioEvent: AudioEvent) {
        if (pitch is PitchResult.Pitch)
            _state.value = getTuning(pitch.hz)
    }

    private fun getTuning(detectedFrequency: Float): Tuning {
        var minDeviation = Int.MAX_VALUE
        var closestNote = ChromaticScale.notes.first()

        ChromaticScale.notes.forEach { note ->
            val deviation = getTuningDeviation(note.frequency, detectedFrequency)
            if (deviation.absoluteValue < minDeviation.absoluteValue) {
                minDeviation = deviation
                closestNote = note
            }
        }

        val deviationResult = TuningDeviationResult.Detected(
            value = minDeviation,
            precision = TuningDeviationPrecision.fromDeviation(
                deviation = minDeviation,
                offset = 2
            )
        )

        return Tuning(closestNote, detectedFrequency, deviationResult)
    }

    private fun getTuningDeviation(standardFrequency: Float, detectedFrequency: Float) =
        try {
            (1200 * log2(detectedFrequency / standardFrequency)).roundToInt()
        } catch (ex: Exception) {
            0
        }
}

