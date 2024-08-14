package tuner

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
import kotlinx.coroutines.flow.MutableStateFlow
import kotlinx.coroutines.flow.asStateFlow
import kotlinx.coroutines.withContext
import kotlinx.datetime.Clock
import org.koin.core.component.KoinComponent
import utils.logError
import kotlin.math.absoluteValue
import kotlin.math.log2
import kotlin.math.roundToInt


class Tuner(private val audioManager: AudioManager) :
    PitchProcessor.DetectedPitchHandler, KoinComponent {


    private val _state by lazy { MutableStateFlow<Tuning>(Tuning()) }
    val state by lazy { _state.asStateFlow() }

    private var audioDispatcher: AudioDispatcher? = null
    private var pitchProcessor: PitchProcessor? = null


    suspend fun startRecording() {
        withContext(IODispatcher) {
            runCatching {
                val bufferSize = audioManager.getBufferSize()
                co.touchlab.kermit.Logger.i("Start Recording ${bufferSize}")
                pitchProcessor = PitchProcessor(
                    AudioConfig.SAMPLE_RATE.toFloat(),
                    bufferSize,
                    this@Tuner,
                    algorithm = PitchProcessor.PitchEstimationAlgorithm.YIN
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
        if (pitch is PitchResult.Pitch) {
            if (pitch.pitched)
                co.touchlab.kermit.Logger.i("pitch ${pitch} ${Clock.System.now()}")
            getTuning(pitch.hz)?.also {
                _state.value = it
            }
        }
    }

    private fun getTuning(detectedFrequency: Float): Tuning? {
        var minDeviation = Int.MAX_VALUE
        var closestNote = ChromaticScale.notes.first()

        ChromaticScale.notes.forEach { note ->
            val deviation = getTuningDeviation(note.frequency, detectedFrequency)
            if (deviation.absoluteValue < minDeviation.absoluteValue) {
                minDeviation = deviation
                closestNote = note
            }
        }
        if (minDeviation == 0) return null

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

